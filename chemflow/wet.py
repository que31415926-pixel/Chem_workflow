"""Produce only the next actionable comparison. Evidence flags are reviewed by a chemist."""

from .network import ROOT, read_json, write_table


def handoff(project, state_path, dest):
    plan = read_json(ROOT / "projects" / project / "wet_plan.json")
    state = read_json(state_path or (ROOT / "projects" / project / "state.json"))
    fields = plan["fields"]
    if set(state) != set(fields) or any(
        v not in ["unknown", *fields[k]] for k, v in state.items()
    ):
        raise ValueError("Invalid state keys/values; see wet_plan.json")
    stage = None
    action = None
    if state["controls"] == "fail":
        action = "修复阳性、标品或回收，重新验证装置及同批对照。"
    elif state["controls"] == "unknown":
        stage = plan["stages"][0]
    elif project == "fe" and state["salt"] in ("unknown", "unusable"):
        action = "先定量实际可递送的盐，核对暗置损失并检查有同骨架一级对照的改环前体。原盐不可投足时铁端记为不可评价。"
    elif project == "co" and state["balance"] == "fail":
        action = "补查液相/气相校准、硼和氢收支、未归属峰与回收，再作路线判断。"
    else:
        outcome = state["parent"] if project == "co" else state["target"]
        if outcome == "uncertain":
            action = "补独立日配对或改进定量限，保留原始区间。"
        elif outcome in ("works", "meets"):
            if project == "co" and state["balance"] != "pass":
                action = "先闭合氢/硼收支并复核目标产物和候选基准。"
            elif state["heldout"] == "unknown":
                action = "暂停开发新配体，在未参与选择的一种真实底物或条件上冻结预测并同批验证。"
            elif state["heldout"] == "fail":
                action = "限定结论到已验证底物，记录留出失败的适用边界。"
            else:
                action = "整理同批对照、独立验证与不确定度，结束这一轮。"
        elif (
            project == "co" and outcome == "shortfall" and state["balance"] == "unknown"
        ):
            action = "在选干预前完成目标格的芳基硼酸酯、硼、氢及空白收支。"
        else:
            eligible = [
                s
                for s in plan["stages"]
                if all(state[k] == v for k, v in s["when"].items())
            ]
            stage = eligible[-1]
    dest.mkdir(parents=True, exist_ok=True)
    if stage:
        rows = [
            dict(phase=stage["id"], question=stage["question"], **r)
            for r in stage["rows"]
        ]
        if project == "fe" and state["salt"] == "modified_only":
            for r in rows:
                r["comparison"] += "；本轮仅改环系列。原三苯环的目标配对记不可评价。"
    else:
        rows = [
            dict(
                phase="review",
                question="根据实际证据更新判断",
                arm=action,
                readout="记录证据文件、批号、定量方式和局限",
                comparison="同批参照",
            )
        ]
    write_table(dest / "next_batch.csv", rows[0].keys(), rows)
    cols = [
        "project",
        "arm",
        "batch",
        "vial",
        "time_h",
        "light_or_headspace",
        "salt_scaffold_or_catalyst",
        "analyte",
        "value",
        "unit",
        "method",
        "LOD",
        "LOQ",
        "recovery",
        "qc",
        "evidence_file",
        "note",
    ]
    write_table(dest / "wet_records_blank.csv", cols, [])
    return rows


def paired(path, delta=15.0, practical=30.0):
    """Bound paired means / t intervals using independent-day lower and upper limits.
    For censored intervals use bounds on sample SD, not substitution by zero/LOQ.
    """
    import csv
    import math

    from scipy.stats import t

    with open(path, encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    days = set()
    measures = {"target": [], "reference_minus_target": []}
    for r in rows:
        if r["qc"] != "pass":
            raise ValueError(
                "Resolve failed/unreviewed QC before using paired endpoints"
            )
        day = r["independent_block"].strip()
        if not day or day in days:
            raise ValueError("Each independently paired day/block appears once")
        days.add(day)
        tlo, thi, rlo, rhi = [
            float(r[k])
            for k in ("target_lo_pct", "target_hi_pct", "ref_lo_pct", "ref_hi_pct")
        ]
        if not (0 <= tlo <= thi <= 100 and 0 <= rlo <= rhi <= 100):
            raise ValueError("Invalid yield interval")
        measures["target"].append((tlo, thi))
        measures["reference_minus_target"].append((rlo - thi, rhi - tlo))
    n = len(days)
    if n < 3:
        raise ValueError(
            "Project decision rule requires at least three independent matched blocks"
        )
    result = {"n": n, "practical_target_pct": practical, "paired_gap_points": delta}
    crit = t.ppf(0.975, n - 1)
    if not (0 <= practical <= 100 and 0 <= delta <= 100):
        raise ValueError("Thresholds must be in 0..100")
    for name, values in measures.items():
        # Identity: sample variance = sum_{i<j}(x_i-x_j)^2 / (n*(n-1)).
        d2 = sum(
            max(abs(a - d), abs(b - c)) ** 2
            for i, (a, b) in enumerate(values)
            for c, d in values[i + 1 :]
        )
        half = crit * math.sqrt(d2 / (n * n * (n - 1)))
        low = sum(a for a, b in values) / n
        high = sum(b for a, b in values) / n
        result[name] = {
            "mean_bounds": [low, high],
            "conservative_95pct_t_envelope": [low - half, high + half],
        }
    lo, hi = result["target"]["conservative_95pct_t_envelope"]
    gap = result["reference_minus_target"]["conservative_95pct_t_envelope"][0]
    result["target_state"] = (
        "works"
        if lo >= practical
        else ("low" if hi < practical and gap > delta else "uncertain")
    )
    return result
