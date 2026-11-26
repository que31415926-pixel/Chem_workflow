"""Select the next wet-lab batch using ordered rules declared by each project."""

from .network import write_table
from .projects import folder, read_json


def checked_plan(project):
    plan = read_json(folder(project) / "wet_plan.json")
    fields = plan["fields"]
    stages = {s["id"]: s for s in plan["stages"]}
    if len(stages) != len(plan["stages"]):
        raise ValueError("Duplicate wet stage IDs")
    for stage in stages.values():
        if (
            not stage.get("question")
            or not stage.get("rows")
            or any(set(r) != {"arm", "readout", "comparison"} for r in stage["rows"])
        ):
            raise ValueError(f"Incomplete wet stage: {stage['id']}")
    rules = plan["rules"]
    if not rules or rules[-1]["when"] or any(not r["when"] for r in rules[:-1]):
        raise ValueError(
            "Wet rules use first match; put a single unconditional fallback last"
        )
    for rule in rules:
        if ("stage" in rule) == ("review" in rule) or (
            "stage" in rule and rule["stage"] not in stages
        ):
            raise ValueError("Each wet rule requires a known stage or review message")
    for entry in [*rules, *plan.get("notes", [])]:
        for k, v in entry["when"].items():
            vals = v if isinstance(v, list) else [v]
            if k not in fields or not vals or set(vals) - {"unknown", *fields[k]}:
                raise ValueError(f"Invalid wet condition {k}={v!r}")
    return plan, stages


def handoff(project, state_path, dest):
    plan, stages = checked_plan(project)
    state = read_json(state_path or folder(project) / "state.json")
    fields = plan["fields"]
    if set(state) != set(fields) or any(
        v not in ["unknown", *fields[k]] for k, v in state.items()
    ):
        raise ValueError("Invalid state keys/values; see wet_plan.json")

    def matches(when):
        return all(
            state[k] in (v if isinstance(v, list) else [v]) for k, v in when.items()
        )

    selected = next(rule for rule in plan["rules"] if matches(rule["when"]))
    if "stage" in selected:
        stage = stages[selected["stage"]]
        rows = [
            dict(phase=stage["id"], question=stage["question"], **r)
            for r in stage["rows"]
        ]
        suffix = "".join(
            n["comparison"] for n in plan.get("notes", []) if matches(n["when"])
        )
        for row in rows:
            row["comparison"] += suffix
    else:
        rows = [
            dict(
                phase="review",
                question="根据实际证据更新判断",
                arm=selected["review"],
                readout="记录证据文件、批号、定量方式和局限",
                comparison="同批参照",
            )
        ]
    dest.mkdir(parents=True, exist_ok=True)
    write_table(dest / "next_batch.csv", rows[0].keys(), rows)
    write_table(
        dest / "wet_records_blank.csv",
        plan.get(
            "record_columns",
            [
                "batch",
                "arm",
                "time_h",
                "analyte",
                "value",
                "unit",
                "qc",
                "evidence_file",
            ],
        ),
        [],
    )
    return rows


def paired(path, delta=15.0, practical=30.0, min_blocks=3):
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
    if min_blocks < 2 or n < min_blocks:
        raise ValueError(
            f"Require at least {max(2, min_blocks)} independent matched blocks"
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
