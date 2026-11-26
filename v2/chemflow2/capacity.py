"""Estimate one frozen plan from verified, explicitly comparable pilot receipts."""

import hashlib
import json
import math
from pathlib import Path
import statistics

from .controller import load_nodes
from .protocol import _mismatches, _order, load_plan


def _number(value):
    return type(value) in (float, int) and math.isfinite(value) and value >= 0


def _history(path, task_ids):
    data = json.loads(path.read_bytes())
    samples = data["samples"]
    versions = data.get("engine_versions", {})
    if not isinstance(samples, list) or not isinstance(versions, dict):
        raise ValueError("Expected history.samples array and optional engine_versions object")
    result = []
    seen = set()
    for item in samples:
        source = (path.parent / item["receipt"]).resolve()
        if not source.is_relative_to(path.parent) or Path(item["receipt"]).is_absolute():
            raise ValueError("Receipt must be below the history directory")
        raw = source.read_bytes()
        checksum = hashlib.sha256(raw).hexdigest()
        if checksum != item["sha256"] or checksum in seen:
            raise ValueError("Receipt digest mismatch or duplicate receipt")
        seen.add(checksum)
        selection = item.get("comparable_to", [])
        if (not isinstance(selection, list) or not all(isinstance(s, str) and s in task_ids for s in selection)
                or selection and not (isinstance(item.get("rationale"), str) and item["rationale"].strip())):
            raise ValueError("Comparable target ids require a written rationale")
        receipt = json.loads(raw)
        if not isinstance(receipt, dict):
            raise ValueError("Expected a worker receipt object")
        result.append((item, receipt))
    return result, versions


def _comparable(task, node, source_node, receipt, review, version):
    if not source_node or not isinstance(version, str) or not version:
        return False
    recorded = receipt.get("node_capacity_class") or review.get("capacity_class")
    if node.get("capacity_class"):
        if recorded != node["capacity_class"]:
            return False
    elif source_node["name"] != node["name"] or review.get("same_hardware_confirmed") is not True:
        return False
    provenance = receipt.get("provenance", {})
    actual_version = receipt.get("observations", {}).get("engine_version") or provenance.get("engine_version")
    if (provenance.get("engine"), provenance.get("stage"), actual_version) != (
            task["engine"], task["stage"], version):
        return False
    if node.get("arch") and provenance.get("architecture") != node["arch"]:
        return False
    if provenance.get("task_id") != task["id"] and task["id"] not in review.get("comparable_to", []):
        return False
    reservations = provenance.get("reservations", {})
    cores = reservations.get("cores")
    metrics = receipt.get("metrics", {})
    wall = metrics.get("wall_seconds", provenance.get("wall_seconds"))
    if not _number(wall) or wall <= 0:
        return False
    if cores is None and _number(metrics.get("requested_core_seconds")):
        # Stored values are rounded to microseconds; avoid guessing for very short jobs.
        return wall >= .001 and math.isclose(metrics["requested_core_seconds"], wall * task["resources"]["cores"],
                                               rel_tol=1e-5, abs_tol=2e-6)
    return cores == task["resources"]["cores"]


def _peaks(samples, field):
    values = [r.get("metrics", {}).get(field) for _, r in samples]
    values = [v for v in values if _number(v)]
    return max(values, default=None)


def _profile(task, node, samples, nodes, version):
    reasons = _mismatches(task, node)
    comparable = [(review, receipt) for review, receipt in samples
                  if _comparable(task, node, nodes.get(receipt.get("node")), receipt, review, version)]
    passed = [(review, r) for review, r in comparable if r.get("ok") is True and r.get("status") == "succeeded"]
    walls = [r.get("metrics", {}).get("wall_seconds", r["provenance"].get("wall_seconds")) for _, r in passed]
    observations = {key: _peaks(comparable, key) for key in
                    ("job_peak_memory_bytes", "child_max_rss_bytes", "scratch_peak_bytes",
                     "scratch_regular_file_bytes_at_finish", "scratch_allocated_sampled_max_bytes")}
    deficits = {k: max(0, task["resources"][k] - node[k]) for k in ("cores", "memory_mb", "scratch_mb")}
    for quantity, fields in (("memory_mb", ("job_peak_memory_bytes", "child_max_rss_bytes")),
                             ("scratch_mb", ("scratch_peak_bytes", "scratch_regular_file_bytes_at_finish", "scratch_allocated_sampled_max_bytes"))):
        bound = max((observations[k] for k in fields if observations[k] is not None), default=None)
        if bound is not None and bound > task["resources"][quantity] * 2**20:
            lower = math.ceil(bound/2**20)
            deficits[quantity] = max(deficits[quantity], lower - node[quantity])
            reasons.append(f"observed {quantity} lower bound {lower} exceeds task reservation")
    if not version:
        reasons.append("target engine version unspecified")
    if not passed:
        reasons.append("no comparable successful receipt")
    return {"node": node["name"], "deficits": deficits, "reasons": reasons,
            "matched_successful": len(passed), "matched_unsuccessful": len(comparable)-len(passed),
            "observed_envelopes": observations, "sample_sha256": [i["sha256"] for i, _ in comparable],
            "median_success_wall_seconds": statistics.median(walls) if walls else None}


def estimate(plan_path, history_path, nodes_path):
    """Return a conditional throughput lower bound, never a completion or spending promise."""
    plan = load_plan(Path(plan_path))
    nodes = load_nodes(Path(nodes_path))
    by_node = {n["name"]: n for n in nodes}
    by_task = {t["id"]: t for t in plan["tasks"]}
    samples, versions = _history(Path(history_path).resolve(), set(by_task))
    rows = {}
    finishes = {}
    allocation = {}
    for node in nodes:
        remaining = node.get("remaining_reserved_core_hours")
        if remaining is not None and not _number(remaining):
            raise ValueError("remaining_reserved_core_hours must be finite and nonnegative")
        allocation[node["name"]] = {"assigned_tasks": 0, "assigned_wall_seconds": 0.,
                                    "assigned_reserved_core_hours": 0., "remaining_reserved_core_hours": remaining}
    for key in _order(plan["tasks"]):
        task = by_task[key]
        version = task.get("engine_version") or versions.get(task["engine"])
        profiles = [_profile(task, n, samples, by_node, version) for n in nodes]
        viable = [p for p in profiles if not p["reasons"]]
        choice = min(viable, key=lambda p: p["median_success_wall_seconds"]) if viable else None
        rows[key] = {"engine": task["engine"], "stage": task["stage"], "estimated_engine_version": version,
                     "assigned_node": choice["node"] if choice else None, "profiles": profiles}
        if choice:
            name = choice["node"]
            duration = choice["median_success_wall_seconds"]
            allocation[name]["assigned_tasks"] += 1
            allocation[name]["assigned_wall_seconds"] += duration
            allocation[name]["assigned_reserved_core_hours"] += duration * task["resources"]["cores"] / 3600
            if all(dep in finishes for dep in task["deps"]):
                finishes[key] = duration + max((finishes[dep] for dep in task["deps"]), default=0.)
    unknown = [k for k, row in rows.items() if row["assigned_node"] is None]
    limits = []
    for name, summary in allocation.items():
        hours = summary["assigned_reserved_core_hours"]
        remaining = summary["remaining_reserved_core_hours"]
        summary["assigned_quota_shortfall_core_hours"] = None if remaining is None else max(0., hours - remaining)
        node = by_node[name]
        limits.extend([3600*hours / node["cores"], summary["assigned_wall_seconds"]/node.get("concurrency", 1)])
    critical = None if unknown else max(finishes.values())
    return {"study": plan["name"], "status": "incomplete" if unknown else "estimated",
            "unknown_tasks": unknown, "task_profiles": rows, "node_allocation": allocation,
            "critical_path_seconds": critical,
            "conditional_min_elapsed_seconds": None if unknown else max([critical, *limits]),
            "total_reserved_core_hours": None if unknown else sum(x["assigned_reserved_core_hours"] for x in allocation.values()),
            "queue_wait_seconds": None, "actual_cost": None,
            "assumptions": [
                "Fastest measured compatible median per task, fixed allocation, simultaneous node availability.",
                "Bound applies to this assignment and these durations. Contention, retries, transfer, queueing and pauses add time.",
                "Node quota is in operator-confirmed reserved-core hours for the same period. Unassigned tasks are excluded from node totals.",
                "Child RSS and scratch at finish are observed lower bounds. Null aggregate peaks are unknown.",
                "Receipts are user-supplied measurements; their hashes bind the files, not the operator's assertion of comparability."]}
