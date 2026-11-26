"""Staged scientific analyses with explicit dependency bindings and evidence receipts."""

import hashlib
import json
from pathlib import Path

from . import mechanism, science
from .protocol import _relative, _unique_keys


def _reject_constant(value):
    raise ValueError(f"invalid JSON number: {value}")


def _json_path(folder, name):
    name = _relative(name, "analysis filename")
    path = folder / name
    if path.suffix != ".json" or path.resolve() != path:
        raise ValueError(f"analysis requires a task-local, non-symlink JSON path: {name!r}")
    return path


def _read(path):
    raw = path.read_bytes()
    return (json.loads(raw, object_pairs_hook=_unique_keys, parse_constant=_reject_constant),
            hashlib.sha256(raw).hexdigest())


def _bindings(value, inherited, folder, cache, evidence):
    if isinstance(value, list):
        return [_bindings(v, inherited, folder, cache, evidence) for v in value]
    if not isinstance(value, dict):
        return value
    if "$artifact" not in value:
        return {k: _bindings(v, inherited, folder, cache, evidence) for k, v in value.items()}
    if set(value) != {"$artifact", "path"} or not isinstance(value["path"], list):
        raise ValueError("a binding must contain exactly $artifact and a path array")
    name = value["$artifact"]
    if not isinstance(name, str) or name not in inherited:
        raise ValueError(f"undeclared dependency binding: {name!r}")
    if name not in cache:
        data, digest = _read(_json_path(folder, name))
        cache[name] = data
        evidence[name] = {"source": inherited[name], "sha256": digest}
    item = cache[name]
    for index in value["path"]:
        if isinstance(item, dict) and type(index) is str and index in item:
            item = item[index]
        elif isinstance(item, list) and type(index) is int and 0 <= index < len(item):
            item = item[index]
        else:
            raise ValueError(f"{name}: no JSON member at path {value['path']!r}")
    # A value taken from a dependency is data, never another executable binding.
    return item


def _disposition(stage, report):
    operation = report["operation"]
    if stage == "science":
        status = report["status"]
        if status == "pass":
            return "pass" if operation in ("gate", "direction") else "not_applicable"
        if status not in ("fail", "missing"):
            raise ValueError(f"unknown science status: {status}")
        return status
    if operation in {"fit_anchor_threshold", "score_anchor_threshold"}:
        return report["status"]
    if operation == "interval_kinetics":
        return "pass" if report["status"] == "compatible" else "indeterminate"
    if operation in ("isotope_observables", "eyring", "freeze_prediction", "solution_cycle"):
        return "not_applicable"
    if operation in ("irc_endpoints", "salt_quantity"):
        return report["status"]
    if operation == "paired_barriers":
        comparison = {r["status"] for r in report["comparisons"].values()}
        paths = {r["status"] for r in report["paths"].values()}
        if "disagreement" in comparison:
            return "fail"
        if "missing" in paths:
            return "missing"
        if "unqualified" in paths or comparison & {"within_band", "insufficient_evidence"}:
            return "indeterminate"
        if comparison == {"agreement"} and paths == {"qualified"}:
            return "pass"
    if operation == "score_holdout":
        statuses = {r["status"] for r in report["targets"].values()}
        if "outside" in statuses:
            return "fail"
        if "awaiting_observation" in statuses:
            return "missing"
        if "overlap" in statuses:
            return "indeterminate"
        if statuses == {"inside"}:
            return "pass"
    raise ValueError(f"unrecognized result status for {operation}")


def run_analysis(task, folder):
    """Write a declared report, including when evidence falsifies a predeclared test."""
    folder = Path(folder).absolute()
    if folder.resolve() != folder or not folder.is_dir():
        raise ValueError("analysis requires a private staged directory")
    stage = task.get("stage")
    if task.get("engine") != "analysis" or stage not in ("science", "mechanism"):
        raise ValueError("analysis supports only science and mechanism stages")
    outputs = task.get("outputs")
    if not isinstance(outputs, list) or len(outputs) != 1:
        raise ValueError("analysis requires exactly one declared JSON output")
    output = _json_path(folder, outputs[0])
    primary = _json_path(folder, task.get("input"))
    inherited = task.get("from_deps", {})
    if not isinstance(inherited, dict) or output.exists() or primary.name == "" or output == primary:
        raise ValueError("invalid analysis input or output declarations")
    template, primary_hash = _read(primary)
    evidence = {}
    data = _bindings(template, inherited, folder, {}, evidence)
    evaluator = science.evaluate if stage == "science" else mechanism.evaluate
    report = evaluator(data)
    disposition = _disposition(stage, report)
    resolved = json.dumps(data, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)
    receipt = {"primary_sha256": primary_hash, "dependency_artifacts": evidence,
               "resolved_input_sha256": hashlib.sha256(resolved.encode()).hexdigest(),
               "scientific_gate": disposition}
    report["analysis_receipt"] = receipt
    payload = json.dumps(report, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as stream:
        stream.write(payload)
    observations = {"type": "analysis", "operation": report["operation"], "report": outputs[0], **receipt}
    failure = None if disposition in ("pass", "not_applicable") else "scientific_" + disposition
    return observations, failure
