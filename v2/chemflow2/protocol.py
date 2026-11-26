"""Load a ChemFlow protocol and match its tasks to declared node capabilities."""

import hashlib
import json
import re
from collections import deque
from pathlib import Path


_IDENTIFIER = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
_TRANSPORT_FILES = {"chemflow-task.json", "chemflow-result.json", "__chemflow_transfer__.json",
                    "chemflow-worker.out.log", "chemflow-worker.err.log", "protocol.json"}
_RESOURCES = ("cores", "memory_mb", "scratch_mb")


def _mapping(value, label):
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    return value


def _string(value, label):
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise ValueError(f"{label} must be a nonempty string without surrounding whitespace")
    return value


def _identifier(value, label):
    if not _IDENTIFIER.fullmatch(_string(value, label)):
        raise ValueError(f"invalid {label}: {value!r}")
    return value


def _strings(value, label):
    if not isinstance(value, list):
        raise ValueError(f"{label} must be a list")
    seen = set()
    for entry in value:
        _string(entry, label)
        if entry in seen:
            raise ValueError(f"duplicate {label}: {entry}")
        seen.add(entry)
    return value


def _resources(value, label):
    _mapping(value, label)
    for key in _RESOURCES:
        quantity = value.get(key)
        minimum = 0 if key == "scratch_mb" else 1
        if type(quantity) is not int or quantity < minimum:
            raise ValueError(f"{label}.{key} must be an integer >= {minimum}")


def _relative(value, label):
    _string(value, label)
    parts = value.split("/")
    if (any(p in ("", ".", "..") or p.startswith("-") for p in parts)
            or "\\" in value or ":" in value or any(ord(c) < 32 or ord(c) == 127 for c in value)):
        raise ValueError(f"unsafe relative path in {label}: {value!r}")
    first = parts[0].casefold()
    if first in _TRANSPORT_FILES or (first == "orca.out" and (len(parts) > 1 or value != "orca.out")):
        raise ValueError(f"reserved path in {label}: {value!r}")
    return value


def _static_file(root, filename, label):
    path = root / filename
    if not path.is_file() or path.resolve() != path:
        raise ValueError(f"{label} is missing, linked, or outside the protocol directory: {filename}")
    return path


def _digest(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _staging_conflicts(paths, label):
    files = set()
    directories = set()
    for path in paths:
        # A Mac client may stage on a case-insensitive filesystem.
        parts = path.casefold().split("/")
        parents = {"/".join(parts[:i]) for i in range(1, len(parts))}
        normalized = "/".join(parts)
        if normalized in files or normalized in directories or parents & files:
            raise ValueError(f"staging path collision in {label}: {path}")
        files.add(normalized)
        directories.update(parents)


def _unique_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _order(tasks):
    by_id = {task["id"]: task for task in tasks}
    counts = {key: len(task["deps"]) for key, task in by_id.items()}
    children = {key: [] for key in by_id}
    for task in tasks:
        for dep in task["deps"]:
            if dep not in by_id:
                raise ValueError(f"{task['id']} depends on unknown task {dep}")
            children[dep].append(task["id"])
    ready = deque(key for key, count in counts.items() if count == 0)
    ordered = []
    while ready:
        key = ready.popleft()
        ordered.append(key)
        for child in children[key]:
            counts[child] -= 1
            if counts[child] == 0:
                ready.append(child)
    if len(ordered) != len(tasks):
        raise ValueError("cycle in protocol dependencies")
    return ordered


def static_sources(task):
    """Map each task-local input name to its source within the frozen bundle."""
    names = list(task.get("files", []))
    if task.get("input") is not None:
        names.append(task["input"])
    aliases = _mapping(task.get("static_sources", {}), "static_sources")
    if set(aliases) - set(names):
        raise ValueError("static_sources must refer to declared task-local inputs")
    return {name: _relative(aliases.get(name, name), "static source") for name in names}


def load_plan(path: Path) -> dict:
    path = Path(path)
    if path.is_symlink():
        raise ValueError("protocol file must not be a symlink")
    path = path.resolve()
    plan = _mapping(json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_unique_keys), "plan")
    _string(plan.get("name"), "plan.name")
    tasks = plan.get("tasks")
    if not isinstance(tasks, list) or not tasks:
        raise ValueError("plan.tasks must be a nonempty list")

    ids = set()
    all_static = set()
    for task in tasks:
        _mapping(task, "task")
        key = _identifier(task.get("id"), "task id")
        if key.casefold() in ids:
            raise ValueError(f"duplicate task id: {key}")
        ids.add(key.casefold())
        if type(task.get("priority", 0)) is not int:
            raise ValueError(f"{key}.priority must be an integer")
        engine = _identifier(task.get("engine"), f"{key}.engine")
        stage = _identifier(task.get("stage"), f"{key}.stage")
        _resources(task.get("resources"), f"{key}.resources")
        _strings(task.setdefault("requires", []), f"{key}.requires")
        _strings(task.setdefault("deps", []), f"{key}.deps")
        files = _strings(task.setdefault("files", []), f"{key}.files")
        outputs = _strings(task.get("outputs"), f"{key}.outputs")
        if not outputs:
            raise ValueError(f"{key}.outputs must declare at least one output")
        inherited = _mapping(task.setdefault("from_deps", {}), f"{key}.from_deps")
        for name in inherited:
            _relative(name, f"{key}.from_deps")
        for name in outputs:
            _relative(name, f"{key}.outputs")

        static = list(files)
        if task.get("input") is not None:
            static.append(_relative(task["input"], f"{key}.input"))
        elif not (stage == "probe" and engine in {"probe", "command"}):
            raise ValueError(f"{key}.input must specify the primary input file")
        sources = static_sources(task)
        for name in static:
            _relative(name, f"{key}.inputs")
            _static_file(path.parent, sources[name], key)
        all_static.update(sources.values())
        staged_inputs = [*static, *inherited]
        if "orca.out" in staged_inputs:
            raise ValueError(f"{key}: orca.out is reserved for engine output")
        _staging_conflicts([*staged_inputs, *outputs], key)

    attachments = _strings(plan.get("attachments", []), "attachments")
    _staging_conflicts(sorted(all_static | set(attachments)), "protocol inputs and attachments")
    hashes = {}
    for name in attachments:
        _relative(name, "attachment")
        hashes[name] = _digest(_static_file(path.parent, name, "attachment"))
    if "attachment_sha256" in plan and hashes != _mapping(plan["attachment_sha256"], "attachment_sha256"):
        raise ValueError("attachment SHA-256 inventory does not match the files")
    if attachments:
        plan["attachment_sha256"] = hashes
    _order(tasks)
    by_id = {task["id"]: task for task in tasks}
    for task in tasks:
        for dest, source in task["from_deps"].items():
            dep, separator, output = _string(source, "from_deps source").partition("/")
            if not separator or dep not in task["deps"] or output not in by_id[dep]["outputs"]:
                raise ValueError(f"{task['id']}: {dest} refers to undeclared dependency output {source!r}")
    return plan


def validate_node(node: dict) -> dict:
    _mapping(node, "node")
    name = _identifier(node.get("name"), "node.name")
    _identifier(node.get("host"), f"{name}.host")
    _resources(node, name)
    _strings(node.get("engines"), f"{name}.engines")
    _strings(node.get("capabilities", []), f"{name}.capabilities")
    concurrency = node.get("concurrency", 1)
    if type(concurrency) is not int or concurrency < 1:
        raise ValueError(f"{name}.concurrency must be a positive integer")
    return node


def _mismatches(task, node):
    reasons = []
    if task["engine"] not in node["engines"]:
        reasons.append(f"engine {task['engine']} unavailable")
    for capability in task.get("requires", []):
        if capability not in node.get("capabilities", []):
            reasons.append(f"missing capability {capability}")
    for resource in _RESOURCES:
        needed = task["resources"][resource]
        if needed > node[resource]:
            reasons.append(f"requires {needed} {resource}, node has {node[resource]}")
    return reasons


def matches(task: dict, node: dict) -> bool:
    """Check one task against one node's total capacity, without using a substitute engine."""
    validate_node(node)
    return not _mismatches(task, node)


def explain_plan(plan: dict, nodes: list) -> dict:
    """Report static fit and propagate a missing dependency's effect down the graph."""
    for node in nodes:
        validate_node(node)
    _strings([n["name"] for n in nodes], "node names")
    tasks = plan["tasks"]
    by_id = {task["id"]: task for task in tasks}
    rows = {}
    blocked = set()
    for key in _order(tasks):
        task = by_id[key]
        blockers = {n["name"]: _mismatches(task, n) for n in nodes}
        eligible = [name for name, reasons in blockers.items() if not reasons]
        blocked_deps = [dep for dep in task["deps"] if dep in blocked]
        if not eligible or blocked_deps:
            blocked.add(key)
        rows[key] = {"id": key, "engine": task["engine"], "stage": task["stage"],
                     "eligible_nodes": eligible,
                     "blockers": {name: reasons for name, reasons in blockers.items() if reasons},
                     "blocked_by_dependencies": blocked_deps}
    return {"name": plan["name"], "task_count": len(tasks), "runnable": len(tasks) - len(blocked),
            "blocked": [key for key in by_id if key in blocked],
            "tasks": [rows[key] for key in by_id]}
