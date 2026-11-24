"""Resolve project directories and their versioned, data-only presets."""

import json
import math
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IDENT = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z")


def read_json(path):
    with Path(path).open(encoding="utf-8-sig") as f:
        try:
            return json.load(f)
        except json.JSONDecodeError as exc:
            raise ValueError(f"{path}: {exc}") from exc


def write_json(path, value):
    Path(path).write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def projects_dir():
    return (
        Path(os.environ.get("CHEMFLOW_PROJECTS_DIR", ROOT / "projects"))
        .expanduser()
        .resolve()
    )


def identifier(value):
    if not isinstance(value, str) or not IDENT.fullmatch(value):
        raise ValueError(
            f"Invalid identifier {value!r}; use letters, digits and underscores, starting with a letter or underscore"
        )
    return value


def folder(project):
    return projects_dir() / identifier(project)


def discover():
    return sorted(
        p.parent.name
        for p in projects_dir().glob("*/project.json")
        if IDENT.fullmatch(p.parent.name)
    )


# Short entries preserve the original output layout. Objects name individual steps.
STEP_FIELDS = {
    "build-xyz": {"species", "seed"},
    "simulate": {"calibrated"},
    "handoff": {"state"},
    "prepare": {"jobs"},
    "fit": {"data"},
    "balance": {"data"},
    "conformers": {"species", "cores", "quick", "exe", "xyz", "frame"},
    "run-orca": {"job", "exe"},
    "audit": {"job"},
    "thermo": {"recipe"},
}
DEFAULT_OUTPUT = {
    "simulate": "scenarios",
    "handoff": "wet",
    "prepare": "orca_input",
    "run-orca": "orca_run",
    "build-xyz": "structures",
}
REQUIRED_INPUT = {
    "fit": "data",
    "balance": "data",
    "conformers": "species",
    "run-orca": "job",
    "audit": "job",
    "thermo": "recipe",
}


def steps(info):
    records = info.get("workflow")
    if not isinstance(records, list) or not records:
        raise ValueError("workflow must be a nonempty ordered list")
    seen = set()
    for item in records:
        row = {"step": item} if isinstance(item, str) else item
        if not isinstance(row, dict):
            raise ValueError("Workflow steps are strings or objects")
        kind = row.get("step")
        if kind not in STEP_FIELDS or set(row) - {"step", "id", *STEP_FIELDS[kind]}:
            raise ValueError(f"Unknown workflow step or option: {row!r}")
        name = identifier(row.get("id", DEFAULT_OUTPUT.get(kind, kind)))
        if name in seen:
            raise ValueError(f"Output id reused: {name}")
        if kind in REQUIRED_INPUT and not row.get(REQUIRED_INPUT[kind]):
            raise ValueError(f"{name}: missing {REQUIRED_INPUT[kind]}")
        for key in ("data", "state", "calibrated", "job", "recipe", "xyz"):
            if key in row:
                check_input(row[key], seen)
        if "jobs" in row:
            if not isinstance(row["jobs"], list) or not row["jobs"]:
                raise ValueError(f"{name}: jobs must be a nonempty list")
            for j in row["jobs"]:
                if j.get("xyz"):
                    check_input(j["xyz"], seen)
        seen.add(name)
        yield {**row, "step": kind, "id": name}


def check_input(ref, seen):
    if not isinstance(ref, str) or not ref or "\\" in ref:
        raise ValueError(f"Use a project-relative path or @earlier_step/file: {ref!r}")
    if ref.startswith("@"):
        label, _, tail = ref[1:].partition("/")
        if label not in seen:
            raise ValueError(f"Reference to missing or later step: {ref}")
        path = Path(tail)
    else:
        path = Path(ref)
    if path.is_absolute() or path.drive or ".." in path.parts:
        raise ValueError(
            f"Input path must stay in its project or step directory: {ref}"
        )


def manifest(project):
    path = folder(project) / "project.json"
    if not path.is_file():
        raise ValueError(
            f"Unknown project {project!r} under {projects_dir()}. See 'projects' or 'init'."
        )
    info = read_json(path)
    if info.get("schema_version") != 1:
        raise ValueError(f"{path}: unsupported schema version")
    list(steps(info))
    return info


def finite(value):
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(value)
    )
