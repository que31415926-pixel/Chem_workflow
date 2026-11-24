"""Execute the ordered project stages. @id/file links a consumer to earlier output."""

import argparse
import shutil
import tempfile
from collections import Counter
from pathlib import Path

from .analysis import fit
from .material import account, specification
from .network import Network, run_scenarios
from .projects import (
    check_input,
    folder,
    identifier,
    manifest,
    read_json,
    steps,
    write_json,
)
from .qm import (
    ELEMENTS,
    conformers,
    geometry,
    get_species,
    orca_audit,
    prepare,
    qm_settings,
    run_orca,
)
from .structures import build, smiles_molecule
from .thermo import reaction_energy
from .wet import handoff


def resolve(ref, project, outputs):
    check_input(ref, outputs)
    if ref.startswith("@"):
        name, _, tail = ref[1:].partition("/")
        return outputs[name] / tail
    return folder(project) / ref


def jobs(project, stage, outputs, check_files):
    cfg = qm_settings(project)
    tasks = stage.get("jobs", cfg.get("jobs", []))
    if not tasks:
        raise ValueError("prepare requires at least one job")
    seen = set()
    for i, job in enumerate(tasks):
        label = identifier(job.get("id", f"job_{i + 1}"))
        if label in seen:
            raise ValueError(f"Duplicate job name: {label}")
        seen.add(label)
        settings = {
            "xyz": None,
            "frame": 1,
            "mult": None,
            "cores": 1,
            "memory": 1500,
            "stage": "min",
            "method": "screen",
            "solvation": "default",
            "temperature": None,
            "reference_solvent": None,
        }
        if set(job) - {"id", "species", *settings}:
            raise ValueError(f"Unsupported setting in {label}")
        settings.update({k: v for k, v in job.items() if k != "id"})
        if (
            settings["stage"] not in {"min", "ts", "sp", "irc"}
            or "species" not in settings
        ):
            raise ValueError(f"Invalid QM stage or missing species in {label}")
        rec = get_species(project, settings["species"])
        if (
            settings["mult"] is not None
            and settings["mult"] not in rec["multiplicities"]
        ):
            raise ValueError(f"{label}: unregistered spin multiplicity")
        deferred = False
        if settings["xyz"]:
            ref = settings["xyz"]
            settings["xyz"] = resolve(ref, project, outputs)
            deferred = ref.startswith("@") and not check_files
        yield label, argparse.Namespace(project=project, **settings), deferred


def run_project(project, dest, state_path=None, calibrated=None):
    info = manifest(project)
    plan = list(steps(info))
    kinds = {s["step"] for s in plan}
    if state_path and "handoff" not in kinds or calibrated and "simulate" not in kinds:
        raise ValueError("Override has no corresponding configured step")
    outputs = {}
    for stage in plan:
        kind = stage["step"]
        target = dest / stage["id"]

        def source(key):
            return resolve(stage[key], project, outputs) if key in stage else None

        if kind == "build-xyz":
            build(project, target, stage.get("species"), stage.get("seed", 137))
        elif kind == "simulate":
            run_scenarios(Network(project, calibrated or source("calibrated")), target)
        elif kind == "handoff":
            handoff(project, state_path or source("state"), target)
        elif kind == "fit":
            fit(Network(project), source("data"), target)
        elif kind == "balance":
            account(project, source("data"), target)
        elif kind == "conformers":
            args = argparse.Namespace(
                project=project,
                species=stage["species"],
                cores=stage.get("cores", 1),
                quick=stage.get("quick", False),
                exe=stage.get("exe", "crest"),
                xyz=source("xyz"),
                frame=stage.get("frame", 1),
            )
            conformers(args, target)
        elif kind == "prepare":
            tasks = list(jobs(project, stage, outputs, True))
            for label, params, _ in tasks:
                prepare(params, target / label if len(tasks) > 1 else target)
        elif kind == "run-orca":
            target.mkdir(parents=True)
            for name in ("geom.xyz", "job.inp", "manifest.json"):
                shutil.copy2(source("job") / name, target / name)
            result = run_orca(target, stage.get("exe", "orca"))
            if result["status"] == "incomplete_or_review":
                raise ValueError(f"{target}: ORCA output needs review")
        elif kind in ("audit", "thermo"):
            result = (
                orca_audit(source("job"))
                if kind == "audit"
                else reaction_energy(source("recipe"))
            )
            target.mkdir(parents=True)
            write_json(target / "result.json", result)
            if result.get("status") == "incomplete_or_review":
                raise ValueError(f"{target}: ORCA output needs review")
        outputs[stage["id"]] = target


def validate(project):
    info = manifest(project)
    p = folder(project)
    plan = list(steps(info))
    kinds = {s["step"] for s in plan}
    if kinds & {"simulate", "fit"} or (p / "network.json").is_file():
        Network(project)
    if kinds & {"prepare", "conformers"} or "qm" in info:
        qm_settings(project)
        for name in read_json(p / "species.json"):
            rec = get_species(project, name)
            if "composition" in rec:
                comp = rec["composition"]
            elif (
                rec.get("smiles")
                and not (p / "xyz" / (rec.get("seed", name) + ".xyz")).is_file()
            ):
                comp = Counter(a.GetSymbol() for a in smiles_molecule(rec).GetAtoms())
            else:
                comp = Counter(geometry(project, name)[1])
            if not comp or any(
                e not in ELEMENTS or type(n) is not int or n < 1
                for e, n in comp.items()
            ):
                raise ValueError(f"{name}: invalid composition")
            total = sum(ELEMENTS[e] * n for e, n in comp.items())
            if any(
                (total - rec["charge"] - (m - 1)) % 2 for m in rec["multiplicities"]
            ):
                raise ValueError(f"{name}: spin parity mismatch")
    if "balance" in kinds or (p / "material.json").is_file():
        specification(project)
    if "paired" in info:
        vals = info["paired"]
        if (
            not 0 <= vals["practical"] <= 100
            or not 0 <= vals["gap"] <= 100
            or type(vals["min_blocks"]) is not int
            or vals["min_blocks"] < 2
        ):
            raise ValueError("Invalid paired endpoint thresholds")
    with tempfile.TemporaryDirectory() as tmp:
        outputs = {}
        if "handoff" in kinds or (p / "wet_plan.json").is_file():
            handoff(project, None, Path(tmp) / "initial")
        for stage in plan:
            for key in ("data", "state", "calibrated", "job", "recipe", "xyz"):
                if (
                    key in stage
                    and not stage[key].startswith("@")
                    and not resolve(stage[key], project, outputs).exists()
                ):
                    raise ValueError(f"{stage['id']}: missing input {stage[key]}")
            if stage["step"] == "build-xyz":
                records = read_json(p / "species.json")
                names = stage.get(
                    "species", [k for k, v in records.items() if v.get("smiles")]
                )
                if (
                    not isinstance(names, list)
                    or not names
                    or any(
                        n not in records or not records[n].get("smiles") for n in names
                    )
                ):
                    raise ValueError("build-xyz requires species with SMILES")
                for n in names:
                    smiles_molecule(get_species(project, n))
            if stage["step"] == "conformers":
                get_species(project, stage["species"])
            if stage["step"] == "prepare":
                for label, params, deferred in jobs(project, stage, outputs, False):
                    if not deferred:
                        prepare(params, Path(tmp) / stage["id"] / label)
            outputs[stage["id"]] = Path(tmp) / stage["id"]
    return info
