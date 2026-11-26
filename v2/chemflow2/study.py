"""Compile approved structure/method matrices and inspect the fixed v4 source ledger."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile

from .protocol import (_digest, _identifier, _mapping, _relative, _resources, _string, _strings,
                       _static_file, _unique_keys, load_plan)
from .worker import Rejected, validate_deck

_ELEMENTS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
             "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce "
             "Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "
             "Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl "
             "Mc Lv Ts Og").split()
_PLACEHOLDER = re.compile(r"\{\{([^{}]+)\}\}")
_COORDS = re.compile(r"^\s*\*\s*xyzfile\s+\{\{charge\}\}\s+\{\{multiplicity\}\}\s+\{\{xyz\}\}\s*$", re.I)


def _hash(data):
    return hashlib.sha256(data).hexdigest()


def _read_config(data):
    return _mapping(json.loads(data.decode("utf-8"), object_pairs_hook=_unique_keys),
                    "study configuration")


def paper_status(root=None):
    from .paper_catalog import paper_status as read_catalog
    return read_catalog(root)


def _indexed(items, label):
    if not isinstance(items, list) or not items:
        raise ValueError(f"{label} must be a nonempty list")
    result = {}
    for item in items:
        _mapping(item, label)
        key = _identifier(item.get("id"), label)
        if key in result:
            raise ValueError(f"duplicate {label}: {key}")
        result[key] = item
    return result


def _asset(root, filename):
    return _static_file(root, _relative(filename, "asset path"), "study asset").read_bytes()


def _xyz(data):
    lines = data.decode("utf-8").splitlines()
    try:
        count = int(lines[0])
        atoms = [line.split() for line in lines[2:] if line.strip()]
        if count < 1 or len(atoms) != count:
            raise ValueError()
        charge = 0
        for atom in atoms:
            if len(atom) != 4 or not all(math.isfinite(float(v)) for v in atom[1:]):
                raise ValueError()
            charge += 1 + _ELEMENTS.index(atom[0])
        return charge
    except (ValueError, IndexError) as exc:
        raise ValueError("expected a single XYZ geometry with element symbols and finite coordinates") from exc


def _methods(cfg, root):
    methods = _indexed(cfg.get("methods"), "method")
    for name, method in methods.items():
        stage = method.get("stage")
        geometry = method.get("geometry_from")
        if stage == "opt" and geometry is not None:
            raise ValueError(f"{name}: opt must consume the approved seed")
        if stage in {"freq", "sp"} and (geometry not in methods or methods[geometry].get("stage") != "opt"):
            raise ValueError(f"{name}: freq/sp must reference an explicit opt method")
        if stage not in {"opt", "freq", "sp"}:
            raise ValueError(f"unsupported matrix stage {stage!r}, supply a direct plan for other stages")
        _resources(method.get("resources"), name)
        _strings(method.get("requires", []), f"{name}.requires")
        template = _asset(root, method.get("template")).decode("utf-8")
        marks = _PLACEHOLDER.findall(template)
        directives = [line.split("#", 1)[0] for line in template.splitlines()]
        coords = [line for line in directives if re.match(r"^\s*\*\s*xyz", line, flags=re.I)]
        if (sorted(marks) != ["charge", "multiplicity", "xyz"] or len(coords) != 1
                or not _COORDS.fullmatch(coords[0])):
            raise ValueError(f"{name}: require one '* xyzfile {{{{charge}}}} {{{{multiplicity}}}} {{{{xyz}}}}' line")
        method["_template"] = template
    return sorted(methods.values(), key=lambda m: m["stage"] != "opt")


def _write_once(path, data):
    if path.exists():
        if path.read_bytes() != data:
            raise ValueError(f"content-address collision at {path.name}")
    else:
        path.write_bytes(data)


def _compile(config_path, output):
    config_bytes = config_path.read_bytes()
    cfg = _read_config(config_bytes)
    name = _string(cfg.get("name"), "study name")
    scope = cfg.get("scope")
    if scope not in {"engineering_calibration", "paper_inputs"}:
        raise ValueError("scope must be engineering_calibration or paper_inputs")
    if "claim_ids" in cfg:
        claim_ids = set(_strings(cfg["claim_ids"], "claim_ids"))
    elif scope == "paper_inputs":
        try:
            from .paper_catalog import _EXPECTED_CLAIMS as claim_ids
        except ImportError as exc:
            raise ValueError("paper_inputs requires an explicit claim_ids catalog") from exc
    else:
        claim_ids = set()
    methods = _methods(cfg, config_path.parent)
    species = _indexed(cfg.get("species"), "species")
    tasks = {}
    records = []
    archive = output / "study-source"
    archive.mkdir()
    _write_once(archive / config_path.name, config_bytes)
    for method in methods:
        template_dest = archive / method["template"]
        template_dest.parent.mkdir(parents=True, exist_ok=True)
        _write_once(template_dest, method["_template"].encode())
    check_dir = output.parent / "deck-check"
    check_dir.mkdir()
    for species_id, molecule in species.items():
        claims = _strings(molecule.get("claims", []), "claims")
        if (scope == "engineering_calibration" and claims) or (
                scope == "paper_inputs" and (not claims or not set(claims) <= claim_ids)):
            raise ValueError(f"invalid claim links for {scope}: {species_id}")
        states = _indexed(molecule.get("states"), f"{species_id} states")
        conformers = _indexed(molecule.get("conformers"), f"{species_id} conformers")
        for conf, info in conformers.items():
            reviewer = _string(info.get("approved_by"), "approved_by")
            source = _string(info.get("source"), "source")
            xyz = _asset(config_path.parent, info.get("xyz"))
            nuclei = _xyz(xyz)
            xyz_digest = _hash(xyz)
            xyz_dest = archive / info["xyz"]
            xyz_dest.parent.mkdir(parents=True, exist_ok=True)
            _write_once(xyz_dest, xyz)
            seed_name = f"seed-{xyz_digest}.xyz"
            _write_once(output / seed_name, xyz)
            for state, spec in states.items():
                q, mult = spec.get("charge"), spec.get("multiplicity")
                if type(q) is not int or type(mult) is not int or mult <= 0:
                    raise ValueError(f"illegal charge or multiplicity for {species_id}/{state}")
                electrons = nuclei - q
                if electrons < mult - 1 or (electrons - mult + 1) % 2:
                    raise ValueError(f"electron count and multiplicity disagree: {species_id}/{state}")
                producers = {}
                for method in methods:
                    stage = method["stage"]
                    parent = producers.get(method.get("geometry_from"))
                    geom = "optimized.xyz" if parent else seed_name
                    subs = {"charge": str(q), "multiplicity": str(mult), "xyz": geom}
                    deck = _PLACEHOLDER.sub(lambda m: subs[m[1]], method["_template"]).encode()
                    requires = sorted(method.get("requires", []))
                    resources = method["resources"]
                    fingerprint = {"deck": _hash(deck), "parent": parent, "stage": stage,
                                   "requires": requires, "resources": resources}
                    digest = _hash(json.dumps(fingerprint, sort_keys=True).encode())
                    task_id = stage + "-" + digest
                    deck_name = task_id + ".inp"
                    if task_id not in tasks:
                        # Only validate the frozen deck here. A child will use the real
                        # optimized geometry supplied by its parent at execution time.
                        (check_dir / geom).write_bytes(xyz)
                        (check_dir / deck_name).write_bytes(deck)
                        try:
                            validate_deck(check_dir / deck_name, stage, resources)
                        except Rejected as exc:
                            raise ValueError(f"{method['id']}: {exc}") from exc
                    _write_once(output / deck_name, deck)
                    task = {"id": task_id, "engine": "orca", "stage": stage, "study_scope": scope, "input": deck_name,
                            "files": [] if parent else [seed_name],
                            "outputs": ["orca.out"], "deps": [parent] if parent else [],
                            "resources": resources, "requires": requires}
                    if scope == "paper_inputs":
                        task["observations_output"] = "values.json"
                        task["outputs"].append("values.json")
                    if parent:
                        task["from_deps"] = {geom: f"{parent}/{parent}.xyz"}
                    if stage == "opt":
                        task["outputs"].append(task_id + ".xyz")
                    if stage == "freq":
                        task["outputs"].append(task_id + ".hess")
                    tasks[task_id] = task
                    producers[method["id"]] = task_id
                    records.append([species_id, state, q, mult, conf, source, reviewer, xyz_digest,
                                    method["id"], stage, _hash(method["_template"].encode()),
                                    task_id, parent or "", ",".join(claims)])
    with (output / "claims-map.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["species", "state", "charge", "multiplicity", "conformer", "source",
                         "approved_by", "seed_sha256", "method", "stage", "template_sha256",
                         "task_id", "geometry_from_task", "intended_claim_ids"])
        writer.writerows(records)
    evidence = cfg.get("evidence_plan")
    reference = {}
    if evidence is not None:
        _mapping(evidence, "evidence_plan")
        if scope != "paper_inputs":
            raise ValueError("scientific evidence plan requires paper_inputs scope")
        filename = _relative(evidence.get("path"), "evidence_plan")
        data = _asset(config_path.parent, filename)
        if _hash(data) != evidence.get("sha256"):
            raise ValueError("evidence_plan digest differs from the approved source")
        evidence_dest = archive / filename
        evidence_dest.parent.mkdir(parents=True, exist_ok=True)
        _write_once(evidence_dest, data)
        reference = {"evidence_plan": {"path": evidence_dest.relative_to(output).as_posix(),
                                       "sha256": _hash(data)}}
    if "conformer_review" in cfg:
        from .conformers import archive_review_bundle
        archive_review_bundle(config_path.parent, _read_config(config_bytes), archive)
    attachments = {"claims-map.csv", *(p.relative_to(output).as_posix()
                                       for p in archive.rglob("*") if p.is_file())}
    hashes = {name: _digest(output / name) for name in sorted(attachments)}
    plan = {"name": name, "study_scope": scope, "attachments": sorted(attachments),
            "attachment_sha256": hashes, **reference, "study_config_sha256": _hash(config_bytes),
            "claims_map_sha256": hashes["claims-map.csv"],
            "tasks": list(tasks.values())}
    (output / "plan.json").write_text(json.dumps(plan, ensure_ascii=False, indent=2), encoding="utf-8")
    load_plan(output / "plan.json")
    return {"tasks": len(tasks), "matrix_rows": len(records), "scope": scope}


def compile_study(config_path: Path, dest: Path) -> dict:
    """Validate all cells, then publish a complete immutable input bundle in a new directory."""
    dest = Path(dest).absolute()
    if dest.exists() or dest.is_symlink():
        raise ValueError(f"destination already exists: {dest}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".chemflow-study-", dir=dest.parent) as temp:
        output = Path(temp) / "bundle"
        output.mkdir()
        info = _compile(Path(config_path).resolve(), output)
        os.rename(output, dest)
    return {**info, "plan": str(dest / "plan.json"), "claims_map": str(dest / "claims-map.csv")}
