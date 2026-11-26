"""Turn a fully reviewed CREST ensemble into immutable study inputs."""
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile

from .protocol import _digest, _identifier, _relative, _static_file, _string, _unique_keys
from .study import _xyz
from .worker import NUMBER, Rejected, number


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_unique_keys)


def _write(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False) + "\n",
                    encoding="utf-8")


def _raw_frames(path):
    with path.open("rb") as stream:
        while True:
            header = stream.readline()
            if not header:
                return
            if not header.strip():
                continue
            if not header.strip().isdigit() or int(header) < 1:
                raise ValueError("invalid ensemble XYZ header")
            rows = [header, stream.readline()]
            for _ in range(int(header)):
                row = stream.readline()
                if not row:
                    raise ValueError("incomplete ensemble XYZ frame")
                rows.append(row)
            data = b"".join(rows)
            _xyz(data)  # Identical single-frame acceptance as compile-study.
            yield data, [line.split()[0].decode("ascii") for line in rows[2:]]


def _energies(path):
    values = []
    for row in path.read_text(encoding="utf-8").splitlines():
        cols = row.split()
        if not cols:
            continue
        if len(cols) != 2 or cols[0] != str(len(values) + 1) or not re.fullmatch(NUMBER, cols[1]):
            raise ValueError("expected sequential one-based CREST energy rows")
        try:
            values.append(number(cols[1]))
        except (Rejected, ValueError) as exc:
            raise ValueError("invalid finite CREST energy") from exc
    if not values or values[0] != 0 or any(a > b for a, b in zip(values, values[1:])):
        raise ValueError("expected energy-ranked CREST conformers relative to frame one")
    return values


def _inventory(ensemble, energies, output=None):
    for path in (ensemble, energies):
        if not path.is_file() or path.is_symlink():
            raise ValueError("expected regular CREST ensemble and energy files")
    values = _energies(energies)
    reference = None
    frames = []
    for index, (data, elements) in enumerate(_raw_frames(ensemble), 1):
        if reference is None:
            reference = elements
        if elements != reference or index > len(values):
            raise ValueError("CREST atom order, identity, or number of energies does not match")
        digest = hashlib.sha256(data).hexdigest()
        record = {"frame": index, "xyz_sha256": digest, "relative_energy": values[index - 1]}
        if output is not None:
            filename = f"audit/frames/frame-{index:06d}-{digest[:16]}.xyz"
            (output / filename).write_bytes(data)
            record["xyz"] = filename
        frames.append(record)
    if len(frames) != len(values):
        raise ValueError("CREST geometry and energy counts differ")
    return {"ensemble_sha256": _digest(ensemble), "energies_sha256": _digest(energies),
            "energy_unit": "kcal/mol", "atom_symbols": reference, "frames": frames}


def inventory(ensemble_path, energies_path):
    """Return the complete numbered and hashed candidate inventory for review."""
    return _inventory(Path(ensemble_path), Path(energies_path))


def _read_approvals(selection, counts):
    rows = selection.get("frames")
    if not isinstance(rows, list):
        raise ValueError("a review decision is required for every frame")
    by_frame, names = {}, set()
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("invalid conformer review row")
        i = row.get("frame")
        if type(i) is not int or i in by_frame or not 1 <= i <= counts:
            raise ValueError("unknown or duplicate reviewed frame")
        _string(row.get("reason"), "decision reason")
        status = row.get("decision")
        if status not in ("include", "exclude"):
            raise ValueError("each frame needs an include or exclude decision")
        if status == "include":
            name = _identifier(row.get("id"), "approved conformer id")
            if name.casefold() in names:
                raise ValueError("duplicate included conformer identifier")
            names.add(name.casefold())
            for field in ("degeneracy", "rotational_symmetry_number"):
                if type(row.get(field)) is not int or row[field] < 1:
                    raise ValueError("included frames require an explicit positive " + field)
        by_frame[i] = row
    if set(by_frame) != set(range(1, counts + 1)):
        raise ValueError("review leaves CREST candidates unaccounted for")
    if not names:
        raise ValueError("there are no approved candidates for a compute matrix")
    return by_frame


def _template_assets(cfg):
    files = {m["template"] for m in cfg["methods"]}
    for species in cfg["species"]:
        files.update(c["xyz"] for c in species.get("conformers", []))
    if cfg.get("evidence_plan"):
        files.add(cfg["evidence_plan"]["path"])
    return files


def _safe_copy(root, name, dest):
    source = _static_file(root, _relative(name, "study asset"), "study asset")
    target = dest / name
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if _digest(source) != _digest(target):
            raise ValueError("conflicting frozen asset " + name)
    else:
        shutil.copyfile(source, target)
    return target


def _prepare(ensemble, energies, selection_path, out):
    audit = out / "audit"
    (audit / "frames").mkdir(parents=True)
    for source, label in ((ensemble, "crest_conformers.xyz"), (energies, "crest.energies"),
                          (selection_path, "selection.json")):
        if source.is_symlink() or not source.is_file():
            raise ValueError(f"invalid input {source}")
        shutil.copyfile(source, audit / label)
    selection = _load(audit / "selection.json")
    if not isinstance(selection, dict):
        raise ValueError("expected a reviewed JSON object")
    reviewer = _string(selection.get("reviewer"), "reviewer")
    rationale = _string(selection.get("rationale"), "selection rationale")
    species_id = _identifier(selection.get("species_id"), "species_id")
    template = _static_file(selection_path.parent.resolve(), _relative(selection.get("study_template"), "study_template"),
                            "study_template")
    shutil.copyfile(template, audit / "study-template.json")
    if selection.get("study_template_sha256") != _digest(audit / "study-template.json"):
        raise ValueError("study template differs from its review commitment")
    cfg = _load(audit / "study-template.json")
    if "conformer_review" in cfg:
        raise ValueError("provide the base study before this ensemble was reviewed")
    matches = [s for s in cfg.get("species", []) if s["id"] == species_id]
    if len(matches) != 1 or matches[0].get("conformers"):
        raise ValueError("target species must occur once with an empty conformer slot")
    # Optional other approved species and method templates retain their names.
    asset_hashes = {}
    for name in sorted(_template_assets(cfg)):
        if name.casefold() == "study.json" or name.split("/")[0].casefold() == "audit":
            raise ValueError("study asset collides with review bundle " + name)
        asset_hashes[name] = _digest(_safe_copy(template.parent, name, out))

    catalog = _inventory(audit / "crest_conformers.xyz", audit / "crest.energies", out)
    for key in ("ensemble_sha256", "energies_sha256"):
        if selection.get(key) != catalog[key]:
            raise ValueError("review source commitment differs: " + key)
    approvals = _read_approvals(selection, len(catalog["frames"]))
    approved = []
    for frame in catalog["frames"]:
        decision = approvals[frame["frame"]]
        if decision.get("xyz_sha256") != frame["xyz_sha256"]:
            raise ValueError(f"reviewed frame differs from geometry at {frame['frame']}")
        frame.update(decision=decision["decision"], reason=decision["reason"])
        if frame["decision"] == "include":
            for name in ("id", "degeneracy", "rotational_symmetry_number"):
                frame[name] = decision[name]
            approved.append({"id": frame["id"], "xyz": frame["xyz"], "approved_by": reviewer,
                             "source": f"CREST SHA256 {catalog['ensemble_sha256']} frame {frame['frame']}",
                             "degeneracy": frame["degeneracy"],
                             "rotational_symmetry_number": frame["rotational_symmetry_number"]})
    matches[0]["conformers"] = approved
    # The ledger commits every original byte and every extracted frame, including rejected ones.
    frozen_files = {p.relative_to(out).as_posix(): _digest(p) for p in audit.rglob("*") if p.is_file()}
    ledger = {**catalog, "schema_version": 1, "species_id": species_id, "reviewer": reviewer,
              "rationale": rationale, "files": frozen_files, "template_assets": asset_hashes}
    ledger_file = audit / "ledger.json"
    _write(ledger_file, ledger)
    cfg["conformer_review"] = {"path": "audit/ledger.json", "sha256": _digest(ledger_file)}
    _write(out / "study.json", cfg)
    return {"total_frames": len(approvals), "included": len(approved),
            "excluded": len(approvals) - len(approved), "ledger_sha256": _digest(ledger_file)}


def prepare(ensemble_path, energies_path, selection_path, dest):
    """Publish all reviewed frames and an approved, directly compilable study configuration."""
    dest = Path(dest).absolute()
    if dest.is_symlink() or dest.exists():
        raise ValueError("destination already exists")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=dest.parent) as tmp:
        output = Path(tmp) / "review"
        output.mkdir()
        result = _prepare(Path(ensemble_path), Path(energies_path), Path(selection_path), output)
        os.rename(output, dest)
    return {**result, "study": str(dest / "study.json"), "ledger": str(dest / "audit/ledger.json")}


def archive_review_bundle(root, cfg, archive):
    """Recheck review hashes, decisions and matrix rows before exporting a submitted plan."""
    ref = cfg["conformer_review"]
    name = _relative(ref["path"], "conformer review")
    if name != "audit/ledger.json":
        raise ValueError("unexpected conformer review path")
    ledger_path = _safe_copy(root, name, archive)
    if _digest(ledger_path) != ref["sha256"]:
        raise ValueError("conformer review ledger changed")
    ledger = _load(ledger_path)
    if ledger.get("schema_version") != 1:
        raise ValueError("unsupported conformer review schema")
    for filename, digest in {**ledger["files"], **ledger["template_assets"]}.items():
        path = _safe_copy(root, filename, archive)
        if _digest(path) != digest:
            raise ValueError("changed reviewed asset " + filename)
    template = _load(archive / "audit/study-template.json")
    matches = [s for s in template["species"] if s["id"] == ledger["species_id"]]
    if len(matches) != 1 or matches[0].get("conformers"):
        raise ValueError("unexpected conformer review template")
    matches[0]["conformers"] = [
        {"id": f["id"], "xyz": f["xyz"], "approved_by": ledger["reviewer"],
         "source": f"CREST SHA256 {ledger['ensemble_sha256']} frame {f['frame']}",
         "degeneracy": f["degeneracy"], "rotational_symmetry_number": f["rotational_symmetry_number"]}
        for f in ledger["frames"] if f["decision"] == "include"
    ]
    template["conformer_review"] = ref
    if cfg != template:
        raise ValueError("study inputs differ from the reviewed candidate matrix or original template")
