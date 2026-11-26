"""Prepare reviewed ion-cluster poses and expand solvent-specific CREST searches."""
import argparse
import json
import math
from pathlib import Path
import re
import shutil
import tempfile

from .pathway import _asset, _positive, _write
from .protocol import _digest, _identifier, _resources, _string, _unique_keys, load_plan
from .study import _xyz


def _vector(value, n):
    if (not isinstance(value, list) or len(value) != n or
            any(type(x) not in (float, int) or not math.isfinite(x) for x in value)):
        raise ValueError(f"expected finite {n}-component vector")
    return value


def _pose_atoms(component, pose, root):
    data = _asset(root, component["xyz"]).read_bytes()
    _xyz(data)
    lines = [line.split() for line in data.decode().splitlines()[2:] if line.strip()]
    qw, qx, qy, qz = _vector(pose["quaternion_wxyz"], 4)
    if not math.isclose(qw*qw + qx*qx + qy*qy + qz*qz, 1, abs_tol=1e-7):
        raise ValueError("rotation quaternion must be normalized")
    shift = _vector(pose["translation_angstrom"], 3)
    rotation = ((1-2*(qy*qy+qz*qz), 2*(qx*qy-qw*qz), 2*(qx*qz+qw*qy)),
                (2*(qx*qy+qw*qz), 1-2*(qx*qx+qz*qz), 2*(qy*qz-qw*qx)),
                (2*(qx*qz-qw*qy), 2*(qy*qz+qw*qx), 1-2*(qx*qx+qy*qy)))
    atoms = []
    for element, *values in lines:
        src = list(map(float, values))
        coord = [d + sum(a*b for a, b in zip(row, src)) for row, d in zip(rotation, shift)]
        atoms.append((element, coord))
    return atoms


def _assemble(pose, components, root):
    atoms, placed, charge = [], [], 0
    roles = []
    for placement in pose["placements"]:
        record = components[placement["component"]]
        _string(record["source"], "component source")
        fragment = _pose_atoms(record, placement, root)
        if type(record["charge"]) is not int or record["role"] not in {"cation", "anion", "solvent"}:
            raise ValueError("each fragment requires charge and a recognized role")
        required_sign = {"cation": 1, "anion": -1, "solvent": 0}[record["role"]]
        if (record["charge"] > 0) - (record["charge"] < 0) != required_sign:
            raise ValueError("charge conflicts with fragment role")
        if any(math.dist(a, b) < 0.5 for _, a in atoms for _, b in fragment):
            raise ValueError("interfragment atomic overlap in starting pose")
        placed.append(fragment)
        roles.append(record["role"])
        atoms.extend(fragment)
        charge += record["charge"]
    if roles.count("cation") != 1 or roles.count("anion") != 1:
        raise ValueError("pose needs one cation and one anion")
    family = pose["family"]
    if family not in {"contact", "solvent_separated"} or (family == "solvent_separated" and "solvent" not in roles):
        raise ValueError("solvent-separated starting geometry needs explicit solvent")
    def anchor(role):
        index = pose["anchor_atom_indices"][role]
        if type(index) is not int or index < 0 or index >= len(placed[roles.index(role)]):
            raise ValueError(f"invalid {role} anchor index")
        return placed[roles.index(role)][index][1]
    left, right = anchor("cation"), anchor("anion")
    distance = math.dist(left, right)
    lo, hi = _vector(pose["anchor_distance_range_angstrom"], 2)
    if not 0 < lo <= distance <= hi:
        raise ValueError("constructed ion anchors fall outside their reviewed interval")
    if family == "solvent_separated":
        # Require an explicit solvent atom in the intervening region. The
        # orientation/chemical assignment remains part of the frozen pose review.
        if not any(max(math.dist(left, x), math.dist(right, x)) < distance
                   for role, frag in zip(roles, placed) if role == "solvent" for _, x in frag):
            raise ValueError("explicit solvent does not intervene between the reviewed anchors")
    return atoms, charge, distance


def prepare(config_path, destination):
    path, dest = Path(config_path).resolve(), Path(destination).absolute()
    config = json.loads(path.read_text(), object_pairs_hook=_unique_keys)
    _string(config["approved_by"], "pose reviewer")
    _resources(config["resources"], "CREST resources")
    if dest.exists() or dest.is_symlink():
        raise ValueError("destination must be new")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=dest.parent) as tmp:
        out = Path(tmp) / "poses"
        out.mkdir()
        tasks, entries = [], []
        shutil.copyfile(path, out / "poses-source.json")
        for pose in config["poses"]:
            key = _identifier(pose["id"], "pose id")
            atoms, charge, distance = _assemble(pose, config["components"], path.parent)
            filename = key + ".xyz"
            seed = out / filename
            xyz = [str(len(atoms)), key] + [f"{symbol} " + " ".join(f"{x:.12f}" for x in xyz)
                                          for symbol, xyz in atoms]
            seed.write_text("\n".join(xyz) + "\n")
            electrons = _xyz(seed.read_bytes()) - charge
            mult = pose["multiplicity"]
            if type(mult) is not int or mult < 1 or electrons < mult - 1 or (electrons - mult + 1) % 2:
                raise ValueError("invalid ion-pair electronic state")
            settings = dict(config["crest_settings"], charge=charge, uhf=mult - 1)
            if (settings.get("method") not in {"gfn1", "gfn2", "gfnff"} or
                    settings.get("solvent_model") not in {"alpb", "gbsa"} or not
                    re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*", settings.get("solvent", ""))):
                raise ValueError("an explicit supported CREST model and solution environment are required")
            _positive(settings["energy_window_kcal_mol"], "energy window")
            tasks.append({"id": key, "engine": "crest", "stage": "conformers", "input": filename,
                          "engine_version": config["crest_version"], "settings": settings,
                          "resources": config["resources"], "requires": config.get("requires", []),
                          "outputs": ["crest.out", "crest_conformers.xyz", "crest.energies"]})
            entries.append({"id": key, "seed_sha256": _digest(seed), "starting_family": pose["family"],
                            "anchor_distance_angstrom": distance, "charge": charge, "multiplicity": mult,
                            "component_sha256": [_digest(_asset(path.parent, config["components"][p["component"]]["xyz"]))
                                                 for p in pose["placements"]]})
        _write(out / "seeds.json", entries)
        _write(out / "plan.json", {"name": config["id"], "attachments": ["poses-source.json", "seeds.json"],
                                    "tasks": tasks})
        load_plan(out / "plan.json")
        out.rename(dest)
    return {"plan": str(dest / "plan.json"), "poses": len(tasks)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config")
    parser.add_argument("destination")
    args = parser.parse_args()
    print(json.dumps(prepare(args.config, args.destination)))
