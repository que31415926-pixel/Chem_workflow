"""Compile a reviewed TS guess to ORCA stationary-point, frequency and IRC tasks."""
import argparse
import json
import math
from pathlib import Path
import re
import shutil
import tempfile

from .protocol import (_digest, _identifier, _relative, _resources, _static_file, _string,
                       _strings, _unique_keys, load_plan)
from .study import _xyz


def _write(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf8")


def _asset(root, filename):
    return _static_file(root, _relative(filename, "source asset"), "source asset")


def _positive(value, field, zero=False):
    if type(value) not in (int, float) or not math.isfinite(value) or (value < 0 if zero else value <= 0):
        raise ValueError(f"invalid {field}")
    return value


def check_stationary_observations(task, obs, folder):
    from .worker import Rejected, task_path, xyz_atoms

    if task.get("endpoint_signatures") is not None:
        geometry = folder / Path(task["input"]).with_suffix(".xyz")
        xyz_atoms(geometry)
        coordinates = [list(map(float, row.split()[1:])) for row in geometry.read_text().splitlines()[2:] if row.strip()]
        matches = []
        for identity, intervals in task["endpoint_signatures"].items():
            if all(low <= math.dist(coordinates[i], coordinates[j]) <= high for i, j, low, high in intervals):
                matches.append(identity)
        obs["endpoint_candidates"] = matches
        obs["endpoint_geometry_sha256"] = _digest(geometry)
    point = task.get("stationary_point")
    if point is not None:
        if task["stage"] != "freq" or type(point["imaginary_count"]) is not int or point["imaginary_count"] not in (0, 1):
            raise Rejected("input", "invalid stationary-point frequency contract")
        cutoff = _positive(point["tolerance_cm1"], "imaginary tolerance", zero=True)
        count = xyz_atoms(task_path(folder, point["xyz"]))
        modes = obs.get("frequencies_cm-1", [])
        if len(modes) != 3 * count:
            raise Rejected("validation", "Hessian dimension differs from staged structure")
        obs["imaginary_modes_below_tolerance"] = sum(mode < -cutoff for mode in modes)
        if obs["imaginary_modes_below_tolerance"] != point["imaginary_count"]:
            raise Rejected("validation", "imaginary-mode count disagrees with target stationary point")
    state = task.get("state_diagnostics")
    if state is not None:
        if state.get("require_stability") is True and obs.get("wavefunction_stable") is not True:
            raise Rejected("validation", "stable SCF solution not confirmed by stability analysis")
        if "max_s2_deviation" in state:
            mult = state["multiplicity"]
            if type(mult) is not int or mult < 1:
                raise Rejected("input", "invalid spin multiplicity")
            limit = _positive(state["max_s2_deviation"], "spin diagnostic tolerance", zero=True)
            if "spin_squared" not in obs:
                raise Rejected("validation", "spin-squared expectation value missing")
            obs["spin_squared_deviation"] = obs["spin_squared"] - (mult ** 2 - 1) / 4
            if abs(obs["spin_squared_deviation"]) > limit:
                raise Rejected("validation", "spin-squared expectation exceeds reviewed tolerance")


def _emit_chain(config, method, output, guess, branch):
    resources = method["resources"]
    _resources(resources, "method resources")
    cores = resources["cores"]
    maxcore = resources["memory_mb"] * 3 // (4 * cores)
    if maxcore < 1:
        raise ValueError("insufficient per-rank memory")
    keywords = _string(method["keywords"], "method keywords")
    if not re.fullmatch(r"[A-Za-z0-9()=+*,./_ -]+", keywords):
        raise ValueError("method keywords must be one simple ORCA keyword line")
    if set(keywords.lower().split()) & {"opt", "optts", "irc", "freq", "numfreq", "moread", "anfreq", "scants", "neb-ts"}:
        raise ValueError("supply electronic method keywords separately from task stages")
    requires = _strings(method.get("requires", []), "requires")
    q, mult = config["charge"], config["multiplicity"]
    label = _identifier(method["id"], "method id") + "-" + branch
    tasks = []

    def emit(name, stage, directive, geometry, orbitals=None, extra=None, block="", products=()):
        task_id = label + "-" + name
        deck_name = task_id + ".inp"
        bindings = dict(extra or {})
        files = []
        if geometry:
            bindings["structure.xyz"] = geometry
        else:
            files.append("seed.xyz")
        if orbitals:
            bindings["start.gbw"] = orbitals
        deck = f"! {keywords} {directive}\n%pal nprocs {cores} end\n%maxcore {maxcore}\n"
        if name == "ts" and guess:
            files.append(guess)
            deck += f'! MORead\n%moinp "{guess}"\n'
        elif orbitals:
            deck += '! MORead\n%moinp "start.gbw"\n'
        xyz = "structure.xyz" if geometry else "seed.xyz"
        deck += f"{block}* xyzfile {q} {mult} {xyz}\n"
        (output / deck_name).write_text(deck, encoding="utf8")
        task = {"id": task_id, "engine": "orca", "stage": stage, "resources": resources,
                "requires": requires, "files": files, "input": deck_name, "from_deps": bindings,
                "deps": sorted({x.split("/")[0] for x in bindings.values()}),
                "observations_output": "values.json",
                "outputs": ["orca.out", "values.json", *(task_id + x for x in products)]}
        if stage == "freq":
            task["stationary_point"] = {"imaginary_count": int(name == "ts-freq"),
                "tolerance_cm1": config["imaginary_tolerance_cm1"], "xyz": xyz}
        stability_job = name.endswith("-stability")
        if stability_job or (stage in {"ts", "opt"} and "max_s2_deviation" in config):
            task["state_diagnostics"] = {"multiplicity": mult, "require_stability": stability_job}
            if "max_s2_deviation" in config:
                task["state_diagnostics"]["max_s2_deviation"] = config["max_s2_deviation"]
        if "engine_version" in method:
            task["engine_version"] = _string(method["engine_version"], "engine version")
        tasks.append(task)
        return task_id

    def artifact(ident, ext):
        return f"{ident}/{ident}{ext}"

    search = config.get("scan")
    if search:
        atoms, limits = search["bond_indices"], search["range_angstrom"]
        if (len(atoms) != 2 or len(set(atoms)) != 2 or
                any(type(i) is not int or i < 0 or i >= config["_n_atoms"] for i in atoms)
                or len(limits) != 2 or limits[0] == limits[1]
                or type(search["points"]) is not int or search["points"] < 3):
            raise ValueError("invalid one-dimensional ScanTS bond scan")
        for v in limits:
            _positive(v, "scan distance")
        block = (f"%geom Scan\n B {atoms[0]} {atoms[1]} = {limits[0]}, {limits[1]}, "
                 f"{search['points']}\n end\n end\n")
        directive = "ScanTS"
    else:
        block, directive = "%geom Calc_Hess true end\n", "OptTS"
    ts = emit("ts", "ts", directive, None, block=block, products=(".xyz", ".gbw"))

    def finish_stationary(parent, name):
        geom, orb = artifact(parent, ".xyz"), artifact(parent, ".gbw")
        if method.get("stability", False):
            diagnostic = emit(name + "-stability", "sp", "SP", geom, orb,
                              block="%scf StabPerform true StabRestartUHFifUnstable false end\n")
            prerequisite = {"stable.json": diagnostic + "/values.json"}
        else:
            prerequisite = {}
        freq_keyword = method.get("frequency_keyword", "Freq")
        if freq_keyword not in {"Freq", "NumFreq", "AnFreq"}:
            raise ValueError("unknown frequency method")
        return emit(name + "-freq", "freq", freq_keyword, geom, orb, extra=prerequisite, products=(".hess",))

    frequency = finish_stationary(ts, "ts")
    irc = emit("irc", "irc", "IRC", artifact(ts, ".xyz"), artifact(ts, ".gbw"),
               {"saddle.hess": artifact(frequency, ".hess")},
               '%irc Direction both InitHess read Hess_Filename "saddle.hess" end\n',
               ("_IRC_F.xyz", "_IRC_B.xyz", "_IRC_Full_trj.xyz"))
    endpoint_inputs = {}
    endpoint_producers = []
    for suffix in ("F", "B"):
        # ORCA can finish either direction last. Use the TS orbitals as a common
        # starting guess instead of silently assigning that last solution to both ends.
        endpoint = emit("end-" + suffix, "opt", "Opt", artifact(irc, "_IRC_" + suffix + ".xyz"),
                        artifact(ts, ".gbw"), products=(".xyz", ".gbw"))
        if "endpoint_signatures" in config:
            tasks[-1]["endpoint_signatures"] = config["endpoint_signatures"]
        minimum_freq = finish_stationary(endpoint, "end-" + suffix)
        endpoint_inputs[suffix + ".json"] = endpoint + "/values.json"
        endpoint_producers.extend([endpoint, minimum_freq])
    if "endpoint_signatures" in config:
        report_id = label + "-match-ends"
        input_name = report_id + ".json"
        _write(output / input_name, {"operation": "irc_endpoints", "expected": list(config["endpoint_signatures"]),
               "branches": {side: {"$artifact": side + ".json", "path": ["observations"]} for side in ("F", "B")}})
        tasks.append({"id": report_id, "engine": "analysis", "stage": "mechanism", "input": input_name,
                      "from_deps": endpoint_inputs, "deps": endpoint_producers, "outputs": ["endpoints.json"],
                      "resources": {"cores": 1, "memory_mb": 256, "scratch_mb": 1}})
    return tasks


def match_endpoints(data):
    expected = set(data["expected"])
    branches = data["branches"]
    if len(expected) != 2 or set(branches) != {"F", "B"}:
        raise ValueError("two expected structures and both IRC branches are required")
    observed = {}
    evidence = {}
    for side, record in branches.items():
        candidates = record["endpoint_candidates"]
        evidence[side] = record["endpoint_geometry_sha256"]
        if not re.fullmatch(r"[0-9a-f]{64}", evidence[side]):
            raise ValueError("missing endpoint geometry digest")
        observed[side] = candidates[0] if len(candidates) == 1 else None
    confirmed = set(observed.values()) == expected and len(set(evidence.values())) == 2
    return {"operation": "irc_endpoints", "branches": observed, "geometry_sha256": evidence,
            "status": "pass" if confirmed else "indeterminate"}


def compile_pathway(config_path, destination):
    path, dest = Path(config_path).resolve(), Path(destination).absolute()
    config = json.loads(path.read_text(), object_pairs_hook=_unique_keys)
    for field in ("id", "approved_by", "source", "coordination"):
        _string(config[field], field)
    _positive(config["imaginary_tolerance_cm1"], "imaginary_tolerance_cm1", zero=True)
    if "max_s2_deviation" in config:
        _positive(config["max_s2_deviation"], "max_s2_deviation", zero=True)
    seed = _asset(path.parent, config["xyz"])
    xyz = seed.read_bytes()
    charge, mult = config["charge"], config["multiplicity"]
    if type(charge) is not int or type(mult) is not int or mult < 1:
        raise ValueError("invalid registered electronic state")
    electrons = _xyz(xyz) - charge
    if electrons < mult - 1 or (electrons - mult + 1) % 2:
        raise ValueError("seed and electronic state disagree")
    config["_n_atoms"] = int(xyz.splitlines()[0])
    if "endpoint_signatures" in config:
        signatures = config["endpoint_signatures"]
        if not isinstance(signatures, dict) or len(signatures) != 2:
            raise ValueError("declare two endpoint identities")
        for ident, intervals in signatures.items():
            _identifier(ident, "endpoint identity")
            if not isinstance(intervals, list) or not intervals:
                raise ValueError("each endpoint needs at least one reviewed distance interval")
            for i, j, low, high in intervals:
                if (any(type(k) is not int or not 0 <= k < config["_n_atoms"] for k in (i, j))
                        or i == j or not 0 <= _positive(low, "minimum distance", zero=True) < _positive(high, "maximum distance")):
                    raise ValueError("invalid endpoint atom index or distance interval")
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        raise ValueError("destination must be new")
    with tempfile.TemporaryDirectory(dir=dest.parent) as tmp:
        out = Path(tmp) / "bundle"
        out.mkdir()
        shutil.copyfile(path, out / "pathway-source.json")
        (out / "seed.xyz").write_bytes(xyz)
        tasks = []
        for method in config["methods"]:
            if type(method.get("stability", False)) is not bool:
                raise ValueError("stability must be boolean")
            guesses = _strings(method.get("guess_files", []), "orbital guesses")
            branches = [(None, "default")]
            hashes = set()
            for filename in guesses:
                src = _asset(path.parent, filename)
                digest = _digest(src)
                if not src.stat().st_size or digest in hashes:
                    raise ValueError("empty or duplicated orbital guess")
                hashes.add(digest)
                name = "guess-" + digest + ".gbw"
                shutil.copyfile(src, out / name)
                branches.append((name, "guess-" + digest[:16]))
            for guess, branch in branches:
                tasks.extend(_emit_chain(config, method, out, guess, branch))
        _write(out / "plan.json", {"name": config["id"], "attachments": ["pathway-source.json"],
            "attachment_sha256": {"pathway-source.json": _digest(out / "pathway-source.json")},
            "tasks": tasks})
        load_plan(out / "plan.json")
        out.rename(dest)
    return {"tasks": len(tasks), "plan": str(dest / "plan.json")}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config")
    parser.add_argument("destination")
    args = parser.parse_args()
    print(json.dumps(compile_pathway(args.config, args.destination)))
