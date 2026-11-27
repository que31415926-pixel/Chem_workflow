"""Frozen CREST 3 conformer-search adapter with checked geometry and energy files."""
import math
from pathlib import Path
import re

from .worker import (NUMBER, Rejected, configured_binary, launch, number, sha256, task_path,
                     xyz_atoms)

ORDER_TOLERANCE_KCAL = 0.0015  # crest.energies resolution (0.001 kcal/mol) plus float noise


def _xyz_frames(path):
    with path.open(encoding="utf-8") as stream:
        while True:
            header = stream.readline()
            if not header:
                return
            if not header.strip():
                continue
            if not header.strip().isdigit() or not 1 <= int(header) <= 1_000_000 or not stream.readline():
                raise Rejected("validation", "CREST produced an invalid XYZ header")
            elements = []
            for _ in range(int(header)):
                parts = stream.readline().split()
                if len(parts) != 4 or not re.fullmatch(r"[A-Za-z]{1,3}", parts[0]):
                    raise Rejected("validation", "Invalid CREST atom record")
                for value in parts[1:]:
                    if not re.fullmatch(NUMBER, value):
                        raise Rejected("validation", "Invalid CREST coordinate")
                    number(value)
                elements.append(parts[0].capitalize())
            yield elements


def run_crest(task, folder, provenance):
    xyz = task_path(folder, task["input"])
    if xyz.suffix.lower() != ".xyz":
        raise Rejected("input", "CREST expects a single-frame XYZ")
    xyz_atoms(xyz)
    required = {"crest_conformers.xyz", "crest.energies"}
    if not required.issubset(task["outputs"]):
        raise Rejected("input", "Declare both CREST conformer XYZ and relative-energy files")
    settings = task.get("settings", {})
    fields = {"method", "charge", "uhf", "solvent_model", "energy_window_kcal_mol"}
    if not isinstance(settings, dict) or not fields.issubset(settings) or settings.keys() - fields - {"solvent"}:
        raise Rejected("input", "Specify only the supported frozen CREST settings")
    method, chrg, uhf = (settings[key] for key in ("method", "charge", "uhf"))
    if method not in ("gfn1", "gfn2", "gfnff") or type(chrg) is not int or type(uhf) is not int or uhf < 0:
        raise Rejected("input", "Invalid CREST method, charge, or unpaired-electron count")
    window = settings["energy_window_kcal_mol"]
    if type(window) not in (int, float) or not math.isfinite(window) or window <= 0:
        raise Rejected("input", "Specify a positive finite CREST energy window")
    solvent, model = settings.get("solvent"), settings["solvent_model"]
    if not ((model == "gas" and solvent is None) or
            (model in ("alpb", "gbsa") and isinstance(solvent, str) and re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*", solvent))):
        raise Rejected("input", "Specify gas or a named solvent and a supported model")
    required_version = task.get("engine_version")
    if not isinstance(required_version, str) or not re.fullmatch(r"\d+(?:\.\d+){1,3}(?:[-+][A-Za-z0-9.-]+)?", required_version):
        raise Rejected("input", "Pin a CREST release in engine_version")
    executable = configured_binary("crest")
    threads = task["resources"]["cores"]
    args = [str(executable), xyz.relative_to(folder).as_posix(), "--" + method, "--chrg", str(chrg),
            "--uhf", str(uhf), "--ewin", str(window), "-T", str(threads)]
    if model != "gas":
        args.extend(["--" + model, solvent])
    log = folder / "crest.out"
    provenance.update(executable=str(executable), executable_sha256=sha256(executable),
                      argv=args, reservations=task["resources"])
    code = launch(args, folder, log, threads, task.get("timeout_s", 7 * 86400))
    provenance["exit_code"] = code
    if code:
        raise Rejected("engine", f"CREST exited with status {code}")
    normal, version = False, None
    with log.open(errors="replace") as stream:
        for line in stream:
            if line.strip().lower() == "crest terminated normally.":
                normal = True
            match = re.match(r"\s*(?:CREST\s+)?Version\s*:?\s*v?(\d+(?:\.\d+){1,3}(?:[-+][A-Za-z0-9.-]+)?)\b", line, re.I)
            if match and version is None:
                version = match[1]
    provenance["engine_version"] = version
    if not normal:
        raise Rejected("engine", "CREST normal termination marker absent")
    if version != required_version:
        raise Rejected("capability", f"CREST version {version!r} differs from {required_version}")
    energies = []
    for line in task_path(folder, "crest.energies").read_text().splitlines():
        row = line.split()
        if not row:
            continue
        if (len(row) != 2 or row[0] != str(len(energies) + 1)
                or not re.fullmatch(NUMBER, row[1])):
            raise Rejected("validation", "Invalid index or energy in crest.energies")
        energies.append(number(row[1]))
    frames = 0
    expected = next(_xyz_frames(xyz))
    for elements in _xyz_frames(task_path(folder, "crest_conformers.xyz")):
        if elements != expected:
            raise Rejected("validation", "CREST atom count, order, or element changed between frames")
        frames += 1
    if not frames or frames != len(energies):
        raise Rejected("validation", "CREST ensemble and energy counts differ")
    if energies[0] != 0 or min(energies) < -ORDER_TOLERANCE_KCAL:
        raise Rejected("validation", "Expected CREST energies relative to the lowest frame 1")
    # crest.energies is printed to 0.001 kcal/mol and CREGEN groups near-degenerate rotamers, so adjacent
    # values can be out of order by rounding. Keep the ensemble; record (do not hide) any disorder.
    disorder = [i + 2 for i, (a, b) in enumerate(zip(energies, energies[1:])) if a > b]
    worst = max((a - b for a, b in zip(energies, energies[1:])), default=0)
    result = {"normal_termination": True, "engine_version": version, "conformers": frames,
              "n_atoms": len(expected), "relative_energies_kcal_mol": energies,
              "scientific_gate": "not_evaluated"}
    if disorder:
        result["energy_order"] = {"out_of_order_frames": disorder, "max_inversion_kcal_mol": worst,
                                  "within_rounding": worst <= ORDER_TOLERANCE_KCAL}
    return result
