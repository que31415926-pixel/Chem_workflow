"""Prepare and audit one ORCA job, or run an isolated CREST search."""

import hashlib
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path

from .projects import finite, folder, identifier, manifest, read_json

SYMBOLS = """H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og""".split()
ELEMENTS = {s: i for i, s in enumerate(SYMBOLS, 1)}
ORGANIC_ELEMENTS = set("H B C N O F Si P S Cl Ge As Se Br Sn Te I".split())
METHODS = {
    "literature": "wB97X-D3 RIJCOSX def2-TZVP def2/J",
    "screen": "r2SCAN-3c",
    "check": "PBE0 D4 def2-TZVP def2/J RIJCOSX",
}
JOBS = {"min": "TightOpt Freq", "ts": "OptTS Freq", "sp": "SP", "irc": "IRC"}


def solvent_block(preset):
    if not isinstance(preset, dict):
        raise ValueError("Solvent preset must be an object")
    modes = {"gas", "cpcm", "epsilon"}.intersection(preset)
    if len(modes) != 1:
        raise ValueError("Solvent preset requires exactly one of gas, cpcm, or epsilon")
    if "gas" in preset:
        if preset["gas"] is not True:
            raise ValueError("Use gas=true for a gas-phase calculation")
        return "", "", preset.get("label", "gas")
    if "cpcm" in preset:
        name = preset["cpcm"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", name):
            raise ValueError("Invalid named ORCA solvent")
        return f"CPCM({name})", "", preset.get("label", name)
    eps, n = preset["epsilon"], preset["refrac"]
    if not finite(eps) or not finite(n) or eps < 1 or n < 1:
        raise ValueError("Invalid dielectric constant or refractive index")
    return (
        "CPCM",
        f"%cpcm epsilon {float(eps)} refrac {n:.2f} end\n"
        if round(n, 2) == n
        else f"%cpcm epsilon {eps} refrac {n} end\n",
        preset.get("label", f"epsilon={eps}, n={n}"),
    )


def qm_settings(project):
    cfg = manifest(project)["qm"]
    if (
        not finite(cfg["temperature_K"])
        or cfg["temperature_K"] <= 0
        or "default" not in cfg["solvents"]
    ):
        raise ValueError(
            "QM requires positive temperature_K and a default solvent preset"
        )
    for p in cfg["solvents"].values():
        solvent_block(p)
    return cfg


def get_species(project, name):
    entries = read_json(folder(project) / "species.json")
    if name not in entries:
        raise ValueError("Unknown species. Options: " + ", ".join(entries))
    rec = entries[name]
    identifier(rec.get("seed", name))
    if (
        type(rec["charge"]) is not int
        or not rec["multiplicities"]
        or any(type(m) is not int or m < 1 for m in rec["multiplicities"])
    ):
        raise ValueError(f"{name}: invalid charge or multiplicities")
    return rec


def read_frames(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        count = int(lines[i])
        if count <= 0 or i + 2 + count > len(lines):
            raise ValueError("Invalid or truncated XYZ frame")
        frame = lines[i : i + 2 + count]
        elements = []
        for line in frame[2:]:
            fields = line.split()
            if (
                len(fields) != 4
                or fields[0] not in ELEMENTS
                or not all(math.isfinite(float(x)) for x in fields[1:])
            ):
                raise ValueError("Unsupported XYZ row: " + line)
            elements.append(fields[0])
        yield elements, "\n".join(frame) + "\n"
        i += count + 2


def geometry(project, name, xyz=None, frame=1):
    rec = get_species(project, name)
    if "composition" in rec:
        if not xyz:
            raise ValueError("Import literature coordinates first and supply --xyz")
        expected = [el for el, n in rec["composition"].items() for _ in range(n)]
    else:
        path = folder(project) / "xyz" / (rec.get("seed", name) + ".xyz")
        if path.is_file():
            expected, default = next(read_frames(path))
        elif rec.get("smiles") and xyz:
            from .structures import smiles_molecule

            expected = [a.GetSymbol() for a in smiles_molecule(rec).GetAtoms()]
        else:
            raise ValueError(
                f"{name}: missing seed XYZ; use build-xyz and pass its result with --xyz"
            )
    if not xyz:
        if frame != 1:
            raise ValueError("Multiple frames require --xyz ensemble.xyz")
        return rec, expected, default
    from collections import Counter

    if frame < 1:
        raise ValueError("Frame numbers start at 1")
    for n, (els, body) in enumerate(read_frames(xyz), 1):
        if n == frame:
            if Counter(els) != Counter(expected):
                raise ValueError("XYZ composition differs from the registered species")
            return rec, els, body
    raise ValueError("Selected frame missing")


def prepare(args, dest):
    rec, els, body = geometry(args.project, args.species, args.xyz, args.frame)
    charge = rec["charge"]
    mult = rec["multiplicities"][0] if args.mult is None else args.mult
    if (
        mult not in rec["multiplicities"]
        or (sum(ELEMENTS[e] for e in els) - charge - (mult - 1)) % 2
    ):
        raise ValueError("Unregistered spin or charge/multiplicity parity mismatch")
    if args.stage in ("ts", "irc") and not args.xyz:
        raise ValueError(
            "Supply the reaction-complex TS guess or optimized saddle with --xyz"
        )
    if args.cores < 1 or args.memory < 100:
        raise ValueError("Check cores and memory (MB per core)")
    cfg = qm_settings(args.project)
    temp = (
        cfg["temperature_K"]
        if getattr(args, "temperature", None) is None
        else args.temperature
    )
    if not finite(temp) or temp <= 0:
        raise ValueError("Temperature must be positive in K")
    name = getattr(args, "reference_solvent", None)
    if name:
        if args.solvation != "default":
            raise ValueError(
                "Choose either reference solvent or a project solvent preset"
            )
        preset = {"gas": True} if name == "gas" else {"cpcm": name}
        preset["label"] = f"Reference medium: {name}"
    else:
        if args.solvation not in cfg["solvents"]:
            raise ValueError("Available solvent presets: " + ", ".join(cfg["solvents"]))
        preset = cfg["solvents"][args.solvation]
    solvent, block, note = solvent_block(preset)
    methods = {**METHODS, **cfg.get("methods", {})}
    if args.method not in methods:
        raise ValueError("Available method names: " + ", ".join(methods))
    geom = "%geom Calc_Hess true end\n" if args.stage == "ts" else ""
    freq = f"%freq Temp {temp:.2f} end\n" if args.stage in ("min", "ts") else ""
    parallel = f"%pal nprocs {args.cores} end\n" if args.cores > 1 else ""
    header = f"! {methods[args.method]} TightSCF {JOBS[args.stage]} {solvent}\n"
    text = header + parallel + f"%maxcore {args.memory}\n" + block + geom + freq
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "geom.xyz").write_text(body)
    (dest / "job.inp").write_text(text + f"* xyzfile {charge} {mult} geom.xyz\n")
    meta = dict(
        project=args.project,
        species=args.species,
        note=rec.get("note", ""),
        stage=args.stage,
        method=methods[args.method],
        temperature_K=temp,
        solvent=note,
        charge=charge,
        multiplicity=mult,
        frame=args.frame,
        input_sha256=hashlib.sha256((dest / "job.inp").read_bytes()).hexdigest(),
        geometry_sha256=hashlib.sha256((dest / "geom.xyz").read_bytes()).hexdigest(),
    )
    (dest / "manifest.json").write_text(json.dumps(meta, indent=2) + "\n")


def conformers(args, dest):
    rec, els, body = geometry(
        args.project,
        args.species,
        getattr(args, "xyz", None),
        getattr(args, "frame", 1),
    )
    if not rec.get("crest_enabled", set(els) <= ORGANIC_ELEMENTS):
        raise ValueError(
            "Use reviewed coordination/spin alternatives for metals; this preset is organic only"
        )
    exe = shutil.which(args.exe)
    if not exe:
        raise ValueError("CREST executable missing: " + args.exe)
    if args.cores < 1:
        raise ValueError("--cores must be positive")
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "start.xyz").write_text(body)
    solvent = qm_settings(args.project).get("crest_solvent")
    command = [
        str(Path(exe).resolve()),
        "start.xyz",
        "--gfn2",
        "--chrg",
        str(rec["charge"]),
        "--uhf",
        str(rec["multiplicities"][0] - 1),
        *(["--alpb", solvent] if solvent else []),
        "-T",
        str(args.cores),
    ]
    if args.quick:
        command.append("--quick")
    (dest / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    with (dest / "crest.log").open("w") as out:
        subprocess.run(
            command,
            cwd=dest,
            stdout=out,
            stderr=subprocess.STDOUT,
            check=True,
            env={
                **os.environ,
                "OMP_NUM_THREADS": str(args.cores),
                "OPENBLAS_NUM_THREADS": "1",
            },
        )
    path = dest / "crest_conformers.xyz"
    if (
        "CREST terminated normally"
        not in (dest / "crest.log").read_text(encoding="utf-8", errors="replace")
        or not path.is_file()
        or not list(read_frames(path))
    ):
        raise ValueError("No parseable CREST ensemble: inspect crest.log")


def orca_audit(folder):
    meta = read_json(folder / "manifest.json")
    for key, file in [("geometry_sha256", "geom.xyz"), ("input_sha256", "job.inp")]:
        if hashlib.sha256((folder / file).read_bytes()).hexdigest() != meta[key]:
            raise ValueError(
                "Prepared inputs changed; regenerate a job with the reviewed coordinates/method"
            )
    data = (folder / "orca.out").read_text(errors="replace")
    # Each directory is for one ORCA invocation. Convergence of the SCF alone is insufficient.
    normal = "ORCA TERMINATED NORMALLY" in data
    opt = "THE OPTIMIZATION HAS CONVERGED" in data
    freq_section = re.findall(
        r"VIBRATIONAL FREQUENCIES(.*?)(?:NORMAL MODES|THERMOCHEMISTRY|$)", data, re.S
    )
    freqs = (
        re.findall(
            r"^\s*\d+:\s*([-+]?\d+(?:\.\d+)?)\s*cm\*\*-1", freq_section[-1], re.M
        )
        if freq_section
        else []
    )
    imaginary = [float(v) for v in freqs if float(v) < -0.1]

    def last(pattern):
        hits = re.findall(pattern, data, re.M)
        return float(hits[-1].replace("D", "E")) if hits else None

    en = last(r"FINAL SINGLE POINT ENERGY\s+([-+.\dEDed]+)")
    gibbs = last(r"Final Gibbs free energy\s+.*?([-+]\d+\.\d+)\s+Eh")
    enthalpy = last(r"Total enthalpy\s+.*?([-+]\d+\.\d+)\s+Eh")
    temp = last(r"^\s*Temperature\s+\.*\s*([\d.]+)\s*K")
    pressure = last(r"^\s*Pressure\s+\.*\s*([\d.]+)\s*atm")
    stage = meta["stage"]
    ok = normal and en is not None
    if stage in ("min", "ts"):
        ok = (
            ok
            and opt
            and bool(freqs)
            and gibbs is not None
            and temp is not None
            and abs(temp - meta["temperature_K"]) <= 0.05
            and pressure is not None
            and len(imaginary) == int(stage == "ts")
        )
    status = "numerical_checks_passed" if ok else "incomplete_or_review"
    if stage == "irc" and ok:
        status = "review_both_endpoints"
    return {
        **meta,
        "status": status,
        "normal_termination": normal,
        "opt_converged": opt,
        "reported_temperature_K": temp,
        "reported_pressure_atm": pressure,
        "electronic_Eh": en,
        "thermal_G_Eh": gibbs,
        "enthalpy_Eh": enthalpy,
        "imaginary_cm1": imaginary,
        "orca_output_sha256": hashlib.sha256(
            (folder / "orca.out").read_bytes()
        ).hexdigest(),
        "interpretation": "Review spin, structure and TS mode/IRC. Energies are not yields.",
    }


def run_orca(job, exe_name):
    exe = shutil.which(exe_name)
    if not exe:
        raise ValueError("Install licensed ORCA; supply --exe /full/path/to/orca")
    manifest = job / "manifest.json"
    if not manifest.is_file() or not (job / "job.inp").is_file():
        raise ValueError("Use a prepared directory")
    with (job / "orca.out").open("x") as out:
        subprocess.run(
            [str(Path(exe).resolve()), "job.inp"],
            cwd=job,
            stdout=out,
            stderr=subprocess.STDOUT,
            check=True,
        )
    return orca_audit(job)
