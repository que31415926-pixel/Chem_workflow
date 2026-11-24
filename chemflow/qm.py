"""Prepare and audit one ORCA job, or run an isolated CREST search."""

import hashlib
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path

from .network import ROOT, read_json

ELEMENTS = {
    "H": 1,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Co": 27,
    "Fe": 26,
    "Br": 35,
}
METHODS = {
    "literature": "wB97X-D3 RIJCOSX def2-TZVP def2/J",
    "screen": "r2SCAN-3c",
    "check": "PBE0 D4 def2-TZVP def2/J RIJCOSX",
}
JOBS = {"min": "TightOpt Freq", "ts": "OptTS Freq", "sp": "SP", "irc": "IRC"}


def get_species(project, name):
    entries = read_json(ROOT / "projects" / project / "species.json")
    if name not in entries:
        raise ValueError("Unknown species. Options: " + ", ".join(entries))
    return entries[name]


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
        path = ROOT / "projects" / project / "xyz" / (rec.get("seed", name) + ".xyz")
        expected, default = next(read_frames(path))
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
    mult = args.mult or rec["multiplicities"][0]
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
    temp = (
        args.temperature
        if getattr(args, "temperature", None) is not None
        else (323.15 if args.project == "co" else 298.15)
    )
    if not math.isfinite(temp) or temp <= 0:
        raise ValueError("Temperature must be positive in K")
    if getattr(args, "reference_solvent", None):
        if args.solvation != "default":
            raise ValueError("Choose a reference solvent or dielectric sensitivity")
        name = args.reference_solvent
        solvent = "" if name == "gas" else f"CPCM({name})"
        block = ""
        note = f"Reference medium: {name}"
    elif args.project == "co":
        if args.solvation != "default":
            raise ValueError("Fe-only dielectric sensitivity option")
        solvent = "CPCM(THF)"
        block = ""
        note = "THF"
    else:
        eps = {"low": 22.0, "default": 30.0, "high": 38.0}[args.solvation]
        solvent = "CPCM"
        block = f"%cpcm epsilon {eps:.1f} refrac 1.40 end\n"
        note = f"DMA/iPrOH continuum sensitivity eps={eps}; assumed, n=1.40"
    geom = "%geom Calc_Hess true end\n" if args.stage == "ts" else ""
    freq = f"%freq Temp {temp:.2f} end\n" if args.stage in ("min", "ts") else ""
    parallel = f"%pal nprocs {args.cores} end\n" if args.cores > 1 else ""
    header = f"! {METHODS[args.method]} TightSCF {JOBS[args.stage]} {solvent}\n"
    text = header + parallel + f"%maxcore {args.memory}\n" + block + geom + freq
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "geom.xyz").write_text(body)
    (dest / "job.inp").write_text(text + f"* xyzfile {charge} {mult} geom.xyz\n")
    meta = dict(
        project=args.project,
        species=args.species,
        note=rec["note"],
        stage=args.stage,
        method=METHODS[args.method],
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
    rec, els, body = geometry(args.project, args.species)
    if "Fe" in els or "Co" in els:
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
    solvent = "thf" if args.project == "co" else "dmso"
    command = [
        str(Path(exe).resolve()),
        "start.xyz",
        "--gfn2",
        "--chrg",
        str(rec["charge"]),
        "--uhf",
        str(rec["multiplicities"][0] - 1),
        "--alpb",
        solvent,
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
