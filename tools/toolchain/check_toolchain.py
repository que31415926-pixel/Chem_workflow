#!/usr/bin/env python3
"""Check that this device's chemistry toolchain matches toolchain.lock.json.

Usage:
    python3 tools/toolchain/check_toolchain.py --orca /ABS/PATH/orca [--no-orca-jobs] [--json out.json]

Checks (standard library only, about 1 minute):
  ORCA    absolute path, version printed by a 1-process water HF job == lock, a 2-process (MPI) job
          terminates normally, binary SHA-256 == lock on Linux x86-64 (recorded elsewhere);
  MPI     mpirun on PATH and its Open MPI version == lock (WARN on macOS, where ORCA may bundle MPI);
  CREST / xTB / tblite / Python packages  versions == lock (WARN if missing: not every device needs them).
Writes toolchain_receipt.json (or --json) so that results can cite which device and versions made them.
Exit code 1 if any FAIL.
"""
import argparse
import hashlib
import importlib.metadata as md
import json
import os
import pathlib
import platform
import re
import shutil
import subprocess
import sys
import tempfile

REPO = pathlib.Path(__file__).resolve().parents[2]
LOCK = json.loads((REPO / "toolchain.lock.json").read_text())
WATER = "O 0 0 0.117\nH 0 0.757 -0.467\nH 0 -0.757 -0.467\n"
PYPKG = {"rdkit": "rdkit", "numpy": "numpy", "scipy": "scipy",
         "python-libsbml": "python-libsbml", "google-cloud-storage": "google-cloud-storage"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd, **kw):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=kw.pop("timeout", 60), **kw)
        return r.returncode, r.stdout + r.stderr
    except (OSError, subprocess.TimeoutExpired) as e:
        return None, str(e)


def orca_job(orca, nprocs):
    d = pathlib.Path(tempfile.mkdtemp(prefix="cf_tc_"))
    pal = "" if nprocs == 1 else f"%pal nprocs {nprocs} end\n"
    (d / "w.inp").write_text("! HF def2-SVP\n" + pal + "* xyz 0 1\n" + WATER + "*\n")
    try:
        _, out = run([orca, "w.inp"], cwd=d, timeout=600)
    finally:
        shutil.rmtree(d, ignore_errors=True)
    m = re.search(r"Program Version\s+(\S+)", out)
    return "ORCA TERMINATED NORMALLY" in out, (m.group(1) if m else None), out


def first_version(text):
    m = re.search(r"(\d+\.\d+(?:\.\d+)?)", text or "")
    return m.group(1) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orca", help="absolute path of the ORCA 6.1.0 binary")
    ap.add_argument("--no-orca-jobs", action="store_true", help="skip the two water test jobs")
    ap.add_argument("--json", default="toolchain_receipt.json")
    a = ap.parse_args()
    rows, fails = [], 0
    plat = f"{platform.system().lower()}-{platform.machine()}"

    def line(tag, item, text):
        nonlocal fails
        fails += tag == "FAIL"
        rows.append({"status": tag, "item": item, "detail": text})
        print(f"{tag:5s} {item:22s} {text}", flush=True)

    line("INFO", "host", f"{platform.node()} | {plat} | python {platform.python_version()}")
    want = LOCK["orca"]["version"]
    lin = LOCK["orca"]["platforms"]["linux-x86_64"]
    if not a.orca:
        line("WARN", "orca", "no --orca given; ORCA not checked (required on any device that runs DFT)")
    elif not (os.path.isabs(a.orca) and os.path.isfile(a.orca) and os.access(a.orca, os.X_OK)):
        line("FAIL", "orca", f"not an absolute path to an executable: {a.orca}")
    else:
        digest = sha256(a.orca)
        if plat in ("linux-x86_64", "linux-amd64"):
            line("PASS" if digest == lin["orca_binary_sha256"] else "FAIL", "orca binary sha256",
                 digest + ("" if digest == lin["orca_binary_sha256"] else f" != lock {lin['orca_binary_sha256']}"))
        else:
            line("INFO", "orca binary sha256", f"{digest} (platform build; recorded, not locked)")
        if not a.no_orca_jobs:
            ok1, ver, out1 = orca_job(a.orca, 1)
            line("PASS" if ver == want else "FAIL", "orca version", f"{ver} (lock {want})")
            line("PASS" if ok1 else "FAIL", "orca 1-process job",
                 "terminated normally" if ok1 else " | ".join(out1.strip().splitlines()[-6:]))
            ok2, _, out2 = orca_job(a.orca, 2)
            line("PASS" if ok2 else "FAIL", "orca 2-process (MPI)",
                 "terminated normally" if ok2 else
                 "failed; put Open MPI on PATH/LD_LIBRARY_PATH. " + " | ".join(out2.strip().splitlines()[-6:]))

    mpirun = shutil.which("mpirun")
    wmpi = LOCK["openmpi"]["version"]
    if not mpirun:
        line("WARN" if not a.orca else "FAIL", "mpirun", "not on PATH")
    else:
        _, out = run([mpirun, "--version"])
        v = first_version(out)
        tag = "PASS" if v == wmpi else ("WARN" if platform.system() == "Darwin" else "FAIL")
        line(tag, "open mpi", f"{v} at {mpirun} (lock {wmpi})")

    cf = LOCK["conda_forge"]
    for tool, flag in (("crest", "--version"), ("xtb", "--version")):
        exe = shutil.which(tool)
        if not exe:
            line("WARN", tool, f"not on PATH (lock {cf[tool]}; needed only on devices that run conformer searches)")
            continue
        _, out = run([exe, flag])
        m = re.search(r"[Vv]ersion\s+(\d+\.\d+(?:\.\d+)?)", out)
        v = m.group(1) if m else first_version(out)
        line("PASS" if v == cf[tool] else "FAIL", tool, f"{v} (lock {cf[tool]})")

    pyv = platform.python_version()
    line("PASS" if pyv == cf["python"] else "WARN", "python", f"{pyv} (lock {cf['python']}; minor 3.12 required)")
    if not pyv.startswith("3.12."):
        line("FAIL", "python minor", "ChemFlow is locked to Python 3.12")
    for key, dist in PYPKG.items():
        try:
            v = md.version(dist)
        except md.PackageNotFoundError:
            line("WARN", key, f"not installed (lock {cf[key]})")
            continue
        line("PASS" if v == cf[key] else "FAIL", key, f"{v} (lock {cf[key]})")

    verdict = "MATCHES toolchain.lock.json" if not fails else f"{fails} FAIL line(s): fix before producing results"
    line("RESULT", "toolchain", verdict)
    pathlib.Path(a.json).write_text(json.dumps({"lock_locked_on": LOCK["locked_on"], "platform": plat,
                                                "checks": rows}, indent=2) + "\n")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
