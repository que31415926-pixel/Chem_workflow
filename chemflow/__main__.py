"""Run with python -m chemflow from the repository root."""

import argparse
import importlib.metadata
import json
import os
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path

from .analysis import fit
from .material import account
from .network import Network, run_scenarios
from .projects import ROOT, discover, folder, manifest, read_json, steps, write_json
from .qm import conformers, orca_audit, prepare, run_orca
from .runner import run_project, validate
from .structures import build
from .thermo import reaction_energy
from .wet import handoff, paired


@contextmanager
def output_dir(path):
    path = Path(path)
    if path.exists():
        raise ValueError(f"Output already exists; select a new --out: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix=path.name + ".partial-", dir=path.parent))
    try:
        yield tmp
        if path.exists():
            raise ValueError("Output appeared during run; choose a new path")
        tmp.rename(path)
    except BaseException:
        print(f"Partial logs retained in {tmp}", file=sys.stderr)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--projects-dir", type=Path, help="Directory containing project folders"
    )
    commands = parser.add_subparsers(dest="cmd", required=True)
    commands.add_parser(
        "doctor", help="Show installed Python packages and external executables"
    )

    def command(name, help, project=True, dest=True):
        cmd = commands.add_parser(name, help=help)
        if project:
            cmd.add_argument("project", help="Project folder name; see projects")
        if dest:
            cmd.add_argument("--out", required=True, type=Path)
        return cmd

    command(
        "demo",
        "Run the declared workflow of each discovered project",
        False,
    )
    commands.add_parser("projects", help="List discovered project presets")
    init = commands.add_parser(
        "init", help="Create a self-contained preset from a template"
    )
    init.add_argument("name")
    init.add_argument(
        "--from",
        dest="source",
        help="Clone an existing project instead of the generic starter",
    )
    init.add_argument("--title")
    command("plan", "Show the ordered steps and input references", dest=False)
    command(
        "validate",
        "Validate configuration, cross-references and available input geometries",
        dest=False,
    )
    runner = command("run", "Run just this project's declared workflow")
    runner.add_argument("--state", type=Path)
    runner.add_argument("--calibrated", type=Path)
    command("list", "List registered geometry seeds", dest=False)
    b = command("build-xyz", "Build organic XYZs from configured SMILES (RDKit)")
    b.add_argument("--species", nargs="+")
    b.add_argument("--seed", type=int, default=137)
    c = command("conformers", "Run CREST on an organic species")
    c.add_argument("species")
    c.add_argument("--xyz", type=Path)
    c.add_argument("--frame", type=int, default=1)
    c.add_argument("--exe", default="crest")
    c.add_argument("--cores", type=int, default=1)
    c.add_argument("--quick", action="store_true")
    q = command(
        "prepare", "Generate a single ORCA input plus its geometry and run manifest"
    )
    q.add_argument("species")
    q.add_argument("--stage", choices=["min", "ts", "sp", "irc"], default="min")
    q.add_argument("--xyz", type=Path)
    q.add_argument("--frame", type=int, default=1)
    q.add_argument("--mult", type=int)
    q.add_argument("--cores", type=int, default=1)
    q.add_argument("--memory", type=int, default=1500, help="MB per core")
    q.add_argument("--temperature", type=float, help="Thermochemistry temperature in K")
    q.add_argument(
        "--reference-solvent",
        metavar="NAME",
        help="Reference calculations only",
    )
    q.add_argument("--method", default="screen")
    q.add_argument("--solvation", default="default")
    o = commands.add_parser(
        "run-orca", help="Run locally installed ORCA in a prepared job"
    )
    o.add_argument("job", type=Path)
    o.add_argument("--exe", default="orca")
    a = commands.add_parser("audit", help="Read orca.out produced by a prepared job")
    a.add_argument("job", type=Path)
    t = commands.add_parser(
        "thermo", help="Read a JSON recipe for a balanced solution reaction"
    )
    t.add_argument("recipe", type=Path)
    s = command(
        "simulate", "Export SciPy trajectories, discrimination scores and matching SBML"
    )
    s.add_argument("--calibrated", type=Path)
    f = command("fit", "Fit the declared rate parameters; reserve an entire condition")
    f.add_argument("--data", type=Path, required=True)
    m = command(
        "balance", "Compute molar atom/fragment balances from calibrated wet data"
    )
    m.add_argument("--data", type=Path, required=True)
    h = command("handoff", "Generate the next conditional wet-lab batch")
    h.add_argument("--state", type=Path)
    for name in ("paired", "paired-fe"):
        pp = commands.add_parser(
            name,
            help="Matched percent-yield decisions using project-specific thresholds",
        )
        if name == "paired":
            pp.add_argument("project")
        pp.add_argument("data", type=Path)
        pp.add_argument("--practical", type=float)
        pp.add_argument("--gap", type=float)
        pp.add_argument("--min-blocks", type=int)
    args = parser.parse_args()
    if args.projects_dir:
        os.environ["CHEMFLOW_PROJECTS_DIR"] = str(args.projects_dir.resolve())
    if args.cmd == "projects":
        for name in discover():
            print(name, manifest(name).get("title", name), sep="  ")
        return
    if args.cmd == "init":
        target = folder(args.name)
        if args.source:
            manifest(args.source)
            source = folder(args.source)
        else:
            source = ROOT / "templates" / "starter"
        with output_dir(target) as dest:
            shutil.copytree(source, dest, dirs_exist_ok=True)
            info = read_json(dest / "project.json")
            info["title"] = args.title or args.name
            write_json(dest / "project.json", info)
        print(target)
        return
    if args.cmd == "validate":
        validate(args.project)
        print(args.project, "OK")
        return
    if args.cmd == "doctor":
        for name in ("numpy", "scipy", "python-libsbml", "rdkit"):
            try:
                print(name, importlib.metadata.version(name))
            except importlib.metadata.PackageNotFoundError:
                print(name, "unavailable")
        for exe in ("crest", "xtb", "orca", "CopasiSE"):
            print(exe, shutil.which(exe) or "not on PATH")
        return
    if args.cmd == "list":
        print((folder(args.project) / "species.json").read_text(encoding="utf-8"))
        return
    if args.cmd == "plan":
        print(
            json.dumps(
                list(steps(manifest(args.project))), indent=2, ensure_ascii=False
            )
        )
        return
    if args.cmd in ("audit", "run-orca"):
        result = (
            run_orca(args.job, args.exe)
            if args.cmd == "run-orca"
            else orca_audit(args.job)
        )
        print(json.dumps(result, indent=2, ensure_ascii=False))
        if result["status"] == "incomplete_or_review":
            sys.exit(2)
        return
    if args.cmd == "thermo":
        print(json.dumps(reaction_energy(args.recipe), indent=2))
        return
    if args.cmd in ("paired", "paired-fe"):
        project = "fe" if args.cmd == "paired-fe" else args.project  # legacy CLI alias
        config = manifest(project).get("paired")
        if config is None:
            raise ValueError("This project has no paired percent-yield decision policy")
        opts = {
            k: (getattr(args, k) if getattr(args, k) is not None else config[k])
            for k in ("gap", "practical", "min_blocks")
        }
        print(
            json.dumps(
                paired(args.data, opts["gap"], opts["practical"], opts["min_blocks"]),
                indent=2,
            )
        )
        return
    with output_dir(args.out) as dst:
        if args.cmd == "demo":
            presets = discover()
            if not presets:
                raise ValueError("No projects found. Use init to create one.")
            for project in presets:
                run_project(project, dst / project)
        elif args.cmd == "run":
            run_project(args.project, dst, args.state, args.calibrated)
        elif args.cmd == "prepare":
            prepare(args, dst)
        elif args.cmd == "build-xyz":
            build(args.project, dst, args.species, args.seed)
        elif args.cmd == "conformers":
            conformers(args, dst)
        elif args.cmd == "simulate":
            run_scenarios(Network(args.project, args.calibrated), dst)
        elif args.cmd == "fit":
            fit(Network(args.project), args.data, dst)
        elif args.cmd == "balance":
            account(args.project, args.data, dst)
        elif args.cmd == "handoff":
            handoff(args.project, args.state, dst)
    print(args.out)


if __name__ == "__main__":
    try:
        main()
    except (
        ValueError,
        OSError,
        KeyError,
        RuntimeError,
        subprocess.CalledProcessError,
    ) as exc:
        sys.exit(f"chemflow: {exc}")
