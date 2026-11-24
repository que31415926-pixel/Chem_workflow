"""Run with python -m chemflow from the repository root."""

import argparse
import importlib.metadata
import json
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path

from .analysis import fit
from .material import account
from .network import PROJECTS, ROOT, Network, run_scenarios
from .qm import conformers, orca_audit, prepare
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
    commands = parser.add_subparsers(dest="cmd", required=True)
    commands.add_parser(
        "doctor", help="Show installed Python packages and external executables"
    )

    def command(name, help, project=True, dest=True):
        cmd = commands.add_parser(name, help=help)
        if project:
            cmd.add_argument("project", choices=PROJECTS)
        if dest:
            cmd.add_argument("--out", required=True, type=Path)
        return cmd

    command(
        "demo",
        "Build both baseline handoffs, scenario tables, SBML and metal input demos",
        False,
    )
    command("list", "List registered geometry seeds", dest=False)
    c = command("conformers", "Run CREST on an organic species")
    c.add_argument("species")
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
        choices=["gas", "benzene", "THF"],
        help="Reference calculations only",
    )
    q.add_argument(
        "--method", choices=["screen", "check", "literature"], default="screen"
    )
    q.add_argument("--solvation", choices=["default", "low", "high"], default="default")
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
    f = command("fit", "Fit one scalar per scenario; reserve an entire condition")
    f.add_argument("--data", type=Path, required=True)
    m = command(
        "balance", "Compute molar atom/fragment balances from calibrated wet data"
    )
    m.add_argument("--data", type=Path, required=True)
    h = command("handoff", "Generate the next conditional wet-lab batch")
    h.add_argument("--state", type=Path)
    pp = commands.add_parser(
        "paired-fe", help="Classify matched measured yields with interval/LOQ support"
    )
    pp.add_argument("data", type=Path)
    pp.add_argument("--practical", type=float, default=30.0)
    pp.add_argument("--gap", type=float, default=15.0)
    args = parser.parse_args()
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
        print(
            (ROOT / "projects" / args.project / "species.json").read_text(
                encoding="utf-8"
            )
        )
        return
    if args.cmd == "run-orca":
        exe = shutil.which(args.exe)
        if not exe:
            raise ValueError("Install licensed ORCA; supply --exe /full/path/to/orca")
        manifest = args.job / "manifest.json"
        if not manifest.is_file() or not (args.job / "job.inp").is_file():
            raise ValueError("Use a prepared directory")
        with (args.job / "orca.out").open("x") as out:
            subprocess.run(
                [str(Path(exe).resolve()), "job.inp"],
                cwd=args.job,
                stdout=out,
                stderr=subprocess.STDOUT,
                check=True,
            )
    if args.cmd in ("audit", "run-orca"):
        result = orca_audit(args.job)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        if result["status"] == "incomplete_or_review":
            sys.exit(2)
        return
    if args.cmd == "thermo":
        print(json.dumps(reaction_energy(args.recipe), indent=2))
        return
    if args.cmd == "paired-fe":
        print(json.dumps(paired(args.data, args.gap, args.practical), indent=2))
        return
    with output_dir(args.out) as dst:
        if args.cmd == "demo":
            for project, metal in [("co", "co_me"), ("fe", "fe_oep_cl")]:
                net = Network(project)
                run_scenarios(net, dst / project / "scenarios")
                handoff(project, None, dst / project / "wet")
                # Explicit demonstration spin. Research runs compare all plausible multiplicities.
                demo = argparse.Namespace(
                    project=project,
                    species=metal,
                    xyz=None,
                    frame=1,
                    mult=1 if project == "co" else 6,
                    cores=1,
                    memory=1500,
                    stage="min",
                    method="screen",
                    solvation="default",
                )
                prepare(demo, dst / project / "orca_input")
        elif args.cmd == "prepare":
            prepare(args, dst)
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
