"""ChemFlow 2: plan locally, execute on existing remote compute."""
import argparse
import importlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import uuid

from . import __version__
from . import controller as c
from .protocol import load_plan, explain_plan
from .transport import client_call, client_stream, TransportError, _unpack


COMMAND_MODULES = {
    "batch": ("batch", "plan and manage cloud executions"),
    "budget": ("budget", "estimate and reconcile a research budget"),
    "cloud-setup": ("cloud_setup", "configure private cloud resources"),
    "qualify": ("qualify", "verify a cloud calibration and register its capabilities"),
    "research": ("research", "inspect and compile a research pack"),
}


def show(data):
    print(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False))


def parser():
    cli = argparse.ArgumentParser(prog="chemflow2",description=__doc__)
    cli.add_argument("--version",action="version",version=__version__)
    sub = cli.add_subparsers(dest="cmd", required=True)
    for name in ("validate","plan","submit","slurm-export"):
        p = sub.add_parser(name)
        p.add_argument("plan",type=Path)
        if name=="plan":
            p.add_argument("--nodes", type=Path)
        if name=="submit":
            p.add_argument("--controller",required=True)
            p.add_argument("--id")
        if name=="slurm-export":
            p.add_argument("destination",type=Path)
    for name in ("status","cancel","report","sync"):
        p=sub.add_parser(name)
        p.add_argument("id")
        p.add_argument("--controller",required=True)
        if name=="sync":
            p.add_argument("destination",type=Path)
    analyses = ["science", "mechanism", "compile-study", "compile-pathway", "ion-pairs"]
    if importlib.util.find_spec("chemflow2.claims"):
        analyses.append("claims")
    for name in analyses:
        p=sub.add_parser(name)
        if name=="claims":
            p.add_argument("run_dir",type=Path)
        p.add_argument("input",type=Path)
        if name in {"compile-study", "compile-pathway", "ion-pairs"}:
            p.add_argument("destination",type=Path)
    p=sub.add_parser("conformer-inventory")
    p.add_argument("ensemble",type=Path)
    p.add_argument("energies",type=Path)
    p=sub.add_parser("prepare-conformers")
    for item in ("ensemble","energies","selection","destination"):
        p.add_argument(item,type=Path)
    p=sub.add_parser("estimate")
    for item in ("plan","history","nodes"):
        p.add_argument(item,type=Path)
    if importlib.util.find_spec("chemflow2.salts"):
        p=sub.add_parser("salt-inventory")
        p.add_argument("--prepare",type=Path)
        p.add_argument("--rdkit",action="store_true")
    for name, (_, help_text) in COMMAND_MODULES.items():
        sub.add_parser(name, help=help_text)
    if importlib.util.find_spec("chemflow2.paper_catalog"):
        sub.add_parser("paper-status")
    remote=sub.add_parser("remote").add_subparsers(dest="action",required=True)
    add=remote.add_parser("add",help="register one verified node on an existing controller")
    add.add_argument("profile",type=Path)
    add.add_argument("--controller",required=True)
    inspect=remote.add_parser("doctor")
    inspect.add_argument("host")
    p=sub.add_parser("doctor")
    p.add_argument("--host",help="existing SSH host alias; omitted shows this host")
    p=sub.add_parser("serve",help="start a controller on an existing remote host")
    p.add_argument("--nodes",type=Path)
    p.add_argument("--once",action="store_true",help="drain currently available jobs then exit")
    p=sub.add_parser("serve-stdio",help=argparse.SUPPRESS)
    p.add_argument("verb",choices=["submit","status","report","cancel","sync","doctor","register"])
    p.add_argument("id")
    sub.add_parser("worker",help="execute one staged task in the current directory")
    return cli


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    if args and args[0] in COMMAND_MODULES:
        module, _ = COMMAND_MODULES[args[0]]
        return importlib.import_module("." + module, __package__).main(args[1:])
    a=parser().parse_args(args)
    if a.cmd=="worker":
        from .worker import run_task
        receipt=run_task(c.read(Path("chemflow-task.json")),Path.cwd())
        c.atomic(Path("chemflow-result.json"),receipt)
        return 0 if receipt["ok"] else 1
    if a.cmd in {"validate","plan"}:
        plan=load_plan(a.plan)
        show({"name":plan["name"],"tasks":len(plan["tasks"]),"status":"valid"}
              if a.cmd=="validate" else explain_plan(plan, c.load_nodes(a.nodes) if a.nodes else []))
    elif a.cmd=="serve":
        root=c.root_dir()
        path=a.nodes or root/"nodes.json"
        c.serve(root,c.load_nodes(path),once=a.once,nodes_path=path)
    elif a.cmd=="remote":
        if a.action=="doctor":
            show(json.loads(client_call(a.host,["serve-stdio","doctor","node"])))
        else:
            node=c.read(a.profile)
            if isinstance(node,list) and len(node)==1:
                node=node[0]
            c.validate_nodes([node])
            show(json.loads(client_call(a.controller,["serve-stdio","register","node"],
                                        stdin=json.dumps(node).encode())))
    elif a.cmd=="doctor":
        show(json.loads(client_call(a.host, ["serve-stdio","doctor","node"])) if a.host else c.doctor())
    elif a.cmd=="submit":
        run_id=a.id or uuid.uuid4().hex
        c.run_path(Path("."),run_id)
        with tempfile.TemporaryFile() as wire,tempfile.TemporaryFile() as reply:
            c.bundle(a.plan.resolve(),wire)
            wire.seek(0)
            client_stream(a.controller,["serve-stdio","submit",run_id],stdin=wire,stdout=reply)
            reply.seek(0)
            show(json.load(reply))
    elif a.cmd=="sync":
        dst=a.destination.absolute()
        if dst.exists():
            raise ValueError("sync destination must be new")
        dst.parent.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(dir=dst.parent) as tmp,tempfile.TemporaryFile() as wire:
            src=Path(tmp)/"verified"
            src.mkdir()
            client_stream(a.controller,["serve-stdio","sync",a.id],stdout=wire)
            wire.seek(0)
            checks=_unpack(wire,src)
            if c.read(src/"state.json")["id"] != a.id:
                raise ValueError("unexpected run in response")
            src.rename(dst)
            show({"id":a.id,"files_verified":len(checks),"destination":str(dst)})
    elif a.cmd in {"status","report","cancel"}:
        show(json.loads(client_call(a.controller,["serve-stdio",a.cmd,a.id])))
    elif a.cmd=="serve-stdio":
        data=c.service(c.root_dir(),a.verb,a.id,sys.stdin.buffer,sys.stdout.buffer)
        if data is not None:
            show(data)
    elif a.cmd=="estimate":
        from .capacity import estimate
        show(estimate(a.plan,a.history,a.nodes))
    elif a.cmd=="salt-inventory":
        from .salts import inventory
        show(inventory(a.prepare,verify_rdkit=a.rdkit))
    elif a.cmd in {"conformer-inventory","prepare-conformers"}:
        from .conformers import inventory,prepare
        show(inventory(a.ensemble,a.energies) if a.cmd=="conformer-inventory"
             else prepare(a.ensemble,a.energies,a.selection,a.destination))
    elif a.cmd=="compile-pathway":
        from .pathway import compile_pathway
        show(compile_pathway(a.input,a.destination))
    elif a.cmd=="ion-pairs":
        from .ion_pairs import prepare
        show(prepare(a.input,a.destination))
    elif a.cmd=="slurm-export":
        from .slurm import export
        show(export(a.plan,a.destination))
    elif a.cmd in {"compile-study","paper-status"}:
        from .study import compile_study,paper_status
        show(paper_status() if a.cmd=="paper-status" else compile_study(a.input,a.destination))
    elif a.cmd=="claims":
        from .claims import assess
        show(assess(a.run_dir,a.input))
    else:
        if a.cmd=="science":
            from .science import evaluate
        else:
            from .mechanism import evaluate
        show(evaluate(c.read(a.input)))
    return 0


def entrypoint():
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        sys.exit(f"chemflow2: {exc}")


if __name__ == "__main__":
    entrypoint()
