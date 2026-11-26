"""Portable submission of new research DAGs to an existing Slurm allocation."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import uuid
from datetime import datetime, timezone
from graphlib import TopologicalSorter
from pathlib import Path

from .protocol import _relative, load_plan, matches, static_sources


def _hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write(path, value):
    temp = path.with_name("." + path.name + "." + uuid.uuid4().hex)
    try:
        temp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def _read(path):
    return json.loads(path.read_text())


def _declared_copy(source, dest, digest):
    if source.is_symlink() or not source.is_file():
        raise ValueError(f"Missing or non-regular input {source}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, dest)
    if _hash(dest) != digest:
        raise ValueError(f"SHA256 mismatch for {source}")


def export(plan_path: Path, dest: Path) -> dict:
    """Package the validated inputs and the exact installed executor for transfer."""
    plan_path = Path(plan_path).resolve()
    dest = Path(dest).absolute()
    plan = load_plan(plan_path)
    files = {f for task in plan["tasks"] for f in static_sources(task).values()}
    for name in plan.get("attachments", []):
        files.add(_relative(name, "plan attachment"))
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.mkdir(mode=0o700)
    try:
        inputs = dest / "inputs"
        inputs.mkdir()
        _write(dest / "protocol.json", plan)
        for name in sorted(files):
            src = plan_path.parent / name
            if not src.resolve().is_relative_to(plan_path.parent):
                raise ValueError(f"Input outside protocol directory: {name}")
            _declared_copy(src, inputs / name, _hash(src))
        package = dest / "chemflow2"
        package.mkdir()
        for source in Path(__file__).parent.glob("*.py"):
            shutil.copyfile(source, package / source.name)
        if any(t["engine"] in {"upstream", "chemflow"} for t in plan["tasks"]):
            from .upstream import VENDOR
            if not (VENDOR / "chemflow" / "__main__.py").is_file():
                raise ValueError("Upstream task requires the bundled ChemFlow v1 sources")
            shutil.copytree(VENDOR, dest / "vendor/chemflow-v1",
                            ignore=shutil.ignore_patterns("__pycache__", "*.pyc", ".git"))
        (dest / "submit.py").write_text(
            "from pathlib import Path\nfrom chemflow2.slurm import submit_cli\n"
            "if __name__ == '__main__':\n    submit_cli(Path(__file__).resolve().parent)\n"
        )
        (dest / "collect.py").write_text(
            "import json, sys\nfrom pathlib import Path\nfrom chemflow2.slurm import collect\n"
            "if __name__ == '__main__':\n"
            "    report = collect(Path(__file__).resolve().parent,sys.argv[1])\n"
            "    print(json.dumps(report,indent=2))\n"
            "    sys.exit(0 if report['all_tasks_succeeded'] else 1)\n"
        )
        # sbatch copies the shell file to a spool location. Use the absolute shared
        # root passed as argv. Compute nodes need no network or SSH credentials.
        (dest / "stage.sh").write_text(
            '#!/bin/sh\nset -eu\ncd "$1"\nexec python3 -m chemflow2.slurm "$2" "$3"\n'
        )
        checks = {str(p.relative_to(dest)): _hash(p) for p in dest.rglob("*") if p.is_file()}
        _write(dest / "checksums.json", checks)
        return {"bundle": str(dest), "tasks": len(plan["tasks"]), "inputs": len(files),
                "protocol_sha256": checks["protocol.json"]}
    except BaseException:
        shutil.rmtree(dest)
        raise


def _verify(root, names, checks):
    for name in names:
        path = root / name
        if name not in checks or not path.is_file() or path.is_symlink() or _hash(path) != checks[name]:
            raise ValueError(f"Changed or missing bundle file: {name}")


def _validate_site(site, tasks):
    # 'host' is unused by Slurm, but shared resource admission uses the same
    # node schema as SSH. A site file is supplied by the operator.
    limits = dict(site)
    limits.setdefault("host", "slurm")
    limits.setdefault("name", "slurm")
    for task in tasks:
        if not matches(task, limits):
            raise ValueError(f"Site cannot admit task {task['id']}")
    if any(t["engine"] == "orca" for t in tasks):
        if not site.get("architecture") or not re.fullmatch("[0-9a-f]{64}", str(site.get("orca_sha256", ""))):
            raise ValueError("ORCA site profile requires architecture and orca_sha256")
    minutes = site.get("wall_minutes")
    if isinstance(minutes, bool) or not isinstance(minutes, int) or minutes < 1:
        raise ValueError("Site requires a positive integer wall_minutes")


def submit(root, site, partition=None, check=False):
    """Run on an allowed submit host, with bundle on a shared filesystem."""
    root = root.resolve()
    checks = _read(root / "checksums.json")
    _verify(root, checks, checks)
    tasks = _read(root / "protocol.json")["tasks"]
    _validate_site(site, tasks)
    if partition and not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", partition):
        raise ValueError("Invalid Slurm partition")
    ordered = list(TopologicalSorter({t["id"]: t.get("deps", []) for t in tasks}).static_order())
    if check:
        return {"status": "admissible", "ordered_tasks": ordered}
    run_id = uuid.uuid4().hex
    run = root / "runs" / run_id
    run.mkdir(parents=True, mode=0o700)
    (run / "tasks").mkdir()
    _write(run / "site.json", site)
    _write(run / "frozen.json", {"checksums_sha256": _hash(root / "checksums.json"),
                                  "site_sha256": _hash(run / "site.json")})
    jobs = {}
    state = {"run": run_id, "jobs": jobs, "state": "submitting"}
    tasks_by_id = {t["id"]: t for t in tasks}
    try:
        for task_id in ordered:
            task = tasks_by_id[task_id]
            resources = task["resources"]
            mpi = task["engine"] == "orca"
            ranks, threads = (resources["cores"], 1) if mpi else (1, resources["cores"])
            args = ["sbatch", "--parsable", "--nodes=1", f"--cpus-per-task={threads}",
                    f"--ntasks={ranks}", f"--mem={resources['memory_mb']}M",
                    "--time=" + str(site["wall_minutes"]), "--kill-on-invalid-dep=yes",
                    "--job-name=cf-" + task_id, "--output=" + str(run / "%j.log")]
            if site.get("reserve_tmp", False):
                args.append(f"--tmp={resources['scratch_mb']}M")
            if partition:
                args.append("--partition=" + partition)
            if task.get("deps"):
                args.append("--dependency=afterok:" + ":".join(jobs[d] for d in task["deps"]))
            reply = subprocess.run(args + [str(root / "stage.sh"), str(root), run_id, task_id],
                                   text=True, capture_output=True, check=True)
            job = reply.stdout.strip().split(";", 1)[0]
            if not re.fullmatch("[0-9]+", job):
                raise ValueError(f"Unexpected sbatch response for {task_id}: {reply.stdout!r}")
            jobs[task_id] = job
            _write(run / "submission.json", state)
        state["state"] = "submitted"
    except (OSError, ValueError, subprocess.SubprocessError, KeyboardInterrupt) as exc:
        state["state"] = "submission_incomplete"
        state["error"] = str(exc)
        if isinstance(exc, subprocess.CalledProcessError):
            state["stderr"] = exc.stderr
        raise RuntimeError(f"Incomplete submission. Inspect {run}/submission.json and squeue") from exc
    finally:
        _write(run / "submission.json", state)
    return state


def _check_allocation(task, site):
    resources = task["resources"]
    if site.get("architecture") and platform.machine() != site["architecture"]:
        raise ValueError("Allocated node architecture differs from site profile")
    if int(os.environ["SLURM_NTASKS"]) * int(os.environ.get("SLURM_CPUS_PER_TASK", 1)) < resources["cores"]:
        raise ValueError("Allocated CPU count below request")
    if int(os.environ["SLURM_MEM_PER_NODE"]) < resources["memory_mb"]:
        raise ValueError("Allocated memory below request")
    if task["engine"] == "orca":
        binary = shutil.which(os.environ.get("CHEMFLOW_ORCA", "orca"))
        if not binary or _hash(Path(binary)) != site["orca_sha256"]:
            raise ValueError("ORCA executable differs from approved compute node binary")


def run_stage(root: Path, run_id: str, task_id: str):
    """Run inside a single allocation and publish artifacts on shared storage."""
    if not re.fullmatch("[0-9a-f]{32}", run_id):
        raise ValueError("Invalid run ID")
    root = root.resolve()
    run = root / "runs" / run_id
    frozen = _read(run / "frozen.json")
    if (_hash(root / "checksums.json") != frozen["checksums_sha256"]
            or _hash(run / "site.json") != frozen["site_sha256"]):
        raise ValueError("Frozen admission data changed")
    checks = _read(root / "checksums.json")
    _verify(root, ["protocol.json"] + [k for k in checks if k.startswith(("chemflow2/", "vendor/"))], checks)
    task = next((t for t in _read(root / "protocol.json")["tasks"] if t["id"] == task_id), None)
    if task is None:
        raise ValueError("Unknown task")
    final = run / "tasks" / task_id
    if final.exists():
        raise ValueError(f"Task already published: {task_id}")
    # Slurm requeue is a new start; do not merge half-written artifacts.
    staging = Path(tempfile.mkdtemp(prefix=".publish-", dir=final.parent))
    result = {"ok": False, "status": "failed", "failure_kind": "slurm_stage",
              "error": None, "artifacts": {}}
    try:
        site = _read(run / "site.json")
        _validate_site(site, [task])
        _check_allocation(task, site)
        local = Path(os.environ.get("SLURM_TMPDIR") or os.environ.get("TMPDIR") or tempfile.gettempdir())
        if shutil.disk_usage(local).free // 2**20 < task["resources"]["scratch_mb"]:
            raise ValueError("Insufficient free scratch on compute node")
        with tempfile.TemporaryDirectory(prefix="chemflow-", dir=local) as folder:
            work = Path(folder)
            for name, source in static_sources(task).items():
                _declared_copy(root / "inputs" / source, work / name, checks["inputs/" + source])
            for dest, ref in task.get("from_deps", {}).items():
                parent, name = ref.split("/", 1)
                parent_dir = run / "tasks" / parent
                parent_result = _read(parent_dir / "chemflow-result.json")
                if parent_result.get("ok") is not True:
                    raise ValueError(f"Dependency {parent} did not succeed")
                _declared_copy(parent_dir / name, work / dest, parent_result["artifacts"][name])
            from .worker import run_task
            result = run_task(task, work)
            for name, digest in result["artifacts"].items():
                relative = Path(name)
                if relative.is_absolute() or ".." in relative.parts or name == "chemflow-result.json":
                    raise ValueError("Invalid artifact name")
                _declared_copy(work / name, staging / name, digest)
            if result["ok"] and not set(task["outputs"]).issubset(result["artifacts"]):
                raise ValueError("Incomplete worker publication")
    except Exception as exc:
        result = {"ok": False, "status": "failed", "failure_kind": "slurm_stage",
                  "error": str(exc), "artifacts": {}}
    result["slurm"] = {"job_id": os.environ.get("SLURM_JOB_ID"),
                       "node": platform.node(), "finished_utc": datetime.now(timezone.utc).isoformat(),
                       "protocol_sha256": checks["protocol.json"]}
    _write(staging / "chemflow-result.json", result)
    staging.rename(final)
    return result


def collect(root, run_id):
    """Verify a returned run and summarize software completion per task."""
    if not re.fullmatch("[0-9a-f]{32}", run_id):
        raise ValueError("Invalid run ID")
    run = root / "runs" / run_id
    frozen = _read(run / "frozen.json")
    if (_hash(root / "checksums.json") != frozen["checksums_sha256"]
            or _hash(run / "site.json") != frozen["site_sha256"]):
        raise ValueError("Frozen admission data changed")
    checks = _read(root / "checksums.json")
    _verify(root, checks, checks)
    statuses = {}
    for task in _read(root / "protocol.json")["tasks"]:
        task_dir = run / "tasks" / task["id"]
        receipt = task_dir / "chemflow-result.json"
        if not receipt.exists():
            statuses[task["id"]] = "no_receipt"
            continue
        try:
            result = _read(receipt)
            artifacts = result["artifacts"]
            for name, digest in artifacts.items():
                relative = Path(name)
                if relative.is_absolute() or ".." in relative.parts:
                    raise ValueError("Invalid artifact path")
                source = task_dir / name
                if source.is_symlink() or not source.is_file() or _hash(source) != digest:
                    raise ValueError("Artifact missing or changed")
            if result["ok"] is True:
                if not set(task["outputs"]).issubset(artifacts):
                    raise ValueError("Missing declared outputs")
                statuses[task["id"]] = "succeeded"
            else:
                statuses[task["id"]] = "failed"
        except (OSError, KeyError, ValueError, AttributeError, TypeError):
            statuses[task["id"]] = "invalid_receipt"
    return {"run": run_id, "tasks": statuses,
            "all_tasks_succeeded": all(v == "succeeded" for v in statuses.values())}


def submit_cli(root):
    parser = argparse.ArgumentParser()
    parser.add_argument("site", type=Path)
    parser.add_argument("--partition")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    print(json.dumps(submit(root, _read(args.site), args.partition, args.check), indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Expected run ID and task ID")
    sys.exit(0 if run_stage(Path.cwd(), sys.argv[1], sys.argv[2])["ok"] else 1)
