"""Single-writer remote coordinator. Only worker hosts execute research tasks."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from datetime import datetime, timezone
import heapq
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time

from .protocol import _order, load_plan, matches, static_sources, validate_node
from .transport import execute, _pack, _unpack

TERMINAL = {"finished", "failed", "blocked", "cancelled"}
RESOURCES = ("cores", "memory_mb", "scratch_mb")


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def atomic(path, data):
    raw = json.dumps(data, ensure_ascii=False, allow_nan=False, indent=2) + "\n"
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False, encoding="utf-8") as f:
        try:
            f.write(raw)
            f.flush()
            os.fsync(f.fileno())
        except BaseException:
            Path(f.name).unlink(missing_ok=True)
            raise
    os.replace(f.name, path)


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def root_dir():
    return Path(os.environ.get("CHEMFLOW_HOME", "~/.chemflow2")).expanduser().resolve()


def run_path(root, name):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,79}", name):
        raise ValueError("invalid run identifier")
    return root / "runs" / name


def validate_nodes(items):
    if not isinstance(items, list):
        raise ValueError("nodes.json must contain a JSON array")
    by_name = {}
    hosts = set()
    for n in items:
        validate_node(n)
        place = n.get("placement","ssh")
        if place not in {"ssh","controller"}:
            raise ValueError("placement must be ssh or controller")
        if place == "controller" and any(x.get("placement")=="controller" for x in by_name.values()):
            raise ValueError("declare the controller compute capacity only once")
        if n["name"] in by_name or n["host"] in hosts:
            raise ValueError("each SSH host alias and node name must occur once")
        by_name[n["name"]] = n
        hosts.add(n["host"])
    return list(by_name.values())


def load_nodes(path):
    return validate_nodes(read(path))


def register_node(root, node):
    """Add one explicitly qualified host. Concurrent registrations never overwrite another node."""
    import fcntl
    root.mkdir(parents=True, mode=0o700, exist_ok=True)
    validate_node(node)
    path = root/"nodes.json"
    with (root/"nodes.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        nodes = load_nodes(path) if path.exists() else []
        if node in nodes:
            return {"name":node["name"], "status":"already_registered"}
        nodes.append(node)
        atomic(path,validate_nodes(nodes))
    return {"name":node["name"], "status":"registered"}


def declared_files(plan):
    static = {f for t in plan["tasks"] for f in static_sources(t).values()}
    static.update(plan.get("attachments", []))
    return static


def bundle(path, out):
    """Ship only the referenced inputs; config files elsewhere stay on the laptop."""
    plan = load_plan(path)
    static = declared_files(plan)
    if "protocol.json" in static:
        raise ValueError("protocol.json is reserved for the manifest")
    with tempfile.TemporaryDirectory() as tmp:
        stage = Path(tmp)
        atomic(stage / "protocol.json", plan)
        for name in static:
            dst = stage / name
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path.parent / name, dst)
        _pack(stage, {"protocol.json", *static}, out)


def accept(root, run_id, stream):
    dst = run_path(root, run_id)
    dst.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    with tempfile.TemporaryDirectory(dir=dst.parent, prefix=".incoming-") as tmp:
        staged = Path(tmp) / "run"
        src = staged / "source"
        src.mkdir(parents=True)
        checksums = _unpack(stream, src)
        plan = load_plan(src / "protocol.json")
        if set(checksums) != {"protocol.json", *declared_files(plan)}:
            raise ValueError("bundle contents differ from the declared inputs")
        (staged / "jobs").mkdir()
        atomic(staged / "source-sha256.json", checksums)
        atomic(staged / "state.json", {"id": run_id, "study": plan["name"], "status": "queued",
               "submitted_utc": now(), "tasks": {t["id"]: {"status": "pending"} for t in plan["tasks"]}})
        if dst.exists():
            raise ValueError("run identifier already exists; choose a new one")
        staged.rename(dst)
    return {"id": run_id, "status": "queued"}


def prepare(run, task):
    path = run / "jobs" / task["id"]
    path.mkdir()
    src = run / "source"
    inputs = read(run / "source-sha256.json")
    refs = [(name, src / source, inputs[source]) for name, source in static_sources(task).items()]
    for name, ref in task["from_deps"].items():
        parent, artifact = ref.split("/", 1)
        parent_dir = run / "jobs" / parent
        result = read(parent_dir / "chemflow-result.json")
        if result.get("ok") is not True:
            raise ValueError(f"dependency {parent} did not succeed")
        refs.append((name, parent_dir / artifact, result["artifacts"][artifact]))
    for name, source, expected in refs:
        dest = path / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, dest)
        if digest(dest) != expected:
            raise ValueError(f"checksum mismatch for input {name}")
    return path


def attempt(run, task, node, executor):
    try:
        path = prepare(run, task)
        result = executor(node, task, path)
        if result.get("ok") is True:
            for name in task["outputs"]:
                if result["artifacts"].get(name) != digest(path / name):
                    raise ValueError(f"missing or changed worker artifact {name}")
    except Exception as exc:
        result = {"ok": False, "status":"failed", "failure_kind":"dispatch_or_integrity",
                  "error": str(exc), "artifacts": {}, "observations": {}, "provenance": {}}
    result.update(node=node["name"], node_capacity_class=node.get("capacity_class"), controller_finished_utc=now())
    path = run / "jobs" / task["id"]
    path.mkdir(exist_ok=True)
    atomic(path / "chemflow-result.json", result)
    return result


@contextmanager
def controller_lock(root):
    import fcntl
    root.mkdir(parents=True, exist_ok=True)
    with (root / "controller.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError("another controller is already active") from exc
        yield


def schedule_order(tasks):
    """Prefer ready ancestors of explicit priorities, then the longest remaining chain."""
    topo = _order(tasks)
    task_by_id = {t["id"]: t for t in tasks}
    children = {name: [] for name in topo}
    remaining = {}
    inherited = {}
    depth = {}
    for task in tasks:
        key = task["id"]
        remaining[key] = len(task["deps"])
        for parent in task["deps"]:
            children[parent].append(key)
    for key in reversed(topo):
        inherited[key] = max([task_by_id[key].get("priority", 0), *[inherited[k] for k in children[key]]])
        depth[key] = 1 + max((depth[k] for k in children[key]), default=0)
    index = {key:i for i,key in enumerate(topo)}
    def score(key):
        return (-inherited[key], -depth[key], index[key], key)
    ready = [score(key) for key in topo if not remaining[key]]
    heapq.heapify(ready)
    order = []
    while ready:
        *_, key = heapq.heappop(ready)
        order.append(key)
        for child in children[key]:
            remaining[child] -= 1
            if not remaining[child]:
                heapq.heappush(ready, score(child))
    return order


def registered_worker(node, task, path):
    if node.get("placement", "ssh") == "ssh":
        return execute(node, task, path)
    if node.get("placement") != "controller" or platform.system() != "Linux":
        raise ValueError("controller workers require explicit placement on a Linux controller")
    # This process runs on the remote Linux service host. The client never calls it.
    atomic(path/"chemflow-task.json", task)
    env = os.environ.copy()
    checkout = str(Path(__file__).resolve().parents[1])
    env["PYTHONPATH"] = checkout + os.pathsep + env.get("PYTHONPATH", "")
    with (path/"chemflow-worker.out.log").open("xb") as out, (path/"chemflow-worker.err.log").open("xb") as err:
        done = subprocess.run([sys.executable, "-m", "chemflow2", "worker"], cwd=path, env=env,
                              stdin=subprocess.DEVNULL, stdout=out, stderr=err, check=False)
    receipt = read(path/"chemflow-result.json")
    if done.returncode != (0 if receipt.get("ok") else 1):
        raise ValueError("worker exit status disagrees with the receipt")
    # Keep registered evidence and frozen inputs, discard out-of-scope scratch.
    keep = {Path(n) for n in receipt["artifacts"]} | {Path(n) for n in task["files"]} | {
        Path(n) for n in task["from_deps"]} | {Path(n) for n in
            ("chemflow-task.json","chemflow-result.json","chemflow-worker.out.log","chemflow-worker.err.log")}
    if task.get("input"):
        keep.add(Path(task["input"]))
    for member in path.rglob("*"):
        if (member.is_file() or member.is_symlink()) and member.relative_to(path) not in keep:
            member.unlink()
    return receipt


def serve(root, nodes, *, once=False, interval=1.0, executor=registered_worker, nodes_path=None):
    """Schedule across all new runs with one shared capacity ledger."""
    limit = 64
    revision = None
    active, cache = {}, {}
    runs = root / "runs"
    runs.mkdir(parents=True, exist_ok=True)
    with controller_lock(root), ThreadPoolExecutor(max_workers=limit) as pool:
        while True:
            if nodes_path is not None:
                info = nodes_path.stat()
                current = (info.st_ino, info.st_mtime_ns)
                if current != revision:
                    nodes = load_nodes(nodes_path)
                    revision = current
            for run in sorted(runs.iterdir()):
                if run.name.startswith(".") or not (run / "state.json").is_file():
                    continue
                try:
                    state = read(run / "state.json")
                except (OSError, ValueError) as exc:
                    # Preserve the corrupted journal for inspection; other runs can proceed.
                    corrupted = run / "state.unreadable.json"
                    if (run/"state.json").is_file() and not corrupted.exists():
                        shutil.copyfile(run/"state.json", corrupted)
                    atomic(run/"state.json", {"id":run.name, "status":"failed", "tasks":{},
                                                "error":f"invalid state file: {exc}"})
                    continue
                if state["status"] in TERMINAL:
                    continue
                old = json.dumps(state, sort_keys=True)
                if run not in cache:
                    source = run/"source"/"protocol.json"
                    try:
                        expected = read(run/"source-sha256.json")
                        if digest(source) != expected["protocol.json"]:
                            raise ValueError("frozen protocol checksum mismatch")
                        plan = load_plan(source)
                        for name in plan.get("attachments", []):
                            if digest(source.parent/name) != expected[name]:
                                raise ValueError(f"study attachment checksum mismatch: {name}")
                    except (OSError, ValueError, KeyError) as exc:
                        state.update(status="failed", error=f"frozen input validation failed: {exc}")
                        atomic(run/"state.json", state)
                        continue
                    cache[run] = plan, {x["id"]: x for x in plan["tasks"]}, schedule_order(plan["tasks"])
                plan, tasks, order = cache[run]
                rows = state["tasks"]
                cancel = (run / "cancel.request").exists()
                for key, row in rows.items():
                    marker = run, key
                    if row["status"] == "running":
                        if marker not in active:
                            row.update(status="failed", error="controller interrupted; remote attempt needs inspection")
                        elif active[marker][0].done():
                            receipt = active.pop(marker)[0].result()
                            row.update(status="succeeded" if receipt["ok"] else "failed",
                                       failure_kind=receipt.get("failure_kind"), error=receipt.get("error"))
                    if row["status"] == "pending" and cancel:
                        row.update(status="cancelled")
                # Traverse the DAG so a blocked parent propagates in the same pass.
                for key in order:
                    row, task = rows[key], tasks[key]
                    if row["status"] != "pending":
                        continue
                    deps = [rows[d]["status"] for d in task["deps"]]
                    if any(s in {"failed", "blocked", "cancelled"} for s in deps):
                        row.update(status="blocked", error="incomplete dependency")
                        continue
                    if any(s != "succeeded" for s in deps):
                        continue
                    eligible = [n for n in nodes if matches(task, n)]
                    if not eligible:
                        row.update(status="blocked", error="no eligible node; see plan")
                        continue
                    for node in eligible:
                        reserved = [r for _, name, r in active.values() if name == node["name"]]
                        if len(reserved) >= node.get("concurrency", 1) or len(active) >= limit:
                            continue
                        if any(task["resources"][k]+sum(r[k] for r in reserved)>node[k] for k in RESOURCES):
                            continue
                        row.update(status="running", node=node["name"], submitted_utc=now())
                        # Persist dispatch before starting; never replay an ambiguous attempt.
                        state["status"] = "running"
                        atomic(run/"state.json", state)
                        active[run,key] = (pool.submit(attempt, run, task, node, executor),
                                           node["name"],task["resources"])
                        break
                statuses = {r["status"] for r in rows.values()}
                if statuses & {"running", "pending"}:
                    state["status"] = "draining" if cancel else "running"
                else:
                    state["status"] = ("cancelled" if cancel else "failed" if "failed" in statuses
                                       else "blocked" if "blocked" in statuses else "finished")
                    state["finished_utc"] = now()
                    cache.pop(run,None)
                if json.dumps(state, sort_keys=True) != old:
                    atomic(run/"state.json", state)
            if once and not active:
                remaining = [p for p in runs.glob("*/state.json") if not p.parent.name.startswith(".")
                             and read(p).get("status") not in TERMINAL]
                if not remaining:
                    return
            time.sleep(interval)


def service(root, verb, identifier, instream, outstream):
    if verb == "doctor":
        return doctor()
    if verb == "register":
        return register_node(root,json.load(instream))
    if verb == "submit":
        return accept(root, identifier, instream)
    run = run_path(root, identifier)
    state = read(run/"state.json")
    if verb == "cancel":
        if state["status"] not in TERMINAL:
            (run/"cancel.request").touch(exist_ok=True)
            return {"id":identifier, "status":"drain_requested", "active_attempts":"finish"}
    if verb in {"status", "cancel"}:
        return state
    if verb == "report":
        receipts = {p.parent.name:read(p) for p in sorted((run/"jobs").glob("*/chemflow-result.json"))}
        return {"run":state,"tasks":receipts, "source_sha256":read(run/"source-sha256.json"),
                "scientific_release":"requires independent claims assessment"}
    if verb == "sync":
        if state["status"] not in TERMINAL:
            raise ValueError("sync requires a terminal run to give a consistent snapshot")
        _pack(run, [p.relative_to(run).as_posix() for p in run.rglob("*") if p.is_file()],outstream)
        return None
    raise ValueError("unknown controller operation")


def doctor():
    # Inventory reports discovery, not method/architecture validation or cloud entitlement.
    cpu = len(os.sched_getaffinity(0)) if hasattr(os,"sched_getaffinity") else os.cpu_count()
    mem = None
    if Path("/proc/meminfo").is_file():
        fields = dict(line.split(":",1) for line in Path("/proc/meminfo").read_text().splitlines())
        mem = int(fields["MemAvailable"].split()[0]) // 1024
    cg = Path("/sys/fs/cgroup")
    if (cg/"memory.max").is_file():
        ceiling = (cg/"memory.max").read_text().strip()
        if ceiling.isdecimal():
            headroom = max(0, int(ceiling)-int((cg/"memory.current").read_text()))
            mem = min(mem,headroom//2**20) if mem is not None else headroom//2**20
    if (cg/"cpu.max").is_file():
        quota, period = (cg/"cpu.max").read_text().split()
        if quota.isdecimal():
            cpu = min(cpu,max(1,int(quota)//int(period)))
    binaries = {}
    for name in ("orca","crest","xtb","sbatch"):
        path = os.environ.get("CHEMFLOW_"+name.upper()) or shutil.which(name)
        if path and Path(path).is_file() and os.access(path,os.X_OK):
            binaries[name]={"path":str(Path(path).resolve()),"sha256":digest(Path(path))}
    return {"platform":platform.system(),"arch":platform.machine(),"available_cores":cpu,
            "available_memory_mb":mem,"scratch_free_mb":shutil.disk_usage(tempfile.gettempdir()).free//2**20,
            "binaries":binaries,"method_qualification":"requires representative run"}
