"""Google Workflows/Batch command line and legacy offline dev2 fixtures."""
import argparse
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile
import uuid

from . import controller as c
from .protocol import load_plan, matches
from .transport import _unpack

MACHINE = "e2-highmem-8"


def check_profile(plan, profile):
    node = dict(profile["node"], placement="controller")
    c.validate_nodes([node])
    if node["cores"] > 8 or node["memory_mb"] > 60 * 1024:
        raise ValueError("profile exceeds e2-highmem-8 usable reservation")
    if node.get("concurrency", 1) != 1:
        raise ValueError("this backend executes one scientific task at a time")
    if any(not matches(t, node) for t in plan["tasks"]):
        raise ValueError("at least one task exceeds the registered image capabilities or capacity")
    return node


def job_spec(plan, profile, run, input_hash):
    node = check_profile(plan, profile)
    image = profile["image"]
    if not re.fullmatch(r"[a-z0-9./_-]+@sha256:[0-9a-f]{64}", image):
        raise ValueError("image must use an immutable sha256 digest")
    if not re.fullmatch(r"[a-z0-9][a-z0-9.-]{1,61}[a-z0-9]", profile["bucket"]):
        raise ValueError("invalid bucket")
    hours = profile["max_hours"]
    if type(hours) is not int or not 1 <= hours <= 168:
        raise ValueError("max_hours must be 1..168 for this backend")
    volume = "/mnt/disks/chemflow"
    args = ["-m","chemflow2.batch","execute",volume,run,input_hash,json.dumps(node)]
    task = {"computeResource":{"cpuMilli":node["cores"]*1000,"memoryMib":node["memory_mb"]},
            "maxRunDuration":f"{hours*3600}s","maxRetryCount":0,
            "volumes":[{"gcs":{"remotePath":profile["bucket"]},"mountPath":volume}],
            "runnables":[{"container":{"imageUri":image,"entrypoint":"python3","commands":args,
                                       "volumes":[f"{volume}:{volume}:rw"]}}]}
    # budget additional disk for VM OS, private image, receipts and retained inputs
    policy = {"machineType":MACHINE,"bootDisk":{"type":"pd-balanced",
              "sizeGb":math.ceil(node["scratch_mb"]*2**20/10**9)+profile["overhead_gb"]}}
    if type(profile["overhead_gb"]) is not int or profile["overhead_gb"] < 32:
        raise ValueError("overhead_gb must reserve at least 32 GB beyond scratch")
    return {"taskGroups":[{"taskCount":1,"taskSpec":task}],
            "allocationPolicy":{"instances":[{"policy":policy}],
                                "serviceAccount":{"email":profile["service_account"]}},
            "logsPolicy":{"destination":"CLOUD_LOGGING"}}


def execute(mount, run, expected, node):
    if not re.fullmatch(r"cf-[a-f0-9]{32}",run):
        raise ValueError("invalid run ID")
    remote = mount/"chemflow"/run
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp).resolve()
        frozen = root/"input.tgz"
        import shutil
        shutil.copyfile(remote/"input.tgz",frozen)
        if c.digest(frozen) != expected:
            raise ValueError("uploaded study checksum mismatch")
        with frozen.open("rb") as src:
            c.accept(root,run,src)
        c.serve(root,[node],once=True,interval=.1)
        report = c.service(root,"report",run,None,None)
        with (remote/"result.tgz").open("wb") as dst:
            c.service(root,"sync",run,None,dst)
        return report["run"]["status"] == "finished"


def gcloud(profile, *args, as_json=False):
    command = ["gcloud", "--quiet", "--project=" + profile["project"], *args]
    if as_json:
        output = subprocess.check_output([*command, "--format=json"], text=True)
        return json.loads(output)
    subprocess.run(command, check=True)


def quota_profile(profile):
    regional = gcloud(profile, "compute", "regions", "describe", profile["region"], as_json=True)
    global_limits = gcloud(profile, "compute", "project-info", "describe", as_json=True)
    def headroom(records, metric):
        rows = [q for q in records.get("quotas", []) if q["metric"] == metric]
        if len(rows) != 1:
            raise ValueError("cannot verify quota for " + metric)
        return max(0, math.floor(float(rows[0]["limit"]) - float(rows[0]["usage"])))
    cpus = headroom(regional, "CPUS")
    if any(q["metric"] == "CPUS_ALL_REGIONS" for q in global_limits.get("quotas", [])):
        cpus = min(cpus, headroom(global_limits, "CPUS_ALL_REGIONS"))
    disk = headroom(regional, "SSD_TOTAL_GB")
    return dict(profile, regional_vcpu_limit=min(cpus, profile.get("regional_vcpu_limit", cpus)),
                regional_disk_limit_gb=min(disk, profile.get("regional_disk_limit_gb", disk)))


def lifecycle(profile, run):
    location = "--location=" + profile["region"]
    workflows = gcloud(profile, "workflows", "list", location, as_json=True)
    present = any(w["name"].rsplit("/", 1)[-1] == run for w in workflows)
    executions = gcloud(profile, "workflows", "executions", "list", run, location, as_json=True) if present else []
    jobs = gcloud(profile, "batch", "jobs", "list", location, "--filter=labels.chemflow-run=" + run, as_json=True)
    jobs = [job for job in jobs if job.get("labels", {}).get("chemflow-run") == run
            and re.fullmatch(re.escape(run) + r"-a\d+", job["name"].rsplit("/", 1)[-1])]
    return {"workflow": present, "executions": executions, "jobs": jobs}


def manage(profile, run, verb, delete_evidence=False):
    report = lifecycle(profile, run)
    if verb == "status":
        print(json.dumps(report, indent=2))
        return
    location = "--location=" + profile["region"]
    active = [e for e in report["executions"] if e["state"] not in {"FAILED", "SUCCEEDED", "CANCELLED"}]
    if verb == "cancel":
        for execution in active:
            gcloud(profile, "workflows", "executions", "cancel", execution["name"], location, "--workflow=" + run)
        report = lifecycle(profile, run)
    for job in report["jobs"]:
        terminal = job["status"]["state"] in {"FAILED", "SUCCEEDED", "CANCELLED"}
        if verb == "cancel" and not terminal:
            gcloud(profile, "batch", "jobs", "cancel", job["name"].rsplit("/", 1)[-1], location)
        elif verb == "cleanup" and not terminal:
            raise ValueError("cancel running allocations before cleanup")
    if verb == "cleanup":
        if active:
            raise ValueError("cancel the active workflow before cleanup")
        for job in report["jobs"]:
            gcloud(profile, "batch", "jobs", "delete", job["name"].rsplit("/", 1)[-1], location)
        if report["workflow"]:
            gcloud(profile, "workflows", "delete", run, location)
        if delete_evidence:
            gcloud(profile, "storage", "rm", "--recursive", f"gs://{profile['bucket']}/chemflow/{run}/")


def collect(downloads, destination, run):
    """Rebuild a source-verifiable run directory, including any partial evidence."""
    import shutil
    from .batch_worker import accept_source, restore

    class LocalStore:
        def get(self, name, stream):
            with (downloads / name).open("rb") as src:
                shutil.copyfileobj(src, stream)

    meta = c.read(downloads / "submission.json")
    if meta["run"] != run:
        raise ValueError("run identity mismatch")
    store = LocalStore()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        raise ValueError("destination already exists")
    with tempfile.TemporaryDirectory(dir=destination.parent) as tmp:
        work = accept_source(store, Path(tmp), run, meta["input_hash"])
        state = c.read(work / "state.json")
        commits = {}
        for key in state["tasks"]:
            marker = downloads / "commits" / (key + ".json")
            status = "unpublished"
            if marker.exists():
                raw = marker.read_bytes()
                result = restore(store, work / "jobs" / key, key, run, meta["input_hash"], raw)
                commits[key] = json.loads(raw)
                status = result["status"]
            state["tasks"][key] = {"status": status}
        from .protocol import _order
        tasks = {t["id"]: t for t in load_plan(work / "source/protocol.json")["tasks"]}
        for key in _order(list(tasks.values())):
            unmet = [dep for dep in tasks[key]["deps"] if state["tasks"][dep]["status"] in ("failed", "blocked")]
            if state["tasks"][key]["status"] == "unpublished" and unmet:
                state["tasks"][key] = {"status": "blocked", "blocked_by": unmet}
        statuses = {row["status"] for row in state["tasks"].values()}
        state["status"] = "finished" if statuses == {"succeeded"} else (
            "failed" if "failed" in statuses else "incomplete")
        c.atomic(work / "state.json", state)
        c.atomic(work / "cloud-commits.json", commits)
        shutil.copyfile(downloads / "submission.json", work / "cloud-submission.json")
        if (downloads / "cloud-status.json").is_file():
            shutil.copyfile(downloads / "cloud-status.json", work / "cloud-status.json")
        work.rename(destination)


def main(argv=None):
    from .batch_dag import encoded, schedule, workflow
    cli = argparse.ArgumentParser(prog="chemflow2 batch")
    subs = cli.add_subparsers(dest="verb", required=True)
    for name in ("plan", "submit", "status", "fetch", "cancel", "cleanup"):
        p = subs.add_parser(name)
        p.add_argument("profile", type=Path)
        if name in ("plan", "submit"):
            p.add_argument("protocol", type=Path)
        else:
            p.add_argument("run")
        if name == "fetch":
            p.add_argument("destination", type=Path)
        if name == "cleanup":
            p.add_argument("--delete-evidence", action="store_true",
                           help="permanently remove all cloud objects for this run")
    # Retain the dev2 private image entry point. New submissions never use it.
    p = subs.add_parser("execute", help=argparse.SUPPRESS)
    p.add_argument("mount", type=Path)
    for key in ("run", "expected", "node"):
        p.add_argument(key)
    args = cli.parse_args(argv)
    if args.verb == "execute":
        return 0 if execute(args.mount, args.run, args.expected, json.loads(args.node)) else 1
    profile = c.read(args.profile)
    run = getattr(args, "run", None) or "cf-" + uuid.uuid4().hex
    if not re.fullmatch(r"cf-[a-f0-9]{32}", run):
        raise ValueError("invalid run ID")
    if args.verb in {"status", "cancel", "cleanup"}:
        manage(profile, run, args.verb, getattr(args, "delete_evidence", False))
        return 0
    prefix = f"gs://{profile['bucket']}/chemflow/{run}/"
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        if args.verb == "fetch":
            if args.destination.exists():
                raise ValueError("destination already exists")
            report = lifecycle(profile, run)
            if any(e["state"] not in {"FAILED", "SUCCEEDED", "CANCELLED"} for e in report["executions"]):
                raise ValueError("fetch requires a terminal workflow")
            if any(job["status"]["state"] not in {"FAILED", "SUCCEEDED", "CANCELLED"} for job in report["jobs"]):
                raise ValueError("fetch requires terminal Batch jobs; cancel active jobs after a workflow failure")
            gcloud(profile, "storage", "rsync", "--recursive", prefix, str(work))
            c.atomic(work / "cloud-status.json", report)
            collect(work, args.destination, run)
            return 0
        submitted_plan_sha256 = c.digest(args.protocol)
        if profile.get("calibration_plan_sha256") is not None and submitted_plan_sha256 != profile["calibration_plan_sha256"]:
            raise ValueError("profile is restricted to its frozen calibration protocol")
        source = work / "input.tgz"
        with source.open("wb") as out:
            c.bundle(args.protocol.resolve(), out)
        if c.digest(args.protocol) != submitted_plan_sha256:
            raise ValueError("protocol changed during freezing")
        input_hash = c.digest(source)
        if args.verb == "submit":
            if (profile["image"].endswith("0" * 64) or profile["project"] == "YOUR_PROJECT"
                    or not profile.get("workflow_service_account")):
                raise ValueError("complete cloud setup before submitting")
            profile = quota_profile(profile)
        with source.open("rb") as frozen:
            c.accept(work, run, frozen)
        frozen_plan = c.run_path(work, run) / "source/protocol.json"
        allocations = schedule(load_plan(frozen_plan), profile, run, input_hash)
        rows = [{k: v for k, v in a.items() if k not in {"spec", "piece"}} for a in allocations]
        if args.verb == "plan":
            print(json.dumps({"backend": "workflows-batch", "allocations": rows}, indent=2))
            return 0
        from .budget import reserve_dispatch
        print(f"Run {run}; storage {prefix}", flush=True)
        reservation = reserve_dispatch(profile, args.profile, run, rows)
        def upload(local, name, content_type):
            gcloud(profile, "storage", "cp", "--if-generation-match=0", "--content-type=" + content_type,
                   str(local), prefix + name)
        meta = {"run": run, "input_hash": input_hash, "allocations": rows, "reservation": reservation,
                "submitted_plan_sha256": submitted_plan_sha256,
                **{key: profile[key] for key in ("project", "region", "bucket", "image")}}
        c.atomic(work / "submission.json", meta)
        upload(work / "submission.json", "submission.json", "application/json")
        upload(source, "input.tgz", "application/gzip")
        for i, allocation in enumerate(allocations):
            for part, folder in (("piece", "pieces"), ("spec", "specs")):
                data = work / "control.json"
                obj = allocation[part] if part == "piece" else {
                    "spec": allocation["spec"], "gates": allocation["gates"]}
                data.write_bytes(encoded(obj))
                upload(data, f"{folder}/{i}.json", "application/json")
        definition = work / "workflow.json"
        c.atomic(definition, workflow(profile, run, len(rows), reservation["valid_until"]))
        location = "--location=" + profile["region"]
        gcloud(profile, "workflows", "deploy", run, location, "--source=" + str(definition),
               "--service-account=" + profile["workflow_service_account"])
        gcloud(profile, "workflows", "execute", run, location)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
