"""Promote only method tags actually exercised by a reviewed private Batch calibration."""
import argparse
import json
import re
import tempfile
from pathlib import Path

from . import batch, controller
from .protocol import load_plan
from .transport import _result
from .worker import read_orca_log, verify_orca


def assess(profile, evidence, run, reviewer, orca_hash):
    image = profile["image"]
    if not re.fullmatch(r"[a-z0-9./_-]+@sha256:[0-9a-f]{64}", image):
        raise ValueError("immutable BYOL image required")
    if not reviewer.strip() or not re.fullmatch(r"[a-f0-9]{64}", orca_hash) or len(set(orca_hash)) == 1:
        raise ValueError("identify the scientific reviewer and independently verified ORCA executable hash")
    meta = controller.read(evidence/"cloud-submission.json")
    if (meta["run"] != run or any(meta[k] != profile[k] for k in ("project", "region", "bucket", "image")) or
            not profile.get("calibration_plan_sha256") or meta["submitted_plan_sha256"] != profile["calibration_plan_sha256"]):
        raise ValueError("calibration identity or submitted plan digest differs")
    state = controller.read(evidence/"state.json")
    status = controller.read(evidence/"cloud-status.json")
    if (state["id"] != run or state["status"] != "finished" or not status["workflow"] or
            len(status["executions"]) != 1 or status["executions"][0]["state"] != "SUCCEEDED"):
        raise ValueError("real successful cloud workflow required before cleanup")
    allocations = {a["id"]: a["task_ids"] for a in meta["allocations"]}
    observed = {}
    for job in status["jobs"]:
        key = job["name"].split("/")[-1]
        groups = job["taskGroups"]
        if (job["status"]["state"] != "SUCCEEDED" or job["labels"].get("chemflow-run") != run
                or not job.get("uid") or len(groups) != 1 or key in observed
                or job["allocationPolicy"]["serviceAccount"]["email"] != profile["service_account"]):
            raise ValueError("Batch allocation identity or status differs")
        containers = groups[0]["taskSpec"]["runnables"]
        if len(containers) != 1:
            raise ValueError("expected exactly one production worker runnable")
        container = containers[0]["container"]
        if (container["imageUri"] != image or container["entrypoint"] != "python3"
                or container["commands"][:4] != ["-m", "chemflow2.batch_worker", profile["bucket"], run]):
            raise ValueError("job did not execute the expected private worker")
        observed[key] = job["uid"]
    if not allocations or allocations.keys() != observed.keys():
        raise ValueError("incomplete set of successful Batch allocations")
    plan = load_plan(evidence/"source/protocol.json")
    ids = [task for group in allocations.values() for task in group]
    if (plan.get("study_scope") != "engineering_calibration" or len(set(ids)) != len(ids)
            or set(ids) != {t["id"] for t in plan["tasks"]}):
        raise ValueError("incomplete frozen engineering calibration")
    records = {}
    artifacts = []
    for task in plan["tasks"]:
        key = task["id"]
        if task["engine"] != "orca" or task.get("study_scope") != "engineering_calibration":
            raise ValueError("only representative ORCA calibration tasks can qualify methods")
        folder = evidence/"jobs"/key
        receipt = _result(folder)
        prov = receipt["provenance"]
        log = read_orca_log(folder/"orca.out")
        deck = folder/task["input"]
        if (not receipt["ok"] or prov.get("task_id") != key or prov.get("engine") != "orca" or
                prov.get("stage") != task["stage"] or prov.get("platform") != "Linux" or
                prov.get("architecture") != "x86_64" or prov.get("executable") != "/opt/orca/orca" or
                prov.get("executable_sha256") != orca_hash or prov.get("exit_code") != 0 or
                prov.get("orca_processes") != task["resources"]["cores"] or
                prov.get("input_sha256", {}).get(task["input"]) != controller.digest(deck) or
                not re.fullmatch(r"6\.\d+\.\d+", log.get("engine_version", "")) or
                not receipt["metrics"].get("wall_seconds", 0) > 0 or "orca.out" not in receipt["artifacts"]):
            raise ValueError("unverified ORCA executable, task provenance, runtime or log")
        verify_orca(task["stage"], log, deck, folder, task["outputs"])
        for tag in task["requires"]:
            records.setdefault(tag, []).append(key)
        artifacts.append({"id":key, "stage":task["stage"], "orca_version":log["engine_version"],
                          "log_sha256":receipt["artifacts"]["orca.out"],
                          "receipt_sha256":controller.digest(folder/"chemflow-result.json")})
    if not records or set(records) != set(profile["node"]["capabilities"]):
        raise ValueError("tags must match the prospective calibration profile exactly")
    qualified = dict(profile)
    qualified.pop("calibration_plan_sha256")
    qualified["qualification"] = {"reviewer":reviewer.strip(), "reviewed_at":controller.now(),
        "run":run, "project":meta["project"], "image":image, "input_sha256":meta["input_hash"],
        "submitted_plan_sha256":meta["submitted_plan_sha256"], "orca_sha256":orca_hash,
        "cloud_job_uids":observed, "tasks_by_capability":records, "evidence":artifacts}
    return qualified


def main(argv=None):
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("profile", type=Path)
    cli.add_argument("run")
    cli.add_argument("destination", type=Path)
    cli.add_argument("--reviewer", required=True)
    cli.add_argument("--authorized-orca-sha256", required=True,
                     help="SHA256 of the ORCA executable obtained from the authorized distribution")
    args = cli.parse_args(argv)
    if args.destination.exists():
        raise ValueError("provide a new output path")
    profile = controller.read(args.profile)
    with tempfile.TemporaryDirectory() as temp:
        evidence = Path(temp)/"evidence"
        # Never qualify caller-supplied fixture directories. Fetch via the signed-in account.
        batch.main(["fetch", str(args.profile), args.run, str(evidence)])
        qualified = assess(profile, evidence, args.run, args.reviewer, args.authorized_orca_sha256)
    with args.destination.open("x") as output:
        json.dump(qualified, output, indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
