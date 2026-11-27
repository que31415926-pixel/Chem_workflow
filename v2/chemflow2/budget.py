"""Conservative USD admission for an isolated, reviewed trial project.

The account owner supplies dated provider evidence. A generation-conditional
GCS ledger holds reservations across independent laptops and workflow runs.
"""
import argparse
import csv
from datetime import datetime, timedelta, timezone
from decimal import Decimal, ROUND_CEILING
import hashlib
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

EXTRAS = ("storage", "image", "logging", "network", "operations", "failure", "cleanup")
CENT = Decimal("0.01")


def now():
    return datetime.now(timezone.utc)


def timestamp(value):
    if not isinstance(value, str):
        raise ValueError("missing timestamp for budget validity")
    result = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if result.tzinfo is None:
        raise ValueError("budget timestamps require a timezone")
    return result.astimezone(timezone.utc)


def money(value):
    try:
        number = Decimal(str(value))
    except ArithmeticError as exc:
        raise ValueError("invalid budget value") from exc
    if not number.is_finite() or number < 0 or number > 10 ** 12:
        raise ValueError("invalid negative, nonfinite or out-of-range budget value")
    return number


def usd(value):
    return str(value.quantize(CENT, rounding=ROUND_CEILING))


def evidence(directory, record):
    path = Path(record["file"])
    absolute = (directory / path).resolve()
    if path.is_absolute() or not absolute.is_relative_to(directory.resolve()):
        raise ValueError("evidence must be local to its policy directory")
    data = absolute.read_bytes()
    if not data or hashlib.sha256(data).hexdigest() != record["sha256"]:
        raise ValueError("missing or changed provider evidence")
    return data


def policy_for(profile, profile_path):
    path = Path(profile_path).parent / profile["budget_policy"]
    policy = json.loads(path.read_bytes(), parse_float=Decimal)
    if any(policy[k] != profile[k] for k in ("project", "bucket", "region")):
        raise ValueError("budget policy must match project, canonical bucket, and region")
    if policy["currency"] != "USD":
        raise ValueError("USD tariff and credits required")
    for k in ("provider_blocks_paid_usage", "trial_eligibility_confirmed",
              "dedicated_account", "billing_current", "retention_and_cleanup_reviewed"):
        if policy.get(k) is not True:
            raise ValueError("unconfirmed budget prerequisite: " + k)
    for key in ("account", "tariff", "other_cost_bounds"):
        evidence(path.parent, policy["evidence"][key])
    if not timestamp(policy["observed_at"]) <= now() < timestamp(policy["valid_until"]):
        raise ValueError("stale or future billing evidence")
    return policy


def quote(policy, allocations):
    if not allocations:
        raise ValueError("cannot quote an empty dispatch")
    overhead = money(policy["vm_overhead_hours"])
    disk_rate = money(policy["disk_usd_per_gb_hour"])
    if overhead <= 0 or disk_rate <= 0:
        raise ValueError("disk tariff and provisioning/teardown allowance required")
    totals = {"vm": Decimal(0), "disk": Decimal(0)}
    details = []
    seen = set()
    longest = Decimal(0)
    for row in allocations:
        ids = row["task_ids"]
        if (not isinstance(ids, list) or not ids or not all(isinstance(i,str) and i for i in ids)
                or len(set(ids)) != len(ids) or seen.intersection(ids)):
            raise ValueError("each task must be allocated exactly once")
        seen.update(ids)
        n = row["max_instances"]
        hours = money(row["max_hours"])
        disk = money(row["disk_gb"])
        machine = row["machine_type"]
        rate = money(policy["vm_usd_per_hour"][machine])
        if type(n) is not int or n < 1 or min(hours, disk, rate) <= 0:
            raise ValueError("positive instance count, runtime, disk and tariff required")
        longest = max(longest, hours+overhead)
        instance_hours = n * (hours + overhead)
        vm = instance_hours * rate
        drive = instance_hours * disk * disk_rate
        totals["vm"] += vm
        totals["disk"] += drive
        details.append({"id":row["id"], "task_ids":ids, "machine_type":machine,
                        "max_instances":n, "disk_gb":str(disk), "hours_per_slot":str(hours),
                        "vm_and_disk_upper_usd":usd(vm + drive)})
    for k in EXTRAS:
        totals[k] = money(policy["ancillary_upper_usd"][k])
    reserve = sum(totals.values())
    # End date also covers post-job evidence retention, deletion and invoicing tail.
    tail = money(policy["tail_hours"])
    deadline = timestamp(policy["valid_until"]) - timedelta(seconds=float((longest + tail) * 3600))
    return {"allocations":details, "itemized_usd":{k:usd(v) for k,v in totals.items()},
            "upper_usd":usd(reserve), "valid_until":deadline.isoformat().replace("+00:00","Z")}


def admission(policy, ledger, run_id, allocation_quote):
    if timestamp(allocation_quote["valid_until"]) <= now():
        raise ValueError("insufficient verified credit lifetime for a dispatch")
    if ledger["project"] != policy["project"] or ledger["currency"] != "USD":
        raise ValueError("unexpected ledger identity")
    holds = ledger["holds"]
    if run_id in holds:
        if holds[run_id] != allocation_quote:
            raise ValueError("run identifier already reserved for different inputs")
        return
    available = min(money(policy["credit_remaining_usd"]), money(policy["spend_limit_usd"]))
    outstanding = sum((money(v["upper_usd"]) for v in holds.values()), Decimal(0))
    deficit = outstanding + money(allocation_quote["upper_usd"]) - available
    if deficit > 0:
        ids = [i for row in allocation_quote["allocations"] for i in row["task_ids"]]
        raise ValueError(f"budget shortfall USD {usd(deficit)} for dispatch of {', '.join(ids)}")
    holds[run_id] = allocation_quote


def reserve_dispatch(profile, profile_path, run_id, allocations):
    policy = policy_for(profile, profile_path)
    bound = quote(policy, allocations)
    obj = "gs://" + profile["bucket"] + "/chemflow-budget/v1/ledger.json"
    cli = ["gcloud", "--project=" + profile["project"], "storage"]
    with TemporaryDirectory() as scratch:
        local = Path(scratch)/"ledger.json"
        # The create-only operation cannot overwrite an existing record.
        initial = {"project":profile["project"], "currency":"USD", "holds":{}}
        admission(policy, dict(initial, holds={}), run_id, bound)
        local.write_text(json.dumps(initial))
        subprocess.run([*cli,"cp","--if-generation-match=0",str(local),obj],
                       capture_output=True, check=False)
        gen = subprocess.check_output([*cli,"objects","describe",obj,"--format=value(generation)"],
                                      text=True).strip()
        if not gen.isdecimal() or int(gen) <= 0:
            raise ValueError("cannot establish current shared budget state")
        subprocess.run([*cli,"cp",obj+"#"+gen,str(local)],check=True)
        ledger = json.loads(local.read_bytes())
        admission(policy, ledger, run_id, bound)
        local.write_text(json.dumps(ledger))
        # A race or unavailable billing state aborts submit; no unlocked fallback.
        subprocess.run([*cli,"cp","--if-generation-match="+gen,str(local),obj],check=True)
    return bound


def reconcile(path):
    """Reconcile immutable worker receipts with reviewed attributable billing rows."""
    report = json.loads(path.read_bytes())
    root = path.parent
    export = evidence(root, report["billing_export"])
    rows = list(csv.DictReader(export.decode("utf-8-sig").splitlines()))
    total = Decimal(0)
    credited = Decimal(0)
    seen = set()
    for row in rows:
        if (row["run_id"] != report["run_id"] or row["project"] != report["project"] or
                row["currency"] != "USD" or not row["line_id"] or row["line_id"] in seen):
            raise ValueError("mismatched billing scope or duplicate line")
        seen.add(row["line_id"])
        total += money(row["gross_usd"])
        credited += money(row["credit_usd"])
    if not rows or credited > total:
        raise ValueError("missing or inconsistent billing rows")
    tasks = []
    receipt_ids = set()
    for ref in report["receipts"]:
        receipt = json.loads(evidence(root, ref))
        task_id = receipt["task_id"]
        if receipt.get("run_id") != report["run_id"] or task_id in receipt_ids:
            raise ValueError("duplicate task or receipt from a different cloud run")
        receipt_ids.add(task_id)
        metrics = receipt.get("metrics")
        provenance = receipt["provenance"]
        tasks.append({"task":task_id, "stage":provenance.get("stage"),
                      "engine_version":receipt["observations"].get("engine_version"),
                      "status":receipt["status"], "failure_kind":receipt.get("failure_kind"),
                      "metrics":metrics,"sha256":ref["sha256"]})
    if not tasks:
        raise ValueError("at least one attributable receipt required")
    return {"project":report["project"],"run_id":report["run_id"],"receipts":tasks,
            "billing_export_sha256":report["billing_export"]["sha256"],
            "reported_gross_usd":str(total),"reported_credit_usd":str(credited),
            "reported_net_usd":str(total - credited), "billing_lines":len(rows)}


def settle(ledger, settlements, approved_by, evidence_note, finished, when=None):
    """Replace finished runs' upper-bound holds by audited settled amounts. Pure; returns (ledger, record).

    Holds are never deleted: each keeps its original bound and an audit block. A run still in progress,
    an unknown run, a second settlement, or a settled amount above the original bound is refused.
    """
    if not approved_by.strip() or not evidence_note.strip():
        raise ValueError("settlement needs the approving account holder and the billing evidence")
    holds = ledger["holds"]
    stamp = when or now().isoformat().replace("+00:00", "Z")
    changes = []
    for run, amount in settlements.items():
        hold = holds.get(run)
        if hold is None:
            raise ValueError("no hold for " + run)
        if hold.get("settled"):
            raise ValueError(run + " is already settled")
        if run not in finished:
            raise ValueError(run + " is not confirmed finished; refusing to settle a live run")
        value = money(amount)
        original = money(hold["upper_usd"])
        if value <= 0 or value > original:
            raise ValueError(f"settled amount for {run} must be positive and not exceed {usd(original)}")
        hold.update({"original_upper_usd": hold["upper_usd"], "upper_usd": usd(value), "settled": True,
                     "settlement_audit": {"approved_by": approved_by, "evidence": evidence_note,
                                          "settled_utc": stamp}})
        changes.append({"run": run, "original_upper_usd": usd(original), "settled_upper_usd": usd(value)})
    total = sum((money(v["upper_usd"]) for v in holds.values()), Decimal(0))
    return ledger, {"changes": changes, "holds_total_after_usd": usd(total), "approved_by": approved_by,
                    "evidence": evidence_note, "settled_utc": stamp}


def finished_runs(profile, runs):
    """Runs whose latest workflow execution is no longer ACTIVE (cloud builds count as finished)."""
    done = set()
    for run in runs:
        if run.startswith("build-"):
            done.add(run)
            continue
        state = subprocess.run(["gcloud", "--project=" + profile["project"], "workflows", "executions", "list",
                                run, "--location=" + profile["region"], "--limit=1", "--format=value(state)"],
                               capture_output=True, text=True)
        if state.returncode == 0 and state.stdout.strip() in {"SUCCEEDED", "FAILED", "CANCELLED"}:
            done.add(run)
        elif state.returncode != 0 and "NOT_FOUND" in state.stderr:
            done.add(run)  # workflow never deployed: nothing can still be running
    return done


def settle_in_bucket(profile, settlements, approved_by, evidence_note, backup_dir):
    obj = "gs://" + profile["bucket"] + "/chemflow-budget/v1/ledger.json"
    cli = ["gcloud", "--project=" + profile["project"], "storage"]
    backup_dir = Path(backup_dir)
    backup_dir.mkdir(parents=True, exist_ok=True)
    gen = subprocess.check_output([*cli, "objects", "describe", obj, "--format=value(generation)"],
                                  text=True).strip()
    if not gen.isdecimal():
        raise ValueError("cannot establish current shared budget state")
    backup = backup_dir / f"ledger_gen{gen}.json"
    subprocess.run([*cli, "cp", obj + "#" + gen, str(backup)], check=True)
    ledger = json.loads(backup.read_bytes())
    ledger, record = settle(ledger, settlements, approved_by, evidence_note,
                            finished_runs(profile, settlements))
    new = backup_dir / f"ledger_after_gen{gen}.json"
    new.write_text(json.dumps(ledger))
    subprocess.run([*cli, "cp", "--if-generation-match=" + gen, str(new), obj], check=True)
    record.update(ledger_generation_before=gen, backup=str(backup))
    (backup_dir / f"settlement_gen{gen}.json").write_text(json.dumps(record, indent=2, ensure_ascii=False))
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    verbs = parser.add_subparsers(dest="verb", required=True)
    q = verbs.add_parser("quote")
    q.add_argument("profile",type=Path)
    q.add_argument("allocations",type=Path)
    r = verbs.add_parser("reconcile")
    r.add_argument("manifest",type=Path)
    st = verbs.add_parser("settle", help="settle finished runs' holds after the account holder audits billing")
    st.add_argument("profile", type=Path)
    st.add_argument("settlement", nargs="+", help="RUN=USD, e.g. cf-0123...=3.00 (finished runs only)")
    st.add_argument("--approved-by", required=True)
    st.add_argument("--evidence", required=True, help="what the holder checked, e.g. console credit and date")
    st.add_argument("--backup-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.verb == "settle":
        pairs = dict(item.split("=", 1) for item in args.settlement)
        profile = json.loads(args.profile.read_bytes())
        print(json.dumps(settle_in_bucket(profile, pairs, args.approved_by, args.evidence, args.backup_dir),
                         indent=2, ensure_ascii=False))
        return
    if args.verb == "quote":
        profile = json.loads(args.profile.read_bytes())
        allocations = json.loads(args.allocations.read_bytes())
        if isinstance(allocations, dict):
            allocations = allocations["allocations"]
        policy = policy_for(profile, args.profile)
        result = quote(policy, allocations)
        available = min(money(policy["credit_remaining_usd"]), money(policy["spend_limit_usd"]))
        result["shortfall_before_existing_holds_usd"] = usd(max(Decimal(0), money(result["upper_usd"]) - available))
        result["dispatch_deadline_passed"] = timestamp(result["valid_until"]) <= now()
    else:
        result = reconcile(args.manifest)
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
