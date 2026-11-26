import io
import json
from datetime import datetime, timezone
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from chemflow2 import batch, batch_dag, batch_worker, controller
from chemflow2.protocol import load_plan

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "examples/probe/protocol.json"
RUN = "cf-" + "a" * 32


class LocalBucket:
    def __init__(self, directory):
        self.directory = directory

    def create(self, name, stream):
        path = self.directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        try:
            with path.open("xb") as dest:
                shutil.copyfileobj(stream, dest)
        except FileExistsError:
            return False
        return True

    def get(self, name, stream):
        with (self.directory / name).open("rb") as src:
            shutil.copyfileobj(src, stream)

    def read(self, name):
        path = self.directory / name
        return path.read_bytes() if path.exists() else None


class BatchDagTests(unittest.TestCase):
    def setUp(self):
        self.profile = controller.read(ROOT / "examples/batch-profile.example.json")
        self.plan = load_plan(SOURCE)

    def test_dimensions_packing_gates_and_quota(self):
        tasks = self.plan["tasks"]
        for i in range(5):
            tasks.append(dict(tasks[0], id="cheap-" + str(i)))
        tasks.append(dict(tasks[0], id="large", resources={"cores": 4, "memory_mb": 32000, "scratch_mb": 10240}))
        self.profile.update(tasks_per_slot=2, max_parallel_vms=2)
        self.profile["node"]["engines"].append("analysis")
        tasks.append(dict(tasks[0], id="gate", engine="analysis", deps=["consume"]))
        tasks.append(dict(tasks[0], id="gated", deps=["gate"]))
        tasks.append(dict(tasks[0], id="child-of-gated", deps=["gated"]))
        rows = batch_dag.schedule(self.plan, self.profile, RUN, "a" * 64)
        lookup = {k: row for row in rows for k in row["task_ids"]}
        self.assertEqual(len(lookup), len(tasks))
        self.assertLess(lookup["produce"]["level"], lookup["consume"]["level"])
        self.assertNotEqual(lookup["large"]["machine_type"], lookup["produce"]["machine_type"])
        self.assertEqual(lookup["child-of-gated"]["gates"], ["gate"])
        self.assertTrue(any(row["max_instances"] == 2 for row in rows))
        self.assertTrue(all(row["max_instances"] <= 2 for row in rows))
        self.assertTrue(all(row["spec"]["taskGroups"][0]["taskCountPerNode"] == 1 for row in rows))
        self.profile["node"]["memory_mb"] = 12000
        with self.assertRaisesRegex(ValueError, "large"):
            batch_dag.schedule(self.plan, self.profile, RUN, "a" * 64)

    def test_isolated_slot_publication_and_partial_fetch(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            store = LocalBucket(base)
            wire = io.BytesIO()
            controller.bundle(SOURCE, wire)
            raw = wire.getvalue()
            digest = batch_dag.sha(raw)
            store.create("input.tgz", io.BytesIO(raw))
            controller.atomic(base / "submission.json", {"run": RUN, "input_hash": digest})
            slots = batch_dag.schedule(self.plan, self.profile, RUN, digest)
            with self.assertRaisesRegex(ValueError, "no committed"):
                batch_worker.run_slot(store, RUN, slots[1]["piece"], 0)
            with self.assertRaisesRegex(ValueError, "authorization expired"):
                batch_worker.run_slot(store, RUN, slots[0]["piece"], 0, datetime(2000, 1, 1, tzinfo=timezone.utc))
            self.assertTrue(batch_worker.run_slot(store, RUN, slots[0]["piece"], 0))
            commit = store.read("commits/produce.json")
            self.assertTrue(batch_worker.run_slot(store, RUN, slots[0]["piece"], 0))
            self.assertEqual(commit, store.read("commits/produce.json"))
            batch.collect(base, base / "partial", RUN)
            self.assertEqual(controller.read(base / "partial/state.json")["status"], "incomplete")
            self.assertTrue(batch_worker.run_slot(store, RUN, slots[1]["piece"], 0))
            batch.collect(base, base / "complete", RUN)
            self.assertEqual(controller.read(base / "complete/state.json")["status"], "finished")
            self.assertEqual(controller.digest(base / "complete/jobs/produce/sample.json"),
                             controller.digest(base / "complete/jobs/consume/parent.json"))
            attempts = list((base / "attempts/produce").iterdir())
            self.assertEqual(len(attempts), 1)
            attempts[0].write_bytes(b"altered")
            with self.assertRaisesRegex(ValueError, "checksum"):
                batch.collect(base, base / "broken", RUN)
            self.assertFalse((base / "broken").exists())

    def test_scientific_failure_commits_and_blocks_consumer(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            store = LocalBucket(root)
            with (root / "input.tgz").open("wb") as out:
                controller.bundle(SOURCE, out)
            digest = controller.digest(root / "input.tgz")
            slots = batch_dag.schedule(self.plan, self.profile, RUN, digest)
            def negative(node, task, path):
                return {"ok": False, "status": "failed", "failure_kind": "scientific_fail", "artifacts": {}}
            with patch.object(controller, "registered_worker", negative):
                self.assertTrue(batch_worker.run_slot(store, RUN, slots[0]["piece"], 0))
            self.assertFalse(json.loads(store.read("commits/produce.json"))["ok"])
            with self.assertRaisesRegex(ValueError, "required task"):
                batch_worker.run_slot(store, RUN, slots[1]["piece"], 0)
            controller.atomic(root / "submission.json", {"run": RUN, "input_hash": digest})
            batch.collect(root, root / "result", RUN)
            self.assertEqual(controller.read(root / "result/state.json")["tasks"]["consume"]["status"], "blocked")

    def test_workflow_guards_precede_billable_job(self):
        wf = batch_dag.workflow(self.profile, RUN, 3, "2030-01-01T00:00:00Z")
        loop = wf["main"]["steps"][0]["allocate"]["for"]
        self.assertEqual(loop["range"], [0, 2])
        steps = [next(iter(step)) for step in loop["steps"]]
        self.assertLess(steps.index("guards"), steps.index("submit"))
        self.assertLess(steps.index("checkAuthorization"), steps.index("submit"))

    def test_submit_uploads_frozen_inputs_control_and_workflow_only_after_admission(self):
        with tempfile.TemporaryDirectory() as tmp:
            profile = Path(tmp) / "profile.json"
            self.profile.update(project="research-project", workflow_service_account="wf@research-project",
                                image="example.io/private/worker@sha256:" + "b" * 64)
            controller.atomic(profile, self.profile)
            controls = {}
            commands = []
            def cli(config, *args, as_json=False):
                commands.append(args)
                if args[:3] == ("compute", "regions", "describe"):
                    return {"quotas": [{"metric": "SSD_TOTAL_GB", "limit": 256, "usage": 5},
                                       {"metric": "CPUS", "limit": 8, "usage": 2}]}
                if args[:3] == ("compute", "project-info", "describe"):
                    return {"quotas": [{"metric": "CPUS_ALL_REGIONS", "limit": 4, "usage": 1}]}
                if args[:2] == ("storage", "cp"):
                    self.assertIn("--if-generation-match=0", args)
                    controls[args[-1].split("/chemflow/")[1].split("/", 1)[1]] = Path(args[-2]).read_bytes()
            from contextlib import redirect_stdout
            with patch.object(batch, "gcloud", side_effect=cli), redirect_stdout(io.StringIO()):
                with patch("chemflow2.budget.reserve_dispatch", side_effect=ValueError("insufficient credit")):
                    with self.assertRaisesRegex(ValueError, "credit"):
                        batch.main(["submit", str(profile), str(SOURCE)])
                self.assertFalse(controls)
                with patch("chemflow2.budget.reserve_dispatch", return_value={"valid_until": "2030-01-01T00:00:00Z"}):
                    self.assertEqual(batch.main(["submit", str(profile), str(SOURCE)]), 0)
            self.assertEqual(commands[-2][:2], ("workflows", "deploy"))
            self.assertEqual(commands[-1][:2], ("workflows", "execute"))
            for index in (0, 1):
                spec = json.loads(controls[f"specs/{index}.json"])["spec"]
                task = spec["taskGroups"][0]["taskSpec"]
                self.assertEqual(task["runnables"][0]["container"]["commands"][-1],
                                 batch_dag.sha(controls[f"pieces/{index}.json"]))
            self.assertEqual(batch_dag.sha(controls["input.tgz"]),
                             json.loads(controls["submission.json"])["input_hash"])

    def test_regional_limits_apply_to_whole_vm_and_disk(self):
        for i in range(5):
            self.plan["tasks"].append(dict(self.plan["tasks"][0], id="probe-" + str(i)))
        self.profile.update(max_parallel_vms=8, tasks_per_slot=1, regional_vcpu_limit=6, regional_disk_limit_gb=160,
                            license_parallelism={"probe": 1})
        rows = batch_dag.schedule(self.plan, self.profile, RUN, "a"*64)
        self.assertTrue(all(r["max_instances"] == 1 for r in rows))
        del self.profile["license_parallelism"]
        rows = batch_dag.schedule(self.plan, self.profile, RUN, "a"*64)
        self.assertEqual(max(r["max_instances"] for r in rows), 2)
        self.profile["regional_vcpu_limit"] = 1
        with self.assertRaisesRegex(ValueError, "regional quota"):
            batch_dag.schedule(self.plan, self.profile, RUN, "a"*64)

    def test_cleanup_preserves_evidence_and_refuses_active_jobs(self):
        info = {"workflow": True, "executions": [],
                "jobs": [{"name": RUN + "-a0", "status": {"state": "RUNNING"}}]}
        with patch.object(batch, "lifecycle", return_value=info), patch.object(batch, "gcloud") as cli:
            with self.assertRaisesRegex(ValueError, "cancel running"):
                batch.manage(self.profile, RUN, "cleanup", True)
            cli.assert_not_called()
            info["jobs"][0]["status"]["state"] = "CANCELLED"
            batch.manage(self.profile, RUN, "cleanup", False)
            self.assertFalse(any(c.args[1] == "storage" for c in cli.call_args_list))
            batch.manage(self.profile, RUN, "cleanup", True)
            self.assertEqual(cli.call_args.args[-1], "gs://" + self.profile["bucket"] + "/chemflow/" + RUN + "/")

    def test_order_edges_use_only_verified_commits_and_data_edges_restore_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            self.plan["tasks"].append(dict(self.plan["tasks"][0], id="barrier", deps=["produce"]))
            origin = base / "plan.json"
            controller.atomic(origin, self.plan)
            store = LocalBucket(base)
            with (base / "input.tgz").open("wb") as wire:
                controller.bundle(origin, wire)
            digest = controller.digest(base / "input.tgz")
            self.profile["tasks_per_slot"] = 1
            pieces = {row["task_ids"][0]: row["piece"] for row in
                      batch_dag.schedule(load_plan(origin), self.profile, RUN, digest)}
            batch_worker.run_slot(store, RUN, pieces["produce"], 0)
            with patch.object(store, "get", wraps=store.get) as download:
                batch_worker.run_slot(store, RUN, pieces["barrier"], 0)
            self.assertEqual([call.args[0] for call in download.call_args_list], ["input.tgz"])
            marker = json.loads(store.read("commits/produce.json"))
            (base / "attempts/produce" / (marker["attempt"] + ".tgz")).write_bytes(b"corrupt")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                batch_worker.run_slot(store, RUN, pieces["consume"], 0)
            (base / "commits/barrier.json").unlink()
            marker["task_id"] = "unrelated"
            (base / "commits/produce.json").write_text(json.dumps(marker))
            with self.assertRaisesRegex(ValueError, "identity mismatch"):
                batch_worker.run_slot(store, RUN, pieces["barrier"], 0)
