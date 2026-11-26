import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from chemflow2.capacity import estimate


class CapacityTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.plan = json.loads((Path(__file__).resolve().parents[1] / "examples/probe/protocol.json").read_text())
        self.nodes = [{"name":"test", "host":"test", "engines":["probe"], "cores":1, "memory_mb":64,
                       "scratch_mb":32, "capacity_class":"linux-test-1", "remaining_reserved_core_hours":1}]
        self.history = {"samples":[], "engine_versions":{"probe":"test-only"}}
        self.receipts = []
        for task, wall in zip(self.plan["tasks"], [3600, 1800]):
            self.receipts.append({"node":"test", "ok":True, "status":"succeeded", "node_capacity_class":"linux-test-1",
                                  "provenance":{"engine":"probe", "stage":"probe", "task_id":task["id"],
                                                "reservations":{"cores":1}}, "observations":{"engine_version":"test-only"},
                                  "metrics":{"wall_seconds":wall, "scratch_peak_bytes":None, "child_max_rss_bytes":1024}})
            self.history["samples"].append({})

    def calculation(self):
        for i, (row, record) in enumerate(zip(self.history["samples"], self.receipts)):
            p = self.root / f"receipt-{i}.json"
            p.write_text(json.dumps(record))
            row.update(receipt=p.name, sha256=hashlib.sha256(p.read_bytes()).hexdigest())
        paths = []
        for name, obj in (("plan", self.plan), ("hist", self.history), ("nodes",self.nodes)):
            path = self.root / (name + ".json")
            path.write_text(json.dumps(obj))
            paths.append(path)
        return estimate(*paths)

    def test_measured_chain_and_explicit_quota_shortfall(self):
        result = self.calculation()
        self.assertEqual(result["status"], "estimated")
        self.assertEqual(result["critical_path_seconds"], 5400)
        self.assertEqual(result["conditional_min_elapsed_seconds"], 5400)
        self.assertEqual(result["total_reserved_core_hours"], 1.5)
        self.assertEqual(result["node_allocation"]["test"]["assigned_quota_shortfall_core_hours"], .5)
        self.assertIsNone(result["task_profiles"]["consume"]["profiles"][0]["observed_envelopes"]["scratch_peak_bytes"])
        self.assertIsNone(result["actual_cost"])

    def test_absent_failed_or_mismatched_measurements_remain_unknown(self):
        for mutate in (lambda r:r.update(ok=False, status="failed"), lambda r:r["observations"].update(engine_version="different"),
                       lambda r:r.update(node_capacity_class="older-shape"), lambda r:r["provenance"]["reservations"].update(cores=4)):
            previous = json.loads(json.dumps(self.receipts[0]))
            mutate(self.receipts[0])
            result = self.calculation()
            self.assertEqual(result["unknown_tasks"], ["produce"])
            self.assertIsNone(result["critical_path_seconds"])
            self.assertIsNone(result["total_reserved_core_hours"])
            self.receipts[0] = previous

    def test_observed_memory_envelope_rejects_underprovisioned_node(self):
        self.receipts[0]["metrics"]["child_max_rss_bytes"] = 128 * 2**20
        row = self.calculation()["task_profiles"]["produce"]["profiles"][0]
        self.assertEqual(row["deficits"]["memory_mb"], 64)
        self.assertIsNone(row["observed_envelopes"]["job_peak_memory_bytes"])
        self.assertTrue(any("exceeds task reservation" in reason for reason in row["reasons"]))

    def test_operator_confirmation_of_comparable_task_and_unchanged_host(self):
        self.nodes[0].pop("capacity_class")
        self.receipts[0]["provenance"]["task_id"] = "pilot"
        for row in self.history["samples"]:
            row["same_hardware_confirmed"] = True
        result = self.calculation()
        self.assertEqual(result["unknown_tasks"], ["produce"])
        self.history["samples"][0].update(comparable_to=["produce"], rationale="Reviewed identical fixture and resources")
        self.assertEqual(self.calculation()["unknown_tasks"], [])

    def test_damaged_evidence_refused_before_estimation(self):
        self.calculation()
        (self.root / "receipt-0.json").write_text("modified")
        with self.assertRaisesRegex(ValueError, "digest"):
            estimate(self.root/"plan.json", self.root/"hist.json", self.root/"nodes.json")
