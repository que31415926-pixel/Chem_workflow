"""Regression tests for the 2026-11-26 K2 night: rules that were stricter than the science needed."""
import json
import os
from pathlib import Path
import tempfile
import unittest

from chemflow2 import batch_dag, budget, transport
from chemflow2.worker import run_task

from test_batch_dag import BatchDagTests, RUN
from test_worker_adapters import ExtraAdapters


class SymlinkedTempRoot(unittest.TestCase):
    def test_symlinked_root_is_accepted_but_inner_symlink_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            real = Path(tmp) / "private" / "work"
            real.mkdir(parents=True)
            alias = Path(tmp) / "var"
            alias.symlink_to(Path(tmp) / "private")
            (real / "protocol.json").write_text("{}")
            self.assertEqual(transport._file(alias / "work", "protocol.json"), alias / "work" / "protocol.json")
            (real / "link.json").symlink_to(real / "protocol.json")
            with self.assertRaisesRegex(transport.TransportError, "Symlink"):
                transport._file(alias / "work", "link.json")


class CrestEnergyRounding(ExtraAdapters):
    def test_rounding_inversion_is_kept_and_recorded(self):
        self.set_crest(energies="1 0.000\n2 -0.001\n")
        report = run_task(self.crest_task(), self.case)
        self.assertTrue(report["ok"], report)
        order = report["observations"]["energy_order"]
        self.assertEqual(order["out_of_order_frames"], [2])
        self.assertTrue(order["within_rounding"])

    def test_monotonic_output_has_no_flag(self):
        self.set_crest(energies="1 0.000\n2 0.000\n")
        report = run_task(self.crest_task(), self.case)
        self.assertTrue(report["ok"], report)
        self.assertNotIn("energy_order", report["observations"])

    def test_negative_relative_energy_beyond_rounding_still_rejected(self):
        self.set_crest(energies="1 0.000\n2 -0.500\n")
        report = run_task(self.crest_task(), self.case)
        self.assertEqual(report["failure_kind"], "validation")


class WorkflowAndQuota(BatchDagTests):
    def test_failed_job_does_not_stop_independent_allocations(self):
        wf = json.dumps(batch_dag.workflow(self.profile, RUN, 3, "2030-01-01T00:00:00Z"))
        self.assertIn('job.status.state == \\"SUCCEEDED\\" or job.status.state == \\"FAILED\\"', wf)
        self.assertIn("cancelled or deleted", wf)
        self.assertNotIn("Batch job failed;", wf)

    def test_failed_job_error_from_connector_is_caught(self):
        wf = batch_dag.workflow(self.profile, RUN, 2, "2030-01-01T00:00:00Z")
        steps = wf["main"]["steps"][0]["allocate"]["for"]["steps"]
        submit = next(s["submit"] for s in steps if "submit" in s)
        self.assertIn("try", submit)
        handler = submit["except"]["steps"]
        self.assertIn("OperationError", json.dumps(handler))
        self.assertEqual(handler[0]["failedJob"]["switch"][0]["next"], "inspect")
        self.assertIn("rethrow", handler[1])

    def test_jobs_let_open_mpi_use_all_vcpus(self):
        rows = batch_dag.schedule(self.plan, self.profile, RUN, "a" * 64)
        for row in rows:
            env = row["spec"]["taskGroups"][0]["taskSpec"]["environment"]["variables"]
            self.assertEqual(env["OMPI_MCA_hwloc_base_use_hwthreads_as_cpus"], "1")

    def test_e2_family_quota_caps_parallel_vms(self):
        for i in range(5):
            self.plan["tasks"].append(dict(self.plan["tasks"][0], id="probe-" + str(i)))
        self.profile.update(max_parallel_vms=8, tasks_per_slot=1, regional_vcpu_limit=200,
                            regional_disk_limit_gb=4000)
        unlimited = max(r["max_instances"] for r in batch_dag.schedule(self.plan, self.profile, RUN, "a" * 64))
        self.profile["regional_e2_vcpu_limit"] = 4
        rows = batch_dag.schedule(self.plan, self.profile, RUN, "a" * 64)
        cores = {r["machine_type"]: next(c for n, c, _ in batch_dag.MACHINES if n == r["machine_type"]) for r in rows}
        self.assertTrue(all(r["max_instances"] * cores[r["machine_type"]] <= 4 for r in rows))
        self.assertGreater(unlimited, max(r["max_instances"] for r in rows))


class BudgetSettlement(unittest.TestCase):
    def ledger(self):
        return {"project": "p", "currency": "USD", "holds": {
            "build-1": {"upper_usd": "31.19"}, "cf-done": {"upper_usd": "37.37"},
            "cf-live": {"upper_usd": "61.44"}}}

    def test_settles_finished_runs_and_keeps_audit(self):
        ledger, record = budget.settle(self.ledger(), {"build-1": "3.00", "cf-done": "3"}, "owner",
                                       "console credit 298.71", finished={"build-1", "cf-done"},
                                       when="2026-11-27T00:00:00Z")
        self.assertEqual(ledger["holds"]["cf-done"]["upper_usd"], "3.00")
        self.assertEqual(ledger["holds"]["cf-done"]["original_upper_usd"], "37.37")
        self.assertTrue(ledger["holds"]["build-1"]["settled"])
        self.assertEqual(ledger["holds"]["cf-live"]["upper_usd"], "61.44")
        self.assertEqual(record["holds_total_after_usd"], "67.44")
        self.assertEqual(len(ledger["holds"]), 3)

    def test_refuses_live_unknown_repeat_and_inflated(self):
        with self.assertRaisesRegex(ValueError, "not confirmed finished"):
            budget.settle(self.ledger(), {"cf-live": "1"}, "owner", "e", finished=set())
        with self.assertRaisesRegex(ValueError, "no hold"):
            budget.settle(self.ledger(), {"cf-x": "1"}, "owner", "e", finished={"cf-x"})
        with self.assertRaisesRegex(ValueError, "not exceed"):
            budget.settle(self.ledger(), {"cf-done": "40"}, "owner", "e", finished={"cf-done"})
        ledger, _ = budget.settle(self.ledger(), {"cf-done": "1"}, "owner", "e", finished={"cf-done"})
        with self.assertRaisesRegex(ValueError, "already settled"):
            budget.settle(ledger, {"cf-done": "1"}, "owner", "e", finished={"cf-done"})
        with self.assertRaisesRegex(ValueError, "approving"):
            budget.settle(self.ledger(), {"cf-done": "1"}, " ", "e", finished={"cf-done"})


del BatchDagTests, ExtraAdapters  # do not re-run the imported suites from this module

if __name__ == "__main__":
    unittest.main()
