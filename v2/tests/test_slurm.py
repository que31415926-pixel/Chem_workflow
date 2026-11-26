import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from chemflow2.slurm import collect, export, run_stage, submit


class SlurmTests(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.root = Path(self.scratch.name)
        self.bundle = self.root / "portable bundle"
        self.site = {"name": "test", "host": "slurm", "cores": 1, "memory_mb": 50,
                     "scratch_mb": 1, "engines": ["probe"], "capabilities": [], "wall_minutes": 5}
        self.plan = self.root / "source.json"
        self.tasks = [
            {"id": "down", "engine": "probe", "stage": "probe", "deps": ["up"],
             "from_deps": {"observed.json": "up/initial.json"}, "outputs": ["final.json"]},
            {"id": "up", "engine": "probe", "stage": "probe", "files": ["seed.dat"],
             "outputs": ["initial.json"]},
        ]
        for task in self.tasks:
            task["resources"] = {"cores": 1, "memory_mb": 40, "scratch_mb": 1}
        self.plan.write_text(json.dumps({"name": "offline smoke", "tasks": self.tasks}))
        (self.root / "seed.dat").write_text("portable input")
        export(self.plan, self.bundle)

    def fake_submit(self):
        with patch("chemflow2.slurm.subprocess.run", side_effect=[
            subprocess.CompletedProcess([], 0, "1100\n"), subprocess.CompletedProcess([], 0, "1101;local\n")
        ]) as sbatch:
            state = submit(self.bundle, self.site, partition="cpu")
        first, second = (call.args[0] for call in sbatch.call_args_list)
        self.assertNotIn("--dependency=afterok:1100", first)
        self.assertIn("--dependency=afterok:1100", second)
        self.assertIn("--mem=40M", first)
        self.assertIn("--ntasks=1", first)
        return state["run"]

    def test_stages_on_new_compute_processes_and_verified_retrieval(self):
        run_id = self.fake_submit()
        env = dict(os.environ, SLURM_NTASKS="1", SLURM_MEM_PER_NODE="40", SLURM_TMPDIR=str(self.root))
        for task in ("up", "down"):
            result = subprocess.run(["sh", str(self.bundle / "stage.sh"), str(self.bundle), run_id, task],
                                    cwd=self.root, env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr +
                             (self.bundle / "runs" / run_id / "tasks" / task / "chemflow-result.json").read_text())
        report = collect(self.bundle, run_id)
        self.assertTrue(report["all_tasks_succeeded"])
        parent = self.bundle / "runs" / run_id / "tasks/up/initial.json"
        parent.write_text("changed")
        self.assertEqual(collect(self.bundle, run_id)["tasks"]["up"], "invalid_receipt")

    def test_attached_frozen_evidence_plan_survives_export(self):
        name = "frozen-claims.json"
        (self.root / name).write_text('{"locked":true}')
        plan = json.loads(self.plan.read_text())
        plan["attachments"] = [name]
        self.plan.write_text(json.dumps(plan))
        target = self.root / "with-attachments"
        export(self.plan, target)
        self.assertEqual((target / "inputs" / name).read_text(), '{"locked":true}')
        self.assertEqual(submit(target, self.site, check=True)["status"], "admissible")
        (target / "inputs" / name).write_text('{"locked":false}')
        with self.assertRaisesRegex(ValueError, "Changed"):
            submit(target, self.site, check=True)

    def test_rejects_resources_before_submission(self):
        self.site["memory_mb"] = 35
        with patch("chemflow2.slurm.subprocess.run") as sbatch, self.assertRaisesRegex(ValueError, "cannot admit"):
            submit(self.bundle, self.site)
        sbatch.assert_not_called()

    def test_tampered_input_blocks_submission(self):
        (self.bundle / "inputs/seed.dat").write_text("modified")
        with self.assertRaisesRegex(ValueError, "Changed"):
            submit(self.bundle, self.site, check=True)

    def test_wrong_runtime_allocation(self):
        run_id = self.fake_submit()
        with patch.dict(os.environ, {"SLURM_NTASKS": "1", "SLURM_MEM_PER_NODE": "20"}):
            result = run_stage(self.bundle, run_id, "up")
        self.assertFalse(result["ok"])
        self.assertIn("Allocated memory", result["error"])
        self.assertEqual(collect(self.bundle, run_id)["tasks"], {"up": "failed", "down": "no_receipt"})

    def test_parent_hash_checked_before_using_output(self):
        run_id = self.fake_submit()
        with patch.dict(os.environ, {"SLURM_NTASKS": "1", "SLURM_MEM_PER_NODE": "50", "TMPDIR": str(self.root)}):
            self.assertTrue(run_stage(self.bundle, run_id, "up")["ok"])
            (self.bundle / "runs" / run_id / "tasks/up/initial.json").write_text("modified")
            result = run_stage(self.bundle, run_id, "down")
        self.assertFalse(result["ok"])
        self.assertIn("SHA256 mismatch", result["error"])

    def test_partial_submission_records_accepted_job_id(self):
        with patch("chemflow2.slurm.subprocess.run", side_effect=[
            subprocess.CompletedProcess([], 0, "1150"), subprocess.CalledProcessError(1, "sbatch", stderr="limit")
        ]), self.assertRaisesRegex(RuntimeError, "Incomplete submission"):
            submit(self.bundle, self.site)
        receipt = json.loads(next((self.bundle / "runs").glob("*/submission.json")).read_text())
        self.assertEqual(receipt["jobs"], {"up": "1150"})
        self.assertEqual(receipt["state"], "submission_incomplete")


if __name__ == "__main__":
    unittest.main()
