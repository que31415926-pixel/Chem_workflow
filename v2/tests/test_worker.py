"""Run only inert stand-ins for the scientific executable."""
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from chemflow2.worker import run_task, sha256


class WorkerTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.work = self.root / "attempt"
        self.work.mkdir()
        self.binary = self.root / "fixture"
        mock_env = patch.dict(os.environ, {"CHEMFLOW_ORCA": str(self.binary)})
        mock_env.start()
        self.addCleanup(mock_env.stop)

    def fixture(self, log, files=None, code=0):
        self.binary.write_text("#!/usr/bin/env python3\nfrom pathlib import Path\n"
                               + "print(" + repr(log) + ")\n"
                               + "".join("Path(" + repr(name) + ").write_text(" + repr(value) + ")\n"
                                         for name, value in (files or {}).items())
                               + f"raise SystemExit({code})\n")
        self.binary.chmod(0o700)

    def task(self, stage, outputs=None, deck=None):
        inp = self.work / "job.inp"
        if deck is None:
            keyword = {"opt": "Opt", "freq": "Freq", "sp": "SP", "ts": "OptTS", "irc": "IRC"}[stage]
            deck = "! HF def2-SVP " + keyword + "\n%maxcore 128\n* xyz 0 1\nH 0 0 0\nH 0 0 1\n*\n"
        inp.write_text(deck)
        return dict(id="fixture", stage=stage, engine="orca", input=inp.name,
                    outputs=outputs or ["orca.out"], resources=dict(cores=1, memory_mb=2048, scratch_mb=0))

    def test_probe_records_dependency_hash_and_never_executes_a_command(self):
        (self.work / "upstream.xyz").write_text("example payload")
        task = dict(id="probe", stage="probe", engine="command", outputs=["finished.json"])
        result = run_task(task, self.work)
        self.assertTrue(result["ok"], result)
        payload = json.loads((self.work / "finished.json").read_text())
        self.assertEqual(payload["inputs"]["upstream.xyz"], sha256(self.work / "upstream.xyz"))
        self.assertEqual(result["artifacts"]["finished.json"], sha256(self.work / "finished.json"))
        self.assertFalse(run_task(dict(task, outputs=["malicious"], argv=["not-allowed"]), self.work)["ok"])

    def test_stale_and_escaping_paths_are_never_published(self):
        (self.work / "input.xyz").write_text("input")
        for output in ("../escaped", "input.xyz", "chemflow-task.json"):
            res = run_task(dict(id="p", stage="probe", engine="probe", outputs=[output]), self.work)
            self.assertEqual(res["failure_kind"], "input")
            self.assertFalse(res["artifacts"])
        self.assertFalse((self.root / "escaped").exists())
        (self.work / "outside").symlink_to(self.root)
        res = run_task(dict(id="p", stage="probe", engine="probe", outputs=["ok"]), self.work)
        self.assertEqual(res["failure_kind"], "input")

    def test_opt_and_sp(self):
        self.fixture("Program Version 6.1.1\nTHE OPTIMIZATION HAS CONVERGED\n"
                     "FINAL SINGLE POINT ENERGY -1.1\n****ORCA TERMINATED NORMALLY****",
                     {"job.xyz": "2\nH2\nH 0 0 0\nH 0 0 0.7\n"})
        task = self.task("opt", ["job.xyz"])
        res = run_task(task, self.work)
        self.assertTrue(res["ok"], res)
        self.assertEqual(set(res["artifacts"]), {"job.xyz", "orca.out"})
        self.assertEqual(res["observations"]["n_atoms_geometry"], 2)
        metrics = res["metrics"]
        self.assertEqual(metrics["requested_core_seconds"], metrics["wall_seconds"])
        self.assertGreaterEqual(metrics["child_system_cpu_seconds"], 0)
        self.assertGreaterEqual(metrics["scratch_regular_file_bytes_at_finish"], metrics["published_artifact_bytes"])
        self.assertIsNone(metrics["actual_cost"])
        self.assertIsNone(metrics["scratch_peak_bytes"])
        self.assertEqual(res["observations"]["engine_version"], "6.1.1")
        self.assertEqual(res["observations"]["scientific_gate"], "not_evaluated")
        self.assertEqual(res["provenance"]["input_sha256"]["job.inp"], sha256(self.work / "job.inp"))

    def test_sp_checks_version_and_finite_energy(self):
        self.fixture("Program Version 6.1.1\nFINAL SINGLE POINT ENERGY -1.1D+00\n****ORCA TERMINATED NORMALLY****")
        task = self.task("sp")
        task["engine_version"] = "6.1.1"
        self.assertTrue(run_task(task, self.work)["ok"])
        (self.work / "orca.out").unlink()
        task["engine_version"] = "6.0"
        self.assertEqual(run_task(task, self.work)["failure_kind"], "capability")

    def test_converged_opt_log_cannot_qualify_an_incomplete_frequency(self):
        self.fixture("THE OPTIMIZATION HAS CONVERGED\nVIBRATIONAL FREQUENCIES\n0: 0.00 cm**-1\n"
                     "****ORCA TERMINATED NORMALLY****")
        res = run_task(self.task("freq"), self.work)
        self.assertEqual(res["failure_kind"], "validation")
        self.assertEqual(set(res["artifacts"]), {"orca.out"})

    def test_spectrum_reports_signed_modes_without_preempting_scientific_gate(self):
        lines = ["Number of atoms ... 2", "VIBRATIONAL FREQUENCIES"]
        lines.extend(f"{ix}: {value} cm**-1" for ix, value in enumerate([0, 0, 0, 0, 0, -35.1]))
        lines += ["NORMAL MODES", "Final Gibbs free energy ... -1.07 Eh", "****ORCA TERMINATED NORMALLY****"]
        self.fixture("\n".join(lines))
        res = run_task(self.task("freq"), self.work)
        self.assertTrue(res["ok"], res)
        self.assertEqual(res["observations"]["negative_frequencies_cm-1"], [-35.1])
        self.assertEqual(res["observations"]["gibbs_free_energy_eh"], -1.07)

    def test_irc_requires_endpoints(self):
        self.fixture("IRC PATH SUMMARY\n****ORCA TERMINATED NORMALLY****")
        res = run_task(self.task("irc"), self.work)
        self.assertEqual(res["failure_kind"], "validation")

    def test_normal_exit_alone_and_input_echo_cannot_pass_sp(self):
        self.fixture("|  1> # FINAL SINGLE POINT ENERGY -1.2\n****ORCA TERMINATED NORMALLY****")
        res = run_task(self.task("sp"), self.work)
        self.assertEqual(res["failure_kind"], "validation")
        self.assertFalse(res["ok"])

    def test_stage_mismatch_and_memory_overcommit_fail_before_launch(self):
        task = self.task("freq", deck="! HF Opt\n")
        self.assertEqual(run_task(task, self.work)["failure_kind"], "input")
        task = self.task("sp", deck="! HF PAL4\n%maxcore 1000\n")
        task["resources"].update(cores=4, memory_mb=4096)
        self.assertEqual(run_task(task, self.work)["failure_kind"], "resources")
        self.assertFalse((self.work / "orca.out").exists())

    def test_failed_executable_records_diagnostic_log_without_outputs(self):
        self.fixture("XC-Hessian failed", {"job.xyz": "partial"}, code=11)
        res = run_task(self.task("opt", ["job.xyz"]), self.work)
        self.assertEqual(res["failure_kind"], "engine")
        self.assertEqual(res["provenance"]["exit_code"], 11)
        self.assertEqual(set(res["artifacts"]), {"orca.out"})
