"""New-project acceptance tests; no project names need to be registered in the engine."""

import csv
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from chemflow.network import Network
from chemflow.projects import ROOT, discover, read_json, write_json
from chemflow.runner import run_project, validate
from chemflow.wet import handoff


class ReuseTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.registry = self.root / "presets"
        self.name = "isomer_lab"
        self.project = self.registry / self.name
        shutil.copytree(ROOT / "templates" / "starter", self.project)
        env = patch.dict(os.environ, {"CHEMFLOW_PROJECTS_DIR": str(self.registry)})
        env.start()
        self.addCleanup(env.stop)

    def test_linked_nonmetal_project(self):
        cfg = read_json(self.project / "project.json")
        cfg["qm"]["temperature_K"] = 312.0
        cfg["qm"]["methods"] = {"custom": "PBE0 def2-SVP"}
        cfg["qm"]["jobs"][0]["method"] = "custom"
        write_json(self.project / "project.json", cfg)
        self.assertEqual(discover(), [self.name])
        validate(self.name)
        out = self.root / "out"
        run_project(self.name, out)
        self.assertIn("312.00", (out / "orca_input/job.inp").read_text())
        self.assertIn("PBE0 def2-SVP", (out / "orca_input/job.inp").read_text())
        self.assertEqual(
            read_json(out / "forecast/provenance.json")["origin"], ["synthetic_fixture"]
        )
        with (out / "closure/material_balance.csv").open() as f:
            self.assertAlmostEqual(
                float(next(csv.DictReader(f))["alkene_recovery_pct"]), 99.0
            )
        rows = handoff(
            self.name, self.project / "examples/state_ready.json", self.root / "updated"
        )
        self.assertEqual(rows[0]["phase"], "paired_start")

    def test_small_workflow_and_invalid_dependencies(self):
        cfg = {"schema_version": 1, "workflow": ["handoff"]}
        write_json(self.project / "project.json", cfg)
        for fn in ["network.json", "material.json", "species.json"]:
            (self.project / fn).unlink()
        validate(self.name)
        out = self.root / "wet_only"
        run_project(self.name, out)
        self.assertEqual([p.name for p in out.iterdir()], ["wet"])
        for bad in ["@later/observations.csv", "../outside.csv"]:
            cfg["workflow"] = [{"step": "fit", "data": bad}]
            write_json(self.project / "project.json", cfg)
            with self.assertRaises(ValueError):
                validate(self.name)

    def test_one_scenario_and_multiple_rates(self):
        p = self.project / "network.json"
        net = read_json(p)
        del net["hypotheses"]["fast_reverse"]
        net["hypotheses"]["slow_forward"]["fit"] = ["forward", "reverse"]
        net.pop("panels")
        write_json(p, net)
        out = self.root / "fit_result"
        run_project(self.name, out)
        fitted = read_json(out / "estimate/calibrated.json")["parameters"][
            "slow_forward"
        ]
        self.assertAlmostEqual(fitted["forward"], 0.09, places=5)
        self.assertAlmostEqual(fitted["reverse"], 0.08, places=5)
        with (out / "forecast/diagnostics.csv").open() as f:
            self.assertEqual(list(csv.DictReader(f)), [])
        net["conditions"]["reference"]["initial_M"] = {"E": 0.2}
        write_json(p, net)
        with self.assertRaisesRegex(ValueError, "different network"):
            Network(self.name, out / "estimate/calibrated.json")

    def test_cli_scaffold_discovery_and_legacy_alias(self):
        def cli(*argv):
            return subprocess.run(
                [sys.executable, "-m", "chemflow", *argv],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )

        cli("init", "another_test")
        self.assertIn("another_test", cli("projects").stdout)
        out = self.root / "cli"
        cli("run", "another_test", "--out", str(out))
        self.assertTrue((out / "forecast/trajectories_M.csv").is_file())
        args = ("--projects-dir", str(ROOT / "projects"))
        data = str(ROOT / "projects/fe/examples/synthetic_paired.csv")
        self.assertEqual(
            cli(*args, "paired", "fe", data).stdout,
            cli(*args, "paired-fe", data).stdout,
        )


if __name__ == "__main__":
    unittest.main()
