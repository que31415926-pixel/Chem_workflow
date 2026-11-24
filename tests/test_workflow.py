import csv
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from chemflow.analysis import fit, observations
from chemflow.material import account
from chemflow.network import ROOT, Network
from chemflow.qm import geometry, orca_audit, prepare
from chemflow.thermo import reaction_energy
from chemflow.wet import handoff, paired


def args(project, species, mult=None, xyz=None):
    return SimpleNamespace(
        project=project,
        species=species,
        xyz=xyz,
        frame=1,
        mult=mult,
        cores=1,
        memory=1500,
        stage="min",
        method="screen",
        solvation="default",
    )


def fixture_output(g="-300.001", imag=None):
    mode = imag or 180.0
    return (
        f"THE OPTIMIZATION HAS CONVERGED\nFINAL SINGLE POINT ENERGY -300.102300\n"
        f"VIBRATIONAL FREQUENCIES\n  1:  {mode:.2f} cm**-1\nNORMAL MODES\n"
        f"Temperature ... 298.15 K\nPressure ... 1.00 atm\n"
        f"Final Gibbs free energy ... {g} Eh\nORCA TERMINATED NORMALLY\n"
    )


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.d = Path(self.temp.name)

    def test_conserved_mass_and_no_metal(self):
        for project in ("fe", "co"):
            m = Network(project)
            for h in m.hypotheses:
                for c in m.conditions:
                    y = m.simulate(h, c, m.data["times_h"])
                    totals = m.balances @ y
                    np.testing.assert_allclose(
                        totals, np.repeat(totals[:, :1], y.shape[1], axis=1), atol=2e-8
                    )
                    self.assertGreaterEqual(y.min(), 0)
            if project == "fe":
                y = m.simulate("slow_SH2", "no_Fe", [24.0])
                self.assertEqual(y[m.names.index("FeBz"), 0], 0.0)
                self.assertGreater(y[m.names.index("P"), 0], 0.0)
                self.assertEqual(
                    m.simulate("slow_SH2", "dark", [24.0])[m.names.index("P"), 0], 0.0
                )
            else:
                y = m.simulate("slow_CH", "blank", [48.0])
                self.assertEqual(y[m.names.index("P"), 0], 0.0)
                self.assertGreater(y[m.names.index("Hgas"), 0], 0.0)

    def test_all_registered_seed_spins(self):
        for p in ("co", "fe"):
            reg = json.loads(
                (ROOT / "projects" / p / "species.json").read_text(encoding="utf-8")
            )
            for name, rec in reg.items():
                if "composition" in rec:
                    continue
                for m in rec["multiplicities"]:
                    prepare(args(p, name, m), self.d / (name + str(m)))
        bad = args("fe", "fe_sec_bn", 3)
        with self.assertRaises(ValueError):
            prepare(bad, self.d / "bad")

    def test_ensemble_selection(self):
        path = ROOT / "projects/co/examples/veratrole_crest.xyz"
        a = geometry("co", "veratrole", path, 1)
        b = geometry("co", "veratrole", path, 2)
        self.assertEqual(a[1], b[1])
        self.assertNotEqual(a[2], b[2])
        with self.assertRaises(ValueError):
            geometry("fe", "pivalic_acid", path, 1)

    def test_fit_and_no_leakage(self):
        for p in ("co", "fe"):
            net = Network(p)
            source = ROOT / "projects" / p / "examples/synthetic_observations.csv"
            path = self.d / p
            fit(net, source, path)
            original = "slow_CH" if p == "co" else "poor_supply"
            restored = Network(p, path / "calibrated.json")
            key = net.hypotheses[original]["fit"]
            self.assertAlmostEqual(
                restored.settings(original, "base")[1][key],
                net.settings(original, "base")[1][key] * 1.7,
                places=4,
            )
            with source.open() as f:
                rows = list(csv.DictReader(f))
            rows[-1]["condition"] = "base"
            broken = self.d / (p + ".csv")
            with broken.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaises(ValueError):
                observations(broken, net)

    def test_gates_and_censored_data(self):
        path = ROOT / "projects/fe/state.json"
        state = json.loads(path.read_text(encoding="utf-8"))
        state.update(controls="pass", salt="modified_only")
        user = self.d / "state.json"
        user.write_text(json.dumps(state))
        result = handoff("fe", user, self.d / "handoff")
        self.assertTrue(all("原三苯环" in r["comparison"] for r in result))
        self.assertEqual(result[0]["phase"], "F1")
        result = paired(ROOT / "projects/fe/examples/synthetic_paired.csv")
        self.assertEqual(result["target_state"], "low")
        self.assertGreater(
            result["target"]["mean_bounds"][1], result["target"]["mean_bounds"][0]
        )

    def test_material_missing_vs_zero(self):
        source = ROOT / "projects/co/examples/synthetic_material.csv"
        with source.open() as f:
            row = next(csv.DictReader(f))
        row["b2_left"] = ""
        path = self.d / "measured.csv"
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=row.keys())
            w.writeheader()
            w.writerow(row)
        account("co", path, self.d / "balance")
        with (self.d / "balance" / "material_balance.csv").open() as f:
            result = next(csv.DictReader(f))
        self.assertEqual(result["assigned_B_pct"], "")
        self.assertAlmostEqual(float(result["TON_total_Co"]), 14.0)

    def test_audit_and_balanced_thermochemistry(self):
        for name in ("fe_sec_bn", "fe_oep", "sec_benzyl_rad"):
            out = self.d / name
            prepare(args("fe", name), out)
            (out / "orca.out").write_text(fixture_output())
        out = self.d / "fe_sec_bn"
        self.assertEqual(orca_audit(out)["status"], "numerical_checks_passed")
        recipe = self.d / "rx.json"
        recipe.write_text(
            json.dumps(
                {
                    "reactants": {"fe_sec_bn": 1},
                    "products": {"fe_oep": 1, "sec_benzyl_rad": 1},
                }
            )
        )
        self.assertEqual(reaction_energy(recipe)["delta_n"], 1)
        recipe.write_text(
            json.dumps({"reactants": {"fe_sec_bn": 1}, "products": {"fe_oep": 1}})
        )
        with self.assertRaises(ValueError):
            reaction_energy(recipe)
        (out / "orca.out").write_text(fixture_output(imag=-182))
        self.assertEqual(orca_audit(out)["status"], "incomplete_or_review")
        (out / "orca.out").write_text("SCF CONVERGED")
        self.assertEqual(orca_audit(out)["status"], "incomplete_or_review")


if __name__ == "__main__":
    unittest.main()
