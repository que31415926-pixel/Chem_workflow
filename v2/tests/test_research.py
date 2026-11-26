import io
import tempfile
import unittest
from pathlib import Path

from chemflow2 import controller, research
from chemflow2.protocol import _digest, load_plan
from chemflow2.worker import run_task


class ResearchTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.protocol = self.root/"design.json"
        self.intake = self.root/"intake.json"
        research.write(self.protocol, {"schema": 1, "name": "review exercise", "steps": [
            {"id": "A", "deps": [], "claims": ["surface"], "groups": [
                {"id": "review", "kind": "review", "axes": {"site": ["x", "y"]},
                 "outputs": ["signed"]}]},
            {"id": "B", "deps": ["A"], "claims": ["surface"], "groups": [
                {"id": "check", "kind": "review", "outputs": ["signed"]}]}]})
        self.cells = {}
        for i, cell in enumerate(["A/review/x", "A/review/y", "B/check"]):
            folder = self.root/str(i)
            folder.mkdir()
            (folder/"record.txt").write_text(f"original {i}")
            form = folder/"form.json"
            research.write(form, {"scope": cell, "rationale": "primary record matches the assignment",
                                 "reviewer": "tester", "decision": "pass", "items": [
                                     {"file": "record.txt", "locator": "Table 1",
                                      "sha256": _digest(folder/"record.txt")}]})
            self.cells[cell] = {"source": f"{i}/form.json", "source_sha256": _digest(form),
                                "exports": {"signed": "evaluate/report.json"}}
        self.write_intake()

    def write_intake(self):
        research.write(self.intake, {"protocol_sha256": _digest(self.protocol), "cells": self.cells})

    def test_parallel_cells_with_same_local_filenames_and_dependency_receipts(self):
        dst = self.root/"compiled"
        research.compile_research(self.protocol, self.intake, dst)
        plan = load_plan(dst/"plan.json")
        self.assertEqual(len(plan["tasks"][2]["deps"]), 2)
        wire = io.BytesIO()
        controller.bundle(dst/"plan.json", wire)
        wire.seek(0)
        state = self.root/"state"
        controller.accept(state, "study1", wire)
        run = controller.run_path(state, "study1")
        for i, task in enumerate(plan["tasks"]):
            work = controller.prepare(run, task)
            self.assertEqual((work/"record.txt").read_text(), f"original {i}")
            result = run_task(task, work)
            self.assertTrue(result["ok"], result)
            self.assertEqual(research.read(work/"report.json")["scientific_gate"], "pass")

    def test_missing_cell_never_compiles_as_whole_study(self):
        del self.cells["B/check"]
        self.write_intake()
        with self.assertRaisesRegex(ValueError, "missing"):
            research.compile_research(self.protocol, self.intake, self.root/"all")
        report = research.compile_research(self.protocol, self.intake, self.root/"first", through="A")
        self.assertEqual(report["cells"], 2)

    def test_scope_and_review_evidence(self):
        self.cells["A/review/y"] = self.cells["A/review/x"]
        self.write_intake()
        with self.assertRaisesRegex(ValueError, "explicitly cover"):
            research.compile_research(self.protocol, self.intake, self.root/"compiled")
        form = self.root/"0"/"form.json"
        data = research.read(form)
        (form.parent/"record.txt").write_text("changed")
        with self.assertRaisesRegex(ValueError, "digest mismatch"):
            research.run_review({"input": "form.json", "outputs": ["out.json"]}, form.parent)

    def test_frozen_quantum_matrix_uses_approved_xyz(self):
        xyz = self.root / "h2.xyz"
        xyz.write_text("2\nfixture\nH 0 0 0\nH 0 0 0.7\n")
        registry = self.root / "registry.json"
        research.write(registry, {"h2": {"states": [{"id": "singlet", "charge": 0, "multiplicity": 1}],
                         "conformers": [{"id": "one", "xyz": xyz.name, "sha256": _digest(xyz),
                                         "source": "fixture", "approved_by": "tester"}]}})
        (self.root / "geom.tpl").write_text("! HF Opt CPCM(water)\n%pal nprocs 1 end\n"
            "%maxcore 100\n* xyzfile {{charge}} {{multiplicity}} {{xyz}}\n")
        research.write(self.root / "methods.json", [{"id": "geom", "stage": "opt", "template": "geom.tpl",
                                    "resources": {"cores": 1, "memory_mb": 1000, "scratch_mb": 200}}])
        research.write(self.protocol, {"schema": 1, "name": "matrix", "steps": [{"id": "A", "deps": [],
                      "claims": ["tested"], "groups": [{"id": "qm", "kind": "study", "outputs": ["structure"],
                      "axes": {"species": ["h2"], "model": ["CPCM"]},
                      "method_sets": {"CPCM": "methods.json"}}]}]})
        prepared = self.root / "matrix"
        research.prepare_quantum_matrix(self.protocol, registry, prepared)
        study = prepared / "A/qm/h2/CPCM/study.json"
        res = research.compile_study(study, self.root / "one")
        self.assertEqual(res["tasks"], 1)
        intake = prepared / "intake.json"
        self.assertEqual(research.assemble(self.protocol, prepared, intake)["registered"], 1)
        report = research.compile_research(self.protocol, intake, self.root / "full")
        self.assertEqual(report["tasks"], 2)
        plan = load_plan(self.root / "full/plan.json")
        self.assertEqual(plan["tasks"][-1]["outputs"], ["matrix.json"])
        self.assertEqual(list(plan["tasks"][-1]["from_deps"].values()),
                         [plan["tasks"][0]["id"] + "/values.json"])
        xyz.write_text("2\nchanged\nH 0 0 0\nH 0 0 0.8\n")
        with self.assertRaisesRegex(ValueError, "changed"):
            research.prepare_quantum_matrix(self.protocol, registry, self.root / "other")

    def test_isotopic_observables(self):
        from chemflow2 import mechanism
        row = {"condition_id": "matched", "initial_rate_unit": "M/s", "isotopic_rates": {"H": 3, "D": 1},
               "pressure_series": [{"hydrogen_pressure_bar": 1, "initial_rate": 2},
                                   {"hydrogen_pressure_bar": 2, "initial_rate": 1}],
               "site_deuterium_mol": {"A": .3, "B": .7}, "deuterium_input_mol": 1}
        data = {"operation": "isotope_observables", "explanations": [dict(row,id=n) for n in ("a","b")]}
        result = mechanism.evaluate(data)["explanations"]["a"]
        self.assertAlmostEqual(result["kinetic_isotope_effect"], 3)
        self.assertAlmostEqual(result["hydrogen_orders"][0]["order"], -1)
        data["explanations"][0]["site_deuterium_mol"] = {"A": .3}
        with self.assertRaisesRegex(ValueError, "inventory"):
            mechanism.evaluate(data)
