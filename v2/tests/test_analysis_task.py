import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from chemflow2.analysis_task import _disposition, run_analysis


class AnalysisTaskTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.folder = Path(tmp.name)
        self.task = {"engine": "analysis", "stage": "science", "input": "input.json",
                     "outputs": ["report.json"], "from_deps": {"up.json": "upstream/parsed.json"}}
        self.primary = self.folder / "input.json"
        self.ref = {"$artifact": "up.json", "path": ["list", 0]}
        self.data = {"operation": "gate", "protocol_id": "p", "gate_id": "example-energy",
                     "gate_spec": {"si_difference": ["max_abs","kcal/mol",0.5],
                                   "electronic_difference_eh": ["max_abs","Eh",None]},
                     "gate_method": {"method_id": "m", "electronic_tolerance_eh": 0.000001},
                     "required_pairs": {"si_difference": ["one"], "electronic_difference_eh": ["state"]},
                     "metrics": {"electronic_difference_eh": {"unit": "Eh", "pairs": {"state": {
                         "left": 0.0, "right": 0.0, "left_source_id": "s", "right_source_id": "t"}}},
                         "si_difference": {"unit": "kcal/mol", "pairs": {
                         "one": {"left": self.ref, "right": 0.0, "left_source_id": "x", "right_source_id": "y"}}}}}

    def run_case(self, upstream):
        self.primary.write_text(json.dumps(self.data), encoding="utf-8")
        (self.folder / "up.json").write_text(json.dumps(upstream), encoding="utf-8")
        return run_analysis(self.task, self.folder)

    def test_inherited_binding_source_and_failed_report(self):
        info, failure = self.run_case({"list": [1]})
        self.assertEqual(failure, "scientific_fail")
        report = json.loads((self.folder / "report.json").read_text())
        self.assertEqual(report["status"], "fail")
        raw = (self.folder / "up.json").read_bytes()
        self.assertEqual(info["dependency_artifacts"]["up.json"],
                         {"source": "upstream/parsed.json", "sha256": hashlib.sha256(raw).hexdigest()})
        self.assertEqual(info["resolved_input_sha256"], report["analysis_receipt"]["resolved_input_sha256"])
        with self.assertRaises(ValueError):
            self.run_case({"list": [0.2]})

    def test_pass_and_missing(self):
        info, failure = self.run_case({"list": [0.5]})
        self.assertIsNone(failure)
        self.assertEqual(info["scientific_gate"], "pass")
        (self.folder / "report.json").unlink()
        self.data["metrics"] = {}
        _, failure = self.run_case({})
        self.assertEqual(failure, "scientific_missing")
        self.assertTrue((self.folder / "report.json").is_file())

    def test_reject_undeclared_broken_or_overloaded_binding(self):
        for ref in ({"$artifact": "../outside.json", "path": []}, {"$artifact": "input.json", "path": []},
                    {"$artifact": "up.json", "path": ["list", -1]},
                    {"$artifact": "up.json", "path": ["list", True]},
                    {"$artifact": "up.json", "path": ["list", 10]},
                    {"$artifact": "up.json", "path": ["list", 0], "default": 0}):
            with self.subTest(ref=ref), self.assertRaises(ValueError):
                self.ref.clear()
                self.ref.update(ref)
                self.run_case({"list": [1]})
        self.assertFalse((self.folder / "report.json").exists())

    def test_distinguish_method_evidence_and_prediction(self):
        for statuses, expected in [(["agreement"], "pass"), (["disagreement"], "fail"),
                                  (["within_band"], "indeterminate")]:
            report = {"operation": "paired_barriers", "paths": {"a": {"status": "qualified"}},
                      "comparisons": {str(i): {"status": x} for i, x in enumerate(statuses)}}
            self.assertEqual(_disposition("mechanism", report), expected)
        for path, expected in [("missing", "missing"), ("unqualified", "indeterminate")]:
            self.assertEqual(_disposition("mechanism", {"operation": "paired_barriers",
                             "paths": {"a": {"status": path}},
                             "comparisons": {"c": {"status": "insufficient_evidence"}}}), expected)
        for status, expected in [("inside", "pass"), ("outside", "fail"),
                                 ("awaiting_observation", "missing"), ("overlap", "indeterminate")]:
            self.assertEqual(_disposition("mechanism", {"operation": "score_holdout",
                             "targets": {"a": {"status": status}}}), expected)
        self.assertEqual(_disposition("mechanism", {"operation": "freeze_prediction"}), "not_applicable")


if __name__ == "__main__":
    unittest.main()
