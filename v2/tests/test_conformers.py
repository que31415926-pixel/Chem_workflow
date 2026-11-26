import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from chemflow2.conformers import inventory, prepare
from chemflow2.study import compile_study
from chemflow2.worker import Rejected


def sha(file):
    return hashlib.sha256(file.read_bytes()).hexdigest()


class ConformerSelectionTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.src = self.root / "input"
        shutil.copytree(Path(__file__).resolve().parents[1] / "examples/water-calibration", self.src)
        study = json.loads((self.src / "study.json").read_bytes())
        study["species"][0]["conformers"] = []
        self.tpl = self.src / "study-template.json"
        self.tpl.write_text(json.dumps(study))
        frame1 = (self.src / "water.xyz").read_text()
        frame2 = frame1.replace("0.957200", "0.965000")
        self.xyz = self.src / "ensemble.xyz"
        self.xyz.write_text(frame1 + frame2)
        self.eng = self.src / "crest.energies"
        self.eng.write_text("1 0.000\n2 0.230\n")
        self.out = self.root / "prepared"
        self.file = self.src / "review.json"
        self.inv = inventory(self.xyz, self.eng)
        self.review = {"species_id": "water", "reviewer": "Example reviewer", "rationale": "Fixture review",
                       "ensemble_sha256": sha(self.xyz), "energies_sha256": sha(self.eng),
                       "study_template": self.tpl.name, "study_template_sha256": sha(self.tpl),
                       "frames": [{"frame": r["frame"], "xyz_sha256": r["xyz_sha256"],
                                   "decision": "exclude", "reason": "Reviewed fixture duplicate"}
                                  for r in self.inv["frames"]]}
        self.review["frames"][0].update(id="selected", decision="include", reason="Fixture reference",
                                        degeneracy=1, rotational_symmetry_number=2)

    def execute(self):
        self.file.write_text(json.dumps(self.review))
        return prepare(self.xyz, self.eng, self.file, self.out)

    def rejects(self):
        with self.assertRaises((ValueError, Rejected, OSError)):
            self.execute()
        self.assertFalse(self.out.exists())

    def test_complete_review_preserves_all_frames_and_compiles(self):
        self.assertEqual(self.inv["ensemble_sha256"], sha(self.xyz))
        self.assertEqual(self.inv["energies_sha256"], sha(self.eng))
        self.assertEqual(self.inv["energy_unit"], "kcal/mol")
        self.assertEqual([r["frame"] for r in self.inv["frames"]], [1, 2])
        self.assertNotEqual(self.inv["frames"][0]["xyz_sha256"], self.inv["frames"][1]["xyz_sha256"])
        self.execute()
        cfg = json.loads((self.out / "study.json").read_bytes())
        included = cfg["species"][0]["conformers"]
        self.assertEqual([r["id"] for r in included], ["selected"])
        self.assertEqual(sha(self.out / included[0]["xyz"]), self.inv["frames"][0]["xyz_sha256"])
        ledger = json.loads((self.out / "audit/ledger.json").read_bytes())
        self.assertEqual([r["decision"] for r in ledger["frames"]], ["include", "exclude"])
        for row in ledger["frames"]:
            self.assertEqual(sha(self.out / row["xyz"]), self.inv["frames"][row["frame"]-1]["xyz_sha256"])
        result = compile_study(self.out / "study.json", self.root / "compile")
        self.assertEqual(result["matrix_rows"], 4)
        with self.assertRaises((ValueError, OSError)):
            self.execute()

    def test_incomplete_decision_and_zero_selected_fail_atomically(self):
        excluded = self.review["frames"].pop()
        self.rejects()
        self.review["frames"] = [dict(self.review["frames"][0], decision="exclude"), excluded]
        self.rejects()

    def test_unmatched_digests_invalid_counts_and_unsafe_template(self):
        for row, key, invalid in ((self.review, "ensemble_sha256", "f"*64),
                                  (self.review, "energies_sha256", "f"*64),
                                  (self.review, "study_template_sha256", "f"*64),
                                  (self.review["frames"][1], "xyz_sha256", "f"*64),
                                  (self.review["frames"][0], "degeneracy", 0),
                                  (self.review["frames"][0], "rotational_symmetry_number", True),
                                  (self.review, "study_template", "../escaped.json")):
            with self.subTest(field=key):
                before = row[key]
                row[key] = invalid
                self.rejects()
                row[key] = before

    def test_missing_energy_frame_rejected(self):
        self.eng.write_text("1 0.0\n")
        with self.assertRaises((ValueError, Rejected)):
            inventory(self.xyz, self.eng)
