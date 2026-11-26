import copy
import unittest

from chemflow2.analysis_task import _disposition
from chemflow2.mechanism import _hash, evaluate


class AnchorThresholdTests(unittest.TestCase):
    def setUp(self):
        anchors = {}
        training, heldout = {}, {}
        self.reg = {"score_id": "reviewed scalar", "unit": "kcal/mol", "score_protocol_sha256": "f"*64,
                    "higher_scores": "stable", "threshold_rule": "midpoint_of_strict_class_gap",
                    "limitations": "Synthetic contract fixture with qualitative labels", "anchors": anchors}
        for i, (role, label, score) in enumerate([("train","unstable",1), ("train","stable",3),
                                                ("heldout","stable",4), ("heldout","unstable",0)]):
            key = str(i)
            ident = {"condition_id":"C" + key, "group_id":"g" + key,
                     "source_sha256":"a"*64, "source_locator":"fixture, record " + key}
            observation = {**ident, "id":key, "label":label}
            anchors[key] = {**ident,"observation_sha256":_hash(observation), "split":role}
            (training if role == "train" else heldout)[key] = {
                "observation":observation,
                "score":{**{k: self.reg[k] for k in ("score_id", "score_protocol_sha256", "unit")},
                         "condition_id":ident["condition_id"], "group_id":ident["group_id"],
                         "artifact_sha256":key*64, "lower":score, "upper":score}}
        self.fit = {"operation":"fit_anchor_threshold","registration":self.reg,
                    "registration_sha256":_hash(self.reg),"training":training}
        self.score = {"operation":"score_anchor_threshold","registration":self.reg,
                      "registration_sha256":_hash(self.reg),"heldout":heldout}

    def score_report(self):
        return evaluate({**self.score,"calibration":evaluate(self.fit)})

    def test_training_only_threshold_and_correct_failure_routing(self):
        frozen = evaluate(self.fit)
        self.assertEqual(frozen["fit"]["oriented_threshold"], 2.)
        verified = self.score_report()
        self.assertTrue(verified["rank_prediction_permitted"])
        self.assertEqual(verified["heldout_class_coverage"], ["stable","unstable"])
        self.assertEqual(verified["shared_source_sha256"], ["a"*64])
        self.score["heldout"]["2"]["score"].update(lower=.5, upper=1.5)
        failure = self.score_report()
        self.assertFalse(failure["rank_prediction_permitted"])
        self.assertEqual(_disposition("mechanism",failure), "fail")
        self.assertEqual(evaluate(self.fit)["fit_sha256"],frozen["fit_sha256"])

    def test_tie_uncertainty_and_absence(self):
        self.score["heldout"]["2"]["score"].update(lower=2., upper=2.)
        self.assertEqual(self.score_report()["status"],"indeterminate")
        del self.score["heldout"]["3"]
        self.assertEqual(self.score_report()["status"],"missing")
        self.fit["training"]["0"]["score"].update(lower=0, upper=4)
        self.assertEqual(evaluate(self.fit)["status"],"indeterminate")
        self.fit["training"]["0"]["score"].update(lower=3, upper=3)
        self.assertEqual(evaluate(self.fit)["status"],"fail")

    def test_frozen_split_and_source_and_model(self):
        bad = copy.deepcopy(self.fit)
        bad["training"]["3"] = self.score["heldout"]["3"]
        with self.assertRaisesRegex(ValueError,"outside"):
            evaluate(bad)
        self.reg["anchors"]["2"]["group_id"] = "g1"
        with self.assertRaisesRegex(ValueError,"commitment"):
            evaluate(self.fit)
        self.fit["registration_sha256"] = _hash(self.reg)
        with self.assertRaisesRegex(ValueError,"overlap"):
            evaluate(self.fit)

    def test_tampered_fit_and_source_label(self):
        frozen = evaluate(self.fit)
        frozen["fit"]["oriented_threshold"] = 100.
        with self.assertRaisesRegex(ValueError,"calibration artifact"):
            evaluate({**self.score, "calibration":frozen})
        self.score["heldout"]["2"]["observation"]["label"] = "unstable"
        with self.assertRaisesRegex(ValueError,"frozen source"):
            self.score_report()
