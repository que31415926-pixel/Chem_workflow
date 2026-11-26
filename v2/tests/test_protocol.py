import json
import tempfile
import unittest
from pathlib import Path

from chemflow2.protocol import explain_plan, load_plan, matches, validate_node


class ProtocolTests(unittest.TestCase):
    def setUp(self):
        self.work = tempfile.TemporaryDirectory()
        self.addCleanup(self.work.cleanup)
        self.root = Path(self.work.name)
        (self.root / "metal.inp").write_text("reference input", encoding="utf-8")
        self.size = {"cores": 2, "memory_mb": 200, "scratch_mb": 100}
        self.first = {"id": "first", "engine": "probe", "stage": "probe", "outputs": ["opt.xyz"],
                      "resources": self.size.copy()}
        self.second = {"id": "second", "engine": "orca", "stage": "sp", "input": "metal.inp",
                       "outputs": ["orca.out"], "resources": self.size.copy(), "deps": ["first"],
                       "from_deps": {"geom.xyz": "first/opt.xyz"}, "requires": ["protocol-A"]}
        self.site = {"name": "campus", "host": "cluster-1", **self.size, "engines": ["orca", "probe"],
                     "capabilities": ["protocol-A"]}

    def load(self, *tasks):
        path = self.root / "plan.json"
        path.write_text(json.dumps({"name": "example", "tasks": tasks}), encoding="utf-8")
        return load_plan(path)

    def test_exact_fit_and_unordered_dependencies(self):
        plan = self.load(self.second, self.first)
        self.assertTrue(matches(plan["tasks"][0], self.site))
        explanation = explain_plan(plan, [self.site])
        self.assertEqual(explanation["runnable"], 2)
        self.assertEqual(explanation["tasks"][0]["eligible_nodes"], ["campus"])

    def test_capability_and_resource_checks_do_not_fall_back(self):
        for override in ({"engines": ["other"]}, {"cores": 1}, {"memory_mb": 100},
                         {"scratch_mb": 99}, {"capabilities": ["protocol-B"]}):
            with self.subTest(override=override):
                self.assertFalse(matches(self.second, {**self.site, **override}))

    def test_blockage_in_producer_propagates_to_eligible_consumer(self):
        plan = self.load(self.second, self.first)
        result = explain_plan(plan, [{**self.site, "engines": ["orca"]}])
        self.assertEqual(set(result["blocked"]), {"first", "second"})
        self.assertEqual(result["tasks"][0]["blocked_by_dependencies"], ["first"])
        self.assertIn("campus", result["tasks"][1]["blockers"])

    def test_invalid_dependency_source_and_cycle(self):
        alterations = ({"deps": ["unknown"]}, {"from_deps": {"a": "first/absent"}},
                       {"from_deps": {"a": "other/opt.xyz"}}, {"deps": ["second"]})
        for override in alterations:
            with self.subTest(override=override), self.assertRaises(ValueError):
                self.load(self.first, {**self.second, **override})
        with self.assertRaisesRegex(ValueError, "cycle"):
            self.load({**self.first, "deps": ["second"]}, self.second)

    def test_reject_reserved_colliding_or_unsafe_paths(self):
        changes = ({"outputs": ["a", "a/x"]}, {"outputs": ["A", "a"]}, {"outputs": ["metal.inp"]},
                   {"outputs": ["../elsewhere"]}, {"outputs": ["/elsewhere"]},
                   {"outputs": ["dir//file"]}, {"outputs": ["chemflow-task.json"]},
                   {"outputs": ["chemflow-result.json/x"]}, {"outputs": ["orca.out/nested"]},
                   {"files": ["metal.inp"]}, {"from_deps": {"orca.out": "first/opt.xyz"}},
                   {"from_deps": {"metal.inp": "first/opt.xyz"}})
        for change in changes:
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.load(self.first, {**self.second, **change})

    def test_static_input_must_exist_and_stay_inside_protocol_dir(self):
        with self.assertRaisesRegex(ValueError, "missing"):
            self.load(self.first, {**self.second, "input": "missing.inp"})
        with tempfile.TemporaryDirectory() as outside:
            target = Path(outside, "outside.inp")
            target.write_text("outside", encoding="utf-8")
            (self.root / "linked.inp").symlink_to(target)
            with self.assertRaisesRegex(ValueError, "outside"):
                self.load(self.first, {**self.second, "input": "linked.inp"})

    def test_validate_nodes_and_resource_units(self):
        for override in ({"name": "../a"}, {"host": "-oProxyCommand=invalid"},
                         {"host": "one two"}, {"cores": False}, {"scratch_mb": -1},
                         {"memory_mb": 1.5}, {"concurrency": 0}):
            with self.subTest(override=override), self.assertRaises(ValueError):
                validate_node({**self.site, **override})

    def test_freeze_and_validate_attachments(self):
        asset = self.root / "review.txt"
        asset.write_text("original", encoding="utf-8")
        body = {"name": "attachments", "tasks": [self.first], "attachments": ["review.txt"]}
        path = self.root / "plan.json"
        path.write_text(json.dumps(body))
        frozen = load_plan(path)
        self.assertEqual(set(frozen["attachment_sha256"]), {"review.txt"})
        path.write_text(json.dumps(frozen))
        asset.write_text("changed", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            load_plan(path)
        (self.root / "linked.txt").symlink_to(asset)
        for name in ["linked.txt", "../review.txt", "protocol.json", "Protocol.json/sub"]:
            with self.subTest(name=name), self.assertRaises(ValueError):
                path.write_text(json.dumps({**body, "attachments": [name]}))
                load_plan(path)

    def test_duplicate_keys_and_task_ids(self):
        with self.assertRaisesRegex(ValueError, "duplicate task id"):
            self.load(self.first, {**self.first, "id": "First"})
        path = self.root / "repeat.json"
        path.write_text('{"name":"a","name":"b"}', encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "duplicate JSON key"):
            load_plan(path)
