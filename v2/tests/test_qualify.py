import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from chemflow2 import controller as c
from chemflow2.qualify import assess, main


class QualificationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.run = "cf-" + "1" * 32
        self.image = "private/worker@sha256:" + "a" * 64
        self.binary = "0123456789abcdef" * 4
        self.profile = {"project":"proj-123", "region":"us-central1", "bucket":"trial-123",
                        "service_account":"worker@example.test", "image":self.image,
                        "calibration_plan_sha256":"f"*64, "node":{"capabilities":["sp-tag"]}}
        group = self.run + "-a0"
        meta = {k:self.profile[k] for k in ("project", "region", "bucket", "image")}
        meta.update(run=self.run, submitted_plan_sha256="f"*64, input_hash="c"*64,
                    allocations=[{"id":group, "task_ids":["sp"]}])
        c.atomic(self.root/"cloud-submission.json", meta)
        c.atomic(self.root/"state.json", {"id":self.run,"status":"finished"})
        container = {"imageUri":self.image, "entrypoint":"python3",
                     "commands":["-m","chemflow2.batch_worker","trial-123",self.run]}
        job = {"name":"jobs/"+group, "uid":"cloud-job-uid", "status":{"state":"SUCCEEDED"},
               "labels":{"chemflow-run":self.run}, "allocationPolicy":{"serviceAccount":{"email":"worker@example.test"}},
               "taskGroups":[{"taskSpec":{"runnables":[{"container":container}]}}]}
        c.atomic(self.root/"cloud-status.json",{"workflow":True,"executions":[{"state":"SUCCEEDED"}],"jobs":[job]})
        self.folder = self.root/"jobs/sp"
        self.folder.mkdir(parents=True)
        self.deck = self.folder/"sp.inp"
        self.deck.write_text("! HF def2-SVP SP\n* xyz 0 1\nHe 0 0 0\n*\n")
        source = self.root/"source"
        source.mkdir()
        (source/"sp.inp").write_bytes(self.deck.read_bytes())
        c.atomic(source/"protocol.json", {"name":"unit", "study_scope":"engineering_calibration",
            "tasks":[{"id":"sp","engine":"orca","stage":"sp","study_scope":"engineering_calibration",
                      "input":"sp.inp","outputs":["orca.out"], "requires":["sp-tag"],
                      "resources":{"cores":1,"memory_mb":256,"scratch_mb":100}}]})
        log = self.folder/"orca.out"
        log.write_text("Program Version 6.1.0\nFINAL SINGLE POINT ENERGY -2.8\nORCA TERMINATED NORMALLY\n")
        prov = {"task_id":"sp","stage":"sp","engine":"orca","platform":"Linux","architecture":"x86_64",
                "executable":"/opt/orca/orca","executable_sha256":self.binary,"orca_processes":1,
                "exit_code":0,"input_sha256":{"sp.inp":c.digest(self.deck)}}
        c.atomic(self.folder/"chemflow-result.json", {"ok":True,"status":"succeeded","provenance":prov,
            "metrics":{"wall_seconds":3}, "artifacts":{"orca.out":c.digest(log)}})

    def evaluate(self):
        return assess(self.profile,self.root,self.run,"Reviewer",self.binary)

    def test_capabilities_have_explicit_receipts(self):
        qualified = self.evaluate()
        self.assertEqual(qualified["qualification"]["tasks_by_capability"],{"sp-tag":["sp"]})
        self.assertNotIn("calibration_plan_sha256", qualified)
        self.assertIn("calibration_plan_sha256", self.profile)

    def test_foreign_binary_plan_or_image_fails(self):
        with self.assertRaisesRegex(ValueError,"unverified"):
            assess(self.profile,self.root,self.run,"Reviewer","abcdef0123456789"*4)
        self.profile["calibration_plan_sha256"] = "1"*64
        with self.assertRaisesRegex(ValueError,"digest"):
            self.evaluate()
        self.profile["calibration_plan_sha256"] = "f"*64
        info = c.read(self.root/"cloud-status.json")
        del info["jobs"][0]["uid"]
        c.atomic(self.root/"cloud-status.json",info)
        with self.assertRaisesRegex(ValueError,"allocation"):
            self.evaluate()

    def test_no_promotion_without_authenticated_fetch(self):
        private = self.root/"input.json"
        private.write_text(json.dumps(self.profile))
        destination = self.root/"qualified.json"
        with patch("chemflow2.batch.main",side_effect=RuntimeError("offline")):
            with self.assertRaisesRegex(RuntimeError,"offline"):
                main([str(private),self.run,str(destination),"--reviewer","Reviewer",
                      "--authorized-orca-sha256",self.binary])
        self.assertFalse(destination.exists())
