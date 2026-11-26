import io
import json
import tarfile
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

from chemflow2 import cloud_setup as setup


class BootstrapTests(unittest.TestCase):
    def test_offline_plan_does_not_touch_cloud(self):
        with patch.object(setup, "call", side_effect=AssertionError), redirect_stdout(io.StringIO()) as out:
            self.assertEqual(setup.main(["trial-project"]), 0)
        result = json.loads(out.getvalue())
        self.assertEqual(result["profile"]["bucket"], "trial-project-chemflow-v3")
        self.assertIn("workflowexecutions", result["required_apis"])
        self.assertEqual(result["build_machine"], setup.build_allocation()["machine_type"])

    def test_staging_allowlist_omits_research_secrets_and_projects(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("chemflow2/__init__.py", "chemflow2/private.json", "chemflow/__main__.py",
                         "deploy/byol.Dockerfile", "projects/secret/geometry.xyz", "credentials.json"):
                (root/name).parent.mkdir(parents=True, exist_ok=True)
                (root/name).write_text("fixture")
            installer = root/"orca_6_1_0_linux_x86-64_shared_openmpi416.run"
            installer.write_bytes(b"fixture")
            target = root/"source.tgz"
            setup.stage(root, installer, target)
            with tarfile.open(target, "r:gz") as tar:
                self.assertEqual(set(tar.getnames()), {"Dockerfile", "orca.run", "chemflow2/__init__.py",
                    "vendor/chemflow-v1/chemflow/__main__.py"})

    def test_live_build_fails_before_cloud_if_budget_unverified(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            installer = root/"orca_6_1_0_linux_x86-64_shared_openmpi416.run"
            installer.write_bytes(b"fixture")
            with patch.object(setup, "call", side_effect=AssertionError), self.assertRaises(FileNotFoundError):
                setup.main(["trial-project", "--apply", "--byol-confirmed", "--installer", str(installer),
                            "--profile", str(root/"profile.json"), "--budget-policy", str(root/"absent")])

    def test_source_deleted_when_build_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            profile = setup.settings("trial-project", "us-central1", root/"policy.json")
            calls = []
            def gcloud(project, *args, **kwargs):
                calls.append(args)
                if args[:2] == ("builds", "submit"):
                    raise RuntimeError("build failed")
            with patch.object(setup, "call", side_effect=gcloud), patch.object(setup, "stage"), \
                    patch.object(setup.budget, "reserve_dispatch") as reserve:
                with self.assertRaisesRegex(RuntimeError, "build failed"):
                    setup.build_image(profile, root/"installer", root, root, root/"profile.json")
            self.assertEqual(reserve.call_count, 1)
            self.assertEqual(calls[0][:2], ("storage", "cp"))
            self.assertEqual(calls[-1], ("storage", "rm", calls[0][-1]))

    def test_retry_checks_resource_identity_and_never_adopts_replacements(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            profile = setup.settings("trial-project", "us-central1", p/"policy.json")
            journal = setup.ResourceJournal(profile, p/"journal.json")
            with patch.object(setup, "call") as gcloud, patch.object(journal, "info", return_value="123"):
                journal.create("job", "create")
                journal.create("job", "create")
                self.assertEqual(gcloud.call_count, 1)
            loaded = setup.ResourceJournal(profile, journal.path)
            with patch.object(setup, "call") as gcloud, patch.object(loaded, "info", return_value="456"):
                with self.assertRaisesRegex(ValueError, "replaced"):
                    loaded.remove("job", "delete")
                gcloud.assert_not_called()

    def test_teardown_refuses_active_work_and_preserves_evidence_by_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            profile = setup.settings("trial-project", "us-central1", p/"policy.json")
            journal = setup.ResourceJournal(profile, p/"journal.json")
            journal.data["created"] = {"bucket":"identity"}
            journal.save()
            with patch.object(setup, "call", return_value="run"):
                with self.assertRaisesRegex(ValueError, "running"):
                    setup.teardown(profile, journal, False)
            with patch.object(setup, "call", return_value=""), redirect_stdout(io.StringIO()):
                setup.teardown(profile, journal, False)
            self.assertEqual(journal.data["created"], {"bucket":"identity"})

    def test_quota_limits_use_remaining_regional_cpu_and_ssd(self):
        snapshot = {"quotas":[{"metric":"CPUS", "limit":10, "usage":2},
                              {"metric":"SSD_TOTAL_GB", "limit":250, "usage":18}]}
        self.assertEqual(setup.quota_limits(snapshot), {"regional_vcpu_limit":8, "regional_disk_limit_gb":232})
        snapshot["quotas"][0]["usage"] = 10
        with self.assertRaisesRegex(ValueError, "no remaining"):
            setup.quota_limits(snapshot)
