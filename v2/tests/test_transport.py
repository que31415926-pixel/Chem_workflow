import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile
import tempfile
import unittest
from unittest.mock import patch

from chemflow2 import transport


class TransportTests(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory()
        self.addCleanup(self.workspace.cleanup)
        self.base = Path(self.workspace.name).resolve()
        self.case = self.base / "case"
        self.case.mkdir()
        (self.case / "input.xyz").write_bytes(b"fixture")
        self.task = {"id": "test", "outputs": ["answer.xyz"]}

    def _simulate(self, ok, corrupt=False):
        def remote(command, **kwargs):
            if "stdin" in kwargs and "stdout" in kwargs and "worker-bundle" in command[-1]:
                transport._worker_bundle(kwargs["stdin"], kwargs["stdout"])
                return subprocess.CompletedProcess(command, 0, stderr=b"")
            stage = Path(command[-1])
            self.assertEqual({p.name for p in stage.iterdir()}, {"input.xyz", "chemflow-task.json"})
            self.assertEqual(json.loads((stage / "chemflow-task.json").read_bytes()), self.task)
            (stage / "scratch.tmp").write_bytes(b"excluded")
            (stage / "answer.xyz").write_bytes(b"result")
            kwargs["stderr"].write(b"worker diagnostics")
            digest = hashlib.sha256(b"result").hexdigest()
            result = {"ok": ok, "status": "succeeded" if ok else "failed", "error": None,
                      "failure_kind": None, "artifacts": {"answer.xyz": "wrong" if corrupt else digest},
                      "observations": {}, "provenance": {}}
            (stage / "chemflow-result.json").write_text(json.dumps(result))
            return subprocess.CompletedProcess(command, 0 if ok else 1)
        return remote

    def test_submit_and_fetch_verified_success(self):
        with patch.object(subprocess, "run", side_effect=self._simulate(True)):
            result = transport.execute({"host": "campus-1"}, self.task, self.case)
        self.assertTrue(result["ok"])
        self.assertEqual((self.case / "answer.xyz").read_bytes(), b"result")
        self.assertFalse((self.case / "scratch.tmp").exists())
        self.assertEqual((self.case / "chemflow-worker.err.log").read_bytes(), b"worker diagnostics")

    def test_failed_task_preserves_artifacts_and_diagnostics(self):
        with patch.object(subprocess, "run", side_effect=self._simulate(False)):
            result = transport.execute({"host": "campus-1"}, self.task, self.case)
        self.assertFalse(result["ok"])
        self.assertEqual((self.case / "answer.xyz").read_bytes(), b"result")

    def test_declared_checksum_mismatch_rejects_result(self):
        with patch.object(subprocess, "run", side_effect=self._simulate(True, corrupt=True)):
            result = transport.execute({"host": "campus-1"}, self.task, self.case)
        self.assertFalse(result["ok"])
        self.assertEqual(result["failure_kind"], "transport")
        self.assertFalse((self.case / "answer.xyz").exists())

    def test_wire_manifest_and_paths_are_verified(self):
        for name, link in (("../escape", False), ("symlink", True), ("valid", False)):
            with self.subTest(name=name):
                wire = io.BytesIO()
                with tarfile.open(fileobj=wire, mode="w:gz") as archive:
                    info = tarfile.TarInfo(name)
                    info.type = tarfile.SYMTYPE if link else tarfile.REGTYPE
                    info.linkname = "../escape" if link else ""
                    archive.addfile(info)
                    manifest = json.dumps({name: "bad digest"}).encode()
                    info = tarfile.TarInfo("__chemflow_transfer__.json")
                    info.size = len(manifest)
                    archive.addfile(info, io.BytesIO(manifest))
                wire.seek(0)
                with tempfile.TemporaryDirectory(dir=self.base) as tmp:
                    with self.assertRaises(transport.TransportError):
                        transport._unpack(wire, Path(tmp))
        self.assertFalse((self.base / "escape").exists())

    def test_symlink_input_refused_before_connecting(self):
        (self.case / "linked").symlink_to(self.case / "input.xyz")
        with patch.object(subprocess, "run") as ssh:
            with self.assertRaises(transport.TransportError):
                transport.execute({"host": "campus-1"}, self.task, self.case)
        ssh.assert_not_called()

    def test_alias_and_shell_quoting(self):
        for bad_alias in ("-oProxyCommand=false", "user@node", "a;echo", ""):
            with self.assertRaises(transport.TransportError):
                transport._ssh(bad_alias, ["python3"])
        with patch.object(subprocess, "run", return_value=subprocess.CompletedProcess([], 0, b"ok", b"")) as ssh:
            value = transport.client_call("campus-1", ["serve-stdio", "status", "a b"], stdin=b"bytes")
        self.assertEqual(value, b"ok")
        argv = ssh.call_args.args[0]
        self.assertIn("StrictHostKeyChecking=yes", argv)
        self.assertIn("BatchMode=yes", argv)
        self.assertIn("'a b'", argv[-1])

    def test_ssh_failure_has_diagnostic(self):
        proc = subprocess.CompletedProcess([], 255, b"", b"Host key verification failed")
        with patch.object(subprocess, "run", return_value=proc):
            with self.assertRaisesRegex(transport.TransportError, "Host key verification failed"):
                transport.client_call("campus-1", ["serve-stdio", "status", "run"])
