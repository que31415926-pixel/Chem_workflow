"""Run one staged task over an existing, verified OpenSSH connection."""

import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import shlex
import subprocess
import sys
import tarfile
import tempfile


class TransportError(RuntimeError):
    """Remote execution or evidence transfer could not be verified."""


_MANIFEST = "__chemflow_transfer__.json"
_TASK = "chemflow-task.json"
_RESULT = "chemflow-result.json"
_LOGS = ("chemflow-worker.out.log", "chemflow-worker.err.log")


def _relative(name):
    if not isinstance(name, str) or not name or "\\" in name or "\x00" in name:
        raise TransportError(f"Invalid archive path: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or str(path) != name or name == ".":
        raise TransportError(f"Invalid archive path: {name!r}")
    return path


def _file(root, name):
    relative = _relative(name)
    path = root / relative
    # Resolve the root once so a symlinked temporary directory (macOS: /var -> /private/var)
    # is accepted, while a symlink anywhere below the root is still refused.
    if path.resolve() != Path(root).resolve() / relative:
        raise TransportError(f"Symlink in transfer path: {name!r}")
    return path


class _HashingReader:
    def __init__(self, stream):
        self.stream = stream
        self.sha = hashlib.sha256()

    def read(self, size):
        data = self.stream.read(size)
        self.sha.update(data)
        return data


def _pack(root, names, stream):
    hashes = {}
    with tarfile.open(fileobj=stream, mode="w|gz", compresslevel=1) as archive:
        for name in sorted(names):
            path = _file(root, name)
            if name == _MANIFEST or name in hashes or not path.is_file():
                raise TransportError(f"Invalid transfer member: {name!r}")
            info = tarfile.TarInfo(name)
            info.size = path.stat().st_size
            with path.open("rb") as source:
                reader = _HashingReader(source)
                archive.addfile(info, reader)
                hashes[name] = reader.sha.hexdigest()
        manifest = json.dumps(hashes, sort_keys=True).encode()
        info = tarfile.TarInfo(_MANIFEST)
        info.size = len(manifest)
        archive.addfile(info, io.BytesIO(manifest))


def _unpack(stream, dest):
    """Extract regular files into a new private directory, then verify all bytes."""
    seen = {}
    expected = None
    with tarfile.open(fileobj=stream, mode="r|gz") as archive:
        for member in archive:
            path = _file(dest, member.name)
            if member.type not in (tarfile.REGTYPE, tarfile.AREGTYPE) or expected is not None:
                raise TransportError(f"Invalid transfer member: {member.name!r}")
            with archive.extractfile(member) as source:
                if member.name == _MANIFEST:
                    if member.size > 8 * 1024 * 1024:
                        raise TransportError("Transfer manifest too large")
                    expected = json.load(source)
                    if not isinstance(expected, dict):
                        raise TransportError("Invalid transfer manifest")
                    continue
                path.parent.mkdir(parents=True, exist_ok=True)
                sha = hashlib.sha256()
                with path.open("xb") as out:
                    for data in iter(lambda: source.read(1024 * 1024), b""):
                        out.write(data)
                        sha.update(data)
                seen[member.name] = sha.hexdigest()
    if seen != expected:
        raise TransportError("Transfer manifest missing or checksum mismatch")
    return seen


def _ssh(host, argv):
    if not isinstance(host, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", host):
        raise TransportError("Node host must be an OpenSSH alias")
    return ["ssh", "-T", "-o", "BatchMode=yes", "-o", "StrictHostKeyChecking=yes",
            "-o", "ConnectTimeout=15", "-o", "ServerAliveInterval=30",
            "-o", "ServerAliveCountMax=3", "-o", "ForwardAgent=no",
            "-o", "ClearAllForwardings=yes", "--", host, shlex.join(argv)]


def _check(proc):
    if proc.returncode:
        message = proc.stderr.decode("utf-8", "replace")[-4000:]
        raise TransportError(f"SSH or remote command failed ({proc.returncode}): {message}")


def client_call(host, argv, *, stdin=None):
    """Call the installed remote CLI. Both stdin and the returned stdout are bytes."""
    proc = subprocess.run(_ssh(host, ["python3", "-m", "chemflow2", *argv]),
                          input=b"" if stdin is None else stdin,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _check(proc)
    return proc.stdout


def client_stream(host, argv, *, stdin=None, stdout=None):
    """Exchange a submission or evidence bundle without loading it into memory."""
    proc = subprocess.run(_ssh(host, ["python3", "-m", "chemflow2", *argv]),
                          stdin=subprocess.DEVNULL if stdin is None else stdin,
                          stdout=subprocess.DEVNULL if stdout is None else stdout,
                          stderr=subprocess.PIPE)
    _check(proc)


def _result(root, verified=None):
    path = _file(root, _RESULT)
    result = json.loads(path.read_bytes())
    if (not isinstance(result, dict) or type(result.get("ok")) is not bool
            or result.get("status") != ("succeeded" if result["ok"] else "failed")
            or not isinstance(result.get("artifacts"), dict)):
        raise TransportError("Invalid worker result")
    for name, digest in result["artifacts"].items():
        if name in (_TASK, _RESULT, _MANIFEST, *_LOGS):
            raise TransportError(f"Reserved artifact name: {name}")
        artifact = _file(root, name)
        if not artifact.is_file():
            raise TransportError(f"Missing artifact: {name!r}")
        if verified is None:
            with artifact.open("rb") as stream:
                reader = _HashingReader(stream)
                while reader.read(1024 * 1024):
                    pass
                computed = reader.sha.hexdigest()
        else:
            computed = verified.get(name)
        if computed != digest:
            raise TransportError(f"Artifact missing or checksum differs: {name!r}")
    return result


def _worker_bundle(source, target):
    """Remote endpoint. Its standard output is reserved exclusively for the bundle."""
    with tempfile.TemporaryDirectory(prefix="chemflow-worker-") as tmp:
        private = Path(tmp).resolve()
        stage = private / "case"
        stage.mkdir()
        initial = _unpack(source, stage)
        if _TASK not in initial or any(name in initial for name in (_RESULT, *_LOGS)):
            raise TransportError("Invalid input bundle control files")

        with (private / _LOGS[0]).open("xb") as stdout, (private / _LOGS[1]).open("xb") as stderr:
            # Import the installed application before entering a case-controlled cwd.
            # This also works with a user-site install on Python 3.10.
            launcher = (
                "import os,runpy,sys,chemflow2;"
                "sys.path[:]=[os.path.abspath(p) for p in sys.path];"
                "os.chdir(sys.argv[1]);sys.argv=['chemflow2','worker'];"
                "runpy.run_module('chemflow2',run_name='__main__',alter_sys=True)"
            )
            proc = subprocess.run([sys.executable, "-c", launcher, str(stage)],
                                  stdin=subprocess.DEVNULL, stdout=stdout, stderr=stderr)

        for name in _LOGS:
            os.replace(private / name, stage / name)
        try:
            result = _result(stage)
            if proc.returncode != (0 if result["ok"] else 1):
                raise TransportError(f"Worker status and exit code disagree ({proc.returncode})")
        except (OSError, ValueError, TransportError) as exc:
            result = {"ok": False, "status": "failed", "error": f"Unverified worker result: {exc}",
                      "failure_kind": "transport", "artifacts": {}, "observations": {},
                      "provenance": {"worker_exit_code": proc.returncode}}
            (stage / _RESULT).write_text(json.dumps(result), encoding="utf-8")
        _pack(stage, [*result["artifacts"], _RESULT, *_LOGS], target)


def execute(node, task, folder):
    """Stage one attempt, execute the installed worker, fetch verified output and logs."""
    folder = Path(folder).resolve()
    host = node.get("host") or node.get("ssh_host")
    command = _ssh(host, ["python3", "-m", "chemflow2.transport", "worker-bundle"])
    taskfile = _file(folder, _TASK)
    if taskfile.exists():
        if json.loads(taskfile.read_bytes()) != task:
            raise TransportError("Staged task JSON does not match the requested task")
    else:
        taskfile.write_text(json.dumps(task), encoding="utf-8")
    names = []
    for item in folder.rglob("*"):
        if item.is_symlink() or not item.is_dir():
            names.append(item.relative_to(folder).as_posix())
    if any(name in names for name in (_RESULT, *_LOGS)):
        raise TransportError("Stage each attempt in a fresh folder")

    with tempfile.TemporaryFile() as incoming, tempfile.TemporaryFile() as outgoing:
        _pack(folder, names, outgoing)
        outgoing.seek(0)
        proc = subprocess.run(command, stdin=outgoing, stdout=incoming, stderr=subprocess.PIPE)
        _check(proc)
        incoming.seek(0)
        with tempfile.TemporaryDirectory(prefix="chemflow-receive-", dir=folder.parent) as tmp:
            stage = Path(tmp).resolve()
            received = _unpack(incoming, stage)
            result = _result(stage, received)
            if set(received) != {*result["artifacts"], _RESULT, *_LOGS}:
                raise TransportError("Returned files do not match the worker result")
            for name in received:
                dest = _file(folder, name)
                if dest.exists():
                    raise TransportError(f"Output collides with staged input: {name!r}")
            # Publish the result marker last, after every declared file is in place.
            for name in sorted(received, key=lambda key: key == _RESULT):
                dest = folder / name
                dest.parent.mkdir(parents=True, exist_ok=True)
                os.replace(stage / name, dest)
    return result


if __name__ == "__main__":
    if sys.argv[1:] != ["worker-bundle"]:
        raise SystemExit("Expected: python3 -m chemflow2.transport worker-bundle")
    try:
        _worker_bundle(sys.stdin.buffer, sys.stdout.buffer)
    except (OSError, ValueError, tarfile.TarError, TransportError) as error:
        print(f"Worker transfer failed: {error}", file=sys.stderr)
        raise SystemExit(2)
