"""Run one cloud slot, committing each scientific task separately to private storage."""
import argparse
import hashlib
import io
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import re
import shutil
import tempfile
import uuid

from . import controller as c
from .batch_dag import encoded, sha
from .protocol import load_plan
from .transport import _pack, _result, _unpack


class BucketStore:
    def __init__(self, bucket, run):
        from google.cloud import storage
        self.bucket = storage.Client().bucket(bucket)
        self.prefix = f"chemflow/{run}/"

    def get(self, name, stream):
        self.bucket.blob(self.prefix + name).download_to_file(stream)

    def read(self, name):
        from google.api_core.exceptions import NotFound
        try:
            return self.bucket.blob(self.prefix + name).download_as_bytes()
        except NotFound:
            return None

    def create(self, name, stream):
        from google.api_core.exceptions import PreconditionFailed
        try:
            self.bucket.blob(self.prefix + name).upload_from_file(stream, if_generation_match=0,
                content_type="application/json" if name.endswith(".json") else "application/gzip")
            return True
        except PreconditionFailed:
            return False


def accept_source(store, root, run, input_hash):
    with tempfile.TemporaryFile() as wire:
        store.get("input.tgz", wire)
        wire.seek(0)
        digest = hashlib.sha256()
        for block in iter(lambda: wire.read(2**20), b""):
            digest.update(block)
        if digest.hexdigest() != input_hash:
            raise ValueError("input checksum mismatch")
        wire.seek(0)
        c.accept(root, run, wire)
    return c.run_path(root, run)


def verify_commit(raw, key, run, input_hash):
    if raw is None:
        raise ValueError(f"no committed result for {key}")
    marker = json.loads(raw)
    if (marker["run"] != run or marker["task_id"] != key or marker["input_hash"] != input_hash
            or not re.fullmatch(r"[0-9a-f]{32}", marker["attempt"])
            or not re.fullmatch(r"[0-9a-f]{64}", marker["sha256"]) or type(marker["ok"]) is not bool):
        raise ValueError("published task identity mismatch")
    return marker


def restore(store, directory, key, run, input_hash, raw):
    marker = verify_commit(raw, key, run, input_hash)
    directory.mkdir()
    # Keep compressed transfer data outside the job tree and filesystem namespace.
    with tempfile.TemporaryFile() as wire:
        store.get(f"attempts/{key}/{marker['attempt']}.tgz", wire)
        wire.seek(0)
        digest = hashlib.sha256()
        for block in iter(lambda: wire.read(2**20), b""):
            digest.update(block)
        if digest.hexdigest() != marker["sha256"]:
            raise ValueError("published task checksum mismatch")
        wire.seek(0)
        verified = _unpack(wire, directory)
    result = _result(directory, verified)
    if result.get("run_id") != run or result.get("task_id") != key or result["ok"] != marker["ok"]:
        raise ValueError("published task status mismatch")
    return result


def completed(result):
    # A completed gate may falsify its hypothesis. Its evidence remains available
    # while the cloud workflow skips only tasks that depend on that gate.
    return result["ok"] or result.get("failure_kind") in {
        "scientific_fail", "scientific_missing", "scientific_indeterminate"}


def run_slot(store, run, control, index, deadline=None):
    """A queued VM rechecks the reviewed reservation before each scientific task."""
    if control["run"] != run or not (0 <= index < len(control["slots"])):
        raise ValueError("invalid cloud task slot")
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        input_hash = control["input_hash"]
        work = accept_source(store, root, run, input_hash)
        plan = load_plan(work / "source/protocol.json")
        tasks = {task["id"]: task for task in plan["tasks"]}
        all_ok = True
        checked = {}
        for key in control["slots"][index]:
            task = tasks[key]
            folder = work / "jobs" / key
            receipt = f"commits/{key}.json"
            existing = store.read(receipt)
            if existing is not None:
                restored = restore(store, folder, key, run, input_hash, existing)
                all_ok = completed(restored) and all_ok
                continue
            data_parents = {ref.split("/", 1)[0] for ref in task["from_deps"].values()}
            for dep in task["deps"]:
                if dep not in checked:
                    raw = store.read(f"commits/{dep}.json")
                    if not verify_commit(raw, dep, run, input_hash)["ok"]:
                        raise ValueError(f"required task {dep} failed")
                    checked[dep] = raw
                target = folder.parent / dep
                if dep in data_parents and not target.exists():
                    restore(store, target, dep, run, input_hash, checked[dep])
            if deadline is not None and datetime.now(timezone.utc) >= deadline:
                raise ValueError("cloud dispatch authorization expired while queued")
            result = c.attempt(work, task, control["node"], c.registered_worker)
            result.update(run_id=run, task_id=key)
            c.atomic(folder / "chemflow-result.json", result)
            attempt = uuid.uuid4().hex
            with tempfile.NamedTemporaryFile(dir=root) as wire:
                names = [p.relative_to(folder).as_posix() for p in folder.rglob("*") if p.is_file()]
                _pack(folder, names, wire)
                wire.flush()
                digest = c.digest(Path(wire.name))
                wire.seek(0)
                if not store.create(f"attempts/{key}/{attempt}.tgz", wire):
                    raise ValueError("attempt ID collision")
            marker = {"run": run, "task_id": key, "input_hash": input_hash,
                      "attempt": attempt, "sha256": digest, "ok": result["ok"]}
            if not store.create(receipt, io.BytesIO(encoded(marker))):
                # Concurrent execution cannot replace the winner or publish partial data.
                shutil.rmtree(folder)
                result = restore(store, folder, key, run, input_hash, store.read(receipt))
            all_ok = completed(result) and all_ok
            # Each member in a slot belongs to the same dependency frontier.
            # Parent artifacts can be reused by its siblings; free member scratch now.
            shutil.rmtree(folder)
        return all_ok


def main():
    cli = argparse.ArgumentParser()
    for name in ("bucket", "run", "piece", "digest"):
        cli.add_argument(name)
    args = cli.parse_args()
    if not re.fullmatch(r"cf-[0-9a-f]{32}", args.run):
        raise ValueError("invalid run ID")
    store = BucketStore(args.bucket, args.run)
    raw = store.read(args.piece)
    if raw is None or sha(raw) != args.digest:
        raise ValueError("cloud slot descriptor checksum mismatch")
    control = json.loads(raw)
    meta = json.loads(store.read("submission.json"))
    if meta["run"] != args.run or meta["input_hash"] != control["input_hash"]:
        raise ValueError("run reservation identity mismatch")
    from .budget import timestamp
    deadline = timestamp(meta["reservation"]["valid_until"])
    return 0 if run_slot(store, args.run, control, int(os.environ["BATCH_TASK_INDEX"]), deadline) else 1


if __name__ == "__main__":
    raise SystemExit(main())
