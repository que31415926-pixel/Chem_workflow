import json
from pathlib import Path
import tempfile
import unittest

from chemflow2 import batch, controller as c
from chemflow2.protocol import load_plan
from chemflow2.transport import _unpack

ROOT=Path(__file__).resolve().parents[1]
PROTOCOL=ROOT/"examples/probe/protocol.json"
RUN="cf-"+"f"*32


class BatchTests(unittest.TestCase):
    def test_offline_worker_and_verified_return(self):
        with tempfile.TemporaryDirectory() as tmp:
            mount=Path(tmp); remote=mount/"chemflow"/RUN
            remote.mkdir(parents=True)
            wire=remote/"input.tgz"
            with wire.open("wb") as out:
                c.bundle(PROTOCOL,out)
            node=batch.check_profile(load_plan(PROTOCOL),c.read(ROOT/"examples/batch-profile.example.json"))
            with self.assertRaisesRegex(ValueError,"checksum"):
                batch.execute(mount,RUN,"wrong",node)
            self.assertTrue(batch.execute(mount,RUN,c.digest(wire),node))
            dst=mount/"verified";dst.mkdir()
            with (remote/"result.tgz").open("rb") as stream:
                _unpack(stream,dst)
            self.assertEqual(c.read(dst/"state.json")["status"],"finished")
            self.assertEqual(c.digest(dst/"jobs/produce/sample.json"),c.digest(dst/"jobs/consume/parent.json"))

    def test_spec_and_capacity(self):
        profile=c.read(ROOT/"examples/batch-profile.example.json");plan=load_plan(PROTOCOL)
        spec=batch.job_spec(plan,profile,RUN,"a"*64)
        self.assertEqual(spec["taskGroups"][0]["taskSpec"]["maxRetryCount"],0)
        profile["node"]["engines"]=[]
        with self.assertRaises(ValueError):
            batch.job_spec(plan,profile,RUN,"a"*64)
