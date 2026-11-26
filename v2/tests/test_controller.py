import io
import json
from pathlib import Path
import tempfile
import threading
import time
import unittest

from chemflow2 import controller as c
from chemflow2.transport import _pack
from chemflow2.worker import run_task

NODE = {"name":"n","host":"n","engines":["probe"],"cores":1,"memory_mb":64,
        "scratch_mb":2,"concurrency":4}
PROBE = Path(__file__).resolve().parents[1]/"examples/probe/protocol.json"


def submit(root, rid, protocol=PROBE):
    wire=io.BytesIO()
    c.bundle(protocol,wire)
    wire.seek(0)
    c.accept(root,rid,wire)


def worker(node,task,folder):
    return run_task(task,folder)


class TestController(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)

    def test_scheduling_uses_global_resource_reservations(self):
        for rid in ("one","two"):
            submit(self.root,rid)
        count=peak=0
        lock=threading.Lock()
        def watch(node,task,path):
            nonlocal count,peak
            with lock:
                count+=1
                peak=max(peak,count)
            time.sleep(0.012)
            receipt=run_task(task,path)
            with lock:
                count-=1
            return receipt
        c.serve(self.root,[NODE],once=True,interval=.001,executor=watch)
        self.assertEqual(peak,1)
        for rid in ("one","two"):
            run=c.run_path(self.root,rid)
            self.assertEqual(c.read(run/"state.json")["status"],"finished")
            self.assertEqual(c.digest(run/"jobs/produce/sample.json"),c.digest(run/"jobs/consume/parent.json"))

    def test_missing_machine_blocks_dag(self):
        submit(self.root,"empty")
        c.serve(self.root,[],once=True,interval=0,executor=worker)
        result=c.read(c.run_path(self.root,"empty")/"state.json")
        self.assertEqual(result["status"],"blocked")
        self.assertEqual(result["tasks"]["consume"]["error"],"incomplete dependency")

    def test_one_damaged_protocol_does_not_stop_other_run(self):
        for name in ("bad","ok"):
            submit(self.root,name)
        (c.run_path(self.root,"bad")/"source/protocol.json").unlink()
        c.serve(self.root,[NODE],once=True,interval=.001,executor=worker)
        self.assertEqual(c.read(c.run_path(self.root,"bad")/"state.json")["status"],"failed")
        self.assertEqual(c.read(c.run_path(self.root,"ok")/"state.json")["status"],"finished")

    def test_cancel_stops_queued_tasks(self):
        submit(self.root,"stop")
        c.service(self.root,"cancel","stop",None,None)
        c.serve(self.root,[NODE],once=True,interval=0,executor=worker)
        state=c.read(c.run_path(self.root,"stop")/"state.json")
        self.assertEqual(state["status"],"cancelled")
        self.assertFalse(any((c.run_path(self.root,"stop")/"jobs").iterdir()))

    def test_rejected_extra_file_and_duplicate_id(self):
        submit(self.root,"taken")
        with self.assertRaises(ValueError):
            submit(self.root,"taken")
        stage=self.root/"tmp"
        stage.mkdir()
        (stage/"protocol.json").write_bytes(PROBE.read_bytes())
        (stage/"extra").write_text("undeclared")
        wire=io.BytesIO()
        _pack(stage,["protocol.json","extra"],wire)
        wire.seek(0)
        with self.assertRaises(ValueError):
            c.accept(self.root,"extra",wire)
        self.assertFalse(c.run_path(self.root,"extra").exists())

    def test_dedicated_linux_controller_can_work_without_self_ssh(self):
        import platform
        if platform.system()!="Linux":
            self.skipTest("Linux service entry point")
        submit(self.root,"single")
        c.serve(self.root,[dict(NODE,placement="controller")],once=True,interval=.002)
        run=c.run_path(self.root,"single")
        self.assertEqual(c.read(run/"state.json")["status"],"finished")
        self.assertEqual(c.read(run/"jobs/consume/chemflow-result.json")["observations"]["type"],"transport_probe")

    def test_critical_chain_and_prospective_priorities(self):
        tasks = [{"id":"loose","deps":[]},{"id":"a","deps":[]},{"id":"b","deps":["a"]},
                 {"id":"c","deps":["b"]}]
        self.assertEqual(c.schedule_order(tasks)[0],"a")
        tasks[0]["priority"]=5
        self.assertEqual(c.schedule_order(tasks)[0],"loose")
        tasks[-1]["priority"]=10
        order=c.schedule_order(tasks)
        self.assertEqual(order,["a","b","c","loose"])

    def test_remote_node_registration_is_additive_and_idempotent(self):
        self.assertEqual(c.register_node(self.root,NODE)["status"],"registered")
        self.assertEqual(c.register_node(self.root,NODE)["status"],"already_registered")
        with self.assertRaises(ValueError):
            c.register_node(self.root,dict(NODE,cores=2))
        self.assertEqual(c.load_nodes(self.root/"nodes.json"),[NODE])
        second=dict(NODE,name="other",host="other")
        self.assertEqual(c.register_node(self.root,second)["status"],"registered")
        self.assertEqual(c.load_nodes(self.root/"nodes.json"),[NODE,second])
