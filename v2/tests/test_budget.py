import hashlib
import json
from datetime import timedelta
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
import subprocess

from chemflow2 import budget
from chemflow2.usage import UsageSampler


class BudgetTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.profile = {"project":"p", "bucket":"b", "region":"r", "budget_policy":"budget.json"}
        proof = self.write("proof",b"unit test evidence")
        self.policy = dict(self.profile,currency="USD",observed_at=(budget.now()-timedelta(hours=1)).isoformat(),
                           valid_until=(budget.now()+timedelta(days=10)).isoformat(), vm_overhead_hours="0.5",
                           vm_usd_per_hour={"m":"1.125"},disk_usd_per_gb_hour="0.01", tail_hours=4,
                           ancillary_upper_usd={k:"0.04" for k in budget.EXTRAS},credit_remaining_usd=10,
                           spend_limit_usd=8, evidence={k:proof for k in ["account","tariff","other_cost_bounds"]})
        for key in ["provider_blocks_paid_usage","trial_eligibility_confirmed","billing_current",
                    "dedicated_account","retention_and_cleanup_reviewed"]:
            self.policy[key]=True
        self.alloc = [{"id":"wave","task_ids":["a","b"],"machine_type":"m","max_instances":2,"max_hours":2,"disk_gb":10}]

    def write(self,name,data):
        (self.root/name).write_bytes(data)
        return {"file":name,"sha256":hashlib.sha256(data).hexdigest()}

    def read(self):
        self.write("budget.json",json.dumps(self.policy).encode())
        return budget.policy_for(self.profile,self.root/"profile.json")

    def test_fail_closed_unknown_and_tamper(self):
        self.policy["billing_current"]=False
        with self.assertRaisesRegex(ValueError,"billing_current"): self.read()
        self.policy["billing_current"]=True
        self.write("proof",b"modified")
        with self.assertRaisesRegex(ValueError,"provider evidence"): self.read()

    def test_full_reservation_shortfall_and_idempotence(self):
        policy = self.read()
        quote=budget.quote(policy,self.alloc)
        self.assertEqual(quote["upper_usd"],"6.41")
        ledger={"project":"p","currency":"USD","holds":{}}
        budget.admission(policy,ledger,"run1",quote)
        budget.admission(policy,ledger,"run1",quote)
        with self.assertRaisesRegex(ValueError,"shortfall USD 4.82"):
            budget.admission(policy,ledger,"run2",quote)
        self.assertEqual(set(ledger["holds"]),{"run1"})

    def test_expired_coverage_and_duplicate_tasks(self):
        self.policy["valid_until"]=(budget.now()+timedelta(hours=4)).isoformat()
        policy=self.read()
        with self.assertRaisesRegex(ValueError,"lifetime"):
            budget.admission(policy,{}, "r",budget.quote(policy,self.alloc))
        with self.assertRaisesRegex(ValueError,"exactly once"):
            budget.quote(policy,self.alloc*2)

    def test_scratch_sampling_and_ledger_reconciliation(self):
        sampler=UsageSampler(self.root)
        sampler.start()
        self.write("attempt",b"energy"*3000)
        self.assertGreaterEqual(sampler.finish()["attempt_scratch_allocated_sampled_max_bytes"],8192)
        receipt={"run_id":"r","task_id":"SP","metrics":{"wall_seconds":12},"provenance":{"task_id":"SP","stage":"sp"},
                 "status":"failed","observations":{"engine_version":"6.1.1"}}
        export=b"line_id,project,run_id,currency,gross_usd,credit_usd\nslot,p,r,USD,0.151,0.101\n"
        rec=self.write("r.json",json.dumps(receipt).encode())
        billing=self.write("bill.csv",export)
        path=self.root/"manifest.json"
        path.write_text(json.dumps(dict(run_id="r",project="p",billing_export=billing,receipts=[rec])))
        self.assertEqual(budget.reconcile(path)["reported_net_usd"],"0.050")

    def test_competing_reservation_refuses_dispatch(self):
        self.read()
        ledger = {"project":"p","currency":"USD","holds":{}}
        def transfer(args, **kwargs):
            if "cp" in args and "--if-generation-match=0" in args:
                return subprocess.CompletedProcess(args,1)
            if "cp" in args and args[-2].endswith("#75"):
                Path(args[-1]).write_text(json.dumps(ledger))
                return subprocess.CompletedProcess(args,0)
            if "--if-generation-match=75" in args:
                raise subprocess.CalledProcessError(1,args)
            raise AssertionError(args)
        with patch("chemflow2.budget.subprocess.run",side_effect=transfer), patch(
                "chemflow2.budget.subprocess.check_output",return_value="75"):
            with self.assertRaises(subprocess.CalledProcessError):
                budget.reserve_dispatch(self.profile,self.root/"profile.json","r",self.alloc)
