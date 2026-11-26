"""Synthetic integration tests for worker dispatch; no scientific software is invoked."""
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from chemflow2.worker import run_task


class ExtraAdapters(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.case = self.root / 'case'
        self.case.mkdir()
        self.xyz = '2\nexample\nH 0 0 0\nH 0 0 0.7\n'
        self.limits = dict(cores=2, memory_mb=2048, scratch_mb=1)
        self.binary = self.root / 'fixture'
        self.env = patch.dict(os.environ, {'CHEMFLOW_CREST': str(self.binary), 'CHEMFLOW_ORCA': str(self.binary)})
        self.env.start()
        self.addCleanup(self.env.stop)

    def crest_task(self):
        (self.case / 'source.xyz').write_text(self.xyz)
        return dict(id='sampling', engine='crest', stage='conformers', engine_version='3.0.2',
                    input='source.xyz', outputs=['crest_conformers.xyz', 'crest.energies'],
                    resources=self.limits, settings=dict(charge=0, uhf=0, method='gfn2',
                    energy_window_kcal_mol=6.5, solvent_model='alpb', solvent='thf'))

    def set_crest(self, energies='1 0.000\n2 1.200\n', version='3.0.2'):
        self.binary.write_text("#!/usr/bin/env python3\nfrom pathlib import Path\nimport sys,os\n"
                "assert sys.argv[sys.argv.index('-T')+1]=='2' and os.environ['OMP_NUM_THREADS']=='2'\n"
                "assert sys.argv[sys.argv.index('--alpb')+1]=='thf'\n"
                f"Path('crest_conformers.xyz').write_text({self.xyz*2!r})\n"
                f"Path('crest.energies').write_text({energies!r})\n"
                f"print('  Version {version}, Thu')\nprint('CREST terminated normally.')\n")
        self.binary.chmod(0o700)

    def test_crest(self):
        self.set_crest()
        report = run_task(self.crest_task(), self.case)
        self.assertTrue(report['ok'], report)
        self.assertEqual(report['observations']['relative_energies_kcal_mol'], [0,1.2])
        self.assertEqual(set(report['artifacts']), {'crest.energies','crest.out','crest_conformers.xyz'})
        self.assertEqual(report['provenance']['engine_version'], '3.0.2')

    def test_crest_reject_missing_frame(self):
        self.set_crest(energies='1 0.000\n')
        report = run_task(self.crest_task(), self.case)
        self.assertEqual(report['failure_kind'], 'validation')
        self.assertEqual(set(report['artifacts']), {'crest.out'})

    def test_crest_reject_settings_and_incompatible_binary(self):
        task = self.crest_task()
        task['settings']['argv'] = ['anything']
        self.assertEqual(run_task(task, self.case)['failure_kind'], 'input')
        del task['settings']['argv']
        self.set_crest(version='3.1.0')
        self.assertEqual(run_task(task, self.case)['failure_kind'], 'capability')

    def test_orca_json_for_downstream(self):
        (self.case / 'qm').mkdir()
        (self.case / 'qm' / 'run.inp').write_text('! HF SP\n%maxcore 100\n')
        task = dict(id='single', engine='orca', stage='sp', input='qm/run.inp',
                    observations_output='data/qm.json',outputs=['data/qm.json'],resources=self.limits)
        self.binary.write_text("#!/bin/sh\nprintf '%s\\n' 'Program Version 6.1.1' 'FINAL SINGLE POINT ENERGY -1.0' '****ORCA TERMINATED NORMALLY****'\n")
        self.binary.chmod(0o700)
        result = run_task(task, self.case)
        self.assertTrue(result['ok'],result)
        parsed = json.loads((self.case/'data/qm.json').read_text())
        self.assertEqual(parsed['observations']['electronic_energy_eh'], -1)
        self.assertEqual(parsed['log_sha256'], result['artifacts']['orca.out'])

    def test_missing_science_report_is_retained_and_blocks_dependents(self):
        (self.case/'analysis.json').write_text(json.dumps(dict(operation='direction',protocol_id='example',
                unit='kcal/mol',expected_sign=1,left_id='A',right_id='B',required_methods=['M1'],pairs={})))
        task = dict(id='scientific',engine='analysis',stage='science',input='analysis.json',
                    outputs=['report.json'],from_deps={})
        result = run_task(task, self.case)
        self.assertFalse(result['ok'], result)
        self.assertEqual(result['failure_kind'],'scientific_missing')
        self.assertEqual(set(result['artifacts']),{'report.json'})
        self.assertEqual(json.loads((self.case/'report.json').read_text())['status'],'missing')

    def test_upstream_dispatch(self):
        def fake(task, folder):
            (folder/'products').mkdir()
            (folder/'products'/'test.csv').write_text('value\n1\n')
            (folder/'upstream.out').write_text('fixture\n')
            return {'validation_scope':'fixture'}
        with patch('chemflow2.upstream.run_upstream', side_effect=fake):
            task = dict(id='sim',engine='upstream',stage='simulate',outputs=['products/test.csv'],resources=self.limits)
            result = run_task(task,self.case)
            self.assertTrue(result['ok'],result)
            self.assertEqual(set(result['artifacts']),{'products/test.csv','upstream.out'})
