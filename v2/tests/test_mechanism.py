import copy
import unittest
from datetime import datetime, timedelta, timezone

from chemflow2 import mechanism


def digest(n):
    return format(n, '064x')


def paired():
    paths, specs = [], []
    for i, key in enumerate(('A', 'B')):
        spec = dict(id=key, reactant=f'{key}-R', product=f'{key}-P', method_id='M1', state_id='S', multiplicity=1)
        specs.append(spec)
        row = dict(id=key, irc_endpoints=[spec['product'], spec['reactant']], evidence={
            'frequency': digest(1), 'irc_forward': digest(2), 'irc_backward': digest(3)},
            review={'reviewer': 'reviewer', **{f: {'accepted': True, 'artifact_sha256': digest(j)}
                    for j, f in enumerate(('mode_assignment', 'endpoint_structures', 'electronic_state'), 10)}})
        for kind, energy in [('reactant', 0), ('ts', 10 + 4 * i), ('product', -3)]:
            row[kind] = dict(composition={'H': 2, 'O': 1}, charge=0, multiplicity=1, state_id='S', method_id='M1',
                             condition_id='C', gibbs=energy, gibbs_unit='kcal/mol', energy_sha256=digest(41 + i),
                             species=spec.get(kind, 'saddle'), geometry_type='nonlinear',
                             vibrational_modes_cm1=[-900 if kind == 'ts' else 510, 232, 690])
        paths.append(row)
    protocol = dict(condition_id='C', temperature_K=323.15, imaginary_tolerance_cm1=0, paths=specs,
                    comparisons=[dict(id='A-minus-B', left='A', right='B', expected_sign=-1,
                                      indeterminate_band_kcal_mol=0.5)])
    return dict(operation='paired_barriers', protocol=protocol, protocol_sha256=mechanism._hash(protocol), paths=paths)


def prediction():
    return dict(author='analysis', observable='fraction', unit='1', condition_id='run condition',
                experimental_procedure='measurement-v1', falsification_rule='outside interval',
                protocol_sha256=digest(4), model_sha256=digest(5),
                training=[dict(sample_id='train', group_id='group-1', observation_sha256=digest(22))],
                predictions=[dict(sample_id='test', group_id='group-2', lower=0.4, upper=0.6)])


def scored():
    now = datetime.now(timezone.utc)
    record = mechanism.evaluate(dict(operation='freeze_prediction', specification=prediction()))['frozen']
    record['frozen_at'] = (now - timedelta(days=2)).isoformat()
    obs = dict(sample_id='test', group_id='group-2', unit='1', observable='fraction', condition_id='run condition',
               observed_at=(now - timedelta(days=1)).isoformat(), source_type='experiment', observer='lab',
               observation_sha256=digest(33), value=0.47, absolute_uncertainty=0.01)
    return dict(operation='score_holdout', frozen=record, sha256=mechanism._hash(record), observations=[obs])


class MechanismTests(unittest.TestCase):
    def test_reverse_irc_and_paired_barrier(self):
        data = paired()
        data['paths'][1]['ts'].update(gibbs=14 * 4.184, gibbs_unit='kJ/mol')
        result = mechanism.evaluate(data)
        self.assertTrue(result['all_directions_agree'])
        self.assertEqual(result['paths']['A']['reverse_kcal_mol'], 13)
        self.assertEqual(result['comparisons']['A-minus-B']['left_minus_right_kcal_mol'], -4)

    def test_missing_path_unqualified(self):
        data = paired()
        data['paths'].pop()
        report = mechanism.evaluate(data)
        self.assertEqual(report['paths']['B']['status'], 'missing')
        self.assertFalse(report['all_directions_agree'])

    def test_close_or_inverted_barriers(self):
        data = paired()
        for energy, answer in [(10.5, 'within_band'), (8, 'disagreement')]:
            data['paths'][1]['ts']['gibbs'] = energy
            result = mechanism.evaluate(data)
            self.assertEqual(result['comparisons']['A-minus-B']['status'], answer)

    def test_tamper_detected(self):
        data = paired()
        data['protocol']['comparisons'][0]['expected_sign'] = 1
        with self.assertRaisesRegex(ValueError, 'commitment'):
            mechanism.evaluate(data)

    def test_invalid_ts_surface_mode_and_unreviewed(self):
        base = paired()
        for field, value in [('multiplicity', 3), ('vibrational_modes_cm1', [-800, -60, 300]),
                              ('condition_id', 'different solvent'), ('vibrational_modes_cm1', [-800])]:
            data = copy.deepcopy(base)
            data['paths'][0]['ts'][field] = value
            self.assertEqual(mechanism.evaluate(data)['paths']['A']['status'], 'unqualified')
        data = copy.deepcopy(base)
        data['paths'][0]['review']['electronic_state']['accepted'] = False
        self.assertFalse(mechanism.evaluate(data)['all_directions_agree'])
        base['paths'][0]['irc_endpoints'][0] = 'wrong'
        self.assertFalse(mechanism.evaluate(base)['all_directions_agree'])

    def test_group_split(self):
        data = prediction()
        data['predictions'][0]['group_id'] = 'group-1'
        with self.assertRaisesRegex(ValueError, 'overlapping'):
            mechanism.evaluate(dict(operation='freeze_prediction', specification=data))

    def test_freeze_then_score(self):
        frozen = mechanism.evaluate(dict(operation='freeze_prediction', specification=prediction()))
        self.assertEqual(frozen['sha256'], mechanism._hash(frozen['frozen']))
        data = scored()
        for val, ans in [(0.47, 'inside'), (0.6, 'overlap'), (0.8, 'outside')]:
            data['observations'][0]['value'] = val
            self.assertEqual(mechanism.evaluate(data)['targets']['test']['status'], ans)
        data['observations'] = []
        self.assertFalse(mechanism.evaluate(data)['all_observed'])

    def test_reject_circular_or_incompatible_holdout(self):
        for field, value in [('source_type', 'simulation'), ('observation_sha256', digest(22)),
                             ('observed_at', '2000-01-01T00:00:00Z'), ('unit', '%')]:
            data = scored()
            data['observations'][0][field] = value
            with self.assertRaises(ValueError):
                mechanism.evaluate(data)
        data = scored()
        data['frozen']['specification']['predictions'][0]['upper'] = 1
        with self.assertRaisesRegex(ValueError, 'changed'):
            mechanism.evaluate(data)

    def test_archived_holdout_uses_release_event_instead_of_measurement_date(self):
        data = scored()
        spec = data['frozen']['specification']
        spec.update(holdout_design='archived', custodian='archive owner')
        source = 'citation and exact source table location'
        spec['predictions'][0].update(archived_source_sha256=digest(33), source_locator=source)
        data['sha256'] = mechanism._hash(data['frozen'])
        obs = data['observations'][0]
        obs.update(unmasked_at=obs['observed_at'], observed_at='2004-01-01T00:00:00Z',
                   source_locator=source, released_by='archive owner', unmask_receipt_sha256=digest(48))
        self.assertEqual(mechanism.evaluate(data)['targets']['test']['status'], 'inside')
        obs['released_by'] = 'analyst'
        with self.assertRaisesRegex(ValueError, 'custodian'):
            mechanism.evaluate(data)

    def test_eyring(self):
        data = dict(operation='eyring', model='unimolecular_TST', gibbs_unit='kcal/mol', activation_gibbs=15.0,
                    standard_state='1 M', temperature_K=300, transmission_coefficient=1,
                    barrier_evidence_sha256=digest(1))
        k = mechanism.evaluate(data)['k_per_s']
        self.assertTrue(72.5 < k < 74.5)
        data['model'] = 'bimolecular'
        with self.assertRaises(ValueError):
            mechanism.evaluate(data)
