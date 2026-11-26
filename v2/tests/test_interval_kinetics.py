import importlib.util
import math
import unittest

from chemflow2 import mechanism
from chemflow2.analysis_task import _disposition


def reversible():
    training = {'role': 'training', 'initial_total_M': {'A': 1}, 'yield_interval': [0, 1],
                'observation_sha256': 'a' * 64, 'observable': {'product': 'B', 'limiting_pool': 'mass'}}
    return {'operation': 'interval_kinetics', 'model_source_sha256': 'c' * 64,
            'parameter_grid': {'ts': [22., 23.]}, 'temperature_K': 298.15, 'end_h': 1.,
            'species_gibbs_kcal_mol': {'A': 0, 'B': 0}, 'conserved': {'mass': {'A': 1, 'B': 1}},
            'solver_yield_tolerance': 1e-6,
            'reactions': [{'type': 'thermal', 'id': 'iso', 'lhs': {'A': 1}, 'rhs': {'B': 1},
                           'ts_gibbs_kcal_mol': 'ts'}],
            'conditions': {'fit': training,
                           'blind': {'role': 'holdout', 'initial_total_M': {'A': 0.5},
                                     'observable': training['observable']}}}


@unittest.skipUnless(importlib.util.find_spec('scipy'), 'requires SciPy')
class KineticsTests(unittest.TestCase):
    def test_shared_transition_state_matches_exact_closed_kinetics(self):
        data = reversible()
        result = mechanism.evaluate(data)
        self.assertEqual(result['feasible_points'], 2)
        for record in result['points']:
            k = (3600 * mechanism.BOLTZMANN_OVER_PLANCK * 298.15 *
                 math.exp(-record['parameters']['ts'] / (298.15 * mechanism.GAS_CONSTANT_KCAL)))
            self.assertAlmostEqual(record['yield_fractions']['fit'], -math.expm1(-2*k)/2, places=7)
        self.assertEqual(result['holdouts_excluded_from_selection'], ['blind'])
        self.assertLess(result['training_sensitivity_per_kcal_mol']['ts']['fit'][1], 0)

    def test_endpoint_rank_detects_unresolved_parameters(self):
        data = reversible()
        data['species_gibbs_kcal_mol']['B'] = 'offset'
        data['parameter_grid']['offset'] = [-0.1, 0.1]
        result = mechanism.evaluate(data)
        for record in result['points']:
            self.assertEqual(record['local_identifiability']['rank'], 1)
            self.assertEqual(record['local_identifiability']['unresolved_dimensions'], 1)

    def test_partition_and_detailed_balance(self):
        data = reversible()
        data.update(end_h=100, henry_bar_L_mol={'B': 0.0831446261815324 * 298.15})
        for c in data['conditions'].values():
            c.update(liquid_L=1, headspace_L=1)
        # At equilibrium dissolved B equals dissolved A; half of total B is in headspace.
        result = mechanism.evaluate(data)
        self.assertAlmostEqual(result['yield_envelopes']['fit'][0], 2/3, places=7)

    def test_blinding_mass_balance_and_no_scientific_stop_from_empty_grid(self):
        data = reversible()
        data['conditions']['blind']['yield_interval'] = [0, .48]
        with self.assertRaisesRegex(ValueError, 'holdout'):
            mechanism.evaluate(data)
        del data['conditions']['blind']['yield_interval']
        data['conditions']['fit']['yield_interval'] = [.95, 1]
        result = mechanism.evaluate(data)
        self.assertEqual(_disposition('mechanism', result), 'indeterminate')
        self.assertEqual(result['yield_envelopes'], {})
        del data['conserved']['mass']['B']
        with self.assertRaisesRegex(ValueError, 'conservation'):
            mechanism.evaluate(data)

    def test_photon_competition_and_no_metal_arm(self):
        data = reversible()
        data.update(parameter_grid={}, end_h=24, light={'epsilon_L_mol_cm': {'A': 10, 'M': 100}},
                    species_gibbs_kcal_mol={'A': 0, 'B': 0, 'M': 0},
                    reactions=[{'id': 'excite', 'type': 'photon', 'lhs': {'A': 1}, 'rhs': {'B': 1},
                                'quantum_yield': 1}])
        data['conserved']['metal'] = {'M': 1}
        for c in data['conditions'].values():
            c.update(initial_total_M={'A': 1}, path_cm=1, incident_photons_M_h=.02)
        data['conditions']['fit']['initial_total_M']['M'] = 1
        data['conditions']['blind']['absent_pools'] = ['metal']
        data['phase_contrast'] = {'primary': 'blind', 'secondary': 'fit', 'minimum_difference': .2}
        result = mechanism.evaluate(data)
        self.assertTrue(result['points'][0]['secondary_suppressed'])
        data['conditions']['blind']['initial_total_M']['M'] = .1
        with self.assertRaisesRegex(ValueError, 'absent pool'):
            mechanism.evaluate(data)
