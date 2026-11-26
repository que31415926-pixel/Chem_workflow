"""Numerical fixtures are synthetic, unrelated to the paper cases."""
import math
import unittest

from chemflow2 import composite, science


class CompositeTests(unittest.TestCase):
    def setUp(self):
        self.common = {"unit": "Eh", "geometry_sha256": "1"*64, "medium_sha256": "2"*64,
                       "method_settings_sha256": "3"*64, "charge": 0, "multiplicity": 1,
                       "state_id": "synthetic-singlet"}
        self.model = {"model_id": "fake-model", "parameter_source": "analytic synthetic fixture",
                      "basis_family": "def2-VPP", "scf_model": "exp_sqrt_cardinal", "scf_alpha": 2.5,
                      "correlation_model": "inverse_power", "correlation_beta": 3.0}
        self.cbs = {"protocol_id": "analytical-only", "operation": "two_point_cbs", "model": self.model,
                    "frozen_model_sha256": composite.model_digest(self.model)}
        for key, x, basis, hash_char in [("tz", 3, "def2-TZVPP", "a"), ("qz", 4, "def2-QZVPP", "b")]:
            scf = -10 + 2*math.exp(-2.5*math.sqrt(x))
            corr = -0.25 + 0.75*x**-3
            self.cbs[key] = {**self.common, "task_id": key, "artifact_sha256": hash_char*64,
                             "unit": "Eh", "quantity": "electronic_components", "basis": basis,
                             "cardinal": x, "scf_eh": scf, "correlation_eh": corr, "electronic_eh": scf+corr}
        self.thermo = {"operation": "composite_thermochemistry", "protocol_id": "illustration",
                       "high_level": {**self.common, "quantity": "electronic_energy", "electronic_eh": -10.25,
                                      "task_id": "high", "artifact_sha256": "c"*64},
                       "frequency": {**self.common, "quantity": "qrrho_gibbs_and_electronic",
                                     "method_settings_sha256": "5"*64, "electronic_eh": -10.1,
                                     "qrrho_gibbs_eh": -10.08, "temperature_k": 298.15,
                                     "reported_standard_state": "1_atm", "qrrho_model_id": "toy-qrrho",
                                     "reviewed_by": "unit-test", "task_id": "low", "artifact_sha256": "d"*64},
                       "standard_state_policy": {"source": "1_atm", "target": "1_M",
                                                 "mode": "ideal_gas_pressure_to_1M"}}

    def test_extrapolation_matches_analytical_limit_via_dispatcher(self):
        result = science.evaluate(self.cbs)
        self.assertEqual(result["status"], "pass")
        self.assertAlmostEqual(result["scf_cbs_eh"], -10)
        self.assertAlmostEqual(result["correlation_cbs_eh"], -0.25)
        self.assertAlmostEqual(result["electronic_eh"], -10.25)
        self.assertEqual(result["source_records"]["tz"]["task_id"], "tz")

    def test_rejects_wrong_partition_state_family_and_reused_source(self):
        for changes in ({"correlation_eh": 0}, {"basis": "def2-QZVP"}, {"cardinal": 4.1},
                        {"geometry_sha256": "f"*64}, {"charge": 1}, {"state_id": "other"},
                        {"medium_sha256": "e"*64}, {"method_settings_sha256": "a"*64},
                        {"unit": "kcal/mol"}, {"task_id": "tz"}, {"artifact_sha256": "a"*64}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                science.evaluate({**self.cbs, "qz": {**self.cbs["qz"], **changes}})

    def test_exponents_and_model_must_be_frozen(self):
        for changes in ({"scf_alpha": 7}, {"correlation_beta": 0}, {"scf_model": "other"}):
            modified = {**self.model, **changes}
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                science.evaluate({**self.cbs, "model": modified})
        self.cbs["model"].pop("scf_alpha")
        self.cbs["frozen_model_sha256"] = composite.model_digest(self.model)
        with self.assertRaisesRegex(ValueError, "scf_alpha"):
            science.evaluate(self.cbs)

    def test_composite_and_standard_state_direction(self):
        result = science.evaluate(self.thermo)
        self.assertAlmostEqual(result["thermal_qrrho_correction_eh"], 0.02)
        self.assertAlmostEqual(result["per_species_standard_state_correction_eh"]*627.509474,
                               1.894, delta=0.003)
        self.assertAlmostEqual(result["gibbs_eh"], -10.23 + result["per_species_standard_state_correction_eh"])
        self.thermo["frequency"]["reported_standard_state"] = "1_M"
        with self.assertRaisesRegex(ValueError, "conflicts"):
            science.evaluate(self.thermo)
        self.thermo["standard_state_policy"] = {"source": "1_M", "target": "1_M", "mode": "as_reported"}
        self.assertEqual(science.evaluate(self.thermo)["per_species_standard_state_correction_eh"], 0)

    def test_invalid_thermal_sources(self):
        for changes in ({"quantity": "rrho_uncorrected"}, {"qrrho_model_id": ""},
                        {"reviewed_by": ""}, {"geometry_sha256": "0"*64}, {"temperature_k": 0}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                science.evaluate({**self.thermo, "frequency": {**self.thermo["frequency"], **changes}})
        self.thermo["standard_state_policy"]["extra_correction"] = 1
        with self.assertRaisesRegex(ValueError, "only"):
            science.evaluate(self.thermo)
