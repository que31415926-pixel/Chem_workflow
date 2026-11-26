import copy
import unittest

from chemflow2.solution import salt_quantity


class SaltQuantityTests(unittest.TestCase):
    def setUp(self):
        self.review = {"reviewer": "test", "decision": "pass", "items": [{"sha256": "1"*64}],
                       "qualified_methods": {"A": {"method_settings_sha256": "2"*64, "solvent_model": "test",
                                                   "reference_couple": "reference"}}}
        self.energies = {}
        for ident, charge, composition, energy in [
                ("salt",1,{"C":1,"H":1},5), ("cat",1,{"C":1},8), ("neutral",0,{"H":1},2),
                ("reduced",0,{"C":1,"H":1},2), ("ref-ox",1,{"N":1},3), ("ref-red",0,{"N":1},1)]:
            self.energies[ident] = {"charge":charge, "composition":composition, "gibbs":energy,
                  "solvent_model":"test", "solvent":"test-solvent", "temperature_K":300,
                  "method":"A", "standard_state":"1M", "unit":"kcal/mol", "energy_sha256":"a"*64,
                  "absolute_error_bound":0}
        self.calc={"method_settings_sha256":"2"*64, "salt":"test-salt", "solvent":"test-solvent",
                   "salt_reactant":"salt"}

    def run_case(self, quantity, cycle):
        return salt_quantity({"quantity":quantity, "salt":"test-salt", "solvent":"test-solvent",
                              "method_review":self.review, "by_method":{"A":dict(self.calc, cycle=cycle)}})

    def test_independent_dissociation_and_one_electron_referenced_potential(self):
        ionic = {"kind":"reaction", "species":self.energies,
                 "stoichiometry":{"salt":-1,"cat":1,"neutral":1}}
        rxn=self.run_case("ionic_reaction", ionic)
        self.assertEqual(rxn["by_method"]["A"]["delta_g_kcal_mol"],5)
        potential={"kind":"relative_reduction","species":self.energies,
                   "stoichiometry":{"salt":-1,"reduced":1},
                   "reference":{"id":"reference","stoichiometry":{"ref-ox":-1,"ref-red":1}}}
        reduction=self.run_case("one_electron_reduction",potential)
        self.assertAlmostEqual(reduction["by_method"]["A"]["potential_V"],4184/96485.33212)
        with self.assertRaisesRegex(ValueError,"requested quantity"):
            self.run_case("ionic_reaction",potential)
        with self.assertRaisesRegex(ValueError,"reviewed couple"):
            self.run_case("one_electron_reduction",dict(potential,reference=dict(potential["reference"],id="other")))
        self.review["qualified_methods"]["B"]=copy.deepcopy(self.review["qualified_methods"]["A"])
        with self.assertRaisesRegex(ValueError,"every qualified"):
            self.run_case("one_electron_reduction",potential)

    def test_reduced_fragmentation_requires_correct_charge_and_endpoint_qualification(self):
        from test_mechanism import paired
        fixture = paired()
        spec, path = fixture["protocol"]["paths"][0], fixture["paths"][0]
        method = self.review["qualified_methods"].pop("A")
        method["imaginary_tolerance_cm1"] = 0
        self.review["qualified_methods"]["M1"] = method
        calc = dict(self.calc, specification=spec, path=path, oxidized_assembly_charge=1,
                    condition_id="C", salt_reactant=spec["reactant"])
        data = {"quantity":"reduced_CN_barrier", "salt":"test-salt", "solvent":"test-solvent",
                "method_review":self.review, "by_method":{"M1":calc}}
        self.assertEqual(salt_quantity(data)["by_method"]["M1"]["forward_kcal_mol"], 10)
        path["irc_endpoints"] = ["wrong", "endpoint"]
        self.assertEqual(salt_quantity(data)["status"], "indeterminate")
        calc["oxidized_assembly_charge"] = 0
        with self.assertRaisesRegex(ValueError, "oxidation level"):
            salt_quantity(data)
