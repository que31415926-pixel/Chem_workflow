import json
from pathlib import Path
import tempfile
import unittest

from chemflow2 import mechanism
from chemflow2.ion_pairs import prepare
from chemflow2.pathway import check_stationary_observations, compile_pathway, pathway_gibbs, read_thermochemistry
from chemflow2.protocol import load_plan
from chemflow2.worker import Rejected, read_orca_log, validate_deck


class PathwayTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.inp = self.root / "settings.json"
        self.out = self.root / "compiled"
        self.resources = dict(cores=2, memory_mb=1024, scratch_mb=32)
        self.seed = self.root / "seed.xyz"
        self.seed.write_text("2\nsynthetic geometry\nH 0 0 0\nH 0 0 0.74\n")

    def test_branch_artifacts_and_multiplicity(self):
        (self.root / "seed.gbw").write_bytes(b"fixture")
        spec = dict(id="path", xyz="seed.xyz", approved_by="reviewer", source="contract test",
                    charge=0, multiplicity=1, coordination="declared", imaginary_tolerance_cm1=22,
                    endpoint_signatures={"reactant": [[0, 1, 0.4, 1.0]], "product": [[0, 1, 1.3, 2.0]]},
                    methods=[dict(id="m", keywords="HF def2-SVP", resources=self.resources,
                                  stability=True, guess_files=["seed.gbw"])])
        self.inp.write_text(json.dumps(spec))
        compile_pathway(self.inp, self.out)
        tasks = load_plan(self.out / "plan.json")["tasks"]
        self.assertEqual(len(tasks), 22)
        irc = next(t for t in tasks if t["stage"] == "irc")
        self.assertIn("saddle.hess", irc["from_deps"])
        self.assertIn("Hess_Filename", (self.out / irc["input"]).read_text())
        endings = [t for t in tasks if t["stage"] == "opt" and t["id"].startswith("m-default")]
        self.assertTrue(endings[0]["from_deps"]["structure.xyz"].endswith("_IRC_F.xyz"))
        self.assertTrue(endings[1]["from_deps"]["structure.xyz"].endswith("_IRC_B.xyz"))
        (self.out / "structure.xyz").write_bytes(self.seed.read_bytes())
        for task in tasks:
            if task["engine"] == "orca":
                self.assertEqual(validate_deck(self.out / task["input"], task["stage"], self.resources), 2)
        report = mechanism.evaluate(dict(operation="irc_endpoints", expected=["r", "p"], branches={
            k: dict(endpoint_candidates=[v], endpoint_geometry_sha256=d*64) for k,v,d in (("F","p","a"),("B","r","b"))}))
        self.assertEqual(report["status"], "pass")
        spec["multiplicity"] = 2
        self.inp.write_text(json.dumps(spec))
        with self.assertRaisesRegex(ValueError, "electronic state"):
            compile_pathway(self.inp, self.root / "rejected")

    def test_split_control_siblings_and_verified_composite(self):
        import hashlib
        spec = dict(id="path", xyz="seed.xyz", approved_by="reviewer", source="test",
                    charge=0, multiplicity=1, coordination="declared", imaginary_tolerance_cm1=22,
                    endpoint_signatures={"r": [[0, 1, 0.4, 1.0]], "p": [[0, 1, 1.3, 2.0]]},
                    methods=[dict(id="m", keywords="HF def2-SVP", resources=self.resources,
                                  stage_resources={"sp": dict(cores=1, memory_mb=1024, scratch_mb=20)},
                                  engine_version="6.1.0", single_point=dict(keywords="HF def2-TZVP",
                                      requires=["sp-model"], temperature_K=298.15, pressure_atm=1.0))])
        self.inp.write_text(json.dumps(spec))
        compile_pathway(self.inp, self.out)
        tasks = {t["id"]: t for t in load_plan(self.out/"plan.json")["tasks"]}
        for point in ("ts", "end-F", "end-B"):
            parent = tasks["m-default-" + point]
            freq, sp, join = (tasks[parent["id"] + "-" + suffix] for suffix in ("freq", "sp", "gibbs"))
            self.assertEqual(sp["from_deps"], freq["from_deps"])
            self.assertEqual(sp["deps"], [parent["id"]])
            self.assertEqual(sp["requires"], ["sp-model"])
            self.assertEqual(sp["resources"]["cores"], 1)
            self.assertEqual(freq["resources"]["cores"], 2)
            self.assertIn(sp["id"], join["deps"])
            self.assertIn(freq["id"], join["deps"])
            self.assertIn("m-default-match-ends", join["deps"])
            request = json.loads((self.out/join["input"]).read_bytes())
            geometry = "a" * 64
            records = {"endpoints": {"status": "pass", "geometry_sha256": {"F": geometry, "B": geometry}}}
            for role, task in (("optimization", parent), ("frequency", freq), ("single", sp)):
                sha = hashlib.sha256((self.out/task["input"]).read_bytes()).hexdigest()
                records[role] = {"input_sha256": {task["input"]: sha, "structure.xyz": geometry},
                                 "observations": {"normal_termination": True, "engine_version": "6.1.0"}}
            records["optimization"]["observations"]["optimization_converged"] = True
            records["optimization"]["geometry_sha256"] = {parent["id"]+".xyz": geometry}
            records["frequency"]["observations"].update(imaginary_modes_below_tolerance=int(point=="ts"),
                thermochemistry=dict(temperature_K=298.15, pressure_atm=1, correction_eh=-.02))
            records["single"]["observations"]["electronic_energy_eh"] = -4
            request["records"] = records
            self.assertAlmostEqual(mechanism.evaluate(request)["gibbs_eh"], -4.02)
            records["single"]["input_sha256"]["structure.xyz"] = "b"*64
            with self.assertRaisesRegex(ValueError, "optimized geometry"):
                pathway_gibbs(request)

    def test_thermal_output_uses_one_complete_condition(self):
        path = self.root/"orca.out"
        block = ("THERMOCHEMISTRY AT 298.15K\n Temperature ... 298.15 K\n"
                 "Pressure ... 1.00 atm\n G-E(el) ... -0.020000 Eh -12.55 kcal/mol\n")
        condition = dict(temperature_K=298.15, pressure_atm=1.0)
        path.write_text(block + block.replace("298.15", "300.00").replace("-0.020000","0.010000"))
        self.assertEqual(read_thermochemistry(path, condition)["correction_eh"], -.02)
        path.write_text(block + block)
        with self.assertRaisesRegex(Rejected, "ambiguous"):
            read_thermochemistry(path, condition)
        path.write_text(block.replace("1.00 atm", "1.01 atm"))
        with self.assertRaisesRegex(Rejected, "missing"):
            read_thermochemistry(path, condition)

    def test_numerical_gates(self):
        log = self.root / "orca.out"
        log.write_text("Expectation value of <S**2> : 0.7554\nStability Analysis indicates a STABLE HF/KS wave function.\n")
        obs = read_orca_log(log)
        task = dict(state_diagnostics=dict(require_stability=True, max_s2_deviation=0.02, multiplicity=2))
        check_stationary_observations(task, obs, self.root)
        obs["spin_squared"] = 2.
        with self.assertRaises(Rejected):
            check_stationary_observations(task, obs, self.root)
        task = dict(stage="freq", stationary_point=dict(imaginary_count=1, tolerance_cm1=25, xyz="seed.xyz"))
        obs = {"frequencies_cm-1":[0,0,0,0,0,-800]}
        check_stationary_observations(task, obs, self.root)
        obs["frequencies_cm-1"][0] = -40
        with self.assertRaisesRegex(Rejected, "imaginary-mode"):
            check_stationary_observations(task, obs, self.root)

    def test_ion_pair_assembly_and_explicit_solvent(self):
        components = {}
        for elem, role, charge in (("Li", "cation", 1), ("F", "anion", -1), ("He", "solvent", 0)):
            (self.root / (elem+".xyz")).write_text(f"1\nfixture\n{elem} 0 0 0\n")
            components[elem] = dict(role=role, charge=charge, xyz=elem+".xyz", source="fixture")
        pose = dict(id="placed", family="solvent_separated", multiplicity=1,
                    anchor_atom_indices=dict(cation=0, anion=0), anchor_distance_range_angstrom=[4, 5],
                    placements=[dict(component=elem, quaternion_wxyz=[1,0,0,0],
                                     translation_angstrom=[dx,0,0]) for elem,dx in (("Li",0),("F",4.5),("He",2.2))])
        spec = dict(id="assembly", approved_by="test", components=components, poses=[pose],
                    crest_version="3.0.2", resources=self.resources,
                    crest_settings=dict(method="gfn2", solvent_model="alpb", solvent="water", energy_window_kcal_mol=5))
        self.inp.write_text(json.dumps(spec))
        prepare(self.inp, self.out)
        task = load_plan(self.out / "plan.json")["tasks"][0]
        self.assertEqual(task["settings"]["charge"], 0)
        self.assertEqual(task["settings"]["uhf"], 0)
        pose["placements"][-1]["translation_angstrom"] = [9,0,0]
        self.inp.write_text(json.dumps(spec))
        with self.assertRaisesRegex(ValueError, "intervene"):
            prepare(self.inp, self.root / "rejected")

    def test_relative_potential_sign_balance_and_error(self):
        species = {}
        for key, q, energy in (("ox",1,0),("red",0,-25),("std_ox",1,0),("std_red",0,-2)):
            species[key] = dict(composition={"H":2}, charge=q, gibbs=energy, absolute_error_bound=.1,
                                unit="kcal/mol", method="M", solvent="water", solvent_model="SMD",
                                standard_state="1M", temperature_K=298, energy_sha256="a"*64)
        op = dict(operation="solution_cycle", kind="relative_reduction", species=species,
                  stoichiometry={"ox":-1, "red":1}, reference=dict(id="Fc+/Fc",stoichiometry={"std_ox":-1, "std_red":1}))
        result = mechanism.evaluate(op)
        self.assertAlmostEqual(result["potential_V"],23*4184/96485.33212)
        self.assertLess(result["potential_interval_V"][0], result["potential_V"])
        species["std_red"]["solvent"] = "DMA"
        with self.assertRaisesRegex(ValueError, "different"):
            mechanism.evaluate(op)
