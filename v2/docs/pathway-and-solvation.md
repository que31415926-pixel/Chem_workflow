# Reaction paths, explicit ion pairs and solution thermochemistry

## Compile an ORCA path

Run `python3 -m chemflow2.pathway private/path.json path-bundle`, then submit `path-bundle/plan.json` using the chosen backend. Relative XYZ and orbital-guess paths resolve from the input JSON directory. Use a separate config for each reviewed charge, spin surface, coordination, and seed. The following illustrates the input format. Replace the example method, tolerances, atom indices, endpoints and resource values with those of the approved study.

```json
{
  "id": "reviewed-path",
  "source": "original structure record and atom map",
  "approved_by": "structure reviewer",
  "coordination": "reviewed coordination and electronic-state assignment",
  "xyz": "mapped-ts-seed.xyz",
  "charge": 0,
  "multiplicity": 2,
  "imaginary_tolerance_cm1": 20,
  "max_s2_deviation": 0.15,
  "methods": [{
    "id": "surface-a",
    "keywords": "PBE0 D4 def2-TZVP TightSCF CPCM(THF)",
    "frequency_keyword": "Freq",
    "stability": true,
    "engine_version": "6.1.0",
    "guess_files": ["alternative-orbitals.gbw"],
    "requires": ["orca6-native"],
    "resources": {"cores": 4, "memory_mb": 16384, "scratch_mb": 32768}
  }],
  "endpoint_signatures": {
    "reactant": [[0, 1, 1.0, 1.25], [1, 2, 1.7, 3.5]],
    "product": [[0, 1, 1.8, 3.5], [1, 2, 1.0, 1.25]]
  }
}
```

Each branch runs OptTS from the reviewed coordinates. To use ScanTS supply, for instance, `"scan":{"bond_indices":[1,2],"range_angstrom":[2.4,1.4],"points":16}`. The default orbital guess and each additional supplied GBW are independent branches. All branches preserve method, solvent, charge, multiplicity and atom order through frequency and IRC. The IRC reads the verified TS Hessian and generates both trajectories/endpoints. Both endpoints undergo independent optimizations and frequency checks. An optional stability job reads the converged orbitals at each stationary point and checks the printed stability result. Set the S² tolerance only when appropriate for the chosen reference. Specifically designed broken-symmetry states require case-specific diagnostics and expert review.

Endpoint signatures list [atom_i, atom_j, minimum Å, maximum Å] constraints using zero-based indices. Both minima must uniquely match different predeclared signatures to pass the optional endpoint gate. Reviewers should inspect the imaginary displacement in the Hessian and full trajectory, confirm connectivity and spin-state assignments, and feed signed evidence digests into `paired_barriers`. Distance and spin diagnostics alone do not establish the mechanistic claim. Independent higher-level single points and composite thermochemistry can consume the optimized XYZs. [ORCA IRC](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/irc.html) · [ORCA stability analysis](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/stabilityanalysis.html)

## Compile ion-pair conformer searches

Run `python3 -m chemflow2.ion_pairs private/poses.json ion-bundle`. Input contains `id`, `approved_by`, `resources`, `crest_version`, `crest_settings` and `components` keyed by component ID. Each component has `xyz`, `source`, integer `charge` and `role` of `cation`, `anion` or `solvent`. Set CREST options explicitly, e.g. `{"method":"gfn2","solvent_model":"alpb","solvent":"water","energy_window_kcal_mol":6.0}`; freeze the actual solvent and window in the registered study.

Supply one or more entries in `poses`. Each entry has a unique `id`, `family` (`contact` or `solvent_separated`), `multiplicity`, `anchor_atom_indices` keyed by cation/anion, the reviewed `anchor_distance_range_angstrom`, and a list of `placements`. Each placement takes the form `{"component":"fragment-id","quaternion_wxyz":[1,0,0,0],"translation_angstrom":[0,0,0]}`. Rotations are about the fragment file origin, then translation is applied. Indices in this input refer to the local fragment XYZ. Place at least one explicit intervening solvent molecule in a solvent-separated pose. Each generated XYZ feeds a separate CREST task. `seeds.json` records the initial distances and source hashes. Review the resulting ensembles for ion migration, connectivity, coverage and duplicate conformers before admitting them for DFT refinement. Explicit-solvent clusters of different composition require a balanced solvent chemical-potential cycle.

## Calculate a solution cycle

Pass a JSON file to `python3 -m chemflow2 mechanism cycle.json`, or use an analysis task with stage `mechanism`, operation `solution_cycle` and dependency bindings. Set `kind` to `reaction`, `solvation` or `relative_reduction`. Provide `stoichiometry` as signed integer coefficients keyed by species ID (positive for products). The `species` map stores each species' ensemble Gibbs energy and absolute uncertainty bound, e.g.

```json
{
  "composition": {"C": 1, "H": 4},
  "charge": 0,
  "multiplicity": 1,
  "species_id": "methane",
  "gibbs": -40.1,
  "unit": "hartree",
  "absolute_error_bound": 0.001,
  "energy_sha256": "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
  "method": "registered-composite-and-ensemble",
  "solvent": "water",
  "solvent_model": "SMD",
  "standard_state": "1M",
  "temperature_K": 298.15
}
```

Gas-to-solution transfer requires the same species, spin, method, temperature and concentration standard state on both sides. Apply the experiment's gas-to-solution standard-state correction explicitly when constructing the energy records. For relative reduction potentials, provide both the target reduction and `reference: {"id":"Fc+/Fc","stoichiometry":{"reference_ox":-1,"reference_red":1}}`. Include all referenced species in the same map. The evaluator balances nuclei, infers electron counts from the formal charge change, and subtracts reference-couple free energy per electron under matching conditions. Returned bounds propagate the supplied absolute energy bounds conservatively. Inputs must already use the reviewed conformer ensembles and thermal corrections. Do not interpret these computed potentials as absolute-electrode values.

### Split thermochemistry control

For a separate, prospectively defined control, add `single_point` to a reviewed pathway method. Pin `engine_version`, supply reviewed `endpoint_signatures`, and specify the complete electronic method for the independent single point. Each optimized TS or IRC endpoint geometry feeds a frequency calculation and a single point. Both inherit the same starting orbitals and any requested low-level stability prerequisite. An enabled stability calculation checks its own electronic model. Review the high-level electronic state separately. The original path definition, basis, medium and convergence criteria belong in the frozen method record. Example method fields:

```json
{
  "stage_resources": {
    "freq": {"cores": 8, "memory_mb": 48000, "scratch_mb": 120000},
    "sp": {"cores": 4, "memory_mb": 24000, "scratch_mb": 60000}
  },
  "engine_version": "6.1.0",
  "single_point": {
    "keywords": "PBE0 def2-TZVP TightSCF CPCM(water)",
    "requires": ["reviewed-high-sp"],
    "temperature_K": 298.15,
    "pressure_atm": 1.0
  }
}
```

The reservations above illustrate syntax. Measure representative jobs and freeze comparable settings before selecting production resources. Stage overrides replace the entire reservation for the selected stage; other stages inherit the method reservation. `sp` also covers any stability jobs. The scheduler may serialize eligible siblings if capacity is limited.

Each `*-gibbs` task checks its optimization, frequency and single-point evidence, their input deck hashes, common XYZ hash, actual ORCA version, expected saddle or minimum mode count, matching printed thermochemistry temperature and pressure, and the reviewed end structures reached from the IRC. It then adds the ORCA `G-E(el)` term to the single-point electronic energy (hartree). This yields a single-structure composite at the recorded pressure. Basis extrapolation, ensemble averaging, concentration/standard-state conversion, physical mode-direction assessment and the independent paired-barrier gate remain separate operations. Same-geometry sibling eligibility is a scheduling opportunity; report speedup only from an independently measured paired run.

[ORCA thermochemistry conventions](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/thermochemistry.html)
