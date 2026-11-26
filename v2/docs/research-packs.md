# Private research packs

Run from the engine root. A study is a user-supplied private directory. Its protocol lists named steps, prerequisites, result identifiers, acceptance criteria, failure routes and groups of work. The Cartesian product of each group's axes defines the individual cells. The private file inventory binds the protocol to its frozen sources. Paths below refer to the local workspace.

```sh
python -m chemflow2 research inventory private/study/protocol.json
python -m chemflow2 research compile private/study/protocol.json private/intake.json submission
python -m chemflow2 validate submission/plan.json
```

Use `--through STEP` to compile an early step and all of its declared prerequisites with the reviews and inputs available for that phase. Review original sources privately. If the owner has to choose between conflicting routes in the source protocol, record the exact clause, alternatives, scope, decision and reviewer in the source review. A missing coordinate, primary source, held-out measurement or electronic state is an input blocker. Only a reviewed observation or a completed gate can support a scientific stop.

The intake is a JSON object with `protocol_sha256` from inventory and a `cells` object. Keys are the complete cell identifiers emitted by inventory. Each record supplies:

* `source` and `source_sha256` locate and freeze a private input relative to the intake.
* `exports` maps the required scientific product names to `task-id/declared-output`. Add named outputs for additional spin states, geometries, basis components and alternative methods. All underlying task outputs remain available.
* `bindings` binds a prerequisite cell's named product to a task-local file. For example, `{"task":"evaluate","as":"left.json","cell":"A/measure/x","output":"value"}` connects the producer's current-run bytes. Analysis requests reference a member with `{"$artifact":"left.json","path":["observations","electronic_energy_eh"]}`. All downstream analysis bindings must occur in the recipe.

| Group kind | Source and runtime |
|---|---|
| review | Human handoff with exact `scope` equal to the cell ID, reviewer, rationale, decision (`pass`, `scientific_stop`, `awaiting_input`, `disputed`) and `items`. Each item supplies the local original file, SHA-256, and exact page/table or run locator. The worker validates the original bytes, records the decision and blocks dependent work on an adverse or incomplete review. |
| study | Reviewed Opt/Freq/SP matrix for compile-study. Use paper_inputs, reviewed XYZ, all electronic states and included conformers. Set the method's `functional` to its ORCA keyword for any pack validating a minimum number of DFT methods. Per-task `values.json` carries parsed observables and the ORCA log digest. |
| ion_pairs | Poses, fragment XYZ, charges, solvent molecules, transformations and distance windows consumed by the ion-pair generator. The compiled workers run the CREST searches. Expert review of the resulting conformers precedes downstream production. |
| pathway | Reviewed stationary-point guess, method and spin settings consumed by the pathway compiler. Each configuration may cover multiple methods/guesses; additional independently reviewed sites and states may be supplied as distinct axes/levels in the private protocol before freezing intake. |
| analysis | A science or mechanism request. The protocol fixes the operation, required gate and benchmark rules. The compiled request contains the numerical thresholds. |
| plan | Validated native plan with all referenced decks/recipes and required stages. Use for original SI procedures and multi-task calculations that need method-specific operators. Engineering probe tasks are rejected. |

After compiling, `research-map.json` maps every cell, scientific output, original task ID and published task ID. Cells can use the same local filename for different molecules. The compiler stores the correct immutable bytes for each task. Computational Results remain pending until run outputs, benchmark qualifications, expert interpretation and required independent checks exist. A protocol and input hash provide local integrity and traceability; use your project's authenticated external registration for prospective freeze dates or independent review.

`research assemble PROTOCOL PRIVATE_DIR PRIVATE_DIR/intake.json` inventories actual inputs and computes their digests. Inside PRIVATE_DIR, use one directory per full cell identifier. Filenames by type are review.json, request.json, study.json, poses.json, pathway.json or plan.json. Place bindings and exports in links.json next to the source. For a review with a single report, the exports are generated automatically. The assemble command lists missing cells and leaves them unregistered. The compile command validates coverage again.

For a pathway requiring a multi-state registry, pathway.json declares `registry`, its `registry_sha256` and a `cases` map. Each case entry contains the reviewed configuration's relative `config` path and `sha256`. The registry supplies the reviewer and a cases map keyed by the identical case IDs with coordination, charge and multiplicity. Each configuration has both endpoint signatures and all required electronic methods. A missing site or spin case from the frozen registry prevents compilation. Native product references for this group begin with the case ID followed by a period and the local task ID.

A private pack may ship frozen ORCA method sets for a molecular matrix. `research prepare-quantum PROTOCOL APPROVED_REGISTRY DEST` generates one ordinary study per required species and solvent model. The registry maps each species ID to `states` and `conformers` as in compile-study. Every conformer supplies its approved XYZ path, SHA-256, source and reviewer. Missing species, method components or changed XYZ files prevent generation. This command prepares electronic-structure inputs; benchmark records, method qualification, high-level analyses and evidence bindings still belong to the research intake.

Each compiled molecular cell collects all successful ORCA observations into matrix.json. Its rows identify species, state, charge, spin multiplicity, conformer, method, frozen template, input hashes, produced geometry hashes and source log hash. For generated method matrices, prepare-quantum writes links.json alongside study.json so assemble can register the cell immediately.

A `salt_quantity` analysis binds `method-review.json` from the method qualification review and `computed.json` from the matching state or pathway cell. The method review supplies a `qualified_methods` map keyed by the qualified electronic method, including its settings hash and solvent model. For a relative potential it also specifies the reference couple. For a TS barrier it specifies the imaginary-mode tolerance. Each salt/solvent needs calculations for exactly these approved methods. Independent outputs report ionic dissociation ΔG, ionic activation ΔG‡, one-electron reduction potential relative to the reviewed reference, and post-reduction C–N cleavage ΔG‡. Failure to qualify the modes, branches or electronic state marks the barrier indeterminate. Downstream anchor scoring requires the training-only threshold fit and independent held-out evidence.

A private intake can also freeze an independent audit plan. Add `evidence_plan` with the plan's relative `path` and `sha256` to intake.json before compiling. The file declares `purpose: paper_inputs`, the `paper` matching the study name, named independent `reviewers`, and `obligations` mapping claim or step IDs to `jobs` (published task-ID/output, available in research-map.json) and `sources` (source ID to evidence category). The compiler verifies each job reference and attaches the hashed plan. Independent review and original evidence records are supplied separately after execution. A scoped private audit reports only that paper's obligations; audit each paper to evaluate both.

Prospective predictions refer to the SHA-256 of the original research protocol. After execution the evidence ledger refers to the compiled protocol SHA-256. The independent audit checks both and matches the externally registered prediction to a verified freeze task from the same run.
