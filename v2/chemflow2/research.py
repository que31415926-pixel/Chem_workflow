"""Compile reviewed research matrices into ordinary ChemFlow dependency plans."""
import argparse
import csv
from collections import defaultdict
from itertools import product
import json
from pathlib import Path
import os
import shutil
import tempfile

from .protocol import (_digest, _identifier, _order, _relative, _static_file,
                       _string, _unique_keys, load_plan, static_sources)
from .study import compile_study


def read(path):
    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_unique_keys)


def write(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8")


def asset(root, name):
    return _static_file(root, _relative(name, "research input"), "research input")


def inventory(protocol):
    pack = read(protocol)
    if pack.get("schema") != 1:
        raise ValueError("unsupported research schema")
    for filename, checksum in pack.get("source_sha256", {}).items():
        if _digest(asset(Path(protocol).resolve().parent, filename)) != checksum:
            raise ValueError(f"frozen protocol source changed: {filename}")
    steps = pack["steps"]
    if len({s["id"] for s in steps}) != len(steps):
        raise ValueError("duplicate study step")
    _order(steps)
    rows = {}
    for step in steps:
        step_id = _identifier(step["id"], "step")
        for group in step["groups"]:
            group_id = _identifier(group["id"], "group")
            axes = group.get("axes", {})
            for key, values in axes.items():
                _identifier(key, "axis name")
                if not isinstance(values, list) or not values or len(set(values)) != len(values):
                    raise ValueError("axis must have unique nonempty levels")
                for value in values:
                    _identifier(value, "axis level")
            for values in product(*axes.values()):
                cell = "/".join([step_id, group_id, *values])
                if cell in rows:
                    raise ValueError("overlapping or duplicate research cell")
                rows[cell] = {"step": step_id, "group": group_id, "kind": group["kind"],
                              "axes": dict(zip(axes, values)), "outputs": group["outputs"],
                              "requires_stages": group.get("requires_stages", []),
                              "parameters": group.get("parameters", {}), "required_bindings": group.get("required_bindings", {}),
                              "method_sets": group.get("method_sets"), "operation": group.get("operation"), "method_contract": group.get("method_contract", {}),
                              "gate": group.get("gate"), "gate_spec": group.get("gate_spec"), "targets": group.get("targets", []),
                              "needs": group.get("needs", []), "claims": step["claims"]}
    if not rows:
        raise ValueError("empty research pack")
    return pack, rows


def run_review(task, folder):
    """Accept an externally reviewed decision, verifying each cited staged artifact."""
    folder = Path(folder).resolve()
    form = read(asset(folder, task["input"]))
    decision = form["decision"]
    if decision not in {"pass", "scientific_stop", "awaiting_input", "disputed"}:
        raise ValueError("unknown review decision")
    if not form["items"]:
        raise ValueError("review requires primary source artifacts")
    evidence = []
    for item in form["items"]:
        path = asset(folder, item["file"])
        digest = _digest(path)
        if digest != item["sha256"]:
            raise ValueError(f"review source digest mismatch: {item['file']}")
        evidence.append({"sha256": digest, "locator": _string(item["locator"], "source locator"),
                         "file": item["file"]})
    for field in ("reviewer", "rationale", "scope"):
        _string(form[field], field)
    status = {"pass": "pass", "scientific_stop": "fail",
              "awaiting_input": "missing", "disputed": "indeterminate"}[decision]
    result = {**form, "items": evidence, "scientific_gate": status, "type": "external_review"}
    if len(task["outputs"]) != 1:
        raise ValueError("one review output required")
    output = folder / _relative(task["outputs"][0], "review result")
    with output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, allow_nan=False)
    return result, None if status == "pass" else "scientific_" + status


def check_study_design(row, config, source):
    contract = row["method_contract"]
    if "species" in row["axes"]:
        if {m["id"] for m in config["species"]} != {row["axes"]["species"]}:
            raise ValueError("species in study does not match the matrix cell")
    keywords = []
    for method in config["methods"]:
        text = asset(source.parent, method["template"]).read_text().casefold()
        text = "\n".join(line.split("#", 1)[0] for line in text.splitlines())
        keywords.append(text)
    for terms in contract.get("required_decks", []):
        if not any(all(term.casefold() in deck for term in terms) for deck in keywords):
            raise ValueError(f"missing mandatory method or basis: {terms}")
    minimum = contract.get("min_dft_functionals", 0)
    if minimum:
        frozen = set()
        for method, deck in zip(config["methods"], keywords):
            functional = method.get("functional")
            if functional:
                token = _string(functional, "functional").casefold()
                terms = {word for line in deck.splitlines() if line.strip().startswith("!")
                         for word in line.strip()[1:].split()}
                if token not in terms or token in {"sp", "tightscf", "tightpno", "freq", "opt", "tightopt"}:
                    raise ValueError("functional must match the electronic method in the ORCA keyword line")
                frozen.add(token)
        if len(frozen) < minimum:
            raise ValueError("insufficient distinct reviewed DFT functionals")
    model = row["axes"].get("model")
    if model == "SMD" and not all("smd" in k for k in keywords):
        raise ValueError("matrix expects SMD for all solution jobs")
    if model == "CPCM" and not all("cpcm" in k and "smd" not in k for k in keywords):
        raise ValueError("matrix expects CPCM for all solution jobs")


def compile_pathway_cases(source, work, contract):
    """Expand the reviewed species/spin/site list into disjoint ORCA pathways."""
    from .pathway import compile_pathway
    config = read(source)
    if "cases" not in config:
        if contract.get("require_case_registry"):
            raise ValueError("registered site/state pathway matrix required")
        compile_pathway(source, work / "bundle")
        return
    inventory_file = asset(source.parent, config["registry"])
    if _digest(inventory_file) != config["registry_sha256"]:
        raise ValueError("reviewed pathway registry changed")
    registry = read(inventory_file)
    _string(registry["reviewer"], "pathway registry reviewer")
    expected = registry["cases"]
    if not isinstance(expected, dict) or not expected or set(expected) != set(config["cases"]):
        raise ValueError("pathway cases must cover the reviewed registry exactly")
    out = work / "bundle"
    out.mkdir()
    attachments = ["pathway-registry.json", "pathway-cases.json"]
    shutil.copyfile(inventory_file, out / attachments[0])
    shutil.copyfile(source, out / attachments[1])
    tasks = []
    for key, entry in config["cases"].items():
        _identifier(key, "pathway case")
        spec = asset(source.parent, entry["config"])
        if _digest(spec) != entry["sha256"]:
            raise ValueError("pathway input changed")
        case = read(spec)
        for field in ("coordination", "charge", "multiplicity"):
            if case[field] != expected[key][field]:
                raise ValueError(f"{key}: pathway does not match the approved {field}")
        if not case.get("endpoint_signatures"):
            raise ValueError("registered pathways require two reviewed endpoint identities")
        if len({m["keywords"] for m in case["methods"]}) < contract.get("min_pathway_methods", 1):
            raise ValueError("required independent pathway methods missing")
        local = out / key
        compile_pathway(spec, local)
        plan = load_plan(local / "plan.json")
        attachments += [key + "/" + item for item in plan.get("attachments", [])]
        for task in plan["tasks"]:
            task["static_sources"] = {name: key + "/" + filename
                                      for name, filename in static_sources(task).items()}
            task["id"] = key + "." + task["id"]
            task["deps"] = [key + "." + dep for dep in task["deps"]]
            task["from_deps"] = {name: key + "." + ref for name, ref in task["from_deps"].items()}
            tasks.append(task)
    write(out / "plan.json", {"name": "registered pathways", "tasks": tasks, "attachments": attachments})


def collect_matrix(data):
    cells = data["rows"]
    if not isinstance(cells, list) or not cells:
        raise ValueError("empty electronic-structure results matrix")
    seen = set()
    for record in cells:
        key = tuple(record[x] for x in ("species", "state", "conformer", "method"))
        value = record["result"]
        if key in seen or not value["observations"].get("normal_termination"):
            raise ValueError("duplicate matrix record or unsuccessful electronic calculation")
        _identifier(value["task_id"], "ORCA task")
        seen.add(key)
    return {"operation": "collect_matrix", "status": "pass", "rows": cells}


def add_matrix_collector(bundle):
    with (bundle / "claims-map.csv").open(newline="", encoding="utf-8") as stream:
        entries = list(csv.DictReader(stream))
    plan = read(bundle / "plan.json")
    bindings, records = {}, []
    for row in entries:
        key = row["task_id"]
        local = key + ".json"
        bindings[local] = key + "/values.json"
        records.append({field: row[field] for field in
                        ("species", "state", "charge", "multiplicity", "conformer", "method",
                         "stage", "seed_sha256", "template_sha256")})
        records[-1]["result"] = {"$artifact": local, "path": []}
    write(bundle / "matrix-request.json", {"operation": "collect_matrix", "rows": records})
    plan["tasks"].append({"id": "matrix", "engine": "analysis", "stage": "science",
                          "input": "matrix-request.json", "deps": sorted({r["task_id"] for r in entries}),
                          "outputs": ["matrix.json"], "from_deps": bindings,
                          "resources": {"cores": 1, "memory_mb": 512, "scratch_mb": 64}})
    write(bundle / "plan.json", plan)


def _local_plan(row, record, root, work):
    kind = row["kind"]
    source = asset(root, record["source"])
    if record["source_sha256"] != _digest(source):
        raise ValueError(f"research input changed: {record['source']}")
    if kind == "study":
        cfg = read(source)
        check_study_design(row, cfg, source)
        if cfg.get("scope") != "paper_inputs":
            raise ValueError("study matrices require approved paper_inputs")
        compile_study(source, work / "bundle")
        add_matrix_collector(work / "bundle")
        source = work / "bundle/plan.json"
    elif kind == "ion_pairs":
        from .ion_pairs import prepare
        prepare(source, work / "bundle")
        source = work / "bundle/plan.json"
    elif kind == "pathway":
        compile_pathway_cases(source, work, row["method_contract"])
        source = work / "bundle/plan.json"
    elif kind in {"analysis", "review"}:
        data = read(source)
        content = _digest(source)[:24]
        name = "request-" + content + ".json"
        shutil.copyfile(source, work / name)
        files = []
        if kind == "review":
            if data.get("scope") != record["cell"]:
                raise ValueError(f"review must explicitly cover {record['cell']}")
            for item in data["items"]:
                path = asset(source.parent, item["file"])
                dest = work / item["file"]
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(path, dest)
                files.append(item["file"])
            stage = "review"
        else:
            required = dict(row["parameters"])
            if row["operation"] == "salt_quantity":
                required.update({k: row["axes"][k] for k in ("salt", "solvent")})
                if data.get("method_review") != {"$artifact": "method-review.json", "path": []}:
                    raise ValueError("salt analyses require the original bound method qualification report")
                if "computed.json" not in artifact_names(data.get("by_method")):
                    raise ValueError("salt quantities require current-run computed evidence")
            if "fit.json" in row["required_bindings"] and data.get("calibration") != {"$artifact": "fit.json", "path": []}:
                raise ValueError("held-out scoring requires the complete frozen threshold fit")
            if any(data.get(k) != v for k, v in required.items()):
                raise ValueError("research cell differs from protocol-frozen analysis parameters")
            if data.get("operation") != row["operation"]:
                raise ValueError("analysis operation differs from protocol")
            if row["gate"] and data.get("gate_id", data.get("gate")) != row["gate"]:
                raise ValueError("analysis uses the wrong scientific gate")
            if row["gate"]:
                from .science import gate_rules
                if row["gate_spec"] is not None:
                    if "gate_spec" in data and data["gate_spec"] != row["gate_spec"]:
                        raise ValueError("intake attempts to change the protocol threshold")
                    data["gate_spec"] = row["gate_spec"]
                data["gate_spec"] = gate_rules(data)
                write(work / name, data)
            stage = row["requires_stages"][0]
        task = {"id": "evaluate", "engine": "analysis", "stage": stage, "files": files,
                "input": name, "outputs": ["report.json"],
                "resources": record.get("resources", {"cores": 1, "memory_mb": 512, "scratch_mb": 64})}
        write(work / "plan.json", {"name": record["cell"], "tasks": [task]})
        source = work / "plan.json"
    elif kind != "plan":
        raise ValueError(f"unsupported cell kind {kind}")
    plan = load_plan(source)
    stages = {task["stage"] for task in plan["tasks"]}
    if any(task["engine"] in {"probe", "command"} for task in plan["tasks"]):
        raise ValueError("engineering probes cannot implement a research cell")
    if not set(row["requires_stages"]) <= stages:
        raise ValueError(f"{record['cell']}: missing required stages")
    return source.parent, plan


def artifact_names(obj):
    if isinstance(obj, list):
        return set().union(*(artifact_names(v) for v in obj))
    if isinstance(obj, dict):
        if "$artifact" in obj:
            return {obj["$artifact"]}
        return set().union(*(artifact_names(v) for v in obj.values()))
    return set()


def compile_research(protocol, inputs, destination, through=None):
    """Compose reviewed cells. Source locators and cross-cell data connections are frozen."""
    protocol, inputs = Path(protocol).resolve(), Path(inputs).resolve()
    pack, rows = inventory(protocol)
    intake = read(inputs)
    if intake["protocol_sha256"] != _digest(protocol):
        raise ValueError("input manifest does not match the private protocol")
    records = intake["cells"]
    if through:
        dependencies = {s["id"]: s["deps"] for s in pack["steps"]}
        if through not in dependencies:
            raise ValueError(f"unknown study step: {through}")
        selected, remaining = set(), [through]
        while remaining:
            step = remaining.pop()
            if step not in selected:
                selected.add(step)
                remaining.extend(dependencies[step])
        rows = {cell: row for cell, row in rows.items() if row["step"] in selected}
        records = {key: value for key, value in records.items() if key in rows}
    if set(records) != set(rows):
        raise ValueError(f"research intake coverage differs; missing={sorted(rows.keys()-records.keys())} "
                         f"unexpected={sorted(records.keys()-rows.keys())}")
    dst = Path(destination).absolute()
    if dst.exists() or dst.is_symlink():
        raise ValueError("destination exists")
    dst.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=dst.parent) as tmp:
        root = Path(tmp).resolve()
        out = root / "published"
        out.mkdir()
        plans, exports, inputs_by_cell = {}, {}, {}
        tasks, attachments = [], set()
        step_cells = defaultdict(list)
        def copy(path, name):
            target = out / _relative(name, "published asset")
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists() and _digest(target) != _digest(path):
                raise ValueError(f"conflicting static file {name}; give distinct source files unique names")
            if not target.exists():
                shutil.copyfile(path, target)
        for index, (cell, row) in enumerate(rows.items()):
            record = dict(records[cell], cell=cell)
            work = root / str(index)
            work.mkdir()
            local, plan = _local_plan(row, record, inputs.parent, work)
            ids = {t["id"]: f"c{index}-{t['id']}" for t in plan["tasks"]}
            step_cells[row["step"]].append(cell)
            exports[cell] = {}
            members = {t["id"]: t for t in plan["tasks"]}
            kinds = {t["engine"] for t in members.values()}
            if row["kind"] in {"study", "pathway"} and "orca" not in kinds:
                raise ValueError("quantum chemistry cell requires an ORCA task")
            for label, ref in record["exports"].items():
                _identifier(label, "research product label")
                owner, sep, filename = ref.partition("/")
                if not sep or owner not in ids or filename not in members[owner]["outputs"]:
                    raise ValueError(f"unknown export {cell}/{ref}")
                exports[cell][label] = ids[owner] + "/" + filename
            if not set(row["outputs"]) <= set(exports[cell]):
                raise ValueError(f"{cell}: all declared research outputs must be bound")
            plans[cell] = {}
            for task in plan["tasks"]:
                aliases = {}
                for name, source in static_sources(task).items():
                    file = local / source
                    archive = "research-inputs/" + _digest(file)
                    copy(file, archive)
                    aliases[name] = archive
                task["static_sources"] = aliases
                original = task["id"]
                task["id"] = ids[original]
                task["deps"] = [ids[d] for d in task["deps"]]
                task["from_deps"] = {name: ids[ref.partition("/")[0]] + "/" + ref.partition("/")[2]
                                     for name, ref in task["from_deps"].items()}
                plans[cell][original] = task
                tasks.append(task)
            inputs_by_cell[cell] = record.get("bindings", [])
            for name in plan.get("attachments", []):
                archive = f"research-sources/c{index}/{name}"
                copy(local / name, archive)
                attachments.add(archive)
        # Depend on all sinks of named prerequisite groups. Explicit data edges bind a
        # particular labelled product from a specific cell to a consumer's local path.
        step_deps = {s["id"]: s["deps"] for s in pack["steps"]}
        for cell, row in rows.items():
            members = list(plans[cell].values())
            starts = [t for t in members if not t["deps"]]
            prerequisites = {c for s in step_deps[row["step"]] for c in step_cells[s]}
            grouped = {c for c in step_cells[row["step"]] if rows[c]["group"] in row["needs"]}
            for need in row["needs"]:
                if not any(rows[c]["group"] == need for c in grouped):
                    raise ValueError(f"unknown prerequisite group {need}")
            for prev in prerequisites | grouped:
                other = list(plans[prev].values())
                other_ids = {item["id"] for item in other}
                children_deps = {d for t in other for d in t["deps"] if d in other_ids}
                sinks = [t["id"] for t in other if t["id"] not in children_deps]
                for task in starts:
                    task["deps"] = sorted(set(task["deps"]) | set(sinks))
            binding_rows = inputs_by_cell[cell]
            for name, requirement in row["required_bindings"].items():
                matches = [b for b in binding_rows if b["as"] == name]
                if len(matches) != 1:
                    raise ValueError(f"{cell}: required scientific binding missing: {name}")
                b = matches[0]
                producer = rows[b["cell"]]
                anchors = {k: v for k, v in row["axes"].items() if k in producer["axes"]}
                if (producer["group"] != requirement["group"] or
                        any(producer["axes"][k] != v for k, v in anchors.items()) or
                        any(producer["axes"].get(k) != v for k,v in requirement.get("axes", {}).items())):
                    raise ValueError(f"{cell}: wrong producer for {name}")
            for bind in binding_rows:
                if bind["cell"] not in prerequisites | grouped:
                    raise ValueError(f"{cell}: binding outside the declared prerequisite groups")
                t = plans[cell][bind["task"]]
                source = exports[bind["cell"]][bind["output"]]
                if bind["as"] in t["from_deps"]:
                    raise ValueError("duplicate cross-cell binding")
                t["deps"] = sorted(set(t["deps"]) | {source.partition("/")[0]})
                t["from_deps"][bind["as"]] = source
            if row["kind"] == "analysis":
                refs = artifact_names(read(asset(inputs.parent, records[cell]["source"])))
                bound = {b["as"] for b in inputs_by_cell[cell]}
                if refs != bound or (prerequisites | grouped and not bound):
                    raise ValueError(f"{cell}: analysis bindings must match live artifact references")

        for name in pack.get("source_sha256", {}):
            target = "research-sources/protocol/" + name
            copy(asset(protocol.parent, name), target)
            attachments.add(target)
        for source, name in ((protocol, "research-protocol.json"), (inputs, "research-intake.json")):
            copy(source, name)
            attachments.add(name)
        write(out / "research-map.json", {"cells": rows, "exports": exports, "task_ids": {c: {k: t["id"] for k,t in m.items()} for c,m in plans.items()}})
        attachments.add("research-map.json")
        plan = {"name": pack["name"], "study_scope": "paper_inputs", "tasks": tasks}
        if "evidence_plan" in intake:
            ref = intake["evidence_plan"]
            evidence = asset(inputs.parent, ref["path"])
            if _digest(evidence) != ref["sha256"]:
                raise ValueError("prospective evidence plan digest mismatch")
            scope = read(evidence)
            if scope.get("purpose") != "paper_inputs" or scope.get("paper") != pack["name"]:
                raise ValueError("evidence plan must identify this paper and paper_inputs purpose")
            members_by_id = {t["id"]: t for t in tasks}
            for spec in scope["obligations"].values():
                for job in spec["jobs"]:
                    key, _, filename = job.partition("/")
                    if key not in members_by_id or filename not in members_by_id[key]["outputs"]:
                        raise ValueError(f"evidence obligation references undeclared job output: {job}")
            name = "research-sources/evidence-plan.json"
            copy(evidence, name)
            attachments.add(name)
            plan.update(evidence_plan={"path": name, "sha256": _digest(evidence)},
                        research_paper=pack["name"],
                        research_protocol_sha256=_digest(protocol))
        plan["attachments"] = sorted(attachments)
        write(out / "plan.json", plan)
        load_plan(out / "plan.json")
        os.rename(out, dst)
    return {"cells": len(rows), "tasks": len(tasks), "plan": str(dst / "plan.json")}


def prepare_quantum_matrix(protocol, registrypath, destination):
    """Materialize all built-in method sets against approved per-species conformers."""
    protocol, registry = Path(protocol).resolve(), Path(registrypath).resolve()
    _, rows = inventory(protocol)
    species = read(registry)
    selected = {cell: row for cell, row in rows.items() if row["method_sets"]}
    if not selected:
        raise ValueError("this protocol does not specify a built-in molecular method matrix")
    required = {row["axes"]["species"] for row in selected.values()}
    if required - species.keys():
        raise ValueError(f"missing approved species: {sorted(required - species.keys())}")
    dest = Path(destination).absolute()
    if dest.exists() or dest.is_symlink():
        raise ValueError("destination exists")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=dest.parent) as tmp:
        out = Path(tmp) / "matrix"
        for cell, row in selected.items():
            sid, model = row["axes"]["species"], row["axes"]["model"]
            record = species[sid]
            molecule = dict(record, id=sid, claims=row["claims"])
            target = out / cell
            target.mkdir(parents=True)
            for conf in molecule["conformers"]:
                _string(conf["approved_by"], "conformer reviewer")
                file = asset(registry.parent, conf["xyz"])
                if _digest(file) != conf["sha256"]:
                    raise ValueError(f"{sid}: reviewed geometry changed")
                ref = target / _relative(conf["xyz"], "conformer XYZ")
                ref.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(file, ref)
            methodset = asset(protocol.parent, row["method_sets"][model])
            methods = read(methodset)
            for method in methods:
                name = method["template"]
                src = asset(methodset.parent, name)
                ref = target / _relative(name, "method template")
                ref.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(src, ref)
            study = {"name": cell, "scope": "paper_inputs", "species": [molecule],
                     "claim_ids": row["claims"], "methods": methods}
            write(target / "study.json", study)
            write(target / "links.json", {"exports": {k: "matrix/matrix.json" for k in row["outputs"]}})

            check_study_design(row, study, target / "study.json")
        os.rename(out, dest)
    return {"prepared_cells": len(selected), "species": len(required), "directory": str(dest)}


def assemble(protocol, directory, destination):
    """Discover cell inputs in a private directory without filling missing science."""
    protocol, directory = Path(protocol).resolve(), Path(directory).resolve()
    output = Path(destination).absolute()
    if output.parent != directory or output.exists() or output.is_symlink():
        raise ValueError("intake must be a new file in the private input directory")
    _, rows = inventory(protocol)
    records, missing = {}, []
    for cell, row in rows.items():
        filename = {"review": "review.json", "analysis": "request.json", "study": "study.json",
                    "ion_pairs": "poses.json", "plan": "plan.json", "pathway": "pathway.json"}[row["kind"]]
        source = cell + "/" + filename
        if not (directory / source).is_file():
            missing.append(source)
            continue
        file = asset(directory, source)
        links = cell + "/links.json"
        if (directory / links).exists():
            record = read(asset(directory, links))
            if set(record) - {"exports", "bindings", "resources"}:
                raise ValueError(f"unknown keys in {links}")
        else:
            if row["kind"] not in {"review","analysis"}:
                missing.append(links)
                continue
            record = {"exports": {label: "evaluate/report.json" for label in row["outputs"]}}
        records[cell] = dict(record, source=source, source_sha256=_digest(file))
    write(output, {"protocol_sha256": _digest(protocol), "cells": records})
    return {"intake": str(output), "registered": len(records), "expected": len(rows), "missing": missing}


def main(argv=None):
    parser = argparse.ArgumentParser(prog="chemflow2 research")
    sub = parser.add_subparsers(dest="verb", required=True)
    inv = sub.add_parser("inventory")
    inv.add_argument("protocol", type=Path)
    matrix = sub.add_parser("prepare-quantum")
    matrix.add_argument("protocol", type=Path)
    matrix.add_argument("registry", type=Path)
    matrix.add_argument("destination", type=Path)
    asm = sub.add_parser("assemble")
    asm.add_argument("protocol", type=Path)
    asm.add_argument("directory", type=Path)
    asm.add_argument("destination", type=Path)
    comp = sub.add_parser("compile")
    comp.add_argument("protocol", type=Path)
    comp.add_argument("inputs", type=Path)
    comp.add_argument("destination", type=Path)
    comp.add_argument("--through", help="compile the named step and its complete prerequisite closure")
    args = parser.parse_args(argv)
    if args.verb == "inventory":
        pack, rows = inventory(args.protocol)
        result = {"name": pack["name"], "protocol_sha256": _digest(args.protocol), "cells": rows}
    elif args.verb == "prepare-quantum":
        result = prepare_quantum_matrix(args.protocol, args.registry, args.destination)
    elif args.verb == "assemble":
        result = assemble(args.protocol, args.directory, args.destination)
    else:
        result = compile_research(args.protocol, args.inputs, args.destination, args.through)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
