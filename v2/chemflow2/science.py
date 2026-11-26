"""Numerical analyses on explicit, frozen chemistry inputs. No external engines."""

import math
from fractions import Fraction

GAS_CONSTANT = 8.31446261815324 / 4184.0  # kcal mol^-1 K^-1
FARADAY = 96485.33212331002  # C/mol

# Legacy private packs can provide their original, fixed benchmark rules.
try:
    from .v4_rules import GATES
except ImportError:
    GATES = {}


def gate_rules(data):
    gate = _id(_get(data, "gate_id"))
    registered = GATES.get(gate)
    if "gate_spec" not in data:
        if registered is None:
            raise ValueError(f"unknown gate: {gate}. Supply a frozen gate_spec.")
        return registered
    rules = {}
    for name, row in _mapping(data["gate_spec"]).items():
        if not isinstance(row, (tuple, list)) or len(row) != 3:
            raise ValueError(f"invalid benchmark rule for {name}")
        rule, unit, limit = row
        if ((rule not in ("mae", "max_abs") and not (rule is None and name == "redox_error"))
                or not _id(unit)):
            raise ValueError("invalid benchmark aggregation or unit")
        if limit is not None:
            if _number(limit) < 0:
                raise ValueError("negative benchmark threshold")
        elif name != "electronic_difference_eh":
            raise ValueError("missing benchmark threshold")
        rules[name] = (rule, unit, limit)
    if not rules or (registered is not None and rules != registered):
        raise ValueError("empty gate or mismatch with registered benchmark")
    return rules


def _get(data, key):
    if not isinstance(data, dict) or key not in data:
        raise ValueError(f"required field: {key}")
    return data[key]


def _id(value):
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise ValueError("expected a non-empty identifier without surrounding whitespace")
    return value


def _number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"expected a finite number, got {value!r}")
    return value


def _int(value):
    if type(value) is not int:
        raise ValueError(f"expected an integer, got {value!r}")
    return value


def _positive(value):
    if _number(value) <= 0:
        raise ValueError("expected a positive number")
    return value


def _ids(values):
    if not isinstance(values, list) or not values:
        raise ValueError("expected a non-empty list of distinct identifiers")
    result = [_id(v) for v in values]
    if len(result) != len(set(result)):
        raise ValueError("duplicate identifier")
    return result


def _mapping(value):
    if not isinstance(value, dict):
        raise ValueError("expected an object keyed by identifier")
    for key in value:
        _id(key)
    return value


def _unit(data, target):
    if _get(data, "unit") != target:
        raise ValueError(f"expected unit {target!r}")


def _check_subset(records, expected):
    extra = sorted(set(records) - set(expected))
    if extra:
        raise ValueError(f"undeclared records: {extra}")


def _identity(data, fields):
    return {key: _id(_get(data, key)) for key in fields}


def _ensemble(data):
    metadata = _identity(data, ("species_id", "standard_state_id", "condition_id"))
    _unit(data, "kcal/mol")
    temp = _positive(_get(data, "temperature_K"))
    conformers = _mapping(_get(data, "degeneracies"))
    if not conformers:
        raise ValueError("declare all required conformers and their degeneracies")
    for value in conformers.values():
        _positive(_int(value))
    methods = _ids(_get(data, "methods"))
    energies = _mapping(_get(data, "energies"))
    _check_subset(energies, methods)
    missing = {}
    for method in methods:
        rows = _mapping(energies.get(method, {}))
        _check_subset(rows, conformers)
        for row in rows.values():
            _number(_get(row, "g"))
            _id(_get(row, "source_id"))
        uncovered = sorted(set(conformers) - set(rows))
        if uncovered:
            missing[method] = uncovered
    if missing:
        return {"status": "missing", "missing_by_method": missing}

    result = {}
    rt = GAS_CONSTANT * temp
    for method, rows in energies.items():
        baseline = min(row["g"] for row in rows.values())
        terms = {conf: conformers[conf] * math.exp(-(row["g"] - baseline) / rt)
                 for conf, row in rows.items()}
        total = math.fsum(terms.values())
        result[method] = {"g": baseline - rt * math.log(total),
                          "populations": {conf: w / total for conf, w in terms.items()},
                          "source_ids": {conf: row["source_id"] for conf, row in rows.items()}}
    return dict(metadata, status="pass", unit="kcal/mol", temperature_K=temp, by_method=result)


def _reaction(data):
    _unit(data, "kcal/mol")
    metadata = _identity(data, ("reaction_id", "method_id", "condition_id"))
    temp = _positive(_get(data, "temperature_K"))
    coeffs = _mapping(_get(data, "coefficients"))
    for n in coeffs.values():
        _int(n)
        if n == 0:
            raise ValueError("stoichiometric coefficients must be nonzero")
    if not any(n < 0 for n in coeffs.values()) or not any(n > 0 for n in coeffs.values()):
        raise ValueError("a reaction requires reactants (-) and products (+)")
    states = _mapping(_get(data, "standard_states"))
    if set(states) != set(coeffs):
        raise ValueError("freeze a standard state for each reaction participant")
    for state in states.values():
        _id(state)
    rows = _mapping(_get(data, "species"))
    _check_subset(rows, coeffs)
    atoms, charge, terms = {}, 0, []
    for species, row in rows.items():
        _unit(row, "kcal/mol")
        _id(_get(row, "source_id"))
        if _get(row, "species_id") != species or _get(row, "standard_state_id") != states[species]:
            raise ValueError(f"{species}: incompatible identity or standard state")
        for field in ("method_id", "condition_id"):
            if _get(row, field) != metadata[field]:
                raise ValueError(f"{species}: incompatible {field}")
        if _positive(_get(row, "temperature_K")) != temp:
            raise ValueError(f"{species}: incompatible temperature")
        coefficient = coeffs[species]
        charge += coefficient * _int(_get(row, "charge"))
        formula = _mapping(_get(row, "composition"))
        if not formula:
            raise ValueError(f"{species}: missing composition")
        for atom, count in formula.items():
            atoms[atom] = atoms.get(atom, 0) + coefficient * _positive(_int(count))
        terms.append(coefficient * _number(_get(row, "g")))
    missing = sorted(set(coeffs) - set(rows))
    if missing:
        return {"status": "missing", "missing_species": missing}
    if charge or any(atoms.values()):
        raise ValueError("reaction does not conserve composition or charge")
    return dict(metadata, status="pass", delta_g=math.fsum(terms), unit="kcal/mol",
                temperature_K=temp, standard_states={k: v["standard_state_id"] for k, v in rows.items()},
                source_ids={k: v["source_id"] for k, v in rows.items()})


def _redox(data):
    _unit(data, "kcal/mol")
    metadata = _identity(data, ("analyte_couple_id", "reference_couple_id", "condition_id", "source_id"))
    if metadata["analyte_couple_id"] == metadata["reference_couple_id"]:
        raise ValueError("specify distinct analyte and reference couples")
    if _get(data, "cycle") != "reduce_analyte_oxidize_reference":
        raise ValueError("expected cycle reduce_analyte_oxidize_reference")
    electrons = _positive(_int(_get(data, "electrons_transferred")))
    dg = _number(_get(data, "cycle_delta_g"))
    return dict(metadata, status="pass", potential_vs_reference=-dg * 4184.0 / (electrons * FARADAY),
                unit="V", electrons_transferred=electrons)


def _fraction(value):
    return Fraction(str(_number(value)))


def _method_settings(spec, method, expected_sets):
    required = ["method_id"]
    if "solvent_difference" in spec:
        required.extend(["solvent_models", "reaction_contrasts"])
    if "electronic_difference_eh" in spec:
        required.append("electronic_tolerance_eh")
    if "redox_error" in spec:
        required.append("redox_aggregation")
    missing = [k for k in required if k not in method]
    if missing:
        return missing
    _id(method["method_id"])
    if "electronic_difference_eh" in spec and _number(method["electronic_tolerance_eh"]) < 0:
        raise ValueError("electronic noise tolerance must be nonnegative and expressed in Eh")
    if "redox_error" in spec and method["redox_aggregation"] not in ("mae", "max_abs"):
        raise ValueError("freeze redox_aggregation as mae or max_abs")
    if "solvent_difference" in spec:
        if len(_ids(method["solvent_models"])) != 2:
            raise ValueError("freeze two distinct solvent models in order")
        contrasts = _mapping(method["reaction_contrasts"])
        expected = set(expected_sets.get("solvent_difference", []))
        if set(contrasts) != expected or any(len(_ids(x)) != 2 for x in contrasts.values()):
            raise ValueError("freeze two ordered, distinct reaction IDs for every solvent contrast")
    return []


def _metric_error(metric, name, pair, method):
    _identity(pair, ("left_source_id", "right_source_id"))
    if metric == "solvent_difference":
        reactions = method["reaction_contrasts"][name]
        ddg = []
        for side, model in zip(("left", "right"), method["solvent_models"]):
            values = _mapping(_get(pair, side))
            if _get(values, "model_id") != model:
                raise ValueError("the observed solvent model differs from the frozen order")
            dg = _mapping(_get(values, "reaction_delta_g"))
            if set(dg) != set(reactions):
                raise ValueError("solvent contrast must contain exactly the two declared reaction free energies")
            ddg.append(_fraction(dg[reactions[0]]) - _fraction(dg[reactions[1]]))
        return abs(ddg[0] - ddg[1])
    return abs(_fraction(_get(pair, "left")) - _fraction(_get(pair, "right")))


def _gate(data):
    gate = _id(_get(data, "gate_id"))
    spec = gate_rules(data)
    inputs = _mapping(_get(data, "metrics"))
    expected_sets = _mapping(_get(data, "required_pairs"))
    method = _mapping(data.get("gate_method", {}))
    _check_subset(inputs, spec)
    _check_subset(expected_sets, spec)
    for expected in expected_sets.values():
        _ids(expected)
    missing_method = _method_settings(spec, method, expected_sets)
    missing_sets = sorted(set(spec) - set(expected_sets))
    if missing_method or missing_sets:
        return {"gate_id": gate, "status": "missing", "missing_gate_method_fields": missing_method,
                "missing_required_sets": missing_sets, "metrics": {}}
    details = {}
    for metric, (rule, unit, limit) in spec.items():
        if metric == "redox_error":
            rule = method["redox_aggregation"]
        if metric == "electronic_difference_eh":
            limit = method["electronic_tolerance_eh"]
        expected = expected_sets[metric]
        pairs = {}
        if metric in inputs:
            _unit(inputs[metric], unit)
            pairs = _mapping(_get(inputs[metric], "pairs"))
        _check_subset(pairs, expected)
        errors = [_metric_error(metric, name, pair, method) for name, pair in pairs.items()]
        missing = sorted(set(expected) - set(pairs))
        detail = {"rule": rule, "unit": unit, "limit": limit, "missing": missing}
        if missing:
            detail["status"] = "missing"
        else:
            error = sum(errors) / len(errors) if rule == "mae" else max(errors)
            passed = error <= Fraction(str(limit))
            detail.update(status="pass" if passed else "fail", value=float(error))
        details[metric] = detail
    statuses = [m["status"] for m in details.values()]
    status = "fail" if "fail" in statuses else ("missing" if "missing" in statuses else "pass")
    return {"gate_id": gate, "status": status, "gate_method": method, "metrics": details}


def _direction(data):
    _unit(data, "kcal/mol")
    left, right = (_id(_get(data, k)) for k in ("left_id", "right_id"))
    if left == right:
        raise ValueError("the compared entities must differ")
    target = _int(_get(data, "expected_sign"))
    if target not in (-1, 1):
        raise ValueError("expected_sign must be -1 or 1")
    required = _ids(_get(data, "required_methods"))
    pairs = _mapping(_get(data, "pairs"))
    _check_subset(pairs, required)
    diffs = {}
    for method, pair in pairs.items():
        _identity(pair, ("left_source_id", "right_source_id"))
        diffs[method] = _number(_get(pair, "left")) - _number(_get(pair, "right"))
    missing = sorted(set(required) - set(pairs))
    status = "fail" if any(target * v <= 0 for v in diffs.values()) else ("missing" if missing else "pass")
    return {"status": status, "left_id": left, "right_id": right, "expected_sign": target,
            "unit": "kcal/mol", "differences_left_minus_right": diffs, "missing_methods": missing}


def evaluate(data: dict) -> dict:
    """Evaluate one declared scientific operation; raise ValueError for invalid input."""
    operations = {"ensemble": _ensemble, "reaction": _reaction, "redox": _redox,
                  "gate": _gate, "direction": _direction}
    operation = _id(_get(data, "operation"))
    if operation == "collect_matrix":
        from .research import collect_matrix
        return collect_matrix(data)
    if operation in ("two_point_cbs", "composite_thermochemistry"):
        from .composite import evaluate as composite
        return composite(data)
    if operation not in operations:
        raise ValueError(f"unknown science operation: {operation}")
    protocol = _id(_get(data, "protocol_id"))
    return dict(operations[operation](data), operation=operation, protocol_id=protocol)
