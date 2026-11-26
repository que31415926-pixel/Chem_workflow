"""Stoichiometric solution thermochemistry using matched ensemble Gibbs energies."""
from collections import Counter
import math
import re


def _quantity(value, label):
    if type(value) not in (float, int) or not math.isfinite(value):
        raise ValueError(f"{label} must be finite")
    return value


def _identity(row):
    identifiers = ("method", "standard_state", "solvent", "solvent_model")
    if any(not isinstance(row.get(k), str) or not row[k].strip() for k in identifiers):
        raise ValueError("incomplete solution reference conditions")
    temp = _quantity(row["temperature_K"], "temperature")
    if temp <= 0:
        raise ValueError("positive absolute temperature required")
    return (*[row[k] for k in identifiers], temp)


def _difference(terms, records, *, solvation=False):
    if not isinstance(terms, dict) or len(terms) < 2:
        raise ValueError("a cycle needs distinct reactant and product entries")
    atoms, charge, condition = Counter(), 0, None
    values, errors, phases = [], [], set()
    for key, coeff in terms.items():
        if type(coeff) is not int or coeff == 0:
            raise ValueError("use signed nonzero integer stoichiometric coefficients")
        rec = records[key]
        composition = rec["composition"]
        if (not isinstance(composition, dict) or not composition or
                any(not re.fullmatch(r"[A-Z][a-z]?", k) or type(v) is not int or v < 1 for k, v in composition.items())
                or type(rec["charge"]) is not int):
            raise ValueError("invalid species composition or charge")
        for element, n in composition.items():
            atoms[element] += coeff * n
        charge += coeff * rec["charge"]
        conditions = _identity(rec)
        if solvation:
            phases.add(conditions[2])
            conditions = (conditions[0], conditions[1], conditions[-1])
        if condition is None:
            condition = conditions
        if condition != conditions:
            raise ValueError("cycle legs use different methods, temperatures, solvents or standard states")
        if not re.fullmatch(r"[0-9a-f]{64}", rec["energy_sha256"]):
            raise ValueError("missing Gibbs energy evidence digest")
        factor = {"hartree": 627.5094740631, "kcal/mol": 1, "kJ/mol": 1/4.184}[rec["unit"]]
        energy = _quantity(rec["gibbs"], "Gibbs energy")
        bound = _quantity(rec["absolute_error_bound"], "absolute Gibbs error bound")
        if bound < 0:
            raise ValueError("absolute error bound cannot be negative")
        values.append(coeff * energy * factor)
        errors.append(abs(coeff) * bound * factor)
    if any(atoms.values()):
        raise ValueError("cycle fails elemental balance")
    if solvation and ("gas" not in phases or len(phases) != 2):
        raise ValueError("solvation requires gas and solution legs")
    return math.fsum(values), math.fsum(errors), charge, condition


def evaluate(data):
    kind = data["kind"]
    species = data["species"]
    terms = data["stoichiometry"]
    if kind not in ("reaction", "solvation", "relative_reduction"):
        raise ValueError("unknown solution cycle kind")
    if kind == "solvation":
        if sorted(terms.values()) != [-1, 1]:
            raise ValueError("solvation requires a one-to-one transfer")
        reactant, product = [species[next(k for k, v in terms.items() if v == c)] for c in (-1, 1)]
        if (reactant["solvent"] != "gas" or reactant["species_id"] != product["species_id"]
                or reactant["multiplicity"] != product["multiplicity"]):
            raise ValueError("solvation must transfer the same species/state from gas to solution")
    dg, err, charge, condition = _difference(terms, species, solvation=kind == "solvation")
    result = {"operation": "solution_cycle", "kind": kind, "delta_g_kcal_mol": dg,
              "error_bound_kcal_mol": err}
    used = set(terms)
    if kind == "relative_reduction":
        electrons = -charge
        ref = data["reference"]
        ref_terms = ref["stoichiometry"]
        if set(ref_terms) & used:
            raise ValueError("relative potential needs an independently identified reference couple")
        used.update(ref_terms)
        ref_dg, ref_err, ref_charge, ref_condition = _difference(ref_terms, species)
        reference_id = ref["id"]
        if (not isinstance(reference_id, str) or not reference_id.strip() or electrons <= 0
                or ref_charge >= 0 or ref_condition != condition):
            raise ValueError("target and reference must be reductions under identical conditions")
        # A same-solution isodesmic cycle cancels the chemical potential of the electron.
        faraday_kcal = 96485.33212 / 4184
        potential = (-dg / electrons + ref_dg / (-ref_charge)) / faraday_kcal
        uncertainty = (err / electrons + ref_err / (-ref_charge)) / faraday_kcal
        result.update(reference=reference_id, potential_V=potential,
                      potential_interval_V=[potential - uncertainty, potential + uncertainty])
    elif charge:
        raise ValueError("cycle is not charge balanced")
    result["energy_evidence"] = {key: species[key]["energy_sha256"] for key in sorted(used)}
    return result


def salt_quantity(data):
    """Compute a prospectively qualified method-by-method solution observable."""
    from .mechanism import _digest, _qualify_path

    observable = data["quantity"]
    contract = data["method_review"]
    if contract.get("decision") != "pass" or not contract.get("reviewer") or not contract.get("items"):
        raise ValueError("a recorded method qualification review is required")
    settings = contract["qualified_methods"]
    calculations = data["by_method"]
    if not settings or set(settings) != set(calculations):
        raise ValueError("calculate every qualified method exactly once")
    results = {}
    for method, experiment in calculations.items():
        qualification = settings[method]
        _digest(qualification["method_settings_sha256"], "method settings")
        if qualification["method_settings_sha256"] != experiment["method_settings_sha256"]:
            raise ValueError("calculation differs from the qualified method")
        if experiment["salt"] != data["salt"] or experiment["solvent"] != data["solvent"]:
            raise ValueError("scientific identity differs from the declared salt or solvent")
        if observable in {"ionic_reaction", "one_electron_reduction"}:
            cycle = experiment["cycle"]
            target = cycle["stoichiometry"]
            species = cycle["species"]
            if any(r["method"] != method or r["solvent"] != data["solvent"]
                   or r["solvent_model"] != qualification["solvent_model"] for r in species.values()):
                raise ValueError("reference conditions differ from the qualified method")
            kind = "relative_reduction" if observable == "one_electron_reduction" else "reaction"
            if cycle["kind"] != kind or not any(v < 0 for v in target.values()):
                raise ValueError("thermochemical cycle does not match the requested quantity")
            if experiment["salt_reactant"] not in target or target[experiment["salt_reactant"]] != -1:
                raise ValueError("one salt formula unit must be consumed")
            if observable == "one_electron_reduction":
                if sum(v*species[k]["charge"] for k, v in target.items()) != -1:
                    raise ValueError("salt reduction must transfer exactly one electron")
                if cycle["reference"]["id"] != qualification["reference_couple"]:
                    raise ValueError("potential reference differs from the reviewed couple")
            elif sum(v for v in target.values() if v > 0) < 2:
                raise ValueError("ionic dissociation requires distinct product fragments")
            results[method] = evaluate(cycle)
        elif observable in {"ionic_barrier", "reduced_CN_barrier"}:
            point = experiment["path"]
            spec = experiment["specification"]
            if (spec["method_id"] != method or spec["reactant"] != experiment["salt_reactant"]
                    or point["reactant"]["charge"] != experiment["oxidized_assembly_charge"] -
                    int(observable == "reduced_CN_barrier")):
                raise ValueError("activation path has wrong initial state, oxidation level or method")
            result = _qualify_path(point, spec, experiment["condition_id"],
                                   _quantity(qualification["imaginary_tolerance_cm1"], "imaginary cutoff"))
            results[method] = result
        else:
            raise ValueError("unknown salt observable")
    qualified = all(row.get("status", "qualified") == "qualified" for row in results.values())
    return {"operation": "salt_quantity", "quantity": observable, "salt": data["salt"],
            "solvent": data["solvent"], "by_method": results,
            "status": "pass" if qualified else "indeterminate"}
