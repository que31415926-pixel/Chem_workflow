"""Aggregate only composition-balanced frequency jobs on a common solution reference."""

import math
from collections import Counter
from pathlib import Path

from .network import read_json
from .qm import orca_audit, read_frames


def reaction_energy(recipe):
    recipe = Path(recipe)
    spec = read_json(recipe)
    if not spec["reactants"] or not spec["products"]:
        raise ValueError("Empty reaction side")
    seen = {}
    formula = Counter()
    charge = 0
    delta_n = 0
    delta_eh = 0.0
    enthalpies = []
    all_enthalpies = True
    for side, sign in [("reactants", -1), ("products", 1)]:
        for relative, number in spec[side].items():
            if not isinstance(number, int) or isinstance(number, bool) or number < 1:
                raise ValueError(
                    "Stoichiometric coefficients must be positive integers"
                )
            folder = recipe.parent / relative
            result = orca_audit(folder)
            if (
                result["stage"] not in ("min", "ts")
                or result["status"] != "numerical_checks_passed"
            ):
                raise ValueError(f"Unusable thermochemistry in {relative}")
            if side == "reactants" and result["stage"] != "min":
                raise ValueError("Expected minimum on reactant side")
            elements, _ = next(read_frames(folder / "geom.xyz"))
            for symbol in elements:
                formula[symbol] += sign * number
            for name in ("project", "temperature_K", "method", "solvent"):
                if name in seen and seen[name] != result[name]:
                    raise ValueError("Unlike thermochemistry: " + name)
                seen[name] = result[name]
            if abs(result["temperature_K"] - result["reported_temperature_K"]) > 0.05:
                raise ValueError("Temperature mismatch")
            if abs(result["reported_pressure_atm"] - 1.0) > 0.005:
                raise ValueError("Expected 1 atmosphere RRHO reference")
            charge += sign * number * result["charge"]
            delta_n += sign * number
            delta_eh += sign * number * result["thermal_G_Eh"]
            if result["enthalpy_Eh"] is None:
                all_enthalpies = False
            else:
                enthalpies.append(sign * number * result["enthalpy_Eh"])
    if any(formula.values()) or charge:
        raise ValueError("Stoichiometry/charge not balanced")
    t = seen["temperature_K"]
    correction = delta_n * 0.00198720426 * t * math.log(0.082057366 * t)
    return {
        **seen,
        "delta_n": delta_n,
        "delta_H_1atm_kcal_mol": sum(enthalpies) * 627.509474
        if all_enthalpies
        else None,
        "delta_G_1atm_kcal_mol": delta_eh * 627.509474,
        "correction_1atm_to_1M_kcal_mol": correction,
        "delta_G_1M_kcal_mol": delta_eh * 627.509474 + correction,
        "interpretation": "All solutes; ideal solution. Handle gas H2 and redox references separately. Review conformer, spin and IRC before interpretation.",
    }
