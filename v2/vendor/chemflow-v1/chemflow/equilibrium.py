"""Coupled gas-liquid equilibrium for ArH + HBpin <=> ArBpin + H2.

This is the C1 module of the v4 plan (K2 first line): it solves the
*thermodynamic equilibrium* of a closed vessel with a liquid phase and an
ideal-gas headspace coupled through Henry's law. It is not a kinetic model,
does not contain activity coefficients, and does not describe catalyst
speciation; it returns an upper bound for a finite-time yield.

Conventions
-----------
* Concentrations are mol/L of liquid volume.
* Gases are ideal. R = 0.083144626 bar L mol-1 K-1.
* Henry volatility H = P(H2) / [H2] in bar L mol-1; larger H means a more
  volatile (less soluble) gas.
* Reaction standard states: solutes 1 M, H2 either 1 M (aq) or 1 bar (g).
  ``deltaG_standard`` selects the convention; the default is "H2(aq,1M)".
  For "H2(g,1bar)" the aq standard free energy is obtained as
  dG(aq) = dG(g) + RT ln H, because mu(aq) = mu(g) + RT ln H at phase
  equilibrium with P = H [H2].
* One H2 is produced per aryl boronate formed. Reactions share one borane
  inventory and one H2 pool. Reactions that form regioisomers of the same
  arene must declare the same ``arene_group``; that group then shares one
  arene inventory (A0 - sum of its products).
"""

import copy
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

R_BAR_L = 0.0831446261815324
R_J = 8.31446261815324
STANDARDS = ("H2(aq,1M)", "H2(g,1bar)")
SWEEP_FIELDS = (
    "temperature_K",
    "liquid_volume_L",
    "headspace_volume_L",
    "borane_equivalents",
    "borane_M",
)
TINY = 1e-300


def _number(value, name, positive=False, nonnegative=False):
    try:
        value = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a number") from None
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    if positive and value <= 0:
        raise ValueError(f"{name} must be positive")
    if nonnegative and value < 0:
        raise ValueError(f"{name} must be non-negative")
    return value


def henry_volatility(henry, temperature_K):
    """H(T) in bar L mol-1 from a reference value and optional enthalpy."""
    if not isinstance(henry, dict) or "value_bar_L_per_mol" not in henry:
        raise ValueError("recipe.henry must contain value_bar_L_per_mol")
    h_ref = _number(henry["value_bar_L_per_mol"], "henry.value_bar_L_per_mol", positive=True)
    t_ref = _number(
        henry.get("reference_temperature_K", temperature_K),
        "henry.reference_temperature_K",
        positive=True,
    )
    enthalpy = _number(henry.get("enthalpy_kJ_per_mol", 0.0), "henry.enthalpy_kJ_per_mol")
    if enthalpy == 0.0:
        return h_ref
    # van't Hoff for dissolution: ln(H/H_ref) = dH_sol/R * (1/T - 1/T_ref)
    return h_ref * math.exp(enthalpy * 1000.0 / R_J * (1.0 / temperature_K - 1.0 / t_ref))


def gas_liquid_factor(henry, temperature_K, liquid_volume_L, headspace_volume_L):
    """Return (gamma, H) with [H2]_liquid = total_H2_M / gamma."""
    h = henry_volatility(henry, temperature_K)
    gamma = 1.0 + h * headspace_volume_L / (R_BAR_L * temperature_K * liquid_volume_L)
    return gamma, h


def _reaction_thermo(reaction, temperature_K, h):
    if "id" not in reaction:
        raise ValueError("each reaction needs an id")
    dg_ref = _number(reaction["deltaG_kJ_per_mol"], f"{reaction['id']}.deltaG_kJ_per_mol")
    standard = reaction.get("deltaG_standard", "H2(aq,1M)")
    if standard not in STANDARDS:
        raise ValueError(f"{reaction['id']}.deltaG_standard must be one of {STANDARDS}")
    t_ref = _number(
        reaction.get("deltaG_reference_temperature_K", temperature_K),
        f"{reaction['id']}.deltaG_reference_temperature_K",
        positive=True,
    )
    enthalpy = reaction.get("deltaH_kJ_per_mol")
    if enthalpy is None:
        dg = dg_ref
    else:
        if "deltaG_reference_temperature_K" not in reaction:
            raise ValueError(
                f"{reaction['id']}: deltaH_kJ_per_mol requires deltaG_reference_temperature_K"
            )
        enthalpy = _number(enthalpy, f"{reaction['id']}.deltaH_kJ_per_mol")
        dg = temperature_K * (dg_ref / t_ref + enthalpy * (1.0 / temperature_K - 1.0 / t_ref))
    if standard == "H2(g,1bar)":
        dg += R_J * temperature_K * math.log(h) / 1000.0
    k = math.exp(-dg * 1000.0 / (R_J * temperature_K))
    return dg, k


def _single_root(a, b, k, gamma):
    """One reaction: solve K = x*(x/gamma)/((a-x)(b-x)) on [0, min(a,b)]."""
    if b <= 0.0 or a <= 0.0:
        return 0.0
    upper = min(a, b)

    def residual(x):
        return k * gamma * (a - x) * (b - x) - x * x

    if residual(upper) >= 0.0:
        return upper
    return brentq(residual, 0.0, upper, xtol=1e-300, rtol=1e-15, maxiter=5000)


def _decimal_root(group_total, ksum, b0, gamma):
    """High-precision bisection for the total extent.

    Double precision cannot resolve sum(S_g) - S once borane is nearly
    exhausted; Decimal arithmetic removes that cancellation.
    """
    from decimal import Decimal, localcontext

    with localcontext() as ctx:
        ctx.prec = 60
        arene = [Decimal.from_float(float(x)) for x in group_total]
        constants = [Decimal.from_float(float(x)) for x in ksum]
        g = Decimal.from_float(float(gamma))
        b = Decimal.from_float(float(b0))

        def residual(total):
            borane = b - total
            out = Decimal(0)
            for a, ksum_g in zip(arene, constants):
                c = g * ksum_g * borane / total
                out += a * c / (Decimal(1) + c)
            return out - total

        lo = Decimal("1e-300")
        hi = b
        if residual(lo) <= 0:
            raise ValueError("equilibrium solve failed to bracket the total extent")
        if residual(hi) >= 0:
            return float(hi), True
        for _ in range(320):
            mid = (lo + hi) / 2
            if residual(mid) > 0:
                lo = mid
            else:
                hi = mid
        return float((lo + hi) / 2), False


def _numeric_roots(group_total, group_index, b0, k, gamma):
    """Return per-reaction extents; reactions in a group share the arene pool.

    Within a group x_i = S_g K_i / Ksum_g. For a trial total extent S the
    group equation S_g S = gamma Ksum_g (A_g - S_g)(B0 - S) gives
    S_g = C A_g / (1 + C) with C = gamma Ksum_g (B0 - S) / S, so the only
    remaining unknown is the scalar S. A bounded one-dimensional root solve is
    exact and stable even when borane or arene is nearly exhausted.
    """
    groups = len(group_total)
    ksum = np.zeros(groups)
    for i, g in enumerate(group_index):
        ksum[g] += k[i]
    if b0 <= 0.0:
        return np.zeros(len(k))

    def group_extents(total):
        borane = b0 - total
        out = np.empty(groups)
        for g in range(groups):
            log_c = (
                math.log(gamma)
                + math.log(ksum[g])
                + math.log(borane)
                - math.log(total)
            )
            if log_c > 700.0:
                out[g] = group_total[g]
            elif log_c < -700.0:
                out[g] = 0.0
            else:
                c = math.exp(log_c)
                if c < 1.0:
                    # stable for small C; avoids A - A/(1+C) cancellation
                    out[g] = group_total[g] * c / (1.0 + c)
                else:
                    out[g] = group_total[g] - group_total[g] / (1.0 + c)
        return out

    def residual(total):
        return float(group_extents(total).sum() - total)

    if groups == 1:
        root = _single_root(float(group_total[0]), b0, float(ksum[0]), gamma)
    else:
        root, at_boundary = _decimal_root(group_total, ksum, b0, gamma)
        if at_boundary or not math.isfinite(root) or root >= b0:
            root = float(np.nextafter(b0, 0.0))
    s = group_extents(root)
    extents = np.zeros(len(k))
    for i, g in enumerate(group_index):
        extents[i] = s[g] * k[i] / ksum[g]
    # The scalar root fixes the total extent; rescaling removes the last few
    # ulps of cancellation in sum(S_g) - S while preserving each group ratio.
    obtained = float(extents.sum())
    if obtained > 0.0:
        extents *= root / obtained
    return extents


def _check_equilibrium(x, group_total, group_index, b0, gamma, k):
    """Reject solutions that violate equilibrium or stoichiometry.

    Returns the group indices accepted at a stoichiometric boundary, where the
    logarithmic residual is not representable in double precision.
    """
    total = float(x.sum())
    if total > b0 + 1e-9 * max(b0, 1.0):
        raise ValueError("equilibrium solve consumed more borane than supplied")
    per_group = np.zeros(len(group_total))
    for i, g in enumerate(group_index):
        per_group[g] += x[i]
    for g, used in enumerate(per_group):
        if used > group_total[g] + 1e-9 * max(group_total[g], 1.0):
            raise ValueError("regioisomer channel consumed more arene than the shared pool")
    if total == 0.0:
        return []
    borane = b0 - total
    h2 = total / gamma
    if borane <= b0 * 1e-6:
        # Borane is exhausted to within double-precision resolution. Every
        # channel is then at the stoichiometric boundary; individual log
        # residuals are not representable or meaningful.
        return sorted({int(g) for g, used in enumerate(per_group) if used > 0.0})
    group_kmax = np.zeros(len(group_total))
    for i, g in enumerate(group_index):
        group_kmax[g] = max(group_kmax[g], k[i])
    limited = set()
    for i, g in enumerate(group_index):
        arene = group_total[g] - per_group[g]
        at_group_bound = per_group[g] >= group_total[g] * (1.0 - 1e-6)
        at_borane_bound = borane <= b0 * 1e-6
        boundary_complete = (at_group_bound or at_borane_bound) and group_kmax[g] >= 1e4
        if arene <= 0.0 or borane <= 0.0:
            if boundary_complete:
                limited.add(int(g))
                continue
            raise ValueError("equilibrium state is outside the stoichiometric domain")
        q = x[i] * h2 / (arene * borane)
        error = abs(math.log(max(q, TINY)) - math.log(k[i]))
        if error > 1e-5:
            if boundary_complete:
                limited.add(int(g))
                continue
            raise ValueError(
                f"equilibrium residual {error:.3g} for reaction {i}; "
                "check bounds, K or starting values"
            )
    return sorted(limited)


def _core(recipe):
    temperature = _number(recipe.get("temperature_K"), "temperature_K", positive=True)
    liquid_volume = _number(recipe.get("liquid_volume_L"), "liquid_volume_L", positive=True)
    headspace = _number(recipe.get("headspace_volume_L", 0.0), "headspace_volume_L", nonnegative=True)
    henry = recipe.get("henry")
    gamma, h = gas_liquid_factor(henry, temperature, liquid_volume, headspace)

    reactions = recipe.get("reactions")
    if not isinstance(reactions, list) or not reactions:
        raise ValueError("recipe.reactions must be a non-empty list")
    ids = [r.get("id") for r in reactions]
    if len(set(ids)) != len(ids):
        raise ValueError("reaction ids must be unique")

    group_ids = []
    group_total = {}
    group_index = []
    arene = []
    for reaction in reactions:
        a = _number(reaction.get("arene_M"), f"{reaction['id']}.arene_M", positive=True)
        group = reaction.get("arene_group", reaction["id"])
        if group not in group_total:
            group_total[group] = a
            group_ids.append(group)
        elif abs(group_total[group] - a) > 1e-12 * max(a, 1.0):
            raise ValueError(
                f"{reaction['id']}: reactions in arene_group '{group}' must declare the same arene_M"
            )
        group_index.append(group_ids.index(group))
        arene.append(a)
    group_total = np.array([group_total[g] for g in group_ids])
    group_index = np.array(group_index, dtype=int)

    dgs, ks = [], []
    for reaction in reactions:
        dg, k = _reaction_thermo(reaction, temperature, h)
        dgs.append(dg)
        ks.append(k)
    ks = np.array(ks)

    if "borane_M" in recipe:
        if "borane_equivalents" in recipe:
            raise ValueError("give either borane_M or borane_equivalents, not both")
        borane = _number(recipe["borane_M"], "borane_M", nonnegative=True)
    else:
        equivalents = _number(recipe.get("borane_equivalents", 1.0), "borane_equivalents", nonnegative=True)
        borane = equivalents * float(group_total.sum())

    if borane <= 0.0:
        extents = np.zeros(len(reactions))
    else:
        extents = _numeric_roots(group_total, group_index, borane, ks, gamma)
    limited = _check_equilibrium(extents, group_total, group_index, borane, gamma, ks)

    total = float(extents.sum())
    h2_liquid = total / gamma
    pressure = h2_liquid * h
    gas_mol = pressure * headspace / (R_BAR_L * temperature)
    dissolved_mol = h2_liquid * liquid_volume
    total_h2_mol = total * liquid_volume
    if total_h2_mol > 0.0:
        balance_error = abs(gas_mol + dissolved_mol - total_h2_mol) / total_h2_mol
    else:
        balance_error = 0.0
    if balance_error > 1e-8:
        raise ValueError("gas-liquid hydrogen balance failed")

    per_group = np.zeros(len(group_total))
    for i, g in enumerate(group_index):
        per_group[g] += extents[i]
    reactions_out = []
    for i, reaction in enumerate(reactions):
        g = group_index[i]
        reactions_out.append(
            {
                "id": reaction["id"],
                "arene_group": group_ids[g],
                "arene_M_initial": float(group_total[g]),
                "arene_M_final": float(max(0.0, group_total[g] - per_group[g])),
                "product_M": float(extents[i]),
                "conversion": float(extents[i] / group_total[g]),
                "deltaG_kJ_per_mol_at_T": float(dgs[i]),
                "deltaG_standard": reaction.get("deltaG_standard", "H2(aq,1M)"),
                "K": float(ks[i]),
            }
        )
    groups_out = [
        {
            "id": group_ids[g],
            "arene_M_initial": float(group_total[g]),
            "arene_M_final": float(max(0.0, group_total[g] - per_group[g])),
            "total_conversion": float(min(1.0, per_group[g] / group_total[g])),
        }
        for g in range(len(group_total))
    ]

    return {
        "temperature_K": temperature,
        "liquid_volume_L": liquid_volume,
        "headspace_volume_L": headspace,
        "henry": {
            "volatility_bar_L_per_mol": h,
            "reference_value_bar_L_per_mol": _number(
                henry["value_bar_L_per_mol"], "henry.value_bar_L_per_mol"
            ),
            "source": henry.get("source", "unspecified"),
        },
        "gas_liquid_factor_gamma": gamma,
        "borane_M_initial": borane,
        "borane_equivalents": borane / float(group_total.sum()),
        "reactions": reactions_out,
        "arene_groups": groups_out,
        "h2": {
            "liquid_M": h2_liquid,
            "partial_pressure_bar": pressure,
            "gas_mol": gas_mol,
            "dissolved_mol": dissolved_mol,
            "total_mol": total_h2_mol,
            "fraction_in_gas": (gas_mol / total_h2_mol) if total_h2_mol > 0.0 else 0.0,
        },
        "total_conversion": total / float(group_total.sum()),
        "stoichiometric_limit": [group_ids[g] for g in limited],
        "data_status": recipe.get("data_status", "unspecified"),
        "assumptions": [
            "ideal gas in the headspace; Henry law with constant volatility",
            "no activity coefficients; concentrations equal activities",
            "one H2 per aryl boronate; all reactions share one borane and one H2 pool",
            "reactions in one arene_group share that arene inventory",
            "equilibrium composition bounds, but does not predict, a finite-time yield",
        ],
    }


def solve(recipe):
    """Validate a recipe and return the equilibrium composition dict."""
    if not isinstance(recipe, dict):
        raise ValueError("recipe must be a JSON object")
    return _core(recipe)


def validate_sweeps(recipe):
    """Check the structure of recipe.sweeps without solving every point."""
    sweeps = recipe.get("sweeps") or {}
    if not isinstance(sweeps, dict):
        raise ValueError("recipe.sweeps must be an object")
    unknown = set(sweeps) - set(SWEEP_FIELDS)
    if unknown:
        raise ValueError(f"unsupported sweep fields: {sorted(unknown)}")
    total = 0
    for field, values in sweeps.items():
        if not isinstance(values, list) or not values:
            raise ValueError(f"sweeps.{field} must be a non-empty list")
        for value in values:
            if not isinstance(value, (int, float)) or isinstance(value, bool):
                raise ValueError(f"sweeps.{field} values must be numbers")
            if not math.isfinite(float(value)):
                raise ValueError(f"sweeps.{field} values must be finite")
            number = float(value)
            if field in ("temperature_K", "liquid_volume_L") and number <= 0:
                raise ValueError(f"sweeps.{field} must be positive")
            if field in ("headspace_volume_L", "borane_equivalents", "borane_M") and number < 0:
                raise ValueError(f"sweeps.{field} must be non-negative")
        total += len(values)
    return total


def sweep(recipe):
    """Return one row per sweep point; ``sweeps`` maps field -> list of values."""
    validate_sweeps(recipe)
    sweeps = recipe.get("sweeps") or {}
    rows = []
    for field, values in sweeps.items():
        for value in values:
            trial = copy.deepcopy(recipe)
            trial.pop("sweeps", None)
            trial[field] = value
            if field == "borane_M":
                trial.pop("borane_equivalents", None)
            elif field == "borane_equivalents":
                trial.pop("borane_M", None)
            try:
                result = _core(trial)
            except ValueError as exc:
                raise ValueError(f"sweep {field}={value}: {exc}") from None
            rows.append(
                {
                    "variable": field,
                    "value": float(value),
                    "total_conversion": result["total_conversion"],
                    "reactions": ";".join(
                        f"{r['id']}:{r['conversion']:.6g}" for r in result["reactions"]
                    ),
                    "partial_pressure_bar": result["h2"]["partial_pressure_bar"],
                    "fraction_in_gas": result["h2"]["fraction_in_gas"],
                    "data_status": result["data_status"],
                }
            )
    return rows


def run_equilibrium(recipe_path, dest):
    """Solve a recipe and write ``summary.json`` plus ``sweeps.csv`` to dest."""
    recipe_path = Path(recipe_path)
    dest = Path(dest)
    recipe = json.loads(recipe_path.read_text(encoding="utf-8"))
    summary = solve(recipe)
    rows = sweep(recipe)
    dest.mkdir(parents=True, exist_ok=True)
    payload = {
        "recipe_path": str(recipe_path),
        "data_status": summary["data_status"],
        "summary": summary,
        "sweeps": rows,
    }
    (dest / "summary.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    with (dest / "sweeps.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "variable",
                "value",
                "total_conversion",
                "reactions",
                "partial_pressure_bar",
                "fraction_in_gas",
                "data_status",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    return summary
