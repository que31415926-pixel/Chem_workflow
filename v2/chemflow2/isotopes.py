"""Compare matched predictions for isotope, hydrogen pressure and deuterium tracing."""
import math


def _number(value, name, positive=False):
    if type(value) not in (int, float) or not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f"invalid {name}")
    return value


def evaluate(data):
    """Evaluate apparent orders by finite differences in a reviewed pressure range."""
    rows = {}
    for entry in data["explanations"]:
        label = entry["id"]
        if not isinstance(label, str) or not label or label in rows:
            raise ValueError("distinct competing mechanism identifiers required")
        condition = entry["condition_id"]
        if not isinstance(condition, str) or not condition:
            raise ValueError("matched condition required")
        rates = entry["isotopic_rates"]
        if set(rates) != {"H", "D"}:
            raise ValueError("matched hydrogen and deuterium rates required")
        rates = {k: _number(v, k, positive=True) for k, v in rates.items()}
        if not isinstance(entry["initial_rate_unit"], str) or not entry["initial_rate_unit"]:
            raise ValueError("declare common rate units")
        points = []
        for point in entry["pressure_series"]:
            p = _number(point["hydrogen_pressure_bar"], "pressure", positive=True)
            r = _number(point["initial_rate"], "initial rate", positive=True)
            points.append((p, r))
        points.sort()
        if len(points) < 2 or len({p for p, _ in points}) != len(points):
            raise ValueError("at least two distinct positive pressures required")
        orders = [{"pressure_bar": [p1, p2], "order": (math.log(r2)-math.log(r1))/(math.log(p2)-math.log(p1))}
                  for (p1, r1), (p2, r2) in zip(points, points[1:])]
        sites = {k: _number(v, "site deuterium mol") for k, v in entry["site_deuterium_mol"].items()}
        inventory = _number(entry["deuterium_input_mol"], "deuterium inventory", positive=True)
        if not sites or any(v < 0 for v in sites.values()) or not math.isclose(
                math.fsum(sites.values()), inventory, rel_tol=1e-7, abs_tol=1e-14):
            raise ValueError("site inventory must include unreacted and unobserved deuterium pools")
        rows[label] = {"condition_id": condition, "initial_rate_unit": entry["initial_rate_unit"],
                       "kinetic_isotope_effect": rates["H"]/rates["D"], "hydrogen_orders": orders,
                       "deuterium_site_fractions": {k: v/inventory for k, v in sites.items()}}
    if len(rows) < 2:
        raise ValueError("at least two competing explanations required")
    return {"operation": "isotope_observables", "explanations": rows}
