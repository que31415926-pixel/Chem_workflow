"""Account for measured material; a blank input is unknown, never zero."""

import csv
import math
from pathlib import Path

from .network import write_table

SCHEMAS = {
    "co": [
        "arene_in",
        "hbpin_in",
        "b2_in",
        "co_atoms_in",
        "target",
        "arene_left",
        "other_aryl_mono",
        "aryl_di",
        "hbpin_left",
        "b2_left",
        "other_b_atoms",
        "h2_gas",
        "h2_liquid",
        "blank_h2_total",
    ],
    "fe": [
        "salt_effective_in",
        "acid_in",
        "target",
        "salt_left",
        "benzyl_dimer",
        "benzyl_side_eq",
        "acid_left",
        "tert_dimer",
        "tert_side_eq",
    ],
}


def account(project, source, dest):
    with Path(source).open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        if not set(["batch", *SCHEMAS[project]]).issubset(reader.fieldnames or []):
            raise ValueError("Missing material-balance columns")
        rows = list(reader)
    if not rows:
        raise ValueError("Material-balance table is empty")
    output = []
    for row in rows:
        v = {k: float(row[k]) if row[k].strip() else None for k in SCHEMAS[project]}
        if any(a is not None and (not math.isfinite(a) or a < 0) for a in v.values()):
            raise ValueError("Amounts must be finite and nonnegative")

        def total(**weights):
            if any(v[k] is None for k in weights):
                return None
            return sum(v[k] * n for k, n in weights.items())

        def ratio(a, b, scale=1.0):
            return scale * a / b if a is not None and b and b > 0 else None

        p = v["target"]
        if project == "co":
            supplied = total(hbpin_in=1, b2_in=2)
            gas = total(h2_gas=1, h2_liquid=1, blank_h2_total=-1)
            result = dict(
                target_yield_pct=ratio(p, v["arene_in"], 100),
                TON_total_Co=ratio(p, v["co_atoms_in"]),
                target_B_efficiency_pct=ratio(p, supplied, 100),
                assigned_B_pct=ratio(
                    total(
                        target=1,
                        other_aryl_mono=1,
                        aryl_di=2,
                        hbpin_left=1,
                        b2_left=2,
                        other_b_atoms=1,
                    ),
                    supplied,
                    100,
                ),
                assigned_arene_pct=ratio(
                    total(target=1, other_aryl_mono=1, aryl_di=1, arene_left=1),
                    v["arene_in"],
                    100,
                ),
                blank_corrected_H2_per_target=ratio(gas, p),
            )
        else:
            result = dict(
                target_yield_on_acid_pct=ratio(p, v["acid_in"], 100),
                assigned_benzyl_pct=ratio(
                    total(target=1, salt_left=1, benzyl_dimer=2, benzyl_side_eq=1),
                    v["salt_effective_in"],
                    100,
                ),
                assigned_tert_pct=ratio(
                    total(target=1, acid_left=1, tert_dimer=2, tert_side_eq=1),
                    v["acid_in"],
                    100,
                ),
            )
        output.append(
            dict(
                batch=row["batch"],
                **result,
                unmeasured=" ".join(k for k, x in v.items() if x is None),
            )
        )
    dest.mkdir(parents=True, exist_ok=True)
    write_table(dest / "material_balance.csv", output[0].keys(), output)
