"""Evaluate declared ratios of weighted sums. Missing measurements propagate."""

import csv
from pathlib import Path

from .network import write_table
from .projects import finite, folder, read_json


def specification(project):
    spec = read_json(folder(project) / "material.json")
    names = spec["fields"]
    metrics = spec["metrics"]
    if (
        not names
        or len(set(names)) != len(names)
        or not metrics
        or set(metrics) & {"batch", "unmeasured"}
    ):
        raise ValueError(
            "Material fields must be unique and metric names distinct from metadata"
        )
    for name, metric in metrics.items():
        for side in ("numerator", "denominator"):
            weights = metric[side]
            if (
                not weights
                or set(weights) - set(names)
                or not all(finite(w) and w != 0 for w in weights.values())
            ):
                raise ValueError(
                    f"{name}.{side}: use defined fields and finite nonzero weights"
                )
        if not finite(metric.get("scale", 1)):
            raise ValueError(f"{name}.scale must be finite")
    return spec


def account(project, source, dest):
    spec = specification(project)
    names = spec["fields"]
    with Path(source).open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        if not set(["batch", *names]).issubset(reader.fieldnames or []):
            raise ValueError("Missing material-balance columns")
        rows = list(reader)
    if not rows:
        raise ValueError("Material-balance table is empty")
    results = []
    for row in rows:
        values = {
            k: float(row[k]) if row[k] and row[k].strip() else None for k in names
        }
        if any(x is not None and (not finite(x) or x < 0) for x in values.values()):
            raise ValueError("Amounts must be finite and nonnegative")

        def total(weights):
            if any(values[n] is None for n in weights):
                return None
            return sum(values[n] * w for n, w in weights.items())

        out = {"batch": row["batch"]}
        for name, metric in spec["metrics"].items():
            a, b = total(metric["numerator"]), total(metric["denominator"])
            out[name] = (
                metric.get("scale", 1) * a / b
                if a is not None and b is not None and b > 0
                else None
            )
        out["unmeasured"] = " ".join(k for k, v in values.items() if v is None)
        results.append(out)
    dest.mkdir(parents=True, exist_ok=True)
    write_table(dest / "material_balance.csv", results[0].keys(), results)
