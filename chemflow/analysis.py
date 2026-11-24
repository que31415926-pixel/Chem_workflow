"""Fit one tractable rate per hypothesis and score independent observed conditions."""

import csv
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from .network import write_table


def observations(path, net):
    with Path(path).open(encoding="utf-8-sig", newline="") as src:
        rows = list(csv.DictReader(src))
    if not rows:
        raise ValueError("Empty observations file")
    keys = set()
    for row in rows:
        if (
            row["split"] not in ("train", "holdout")
            or row["condition"] not in net.conditions
            or row["species"] not in net.names
            or row["kind"] not in ("value", "upper")
            or not row["batch"].strip()
        ):
            raise ValueError(
                "Invalid condition, pool, batch, split or observation kind"
            )
        for col in ("time_h", "value_M", "sigma_M"):
            row[col] = float(row[col])
            if not math.isfinite(row[col]) or row[col] < 0:
                raise ValueError("Bad value in " + col)
        if row["time_h"] > net.data["end_h"] or row["sigma_M"] <= 0:
            raise ValueError("Invalid time or uncertainty")
        key = tuple(row[col] for col in ["batch", "condition", "time_h", "species"])
        if key in keys:
            raise ValueError(
                "Duplicate measurement in one batch; combine technical injections first"
            )
        keys.add(key)
    groups = {
        s: {r["condition"] for r in rows if r["split"] == s}
        for s in ("train", "holdout")
    }
    if not all(groups.values()) or groups["train"] & groups["holdout"]:
        raise ValueError(
            "Reserve at least one whole unused condition; no train/holdout condition overlap"
        )
    return rows


def residuals(net, h, rows, key, value):
    predictions = {}
    for c in sorted({r["condition"] for r in rows}):
        times = sorted({r["time_h"] for r in rows if r["condition"] == c})
        vals = net.simulate(h, c, times, {key: value})
        for j, t in enumerate(times):
            predictions[c, t] = dict(zip(net.names, vals[:, j]))
    result = []
    for r in rows:
        pred = predictions[r["condition"], r["time_h"]][r["species"]]
        err = pred - r["value_M"]
        if r["kind"] == "upper":
            err = max(0.0, err)
        result.append(
            dict(
                **r,
                hypothesis=h,
                predicted_M=pred,
                weighted_residual=err / r["sigma_M"],
            )
        )
    return result


def fit(net, path, dest):
    import hashlib
    import json

    rows = observations(path, net)
    train = [r for r in rows if r["split"] == "train"]
    holdout = [r for r in rows if r["split"] == "holdout"]
    scores = []
    all_rows = []
    params = {}
    for h, definition in net.hypotheses.items():
        key = definition["fit"]
        start = net.settings(h, "base")[1][key]

        def fun(log_factor):
            attempt = start * np.exp(log_factor[0])
            return [
                r["weighted_residual"] for r in residuals(net, h, train, key, attempt)
            ]

        fitted = least_squares(fun, [0.0], bounds=([-np.log(100.0)], [np.log(100.0)]))
        if not fitted.success:
            raise ValueError("Fit failed for " + h)
        value = float(start * np.exp(fitted.x[0]))
        params[h] = {key: value}
        for split, subset in [("train", train), ("holdout", holdout)]:
            pred = residuals(net, h, subset, key, value)
            errors = [r["weighted_residual"] for r in pred]
            scores.append(
                dict(
                    hypothesis=h,
                    fit_parameter=key,
                    fitted_value=value,
                    split=split,
                    n=len(errors),
                    at_search_bound=bool(np.any(fitted.active_mask)),
                    weighted_rms=float(np.sqrt(np.mean(np.square(errors)))),
                )
            )
            all_rows.extend(pred)
    dest.mkdir(parents=True, exist_ok=True)
    write_table(dest / "scores.csv", scores[0].keys(), scores)
    write_table(dest / "residuals.csv", all_rows[0].keys(), all_rows)
    (dest / "calibrated.json").write_text(
        json.dumps(
            {
                "project": net.project,
                "parameters": params,
                "provenance": {
                    "origin": sorted({r.get("origin", "unspecified") for r in rows}),
                    "data_sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                    "parameter_status": "one_param_per_hypothesis",
                },
                "note": "One fitted rate per hypothesis; all other parameters retain their scenario values.",
            },
            indent=2,
        )
        + "\n"
    )
