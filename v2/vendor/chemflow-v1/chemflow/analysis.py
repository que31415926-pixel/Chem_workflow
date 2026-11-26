"""Fit configured rates per hypothesis and score independent observed conditions."""

import csv
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from .network import write_table

# The default finite-difference step (~1e-8 relative) is at the same scale as the
# LSODA residual noise (rtol 3e-8). For the Fe network that produced a spurious
# local minimum (salt_loss 1.037 instead of 1.530 h-1). A 1e-6 relative step
# resolves the gradient with negligible truncation error.
FIT_DIFF_STEP = 1e-6


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


def residuals(net, h, rows, rates):
    predictions = {}
    for c in sorted({r["condition"] for r in rows}):
        times = sorted({r["time_h"] for r in rows if r["condition"] == c})
        vals = net.simulate(h, c, times, rates)
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
    diagnostics = {}
    for h, definition in net.hypotheses.items():
        requested = definition.get("fit")
        if not requested:
            raise ValueError(f"{h}: declare one or more rate IDs to fit")
        keys = [requested] if isinstance(requested, str) else requested
        base = {r["id"]: r["k"] for r in net.rx}
        base.update(definition["rates"])
        starts = np.array([base[k] for k in keys], dtype=float)
        if np.any(starts <= 0):
            raise ValueError(f"{h}: fitted rates must start positive")
        factor = definition.get("fit_factor", 100.0)
        if not math.isfinite(factor) or factor <= 1:
            raise ValueError("fit_factor must exceed 1")

        def rates(x):
            return {k: float(v) for k, v in zip(keys, starts * np.exp(x), strict=True)}

        def fun(x):
            return [r["weighted_residual"] for r in residuals(net, h, train, rates(x))]

        fitted = least_squares(
            fun,
            np.zeros(len(keys)),
            bounds=(-np.log(factor), np.log(factor)),
            diff_step=FIT_DIFF_STEP,
        )
        if not fitted.success:
            raise ValueError("Fit failed for " + h)
        params[h] = rates(fitted.x)
        diagnostics[h] = {
            "nfev": int(fitted.nfev),
            "status": int(fitted.status),
            "optimality": float(fitted.optimality),
            "message": str(fitted.message),
        }
        for split, subset in [("train", train), ("holdout", holdout)]:
            pred = residuals(net, h, subset, params[h])
            errors = [r["weighted_residual"] for r in pred]
            for i, key in enumerate(keys):
                scores.append(
                    dict(
                        hypothesis=h,
                        fit_parameter=key,
                        fitted_value=params[h][key],
                        split=split,
                        n=len(errors),
                        at_search_bound=bool(fitted.active_mask[i]),
                        weighted_rms=float(np.sqrt(np.mean(np.square(errors)))),
                        nfev=int(fitted.nfev),
                        optimality=float(fitted.optimality),
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
                "network_sha256": net.fingerprint,
                "provenance": {
                    "origin": sorted({r.get("origin", "unspecified") for r in rows}),
                    "data_sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                    "parameter_status": "configured_rates_per_hypothesis",
                    "fit_settings": {
                        "optimizer": "scipy.optimize.least_squares",
                        "method": "trf",
                        "diff_step": FIT_DIFF_STEP,
                        "starts": 1,
                    },
                    "fit_diagnostics": diagnostics,
                },
                "note": "Only the configured rates are fitted; other rates retain their scenario values.",
            },
            indent=2,
        )
        + "\n"
    )
