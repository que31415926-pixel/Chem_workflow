"""Explicit mass action in mol/L, hours. One litre is a normalization volume."""

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from .projects import ROOT as ROOT
from .projects import finite, folder, identifier, manifest, read_json


def write_table(path, columns, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def validate_network(data):
    names = set(data["initial_M"])
    if not names or set(data["species_labels"]) != names:
        raise ValueError(
            "Every network species requires an initial concentration and label"
        )
    for key in names:
        identifier(key)
    reserved = {"vial", "k", "time_h", "hypothesis", "condition"}
    if names & reserved:
        raise ValueError("Reserved species name: " + str(names & reserved))
    if any(not finite(v) or v < 0 for v in data["initial_M"].values()):
        raise ValueError("Invalid initial concentrations")
    rx = data["reactions"]
    ids = {identifier(r["id"]) for r in rx}
    if not rx or len(ids) != len(rx):
        raise ValueError("Reaction IDs must be unique")
    for r in rx:
        if not finite(r["k"]) or r["k"] < 0:
            raise ValueError(f"{r['id']}: invalid rate")
        for side in ("lhs", "rhs"):
            if set(r[side]) - names or any(
                type(v) is not int or v <= 0 for v in r[side].values()
            ):
                raise ValueError(f"{r['id']}: invalid species or stoichiometry")
    for law, weights in data["conserved"].items():
        if (
            set(weights) - names
            or not weights
            or any(not finite(v) for v in weights.values())
        ):
            raise ValueError(f"{law}: invalid conservation coefficients")

    def overrides(mapping, allowed):
        if set(mapping) - allowed or any(
            not finite(v) or v < 0 for v in mapping.values()
        ):
            raise ValueError(f"Unknown key or invalid nonnegative value: {mapping}")

    for name, cond in data["conditions"].items():
        identifier(name)
        overrides(cond.get("initial_M", {}), names)
        for key in ("rate_values", "rate_factors"):
            overrides(cond.get(key, {}), ids)
    if not data["hypotheses"] or data["baseline"] not in data["conditions"]:
        raise ValueError("Provide hypotheses and name an existing baseline condition")
    for name, h in data["hypotheses"].items():
        identifier(name)
        overrides(h["rates"], ids)
        keys = h.get("fit", [])
        keys = [keys] if isinstance(keys, str) else keys
        if not isinstance(keys, list) or len(set(keys)) != len(keys) or set(keys) - ids:
            raise ValueError(f"{name}: unknown fit parameter")
    end = data["end_h"]
    times = data["times_h"]
    if (
        not finite(end)
        or end <= 0
        or not times
        or any(not finite(t) or not 0 <= t <= end for t in times)
        or any(a >= b for a, b in zip(times, times[1:]))
    ):
        raise ValueError("Invalid time grid or observation window")
    for p in data.get("panels", []):
        if (
            p["condition"] not in data["conditions"]
            or p.get("baseline", data["baseline"]) not in data["conditions"]
            or p["observable"] not in names
            or not finite(p["time_h"])
            or not 0 <= p["time_h"] <= end
            or not finite(p["sd_M"])
            or p["sd_M"] <= 0
        ):
            raise ValueError(f"Invalid diagnostic panel: {p}")


class Network:
    def __init__(self, project, calibrated=None):
        self.project = project
        manifest(project)
        self.data = read_json(folder(project) / "network.json")
        self.fingerprint = hashlib.sha256(
            json.dumps(self.data, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        validate_network(self.data)
        self.baseline = self.data["baseline"]
        self.names = list(self.data["initial_M"])
        self.rx = self.data["reactions"]
        self.hypotheses = self.data["hypotheses"]
        self.conditions = self.data["conditions"]
        self.provenance = {
            "origin": "illustrative_scenario",
            "parameter_status": "unfitted",
        }
        if calibrated:
            fitted = read_json(calibrated)
            if fitted["project"] != project or set(fitted["parameters"]) != set(
                self.hypotheses
            ):
                raise ValueError(
                    "Calibration belongs to a different project or hypotheses"
                )
            if "network_sha256" not in fitted:
                raise ValueError(
                    "Calibration has no network_sha256; refit with this version"
                )
            self.provenance = fitted.get("provenance")
            if not isinstance(self.provenance, dict):
                raise ValueError("Calibration has no provenance block; refit")
            if not self.provenance.get("origin") or not self.provenance.get("data_sha256"):
                raise ValueError(
                    "Calibration provenance lacks origin/data_sha256; refit"
                )
            if fitted["network_sha256"] != self.fingerprint:
                raise ValueError(
                    "Calibration belongs to different network settings; refit this project"
                )
            self.calibrated = fitted["parameters"]
            for overrides in self.calibrated.values():
                if set(overrides) - {r["id"] for r in self.rx} or not all(
                    finite(v) and v >= 0 for v in overrides.values()
                ):
                    raise ValueError("Unknown or invalid calibrated parameter")
        else:
            self.calibrated = {}
        self.stoich = np.zeros((len(self.names), len(self.rx)))
        self.orders = []
        for j, reaction in enumerate(self.rx):
            self.orders.append(
                [(self.names.index(n), k) for n, k in reaction["lhs"].items()]
            )
            for side, sign in [("lhs", -1), ("rhs", 1)]:
                for n, k in reaction[side].items():
                    self.stoich[self.names.index(n), j] += sign * k
        self.balances = np.array(
            [[v.get(n, 0) for n in self.names] for v in self.data["conserved"].values()]
        )
        self.balances = self.balances.reshape(-1, len(self.names))
        if not np.allclose(self.balances @ self.stoich, 0):
            raise ValueError("Unbalanced reaction in " + project)

    def settings(self, hyp, cond, overrides=None):
        initial = dict(self.data["initial_M"])
        initial.update(self.conditions[cond].get("initial_M", {}))
        rates = {r["id"]: r["k"] for r in self.rx}
        rates.update(self.hypotheses[hyp]["rates"])
        rates.update(self.calibrated.get(hyp, {}))
        rates.update(overrides or {})
        rates.update(self.conditions[cond].get("rate_values", {}))
        for key, multiple in self.conditions[cond].get("rate_factors", {}).items():
            rates[key] *= multiple
        for values in (rates, initial):
            if any(not math.isfinite(x) or x < 0 for x in values.values()):
                raise ValueError("Non-finite or negative setting")
        return initial, rates

    def simulate(self, hyp, cond, times, overrides=None):
        times = np.asarray(times, dtype=float)
        if (
            times.ndim != 1
            or not len(times)
            or not np.isfinite(times).all()
            or np.any(np.diff(times) <= 0)
            or times[0] < 0
            or times[-1] > self.data["end_h"]
        ):
            raise ValueError(
                "Times must increase and lie inside the project observation window"
            )
        init, rates = self.settings(hyp, cond, overrides)
        y0 = np.array([init[n] for n in self.names])
        k = np.array([rates[r["id"]] for r in self.rx])

        def rhs(_, y):
            v = k.copy()
            for j, entries in enumerate(self.orders):
                for i, order in entries:
                    v[j] *= max(y[i], 0.0) ** order
            return self.stoich @ v

        if times[-1] == 0:
            return y0[:, None]
        sol = solve_ivp(
            rhs,
            (0.0, times[-1]),
            y0,
            t_eval=times,
            method="LSODA",
            rtol=3e-8,
            atol=2e-12,
        )
        if not sol.success or not np.isfinite(sol.y).all() or sol.y.min() < -1e-8:
            raise ValueError("Integration failed or negative pool")
        if not np.allclose(
            self.balances @ sol.y, (self.balances @ y0)[:, None], atol=2e-8, rtol=2e-6
        ):
            raise ValueError("Pool balance failed")
        return np.maximum(sol.y, 0.0)

    def export_sbml(self, hyp, cond, path):
        import libsbml as sb

        doc = sb.SBMLDocument(3, 2)
        model = doc.createModel()
        model.setId(self.project + "_" + hyp + "_" + cond)
        model.setName(
            "Scenario model; inspect parameter provenance before interpretation"
        )
        model.setTimeUnits("hour")
        model.setSubstanceUnits("mole")
        model.setExtentUnits("mole")
        model.setVolumeUnits("litre")

        def unit_definition(name, components):
            ud = model.createUnitDefinition()
            ud.setId(name)
            for kind, exponent, multiplier in components:
                unit = ud.createUnit()
                unit.setKind(kind)
                unit.setExponent(exponent)
                unit.setScale(0)
                unit.setMultiplier(multiplier)

        unit_definition("hour", [(sb.UNIT_KIND_SECOND, 1, 3600)])
        for order in {sum(r["lhs"].values()) for r in self.rx}:
            units = [(sb.UNIT_KIND_SECOND, -1, 3600)]
            if order != 1:
                units += [
                    (sb.UNIT_KIND_LITRE, order - 1, 1),
                    (sb.UNIT_KIND_MOLE, 1 - order, 1),
                ]
            unit_definition("rate_order_" + str(order), units)
        comp = model.createCompartment()
        comp.setId("vial")
        comp.setSpatialDimensions(3)
        comp.setSize(1.0)
        comp.setUnits("litre")
        comp.setConstant(True)
        init, rates = self.settings(hyp, cond)
        for name, value in init.items():
            s = model.createSpecies()
            s.setId(name)
            s.setName(self.data["species_labels"][name])
            s.setCompartment("vial")
            s.setSubstanceUnits("mole")
            s.setInitialConcentration(value)
            s.setHasOnlySubstanceUnits(False)
            s.setBoundaryCondition(False)
            s.setConstant(False)
        for record in self.rx:
            rx = model.createReaction()
            rx.setId("rx_" + record["id"])
            rx.setReversible(
                False
            )  # reversible steps are explicit forward/reverse pairs
            for side, factory in [
                ("lhs", rx.createReactant),
                ("rhs", rx.createProduct),
            ]:
                for name, count in record[side].items():
                    term = factory()
                    term.setSpecies(name)
                    term.setStoichiometry(count)
                    term.setConstant(True)
            law = rx.createKineticLaw()
            p = law.createLocalParameter()
            p.setId("k")
            p.setValue(rates[record["id"]])
            p.setUnits("rate_order_" + str(sum(record["lhs"].values())))
            factors = ["vial", "k"] + [f"{n}^{s}" for n, s in record["lhs"].items()]
            law.setMath(sb.parseL3Formula(" * ".join(factors)))
        doc.checkConsistency()
        errors = [
            doc.getError(i).getMessage()
            for i in range(doc.getNumErrors())
            if doc.getError(i).getSeverity() >= sb.LIBSBML_SEV_ERROR
        ]
        if errors:
            raise ValueError("SBML invalid: " + " ".join(errors))
        if not sb.writeSBMLToFile(doc, str(path)):
            raise OSError(f"Could not write {path}")


def run_scenarios(net, dest):
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "provenance.json").write_text(
        json.dumps(
            {
                "project": net.project,
                **net.provenance,
                "note": "Predictions from effective pools, not measured yields. Diagnostic scores assume the specified measurement uncertainty.",
            },
            indent=2,
        )
        + "\n"
    )
    times = net.data["times_h"]
    rows, panels = [], []
    hyps = list(net.hypotheses)
    for h in hyps:
        for c in net.conditions:
            values = net.simulate(h, c, times)
            for i, t in enumerate(times):
                rows.append(
                    dict(
                        hypothesis=h,
                        condition=c,
                        time_h=t,
                        **dict(zip(net.names, values[:, i])),
                    )
                )
            net.export_sbml(h, c, dest / f"{h}__{c}.xml")
    for panel in net.data.get("panels", []):
        c = panel["condition"]
        predictions = []
        for h in hyps:
            obs = panel["observable"]
            t = [panel["time_h"]]
            i = net.names.index(obs)
            perturbed = net.simulate(h, c, t)[i, 0]
            baseline = net.simulate(h, panel.get("baseline", net.baseline), t)[i, 0]
            predictions.append(perturbed - baseline)
        for a in range(len(hyps)):
            for b in range(a + 1, len(hyps)):
                panels.append(
                    dict(
                        condition=c,
                        observable=obs,
                        time_h=t[0],
                        A=hyps[a],
                        B=hyps[b],
                        change_A_M=predictions[a],
                        change_B_M=predictions[b],
                        assumed_sd_of_contrast_M=panel["sd_M"],
                        separation=abs(predictions[a] - predictions[b]) / panel["sd_M"],
                    )
                )
    write_table(
        dest / "trajectories_M.csv",
        ["hypothesis", "condition", "time_h", *net.names],
        rows,
    )
    write_table(
        dest / "diagnostics.csv",
        [
            "condition",
            "observable",
            "time_h",
            "A",
            "B",
            "change_A_M",
            "change_B_M",
            "assumed_sd_of_contrast_M",
            "separation",
        ],
        sorted(panels, key=lambda x: -x["separation"]),
    )
    return len(rows)
