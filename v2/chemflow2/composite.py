"""Component-resolved CBS and composite single-conformer thermochemistry on pinned inputs."""
import hashlib
import json
import math
import re


_BASIS_PAIR = {"def2-VPP": ("def2-TZVPP", "def2-QZVPP"),
               "def2-VP": ("def2-TZVP", "def2-QZVP"),
               "cc-pVXZ": ("cc-pVTZ", "cc-pVQZ"),
               "aug-cc-pVXZ": ("aug-cc-pVTZ", "aug-cc-pVQZ")}
_R = 8.31446261815324       # J mol^-1 K^-1
_HARTREE_J_MOL = 2625499.6394799
_PRESSURE_PA = {"1_atm": 101325.0, "1_bar": 100000.0}


def _text(d, key):
    val = d.get(key)
    if not isinstance(val, str) or not val.strip():
        raise ValueError(f"nonempty {key} required")
    return val


def _sha(d, key):
    val = _text(d, key)
    if not re.fullmatch("[a-f0-9]{64}", val):
        raise ValueError(f"{key} must be a lowercase SHA-256")
    return val


def _object(d, key):
    val = d.get(key)
    if not isinstance(val, dict):
        raise ValueError(f"object required for {key}")
    return val


def _number(d, key, positive=False):
    val = d.get(key)
    if type(val) not in (int, float) or not math.isfinite(val) or (positive and val <= 0):
        raise ValueError(f"invalid numeric {key}")
    return float(val)


def _record(d):
    if d.get("unit") != "Eh":
        raise ValueError("energy unit must be Eh")
    if type(d.get("charge")) is not int or type(d.get("multiplicity")) is not int or d["multiplicity"] < 1:
        raise ValueError("electronic charge and multiplicity are required")
    identity = {key: _sha(d, key) for key in ("geometry_sha256", "medium_sha256")}
    identity.update({key: d[key] for key in ("charge", "multiplicity")})
    identity["state_id"] = _text(d, "state_id")
    source = {"task_id": _text(d, "task_id"), "artifact_sha256": _sha(d, "artifact_sha256"),
              "method_settings_sha256": _sha(d, "method_settings_sha256")}
    return identity, source


def _paired(left, right):
    li, ls = _record(left)
    ri, rs = _record(right)
    if li != ri:
        raise ValueError("inputs differ in geometry, medium, charge, multiplicity or electronic state")
    if ls["task_id"] == rs["task_id"] or ls["artifact_sha256"] == rs["artifact_sha256"]:
        raise ValueError("distinct tasks and source artifacts are required")
    return li, ls, rs


def _limit(smaller, larger, log_ratio):
    # The cardinal 4 tail relative to the 3–4 difference. This form is stable at large exponents.
    delta = -math.expm1(-log_ratio)
    if not delta:
        raise ValueError("CBS exponent prevents a finite solution")
    value = larger + (larger - smaller) * math.exp(-log_ratio) / delta
    if not math.isfinite(value):
        raise ValueError("non-finite CBS solution")
    return value


def model_digest(model):
    """Return the digest to freeze alongside a prospectively selected CBS model."""
    raw = json.dumps(model, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _cbs(data):
    model = _object(data, "model")
    frozen = _sha(data, "frozen_model_sha256")
    if model_digest(model) != frozen:
        raise ValueError("CBS model differs from the frozen definition")
    family = _text(model, "basis_family")
    if family not in _BASIS_PAIR:
        raise ValueError("unrecognised TZ/QZ basis family")
    if model.get("scf_model") != "exp_sqrt_cardinal" or model.get("correlation_model") != "inverse_power":
        raise ValueError("unsupported SCF or correlation extrapolation model")
    _text(model, "model_id")
    _text(model, "parameter_source")
    alpha = _number(model, "scf_alpha", positive=True)
    beta = _number(model, "correlation_beta", positive=True)
    points = [_object(data, key) for key in ("tz", "qz")]
    identity, low, high = _paired(*points)
    if low["method_settings_sha256"] != high["method_settings_sha256"]:
        raise ValueError("CBS requires identical non-basis method settings")
    for point, basis, cardinal in zip(points, _BASIS_PAIR[family], (3, 4)):
        if (point.get("basis") != basis or type(point.get("cardinal")) is not int
                or point["cardinal"] != cardinal or point.get("quantity") != "electronic_components"):
            raise ValueError(f"expected {basis}, cardinal {cardinal}, and electronic_components")
        scf, corr, total = (_number(point, key) for key in ("scf_eh", "correlation_eh", "electronic_eh"))
        # Reconciliation tolerance for printed output (roughly 6e-5 kcal/mol), not a research gate.
        if not math.isclose(scf + corr, total, abs_tol=1e-7, rel_tol=0):
            raise ValueError("reported total does not match SCF plus correlation; inspect the partition")
    scf = _limit(points[0]["scf_eh"], points[1]["scf_eh"], alpha * (2 - math.sqrt(3)))
    corr = _limit(points[0]["correlation_eh"], points[1]["correlation_eh"], beta * math.log(4 / 3))
    if not math.isfinite(scf + corr):
        raise ValueError("non-finite CBS total")
    return {**identity, "quantity": "electronic_energy", "unit": "Eh",
            "electronic_eh": scf + corr, "scf_cbs_eh": scf, "correlation_cbs_eh": corr,
            "model": dict(model), "frozen_model_sha256": frozen, "source_records": {"tz": low, "qz": high}}


def _thermo(data):
    high, low = (_object(data, key) for key in ("high_level", "frequency"))
    identity, high_prov, low_prov = _paired(high, low)
    if high.get("quantity") != "electronic_energy" or low.get("quantity") != "qrrho_gibbs_and_electronic":
        raise ValueError("requires a pure high-level electronic energy and matched qRRHO G/Eel")
    high_e, low_e, low_g = (_number(high, "electronic_eh"), _number(low, "electronic_eh"),
                            _number(low, "qrrho_gibbs_eh"))
    temp = _number(low, "temperature_k", positive=True)
    model = _text(low, "qrrho_model_id")
    reviewed_by = _text(low, "reviewed_by")
    ss = _object(data, "standard_state_policy")
    if set(ss) != {"source", "target", "mode"}:
        raise ValueError("use only source, target, mode in a standard-state policy")
    reported = _text(low, "reported_standard_state")
    target, mode = _text(ss, "target"), _text(ss, "mode")
    if reported not in {*_PRESSURE_PA, "1_M"}:
        raise ValueError("unsupported reported standard state")
    if ss.get("source") != reported:
        raise ValueError("policy source conflicts with the state already applied to the declared qRRHO G")
    if mode == "as_reported" and reported == target:
        standard_correction = 0.0
    elif (mode == "ideal_gas_pressure_to_1M" and reported in _PRESSURE_PA and target == "1_M"):
        ratio = 1000 * _R * temp / _PRESSURE_PA[reported]
        standard_correction = _R * temp * math.log(ratio) / _HARTREE_J_MOL
    else:
        raise ValueError("unsupported, repeated, or inconsistent standard-state correction")
    thermal_correction = low_g - low_e
    value = math.fsum([high_e, thermal_correction, standard_correction])
    if not math.isfinite(value):
        raise ValueError("non-finite composite Gibbs energy")
    return {**identity, "quantity": "single_conformer_gibbs_energy", "unit": "Eh",
            "gibbs_eh": value, "electronic_eh": high_e, "thermal_qrrho_correction_eh": thermal_correction,
            "per_species_standard_state_correction_eh": standard_correction, "temperature_k": temp,
            "standard_state": target, "standard_state_policy": dict(ss),
            "qrrho_model_id": model, "thermochemistry_reviewed_by": reviewed_by,
            "source_records": {"high_level": high_prov, "frequency": low_prov}}


def evaluate(data):
    protocol = _text(data, "protocol_id")
    op = _text(data, "operation")
    if op == "two_point_cbs":
        values = _cbs(data)
    elif op == "composite_thermochemistry":
        values = _thermo(data)
    else:
        raise ValueError(f"unknown composite operation: {op}")
    return {"operation": op, "protocol_id": protocol, "status": "pass",
            "scope": "declared_input_numerical_assembly", **values}
