"""Strictly separated calibration and validation for a preregistered scalar classifier."""

from .mechanism import _digest, _hash, _number, _text


def _registration(data):
    registry = data["registration"]
    if _hash(registry) != _digest(data["registration_sha256"], "registration hash"):
        raise ValueError("preregistered anchor specification differs from commitment")
    _text(registry["score_id"], "score identity")
    _text(registry["unit"], "score unit")
    _digest(registry["score_protocol_sha256"], "calculation and score definition")
    _text(registry["limitations"], "qualitative evidence and differing-condition limitations")
    if registry["higher_scores"] not in {"stable", "unstable"}:
        raise ValueError("freeze the direction before calibration")
    if registry["threshold_rule"] != "midpoint_of_strict_class_gap":
        raise ValueError("unsupported frozen calibration rule")
    rows = registry["anchors"]
    if not isinstance(rows, dict) or len(rows) < 3:
        raise ValueError("an anchor registry with separate calibration and validation sets is required")
    groups, ids = {"train": set(), "heldout": set()}, {"train": set(), "heldout": set()}
    for key, item in rows.items():
        _text(key, "anchor identity")
        role = item["split"]
        if role not in ids:
            raise ValueError("invalid split assignment")
        group = _text(item["group_id"], "compound or equivalent evidence group")
        _text(item["condition_id"], "experimental conditions")
        _text(item["source_locator"], "traceable original citation")
        _digest(item["source_sha256"], "source file")
        _digest(item["observation_sha256"], "frozen observation record")
        groups[role].add(group)
        ids[role].add(key)
    if not all(ids.values()) or abs(len(ids["train"]) - len(ids["heldout"])) > 1:
        raise ValueError("anchor rows must be split into two nonempty near-equal sets")
    if groups["train"] & groups["heldout"]:
        raise ValueError("compound or equivalent evidence groups overlap the split")
    return registry, ids


def _rows(registry, supplied, allowed):
    if not isinstance(supplied, dict) or supplied.keys() - allowed:
        raise ValueError("unexpected records outside this split")
    result = {}
    direction = 1 if registry["higher_scores"] == "stable" else -1
    for key, row in supplied.items():
        frozen = registry["anchors"][key]
        observation, score = row["observation"], row["score"]
        if (_hash(observation) != frozen["observation_sha256"] or
                observation["id"] != key or observation["group_id"] != frozen["group_id"] or
                observation["condition_id"] != frozen["condition_id"] or
                observation["source_sha256"] != frozen["source_sha256"] or
                observation["source_locator"] != frozen["source_locator"]):
            raise ValueError("experimental observation does not match frozen source or group")
        label = observation["label"]
        if label not in {"stable", "unstable"}:
            raise ValueError("stability label must be stable or unstable")
        for field in ("score_id", "score_protocol_sha256", "unit"):
            if score[field] != registry[field]:
                raise ValueError("calculated score definition differs from registration")
        if score["group_id"] != frozen["group_id"] or score["condition_id"] != frozen["condition_id"]:
            raise ValueError("score is for the wrong group or experimental conditions")
        evidence = _digest(score["artifact_sha256"], "computed evidence")
        lo, hi = (_number(score[k], k) for k in ("lower", "upper"))
        if lo > hi:
            raise ValueError("reversed computed-score interval")
        interval = [lo, hi] if direction == 1 else [-hi, -lo]
        result[key] = {"label": label, "oriented_interval": interval,
                       "score_artifact_sha256": evidence,
                       "observation_sha256": frozen["observation_sha256"]}
    return result


def _fit(data):
    if set(data) != {"operation", "registration", "registration_sha256", "training"}:
        raise ValueError("threshold calibration only accepts the registered training set")
    registry, ids = _registration(data)
    training = _rows(registry, data["training"], ids["train"])
    frozen = {"registration_sha256": data["registration_sha256"], "training": training,
              "missing": sorted(ids["train"] - training.keys()), "oriented_threshold": None}
    if frozen["missing"] or {x["label"] for x in training.values()} != {"stable", "unstable"}:
        status = "missing"
    else:
        # Every allowed stable score must be above every allowed unstable score.
        lower = max(x["oriented_interval"][1] for x in training.values() if x["label"] == "unstable")
        upper = min(x["oriented_interval"][0] for x in training.values() if x["label"] == "stable")
        frozen["class_gap_edges"] = [lower, upper]
        if lower < upper:
            threshold = lower / 2 + upper / 2
            status = "pass" if lower < threshold < upper else "indeterminate"
            if status == "pass":
                frozen["oriented_threshold"] = threshold
        else:
            # An overlap of uncertainty bounds alone is insufficient to falsify separation.
            max_unstable_lo = max(x["oriented_interval"][0] for x in training.values() if x["label"] == "unstable")
            min_stable_hi = min(x["oriented_interval"][1] for x in training.values() if x["label"] == "stable")
            status = "fail" if max_unstable_lo >= min_stable_hi else "indeterminate"
    frozen["status"] = status
    return {"operation": "fit_anchor_threshold", "status": status, "fit": frozen, "fit_sha256": _hash(frozen)}


def _score(data):
    registry, ids = _registration(data)
    report = data["calibration"]
    fit = report["fit"]
    if (report["operation"] != "fit_anchor_threshold" or _hash(fit) != report["fit_sha256"]
            or fit["registration_sha256"] != data["registration_sha256"]):
        raise ValueError("score requires the matching frozen calibration artifact")
    scored = _rows(registry, data["heldout"], ids["heldout"])
    results = {}
    status = fit["status"]
    if status == "pass":
        if set(fit["training"]) != ids["train"]:
            raise ValueError("frozen training records do not match the split")
        threshold = _number(fit["oriented_threshold"], "frozen threshold")
        for key, record in scored.items():
            lo, hi = record["oriented_interval"]
            if hi < threshold:
                predicted = "unstable"
            elif lo > threshold:
                predicted = "stable"
            else:
                predicted = None
            outcome = "indeterminate" if predicted is None else "pass" if predicted == record["label"] else "fail"
            results[key] = {**record, "prediction": predicted, "status": outcome}
        dispositions = {item["status"] for item in results.values()}
        status = ("fail" if "fail" in dispositions else "missing" if scored.keys() != ids["heldout"]
                  else "indeterminate" if "indeterminate" in dispositions else "pass")
    train_sources = {registry["anchors"][key]["source_sha256"] for key in ids["train"]}
    heldout_sources = {registry["anchors"][key]["source_sha256"] for key in ids["heldout"]}
    return {"operation": "score_anchor_threshold", "status": status, "fit_sha256": report["fit_sha256"],
            "rank_prediction_permitted": status == "pass",
            "heldout": results, "missing": sorted(ids["heldout"] - scored.keys()),
            "heldout_class_coverage": sorted({r["label"] for r in scored.values()}),
            "shared_source_sha256": sorted(train_sources & heldout_sources), "limitations": registry["limitations"]}


def evaluate(data):
    return _fit(data) if data["operation"] == "fit_anchor_threshold" else _score(data)
