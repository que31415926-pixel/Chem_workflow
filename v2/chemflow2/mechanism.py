"""Evaluate reviewed stationary points, paired pathways and frozen holdouts."""
import hashlib
import json
import math
import re
from datetime import datetime, timezone


GAS_CONSTANT_KCAL = 8.31446261815324 / 4184.0
BOLTZMANN_OVER_PLANCK = 1.380649e-23 / 6.62607015e-34
ENERGY_TO_KCAL = {'kcal/mol': 1.0, 'kJ/mol': 1 / 4.184, 'hartree': 627.5094740631}


def _number(x, label, *, lower=None):
    if type(x) not in (float, int) or not math.isfinite(x) or (lower is not None and x < lower):
        raise ValueError(f'{label}: expected a finite number' + (f' >= {lower}' if lower is not None else ''))
    return x


def _text(x, label):
    if not isinstance(x, str) or not x.strip():
        raise ValueError(f'{label}: non-empty text required')
    return x


def _digest(x, label):
    if not isinstance(x, str) or not re.fullmatch(r'[0-9a-f]{64}', x):
        raise ValueError(f'{label}: expected lowercase SHA-256')
    return x


def _hash(record):
    return hashlib.sha256(json.dumps(record, ensure_ascii=False, sort_keys=True,
                                     separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def _time(value):
    instant = datetime.fromisoformat(_text(value, 'timestamp').replace('Z', '+00:00'))
    if instant.tzinfo is None:
        raise ValueError('timestamps require timezone offsets')
    return instant


def _unique(rows, key):
    if not isinstance(rows, list):
        raise ValueError(f'{key}: list required')
    records = {}
    for row in rows:
        ident = _text(row[key], key)
        if ident in records:
            raise ValueError(f'duplicate {key}: {ident}')
        records[ident] = row
    return records


def _qualify_path(path, spec, condition, cutoff):
    reasons = []
    reviewed = path['review']
    _text(reviewed['reviewer'], 'reviewer')
    for field in ('mode_assignment', 'endpoint_structures', 'electronic_state'):
        finding = reviewed[field]
        _digest(finding['artifact_sha256'], field)
        if finding['accepted'] is not True:
            reasons.append(f'{field} not accepted by reviewer')
    evidence = path['evidence']
    for field in ('frequency', 'irc_forward', 'irc_backward'):
        _digest(evidence[field], field)
    if evidence['irc_forward'] == evidence['irc_backward']:
        reasons.append('IRC evidence must distinguish the two branch records')
    expected = [spec['reactant'], spec['product']]
    if expected[0] == expected[1]:
        raise ValueError('expected a distinct reactant and product for this comparison')
    if sorted(path['irc_endpoints']) != sorted(expected):
        reasons.append('IRC endpoints differ from frozen identities')
    energies = {}
    states = [path['reactant'], path['ts'], path['product']]
    comp = states[0]['composition']
    if not isinstance(comp, dict) or not comp or any(not re.fullmatch('[A-Z][a-z]?', k) or type(v) is not int or v <= 0 for k, v in comp.items()):
        raise ValueError('invalid elemental inventory')
    charge, mult = states[0]['charge'], spec['multiplicity']
    if type(charge) is not int or type(mult) is not int or mult < 1:
        raise ValueError('invalid charge or multiplicity')
    for kind, point in zip(('reactant', 'ts', 'product'), states):
        _digest(point['energy_sha256'], 'energy_sha256')
        if (point['condition_id'] != condition or point['method_id'] != spec['method_id']
                or point['state_id'] != spec['state_id'] or point['multiplicity'] != mult
                or point['charge'] != charge or point['composition'] != comp):
            reasons.append(f'{kind} does not match the condition, method, composition or electronic surface')
        if kind != 'ts' and point['species'] != spec[kind]:
            reasons.append(f'{kind} species differs from frozen identity')
        if not isinstance(point['vibrational_modes_cm1'], list) or not point['vibrational_modes_cm1']:
            raise ValueError(f'{kind}: missing vibrational modes')
        modes = [_number(x, 'vibrational mode') for x in point['vibrational_modes_cm1']]
        geometry_type = point['geometry_type']
        if geometry_type not in ('linear', 'nonlinear'):
            raise ValueError('frozen path only supports polyatomic/diatomic linear or nonlinear structures')
        dof = 3 * sum(comp.values()) - (5 if geometry_type == 'linear' else 6)
        if len(modes) != dof:
            reasons.append(f'{kind} vibrational mode count does not match molecular geometry')
        if sum(x < -cutoff for x in modes) != int(kind == 'ts'):
            reasons.append(f'{kind} has an unexpected number of imaginary modes')
        unit = point['gibbs_unit']
        if unit not in ENERGY_TO_KCAL:
            raise ValueError(f'unsupported Gibbs unit: {unit}')
        energies[kind] = _number(point['gibbs'], 'Gibbs energy') * ENERGY_TO_KCAL[unit]
    if reasons:
        return {'status': 'unqualified', 'reasons': reasons}
    return {'status': 'qualified', 'forward_kcal_mol': energies['ts'] - energies['reactant'],
            'reverse_kcal_mol': energies['ts'] - energies['product'],
            'reaction_kcal_mol': energies['product'] - energies['reactant']}


def _barriers(data):
    protocol = data['protocol']
    commitment = _digest(data['protocol_sha256'], 'protocol_sha256')
    if _hash(protocol) != commitment:
        raise ValueError('protocol differs from its supplied commitment')
    condition = _text(protocol['condition_id'], 'condition_id')
    temp = _number(protocol['temperature_K'], 'temperature_K')
    if temp <= 0:
        raise ValueError('positive absolute temperature required')
    cutoff = _number(protocol['imaginary_tolerance_cm1'], 'imaginary_tolerance_cm1', lower=0)
    specs = _unique(protocol['paths'], 'id')
    pairs = _unique(protocol['comparisons'], 'id')
    if not specs or not pairs:
        raise ValueError('declare the required paths and their directional comparisons')
    supplied = _unique(data['paths'], 'id')
    if set(supplied) - specs.keys():
        raise ValueError('evidence supplied for an undeclared path')
    results = {}
    for ident, spec in specs.items():
        for field in ('method_id', 'state_id', 'reactant', 'product'):
            _text(spec[field], field)
        if ident not in supplied:
            results[ident] = {'status': 'missing'}
            continue
        results[ident] = _qualify_path(supplied[ident], spec, condition, cutoff)
    comparisons = {}
    used = set()
    for ident, pair in pairs.items():
        left, right = pair['left'], pair['right']
        if left not in specs or right not in specs or left == right:
            raise ValueError(f'{ident}: invalid pair of declared paths')
        if specs[left]['method_id'] != specs[right]['method_id']:
            raise ValueError('paired barriers must use the same electronic method')
        expected = pair['expected_sign']
        if type(expected) is not int or expected not in (-1, 1):
            raise ValueError('expected_sign must be +1 or -1 (left minus right)')
        band = _number(pair['indeterminate_band_kcal_mol'], 'indeterminate_band_kcal_mol', lower=0)
        used.update((left, right))
        if any(results[x]['status'] != 'qualified' for x in (left, right)):
            comparisons[ident] = {'status': 'insufficient_evidence'}
            continue
        diff = results[left]['forward_kcal_mol'] - results[right]['forward_kcal_mol']
        verdict = 'within_band' if abs(diff) <= band else 'agreement' if diff * expected > 0 else 'disagreement'
        comparisons[ident] = {'status': verdict, 'left_minus_right_kcal_mol': diff}
    if used != set(specs):
        raise ValueError('each declared path must participate in a comparison')
    return {'operation': 'paired_barriers', 'protocol_sha256': commitment, 'condition_id': condition,
            'paths': results, 'comparisons': comparisons, 'all_directions_agree':
            all(row['status'] == 'agreement' for row in comparisons.values())}


def _validate_prediction(spec):
    for field in ('protocol_sha256', 'model_sha256'):
        _digest(spec[field], field)
    for field in ('author', 'observable', 'unit', 'condition_id', 'experimental_procedure', 'falsification_rule'):
        _text(spec[field], field)
    train, targets = _unique(spec['training'], 'sample_id'), _unique(spec['predictions'], 'sample_id')
    if not targets:
        raise ValueError('at least one held-out prediction required')
    design = spec.get('holdout_design', 'prospective')
    if design not in ('prospective', 'archived'):
        raise ValueError('holdout_design must be prospective or archived')
    if design == 'archived':
        if _text(spec['custodian'], 'custodian') == spec['author']:
            raise ValueError('archived holdout requires a custodian separate from prediction author')
        for target in targets.values():
            _text(target['source_locator'], 'archived source locator')
            _digest(target['archived_source_sha256'], 'archived_source_sha256')
    tg = {_text(row['group_id'], 'group_id') for row in targets.values()}
    for row in train.values():
        _digest(row['observation_sha256'], 'observation_sha256')
        if _text(row['group_id'], 'group_id') in tg:
            raise ValueError('overlapping training and held-out sample groups')
    if train.keys() & targets.keys():
        raise ValueError('overlapping training and held-out sample IDs')
    for row in targets.values():
        if _number(row['lower'], 'lower prediction bound') > _number(row['upper'], 'upper prediction bound'):
            raise ValueError('reversed prediction interval')
    return train, targets


def _freeze(data):
    spec = data['specification']
    if set(data) != {'operation', 'specification'}:
        raise ValueError('freeze takes only a prospective specification')
    _validate_prediction(spec)
    frozen = {'specification': spec, 'frozen_at': datetime.now(timezone.utc).isoformat(), 'schema_version': 1}
    return {'operation': 'freeze_prediction', 'frozen': frozen, 'sha256': _hash(frozen)}


def _score_holdout(data):
    frozen = data['frozen']
    expected = _digest(data['sha256'], 'sha256')
    if _hash(frozen) != expected or frozen['schema_version'] != 1:
        raise ValueError('frozen prediction changed or unsupported schema')
    spec = frozen['specification']
    train, targets = _validate_prediction(spec)
    frozen_at = _time(frozen['frozen_at'])
    observed = _unique(data['observations'], 'sample_id')
    if observed.keys() - targets.keys():
        raise ValueError('undeclared holdout observation')
    used_sources = {r['observation_sha256'] for r in train.values()}
    results = {}
    for ident, target in targets.items():
        if ident not in observed:
            results[ident] = {'status': 'awaiting_observation'}
            continue
        obs = observed[ident]
        digest = _digest(obs['observation_sha256'], 'observation_sha256')
        if digest in used_sources:
            raise ValueError('holdout observation artifact matches training evidence')
        obs_time = _time(obs['observed_at'])
        now = datetime.now(timezone.utc)
        if spec.get('holdout_design', 'prospective') == 'archived':
            unmasked = _time(obs['unmasked_at'])
            if not obs_time <= frozen_at < unmasked <= now:
                raise ValueError('archived observation must predate freeze and be released after freeze')
            if (obs['released_by'] != spec['custodian'] or obs['source_locator'] != target['source_locator']
                    or digest != target['archived_source_sha256']):
                raise ValueError('archived observation differs from frozen source or custodian')
            _digest(obs['unmask_receipt_sha256'], 'unmask_receipt_sha256')
        elif not frozen_at < obs_time <= now:
            raise ValueError('prospective observation must follow freeze and precede scoring time')
        if (obs['group_id'] != target['group_id'] or obs['unit'] != spec['unit']
                or obs['observable'] != spec['observable'] or obs['condition_id'] != spec['condition_id']):
            raise ValueError('observation does not match the frozen group, observable, unit or condition')
        if obs['source_type'] != 'experiment':
            raise ValueError('expected experimental data for this holdout')
        _text(obs['observer'], 'observer')
        value = _number(obs['value'], 'observed value')
        uncertainty = _number(obs['absolute_uncertainty'], 'absolute_uncertainty', lower=0)
        lo, hi = value - uncertainty, value + uncertainty
        if lo > target['upper'] or hi < target['lower']:
            status = 'outside'
        elif lo >= target['lower'] and hi <= target['upper']:
            status = 'inside'
        else:
            status = 'overlap'
        results[ident] = {'status': status, 'observation_sha256': digest, 'observed_interval': [lo, hi]}
    return {'operation': 'score_holdout', 'sha256': expected, 'targets': results,
            'all_observed': not (targets.keys() - observed.keys())}


def _rate(data):
    if data['model'] != 'unimolecular_TST' or data['gibbs_unit'] not in ENERGY_TO_KCAL:
        raise ValueError('requires unimolecular_TST and a supported molar Gibbs unit')
    _digest(data['barrier_evidence_sha256'], 'barrier_evidence_sha256')
    _text(data['standard_state'], 'standard_state')
    temp = _number(data['temperature_K'], 'temperature_K')
    kappa = _number(data['transmission_coefficient'], 'transmission_coefficient')
    if temp <= 0 or kappa <= 0:
        raise ValueError('positive temperature and transmission coefficient required')
    barrier = _number(data['activation_gibbs'], 'activation_gibbs') * ENERGY_TO_KCAL[data['gibbs_unit']]
    log_k = math.log(kappa) + math.log(BOLTZMANN_OVER_PLANCK) + math.log(temp) - barrier / (GAS_CONSTANT_KCAL * temp)
    try:
        k = math.exp(log_k)
    except OverflowError:
        k = None
    return {'operation': 'eyring', 'log_k_per_s': log_k, 'k_per_s': k, 'model': data['model'],
            'standard_state': data['standard_state'], 'barrier_evidence_sha256': data['barrier_evidence_sha256']}


def evaluate(data):
    from . import anchor_threshold, isotopes
    from . import interval_kinetics, solution
    from .pathway import match_endpoints

    functions = {'fit_anchor_threshold': anchor_threshold.evaluate, 'score_anchor_threshold': anchor_threshold.evaluate,
                 'salt_quantity': solution.salt_quantity, 'isotope_observables': isotopes.evaluate, 'solution_cycle': solution.evaluate, 'irc_endpoints': match_endpoints,
                 'interval_kinetics': interval_kinetics.evaluate, 'paired_barriers': _barriers, 'freeze_prediction': _freeze,
                 'score_holdout': _score_holdout, 'eyring': _rate}
    if not isinstance(data, dict) or data.get('operation') not in functions:
        raise ValueError('mechanism operation must be one of ' + ', '.join(functions))
    try:
        return functions[data['operation']](data)
    except KeyError as exc:
        raise ValueError(f'missing field: {exc.args[0]}') from exc
