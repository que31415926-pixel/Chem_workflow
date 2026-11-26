"""Finite-grid inference for reversible solution networks with measured endpoint bounds."""
from itertools import product
import math

from .mechanism import BOLTZMANN_OVER_PLANCK, GAS_CONSTANT_KCAL, _digest, _number, _text


def _positive(value, name):
    value = _number(value, name)
    if value <= 0:
        raise ValueError(f'{name} must be positive')
    return value


def _keys(values, names, field):
    if not isinstance(values, dict) or set(values) - set(names):
        raise ValueError(f'{field} contains unknown species or fields')
    return values


def _grid(data):
    axes = data['parameter_grid']
    if not isinstance(axes, dict):
        raise ValueError('parameter_grid must map names to prespecified values')
    for name, values in axes.items():
        _text(name, 'parameter name')
        if not isinstance(values, list) or not values:
            raise ValueError('each parameter requires a nonempty grid')
        for value in values:
            _number(value, name)
        if any(a >= b for a, b in zip(values, values[1:])):
            raise ValueError('parameter grids must be strictly increasing')
    return axes


def _energy(value, point):
    if isinstance(value, str):
        if value not in point:
            raise ValueError(f'unknown parameter: {value}')
        return point[value]
    return _number(value, 'energy in kcal/mol')


def _system(data, point, condition):
    import numpy as np

    energies = data['species_gibbs_kcal_mol']
    names = list(energies)
    species_index = {name: i for i, name in enumerate(names)}
    g = np.array([_energy(v, point) for v in energies.values()])
    temp = _positive(data['temperature_K'], 'temperature_K')
    initial = _keys(condition['initial_total_M'], names, 'initial_total_M')
    y0 = np.array([_number(initial.get(n, 0), n, lower=0) for n in names], dtype=float)
    partition = np.ones(len(names))
    # Totals are moles per litre of liquid. All mass-action terms use dissolved activity.
    for species, h in _keys(data.get('henry_bar_L_mol', {}), names, 'Henry volatility').items():
        h = _positive(h, 'Henry volatility at simulation temperature')
        headspace = _number(condition['headspace_L'], 'headspace_L', lower=0)
        liquid = _positive(condition['liquid_L'], 'liquid_L')
        partition[species_index[species]] += h * headspace / (0.0831446261815324 * temp * liquid)

    def orders(side):
        if not _keys(side, names, 'reaction side') or any(type(x) is not int or x < 1 for x in side.values()):
            raise ValueError('reaction sides require positive integer stoichiometry')
        idx = np.array([species_index[n] for n in side])
        return idx, np.array(list(side.values()))

    reactions = data['reactions']
    ids = [_text(r['id'], 'reaction id') for r in reactions]
    if not ids or len(set(ids)) != len(ids):
        raise ValueError('reaction ids must be unique and nonempty')
    stoich = np.zeros((len(names), len(ids)))
    thermal, photo, fractions = [], [], {}
    for j, r in enumerate(reactions):
        i, a = orders(r['lhs'])
        k, b = orders(r['rhs'])
        stoich[i, j] -= a
        stoich[k, j] += b
        if r['type'] == 'thermal':
            ts = _energy(r['ts_gibbs_kcal_mol'], point)
            barriers = [ts - float(a @ g[i]), ts - float(b @ g[k])]
            if min(barriers) < 0:
                raise ValueError(f"{r['id']}: negative forward or reverse barrier")
            # Dimensionless 1 M activities, prefactor converted from seconds to hours.
            kf, kr = [3600 * BOLTZMANN_OVER_PLANCK * temp * math.exp(-v / (GAS_CONSTANT_KCAL * temp))
                      for v in barriers]
            thermal.append((j, i, a, k, b, kf, kr))
        elif r['type'] == 'photon':
            if len(i) != 1 or a[0] != 1:
                raise ValueError('photon excitation requires one absorbing molecule')
            q = _number(r['quantum_yield'], 'quantum yield', lower=0)
            fractions[i[0]] = fractions.get(i[0], 0) + q
            photo.append((j, i[0], q))
        else:
            raise ValueError('unknown reaction type')
    if any(q > 1 for q in fractions.values()):
        raise ValueError('photon branches exceed available absorbed photons')
    if photo:
        light = data['light']
        epsilon = _keys(light['epsilon_L_mol_cm'], names, 'extinction coefficients')
        extinction = np.array([_number(epsilon.get(n, 0), n, lower=0) for n in names])
        if any(extinction[i] <= 0 for i in fractions):
            raise ValueError('every photoexcitation requires an absorbing species')
        path = _positive(condition['path_cm'], 'optical path')
        flux = _number(condition['incident_photons_M_h'], 'incident photon flux', lower=0)

    pools = data['conserved']
    if not isinstance(pools, dict) or not pools:
        raise ValueError('at least one material balance required')
    rows = []
    for weights in pools.values():
        _keys(weights, names, 'conservation vector')
        row = [_number(weights.get(n, 0), n, lower=0) for n in names]
        if not any(row):
            raise ValueError('empty conserved pool')
        rows.append(row)
    balance = np.array(rows)
    if not np.allclose(balance @ stoich, 0, atol=1e-12, rtol=0):
        raise ValueError('reaction violates declared pool conservation')
    for pool in condition.get('absent_pools', []):
        if pool not in pools or balance[list(pools).index(pool)] @ y0 != 0:
            raise ValueError(f'nonzero or unknown absent pool: {pool}')

    def derivative(_, y):
        activity = np.maximum(y, 0) / partition
        rates = np.zeros(len(reactions))
        for j, i, a, k, b, kf, kr in thermal:
            rates[j] = kf * np.prod(activity[i] ** a) - kr * np.prod(activity[k] ** b)
        if photo:
            absorption = extinction * activity
            total = np.sum(absorption)
            if total > 0:
                absorbed = -math.expm1(-math.log(10) * path * total) * flux
                for j, i, q in photo:
                    rates[j] = q * absorbed * absorption[i] / total
        return stoich @ rates

    observable = condition['observable']
    idx = species_index[observable['product']]
    weights = balance[list(pools).index(observable['limiting_pool'])]
    denominator = float(weights @ y0)
    if weights[idx] <= 0 or denominator <= 0:
        raise ValueError('yield requires a product and nonzero limiting pool basis')
    return derivative, y0, balance, idx, weights[idx] / denominator


def _endpoint(data, point, condition):
    import numpy as np
    from scipy.integrate import solve_ivp

    derivative, y0, balance, idx, scale = _system(data, point, condition)
    time = _positive(data['end_h'], 'observation time')
    yields = []
    for rtol, atol in [(1e-7, 1e-11), (1e-9, 1e-13)]:
        sol = solve_ivp(derivative, (0, time), y0, method='LSODA', rtol=rtol, atol=atol)
        if not sol.success or not np.isfinite(sol.y).all() or sol.y.min() < -1e-9:
            raise ValueError('integration failed or produced negative concentrations')
        if not np.allclose(balance @ sol.y, (balance @ y0)[:, None], atol=1e-9, rtol=1e-7):
            raise ValueError('material balance drift')
        product_yield = sol.y[idx] * scale
        if product_yield.min() < -1e-8 or product_yield.max() > 1 + 1e-8:
            raise ValueError('yield outside physical range')
        yields.append(float(np.clip(product_yield[-1], 0, 1)))
    if abs(yields[0] - yields[1]) > _positive(data['solver_yield_tolerance'], 'solver tolerance'):
        raise ValueError('endpoint sensitive to solver tolerance')
    return yields[-1]


def _identifiability(axes, training, cache, records, yield_tolerance):
    """Local ranks on the frozen grid; derivatives are with respect to normalized axes."""
    import numpy as np

    varying = [(i, name, grid, {v: j for j, v in enumerate(grid)})
               for i, (name, grid) in enumerate(axes.items()) if len(grid) > 1]
    for record in records:
        values = tuple(record['parameters'][a] for a in axes)
        jacobian = np.zeros((len(training), len(varying)))
        noise = np.zeros_like(jacobian)
        for col, (i, _, grid, index) in enumerate(varying):
            j = index[values[i]]
            low, high = grid[max(0, j - 1)], grid[min(j + 1, len(grid) - 1)]
            left = cache[values[:i] + (low,) + values[i+1:]]
            right = cache[values[:i] + (high,) + values[i+1:]]
            scale = (grid[-1] - grid[0]) / (high - low)
            jacobian[:, col] = [(right[n] - left[n]) * scale for n in training]
            noise[:, col] = 2 * yield_tolerance * scale
        singular = np.linalg.svd(jacobian, compute_uv=False)
        threshold = float(np.linalg.norm(noise))
        rank = int(np.count_nonzero(singular > threshold))
        record['local_identifiability'] = {
            'varying_parameters': [name for _, name, _, _ in varying],
            'normalized_jacobian_singular_values': singular.tolist(),
            'numerical_rank_threshold': threshold,
            'rank': rank, 'unresolved_dimensions': len(varying) - rank}


def evaluate(data):
    """Accept a point only if every training endpoint lies inside its frozen interval."""
    axes = _grid(data)
    _digest(data['model_source_sha256'], 'model source')
    conditions = data['conditions']
    roles = {'training': [], 'holdout': [], 'design': []}
    for name, condition in conditions.items():
        role = condition['role']
        if role not in roles:
            raise ValueError('unknown condition role')
        roles[role].append(name)
        required = {'role', 'initial_total_M', 'observable'}
        extras = {'headspace_L', 'liquid_L', 'absent_pools', 'path_cm', 'incident_photons_M_h'}
        if role == 'training':
            required |= {'yield_interval', 'observation_sha256'}
            _digest(condition['observation_sha256'], 'observation evidence')
            low, high = condition['yield_interval']
            if not 0 <= _number(low, 'lower') <= _number(high, 'upper') <= 1:
                raise ValueError('invalid observed yield interval')
        if set(condition) - (required | extras) or required - set(condition):
            raise ValueError(f'{name}: invalid fields for {role} condition')
    if not roles['training']:
        raise ValueError('at least one training endpoint required')
    referenced = {v for v in data['species_gibbs_kcal_mol'].values() if isinstance(v, str)}
    referenced |= {r['ts_gibbs_kcal_mol'] for r in data['reactions']
                   if r['type'] == 'thermal' and isinstance(r['ts_gibbs_kcal_mol'], str)}
    if referenced != set(axes):
        raise ValueError('every energy parameter requires a used, prespecified grid')
    phase = data.get('phase_contrast')
    if phase:
        if {phase['primary'], phase['secondary']} - set(conditions):
            raise ValueError('phase contrast needs two declared conditions')
        if phase['primary'] == phase['secondary']:
            raise ValueError('phase contrast requires distinct conditions')
        if not 0 < _number(phase['minimum_difference'], 'phase threshold') <= 1:
            raise ValueError('phase threshold must be a positive difference in yield fractions')
    # Derived free-energy coordinates retain the joint binding and substitution dependence.
    coordinates = data.get('energy_coordinates', {})
    for name, weights in coordinates.items():
        _text(name, 'energy coordinate')
        if not _keys(weights, axes, 'energy coordinate terms') or not any(weights.values()):
            raise ValueError('energy coordinate needs nonzero coefficients')
        for coefficient in weights.values():
            _number(coefficient, 'energy coordinate coefficient')
    results, feasible, cache = [], [], {}
    for values in product(*axes.values()):
        point = dict(zip(axes, values))
        predicted = {n: _endpoint(data, point, conditions[n]) for n in roles['training']}
        failures = [n for n, y in predicted.items()
                    if not conditions[n]['yield_interval'][0] <= y <= conditions[n]['yield_interval'][1]]
        record = {'parameters': point, 'incompatible_endpoints': failures}
        if coordinates:
            record['energy_coordinates_kcal_mol'] = {
                n: sum(point[k] * coefficient for k, coefficient in weights.items())
                for n, weights in coordinates.items()}
        cache[values] = predicted
        if not failures:
            predicted.update({n: _endpoint(data, point, conditions[n]) for n in roles['holdout'] + roles['design']})
            record['yield_fractions'] = dict(predicted)
            feasible.append(record)
            if phase:
                delta = predicted[phase['primary']] - predicted[phase['secondary']]
                record['yield_difference'] = delta
                record['secondary_suppressed'] = delta >= phase['minimum_difference']
        results.append(record)
    _identifiability(axes, roles['training'], cache, feasible, data['solver_yield_tolerance'])
    bounds = {n: [min(r['yield_fractions'][n] for r in feasible),
                  max(r['yield_fractions'][n] for r in feasible)] for n in conditions if feasible}
    spreads = {n: bounds[n][1] - bounds[n][0] for n in roles['design'] if feasible}
    support = {a: sorted({r['parameters'][a] for r in feasible}) for a in axes}
    sensitivity = {}
    for i, (axis, grid) in enumerate(axes.items()):
        positions = {v: j for j, v in enumerate(grid)}
        endpoint_slopes = {n: [] for n in roles['training']}
        for values, predicted in cache.items():
            pos = positions[values[i]]
            if pos + 1 >= len(grid):
                continue
            other = values[:i] + (grid[pos+1],) + values[i+1:]
            for n in endpoint_slopes:
                slope = (cache[other][n] - predicted[n]) / (grid[pos+1] - grid[pos])
                endpoint_slopes[n].append(slope)
        sensitivity[axis] = {n: [min(v), max(v)] for n, v in endpoint_slopes.items() if v}
    return {'operation': 'interval_kinetics', 'status': 'compatible' if feasible else 'no_sampled_fit',
            'grid_points': len(results), 'feasible_points': len(feasible), 'points': results,
            'sampled_parameter_support': support, 'yield_envelopes': bounds,
            'training_sensitivity_per_kcal_mol': sensitivity,
            'largest_design_spread': max(spreads, key=spreads.get) if spreads else None,
            'design_spreads': spreads, 'holdouts_excluded_from_selection': roles['holdout']}
