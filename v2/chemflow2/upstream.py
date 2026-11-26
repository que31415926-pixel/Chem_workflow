"""Run the pinned ChemFlow scientific CLI using a fixed set of data-only operations."""
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import subprocess
import sys


STAGES = {'simulate', 'fit', 'balance', 'equilibrium'}
VENDOR = Path(__file__).resolve().parents[1] / 'vendor' / 'chemflow-v1'
if not (VENDOR / 'chemflow' / '__main__.py').is_file():
    VENDOR = Path(__file__).resolve().parent / '_upstream'
LAUNCH = ('import runpy,site,sys; sys.path.insert(0,site.getusersitepackages()); '
          'sys.path.insert(0,sys.argv.pop(1)); '
          "runpy.run_module('chemflow',run_name='__main__')")
IDENTIFIER = re.compile(r'[A-Za-z_][A-Za-z_0-9]*\Z')


def _regular(folder, filename):
    if not isinstance(filename, str) or not filename or '\\' in filename or '\x00' in filename:
        raise ValueError('expected a relative input filename')
    if any(piece in ('', '.', '..') for piece in filename.split('/')) or Path(filename).is_absolute():
        raise ValueError('expected an input within the task directory')
    path = folder / filename
    if not path.is_file() or path.resolve() != path:
        raise ValueError(f'missing or linked input: {filename}')
    return path


def _fingerprint_files(folder, paths):
    digest = hashlib.sha256()
    for file in sorted(paths):
        digest.update(file.relative_to(folder).as_posix().encode('utf-8') + b'\0')
        digest.update(hashlib.sha256(file.read_bytes()).digest())
    return digest.hexdigest()


def run_upstream(task, folder):
    """Execute one allowlisted CLI stage. All input and output paths are in the attempt."""
    folder = Path(folder).resolve()
    stage = task['stage']
    if stage not in STAGES:
        raise ValueError('unrecognized upstream stage')
    spec = json.loads(_regular(folder, task.get('input')).read_text(encoding='utf-8'))
    if not isinstance(spec, dict):
        raise ValueError('upstream task input must be a JSON object')
    allowed = {'recipe'} if stage == 'equilibrium' else {'project', 'project_source'}
    if stage in ('fit', 'balance'):
        allowed.add('data')
    if stage == 'simulate':
        allowed.add('calibrated')
    if set(spec) - allowed:
        raise ValueError('unsupported upstream options: ' + ', '.join(sorted(set(spec) - allowed)))
    outputs = task.get('outputs')
    if not isinstance(outputs, list) or not any(isinstance(x, str) and x.startswith('products/') for x in outputs):
        raise ValueError('declare scientific outputs in products/')
    for name in outputs:
        if not isinstance(name, str) or (name != 'upstream.out' and not name.startswith('products/')):
            raise ValueError('upstream outputs must be within products/ or upstream.out')
        if any(s in ('', '.', '..') for s in name.split('/')) or '\\' in name or '\x00' in name:
            raise ValueError('invalid output path')
    product = folder / 'products'
    if product.exists() or product.is_symlink() or (folder / 'upstream.out').exists():
        raise ValueError('upstream stage requires fresh output locations')
    if not (VENDOR / 'chemflow' / '__main__.py').is_file():
        raise FileNotFoundError('worker needs vendor/chemflow-v1 from the source distribution')
    args = []
    inputs = []
    if stage == 'equilibrium':
        recipe = _regular(folder, spec.get('recipe'))
        args.extend((stage, str(recipe)))
    else:
        project = spec.get('project')
        if not isinstance(project, str) or not IDENTIFIER.fullmatch(project):
            raise ValueError('project must be a ChemFlow identifier')
        source = spec.get('project_source', 'staged')
        if source not in ('vendored', 'staged'):
            raise ValueError('project_source must be vendored or staged')
        root = VENDOR if source == 'vendored' else folder
        manifest = _regular(root, f'projects/{project}/project.json')
        config = 'material.json' if stage == 'balance' else 'network.json'
        model = _regular(root, f'projects/{project}/{config}')
        inputs.extend([manifest, model])
        if stage in ('simulate', 'fit'):
            network = json.loads(model.read_text(encoding='utf-8'))
            for key in ('conditions', 'hypotheses'):
                labels = network.get(key)
                if not isinstance(labels, dict) or not labels or any(not IDENTIFIER.fullmatch(n) for n in labels):
                    raise ValueError(f'invalid labels in {key}')
        args.extend(('--projects-dir', str(root / 'projects'), stage, project))
        if stage in ('fit', 'balance'):
            args.extend(('--data', str(_regular(folder, spec.get('data')))))
        if 'calibrated' in spec:
            args.extend(('--calibrated', str(_regular(folder, spec['calibrated']))))
    args.extend(('--out', str(product)))
    cores = task.get('resources', {}).get('cores')
    timeout = task.get('timeout_s', 7 * 86400)
    if type(cores) is not int or cores < 1 or type(timeout) is not int or timeout < 1:
        raise ValueError('a positive core reservation and timeout are required')
    environment = os.environ.copy()
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
                'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        environment[key] = str(cores)
    environment.pop('CHEMFLOW_PROJECTS_DIR', None)
    code = [sys.executable, '-I', '-c', LAUNCH, str(VENDOR), *args]
    log_path = folder / 'upstream.out'
    with log_path.open('xb') as log:
        try:
            proc = subprocess.run(code, cwd=folder, env=environment, stdin=subprocess.DEVNULL, stdout=log,
                                  stderr=subprocess.STDOUT, check=False, timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError('ChemFlow upstream stage exceeded its time reservation') from exc
    if proc.returncode:
        with log_path.open('rb') as log:
            log.seek(max(0, log_path.stat().st_size - 8192))
            reason = log.read().decode(errors='replace')
        if ('ModuleNotFoundError' in reason or 'No module named' in reason) and 'Traceback' in reason:
            raise FileNotFoundError('upstream dependency missing; see upstream.out')
        raise RuntimeError(f'upstream {stage} exited with code {proc.returncode}; see upstream.out')
    for name in outputs:
        file = _regular(folder, name)
        if file.stat().st_size == 0:
            raise RuntimeError(f'empty declared output: {name}')
    packages = {}
    for package in ('numpy', 'scipy', 'python-libsbml'):
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            packages[package] = None
    observations = {'stage': stage, 'python': sys.version.split()[0], 'packages': packages,
                    'code_sha256': _fingerprint_files(VENDOR, (VENDOR / 'chemflow').glob('*.py')),
                    'scientific_gate': 'requires_review'}
    if inputs:
        observations['project_inputs_sha256'] = _fingerprint_files(root, inputs)
        observations['project_source'] = source
    return observations
