"""Single-attempt local worker for a trusted research group on a registered host.

The controller provisions an empty private directory and stages its frozen inputs.
Only explicitly named artifacts and the ORCA log are published to the controller.
"""
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone

from .protocol import _relative


NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
SINGLE_POINT = re.compile(r"^\s*FINAL SINGLE POINT ENERGY\s+(" + NUMBER + r")(?:\s|$)")
MODE = re.compile(r"^\s*(\d+):\s*(" + NUMBER + r")\s+cm\*\*-1")
CONTROL_FILES = {"chemflow-task.json", "chemflow-result.json", "__chemflow_transfer__.json",
                 "chemflow-worker.err.log", "chemflow-worker.out.log"}
STAGES = {"opt", "freq", "sp", "ts", "irc"}


class Rejected(Exception):
    def __init__(self, kind, reason):
        super().__init__(reason)
        self.kind = kind


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_name(name):
    try:
        _relative(name, "worker path")
    except ValueError as exc:
        raise Rejected("input", str(exc)) from exc
    if Path(name).parts[0].casefold() in CONTROL_FILES:
        raise Rejected("input", "Reserved task-local filename " + repr(name))
    return name


def task_path(root, name):
    path = root / safe_name(name)
    if path.resolve() != path:
        raise Rejected("input", "Symlink in task-local path " + repr(name))
    return path


def number(token):
    value = float(token.upper().replace("D", "E"))
    if not math.isfinite(value):
        raise Rejected("validation", "Expected a finite numeric result")
    return value


def xyz_atoms(path):
    """Validate an entire single-frame ORCA XYZ geometry and return its atom count."""
    with path.open(encoding="utf-8", errors="replace") as stream:
        try:
            count = int(stream.readline())
        except ValueError:
            count = 0
        if not 1 <= count <= 1_000_000 or not stream.readline():
            raise Rejected("validation", f"Invalid XYZ header in {path.name}")
        for _ in range(count):
            row = stream.readline().split()
            if len(row) != 4 or not re.fullmatch(r"(?:[A-Za-z]{1,3}:?|\d+)", row[0]):
                raise Rejected("validation", f"Invalid XYZ atom in {path.name}")
            if not all(re.fullmatch(NUMBER, value) for value in row[1:]):
                raise Rejected("validation", f"Invalid XYZ coordinates in {path.name}")
            for value in row[1:]:
                number(value)
        if any(line.strip() for line in stream):
            raise Rejected("validation", f"Expected a single endpoint in {path.name}")
    return count


def read_orca_log(path):
    obs = {"normal_termination": False, "optimization_converged": False}
    mode_table = False
    seen_modes = []
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            # Omit echoed input: a comment containing a result marker cannot qualify a job.
            if line.lstrip().startswith(("#", "|")):
                continue
            if "ORCA TERMINATED NORMALLY" in line and line.strip().strip("* ").strip() == "ORCA TERMINATED NORMALLY":
                obs["normal_termination"] = True
            if "THE OPTIMIZATION HAS CONVERGED" in line:
                obs["optimization_converged"] = True
            if "THE OPTIMIZATION DID NOT CONVERGE" in line or "THE OPTIMIZATION HAS NOT CONVERGED" in line:
                obs["optimization_converged"] = False
            match = SINGLE_POINT.match(line)
            if match:
                obs["electronic_energy_eh"] = number(match[1])
            match = re.search(r"^\s*Program Version\s+(\S+)", line)
            if match:
                obs["engine_version"] = match[1]
            match = re.search(r"^\s*Number of atoms\s*[:. ]+\s*(\d+)\s*$", line, re.I)
            if match:
                obs["n_atoms"] = int(match[1])
            match = re.search(r"^\s*Final Gibbs free energy\s*[:. ]+(" + NUMBER + r")\s+Eh", line, re.I)
            if match:
                obs["gibbs_free_energy_eh"] = number(match[1])
            match = re.search(r"^\s*T1 diagnostic\s*[:. ]+(" + NUMBER + r")", line, re.I)
            if match:
                obs["t1_diagnostic"] = number(match[1])
            spin = re.match(r"^\s*Expectation value of <S\*\*2>\s*:\s*(" + NUMBER + r")\s*$", line, re.I)
            if spin:
                obs["spin_squared"] = number(spin[1])
            stability = re.search(r"Stability Analysis indicates an? (UNSTABLE|STABLE) (?:HF/KS )?wave\s*function", line, re.I)
            if stability:
                obs["wavefunction_stable"] = stability[1].lower() == "stable"
            if line.strip() == "IRC PATH SUMMARY":
                obs["irc_path_summary"] = True
            if line.strip() == "VIBRATIONAL FREQUENCIES":
                mode_table = True
                seen_modes = []
            if mode_table and line.strip() == "NORMAL MODES":
                mode_table = False
                obs["frequency_table_closed"] = True
            match = MODE.match(line) if mode_table else None
            if match:
                if int(match[1]) != len(seen_modes):
                    raise Rejected("validation", "Missing or repeated frequency index")
                seen_modes.append(number(match[2]))
    if seen_modes:
        obs["frequencies_cm-1"] = seen_modes
        obs["negative_frequencies_cm-1"] = [v for v in seen_modes if v < 0]
    return obs


def verify_orca(stage, obs, deck, folder, outputs):
    if not obs["normal_termination"]:
        raise Rejected("engine", "ORCA normal termination marker is absent")
    stem = deck.relative_to(folder).with_suffix("").as_posix()
    if stage in {"sp", "opt", "ts"} and "electronic_energy_eh" not in obs:
        raise Rejected("validation", "Final single-point energy is absent")
    if stage in {"opt", "ts"}:
        if not obs["optimization_converged"]:
            raise Rejected("validation", "Optimization convergence marker is absent")
        geometry = stem + ".xyz"
        if geometry not in outputs:
            raise Rejected("validation", "Declare the optimized geometry " + geometry)
        obs["n_atoms_geometry"] = xyz_atoms(folder / geometry)
    if stage == "freq":
        modes = obs.get("frequencies_cm-1", [])
        if not modes or not obs.get("frequency_table_closed") or len(modes) % 3:
            raise Rejected("validation", "Vibrational spectrum is incomplete")
        if obs.get("n_atoms") and len(modes) != 3 * obs["n_atoms"]:
            raise Rejected("validation", "Frequency count differs from 3N")
    if stage == "irc":
        if not obs.get("irc_path_summary"):
            raise Rejected("validation", "IRC path summary is absent")
        ends = [stem + suffix for suffix in ("_IRC_F.xyz", "_IRC_B.xyz")]
        if not all(name in outputs for name in ends):
            raise Rejected("validation", "Declare both IRC endpoint XYZ files")
        if xyz_atoms(folder / ends[0]) != xyz_atoms(folder / ends[1]):
            raise Rejected("validation", "IRC endpoint atom counts differ")
    obs["validation_scope"] = "numerical markers and declared artifacts"
    obs["scientific_gate"] = "not_evaluated"


def validate_deck(deck, stage, resources):
    if deck.suffix.lower() != ".inp" or deck.stat().st_size > 20 * 1024 * 1024:
        raise Rejected("input", "Expected a task-local ORCA .inp of at most 20 MiB")
    # Decks originate from trusted protocol authors; they are not a tenant sandbox.
    contents = "\n".join(line.split("#", 1)[0] for line in deck.read_text().splitlines())
    keywords = {word.lower() for line in contents.splitlines() if line.lstrip().startswith("!")
                for word in line.lstrip()[1:].split()}
    stages = set()
    if keywords & {"opt", "tightopt", "verytightopt", "looseopt", "normalopt"}:
        stages.add("opt")
    if keywords & {"optts", "scants", "neb-ts", "zoom-neb-ts"}:
        stages.add("ts")
    if keywords & {"anfreq", "numfreq", "freq"}:
        stages.add("freq")
    if "irc" in keywords:
        stages.add("irc")
    if not keywords or (stages or {"sp"}) != {stage}:
        raise Rejected("input", "Frozen ORCA deck keywords do not match single declared stage")
    if re.search(r"(?im)^\s*(?:\$new_job\b|%compound\b|%base\b)", contents):
        raise Rejected("input", "Only a single ORCA job with its default basename is supported")
    for match in re.finditer(r"(?im)^\s*\*\s*xyzfile\s+[-+]?\d+\s+\d+\s+(.+?)\s*$", contents):
        xyz_name = safe_name(match[1].strip('"'))
        if not task_path(deck.parent, xyz_name).is_file():
            raise Rejected("input", "Missing staged geometry " + xyz_name)
    pal = [int(v) for v in re.findall(r"\bnprocs\s+(\d+)\b", contents, re.I)]
    pal += [int(k[3:]) for k in keywords if re.fullmatch(r"pal\d+", k)]
    if len(set(pal)) > 1:
        raise Rejected("input", "Conflicting ORCA PAL declarations")
    cores = pal[0] if pal else 1
    if cores < 1 or cores > resources["cores"]:
        raise Rejected("resources", "ORCA processes exceed the task CPU reservation")
    maxcores = [int(v) for v in re.findall(r"(?im)^\s*%maxcore\s+(\d+)\b", contents)]
    if len(maxcores) > 1:
        raise Rejected("input", "Multiple ORCA MaxCore declarations")
    maxcore = maxcores[0] if maxcores else 4096
    if maxcore < 1 or 4 * cores * maxcore > 3 * resources["memory_mb"]:
        raise Rejected("resources", "ORCA MaxCore × processes exceeds 75% of reserved memory_mb")
    return cores


def launch(args, cwd, log, threads, timeout):
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    with log.open("xb") as output:
        proc = subprocess.Popen(args, cwd=cwd, stdout=output, stderr=subprocess.STDOUT,
                                stdin=subprocess.DEVNULL, env=env, start_new_session=True)
        try:
            return proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            proc.wait()
            raise Rejected("timeout", "Attempt reached its wall-time limit")


def configured_binary(engine):
    selected = os.environ.get("CHEMFLOW_" + engine.upper()) or shutil.which(engine)
    if not selected or not Path(selected).is_absolute() or not os.access(selected, os.X_OK):
        raise Rejected("capability", f"A registered executable {engine} installation is required")
    return Path(selected).resolve()


def check_resources(task, folder):
    resources = task.get("resources", {})
    for name in ("cores", "memory_mb", "scratch_mb"):
        value = resources.get(name)
        if type(value) is not int or value < (0 if name == "scratch_mb" else 1):
            raise Rejected("input", "Invalid resource " + name)
    if resources["scratch_mb"] * 1024**2 > shutil.disk_usage(folder).free:
        raise Rejected("resources", "Insufficient free local scratch space")
    timeout = task.get("timeout_s", 7 * 86400)
    if type(timeout) is not int or timeout < 1:
        raise Rejected("input", "timeout_s must be a positive integer")
    return resources, timeout


def published(folder, names):
    artifacts = {}
    for name in names:
        path = task_path(folder, name)
        if not path.is_file() or not path.stat().st_size:
            raise Rejected("validation", "Missing or empty regular output: " + name)
        artifacts[name] = sha256(path)
    return artifacts


def child_usage():
    # Available on the Linux/macOS deployment targets. Python-only analysis CPU
    # is measured separately by process_time(), since it has no child process.
    try:
        import resource
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    except (ImportError, OSError):
        return None
    scale = 1024 if sys.platform.startswith("linux") else (1 if sys.platform == "darwin" else None)
    return (usage.ru_utime, usage.ru_stime, usage.ru_maxrss * scale if scale else None)


def finish_metrics(started, worker_cpu, previous, folder, result, task):
    elapsed = time.monotonic() - started
    metrics = {
        "wall_seconds": round(elapsed, 6),
        "worker_process_cpu_seconds": round(time.process_time() - worker_cpu, 6),
        "requested_core_seconds": None, "child_user_cpu_seconds": None,
        "child_system_cpu_seconds": None, "child_max_rss_bytes": None,
        "scratch_regular_file_bytes_at_finish": None, "scratch_peak_bytes": None,
        "published_artifact_bytes": None, "actual_cost": None,
        "accounting_scope": "Worker-process CPU includes Python work. Child CPU is the before/after "
                            "RUSAGE_CHILDREN difference for reaped processes and can include concurrent "
                            "children if this interpreter is shared. Max RSS is the largest recorded child "
                            "high-water mark, not a sum across MPI ranks or a cgroup peak. An earlier "
                            "larger child prevents attribution. Scratch is the logical size of regular "
                            "files at finish, including staged inputs and surviving intermediates. "
                            "No tariff, storage charge or billing meter is configured."
    }
    cores = task.get("resources", {}).get("cores") if isinstance(task, dict) and isinstance(task.get("resources"), dict) else None
    if type(cores) is int and cores > 0:
        metrics["requested_core_seconds"] = round(cores * elapsed, 6)
    after = child_usage()
    if previous is not None and after is not None:
        metrics["child_user_cpu_seconds"] = round(after[0] - previous[0], 6)
        metrics["child_system_cpu_seconds"] = round(after[1] - previous[1], 6)
        if previous[2] is not None and after[2] is not None and after[2] > previous[2]:
            metrics["child_max_rss_bytes"] = int(after[2])
    try:
        if folder.is_dir() and not folder.is_symlink():
            metrics["scratch_regular_file_bytes_at_finish"] = sum(
                p.stat().st_size for p in folder.rglob("*") if not p.is_symlink() and p.is_file())
            metrics["published_artifact_bytes"] = sum(task_path(folder.resolve(), name).stat().st_size
                                                      for name in result["artifacts"])
    except (OSError, Rejected):
        # The execution result is authoritative even if local accounting is unavailable.
        pass
    return metrics


def run_task(task, folder):
    started = time.monotonic()
    worker_cpu = time.process_time()
    baseline_child_usage = child_usage()
    folder = Path(folder)
    result = {"ok": False, "status": "failed", "error": None, "failure_kind": None,
              "observations": {}, "artifacts": {}, "provenance": {
                  "started_utc": datetime.now(timezone.utc).isoformat(), "host": platform.node(),
                  "platform": platform.system(), "architecture": platform.machine()}}
    provenance = result["provenance"]
    log_name = None
    sampler = None
    try:
        if not folder.is_dir() or folder.is_symlink():
            raise Rejected("input", "Expected a private attempt directory")
        if not isinstance(task, dict) or not isinstance(task.get("id"), str) or not task["id"]:
            raise Rejected("input", "A nonempty task id is required")
        folder = folder.resolve()
        if task.get("engine") in {"orca", "crest"} and ("BATCH_TASK_INDEX" in os.environ or
                                                          os.environ.get("CHEMFLOW_MEASURE") == "1"):
            from .usage import UsageSampler
            sampler = UsageSampler(folder)
            sampler.start()
        outputs = task.get("outputs")
        if not isinstance(outputs, list) or not outputs:
            raise Rejected("input", "Expected a nonempty list of output filenames")
        for name in outputs:
            if task_path(folder, name).exists():
                raise Rejected("input", "Output already exists in attempt directory: " + name)
        if len(set(outputs)) != len(outputs):
            raise Rejected("input", "Duplicate output name")
        files = {}
        for path in folder.rglob("*"):
            if path.is_symlink():
                raise Rejected("input", "Symlink in staged task inputs")
            if path.is_dir():
                continue
            if not path.is_file():
                raise Rejected("input", "Only regular staged input files are permitted")
            name = path.relative_to(folder).as_posix()
            if name not in CONTROL_FILES:
                files[safe_name(name)] = sha256(path)
        if {"orca.out", "crest.out", "upstream.out"} & files.keys():
            raise Rejected("input", "An attempt may not pre-stage a captured engine log")
        if task.get("input") is not None and safe_name(task["input"]) not in files:
            raise Rejected("input", "Primary input is missing")
        engine, stage = task.get("engine"), task.get("stage")
        provenance.update(task_id=task["id"], engine=engine, stage=stage, input_sha256=files)
        if engine in {"command", "probe"} and stage == "probe":
            if any(key in task for key in ("command", "argv", "script", "executable")):
                raise Rejected("input", "Probe does not accept a command")
            payload = json.dumps({"type": "transport_probe", "task": task["id"], "inputs": files},
                                 sort_keys=True) + "\n"
            for name in outputs:
                path = task_path(folder, name)
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("x") as output:
                    output.write(payload)
            result["observations"] = {"type": "transport_probe", "scientific_gate": "not_applicable"}
        elif engine == "orca" and stage in STAGES:
            resources, timeout = check_resources(task, folder)
            deck = folder / safe_name(task.get("input"))
            cores = validate_deck(deck, stage, resources)
            executable = configured_binary("orca")
            provenance.update(executable=str(executable), executable_sha256=sha256(executable),
                              orca_processes=cores, reservations=resources)
            log = folder / "orca.out"
            log_name = "orca.out"
            code = launch([str(executable), deck.name], deck.parent, log, 1, timeout)
            provenance["exit_code"] = code
            result["observations"] = read_orca_log(log)
            if code != 0:
                raise Rejected("engine", f"ORCA exited with status {code}")
            observations_out = task.get("observations_output")
            if observations_out is not None and (observations_out not in outputs or
                    not isinstance(observations_out, str) or not observations_out.endswith(".json")):
                raise Rejected("input", "observations_output must name a declared JSON output")
            for name in set(outputs) - {observations_out}:
                path = task_path(folder, name)
                if not path.is_file() or not path.stat().st_size:
                    raise Rejected("validation", "Missing or empty regular output: " + name)
            verify_orca(stage, result["observations"], deck, folder, outputs)
            if {"stationary_point", "state_diagnostics", "endpoint_signatures"} & task.keys():
                from .pathway import check_stationary_observations
                try:
                    check_stationary_observations(task, result["observations"], folder)
                except (KeyError, IndexError, TypeError, ValueError) as exc:
                    raise Rejected("input", f"invalid ORCA diagnostics contract: {exc}") from exc
            if "engine_version" in task and result["observations"].get("engine_version") != task["engine_version"]:
                raise Rejected("capability", "Reported ORCA version differs from the frozen version")
            if observations_out is not None:
                path = task_path(folder, observations_out)
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("x") as out:
                    geometries = {name: sha256(task_path(folder, name)) for name in outputs
                                  if name.lower().endswith(".xyz")}
                    json.dump({"observations": result["observations"], "log_sha256": sha256(log),
                               "geometry_sha256": geometries,
                               "input_sha256": files, "task_id": task["id"]}, out, allow_nan=False)
        elif engine == "crest" and stage == "conformers":
            from .crest import run_crest
            check_resources(task, folder)
            log_name = "crest.out"
            result["observations"] = run_crest(task, folder, provenance)
        elif engine == "analysis" and stage in {"science", "mechanism", "review"}:
            from .analysis_task import run_analysis
            from .research import run_review
            evaluator = run_review if stage == "review" else run_analysis
            try:
                observations, scientific_outcome = evaluator(task, folder)
            except (ValueError, KeyError, TypeError) as exc:
                raise Rejected("input", str(exc)) from exc
            result["observations"] = observations
            if scientific_outcome:
                # Preserve legitimate negative findings and missing-evidence reports,
                # while blocking downstream tasks requiring this scientific gate.
                result["artifacts"] = published(folder, outputs)
                raise Rejected(scientific_outcome, "Scientific evaluator returned " + scientific_outcome)
        elif engine in {"upstream", "chemflow"} and stage in {"simulate", "fit", "balance", "equilibrium"}:
            from .upstream import run_upstream
            log_name = "upstream.out"
            check_resources(task, folder)
            try:
                result["observations"] = run_upstream(task, folder)
            except FileNotFoundError as exc:
                raise Rejected("capability", str(exc)) from exc
            except ValueError as exc:
                raise Rejected("input", str(exc)) from exc
            except RuntimeError as exc:
                raise Rejected("engine", str(exc)) from exc
        else:
            raise Rejected("capability", "Unsupported engine/stage pair")
        result["artifacts"] = published(folder, outputs)
        result.update(ok=True, status="succeeded")
    except Rejected as exc:
        result.update(error=str(exc), failure_kind=exc.kind)
    except (OSError, ValueError, TypeError, AttributeError) as exc:
        result.update(error=str(exc), failure_kind="worker")
    finally:
        if log_name and (folder / log_name).is_file() and not (folder / log_name).is_symlink():
            result["artifacts"][log_name] = sha256(folder / log_name)
        metrics = finish_metrics(started, worker_cpu, baseline_child_usage, folder, result, task)
        if sampler is not None and sampler.thread is not None:
            usage = sampler.finish()
            metrics["allocation_usage"] = usage
            metrics["scratch_allocated_sampled_max_bytes"] = usage["attempt_scratch_allocated_sampled_max_bytes"]
        result["metrics"] = metrics
        provenance["wall_seconds"] = metrics["wall_seconds"]
        provenance["finished_utc"] = datetime.now(timezone.utc).isoformat()
    return result
