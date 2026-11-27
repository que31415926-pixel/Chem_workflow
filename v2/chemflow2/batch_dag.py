"""Schedule topological frontiers as resource-matched, capped Google Batch arrays."""
from collections import defaultdict
import hashlib
import json
import math
import re

from .protocol import _order, matches, validate_node

# Usable reservations deliberately leave headroom for the OS and Batch agent.
MACHINES = [
    ("e2-standard-2", 2, 7 * 1024),
    ("e2-standard-4", 4, 14 * 1024),
    ("e2-highmem-4", 4, 28 * 1024),
    ("e2-standard-8", 8, 28 * 1024),
    ("e2-highmem-8", 8, 60 * 1024),
]


MPI_ENVIRONMENT = {"OMPI_MCA_hwloc_base_use_hwthreads_as_cpus": "1",
                   "OMPI_MCA_rmaps_base_oversubscribe": "1"}

def shm_option(memory_mb):
    """Docker gives containers a 64 MB /dev/shm. Multi-process ORCA allocates MPI shared-memory windows there
    (MPI_Win_allocate_shared in LEANSCF needed ~84 MB for a 42-atom r2SCAN-3c CPCM job), so give the container
    half of the task's memory reservation as tmpfs. tmpfs pages count against memory only when written."""
    return f"--shm-size={max(1024, int(memory_mb) // 2)}m"


def encoded(data):
    return json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def sha(data):
    return hashlib.sha256(data).hexdigest()


def positive(value, name, upper):
    if type(value) is not int or not 1 <= value <= upper:
        raise ValueError(f"{name} must be a positive integer <= {upper}")
    return value


def schedule(plan, profile, run, input_hash):
    if not re.fullmatch(r"cf-[0-9a-f]{32}", run) or not re.fullmatch(r"[0-9a-f]{64}", input_hash):
        raise ValueError("invalid run identity or source hash")
    if not re.fullmatch(r"[a-z0-9./_-]+@sha256:[0-9a-f]{64}", profile["image"]):
        raise ValueError("private worker image must have an immutable digest")
    if not re.fullmatch(r"[a-z0-9][a-z0-9.-]{1,61}[a-z0-9]", profile["bucket"]):
        raise ValueError("invalid bucket")
    hours = positive(profile["max_hours"], "max_hours", 168)
    parallel = positive(profile.get("max_parallel_vms", 1), "max_parallel_vms", 1000)
    configured_pack = profile.get("tasks_per_slot")
    if configured_pack is not None:
        positive(configured_pack, "tasks_per_slot", 100)
    overhead = positive(profile["overhead_gb"], "overhead_gb", 65536)
    if overhead < 32:
        raise ValueError("reserve at least 32 GB beyond scratch for OS, image and evidence")
    node = validate_node(dict(profile["node"], placement="controller"))
    tasks = {t["id"]: t for t in plan["tasks"]}
    depths = {}
    gate_ancestors = {}
    groups = defaultdict(list)
    for key in _order(plan["tasks"]):
        task = tasks[key]
        if not matches(task, node):
            raise ValueError(f"{key} exceeds registered software, licence or resource capacity")
        depths[key] = 1 + max((depths[parent] for parent in task["deps"]), default=-1)
        gates = dict.fromkeys(gate for parent in task["deps"] for gate in gate_ancestors[parent])
        gates.update(dict.fromkeys(parent for parent in task["deps"] if tasks[parent]["engine"] == "analysis"))
        gate_ancestors[key] = tuple(sorted(gates, key=lambda t: (depths[t], t)))
        resources = task["resources"]
        choice = next((machine for machine, cores, ram in MACHINES
                       if resources["cores"] <= cores and resources["memory_mb"] <= ram), None)
        if choice is None:
            raise ValueError(f"no supported machine fits {key}")
        disk_gb = math.ceil(resources["scratch_mb"] * 2**20 / 10**9) + overhead
        lightweight = task["engine"] in {"analysis", "probe"} or (
            type(task.get("timeout_s")) is int and 0 < task["timeout_s"] <= hours * 3600 // 8)
        pack = configured_pack or (8 if lightweight else 1)
        groups[(depths[key], choice, disk_gb, gate_ancestors[key], pack)].append(key)

    allocations = []
    for (level, machine_type, disk_gb, gates, pack), ids in sorted(groups.items()):
        # Split frontiers into bounded waves. maxRunDuration covers one packed slot.
        cores = next(cpu for name, cpu, _ in MACHINES if name == machine_type)
        limit = parallel
        if "regional_vcpu_limit" in profile:
            limit = min(limit, profile["regional_vcpu_limit"] // cores)
        if "regional_e2_vcpu_limit" in profile and machine_type.startswith("e2-"):
            limit = min(limit, profile["regional_e2_vcpu_limit"] // cores)
        if "regional_disk_limit_gb" in profile:
            limit = min(limit, profile["regional_disk_limit_gb"] // disk_gb)
        licensed = profile.get("license_parallelism", {})
        for engine in {tasks[t]["engine"] for t in ids}:
            if engine in licensed:
                limit = min(limit, positive(licensed[engine], "license_parallelism", 1000))
        if type(limit) is not int or limit < 1:
            raise ValueError(f"regional quota cannot fit a single {machine_type} VM with {disk_gb} GB for {ids[0]}")
        for offset in range(0, len(ids), limit * pack):
            members = ids[offset:offset + limit * pack]
            slots = [members[i:i + pack] for i in range(0, len(members), pack)]
            allocation_id = f"{run}-a{len(allocations)}"
            schedule_path = f"pieces/{len(allocations)}.json"
            control = {"run": run, "input_hash": input_hash, "node": node, "slots": slots}
            resources = {k: max(tasks[t]["resources"][k] for t in members)
                         for k in ("cores", "memory_mb")}
            commands = ["-m", "chemflow2.batch_worker", profile["bucket"], run,
                        schedule_path, sha(encoded(control))]
            task = {"computeResource": {"cpuMilli": resources["cores"] * 1000,
                                        "memoryMib": resources["memory_mb"]},
                    # E2 vCPUs are hyperthreads: e2-*-4 has 4 vCPUs on 2 cores. Open MPI counts cores
                    # by default, so `mpirun -np 4` (ORCA %pal nprocs 4) refused to start. The scheduler
                    # never requests more cores than the VM has vCPUs, so count hardware threads.
                    "environment": {"variables": MPI_ENVIRONMENT},
                    "maxRunDuration": f"{hours * 3600}s", "maxRetryCount": 0,
                    "runnables": [{"container": {"imageUri": profile["image"],
                                                  "entrypoint": "python3", "commands": commands,
                                                  "options": shm_option(resources["memory_mb"])}}]}
            policy = {"machineType": machine_type, "bootDisk": {"type": "pd-balanced",
                                                                "sizeGb": disk_gb}}
            spec = {"labels": {"chemflow-run": run},
                    "taskGroups": [{"taskSpec": task, "taskCount": len(slots),
                                    "parallelism": len(slots), "taskCountPerNode": 1}],
                    "allocationPolicy": {"instances": [{"policy": policy}],
                                         "serviceAccount": {"email": profile["service_account"]}},
                    "logsPolicy": {"destination": "CLOUD_LOGGING"}}
            allocations.append({"id": allocation_id, "level": level, "machine_type": machine_type,
                                "max_hours": hours, "disk_gb": disk_gb, "max_instances": len(slots),
                                "task_ids": members, "gates": list(gates), "spec": spec, "piece": control})
    return allocations


def workflow(profile, run, count, valid_until):
    """The source is constant-size. Each invocation reads one immutable job spec."""
    # The caller uploads these JSON objects with application/json content-type.
    from urllib.parse import quote
    prefix = quote("chemflow/" + run + "/specs/", safe="")
    storage = f'https://storage.googleapis.com/storage/v1/b/{profile["bucket"]}/o/{prefix}'
    commits = f'https://storage.googleapis.com/storage/v1/b/{profile["bucket"]}/o/' + quote(
        "chemflow/" + run + "/commits/", safe="")
    parent = f'projects/{profile["project"]}/locations/{profile["region"]}'
    connector = "googleapis.batch.v1.projects.locations.jobs."
    steps = [
        {"spec": {"call": "http.get", "args": {
            "url": '${"' + storage + '" + string(i) + ".json"}',
            "query": {"alt": "media"}, "auth": {"type": "OAuth2"}}, "result": "object"}},
        {"enable": {"assign": [{"eligible": True}]}},
        {"guards": {"for": {"value": "gate", "in": "${object.body.gates}", "steps": [
            {"readGate": {"call": "http.get", "args": {
                "url": '${"' + commits + '" + gate + ".json"}', "query": {"alt": "media"},
                "auth": {"type": "OAuth2"}}, "result": "receipt"}},
            {"evaluate": {"switch": [{"condition": '${receipt.body.ok == false}',
                                      "steps": [{"disable": {"assign": [{"eligible": False}],
                                                              "next": "break"}}]}]}}
        ]}}},
        {"skip": {"switch": [{"condition": "${eligible == false}", "next": "continue"}]}},
        {"checkAuthorization": {"switch": [{"condition": '${sys.now() >= time.parse("' + valid_until + '")}',
                                            "steps": [{"expired": {"raise": "Budget authorization expired"}}]}]}},
        # The connector waits for the job and raises OperationError when it ends FAILED (one task failed).
        # Catch only that, so the inspect/check steps below decide: FAILED -> continue, CANCELLED -> stop.
        {"submit": {"try": {"call": connector + "create", "args": {
            "parent": parent, "jobId": '${"' + run + '-a" + string(i)}',
            "body": "${object.body.spec}", "connector_params": {
                "timeout": 7 * 86400 + 3600,
                "polling_policy": {"initial_delay": 5, "multiplier": 1.5, "max_delay": 120}}}},
                    "except": {"as": "e", "steps": [
                        {"failedJob": {"switch": [{"condition": '${"OperationError" in e.tags}', "next": "inspect"}]}},
                        {"rethrow": {"raise": "${e}"}}]}}},
        {"inspect": {"call": connector + "get",
                     "args": {"name": '${"' + parent + '/jobs/' + run + '-a" + string(i)}'},
                     "result": "job"}},
        {"check": {"switch": [{"condition": '${job.status.state == "SUCCEEDED" or job.status.state == "FAILED"}', "next": "continue"},
                              {"condition": '${job.status.state == "DELETION_IN_PROGRESS" or job.status.state == "CANCELLED"}',
                               "steps": [{"stop": {"raise": "Batch job cancelled or deleted; inspect task evidence"}}]}]}},
        {"pause": {"call": "sys.sleep", "args": {"seconds": 60}, "next": "inspect"}},
    ]
    return {"main": {"steps": [{"allocate": {"for": {"value": "i", "range": [0, count - 1],
                                                     "steps": steps}}}]}}
