# Managed cloud execution

Run commands from the directory containing `chemflow2/`. After completing cloud setup and reviewing the budget policy, use the generated private profile.

```sh
python3 -m chemflow2 batch plan profile.json research/plan.json
python3 -m chemflow2 batch submit profile.json research/plan.json
python3 -m chemflow2 batch status profile.json "$RUN"
python3 -m chemflow2 batch fetch profile.json "$RUN" evidence
```

`submit` prints the run ID before creating resources. Set `RUN` to that ID. Once the cloud workflow execution has started the laptop can disconnect. Each run owns one workflow, one private storage prefix and a sequence of labelled Batch arrays. The CLI freezes only files declared in the protocol. The workflow reads one immutable allocation specification at a time.

The planner combines tasks at the same dependency depth with compatible resource demands and scientific gate prerequisites. Each array uses the smallest supported E2 configuration that meets its tasks' RAM and CPU requirements, with a boot disk sized for scratch plus fixed OS, image and transfer headroom. The current supported CPU configurations range from 2 to 8 vCPUs and offer up to 60 GiB usable memory. A task exceeding measured capabilities, its profile, or one VM's limits prevents submission. Verified private image and licence capabilities belong in the profile.

`max_parallel_vms` caps the slots active per allocation and defaults to one. One slot runs on one VM. By default the planner packs up to eight same-frontier analysis or probe tasks in a slot. It also packs other tasks whose individual timeouts fit within one eighth of the slot deadline. Other tasks receive a dedicated slot. Override with `tasks_per_slot` after calibrating startup, runtime and disk requirements. `max_hours` applies to the entire slot. Allocations run sequentially across the priced waves. On live submission, the client checks remaining regional CPU and SSD disk quotas, plus an all-regions CPU limit when supplied by Compute Engine. Each wave also respects configured per-engine seat caps in `license_parallelism`. The planner searches a fixed E2 menu. Concurrent external usage or provider availability can still delay or reject scheduling. The worker validates its stored budget expiry again immediately before starting each scientific task.

Every finished task uploads a uniquely named archive containing captured logs, receipt, inputs and registered outputs. A create-only commit object publishes its checksum. Later tasks verify their parents' identity, checksum and success before staging declared artifacts. A completed analysis reporting a negative, missing or indeterminate gate publishes its result; the workflow skips its descendants and can continue independent allocations. An engine, integrity or infrastructure failure stops further dispatch. If the service kills an in-flight process before publication, its outputs can be unavailable. Batch/Cloud Logging retains service diagnostics according to the project's retention policy.

`fetch` reconstructs and checks the frozen run directory after the workflow terminates. It also works for a failed run and labels uncommitted tasks as unpublished or blocked. Review `state.json`, individual receipts, cloud job history, scientific reports and original outputs to distinguish research outcomes from execution failures.

```sh
python3 -m chemflow2 batch cancel profile.json "$RUN"
python3 -m chemflow2 batch cleanup profile.json "$RUN"
python3 -m chemflow2 batch cleanup profile.json "$RUN" --delete-evidence
```

Cancellation targets this run's active workflow executions and labelled Batch jobs. Cleanup requires terminal resources. By default it deletes the run workflow and jobs. Fetch first to capture cloud execution history. The explicit deletion flag also removes the private run-prefix objects. Shared images, project-wide billing holds, logging retention, storage soft-delete retention and provider billing delays have separate lifecycles.

The local regression exercises independent workers on the existing probe DAG with real staged inputs, parent artifacts, partial retrieval, corrupt archives, duplicate execution, expired dispatch and gate exits. Production Batch execution, ORCA qualification, resource peaks and billing require a real authorized project.
