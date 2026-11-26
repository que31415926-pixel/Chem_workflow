"""Prepare a private Google Batch / Workflows project for one licensed user."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import tempfile
import uuid

from . import budget
from .protocol import load_plan

ROOT = Path(__file__).resolve().parents[1]
APIS = ("artifactregistry", "batch", "cloudbilling", "cloudbuild", "compute", "iam",
        "logging", "storage", "workflows", "workflowexecutions")
BUILD_MACHINE = "E2_HIGHCPU_8"
BUILD_HOURS = 1
BUILD_DISK_GB = 100
# Locked to the single installer in toolchain.lock.json (ORCA 6.1.0, Open MPI 4.1.8).
INSTALLER_NAME = re.compile(r"orca_6_1_0_linux_x86-64_shared_openmpi418\.run")


def build_allocation():
    return {"id": "private-image-build", "task_ids": ["private-image-build"],
            "max_instances": 1, "max_hours": BUILD_HOURS, "disk_gb": BUILD_DISK_GB,
            "machine_type": BUILD_MACHINE}


def check_installer(installer):
    if (installer.is_symlink() or not installer.is_file() or not installer.stat().st_size
            or not INSTALLER_NAME.fullmatch(installer.name)):
        raise ValueError("supply the locked ORCA 6.1.0 installer orca_6_1_0_linux_x86-64_shared_openmpi418.run (see toolchain.lock.json)")


def settings(project, region, policy_path):
    if not re.fullmatch(r"[a-z][a-z0-9-]{4,28}[a-z0-9]", project):
        raise ValueError("supply a Google Cloud project ID")
    if not re.fullmatch(r"[a-z]+-[a-z]+[1-9][0-9]*", region):
        raise ValueError("supply a Google Cloud region")
    return {"project": project, "region": region, "bucket": f"{project}-chemflow-v3",
            "budget_policy": str(policy_path.resolve()),
            "service_account": f"chemflow-job@{project}.iam.gserviceaccount.com",
            "workflow_service_account": f"chemflow-flow@{project}.iam.gserviceaccount.com"}


def quota_limits(snapshot):
    metrics = {q["metric"]: q for q in snapshot["quotas"]}
    limits = {}
    for name, metric in (("regional_vcpu_limit", "CPUS"), ("regional_disk_limit_gb", "SSD_TOTAL_GB")):
        row = metrics.get(metric)
        if row is None:
            raise ValueError("regional quota unavailable: " + metric)
        remaining = math.floor(float(row["limit"]) - float(row["usage"]))
        if remaining <= 0:
            raise ValueError("no remaining regional quota: " + metric)
        limits[name] = remaining
    return limits


def call(project, *args, output=False):
    cmd = ["gcloud", "--project=" + project, *args]
    result = subprocess.run(cmd, text=True, check=True, stdout=subprocess.PIPE if output else None)
    return result.stdout.strip() if output else None


def bind(project, scope, member, role):
    call(project, *scope, "--member=serviceAccount:" + member, "--role=roles/" + role,
         "--condition=None", "--quiet")


def stage(source, installer, archive):
    """Only executable Python modules and the selected installer can leave the laptop."""
    check_installer(installer)
    code = source / "chemflow2"
    vendored = source / "vendor" / "chemflow-v1" / "chemflow"
    # The repository layout may place the original CLI beside the new engine.
    if not vendored.is_dir():
        vendored = code / "_upstream" / "chemflow"
    if not vendored.is_dir():
        vendored = source / "chemflow"
    if not (vendored / "__main__.py").is_file():
        raise ValueError("upstream chemflow Python package missing")
    dockerfile = source / "deploy" / "byol.Dockerfile"
    if not dockerfile.is_file():
        dockerfile = code / "_deploy" / "byol.Dockerfile"
    members = {dockerfile: "Dockerfile", installer: "orca.run"}
    for directory, prefix in ((code, "chemflow2"), (vendored, "vendor/chemflow-v1/chemflow")):
        for file in directory.glob("*.py"):
            if directory == code and file.stem in {"claims", "v4_rules", "salts", "paper_catalog"}:
                continue
            if file.is_symlink() or not file.is_file():
                raise ValueError(f"unexpected linked source {file}")
            members[file] = prefix + "/" + file.name
    with tarfile.open(archive, "w:gz", compresslevel=1) as bundle:
        for file, name in members.items():
            bundle.add(file, arcname=name, recursive=False)


class ResourceJournal:
    """Track only resources successfully created for this particular installation."""

    def __init__(self, profile, path):
        self.profile = profile
        self.path = path
        self.identity = {k: profile[k] for k in ("project", "region", "bucket")}
        if path.exists():
            self.data = json.loads(path.read_text())
            if self.data["identity"] != self.identity:
                raise ValueError("setup journal belongs to a different location")
        else:
            self.data = {"identity": self.identity, "created": {}}

    def save(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile("w", dir=self.path.parent, delete=False) as temp:
            json.dump(self.data, temp)
        Path(temp.name).replace(self.path)

    def info(self, key):
        p = self.profile
        if key == "bucket":
            args = ["storage", "buckets", "describe", "gs://" + p["bucket"], "--raw"]
        elif key == "repo":
            args = ["artifacts", "repositories", "describe", "chemflow-private", "--location=" + p["region"]]
        else:
            args = ["iam", "service-accounts", "describe", f"chemflow-{key}@{p['project']}.iam.gserviceaccount.com"]
        info = json.loads(call(p["project"], *args, "--format=json", output=True))
        if key == "bucket":
            permissions = info["iamConfiguration"]
            project_number = call(p["project"], "projects", "describe", p["project"],
                                  "--format=value(projectNumber)", output=True)
            if (str(info["projectNumber"]) != project_number or info["location"].lower() != p["region"] or
                    permissions["publicAccessPrevention"] != "enforced" or
                    permissions["uniformBucketLevelAccess"]["enabled"] is not True):
                raise ValueError("private bucket project, region or access controls changed")
            return info["timeCreated"]
        return info["createTime"] if key == "repo" else info["uniqueId"]

    def check(self, key):
        if self.info(key) != self.data["created"][key]:
            raise ValueError("resource was replaced since setup: " + key)

    def create(self, key, *args):
        if key in self.data["created"]:
            self.check(key)
            return
        call(self.profile["project"], *args)
        self.data["created"][key] = self.info(key)
        self.save()

    def remove(self, key, *args):
        if key in self.data["created"]:
            self.check(key)
            call(self.profile["project"], *args, "--quiet")
            del self.data["created"][key]
            self.save()


def teardown(profile, journal, purge):
    project, region = profile["project"], profile["region"]
    if not journal.path.is_file():
        raise ValueError("no resource ownership journal for this profile")
    if call(project, "builds", "list", "--region="+region, "--ongoing", "--format=value(id)", output=True):
        raise ValueError("a regional build is still running")
    # Run-level cleanup removes workflows and Batch job records. Verify this before
    # tearing down identities: an idle workflow may be about to dispatch its next job.
    for service in (["workflows", "list"], ["batch", "jobs", "list"]):
        if call(project, *service, "--location="+region, "--format=value(name)", output=True):
            raise ValueError("clean up regional workflows and Batch jobs before tearing down setup")
    journal.remove("repo", "artifacts", "repositories", "delete", "chemflow-private", "--location="+region)
    for short in ("job", "flow", "build"):
        email = f"chemflow-{short}@{project}.iam.gserviceaccount.com"
        journal.remove(short, "iam", "service-accounts", "delete", email)
    if purge:
        journal.remove("bucket", "storage", "rm", "--recursive", "gs://"+profile["bucket"])
    print(json.dumps({"retained_resources":list(journal.data["created"]), "journal":str(journal.path)}))


def create_resources(profile, journal):
    project, region, bucket = (profile[k] for k in ("project", "region", "bucket"))
    uri = "gs://" + bucket
    call(project, "services", "enable", *(api + ".googleapis.com" for api in APIS))
    bill = json.loads(call(project, "billing", "projects", "describe", project,
                           "--format=json", output=True))
    if bill.get("billingEnabled") is not True:
        raise ValueError("project billing is unavailable; check trial eligibility")
    bindings = json.loads(call(project, "projects", "get-iam-policy", project, "--format=json",
                               output=True))["bindings"]
    if any({"allUsers", "allAuthenticatedUsers"}.intersection(b["members"]) for b in bindings):
        raise ValueError("project has public IAM bindings; use a private dedicated project")
    # All creates are exclusive. An existing similarly named bucket, repository or
    # identity is never silently adopted as a destination for licensed software.
    journal.create("bucket", "storage", "buckets", "create", uri, "--location=" + region,
         "--public-access-prevention", "--uniform-bucket-level-access", "--soft-delete-duration=0")
    journal.create("repo", "artifacts", "repositories", "create", "chemflow-private", "--location=" + region,
         "--repository-format=docker")
    for name in ("job", "flow", "build"):
        journal.create(name, "iam", "service-accounts", "create", "chemflow-" + name)
    with tempfile.TemporaryDirectory() as folder:
        lifecycle = Path(folder) / "rules.json"
        lifecycle.write_text(json.dumps({"rule":[{"action":{"type":"Delete"},
            "condition":{"age":1, "matchesPrefix":["private-build/"]}}]}))
        call(project, "storage", "buckets", "update", uri, "--lifecycle-file="+str(lifecycle))
    repo_iam = json.loads(call(project, "artifacts", "repositories", "get-iam-policy", "chemflow-private",
                               "--location=" + region, "--format=json", output=True))
    if any({"allUsers", "allAuthenticatedUsers"}.intersection(b["members"])
           for b in repo_iam.get("bindings", [])):
        raise ValueError("private image repository has public IAM bindings")
    build_sa = f"chemflow-build@{project}.iam.gserviceaccount.com"
    worker = profile["service_account"]
    flow = profile["workflow_service_account"]
    project_scope = ["projects", "add-iam-policy-binding", project]
    worker_scope = ["iam", "service-accounts", "add-iam-policy-binding", worker]
    repository = ["artifacts", "repositories", "add-iam-policy-binding", "chemflow-private",
                  "--location=" + region]
    for account, role in ((worker, "batch.agentReporter"), (worker, "logging.logWriter"),
                          (flow, "batch.jobsEditor"), (flow, "logging.logWriter"),
                          (build_sa, "logging.logWriter")):
        bind(project, project_scope, account, role)
    bind(project, worker_scope, flow, "iam.serviceAccountUser")
    # Artifact Registry and Storage have no --condition=None flag across all SDK versions.
    for account, role in ((worker, "reader"), (build_sa, "writer")):
        call(project, *repository, "--member=serviceAccount:" + account,
             "--role=roles/artifactregistry." + role)
    for account, role, prefix in ((worker, "objectViewer", "chemflow/"),
                                  (worker, "objectCreator", "chemflow/"),
                                  (build_sa, "objectViewer", "private-build/"),
                                  (flow, "objectViewer", "chemflow/")):
        condition = f"expression=resource.name.startsWith('projects/_/buckets/{bucket}/objects/{prefix}'),title=scoped-objects"
        call(project, "storage", "buckets", "add-iam-policy-binding", uri,
             "--member=serviceAccount:" + account, "--role=roles/storage." + role, "--condition=" + condition)
    return json.loads(call(project, "compute", "regions", "describe", region,
                           "--format=json(quotas,status,zones)", output=True))


def build_image(profile, installer, workspace, root, profile_path):
    project, region, bucket = (profile[k] for k in ("project", "region", "bucket"))
    ident = uuid.uuid4().hex
    archive = workspace / "source.tar.gz"
    stage(root, installer, archive)
    uri = f"gs://{bucket}/private-build/{ident}.tar.gz"
    repository = f"{region}-docker.pkg.dev/{project}/chemflow-private/worker"
    image = f"{repository}:build-{ident}"
    config = {"steps":[{"name":"gcr.io/cloud-builders/docker", "args":["build", "-t", image, "."]}],
              "images":[image], "timeout":str(BUILD_HOURS * 3600) + "s",
              "options":{"logging":"CLOUD_LOGGING_ONLY", "machineType":BUILD_MACHINE,
                         "diskSizeGb":str(BUILD_DISK_GB)}}
    path = workspace / "build.json"
    path.write_text(json.dumps(config))
    budget.reserve_dispatch(profile, profile_path, "build-" + ident, [build_allocation()])
    call(project, "storage", "cp", "--if-generation-match=0", str(archive), uri)
    try:
        sa = f"projects/{project}/serviceAccounts/chemflow-build@{project}.iam.gserviceaccount.com"
        call(project, "builds", "submit", uri, "--region="+region, "--config="+str(path),
             "--service-account="+sa, "--suppress-logs")
    finally:
        call(project, "storage", "rm", uri)
    digest = call(project, "artifacts", "docker", "images", "describe", image,
                  "--format=value(image_summary.digest)", output=True)
    if not re.fullmatch(r"sha256:[a-f0-9]{64}", digest):
        raise ValueError("could not identify immutable image after build")
    return repository + "@" + digest


def main(argv=None):
    cli = argparse.ArgumentParser(prog="chemflow2 cloud-setup", description=__doc__)
    cli.add_argument("project")
    cli.add_argument("--region", default="us-central1")
    cli.add_argument("--profile", type=Path, default=Path("chemflow-cloud.json"))
    cli.add_argument("--budget-policy", type=Path, default=Path("budget-policy.json"))
    cli.add_argument("--installer", type=Path)
    cli.add_argument("--calibration-plan", type=Path)
    cli.add_argument("--apply", action="store_true", help="authorize project setup and a private Cloud Build")
    cli.add_argument("--teardown", action="store_true", help="remove recorded private image repository and identities")
    cli.add_argument("--purge-bucket", action="store_true", help="also delete evidence, budget ledger and bucket")
    cli.add_argument("--byol-confirmed", action="store_true", help="confirm private cloud installation is authorized")
    args = cli.parse_args(argv)
    profile = settings(args.project, args.region, args.budget_policy)
    journal = ResourceJournal(profile, args.profile.with_name(args.profile.name + ".resources.json"))
    if args.purge_bucket and not args.teardown:
        raise ValueError("--purge-bucket requires --teardown")
    if args.teardown:
        if args.apply:
            raise ValueError("choose setup or teardown")
        teardown(profile, journal, args.purge_bucket)
        return 0
    if not args.apply:
        print(json.dumps({"profile":profile, "required_apis":APIS, "build_machine":BUILD_MACHINE,
                          "build_hours":BUILD_HOURS, "build_disk_gb":BUILD_DISK_GB}, indent=2))
        return 0
    if not args.byol_confirmed or not args.installer:
        raise ValueError("live setup requires an authorized private BYOL installer")
    if args.profile.exists():
        raise ValueError("profile must be a new path")
    check_installer(args.installer)
    policy = budget.policy_for(profile, args.profile)
    expected = budget.quote(policy, [build_allocation()])
    budget.admission(policy, {"project":args.project, "currency":"USD", "holds":{}},
                     "setup-preflight", expected)
    tags = []
    if args.calibration_plan:
        plan = load_plan(args.calibration_plan)
        if plan.get("study_scope") != "engineering_calibration" or any(task.get("claims") for task in plan["tasks"]):
            raise ValueError("qualification by calibration requires a claim-free plan")
        tags = sorted({tag for task in plan["tasks"] for tag in task["requires"]})
        profile["calibration_plan_sha256"] = hashlib.sha256(args.calibration_plan.read_bytes()).hexdigest()
    if not shutil.which("gcloud"):
        raise RuntimeError("install Google Cloud CLI or use Google Cloud Shell with these files")
    with tempfile.TemporaryDirectory() as folder:
        area = Path(folder)
        accounts = call(args.project, "auth", "list", "--filter=status:ACTIVE", "--format=value(account)", output=True)
        if not accounts:
            call(args.project, "auth", "login")
        quotas = create_resources(profile, journal)
        profile.update(quota_limits(quotas))
        profile["license_parallelism"] = {"orca": 1}
        profile["image"] = build_image(profile, args.installer, area, ROOT, args.profile)
    profile.update({"node":{"name":"private-batch","host":"controller","cores":8,"memory_mb":61440,
                            "scratch_mb":131072,"engines":["probe","orca","crest","analysis","upstream"],
                            "capabilities":tags}, "max_hours":24, "overhead_gb":64,
                    "max_parallel_vms":1, "regional_quota_snapshot":quotas})
    args.profile.parent.mkdir(parents=True, exist_ok=True)
    with args.profile.open("x") as output:
        json.dump(profile, output, indent=2)
    print(str(args.profile))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
