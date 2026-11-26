"""Opt-in sampling inside a scientific worker's private attempt directory."""
import os
from pathlib import Path
import stat
from threading import Event, Thread


def scratch_bytes(folder):
    total = 0
    for base, dirs, files in os.walk(folder, followlinks=False):
        for name in files:
            try:
                info = (Path(base) / name).lstat()
            except FileNotFoundError:
                continue
            if stat.S_ISREG(info.st_mode):
                total += info.st_blocks * 512
    return total


def memory_cgroup():
    try:
        for line in Path("/proc/self/cgroup").read_text().splitlines():
            if line.startswith("0::"):
                relative = Path(line[3:].lstrip("/"))
                base = Path("/sys/fs/cgroup")
                paths = [] if ".." in relative.parts else [base / relative]
                for path in [*paths, base]:
                    if (path / "memory.current").is_file():
                        return path
    except OSError:
        pass
    return None


class UsageSampler:
    """Collect interval maxima; preserve nulls when the platform lacks a meter."""

    def __init__(self, folder):
        self.folder = folder
        self.group = memory_cgroup()
        self.max_memory = None
        self.max_scratch = None
        self.samples = 0
        self.incomplete = False
        self.stop = Event()
        self.thread = None

    def sample(self):
        for field, read in (
            ("max_memory", lambda: int((self.group/"memory.current").read_text())),
            ("max_scratch", lambda: scratch_bytes(self.folder)),
        ):
            if field == "max_memory" and self.group is None:
                continue
            try:
                reading = read()
                setattr(self, field, max(getattr(self, field) or 0, reading))
            except (OSError, ValueError):
                self.incomplete = True
        self.samples += 1

    def start(self):
        def track():
            while not self.stop.wait(5):
                self.sample()
        self.sample()
        self.thread = Thread(target=track, name="chemflow-usage", daemon=True)
        self.thread.start()

    def finish(self):
        self.stop.set()
        self.thread.join()
        self.sample()
        peak = None
        if self.group:
            try:
                peak = int((self.group/"memory.peak").read_text())
            except (OSError, ValueError):
                pass
        return {"interval_s":5, "samples":self.samples, "sampling_error":self.incomplete,
                "attempt_scratch_allocated_sampled_max_bytes":self.max_scratch,
                "cgroup_memory_current_sampled_max_bytes":self.max_memory,
                "cgroup_memory_peak_lifetime_bytes":peak,
                "memory_scope":"worker cgroup including sibling processes; lifetime peak can predate attempt",
                "scratch_scope":"allocated regular-file blocks under attempt directory; spikes between samples may be missed"}
