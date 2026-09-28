"""Per-job directory layout so that concurrent runs never share a project, lock, or scratch file."""
from __future__ import annotations

import os
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from bbsim.naming import dataset_dir_name, frequency_label


@dataclass(frozen=True)
class JobDirs:
    root: Path
    project: Path
    scratch: Path
    logs: Path

    def dataset_dir(self, dataset_id: str, frequency_ghz: float, ephi: int) -> Path:
        return self.root / dataset_dir_name(dataset_id, frequency_ghz, ephi)


def default_job_id(env: Mapping[str, str] | None = None, now: datetime | None = None) -> str:
    """``slurm<SLURM_JOB_ID>`` inside Slurm, otherwise a UTC timestamp."""
    env = os.environ if env is None else env
    slurm_id = env.get("SLURM_JOB_ID")
    if slurm_id:
        return f"slurm{slurm_id}"
    now = now or datetime.now(timezone.utc)
    return now.strftime("%Y%m%dT%H%M%SZ")


def create_job_dirs(output_root: Path, dataset_id: str, frequency_ghz: float, job_id: str) -> JobDirs:
    """Create ``<root>/<dataset_id>_<label>/job_<id>/{project,scratch,logs}``; refuse if it exists."""
    root = Path(output_root) / f"{dataset_id}_{frequency_label(frequency_ghz)}" / f"job_{job_id}"
    if root.exists():
        raise FileExistsError(f"job directory already exists: {root}")
    root.mkdir(parents=True, exist_ok=False)
    dirs = JobDirs(root=root, project=root / "project", scratch=root / "scratch", logs=root / "logs")
    for d in (dirs.project, dirs.scratch, dirs.logs):
        d.mkdir()
    return dirs
