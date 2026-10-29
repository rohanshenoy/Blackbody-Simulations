"""Provenance recording: file hashes, Slurm context, the code checkout's state, JSON writer."""
from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import socket
import subprocess
import sys
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
# Seconds each git call may take on the checkout (NFS on the cluster); a slower one stops the run (HF-057)
GIT_TIMEOUT_S = 30


def sha256_file(path: Path, chunk_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def slurm_environment(env: Mapping[str, str] | None = None) -> dict[str, str]:
    env = os.environ if env is None else env
    return {k: v for k, v in sorted(env.items()) if k.startswith("SLURM_")}


class EnvironmentStateError(RuntimeError):
    """The code checkout's state could not be read: git could not run, failed or timed out."""


@dataclass(frozen=True)
class EnvironmentState:
    """The code checkout a run or a plan uses (HF-057): its ``HEAD`` commit, whether tracked files differ from it, and
    the installed PyAEDT's version."""

    commit: str
    dirty: bool
    pyaedt_version: str

    def record(self) -> dict[str, Any]:
        return asdict(self)


def _git(repo: Path, *args: str) -> str:
    command = ["git", "-C", str(repo), "--no-optional-locks", *args]
    try:
        run = subprocess.run(command, capture_output=True, text=True, check=True, timeout=GIT_TIMEOUT_S)
    except subprocess.TimeoutExpired as exc:
        raise EnvironmentStateError(f"{' '.join(command)} timed out after {GIT_TIMEOUT_S} s") from exc
    except subprocess.CalledProcessError as exc:
        raise EnvironmentStateError(f"{' '.join(command)} failed with exit {exc.returncode}: "
                                    f"{(exc.stderr or '').strip()}") from exc
    except OSError as exc:
        raise EnvironmentStateError(f"git could not run on {repo}: {exc}") from exc
    return run.stdout


def environment_state(repo: Path = REPO_ROOT) -> EnvironmentState:
    """The state of the checkout ``repo``: ``HEAD``, and dirty when ``git --no-optional-locks status --porcelain
    --untracked-files=no`` lists anything (untracked files never count); the PyAEDT version is
    ``importlib.metadata.version("pyaedt")``. One seam for the planner, the task pre-flight and the manifest (spec
    section 4). Nothing is guessed: a git error or timeout raises ``EnvironmentStateError``, and a missing PyAEDT
    ``importlib.metadata.PackageNotFoundError``."""
    commit = _git(repo, "rev-parse", "--verify", "HEAD").strip()
    if not re.fullmatch(r"[0-9a-f]{40}|[0-9a-f]{64}", commit):
        raise EnvironmentStateError(f"git rev-parse HEAD in {repo} printed {commit!r}, not a commit")
    dirty = bool(_git(repo, "status", "--porcelain", "--untracked-files=no").strip())
    return EnvironmentState(commit=commit, dirty=dirty, pyaedt_version=metadata.version("pyaedt"))


def build_manifest(*, environment: EnvironmentState, **sections: Any) -> dict[str, Any]:
    """Standard provenance header plus caller-supplied sections. ``environment`` is the checkout's state the run read
    once before it started (``environment_state``); the header's ``git_commit`` and ``environment_state`` come from
    it, so no second git call is made at the end of a run and no commit is guessed (HF-057)."""
    manifest: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "hostname": socket.gethostname(),
        "slurm": slurm_environment(),
        "git_commit": environment.commit,
        "environment_state": environment.record(),
    }
    manifest.update(sections)
    return manifest


def write_json(path: Path, data: Any) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True, default=_json_default) + "\n")


def _json_default(obj: Any) -> Any:
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (set, tuple)):
        return list(obj)
    if hasattr(obj, "tolist"):  # numpy scalars and arrays
        return obj.tolist()
    return str(obj)
