"""Provenance recording: file hashes, Slurm and git context, JSON writer."""
from __future__ import annotations

import hashlib
import json
import os
import platform
import socket
import subprocess
import sys
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]


def sha256_file(path: Path, chunk_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def slurm_environment(env: Mapping[str, str] | None = None) -> dict[str, str]:
    env = os.environ if env is None else env
    return {k: v for k, v in sorted(env.items()) if k.startswith("SLURM_")}


def git_commit(repo_root: Path = REPO_ROOT) -> str | None:
    try:
        out = subprocess.run(["git", "-C", str(repo_root), "rev-parse", "HEAD"],
                             capture_output=True, text=True, check=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    sha = out.stdout.strip()
    return sha or None


def build_manifest(**sections: Any) -> dict[str, Any]:
    """Standard provenance header plus caller-supplied sections."""
    manifest: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "hostname": socket.gethostname(),
        "slurm": slurm_environment(),
        "git_commit": git_commit(),
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
