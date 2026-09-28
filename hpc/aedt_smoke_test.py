#!/usr/bin/env python3
"""Start headless AEDT 2025 R2, create a scratch HFSS project, report versions, close.

Run on HPC inside a Slurm allocation after `source /home/rshenoy/BBRSim/bb_env.sh`:
    python hpc/aedt_smoke_test.py --workdir /home/rshenoy/BBRSim/scratch/smoke
This verifies desktop startup and the base license; it does not solve anything.
"""
from __future__ import annotations

import argparse
import inspect
import logging
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bbsim.session import aedt_session, aedt_versions  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", default="2025.2")
    parser.add_argument("--workdir", type=Path, default=None)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    import ansys.aedt.core as core
    print("pyaedt version:", core.__version__)
    print("Hfss.__init__ signature:", inspect.signature(core.Hfss.__init__))

    workdir = args.workdir or Path(tempfile.mkdtemp(prefix="aedt_smoke_"))
    workdir.mkdir(parents=True, exist_ok=True)
    project = workdir / "smoke.aedt"
    if project.exists():
        project.unlink()

    t0 = time.perf_counter()
    with aedt_session(None, design="smoke", version=args.version) as hfss:
        print(f"desktop started in {time.perf_counter() - t0:.1f}s")
        print("versions:", aedt_versions(hfss))
        box = hfss.modeler.create_box([0, 0, 0], [1, 1, 1], name="probe", material="vacuum")
        print("created object:", box.name, "model units:", hfss.modeler.model_units)
        hfss.save_project(str(project))
        print("saved:", project, project.exists())
    print("closed cleanly")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
