"""Own exactly one headless AEDT session per job and always release it."""
from __future__ import annotations

import inspect
import logging
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)


def supported_kwargs(callable_obj: Any, wanted: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    """Keep only keyword arguments that ``callable_obj`` accepts.

    PyAEDT renamed constructor arguments across releases; the HPC version is not
    pinned yet, so unsupported names are dropped with a warning rather than crashing.
    """
    params = inspect.signature(callable_obj).parameters
    if any(p.kind is inspect.Parameter.VAR_KEYWORD for p in params.values()):
        return dict(wanted), []
    kept = {k: v for k, v in wanted.items() if k in params}
    dropped = sorted(set(wanted) - set(kept))
    return kept, dropped


@contextmanager
def aedt_session(project_path: Path | None, design: str | None, version: str = "2025.2",
                 non_graphical: bool = True) -> Iterator[Any]:
    """Start a new headless AEDT desktop, open ``project_path``, and release everything on exit.

    ``project_path=None`` starts with a new empty project (used by the smoke test).
    ``new_desktop=True`` guarantees we never attach to someone else's AEDT process.
    """
    from ansys.aedt.core import Hfss  # lazy: the pure modules must import without PyAEDT

    wanted = dict(
        project=str(project_path) if project_path is not None else None,
        design=design,
        version=version,
        non_graphical=non_graphical,
        new_desktop=True,
        close_on_exit=True,
        remove_lock=True,
    )
    kwargs, dropped = supported_kwargs(Hfss.__init__, wanted)
    if dropped:
        log.warning("This PyAEDT's Hfss() does not accept %s; continuing without them", dropped)
    log.info("Starting AEDT %s non_graphical=%s project=%s design=%s", version, non_graphical, project_path, design)
    hfss = Hfss(**kwargs)
    try:
        yield hfss
    finally:
        try:
            hfss.release_desktop(close_projects=True, close_desktop=True)
            log.info("AEDT session released")
        except Exception:  # noqa: BLE001 - must not mask the original error
            log.exception("release_desktop failed")


def aedt_versions(hfss: Any) -> dict[str, str]:
    import ansys.aedt.core as core

    versions = {"pyaedt": str(getattr(core, "__version__", "unknown"))}
    try:
        versions["aedt"] = str(hfss.odesktop.GetVersion())
    except Exception as exc:  # noqa: BLE001
        versions["aedt"] = f"unknown ({exc})"
    return versions
