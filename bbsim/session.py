"""Own exactly one headless AEDT session per job and always release it."""
from __future__ import annotations

import inspect
import logging
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)

# Without these the session could open a blank project or attach to another user's desktop on the node.
ESSENTIAL_KWARGS = ("project", "version", "non_graphical", "new_desktop")


def require_essential_kwargs(callable_obj: Any, extra: tuple[str, ...] = ()) -> None:
    """Raise if the installed Hfss() does not accept the arguments the headless session depends on."""
    params = inspect.signature(callable_obj).parameters
    if any(p.kind is inspect.Parameter.VAR_KEYWORD for p in params.values()):
        return
    missing = [k for k in (*ESSENTIAL_KWARGS, *extra) if k not in params]
    if missing:
        raise RuntimeError(
            f"installed PyAEDT Hfss() does not accept {missing}; this runner needs a PyAEDT release with the "
            f"project/design/version/new_desktop arguments (PyAEDT >= 0.10). Signature: {inspect.signature(callable_obj)}"
        )


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
    from ansys.aedt.core import Hfss, settings  # lazy: the pure modules must import without PyAEDT

    # On an exception inside one of its own methods PyAEDT releases every desktop it knows (HPC step 4a lost
    # its session to a failed mesh-statistics export). This session owns the desktop; it is released here.
    if hasattr(settings, "release_on_exception"):
        settings.release_on_exception = False
    wanted = dict(
        project=str(project_path) if project_path is not None else None,
        design=design,
        version=version,
        non_graphical=non_graphical,
        new_desktop=True,
        close_on_exit=True,
        remove_lock=True,
    )
    require_essential_kwargs(Hfss.__init__, extra=("design",) if design is not None else ())
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
