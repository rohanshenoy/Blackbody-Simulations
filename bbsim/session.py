"""Own exactly one headless AEDT session per job, under a launch lease of the licence guard, and always release it."""
from __future__ import annotations

import inspect
import logging
import os
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from bbsim import licence_guard

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
                 non_graphical: bool = True, *, cleanup_wait_s: float = 120.0,
                 scope_observer: Any = None) -> Iterator[Any]:
    """Start a new headless AEDT desktop under a launch lease, open ``project_path``, and release everything on exit.

    ``project_path=None`` starts with a new empty project (used by the smoke test).
    ``new_desktop=True`` asks for a desktop of our own, but PyAEDT 1.7.0 lets PYAEDT_DESKTOP_PORT and
    PYAEDT_PROCESS_ID override it (desktop.py), so the job would run in, and on release close, another job's AEDT.
    Every PYAEDT_* variable is therefore removed from this process first, whoever the caller is.

    Every owned AEDT start takes a launch lease from the shared licence guard (OR-045, HF-057): without
    ``BBSIM_LICENCE_GUARD`` naming the shared guard the session is refused before ``Hfss()`` is called, and without
    ``BBSIM_GUARD_RESERVATION`` it takes a direct one-slot reservation of its own. The launch scope is recorded before
    ``Hfss()``, and the desktop's process and its descendants after it. The outer ``finally`` covers a constructor that
    raises as well as an error in the body: it releases the desktop when one exists, then polls the scope every 2 s for
    up to ``cleanup_wait_s`` and releases the lease only once the scope has stopped; otherwise the lease stays held and
    the scope it saw is logged (HF-058(2)). ``scope_observer`` replaces this host's ``/proc`` observer, for tests.
    """
    stale = sorted(k for k in os.environ if k.startswith("PYAEDT_"))
    if stale:
        log.warning("Removing %s from the environment so PyAEDT starts its own desktop", ", ".join(stale))
        for k in stale:
            del os.environ[k]
    identity = {"tool": "bbsim.session.aedt_session",
                "project": str(project_path) if project_path is not None else None,
                "design": design, "version": version, "non_graphical": bool(non_graphical)}
    with licence_guard.launch_lease_from_env(identity, cleanup_wait_s=cleanup_wait_s, observer=scope_observer) as lease:
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
        lease.guard.begin_launch_scope(lease, identity)
        log.info("Starting AEDT %s non_graphical=%s project=%s design=%s under launch lease %s", version,
                 non_graphical, project_path, design, lease)
        hfss = None
        try:
            try:
                hfss = Hfss(**kwargs)
            except BaseException as exc:
                _note(lease, "constructor_raised", error=f"{type(exc).__name__}: {exc}")
                raise
            _record_processes(lease, hfss)
            yield hfss
        finally:
            if hfss is not None:
                _record_processes(lease, hfss)          # solver processes started since construction
                try:
                    hfss.release_desktop(close_projects=True, close_desktop=True)
                    log.info("AEDT session released")
                except Exception as exc:  # noqa: BLE001 - must not mask the original error
                    log.exception("release_desktop failed")
                    _note(lease, "release_desktop_failed", error=f"{type(exc).__name__}: {exc}")


def _note(lease: licence_guard.LaunchLease, event: str, **details: Any) -> None:
    try:
        lease.guard.note_lease(lease, event, **details)
    except Exception:  # noqa: BLE001 - the lease stays held whatever is noted; never mask the session's error
        log.exception("launch lease %s: could not record %s", lease, event)


def _record_processes(lease: licence_guard.LaunchLease, hfss: Any) -> None:
    """Record the desktop's process (``hfss.desktop_class.aedt_process_id``) and its descendants in the lease."""
    try:
        pid = getattr(getattr(hfss, "desktop_class", None), "aedt_process_id", None)
        if isinstance(pid, bool) or not isinstance(pid, int) or pid <= 0:
            _note(lease, "no_desktop_pid", value=repr(pid))
            return
        observer = lease.guard.observer
        desktop = observer.process(pid)
        found = [{"pid": pid, "start": desktop.get("start") if desktop else None, "role": "desktop",
                  **({"comm": desktop.get("comm")} if desktop else {})}]
        found += [{**child, "role": "descendant"} for child in observer.descendants(pid)]
        lease.guard.set_lease_processes(lease, found)
    except Exception:  # noqa: BLE001 - the launch cgroup still decides cleanup; without it, cleanup stays unknown
        log.exception("launch lease %s: the desktop's processes could not be recorded", lease)


def pyaedt_version() -> str:
    """The installed PyAEDT's version string, exactly as ``aedt_versions`` records it; needs no AEDT session, so a
    run can compare it with a column-order reference before copying anything."""
    import ansys.aedt.core as core

    return str(getattr(core, "__version__", "unknown"))


def aedt_versions(hfss: Any) -> dict[str, str]:
    versions = {"pyaedt": pyaedt_version()}
    try:
        versions["aedt"] = str(hfss.odesktop.GetVersion())
    except Exception as exc:  # noqa: BLE001
        versions["aedt"] = f"unknown ({exc})"
    return versions
