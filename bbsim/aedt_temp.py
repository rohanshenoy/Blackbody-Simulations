"""AEDT's temp directory for one run (HF-057): the modes ``off``, ``task`` (or a path template) and ``env``, and the
record every run keeps of it.

PyAEDT passes no batch options, so AEDT takes its TempDirectory from per-user settings, which fall back to
``$HOME/Ansoft/.temp`` (the home quota), and the 2025 R2 settings file is shared by every job of the user on a host.

- ``off``: nothing is changed.
- ``task``: inside the session, ``SetTempDirectory(<job>/aedt_temp)`` and read it back; when the session block ends,
  however it ends, set it to the restore value. Never back to what it was before: with several tasks on one host that
  can be another task's job directory. A path template is task mode in the template's directory, with ``{job}``
  replaced by the job directory.
- ``env``: before the session, ``TMPDIR`` and ``TMP`` are exported as ``<job>/aedt_temp`` for AEDT to inherit, and
  nothing in AEDT is set, so nothing is restored and a kill cannot leave a stale setting. AEDT uses them only when no
  TempDirectory is configured, so what it took is read back at session start.

Every mode logs what AEDT reports at session start and records it as ``previous``. A requested directory that AEDT
does not report back, or a temp-directory call that raises, is a warning, or with ``strict`` an ``AedtTempError``.
The record is ``{mode, requested, previous, reported, restored, restore_target}``.
"""
from __future__ import annotations

import logging
import os
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)

MODES = ("off", "task", "env")
JOB_FIELD = "{job}"          # a path template's placeholder for the job directory
JOB_SUBDIR = "aedt_temp"     # the directory task and env mode make inside the job
ENV_VARS = ("TMPDIR", "TMP")


class AedtTempError(RuntimeError):
    """A strict run whose temp directory is not the requested one, or whose temp-directory call raised."""


@dataclass(frozen=True)
class TempSettings:
    """One run's temp-directory settings. ``restore`` is required by ``task``, a template included: the absolute
    TempDirectory AEDT is set back to at the end, the host's TempDirectory before any plan-mode run, read from AEDT's
    .cfg by the cluster probe. ``strict`` needs ``task``. ``template`` (``task`` only) is an absolute path holding
    ``{job}``. A plan may carry a restore value to an ``env`` task, which ignores it."""

    mode: str = "off"
    restore: str | None = None
    strict: bool = False
    template: str | None = None

    def __post_init__(self) -> None:
        if self.mode not in MODES:
            raise ValueError(f"the AEDT temp mode must be one of {', '.join(MODES)}, got {self.mode!r}")
        if self.template is not None:
            if self.mode != "task":
                raise ValueError(f"a path template is task mode, not {self.mode!r}")
            if JOB_FIELD not in self.template:
                raise ValueError(f"the path template {self.template!r} must hold {JOB_FIELD}")
            if not os.path.isabs(self.template.replace(JOB_FIELD, os.sep + "job")):
                raise ValueError(f"the path template {self.template!r} must give an absolute path")
        if self.mode == "task" and not self.restore:
            raise ValueError("--aedt-temp task needs --aedt-temp-restore, the TempDirectory AEDT is set back to at the "
                             "end: the host's TempDirectory before any plan-mode run, read from AEDT's .cfg by the "
                             "cluster probe; a run never restores what it found, which can be another task's job "
                             "directory")
        if self.restore is not None and not os.path.isabs(self.restore):
            raise ValueError(f"--aedt-temp-restore must be an absolute path, got {self.restore!r}")
        if self.strict and self.mode != "task":
            raise ValueError(f"--aedt-temp-strict needs task mode or a path template, not {self.mode!r}")


def parse_mode(text: str) -> tuple[str, str | None]:
    """``--aedt-temp``'s value as (mode, template): one of ``MODES`` with no template, or a path holding ``{job}``,
    which is task mode in that directory."""
    if JOB_FIELD in text:
        return "task", text
    if text in MODES:
        return text, None
    raise ValueError(f"--aedt-temp must be off, task, env or a path template holding {JOB_FIELD}, got {text!r}")


def requested_directory(settings: TempSettings, job_root: Path) -> str | None:
    """The directory a run asks for: ``<job>/aedt_temp`` in task and env mode, the expanded template, or None (off)."""
    if settings.mode == "off":
        return None
    job = os.path.abspath(job_root)
    if settings.template is not None:
        return os.path.normpath(settings.template.replace(JOB_FIELD, job))
    return os.path.join(job, JOB_SUBDIR)


def _record(settings: TempSettings, job_root: Path) -> dict[str, Any]:
    return {"mode": settings.mode, "requested": requested_directory(settings, job_root), "previous": None,
            "reported": None, "restored": None, "restore_target": settings.restore if settings.mode == "task" else None}


def apply_env(settings: TempSettings, job_root: Path) -> dict[str, Any]:
    """Env mode, before the session: make ``<job>/aedt_temp`` and export it as ``TMPDIR`` and ``TMP`` in this process,
    which AEDT inherits when PyAEDT starts it. Other modes change nothing. Returns the run's record (``around_session``
    fills in what AEDT reports)."""
    record = _record(settings, job_root)
    if settings.mode == "env":
        Path(record["requested"]).mkdir(parents=True, exist_ok=True)
        for name in ENV_VARS:
            os.environ[name] = record["requested"]
        log.info("AEDT temp directory, env mode: TMPDIR and TMP exported as %s for AEDT to inherit",
                 record["requested"])
    return record


@contextmanager
def exported_env(settings: TempSettings, job_root: Path) -> Iterator[dict[str, Any]]:
    """``apply_env`` for the length of the block; this process's ``TMPDIR`` and ``TMP`` are put back afterwards (AEDT
    took them when it started), so a second run in the same process starts from the same environment."""
    saved = {name: os.environ.get(name) for name in ENV_VARS}
    try:
        yield apply_env(settings, job_root)
    finally:
        for name, value in saved.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value


def _same_directory(a: str | None, b: str | None) -> bool:
    """Equal after normalization, or two existing names of one directory (``/home`` is ``/resnick/home`` on the
    cluster)."""
    if a is None or b is None:
        return False
    if os.path.normpath(a) == os.path.normpath(b):
        return True
    try:
        return os.path.samefile(a, b)
    except OSError:
        return False


def _problem(settings: TempSettings, message: str, exc: BaseException | None = None, *, quiet: bool = False) -> None:
    """A warning, or with ``strict`` an ``AedtTempError``; ``quiet`` only warns (another exception is on its way)."""
    if settings.strict and not quiet:
        raise AedtTempError(f"{message} (--aedt-temp-strict)") from exc
    log.warning("AEDT temp directory: %s", message)


def _describe(exc: BaseException) -> str:
    return f"{type(exc).__name__}: {exc}"


@contextmanager
def around_session(hfss: Any, settings: TempSettings, job_root: Path) -> Iterator[dict[str, Any]]:
    """Inside the AEDT session: read and log what AEDT reports as its temp directory (``previous``, every mode); in
    task mode point it at the requested directory and read it back, and when the block ends, whatever ends it (an
    exception, ``SystemExit`` from SIGTERM), set it to ``settings.restore``; in env mode read back what AEDT took from
    TMPDIR. Yields the record; ``restore`` may complete it earlier, before the run writes its manifest."""
    record = _record(settings, job_root)
    previous_error = None
    try:
        record["previous"] = str(hfss.odesktop.GetTempDirectory())
    except Exception as exc:  # noqa: BLE001 - reported below as a warning or a strict failure
        previous_error = exc
    try:
        if settings.mode == "task":
            _point_at_job(hfss, settings, record, previous_error)
        else:
            if previous_error is not None:
                _problem(settings, f"GetTempDirectory raised {_describe(previous_error)}", previous_error)
            if settings.mode == "env":
                record["reported"] = record["previous"]
                log.info("AEDT temp directory at session start: %s (env mode exported TMPDIR and TMP = %s)",
                         record["previous"], record["requested"])
                if record["reported"] is not None and not _same_directory(record["requested"], record["reported"]):
                    _problem(settings, f"AEDT reports {record['reported']!r}, not the exported TMPDIR "
                                       f"{record['requested']!r}: a TempDirectory configured for this host wins over "
                                       f"TMPDIR")
            else:
                log.info("AEDT temp directory at session start: %s (mode off: left as it is)", record["previous"])
        yield record
    except BaseException:
        restore(hfss, settings, record, quiet=True)
        raise
    else:
        restore(hfss, settings, record)


def _point_at_job(hfss: Any, settings: TempSettings, record: dict[str, Any],
                  previous_error: BaseException | None) -> None:
    """Task mode: make the requested directory, ``SetTempDirectory`` it and read it back; log previous, requested and
    reported in one line (the manifest is lost when a run fails), then report each problem in turn."""
    requested = record["requested"]
    problems: list[tuple[str, BaseException | None]] = []
    if previous_error is not None:
        problems.append((f"GetTempDirectory raised {_describe(previous_error)}", previous_error))
    try:
        Path(requested).mkdir(parents=True, exist_ok=True)
        hfss.odesktop.SetTempDirectory(requested)
    except Exception as exc:  # noqa: BLE001 - a warning or a strict failure, below
        problems.append((f"SetTempDirectory({requested}) raised {_describe(exc)}", exc))
    try:
        record["reported"] = str(hfss.odesktop.GetTempDirectory())
    except Exception as exc:  # noqa: BLE001
        problems.append((f"GetTempDirectory after SetTempDirectory raised {_describe(exc)}", exc))
    log.info("AEDT temp directory: previous %s, requested %s, reported %s (mode task; set to %s at the end)",
             record["previous"], requested, record["reported"], settings.restore)
    if record["reported"] is not None and not _same_directory(requested, record["reported"]):
        problems.append((f"AEDT reported {record['reported']!r} after SetTempDirectory, not the requested "
                         f"{requested!r}", None))
    for message, exc in problems:
        _problem(settings, message, exc)


def restore(hfss: Any, settings: TempSettings, record: dict[str, Any], *, quiet: bool = False) -> None:
    """Task mode: set AEDT's temp directory to ``settings.restore`` and record what AEDT reports back as ``restored``.
    It acts once: the runner calls it after the run's last AEDT call, before the manifest is written, so the manifest
    records the outcome, and ``around_session`` calls it again on the way out, which acts only if no read-back was
    recorded. A failure is a warning, or with ``strict`` an ``AedtTempError`` unless ``quiet``."""
    if settings.mode != "task" or record["restored"] is not None:
        return
    target = settings.restore
    try:
        hfss.odesktop.SetTempDirectory(target)
    except Exception as exc:  # noqa: BLE001
        _problem(settings, f"restore: SetTempDirectory({target}) raised {_describe(exc)}; AEDT may still use "
                           f"{record['reported']}", exc, quiet=quiet)
        return
    try:
        record["restored"] = str(hfss.odesktop.GetTempDirectory())
    except Exception as exc:  # noqa: BLE001
        _problem(settings, f"restore: GetTempDirectory after SetTempDirectory({target}) raised {_describe(exc)}", exc,
                 quiet=quiet)
        return
    if _same_directory(target, record["restored"]):
        log.info("AEDT temp directory restored to %s", record["restored"])
    else:
        _problem(settings, f"restore: AEDT reported {record['restored']!r} after SetTempDirectory({target!r})",
                 quiet=quiet)
