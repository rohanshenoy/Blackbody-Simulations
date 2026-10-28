"""Frequency grids of the job array and the exit step along l (OR-048, OR-037, HF-051, HF-057).

A grid file lists one array task per line, ``<frequency GHz> [<step along l, mm>]``:

- the frequency is a positive integer of at most 6 digits, written with digits only (no sign, leading zero or decimal
  point, so ``500.0`` is refused), strictly ascending and never repeated;
- the optional step is the requested exit-field resolution along l, a finite positive decimal in mm; without it the
  config's own ``exit_field.resolution_mm[1]`` applies;
- ``#`` starts a comment, on a line of its own or after the values; blank lines are skipped;
- lines end at LF, with one CR before it allowed (CRLF), and fields are separated by spaces and tabs only. Any other
  character outside printable ASCII is refused, in comments too: a bare CR, form feed, vertical tab, U+2028 or a
  no-break space would otherwise split a line differently here than in an editor, and turn a frequency into the
  previous entry's step.

Every refusal is a ``GridError`` that names its line; ``read_grid`` adds the file. ``max_entries`` defaults to 1001,
Slurm's default MaxArraySize (indices 0 to 1000); a caller may pass a smaller cap.

The step is the requested resolution, not the exported lattice step: the runner fits the step to the exit face's
extent (HF-045), so crack1's request of 0.0742 mm at 1010 GHz exports 10 mm over 135 steps, 0.07407 mm, 136 x 51
points. At every frequency the step along l may not exceed c/(4f) (OR-048, refining OR-037 and HF-051), 0.07420605 mm
at 1010 GHz; ``effective_l_resolution_mm`` applies that rule to the entry's step or, without one, to the config's.
The config's 0.1 mm meets it up to c/(4 x 0.1 mm) = 749.48 GHz, which is where OR-037's "above about 0.75 THz" came
from; c/(4f) is 0.416, 0.234 and 0.132 mm at 180, 320 and 569 GHz. Above 749.48 GHz a grid entry gives the step, which
the grids write as c/(4f) rounded down to 3 significant figures.

No grid generator: the earlier log grid (180 to 18000 GHz, 24 per decade) is superseded by HF-039(10), and its
generator waits for S62. The first licensed batch (OR-047) uses the hand-written files in ``configs/frequency_grids/``.

Standard library only, so ``SPEED_OF_LIGHT_M_PER_S`` repeats ``bbsim.sampling``'s constant instead of importing that
module, which imports numpy; the two must stay equal.
"""
from __future__ import annotations

import hashlib
import math
import numbers
import re
from dataclasses import dataclass
from pathlib import Path

SPEED_OF_LIGHT_M_PER_S = 299792458.0    # equal to bbsim.sampling.SPEED_OF_LIGHT_M_PER_S
MAX_ENTRIES = 1001                      # Slurm's default MaxArraySize

_FREQUENCY = re.compile(r"[1-9][0-9]{0,5}")   # at most 6 digits, so every frequency converts to a finite float
_NOT_ALLOWED = re.compile(r"[^\t\x20-\x7e]")  # anything but tab and printable ASCII
_FIELD_SEPARATOR = re.compile(r"[ \t]+")
_STEP = re.compile(r"(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")
_FORM = "an entry is '<frequency GHz> [<step along l, mm>]'"


class GridError(ValueError):
    """A grid file, a grid entry or an exit step along l is refused."""


@dataclass(frozen=True)
class GridEntry:
    """One array task: its frequency and the requested exit step along l (None: the config's own step)."""

    frequency_ghz: int
    l_resolution_mm: float | None


def parse_grid(text: str, max_entries: int = MAX_ENTRIES) -> list[GridEntry]:
    """The entries of a grid file's text, in order; ``GridError`` names the first line refused."""
    if isinstance(max_entries, bool) or not isinstance(max_entries, int) or max_entries < 1:
        raise ValueError(f"max_entries must be a positive integer, got {max_entries!r}")
    entries: list[GridEntry] = []
    first_line: dict[int, int] = {}
    for number, raw in enumerate(text.split("\n"), start=1):
        where = f"line {number}: "
        line = raw[:-1] if raw.endswith("\r") else raw
        bad = _NOT_ALLOWED.search(line)
        if bad:
            char = bad.group()
            raise GridError(f"{where}character {char!r} (U+{ord(char):04X}) is not allowed: a grid line holds "
                            "printable ASCII and tabs and ends at LF or CRLF")
        fields = [field for field in _FIELD_SEPARATOR.split(line.split("#", 1)[0]) if field]
        if not fields:
            continue
        if len(entries) == max_entries:
            raise GridError(f"{where}more than {max_entries} entries")
        if len(fields) > 2:
            raise GridError(f"{where}{len(fields)} fields {fields!r}; {_FORM}")
        frequency = _frequency(fields[0], where)
        step = _step(fields[1], where) if len(fields) == 2 else None
        if frequency in first_line:
            raise GridError(f"{where}duplicate frequency {frequency} (first on line {first_line[frequency]})")
        if entries and frequency < entries[-1].frequency_ghz:
            raise GridError(f"{where}frequency {frequency} after {entries[-1].frequency_ghz}: frequencies must be "
                            "strictly ascending")
        first_line[frequency] = number
        entries.append(GridEntry(frequency, step))
    if not entries:
        raise GridError(f"no frequency entries; {_FORM}")
    return entries


def read_grid(path: Path, max_entries: int = MAX_ENTRIES) -> tuple[list[GridEntry], str]:
    """A grid file's entries and the sha256 of the very bytes parsed. ``GridError`` names the file and the line."""
    data = Path(path).read_bytes()
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        line = data[:exc.start].count(b"\n") + 1
        raise GridError(f"{path}: line {line}: not UTF-8 text") from None
    try:
        entries = parse_grid(text, max_entries)
    except GridError as exc:
        raise GridError(f"{path}: {exc}") from None
    return entries, hashlib.sha256(data).hexdigest()


def max_l_step_mm(frequency_ghz: float) -> float:
    """c/(4f) in mm: the largest exit step along l allowed at frequency f (OR-048, OR-037, HF-051)."""
    frequency = _finite_positive(frequency_ghz)
    if frequency is None:
        raise GridError(f"frequency {frequency_ghz!r} GHz must be a finite positive number")
    return SPEED_OF_LIGHT_M_PER_S / (4.0 * frequency * 1e9) * 1e3


def effective_l_resolution_mm(entry: GridEntry, config_l_resolution_mm: float) -> float:
    """The exit step along l a task requests: the entry's own, else the config's ``exit_field.resolution_mm[1]``.

    At every frequency it may not exceed c/(4f) (OR-048, OR-037, HF-051); a step that does, whichever supplied it, is
    a ``GridError`` naming the frequency, the step and the limit.
    """
    if entry.l_resolution_mm is None:
        value = _positive_mm(config_l_resolution_mm, "the config's step along l")
        source = "the config's step along l (no step in the grid entry)"
    else:
        value = _positive_mm(entry.l_resolution_mm, "the grid entry's step along l")
        source = "the grid entry's step along l"
    limit = max_l_step_mm(entry.frequency_ghz)
    if value > limit:
        raise GridError(f"{entry.frequency_ghz} GHz: {source} {value!r} mm exceeds c/(4f) = {limit:.7g} mm, the "
                        f"largest step along l allowed (OR-048, OR-037, HF-051); give the grid entry a step of at most "
                        f"{limit:.7g} mm")
    return value


def _frequency(token: str, where: str) -> int:
    if not _FREQUENCY.fullmatch(token):
        raise GridError(f"{where}frequency {_shown(token)} is not a positive integer GHz of at most 6 digits (digits "
                        "only, no sign, leading zero or decimal point)")
    return int(token)


def _step(token: str, where: str) -> float:
    if not _STEP.fullmatch(token):
        raise GridError(f"{where}step along l {_shown(token)} is not a positive decimal number of mm")
    value = float(token)
    if not math.isfinite(value):
        raise GridError(f"{where}step along l {_shown(token)} is not finite")
    if value <= 0:
        raise GridError(f"{where}step along l {_shown(token)} must be positive")
    return value


def _positive_mm(value: object, what: str) -> float:
    number = _finite_positive(value)
    if number is None:
        raise GridError(f"{what} {value!r} must be a finite positive number of mm")
    return number


def _finite_positive(value: object) -> float | None:
    """``value`` as a float when it is a finite positive real (not a bool), else None; an int too large for a float
    is None, not an OverflowError."""
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        return None
    try:
        number = float(value)
    except (OverflowError, ValueError):
        return None
    return number if math.isfinite(number) and number > 0 else None


def _shown(token: str) -> str:
    """A token for a message: its repr, cut after 20 characters when longer than 24."""
    return repr(token) if len(token) <= 24 else f"{token[:20]!r}... ({len(token)} characters)"
