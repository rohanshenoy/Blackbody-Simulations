"""``--verify TREE``: a written tree re-checked against BBRsim's loader rules with the standard library only, so it runs
anywhere, Fir's bare ``python3`` included (HF-057, G4-011, G4-012, S65).

It reads the tree alone: no plan, no job directory, and not the contents of ``tree_manifest.json``. Every rule below is
one BBRsim's reader applies when it loads the tree (its code in parentheses), or the stricter form the assembler's trees
take; a tree that passes is one BBRsim's loader accepts as far as these rules reach:

- the tree holds ``waveguides/``, ``SHA256SUMS`` and ``tree_manifest.json``, and nothing else; every entry under the
  tree is a copy, never a symlink, and nothing is a fifo or another special file;
- every entry under ``waveguides/`` is a directory ``<id>_<token>GHz_Ephi=0`` or ``=1`` or a sidecar file
  ``<id>_<token>GHz.dataset.json``: the token is the text after the last ``_`` and before ``GHz``, which must parse in
  full as C ``strtod`` parses it and be above 0 and finite, and the ID is the rest (``BBRCrackLibrary::Discover``). A
  dataset directory without the frequency (``<id>_Ephi=0``) is refused (BBR011); any other entry fails, a name BBRsim
  skips (a token that does not parse, NaN, 0 or below, a frequency in MHz or THz, ``<id>_Ephi=1``) as not a dataset of
  this tree, an infinite token as one BBRsim lists and then refuses at its sidecar (BBR024), and a stray file or
  directory;
- every stem is a whole dataset: its sidecar beside both Ephi directories, so an orphan or a partial trio fails
  (BBR024 without the sidecar; BBR001 and BBR002 without a directory's CSVs), and no two stems of one ID name the same
  frequency ("500" and "5e2"; BBR011);
- each Ephi directory holds ``far_field.csv`` and ``waveguide.csv``, optionally ``manifest.json``, and nothing else; an
  empty one fails (BBR001, BBR002);
- each sidecar is JSON as BBRsim's parser takes it (no NaN or Infinity, no number beyond a double's range), schema 2.x
  (one schema per tree, and the trees this assembler writes are schema 2, S65), whose ``dataset_id`` and
  ``frequency_label`` are its name's and whose ``frequency_ghz`` is within 0.1 % of the directory frequency
  (``BBRDatasetSidecar::Parse``, BBR024); its ``far_field.columns`` and ``exit_field.columns`` are the schema-2 lists
  (absent BBR024, otherwise BBR025); its ``files`` has exactly the four CSV keys, each recording the file's sha256,
  bytes and rows (the non-blank lines after the header), as BBRsim's validator checks them (BBR024);
- per dataset ID, every sidecar's ``tree_invariant.invariant`` agrees with the lowest frequency's on the blocks
  BBRsim's reader compares (``InvariantDiff``, BBR024; G4-012), through the same projection the assembly uses;
- each CSV's header is the schema-2 column list of its file (one leading UTF-8 byte-order mark and one trailing CR
  dropped; ``BBRHFSSData`` reads the columns by position, BBR013), and the ``Freq`` of its first data row (the first
  non-empty line after the header) is within 0.1 % of the directory frequency (``CheckFrequencyColumn``, BBR009);
- every ``SHA256SUMS`` line is ``<sha256>  <path>`` and verifies, names a file of the tree once, and every file under
  ``waveguides/`` and ``tree_manifest.json`` is listed.

BBRsim's loader checks every row's ``Freq`` and much more of each sidecar and table (frames, modes, the per-key rows);
its validator, with ``sha256sum -c SHA256SUMS``, is the check of record where this one stops. ``verify`` prints one
``FAIL:`` line per problem (an entry refused once is not reported again by the rules that would read it), then
``VERIFY PASS: <n> datasets -> <TREE>; <id>: <f>, ... GHz x <keys> keys`` per ID (exit 0) or ``VERIFY FAIL: <n>
problem(s) in <TREE>: <the first>`` (exit 1); a name that is not UTF-8 is printed backslash-escaped. It writes nothing
into the tree; a report goes only to a new file outside it. A path that is not a readable directory, or a report path
that exists, has no directory or lies inside the tree (by device and inode, so a case variant, a symlink or a bind
mount is caught), is refused (exit 2). It imports the standard library and ``bbsim.tree_invariant`` (standard library)
only.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import stat
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from bbsim.tree_invariant import InvariantError, invariant, invariant_diff, schema_major

WAVEGUIDES, TREE_MANIFEST, SUMS = "waveguides", "tree_manifest.json", "SHA256SUMS"
SIDECAR_SUFFIX = ".dataset.json"
CSV_NAMES = ("far_field.csv", "waveguide.csv")
EPHI_FILES = (*CSV_NAMES, "manifest.json")
# The schema-2 columns, by position, as BBRsim's sidecar parser requires them and BBRHFSSData reads them
COLUMNS = {"far_field.csv": ("Freq", "Ephi", "IWavePhi", "IWaveTheta", "Phi", "Theta", "rEphi_real", "rEphi_imag",
                             "rEtheta_real", "rEtheta_imag"),
           "waveguide.csv": ("Freq", "Ephi", "IWavePhi", "IWaveTheta", "OutgoingPower", "IngoingPower", "X", "Y", "Z",
                             "Ex_real", "Ex_imag", "Ey_real", "Ey_imag", "Ez_real", "Ez_imag", "InterceptFactor")}
FREQ_REL = 1e-3                    # BBR009 and the sidecar's frequency_ghz: 0.1 % of the directory frequency
CHUNK = 1 << 20
HEAD_LIMIT = 1 << 20               # the header and the first data row are read from a CSV's first MiB
REPORT_SCHEMA = "bbsim-verify-report/1"
BOM = b"\xef\xbb\xbf"
_BLANKS = b" \t\x0b\x0c"           # with \r and \n, what bytes.strip() strips
_BLANK_START = re.compile(rb"\n[\n \t\x0b\x0c]")      # a line that starts empty or with a blank
_SEPARATORS = re.compile(rb"[\r\n]")
_EPHI_DIR = re.compile(r"(?P<stem>.*)_Ephi=(?P<ephi>[01])")
_SUM_LINE = re.compile(r"(?P<digest>[0-9a-f]{64})  (?P<path>.+)")
# What C strtod parses in full (C locale): leading whitespace, a sign, then a decimal or hexadecimal number, inf,
# infinity or nan in any case, and nothing after it (Python's float() also takes trailing blanks and underscores)
_STRTOD = re.compile(r"[ \t\n\v\f\r]*[+-]?(?:"
                     r"(?P<hex>0[xX](?:[0-9a-fA-F]+\.?[0-9a-fA-F]*|\.[0-9a-fA-F]+)(?:[pP][+-]?[0-9]+)?)"
                     r"|(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
                     r"|(?i:inf|infinity|nan))")
_UNITS = {"mhz": 1e-3, "ghz": 1.0, "thz": 1e3}
_OTHER_UNIT = re.compile(r"[+-]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?[A-Za-z]+")    # 500MHz, 0.5THz
SKIPPED = ("a name BBRsim skips (it reads <id>_<freq>GHz_Ephi=0 or =1 and <id>_<freq>GHz.dataset.json), not a "
           "dataset of this tree")


def strtod_full(token: str) -> float | None:
    """``float(token)`` when C ``strtod`` would parse all of it, else None."""
    match = _STRTOD.fullmatch(token)
    if match is None:
        return None
    text = token.lstrip(" \t\n\v\f\r")
    return float.fromhex(text) if match["hex"] else float(text)


def parse_frequency_ghz(token: str) -> float | None:
    """A ``Freq`` cell as BBRsim reads it, ``<number><MHz|GHz|THz>`` (whitespace anywhere ignored, the unit in any
    case, the number parsed in full and above 0), in GHz; None when it does not parse."""
    text = re.sub(r"[ \t\n\v\f\r]", "", token)
    mult = _UNITS.get(text[-3:].lower()) if len(text) >= 4 else None
    value = strtod_full(text[:-3]) if mult is not None else None
    return value * mult if value is not None and value > 0 else None


def parse_stem(stem: str, *, ephi0: bool = False) -> tuple[str, str, float] | str:
    """``<id>_<token>GHz`` as (dataset ID, token, frequency in GHz), or what is wrong with it. ``ephi0``: the stem is
    an ``_Ephi=0`` directory's, the one name without a frequency BBRsim refuses (``<id>_Ephi=0``, BBR011); any other
    name without one it skips."""
    if not stem.endswith("GHz"):
        last = stem.rpartition("_")[2]
        if _OTHER_UNIT.fullmatch(last):
            return (f"a frequency in a unit other than GHz ({last!r}): BBRsim reads only <id>_<freq>GHz and skips this "
                    "name, not a dataset of this tree")
        if ephi0:
            return ("a dataset directory without a frequency: BBRsim refuses <id>_Ephi=0 for a placed crack <id> "
                    "(BBR011); a dataset's name is <id>_<freq>GHz_Ephi=N")
        return SKIPPED
    dataset_id, sep, token = stem[:-3].rpartition("_")
    if not sep or not dataset_id:
        return "no dataset ID before the frequency (<id>_<freq>GHz)"
    # BBRCrackLibrary::Discover takes a token strtod parses in full to a value above 0: it skips any other (NaN, 0
    # and below), and lists an infinite one, whose sidecar BBRDatasetSidecar::Load then refuses (BBR024)
    value = strtod_full(token)
    if value is None:
        return (f"the frequency {token!r} does not parse as BBRsim parses one (a number in full): BBRsim skips such a "
                "directory, not a dataset of this tree")
    if not value > 0:
        return (f"the frequency {token!r} is not above 0: BBRsim skips such a directory, not a dataset of this tree")
    if not math.isfinite(value):
        return (f"the frequency {token!r} is infinite: BBRsim lists such a directory as a frequency and then refuses "
                "its sidecar (BBR024)")
    return dataset_id, token, value


def _content(piece: bytes) -> bool:
    """Whether a piece of one line (no CR or LF in it) holds a byte other than blanks."""
    return bool(piece.translate(None, _BLANKS))


def _fields(line: bytes, *, header: bool) -> list[str]:
    """A line's comma-separated fields as BBRsim splits them: one trailing CR dropped, no empty last field (and for the
    header one leading UTF-8 byte-order mark dropped first)."""
    if header and line.startswith(BOM):
        line = line[len(BOM):]
    if line.endswith(b"\r"):
        line = line[:-1]
    parts = line.decode("utf-8", "replace").split(",")
    if parts[-1] == "":
        parts.pop()
    return parts


@dataclass
class Facts:
    """What one read of a file gives: its sha256 and size, and for a CSV its rows (the non-blank lines after the
    header, as BBRsim's validator counts them: ``bytes.splitlines`` and ``strip``), its header's fields and its first
    data row (the first non-empty line after the header, one trailing CR dropped; None when none is found within the
    first MiB, ``head_complete`` telling whether that was the whole file)."""
    sha256: str
    size: int
    rows: int | None = None
    header: list[str] | None = None
    first_row: bytes | None = None
    head_complete: bool = True


def scan_file(path: Path, *, csv: bool, chunk: int = CHUNK) -> Facts:
    """One streamed read of ``path`` (``chunk`` bytes at a time): its ``Facts``. Raises ``OSError``."""
    hashed, size, lines, open_line, head = hashlib.sha256(), 0, 0, False, bytearray()
    with open(path, "rb") as fh:
        while data := fh.read(chunk):
            hashed.update(data)
            size += len(data)
            if not csv:
                continue
            if len(head) < HEAD_LIMIT:
                head += data[:HEAD_LIMIT - len(head)]
            if b"\r" in data or _BLANK_START.search(data):
                parts = _SEPARATORS.split(data)                 # at least two: the chunk holds a separator
                lines += (open_line or _content(parts[0])) + sum(map(_content, parts[1:-1]))
                open_line = _content(parts[-1])
                continue
            # Only LF separates lines here, and every line that starts in this chunk after one starts with a byte
            # that is not blank, so each is a non-blank line: count them, and judge the two ends
            first, last = data.find(b"\n"), data.rfind(b"\n")
            if first < 0:
                open_line = open_line or _content(data)
                continue
            lines += (open_line or _content(data[:first])) + data.count(b"\n", first + 1, last + 1)
            open_line = _content(data[last + 1:])
    if not csv:
        return Facts(hashed.hexdigest(), size)
    lines += open_line
    complete = size <= HEAD_LIMIT
    head_lines = bytes(head).split(b"\n")
    if not complete or head_lines[-1] == b"":      # a cut line may go on; nothing after a last LF is no line
        head_lines.pop()
    header = _fields(head_lines[0], header=True) if head_lines else None
    first = next((r for r in (r[:-1] if r.endswith(b"\r") else r for r in head_lines[1:]) if r), None)
    return Facts(hashed.hexdigest(), size, max(lines - 1, 0), header, first, complete)


def _kind(mode: int) -> str:
    return "a symlink" if stat.S_ISLNK(mode) else "a stray directory" if stat.S_ISDIR(mode) else \
        "a stray file" if stat.S_ISREG(mode) else "neither a regular file nor a directory"


def _one_line(text: object) -> str:
    """One printed line: line breaks (a path or an error may hold one) become spaces, other spacing is kept, and what
    UTF-8 cannot encode (a name that is not UTF-8 arrives from os.listdir as surrogates) is backslash-escaped, so the
    line prints under a strict UTF-8 stdout."""
    return re.sub(r"[\r\n]+", " ", str(text)).encode("utf-8", "backslashreplace").decode("utf-8")


def _say(text: object) -> None:
    print(_one_line(text))


def _number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _no_constant(name: str) -> None:
    raise ValueError(f"{name} is not JSON")


def _int(text: str) -> int:
    """A JSON integer as BBRsim's parser takes it: one beyond a double's range is refused (its parser fails on the
    overflow), so no later rule meets an integer that a float conversion cannot hold."""
    value = int(text)
    try:
        float(value)
    except OverflowError:
        raise ValueError(f"the integer {text[:13]}... ({len(text.lstrip('-'))} digits) overflows a double") from None
    return value


def _float(text: str) -> float:
    value = float(text)
    if math.isinf(value):
        raise ValueError(f"the number {text} overflows a double")
    return value


@dataclass
class _Stem:
    name: str
    dataset_id: str
    token: str
    frequency: float
    dirs: set[int] = field(default_factory=set)
    sidecar: bool = False
    refused: set[str] = field(default_factory=set)      # its parts present but refused (a symlink): not "missing"
    sidecar_data: dict | None = None
    projection: dict | None = None


class _Verify:
    """One verification of ``tree`` (the module docstring's rules), collecting every problem."""

    def __init__(self, tree: Path):
        self.tree, self.wg = tree, tree / WAVEGUIDES
        self.problems: list[str] = []
        self.stems: dict[str, _Stem] = {}
        self.files: dict[str, Path] = {}        # the tree's files, by their path in SHA256SUMS
        self.facts: dict[str, Facts] = {}
        self.refused: set[str] = set()          # paths reported once, which no later rule reports again
        self.listed = 0

    def fail(self, problem: str) -> None:
        self.problems.append(_one_line(problem))

    def _lstat(self, path: Path) -> int | None:
        try:
            return os.lstat(path).st_mode
        except OSError as exc:
            self._refuse(path, f"cannot be read ({exc})")
            return None

    def _list(self, directory: Path) -> list[str] | None:
        try:
            return sorted(os.listdir(directory))
        except OSError as exc:
            self._refuse(directory, f"cannot be listed ({exc})")
            return None

    def _rel(self, path: Path) -> str:
        return path.relative_to(self.tree).as_posix()

    def _refuse(self, path: Path, problem: str) -> None:
        self.fail(f"{path}: {problem}")
        self.refused.add(self._rel(path))

    def run(self) -> None:
        found = set()
        for name in self._list(self.tree) or []:
            path, mode = self.tree / name, self._lstat(self.tree / name)
            if mode is None:
                continue
            expected = {WAVEGUIDES: stat.S_ISDIR, SUMS: stat.S_ISREG, TREE_MANIFEST: stat.S_ISREG}.get(name)
            if stat.S_ISLNK(mode):
                self._refuse(path, "a symlink; a tree holds copies, never links")
            elif expected is None:
                self._refuse(path, f"{_kind(mode)}; a tree holds only {WAVEGUIDES}/, {SUMS} and {TREE_MANIFEST}")
            elif not expected(mode):
                self._refuse(path, "not a directory" if name == WAVEGUIDES else "not a regular file")
            else:
                found.add(name)
        if WAVEGUIDES not in found and WAVEGUIDES not in self.refused:
            self.fail(f"{self.wg}: no {WAVEGUIDES}/ directory")
        for name in (SUMS, TREE_MANIFEST):
            if name not in found and name not in self.refused:
                self.fail(f"{self.tree}: no {name}")
        if TREE_MANIFEST in found:
            self.files[TREE_MANIFEST] = self.tree / TREE_MANIFEST
        if WAVEGUIDES in found:
            self._waveguides()
        self._hash()
        for stem in sorted(self.stems.values(), key=lambda s: s.name):
            self._sidecar(stem)
            self._csvs(stem)
        self._agreement()
        if SUMS in found:
            self._sums()

    # --- the entries ----------------------------------------------------------------------------------------------

    def _stem(self, name: str, parsed: tuple[str, str, float]) -> _Stem:
        return self.stems.setdefault(name, _Stem(name, *parsed))

    def _waveguides(self) -> None:
        names = self._list(self.wg)
        if names == []:
            self.fail(f"{self.wg}: holds no dataset")
        for name in names or []:
            path, mode = self.wg / name, self._lstat(self.wg / name)
            if mode is None:
                continue
            if stat.S_ISLNK(mode) or not (stat.S_ISDIR(mode) or stat.S_ISREG(mode)):
                self._refuse(path, "a symlink; a tree holds copies, never links" if stat.S_ISLNK(mode) else
                             "neither a regular file nor a directory")
                self._note_refused_part(name)
                continue
            if stat.S_ISDIR(mode):
                self._dataset_dir(name, path)
            else:
                self._sidecar_file(name, path)

    def _note_refused_part(self, name: str) -> None:
        """A refused entry that names a dataset's part: the part is there (and reported), not missing."""
        match = _EPHI_DIR.fullmatch(name)
        stem = match["stem"] if match else name[:-len(SIDECAR_SUFFIX)] if name.endswith(SIDECAR_SUFFIX) else None
        parsed = parse_stem(stem) if stem is not None else "no"
        if not isinstance(parsed, str):
            self._stem(stem, parsed).refused.add(name)

    def _dataset_dir(self, name: str, path: Path) -> None:
        match = _EPHI_DIR.fullmatch(name)
        if match is None:
            self._refuse(path, "a stray directory, not <id>_<freq>GHz_Ephi=0 or =1")
            return
        parsed = parse_stem(match["stem"], ephi0=match["ephi"] == "0")
        if isinstance(parsed, str):
            self._refuse(path, parsed)
            return
        self._stem(match["stem"], parsed).dirs.add(int(match["ephi"]))
        entries = self._list(path)
        if entries is None:
            return
        if not entries:
            self.fail(f"{path}: an empty directory; it needs far_field.csv and waveguide.csv (BBR001, BBR002)")
            return
        for entry in entries:
            file, mode = path / entry, self._lstat(path / entry)
            if mode is None:
                continue
            if stat.S_ISLNK(mode):
                self._refuse(file, "a symlink; a tree holds copies, never links")
            elif entry not in EPHI_FILES or not stat.S_ISREG(mode):
                self._refuse(file, f"{_kind(mode)}; an Ephi directory holds far_field.csv, waveguide.csv and "
                                   "manifest.json only")
            else:
                self.files[self._rel(file)] = file
        for csv in CSV_NAMES:
            if csv not in entries:
                self.fail(f"{path}: no {csv} ({'BBR001' if csv == 'far_field.csv' else 'BBR002'})")

    def _sidecar_file(self, name: str, path: Path) -> None:
        if not name.endswith(SIDECAR_SUFFIX):
            self._refuse(path, "a stray file, not <id>_<freq>GHz.dataset.json")
            return
        parsed = parse_stem(name[:-len(SIDECAR_SUFFIX)])
        if isinstance(parsed, str):
            self._refuse(path, parsed)
            return
        self._stem(name[:-len(SIDECAR_SUFFIX)], parsed).sidecar = True
        self.files[self._rel(path)] = path

    def _hash(self) -> None:
        """Every file of the tree read once: hashed, and a CSV's rows, header and first data row taken from those
        bytes; a sidecar is parsed from the bytes hashed."""
        for rel, path in sorted(self.files.items()):
            try:
                if rel.endswith(SIDECAR_SUFFIX):
                    data = path.read_bytes()
                    self.facts[rel] = Facts(hashlib.sha256(data).hexdigest(), len(data))
                    self._parse_sidecar(rel, path, data)
                else:
                    self.facts[rel] = scan_file(path, csv=rel.endswith(".csv"))
            except OSError as exc:
                self._refuse(path, f"cannot be read ({exc})")

    def _parse_sidecar(self, rel: str, path: Path, data: bytes) -> None:
        stem = self.stems[rel[len(WAVEGUIDES) + 1:-len(SIDECAR_SUFFIX)]]
        try:
            sidecar = json.loads(data.decode("utf-8"), parse_constant=_no_constant, parse_int=_int,
                                 parse_float=_float)
        except (UnicodeDecodeError, ValueError) as exc:
            self.fail(f"{path}: not valid JSON for BBRsim's parser ({exc}) (BBR024)")
            return
        except RecursionError:
            self.fail(f"{path}: not valid JSON for BBRsim's parser (nested too deep to parse) (BBR024)")
            return
        if not isinstance(sidecar, dict):
            self.fail(f"{path}: the top level must be a JSON object (BBR024)")
            return
        stem.sidecar_data = sidecar

    # --- the datasets ---------------------------------------------------------------------------------------------

    def _sidecar(self, stem: _Stem) -> None:
        """The stem's trio and the identity, files and projection of its sidecar."""
        missing = [f"{stem.name}_Ephi={e}" for e in (0, 1) if e not in stem.dirs
                   and f"{stem.name}_Ephi={e}" not in stem.refused]
        if not stem.sidecar and f"{stem.name}{SIDECAR_SUFFIX}" not in stem.refused:
            missing.append(f"{stem.name}{SIDECAR_SUFFIX}")
        if missing:
            codes = (["BBR001", "BBR002"] if any("_Ephi=" in m for m in missing) else []) + (
                ["BBR024"] if missing[-1].endswith(SIDECAR_SUFFIX) else [])
            self.fail(f"{self.wg / stem.name}: {'an orphan' if len(missing) == 2 else 'an incomplete dataset'}: "
                      f"{', '.join(missing)} missing; a dataset is {stem.name}{SIDECAR_SUFFIX} beside "
                      f"{stem.name}_Ephi=0/ and {stem.name}_Ephi=1/ ({', '.join(codes)})")
        sidecar, path = stem.sidecar_data, self.wg / f"{stem.name}{SIDECAR_SUFFIX}"
        if sidecar is None:
            return
        try:
            major = schema_major(sidecar)
        except InvariantError as exc:
            self.fail(f"{path}: {exc} (BBR024)")
            stem.sidecar_data = None
            return
        if major != 2:
            self.fail(f"{path}: schema_version {sidecar['schema_version']!r}; a tree of this assembler is schema 2.x "
                      "(S65), and a schema-1 tree is BBRsim's validator's to check")
            stem.sidecar_data = None
            return
        if sidecar.get("dataset_id") != stem.dataset_id:
            self.fail(f"{path}: dataset_id {sidecar.get('dataset_id')!r}, but the name says {stem.dataset_id!r} "
                      "(BBR024)")
        label = f"{stem.token}GHz"
        if sidecar.get("frequency_label") != label:
            self.fail(f"{path}: frequency_label {sidecar.get('frequency_label')!r}, but the name says {label!r} "
                      "(BBR024)")
        frequency = sidecar.get("frequency_ghz")
        if not (_number(frequency) and math.isfinite(frequency)):
            self.fail(f"{path}: frequency_ghz must be a finite number, got {frequency!r} (BBR024)")
        elif abs(frequency - stem.frequency) > FREQ_REL * stem.frequency:
            self.fail(f"{path}: frequency_ghz {frequency!r} is not within 0.1 % of the directory's "
                      f"{stem.frequency:g} GHz (BBR024)")
        self._files_block(stem, path, sidecar.get("files"))
        self._columns(path, sidecar)
        try:
            stem.projection = invariant(sidecar)
        except InvariantError as exc:
            self.fail(f"{path}: {exc} (BBR024)")
        except RecursionError:
            self.fail(f"{path}: a block nested too deep to compare (BBR024)")

    def _columns(self, path: Path, sidecar: dict) -> None:
        """The declared column lists, as BBRsim's sidecar parser requires them: present (BBR024), an array (BBR024), and
        exactly the schema-2 list of their file, strings in the order BBRHFSSData reads them (BBR025). A missing
        ``exit_field`` is the projection's to report."""
        for block, csv in (("far_field", "far_field.csv"), ("exit_field", "waveguide.csv")):
            parent, want = sidecar.get(block), list(COLUMNS[csv])
            if not isinstance(parent, dict):
                if block == "far_field":
                    self.fail(f"{path}: far_field is missing or not an object (BBR024)")
                continue
            if "columns" not in parent:
                self.fail(f"{path}: {block}.columns is missing (BBR024)")
                continue
            columns = parent["columns"]
            if not isinstance(columns, list):
                self.fail(f"{path}: {block}.columns must be an array of column names, got {columns!r} (BBR024)")
            elif not all(isinstance(c, str) for c in columns) or columns != want:
                self.fail(f"{path}: {block}.columns {columns!r} differ from the positional order BBRHFSSData reads, "
                          f"{want!r} (BBR025)")

    def _files_block(self, stem: _Stem, path: Path, files: object) -> None:
        want = {f"{stem.name}_Ephi={e}/{csv}" for e in (0, 1) for csv in CSV_NAMES}
        if not isinstance(files, dict):
            self.fail(f"{path}: files is missing or not an object (BBR024)")
            return
        if set(files) != want:
            self.fail(f"{path}: files has the keys {sorted(files)}, not exactly the four {sorted(want)} (BBR024)")
        for rel in sorted(want & set(files)):
            entry = files[rel]
            if not (isinstance(entry, dict) and isinstance(entry.get("sha256"), str)
                    and all(type(entry.get(k)) is int for k in ("bytes", "rows"))):
                self.fail(f"{path}: files[{rel!r}] must hold sha256 (a string), bytes and rows (integers), got "
                          f"{entry!r} (BBR024)")
                continue
            facts = self.facts.get(f"{WAVEGUIDES}/{rel}")
            if facts is None:                        # missing or refused: reported with its directory
                continue
            for key, got in (("sha256", facts.sha256), ("bytes", facts.size), ("rows", facts.rows)):
                if entry[key] != got:
                    self.fail(f"{path}: files[{rel!r}].{key} is {entry[key]!r}, the file gives {got!r} (BBR024)")

    def _csvs(self, stem: _Stem) -> None:
        """Each CSV's header (when its sidecar is schema 2) and the Freq of its first data row."""
        for e in sorted(stem.dirs):
            for csv in CSV_NAMES:
                rel = f"{WAVEGUIDES}/{stem.name}_Ephi={e}/{csv}"
                facts, path = self.facts.get(rel), self.tree / rel
                if facts is None:
                    continue
                want = list(COLUMNS[csv])
                if stem.sidecar_data is not None and facts.header != want:
                    self.fail(f"{path}: its header {facts.header!r} is not the schema-2 {csv} columns {want!r} "
                              "(BBR013)")
                if facts.first_row is None:
                    self.fail(f"{path}: no data row after its header" + (
                        "" if facts.head_complete else f" in its first {HEAD_LIMIT >> 20} MiB"))
                    continue
                token = facts.first_row.split(b",", 1)[0].decode("latin-1")
                value = parse_frequency_ghz(token)
                if value is None or abs(value - stem.frequency) > FREQ_REL * stem.frequency:
                    self.fail(f"{path}: the Freq {token!r} of its first data row is not within 0.1 % of the "
                              f"directory's {stem.frequency:g} GHz (BBR009)")

    def _agreement(self) -> None:
        """Per dataset ID: no frequency twice (BBR011), and every sidecar's projection the lowest frequency's
        (G4-012)."""
        for dataset_id in sorted({s.dataset_id for s in self.stems.values()}):
            group = sorted((s for s in self.stems.values() if s.dataset_id == dataset_id),
                           key=lambda s: (s.frequency, s.name))
            for a, b in zip(group, group[1:]):
                if a.frequency == b.frequency:
                    self.fail(f"{dataset_id}: {a.name} and {b.name} name the same frequency, {a.frequency:g} GHz "
                              "(BBR011)")
            projected = [s for s in group if s.projection is not None]
            for other in projected[1:]:
                differ = invariant_diff(projected[0].projection, other.projection)
                if differ:
                    self.fail(f"{dataset_id}: the sidecars {projected[0].name}{SIDECAR_SUFFIX} and "
                              f"{other.name}{SIDECAR_SUFFIX} disagree on {', '.join(differ)} (BBR024; G4-012)")

    # --- SHA256SUMS -----------------------------------------------------------------------------------------------

    def _sums(self) -> None:
        path = self.tree / SUMS
        try:
            text = path.read_bytes().decode("utf-8")
        except (OSError, UnicodeDecodeError) as exc:
            self.fail(f"{path}: cannot be read as UTF-8 text ({exc})")
            return
        lines = text.split("\n")
        if lines[-1] == "":
            lines.pop()
        listed: dict[str, int] = {}
        for n, line in enumerate(lines, 1):
            match = _SUM_LINE.fullmatch(line)
            if match is None:
                self.fail(f"{path} line {n} is not '<sha256>  <path>': {line!r}")
                continue
            digest, rel = match["digest"], match["path"]
            if rel in listed:
                self.fail(f"{path} line {n} lists {rel} again (first at line {listed[rel]})")
                continue
            listed[rel] = n
            if rel in self.refused or any(rel.startswith(f"{r}/") for r in self.refused):
                continue
            facts = self.facts.get(rel)
            if facts is None:
                self.fail(f"{path} line {n} names {rel}, which is not a file of the tree")
                continue
            if facts.sha256 != digest:
                self.fail(f"{path} line {n}: {rel} has sha256 {facts.sha256}, the line records {digest}")
        for rel in sorted(set(self.facts) - set(listed), key=str.encode):
            self.fail(f"{self.tree / rel}: not listed in {SUMS}")
        self.listed = len(listed)

    # --- the verdict ----------------------------------------------------------------------------------------------

    def datasets(self) -> dict[str, dict]:
        out = {}
        for dataset_id in sorted({s.dataset_id for s in self.stems.values()}):
            group = sorted((s for s in self.stems.values() if s.dataset_id == dataset_id),
                           key=lambda s: (s.frequency, s.name))
            projected = next((s.projection for s in group if s.projection is not None), None)
            excitation = projected["excitation"] if projected is not None else None
            keys = (len(excitation["incident_phi_deg"]) * len(excitation["incident_theta_deg"])
                    if excitation is not None else None)
            out[dataset_id] = {"frequencies_ghz": [s.frequency for s in group], "stems": [s.name for s in group],
                               "incidence_keys": keys}
        return out


def _refused(tree: Path, problem: str) -> int:
    print(_one_line(f"assemble: refused: --verify {tree}: {problem}"), file=sys.stderr)
    return 2


def _inside(path: Path, tree: Path) -> bool:
    """Whether the directory ``path`` is the tree or lies under it: some ancestor of it, as written or resolved, is the
    tree's directory by (device, inode), so a case variant, a symlink, ``..`` or a bind mount of the tree is caught,
    which no comparison of path strings sees."""
    root = os.stat(tree)
    for candidate in (Path(os.path.abspath(path)), Path(os.path.realpath(path))):
        for ancestor in (candidate, *candidate.parents):
            try:
                st = os.stat(ancestor)
            except OSError:
                continue
            if (st.st_dev, st.st_ino) == (root.st_dev, root.st_ino):
                return True
    return False


def verify(tree: Path | str, *, report_path: Path | str | None = None) -> int:
    """Verify ``tree`` (the module docstring's rules); print the ``FAIL:`` lines and the verdict, and with
    ``report_path`` write them as JSON (``REPORT_SCHEMA``) to that new file. Returns 0 (``VERIFY PASS``), 1 (``VERIFY
    FAIL``, also when the report cannot be written) or 2 (refused: the path is not a readable directory, or the report
    path exists, has no directory or lies inside the tree)."""
    tree = Path(os.path.abspath(tree))
    if not os.path.exists(tree):
        return _refused(tree, "no such directory")
    if not os.path.isdir(tree):
        return _refused(tree, "not a directory")
    try:
        os.listdir(tree)
    except OSError as exc:
        return _refused(tree, f"cannot be read ({exc})")
    target = None
    if report_path is not None:
        target = Path(os.path.abspath(report_path))
        if os.path.lexists(target):
            return _refused(tree, f"--report-path {target} exists; a report is a new file")
        if not target.parent.is_dir():
            return _refused(tree, f"--report-path {target} has no directory")
        if _inside(target.parent, tree):
            return _refused(tree, f"--report-path {target} lies inside the tree, which --verify never writes into")
    _say(f"== bbsim verify: {tree} against BBRsim's loader rules, with the standard library only (HF-057)")
    check = _Verify(tree)
    check.run()
    datasets = check.datasets()
    for problem in check.problems:
        _say(f"FAIL: {problem}")
    if check.problems:
        verdict = f"VERIFY FAIL: {len(check.problems)} problem(s) in {tree}: {check.problems[0]}"
    else:
        _say(f"{SUMS}: {check.listed} files, each verifying; every file under {WAVEGUIDES}/ and {TREE_MANIFEST} "
             "listed")
        verdict = f"VERIFY PASS: {len(check.stems)} datasets -> {tree}; " + "; ".join(
            f"{dataset_id}: {', '.join(f'{f:g}' for f in d['frequencies_ghz'])} GHz x {d['incidence_keys']} keys"
            for dataset_id, d in datasets.items())
    if target is not None:
        record = {"schema": REPORT_SCHEMA, "created_utc": datetime.now(timezone.utc).isoformat(), "tree": str(tree),
                  "passed": not check.problems, "problems": check.problems, "datasets": datasets, "verdict": verdict}
        created = False
        try:
            with open(target, "x", encoding="utf-8") as fh:
                created = True
                fh.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
        except OSError as exc:
            if created:                                # never a half-written report
                target.unlink(missing_ok=True)
            _say(f"the report could not be written to {target} ({exc})")
            _say("VERIFY FAIL: the report could not be written")
            return 1
        _say(f"report: {target}")
    _say(verdict)
    return 1 if check.problems else 0
