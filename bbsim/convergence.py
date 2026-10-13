"""Whether a solve met its convergence target, from AEDT's exported convergence table (CDX-025).

AEDT 2025 R2 writes, per design variation, a header such as (HPC, job_tight1, crack1 Ephi=0)::

    Completed : 9
    Maximum   : 20
    Criterion : Delta Mag Energy
    Target    : 0.005
    Current   : 0.0036596
    Converged : Yes

followed by the pass table. Only an explicit ``Converged : Yes`` with the achieved delta at or below the
configured target qualifies; a missing or unreadable table never does.
"""
from __future__ import annotations

import math
import re

_FIELD = re.compile(r"^\s*(Completed|Maximum|Criterion|Target|Current|Converged)\s*:\s*(.*?)\s*$", re.MULTILINE)
CRITERION = "Delta Mag Energy"


def _number(text: str | None) -> float | None:
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def parse_convergence(text: str) -> dict:
    """The header fields of one variation's convergence table; a field absent or unreadable is None."""
    fields = dict(_FIELD.findall(text or ""))
    completed, maximum = _number(fields.get("Completed")), _number(fields.get("Maximum"))
    return {
        "criterion": fields.get("Criterion"),
        "target": _number(fields.get("Target")),
        "achieved": _number(fields.get("Current")),
        "passes_completed": int(completed) if completed is not None else None,
        "max_passes": int(maximum) if maximum is not None else None,
        "converged": fields.get("Converged"),
    }


def qualify(text: str, target: float) -> dict:
    """The parsed table plus ``qualified`` and the ``reasons`` it does not qualify against ``target``."""
    parsed = parse_convergence(text)
    reasons = []
    if parsed["converged"] is None:
        reasons.append("no 'Converged' line (table missing or unreadable)")
    elif parsed["converged"].lower() != "yes":
        reasons.append(f"AEDT reports Converged : {parsed['converged']}")
    if parsed["achieved"] is None:
        reasons.append("no achieved delta")
    elif parsed["achieved"] > target:
        reasons.append(f"achieved delta {parsed['achieved']:g} above the target {target:g}")
    if parsed["target"] is not None and not math.isclose(parsed["target"], target, rel_tol=1e-9):
        reasons.append(f"table target {parsed['target']:g} differs from the configured {target:g}")
    if parsed["criterion"] is not None and parsed["criterion"] != CRITERION:
        reasons.append(f"criterion {parsed['criterion']!r}, expected {CRITERION!r}")
    return {**parsed, "qualified": not reasons, "reasons": reasons}
