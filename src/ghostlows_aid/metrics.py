"""Clean-run summaries on true glucose, for settings calibration (task 6).

Consensus ranges come from paper 1 (ghostlows.metrics.consensus, Battelino 2019), computed on the
1-min plasma truth as paper 1's calibration does. Faulted-versus-clean outcomes are not computed
anywhere in this package before protocol-v1 (rule 11).
"""

from __future__ import annotations

import numpy as np
from ghostlows.metrics.consensus import TAR_THRESHOLDS, TBR_THRESHOLDS, TIR_RANGE

from .loop import LoopResult


def summarise(r: LoopResult, first_minute: int = 0) -> dict:
    """Time in ranges (% of minutes), mean, nadir, insulin and rescue carbs per day."""
    g = r.truth[first_minute:]
    days = g.size / 1440
    events = r.events[r.events["minute"] >= first_minute]
    lo, hi = TIR_RANGE
    return {
        "days": days,
        "valid": bool(r.valid),
        "tir_pct": 100 * float(np.mean((g >= lo) & (g <= hi))),
        "tbr70_pct": 100 * float(np.mean(g < TBR_THRESHOLDS[0])),
        "tbr54_pct": 100 * float(np.mean(g < TBR_THRESHOLDS[1])),
        "tar180_pct": 100 * float(np.mean(g > TAR_THRESHOLDS[0])),
        "tar250_pct": 100 * float(np.mean(g > TAR_THRESHOLDS[1])),
        "mean": float(np.mean(g)),
        "nadir": float(np.min(g)),
        "insulin_u_per_day": float(r.insulin[first_minute:].sum() / days),
        "rescue_g_per_day": float(events.loc[events["kind"] == "rescue", "value"].sum() / days),
        "rescues_per_day": float((events["kind"] == "rescue").sum() / days),
        "controller_errors": int(r.controller_errors),
    }
