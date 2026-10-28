"""The confirmatory family H1-H2 on `oref0_tb` (memo §9): estimates, bootstrap p-values, Holm's
step-down and the four-way rule (meaningful harm / negligible / opposite direction /
inconclusive).

Exactly as memo §9: for hypothesis k, a_k and b_k are the shares of bootstrap replicates of theta
at or below 0 and at or above 0, each (1 + count) / (B + 1); p_k = min(1, 2 min(a_k, b_k)). Holm
at familywise alpha = 0.05 over the two: the smaller p is compared with 0.025 and, only if it
passes, the larger with 0.05. "Meaningful harm" needs that rejection, a positive estimate and the
estimate >= SESOI, and its CI is reported at the Holm-adjusted level of its step (97.5% or 95%).
"Negligible" if the 90% CI lies within (-SESOI, SESOI) (TOST, each at alpha = 0.05, not adjusted);
"opposite direction" if the Holm test rejects with an estimate <= -SESOI (an unexpected direction,
never reported as harm); otherwise "inconclusive". Co-reported outcomes get 95% CIs, uncorrected
and labelled so.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import effects as ef
from .schema import ARMS, PRIMARY_CONTROLLER

ALPHA = 0.05  # memo §9: familywise alpha (Holm)
TOST_ALPHA = 0.05  # memo §9: each TOST at alpha = 0.05, i.e. the 90% CI
HARM, NEGLIGIBLE, OPPOSITE, INCONCLUSIVE, NOT_ESTIMABLE = (
    "meaningful harm",
    "negligible",
    "opposite direction",
    "inconclusive",
    "not estimable",
)


@dataclass(frozen=True)
class Hypothesis:
    id: str
    arm: str
    outcome: str
    sesoi: float
    co_reported: tuple[str, ...]
    controller: str = PRIMARY_CONTROLLER


HYPOTHESES = (
    # SOURCE: Battelino 2019 (a 5-point change in time in range, about 72 min per day, is
    # clinically meaningful). ASSUMPTION (memo §9): the per-event reading; CR.
    Hypothesis(
        "H1", "pisa_typical", "minutes_above_180", 72.0, ("minutes_below_54", "minutes_below_70")
    ),
    # ASSUMPTION (memo §9): 15 min, about 1% of a day, between the consensus targets of < 1% below
    # 54 and < 4% below 70 mg/dL [Battelino 2019]; CR.
    Hypothesis("H2", "interference_large", "minutes_below_70", 15.0, ("minutes_below_54",)),
)


def bootstrap_p(draws: np.ndarray) -> float:
    """Two-sided bootstrap p-value of memo §9 from the finite replicates of theta."""
    d = np.asarray(draws, float)
    d = d[np.isfinite(d)]
    if d.size == 0:
        return math.nan
    a = (1 + np.count_nonzero(d <= 0)) / (d.size + 1)
    b = (1 + np.count_nonzero(d >= 0)) / (d.size + 1)
    return float(min(1.0, 2 * min(a, b)))


def holm(pvalues: Mapping[str, float], alpha: float = ALPHA) -> pd.DataFrame:
    """Holm's step-down over the hypotheses: per id its step (1 = smallest p; ties by id, NaN
    last), the step's level alpha / (m - step + 1) and whether it is rejected (p <= that level,
    and every earlier step rejected)."""
    m = len(pvalues)
    order = sorted(pvalues, key=lambda k: (not np.isfinite(pvalues[k]), pvalues[k], k))
    rows, still = [], True
    for step, k in enumerate(order, 1):
        level = alpha / (m - step + 1)
        p = pvalues[k]
        rejected = bool(still and np.isfinite(p) and p <= level)
        still = rejected
        rows.append({"id": k, "holm_step": step, "holm_alpha": level, "rejected": rejected})
    return pd.DataFrame(rows).set_index("id")


def decide(theta: float, rejected: bool, ci90: tuple[float, float], sesoi: float) -> str:
    """The four-way rule of memo §9, harm first."""
    lo, hi = ci90
    if not (np.isfinite(theta) and np.isfinite(lo) and np.isfinite(hi)):
        return NOT_ESTIMABLE
    if rejected and theta > 0 and theta >= sesoi:
        return HARM
    if rejected and theta < 0 and theta <= -sesoi:
        return OPPOSITE
    if -sesoi < lo and hi < sesoi:
        return NEGLIGIBLE
    return INCONCLUSIVE


def _rows(effects: pd.DataFrame, h: Hypothesis, outcomes: Sequence[str], boot) -> tuple:
    sel = effects[(effects["controller"] == h.controller) & (effects["arm"] == h.arm)]
    est, draws = ef.estimate(sel, list(outcomes), ef.BY_ARM, boot, keep_draws=True)
    return est.set_index("outcome"), draws


def run(
    effects: pd.DataFrame,
    boot: ef.Bootstrap,
    hypotheses: Sequence[Hypothesis] = HYPOTHESES,
    alpha: float = ALPHA,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(confirmatory table, co-reported table) from the paired effects (invalid rule applied).
    `boot` is the confirmatory bootstrap (memo §10: n_boot = 10,000)."""
    found, co = {}, []
    for h in hypotheses:
        est, draws = _rows(effects, h, [h.outcome, *h.co_reported], boot)
        cols = list(est.index)
        main = (
            draws[:, cols.index(h.outcome)] if h.outcome in cols else np.full(boot.n_boot, np.nan)
        )
        row = est.loc[h.outcome] if h.outcome in cols else None
        found[h.id] = (h, row, main)
        for o in h.co_reported:
            if o in cols:
                r = est.loc[o]
                co.append(
                    {
                        "id": h.id,
                        "controller": h.controller,
                        "arm": h.arm,
                        "outcome": o,
                        "theta": r["theta"],
                        "ci_lo": r["ci_lo"],
                        "ci_hi": r["ci_hi"],
                        "level": ef.LEVEL,
                        "n_patients": r["n_patients"],
                        "label": "co-reported, 95% CI, uncorrected",
                    }
                )
    p = {k: bootstrap_p(v[2]) for k, v in found.items()}
    steps = holm(p, alpha)
    out = []
    for k, (h, row, main) in found.items():
        step = steps.loc[k]
        level = 1 - step["holm_alpha"]
        lo, hi = ef.percentile_ci(main[:, None], level)
        lo90, hi90 = ef.percentile_ci(main[:, None], 1 - 2 * TOST_ALPHA)
        theta = float(row["theta"]) if row is not None else math.nan
        rec = {
            "id": k,
            "controller": h.controller,
            "slot": ARMS[h.arm].slot,
            "arm": h.arm,
            "outcome": h.outcome,
            "sesoi": h.sesoi,
            "theta": theta,
            "p_value": p[k],
            "holm_step": int(step["holm_step"]),
            "holm_alpha": float(step["holm_alpha"]),
            "rejected": bool(step["rejected"]),
            "ci_level": level,
            "ci_lo": float(lo[0]),
            "ci_hi": float(hi[0]),
            "ci90_lo": float(lo90[0]),
            "ci90_hi": float(hi90[0]),
            "n_boot": boot.n_boot,
        }
        rec["decision"] = decide(theta, rec["rejected"], (rec["ci90_lo"], rec["ci90_hi"]), h.sesoi)
        for c in ("n_patients", "n_events", "within_se", "between_sd", "mc_se", "mc_flag"):
            rec[c] = row[c] if row is not None else math.nan
        out.append(rec)
    return pd.DataFrame(out), pd.DataFrame(co)


__all__ = [
    "ALPHA",
    "HARM",
    "HYPOTHESES",
    "INCONCLUSIVE",
    "NEGLIGIBLE",
    "NOT_ESTIMABLE",
    "TOST_ALPHA",
    "Hypothesis",
    "bootstrap_p",
    "decide",
    "holm",
    "run",
]
