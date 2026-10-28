"""C1, the pre-specified primary comparison (memo §9; estimation only): per-day harm of paper 1's
five fault types at the typical level in the alarm layer and the insulin layer, Kendall's tau-b
between the two per-fault rankings, compression's share of each layer's total and the difference
between the layers' shares, with the three sensitivity analyses.

Per-day harm of a fault = per-event effect x its daily rate (`RATES`, the one place for rates;
the person-year harm in `secondary` reads the same dict). Every quantity is a function of the
replicate means of one joint patient bootstrap (one resampling index per replicate shared by all
faults, both layers and every sensitivity analysis), so shares and the tau-b are computed per
replicate. A share is estimated only when its layer's total is positive in at least 97.5% of the
replicates; the layers are *discordant* when the 95% CI for the difference of compression's
shares excludes 0 (both shares estimable).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from ghostlows.faults.suite import FaultConfig

from . import effects as ef

_P1 = FaultConfig()
S_MODERATE = 0.5  # paper 1's moderate severity s = 0.5 (memo §6 arms 3, 6, 8 and §9 C1)

# SOURCE: Lobo 2026 (1,260 PISAs with a sensor minimum below 85 mg/dL in 94,201 h)
LOBO_PER_DAY = 0.321

RATES = {
    # Events per person-day. Paper 1 draws its non-PISA faults at full-severity rates and keeps
    # each with probability s, so the rate at s = 0.5 is half of FaultConfig's.
    # SOURCE: Reddy 2026 via paper 1 (a compression artifact on 86% of nights, one per night).
    "pisa_typical": 0.86,
    # SOURCE: paper 1's headline setting, PISA on every night (memo §9 C1 sensitivity analysis 1).
    "pisa_every_night": 1.0,
    # SOURCE: paper 1 FaultConfig.dropout_per_day x s at s = 0.5 (paper 1 marks the rate
    # ASSUMPTION: no source [R4]).
    "dropout_short": _P1.dropout_per_day * S_MODERATE,
    # SOURCE: paper 1 FaultConfig.stale_per_day x s at s = 0.5 (ASSUMPTION in paper 1).
    "frozen_typical": _P1.stale_per_day * S_MODERATE,
    # SOURCE: paper 1 FaultConfig.spike_per_day x s at s = 0.5 (ASSUMPTION in paper 1).
    "spike_typical": _P1.spike_per_day * S_MODERATE,
    # ASSUMPTION (memo §9 C1): drift as a continuous state, two 12-h windows per day, half in
    # each direction, as paper 1 draws the drift sign with probability 1/2 (faults/suite.py).
    "drift_low": 2 * 0.5,
    "drift_high": 2 * 0.5,
    # SOURCE: Vettoretti 2023, G6 trial: 0.13 to 0.24 gaps per day (the long-gap model, arm 4).
    "dropout_long_low": 0.13,
    "dropout_long_high": 0.24,
    # Compression for the person-year harm (memo §10 Q7), as a fraction of nights with a PISA:
    # paper 1's 0.86 reproduces Lobo 2026's 0.321 detected PISAs per day, so a published detected
    # rate r maps to 0.86 * r / 0.321 nights.
    # ASSUMPTION (memo §10 Q7): the published detection definitions are comparable with Lobo's.
    # SOURCE: Lobo 2026 (0.321 per day) via paper 1's calibration (central).
    "pisa_nights_central": 0.86,
    # SOURCE: LoopEval 2026 (0.03 per person-day, the lowest published).
    "pisa_nights_low": 0.86 * 0.03 / LOBO_PER_DAY,
    # SOURCE: Idi 2025 (about 0.33 per day, the highest published).
    "pisa_nights_high": 0.86 * 0.33 / LOBO_PER_DAY,
}

LAYERS = {  # memo §9 C1: layer -> (controller, outcome)
    "alarm": ("open_loop", "trend_false_alarms"),  # excess false alarms per day
    "insulin": ("oref0_tb", "risk_weighted_min"),  # excess risk-weighted minutes per day
}
COMPRESSION = "compression"
ESTIMABLE_SHARE = 0.975  # memo §9: total positive in at least 97.5% of replicates
TIE_DECIMALS = 9  # ASSUMPTION: per-day harms equal to 1e-9 rank as ties (rounding, not effects)


@dataclass(frozen=True)
class Spec:
    """One C1 analysis: fault -> {arm: events per day}, and the layers."""

    name: str
    faults: Mapping[str, Mapping[str, float]]
    layers: Mapping[str, tuple[str, str]] = field(default_factory=lambda: dict(LAYERS))


def _main_faults() -> dict[str, dict[str, float]]:
    return {
        COMPRESSION: {"pisa_typical": RATES["pisa_typical"]},
        "dropout": {"dropout_short": RATES["dropout_short"]},
        "frozen": {"frozen_typical": RATES["frozen_typical"]},
        "spikes": {"spike_typical": RATES["spike_typical"]},
        "drift": {"drift_low": RATES["drift_low"], "drift_high": RATES["drift_high"]},
    }


def specs() -> list[Spec]:
    """The main analysis and the sensitivity analyses of memo §9 C1 (the long-gap analysis at
    both ends of the G6 range)."""
    main = _main_faults()
    every = {**main, COMPRESSION: {"pisa_typical": RATES["pisa_every_night"]}}
    threshold = {**LAYERS, "alarm": ("open_loop", "threshold_false_alarms")}
    out = [Spec("main", main), Spec("every_night", every), Spec("threshold_alarm", main, threshold)]
    for end in ("low", "high"):
        gaps = {**main, "dropout": {"dropout_long": RATES[f"dropout_long_{end}"]}}
        out.append(Spec(f"long_gaps_{end}", gaps))
    return out


def kendall_tau_b(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Kendall's tau-b over the last axis (vectorised over leading axes); NaN where a ranking is
    all ties or a value is missing."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    i, j = np.triu_indices(x.shape[-1], 1)
    dx, dy = np.sign(x[..., i] - x[..., j]), np.sign(y[..., i] - y[..., j])
    with np.errstate(invalid="ignore", divide="ignore"):
        den = np.sqrt((dx != 0).sum(-1) * (dy != 0).sum(-1))
        tau = (dx * dy).sum(-1) / den
    bad = np.isnan(dx).any(-1) | np.isnan(dy).any(-1) | (den == 0)
    return np.where(bad, np.nan, tau)


class _Theta:
    """Per-event theta and its replicates per (controller, arm, outcome), computed once."""

    def __init__(self, effects: pd.DataFrame, boot: ef.Bootstrap):
        self.effects, self.boot, self.cache = effects, boot, {}

    def __call__(self, controller: str, arm: str, outcome: str) -> tuple[float, np.ndarray, int]:
        key = (controller, arm, outcome)
        if key not in self.cache:
            e = self.effects
            sel = e[(e["controller"] == controller) & (e["arm"] == arm)]
            nan = (np.nan, np.full(self.boot.n_boot, np.nan), 0)
            if sel.empty or outcome not in sel:
                self.cache[key] = nan
            else:
                est, draws = ef.estimate(sel, [outcome], ef.BY_ARM, self.boot, keep_draws=True)
                if est.empty:
                    self.cache[key] = nan
                else:
                    row = est.iloc[0]
                    self.cache[key] = (float(row["theta"]), draws[:, 0], int(row["n_patients"]))
        return self.cache[key]


def _ci(draws: np.ndarray) -> tuple[float, float]:
    lo, hi = ef.percentile_ci(np.asarray(draws, float)[:, None])
    return float(lo[0]), float(hi[0])


def run(
    effects: pd.DataFrame, boot: ef.Bootstrap, analyses: Sequence[Spec] | None = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(per-fault table, summary table) of C1 from the paired effects (invalid rule applied);
    `boot` is the confirmatory bootstrap (memo §10: n_boot = 10,000)."""
    theta = _Theta(effects, boot)
    per_fault, summary = [], []
    for spec in analyses or specs():
        faults = list(spec.faults)
        layer_est, layer_draws = {}, {}
        for layer, (controller, outcome) in spec.layers.items():
            est = np.zeros(len(faults))
            draws = np.zeros((boot.n_boot, len(faults)))
            for k, f in enumerate(faults):
                n_pat = []
                for arm, rate in spec.faults[f].items():
                    t, d, n = theta(controller, arm, outcome)
                    est[k] += rate * t
                    draws[:, k] += rate * d
                    n_pat.append(n)
                lo, hi = _ci(draws[:, k])
                per_fault.append(
                    {
                        "analysis": spec.name,
                        "layer": layer,
                        "controller": controller,
                        "outcome": outcome,
                        "fault": f,
                        "arms": "+".join(spec.faults[f]),
                        "rates_per_day": "+".join(f"{r:g}" for r in spec.faults[f].values()),
                        "per_day": est[k],
                        "ci_lo": lo,
                        "ci_hi": hi,
                        "n_patients": min(n_pat),
                    }
                )
            layer_est[layer], layer_draws[layer] = est, draws
        rec = {"analysis": spec.name, "n_boot": boot.n_boot}
        a, b = list(spec.layers)[:2]
        tie = {k: np.round(v, TIE_DECIMALS) for k, v in {**layer_est}.items()}
        tie_d = {k: np.round(v, TIE_DECIMALS) for k, v in layer_draws.items()}
        rec["tau_b"] = float(kendall_tau_b(tie[a], tie[b]))
        rec["tau_b_lo"], rec["tau_b_hi"] = _ci(kendall_tau_b(tie_d[a], tie_d[b]))
        c = faults.index(COMPRESSION)
        share_draws = {}
        for layer in (a, b):
            total, total_d = layer_est[layer].sum(), layer_draws[layer].sum(1)
            p_pos = float(np.mean(total_d > 0))
            ok = p_pos >= ESTIMABLE_SHARE
            with np.errstate(invalid="ignore", divide="ignore"):
                share = layer_est[layer][c] / total
                sd = layer_draws[layer][:, c] / total_d
            share_draws[layer] = sd if ok else None
            lo, hi = _ci(sd) if ok else (np.nan, np.nan)
            rec.update(
                {
                    f"total_{layer}": total,
                    f"p_total_positive_{layer}": p_pos,
                    f"share_estimable_{layer}": ok,
                    f"share_{layer}": share if ok else np.nan,
                    f"share_{layer}_lo": lo,
                    f"share_{layer}_hi": hi,
                }
            )
        if share_draws[a] is not None and share_draws[b] is not None:
            diff = share_draws[a] - share_draws[b]
            rec["share_diff"] = rec[f"share_{a}"] - rec[f"share_{b}"]
            rec["share_diff_lo"], rec["share_diff_hi"] = _ci(diff)
            rec["discordant"] = bool(rec["share_diff_lo"] > 0 or rec["share_diff_hi"] < 0)
        else:
            rec.update(share_diff=np.nan, share_diff_lo=np.nan, share_diff_hi=np.nan)
            rec["discordant"] = pd.NA
        summary.append(rec)
    pf = pd.DataFrame(per_fault)
    if not pf.empty:  # rank 1 = most per-day harm, within analysis and layer
        pf["rank"] = pf.groupby(["analysis", "layer"])["per_day"].rank(ascending=False)
    out = pd.DataFrame(summary)
    if not out.empty:
        out["discordant"] = out["discordant"].astype("boolean")
    return pf, out


__all__ = [
    "COMPRESSION",
    "ESTIMABLE_SHARE",
    "LAYERS",
    "RATES",
    "S_MODERATE",
    "TIE_DECIMALS",
    "Spec",
    "kendall_tau_b",
    "run",
    "specs",
]
