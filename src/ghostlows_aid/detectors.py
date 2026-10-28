"""Paper 1's causal PISA detectors, run one reading at a time as a gate in the loop.

`BaysalOnline` is paper 1's `ghostlows.retro.baysal.baysal_flags` [Baysal 2014; Cameron 2012]
rewritten as a state machine that takes one reading per call: the same Kalman steps (the same
numpy code), the same rules a-f, in the same order, so the flag of every reading equals the batch
flag of the stream up to and including it (tests/test_detectors.py). The batch function is
causal; the online form only avoids re-reading the whole history at every sample.

`MatchedFilterOnline` is paper 1's causal matched filter (`ghostlows.alarms.detectors.
matched_filter_score`, a causal adaptation of [Idi 2023]) on the last five readings, which is all
the score at a sample depends on; a reading is flagged when the score exceeds 0.5, paper 1's
untuned default threshold.

Both see the stream bound for the controller every 5 min (NaN when a reading is missing) and keep
no randomness. In a branch, a detector is warmed with the trunk's readings before t0: it has been
watching since the start; it gates only from t0 (episodes.py).
"""

from __future__ import annotations

from collections import deque

import numpy as np
from ghostlows.alarms.detectors import TEMPLATE_DIFFS, matched_filter_score
from ghostlows.retro.baysal import (
    BAYSAL_DEFAULT,
    BAYSAL_SETS,
    KF_DROPOUT_RESET_MIN,
    KF_MEASUREMENT_SD,
    KF_RATE_SD_PER_MIN,
    KF_RESET_TARGET,
    KF_RESET_TAU_MIN,
    BaysalParams,
)

MF_THRESHOLD = 0.5  # ghostlows.alarms.learned: "untuned default (suppress when gate > 0.5)"


class BaysalOnline:
    """baysal_flags (ghostlows/retro/baysal.py:53-138), one reading per call."""

    name = "baysal"

    def __init__(self, params: str | BaysalParams = BAYSAL_DEFAULT):
        self.set_name = params if isinstance(params, str) else "custom"
        self.p = BAYSAL_SETS[params] if isinstance(params, str) else params
        self.n = 0  # valid readings so far (the batch index k)
        self.x = None
        self.cov = None
        self.t_prev = None
        self.g_prev = None
        self.roc = deque([np.nan] * 5, maxlen=5)  # roc[k-4] .. roc[k]
        self.curv = deque([np.nan] * 2, maxlen=2)  # curv[k-1], curv[k]
        self.in_pisa, self.n_pisa, self.t_start = False, 0, None

    # the batch code's own predict and update (baysal.py:72-87), unchanged
    @staticmethod
    def _predict(x, cov, dt):
        q = KF_RATE_SD_PER_MIN**2
        f = np.array([[1.0, dt], [0.0, 1.0]])
        qd = q * np.array([[dt**3 / 3, dt**2 / 2], [dt**2 / 2, dt]])
        x = f @ x
        if dt > KF_DROPOUT_RESET_MIN:  # long dropout: relax toward the reset target
            w = 1 - np.exp(-(dt - KF_DROPOUT_RESET_MIN) / KF_RESET_TAU_MIN)
            x = np.array([(1 - w) * x[0] + w * KF_RESET_TARGET, (1 - w) * x[1]])
        return x, f @ cov @ f.T + qd

    @staticmethod
    def _update(x, cov, z):
        r = KF_MEASUREMENT_SD**2
        h = np.array([1.0, 0.0])
        s = h @ cov @ h + r
        k = cov @ h / s
        return x + k * (z - x[0]), cov - np.outer(k, h @ cov)

    def _within(self, a, b):
        band = (1 - self.p.ratio_band, 1 + self.p.ratio_band)
        if b == 0 or np.isnan(a) or np.isnan(b):
            return False
        return band[0] <= a / b <= band[1]

    def flag(self, minute: float, value: float) -> bool:
        if value != value:  # missing readings are skipped, as the batch code drops them
            return False
        t, g, p, k = float(minute), float(value), self.p, self.n
        if k == 0:
            self.x = np.array([g, 0.0])
            self.cov = np.diag([KF_MEASUREMENT_SD**2, 1.0])
            roc = np.nan
        else:
            roc = (g - self.g_prev) / (t - self.t_prev)
        curv = (roc - self.roc[-1]) / (t - self.t_prev) if k >= 2 else np.nan
        self.roc.append(roc)
        self.curv.append(curv)
        r = self.roc  # r[-1] = roc[k], r[-2] = roc[k-1], ...
        if k > 0:
            self.x, self.cov = self._predict(self.x, self.cov, t - self.t_prev)
        estimate = self.x[0]
        if not self.in_pisa:
            rule_a = k >= 2
            rule_b = rule_a and r[-1] < p.g_in
            rule_c = rule_b and (r[-2] >= 0 or r[-1] / r[-2] > p.g_ratio)
            if rule_c:
                self.in_pisa, self.n_pisa, self.t_start = True, 1, t
        else:
            self.n_pisa += 1
            rule_d = self.t_prev - self.t_start >= p.max_pisa_min
            rule_e = t - self.t_prev >= p.max_dropout_min
            w = self._within
            rule_f = r[-1] > p.g_out and (
                estimate < g
                or (self.n_pisa >= 3 and self.curv[-1] < 0 and self.curv[-2] < 0)
                or (k >= 3 and w(r[-4], r[-2]) and w(r[-4], r[-1]))
                or (k >= 4 and w(r[-5], r[-2]) and w(r[-5], r[-1]))
            )
            if rule_d or rule_e or rule_f:
                self.in_pisa = False
        flagged = self.in_pisa
        if not flagged:
            self.x, self.cov = self._update(self.x, self.cov, g)  # attenuated readings unused
        self.t_prev, self.g_prev, self.n = t, g, k + 1
        return bool(flagged)

    def warm(self, minutes, values) -> None:
        for m, v in zip(minutes, values, strict=True):
            self.flag(m, v)


class MatchedFilterOnline:
    """matched_filter_score at the newest reading, flagged above MF_THRESHOLD."""

    name = "matched_filter"

    def __init__(self, threshold: float = MF_THRESHOLD):
        self.threshold = threshold
        self.window = deque(maxlen=len(TEMPLATE_DIFFS) + 1)

    def score(self, value: float) -> float:
        self.window.append(float(value))
        if len(self.window) < self.window.maxlen:
            return 0.0
        return float(matched_filter_score(np.array(self.window))[-1])

    def flag(self, minute: float, value: float) -> bool:
        return self.score(value) > self.threshold

    def warm(self, minutes, values) -> None:
        for v in values:
            self.window.append(float(v))


def make_detector(spec: str):
    """ "baysal" (paper 1's pre-specified "trial" set), "baysal_<set>" or "matched_filter"."""
    if spec == "matched_filter":
        return MatchedFilterOnline()
    if spec == "baysal":
        return BaysalOnline(BAYSAL_DEFAULT)
    if spec.startswith("baysal_") and spec[7:] in BAYSAL_SETS:
        return BaysalOnline(spec[7:])
    raise ValueError(f"unknown detector {spec!r}")
