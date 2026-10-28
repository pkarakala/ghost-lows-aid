"""The online detectors equal paper 1's batch detectors, reading by reading (CGM traces only)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from ghostlows.alarms.detectors import matched_filter_score
from ghostlows.faults.suite import FaultConfig
from ghostlows.retro.baysal import BAYSAL_SETS, baysal_flags

from ghostlows_aid.detectors import MF_THRESHOLD, BaysalOnline, MatchedFilterOnline, make_detector
from ghostlows_aid.faults_online import FaultSchedule


def faulted_trace(seed: int, days: int = 4) -> tuple[np.ndarray, np.ndarray]:
    """A CGM trace with PISAs every night, dropouts, frozen values and spikes (paper 1's suite at
    severity 1), clipped to the device range, NaN where missing."""
    rng = np.random.default_rng(seed)
    m = np.arange(days * 288) * 5.0
    clean = 140 + 60 * np.sin(m / 240.0) + 25 * np.sin(m / 37.0) + rng.normal(0, 4, m.size)
    schedule = FaultSchedule.draw(seed, days, FaultConfig(pisa_fraction=1.0, severity=1.0))
    obs = np.array([schedule.apply(k, float(v))[0] for k, v in enumerate(clean)])
    return m, np.clip(obs, 40, 400)


@pytest.mark.parametrize("seed", [1, 2, 3])
@pytest.mark.parametrize("params", list(BAYSAL_SETS))
def test_baysal_online_equals_batch(params, seed):
    m, y = faulted_trace(seed)
    times = pd.Timestamp("2026-01-05") + pd.to_timedelta(m, unit="min")
    batch = baysal_flags(times, y, BAYSAL_SETS[params])
    det = BaysalOnline(params)
    online = np.array([det.flag(t, v) for t, v in zip(m, y, strict=True)], dtype=float)
    assert np.array_equal(online, batch)
    assert batch.sum() > 0  # the trace does trigger the rules


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_matched_filter_online_equals_batch(seed):
    m, y = faulted_trace(seed)
    batch = matched_filter_score(y)
    det = MatchedFilterOnline()
    online = np.array([det.score(v) for v in y])
    assert np.array_equal(online, batch)
    flags = [MatchedFilterOnline().flag(0, v) for v in y[:3]]
    assert flags == [False, False, False]
    assert (batch > MF_THRESHOLD).sum() > 0


def test_warm_start_equals_watching_from_the_start():
    m, y = faulted_trace(4)
    j = 500
    full = BaysalOnline("trial")
    a = [full.flag(t, v) for t, v in zip(m, y, strict=True)]
    warmed = BaysalOnline("trial")
    warmed.warm(m[:j], y[:j])
    assert [warmed.flag(t, v) for t, v in zip(m[j:], y[j:], strict=True)] == a[j:]
    mf_full, mf_warm = MatchedFilterOnline(), MatchedFilterOnline()
    b = [mf_full.flag(0, v) for v in y]
    mf_warm.warm(m[:j], y[:j])
    assert [mf_warm.flag(0, v) for v in y[j:]] == b[j:]


def test_make_detector():
    assert make_detector("baysal").set_name == "trial"  # paper 1's pre-specified set
    assert make_detector("baysal_cautious").p == BAYSAL_SETS["cautious"]
    assert isinstance(make_detector("matched_filter"), MatchedFilterOnline)
    with pytest.raises(ValueError):
        make_detector("iforest")
