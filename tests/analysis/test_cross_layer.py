"""C1: per-day harm, tau-b, compression's shares, estimability, discordance, joint bootstrap."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.stats import kendalltau

from ghostlows_aid.analysis import cross_layer as cl
from ghostlows_aid.analysis import effects as ef
from ghostlows_aid.analysis.synthetic import PATIENTS, Design, synthetic_outcomes

C1_ARMS = (
    "pisa_typical",
    "dropout_short",
    "dropout_long",
    "frozen_typical",
    "spike_typical",
    "drift_low",
    "drift_high",
)
ALARM = {"pisa_typical": 1.0, "dropout_short": 0.1, "dropout_long": 0.5, "frozen_typical": 0.2}
ALARM |= {"spike_typical": 0.3, "drift_low": 0.05, "drift_high": 0.05}
INSULIN = {"pisa_typical": 10.0, "dropout_short": 20.0, "dropout_long": 50.0}
INSULIN |= {"frozen_typical": 5.0, "spike_typical": 2.0, "drift_low": 30.0, "drift_high": 10.0}


def c1(effects=None, effect_fn=None, n_boot=300, **kw):
    d = Design(
        controllers=("oref0_tb", "open_loop"),
        arms=C1_ARMS,
        n_seeds=2,
        effects=effects or {},
        effect_fn=effect_fn,
        seed=11,
        **kw,
    )
    kept, _ = ef.apply_invalid_rule(synthetic_outcomes(d))
    return cl.run(ef.paired_effects(kept), ef.boot_for(kept, n_boot))


def known_effects():
    e = {("open_loop", a, "trend_false_alarms"): v for a, v in ALARM.items()}
    e |= {("open_loop", a, "threshold_false_alarms"): 2 * v for a, v in ALARM.items()}
    e |= {("oref0_tb", a, "risk_weighted_min"): v for a, v in INSULIN.items()}
    return e


def per_day(theta, faults):
    return np.array([sum(r * theta[a] for a, r in arms.items()) for arms in faults.values()])


@pytest.mark.parametrize("seed", range(4))
def test_tau_b_equals_scipy_with_ties(seed):
    rng = np.random.default_rng(seed)
    x = rng.integers(0, 3, (20, 5)).astype(float)
    y = rng.integers(0, 4, (20, 5)).astype(float)
    mine = cl.kendall_tau_b(x, y)
    for k in range(20):
        ref = kendalltau(x[k], y[k], variant="b").statistic
        assert (np.isnan(ref) and np.isnan(mine[k])) or mine[k] == pytest.approx(ref)


def test_rates_follow_paper1_at_moderate_severity():
    assert cl.RATES["dropout_short"] == 2.0 and cl.RATES["frozen_typical"] == 1.0
    assert cl.RATES["spike_typical"] == 1.0 and cl.RATES["pisa_typical"] == 0.86
    assert cl.RATES["drift_low"] + cl.RATES["drift_high"] == 2.0  # two 12-h windows a day


def test_known_per_day_harm_shares_and_tau():
    pf, summary = c1(known_effects(), zero_noise_arms=True)
    s = summary.set_index("analysis")
    for spec in cl.specs():
        alarm_theta = {
            a: v * (2 if spec.name == "threshold_alarm" else 1) for a, v in ALARM.items()
        }
        a, i = per_day(alarm_theta, spec.faults), per_day(INSULIN, spec.faults)
        got = pf[pf["analysis"] == spec.name].set_index(["layer", "fault"])["per_day"]
        assert got.loc["alarm"].to_numpy() == pytest.approx(a)
        assert got.loc["insulin"].to_numpy() == pytest.approx(i)
        row = s.loc[spec.name]
        assert row["share_alarm"] == pytest.approx(a[0] / a.sum())
        assert row["share_insulin"] == pytest.approx(i[0] / i.sum())
        assert row["tau_b"] == pytest.approx(kendalltau(a, i, variant="b").statistic)
        assert row["discordant"]  # exact effects: a degenerate CI that excludes 0
    gaps = pf[(pf["analysis"] == "long_gaps_high") & (pf["fault"] == "dropout")]
    assert set(gaps["arms"]) == {"dropout_long"} and set(gaps["rates_per_day"]) == {"0.24"}


def test_share_not_estimable_without_positive_total():
    pf, summary = c1({}, zero_noise_arms=True)  # every arm identical to its clean arm
    assert (pf["per_day"] == 0).all()
    row = summary.set_index("analysis").loc["main"]
    assert not row["share_estimable_alarm"] and np.isnan(row["share_alarm"])
    assert row["p_total_positive_alarm"] == 0.0 and pd.isna(row["discordant"])
    assert np.isnan(row["share_diff"])


def test_share_estimability_threshold_with_noise():
    pf, summary = c1({})  # zero effects with seed noise: totals straddle 0
    row = summary.set_index("analysis").loc["main"]
    assert row["p_total_positive_alarm"] < cl.ESTIMABLE_SHARE
    assert not row["share_estimable_alarm"]


def test_one_bootstrap_index_for_both_layers():
    rng = np.random.default_rng(0)
    mult = {(p, a): rng.uniform(0.5, 2.0) for p in PATIENTS for a in C1_ARMS}

    def proportional(f):  # insulin = 10 x alarm, patient by patient: identical shares
        arm, ctrl = f["arm"].iloc[0], f["controller"].iloc[0]
        v = np.array([mult[(p, arm)] for p in f["patient_id"]])
        return {"trend_false_alarms": v} if ctrl == "open_loop" else {"risk_weighted_min": 10 * v}

    _, summary = c1({}, effect_fn=proportional, zero_noise_arms=True)
    row = summary.set_index("analysis").loc["main"]
    assert row["share_alarm_hi"] - row["share_alarm_lo"] > 0.01  # the shares do vary
    assert abs(row["share_diff_lo"]) < 1e-12 and abs(row["share_diff_hi"]) < 1e-12
    assert not pd.isna(row["discordant"]) and not row["discordant"]
    assert row["tau_b_lo"] == pytest.approx(1.0) and row["tau_b"] == pytest.approx(1.0)
