"""The alarm layer reproduces paper 1's own benchmark counts (BENCHMARK.md §1-4)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from ghostlows.alarms import GatedAlarm, ThresholdAlarm, TrendAlarm, view
from ghostlows.alarms.learned import BaysalDetector, MatchedFilterDetector
from ghostlows.bench.events import true_events
from ghostlows.bench.scoring import DEFAULT_WINDOWS, fire_and_score, score_alarm
from ghostlows.faults import FaultConfig, apply_faults
from ghostlows.sim.engine import DEFAULT_START, first_order_lag, sensor_params

from ghostlows_aid.analysis import alarm_layer
from ghostlows_aid.analysis.outcomes import ArmTrace
from ghostlows_aid.analysis.schema import ALARM_PREFIXES

PAPER1 = {
    "trend": TrendAlarm,
    "threshold": ThresholdAlarm,
    "gate_aggressive": lambda: GatedAlarm(BaysalDetector("aggressive")),
    "gate_trial": lambda: GatedAlarm(BaysalDetector("trial")),
    "gate_nominal": lambda: GatedAlarm(BaysalDetector("nominal")),
    "gate_cautious": lambda: GatedAlarm(BaysalDetector("cautious")),
    "gate_mf": lambda: GatedAlarm(MatchedFilterDetector()),
}
QUIET = 180  # minutes of flat glucose at the start of the run


def synthetic_truth(days: int = 2, seed: int = 0) -> np.ndarray:
    """1-min glucose: flat for 3 h, then oscillations and designed lows (level 1 and 2, slow and
    steep), some at night where paper 1 puts its PISAs. Not simulator output."""
    rng = np.random.default_rng(seed)
    m = np.arange(days * 1440, dtype=float)
    g = 150 + 45 * np.sin(2 * np.pi * (m - QUIET) / 310.0) * (m >= QUIET)
    for centre in rng.uniform(QUIET + 120, m[-1] - 120, size=14):
        nadir = rng.choice([45.0, 52.0, 62.0, 66.0])
        width = rng.uniform(12, 35)
        base = np.interp(centre, m, g)
        g -= (base - nadir) * np.exp(-0.5 * ((m - centre) / width) ** 2)
    return np.clip(g, 39.0, 400.0)


def paper1_run(days: int = 2, seed: int = 0):
    """(observed samples, fault instances, 1-min truth) in paper 1's format, faults from paper 1's
    apply_faults with a PISA on every night."""
    from simglucose.sensor.noise_gen import CGMNoise

    truth = synthetic_truth(days, seed)
    isf = first_order_lag(truth, 7.0)
    noise_gen = CGMNoise(sensor_params(), seed=seed + 1)
    k = np.arange(0, truth.size, 5)
    noise = np.array([next(noise_gen) for _ in k])
    times = (DEFAULT_START + pd.to_timedelta(k, unit="min")).astype("datetime64[ns]")
    keys = {"patient_id": "sim:adult#001", "seed": 0, "mix": "clean"}
    samples = pd.DataFrame(
        {
            **keys,
            "time": times,
            "isf_glucose": isf[k],
            "cgm_noise": noise,
            "glucose": np.clip(isf[k] + noise, 40, 400),
            "truth_glucose": truth[k],
        }
    )
    observed, faults = apply_faults(samples, FaultConfig(pisa_fraction=1.0), seed=seed)
    t1 = (DEFAULT_START + pd.to_timedelta(np.arange(truth.size), unit="min")).astype(
        "datetime64[ns]"
    )
    truth_df = pd.DataFrame({**keys, "time": t1, "truth_glucose": truth})
    return observed, faults, truth_df


def as_trace(observed, truth_df, pre_min: int, controller: str = "open_loop") -> ArmTrace:
    truth = truth_df["truth_glucose"].to_numpy(float)
    n5 = truth.size // 5
    return ArmTrace(
        controller=controller,
        patient_id="sim:adult#001",
        seed=0,
        slot="night",
        arm="clean",
        onset_minute=pre_min,
        truth=truth,
        basal_u=np.zeros(truth.size),
        bolus_u=np.zeros(truth.size),
        smb_u=np.zeros(truth.size),
        rescue_g=np.zeros(truth.size),
        cgm_display=observed["cgm_observed"].to_numpy(float),
        cgm_controller=np.full(n5, np.nan),
        iob_u=np.full(n5, np.nan),
        cob_g=np.full(n5, np.nan),
        reasons=[""] * n5,
        scheduled_basal_u_per_min=0.0,
        last_bolus_minute=np.nan,
        pre_min=pre_min,
        window_min=truth.size - pre_min,
        data_source="synthetic",
    )


@pytest.fixture(scope="module")
def run():
    return paper1_run()


def test_whole_run_counts_equal_paper1_score_alarm(run):
    observed, faults, truth = run
    lows = true_events(truth)
    assert (lows["level"] == 2).sum() >= 3 and (lows["level"] == 1).sum() >= 5
    trace = as_trace(observed, truth, pre_min=QUIET)
    mine = alarm_layer.score_traces([trace], ALARM_PREFIXES)
    n_false_total = 0
    for prefix, make in PAPER1.items():
        ref = score_alarm(make(), observed, lows, faults).iloc[0]
        # precondition: nothing fires and no low starts before t0 (the quiet start)
        fired = make().alarms(view(observed, "cgm_observed"))
        assert not len(fired) or fired["time"].min() >= DEFAULT_START + pd.Timedelta(minutes=QUIET)
        got = mine[("clean", prefix)]
        assert got["alarms"] == ref["n_alarms"], prefix
        assert got["false_alarms"] == ref["n_false"], prefix
        for level in (1, 2):
            assert got[f"l{level}_caught"] == ref[f"caught_l{level}"], prefix
            assert got[f"l{level}_caught"] + got[f"l{level}_missed"] == ref[f"events_l{level}"]
        n_false_total += got["false_alarms"]
    assert n_false_total > 0  # PISA nights give false alarms: the comparison is not vacuous


def test_window_counts_equal_paper1_classification_after_t0(run):
    observed, faults, truth = run
    lows = true_events(truth)
    pre = 1380  # t0 at 23:00 on day 1: alarms, lows and PISAs both before and after it
    trace = as_trace(observed, truth, pre_min=pre)
    t0 = DEFAULT_START + pd.Timedelta(minutes=pre)
    mine = alarm_layer.score_traces([trace], ["trend", "threshold", "gate_trial"])
    for prefix in ("trend", "threshold", "gate_trial"):
        alarm = PAPER1[prefix]()
        frame = view(observed, "cgm_observed")
        _, cls, scored = fire_and_score(
            alarm, frame, alarm.default_threshold, lows, faults, observed, DEFAULT_WINDOWS
        )
        after = cls[pd.DatetimeIndex(cls["time"]) >= t0]
        ev = scored[pd.DatetimeIndex(scored["onset"]) >= t0]
        got = mine[("clean", prefix)]
        assert got["alarms"] == len(after)
        assert got["false_alarms"] == int((~after["matched"].astype(bool)).sum())
        for level in (1, 2):
            e = ev[ev["level"] == level]
            assert got[f"l{level}_caught"] == int(e["caught"].sum())
            assert got[f"l{level}_missed"] == int((~e["caught"].astype(bool)).sum())
        assert len(cls) > len(after)  # some alarms fell before t0 and were not counted


def test_arms_scored_together_do_not_mix(run):
    observed, faults, truth = run
    a = as_trace(observed, truth, pre_min=QUIET)
    quiet = observed.assign(cgm_observed=150.0)
    b = as_trace(quiet, truth, pre_min=QUIET)
    b.arm = "pisa_typical"
    both = alarm_layer.score_traces([a, b], ["trend"])
    alone = alarm_layer.score_traces([a], ["trend"])
    assert both[("clean", "trend")] == alone[("clean", "trend")]
    # a display that never falls never alarms: every true low is missed
    assert both[("pisa_typical", "trend")]["alarms"] == 0
    assert both[("pisa_typical", "trend")]["l1_caught"] == 0


def test_gated_alarms_only_for_open_loop_night_by_default(run):
    observed, _, truth = run
    trace = as_trace(observed, truth, pre_min=QUIET)
    cols = alarm_layer.alarm_columns([trace])[0]
    assert not np.isnan(cols["gate_trial_false_alarms"])
    closed = as_trace(observed, truth, pre_min=QUIET, controller="oref0_tb")
    cols = alarm_layer.alarm_columns([closed])[0]
    assert np.isnan(cols["gate_trial_false_alarms"]) and not np.isnan(cols["trend_alarms"])


@pytest.mark.slow
def test_paper1_simulation_counts_equal_score_alarm():
    """On paper 1's own simulator run (simulate, default scenario, PISA every night)."""
    from ghostlows.sim.engine import HypoTreatment, Therapy, simulate
    from ghostlows.sim.run import task_seeds
    from ghostlows.sim.scenario import DEFAULT, make_scenario

    meal, perturb, noise = task_seeds(1, "child#003")
    sim = simulate(
        "child#003",
        make_scenario(DEFAULT, 2, meal, perturb),
        Therapy(),
        noise,
        treatment=HypoTreatment(),
    )
    samples = sim.samples.assign(seed=1, mix="default")
    observed, faults = apply_faults(samples, FaultConfig(pisa_fraction=1.0), seed=0)
    truth = sim.truth.assign(seed=1, mix="default")
    lows = true_events(truth)
    trace = as_trace(observed.assign(mix="clean"), truth, pre_min=15)
    trace.patient_id = "sim:child#003"
    mine = alarm_layer.score_traces([trace], ["trend", "threshold", "gate_trial"])
    start = DEFAULT_START + pd.Timedelta(minutes=15)
    for prefix in ("trend", "threshold", "gate_trial"):
        alarm = PAPER1[prefix]()
        _, cls, scored = fire_and_score(
            alarm,
            view(observed, "cgm_observed"),
            alarm.default_threshold,
            lows,
            faults,
            observed,
            DEFAULT_WINDOWS,
        )
        after = cls[pd.DatetimeIndex(cls["time"]) >= start]
        got = mine[("clean", prefix)]
        assert got["false_alarms"] == int((~after["matched"].astype(bool)).sum())
        ev = scored[pd.DatetimeIndex(scored["onset"]) >= start]
        assert got["l2_caught"] == int(ev[ev["level"] == 2]["caught"].sum())
