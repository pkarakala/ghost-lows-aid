"""Online faults: paper 1's post hoc injection reproduced sample by sample, causally, and never
touching truth; single events and the new bias operators. CGM traces only (no outcomes)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from ghostlows.faults import pisa as pisa_model
from ghostlows.faults.suite import FaultConfig, apply_faults, apply_run

from ghostlows_aid.faults_online import BIAS_SHAPE, SINGLE_KINDS, FaultSchedule, SingleFault
from ghostlows_aid.loop import LoopConfig, Person, Routing, make_loop

DAYS = 7
LIMITS = (40.0, 400.0)
T0 = pd.Timestamp("2026-01-05")
CONFIGS = {
    "pisa": FaultConfig(enabled=frozenset({"pisa"}), pisa_fraction=1.0),
    "drift": FaultConfig(enabled=frozenset({"drift"}), severity=1.0, sensor_days=2.0),
    "spike": FaultConfig(enabled=frozenset({"spike"}), severity=1.0, spike_per_day=6.0),
    "stale": FaultConfig(enabled=frozenset({"stale"}), severity=1.0, stale_per_day=6.0),
    "dropout": FaultConfig(enabled=frozenset({"dropout"}), severity=1.0, dropout_per_day=8.0),
    "mixed": FaultConfig(pisa_fraction=0.5),  # paper 1's suite at moderate severity
    "mixed_dense": FaultConfig(
        pisa_fraction=1.0,
        severity=1.0,
        sensor_days=3.0,
        dropout_per_day=8.0,
        stale_per_day=6.0,
        spike_per_day=6.0,
    ),
}
SEEDS = [0, 1, 7, 90001]


def clean_trace(days: int = DAYS, seed: int = 0) -> pd.DataFrame:
    """A clean CGM trace in paper 1's columns; it crosses both ends of the device range."""
    rng = np.random.default_rng(seed)
    m = np.arange(days * 288) * 5.0
    isf = 210 + 190 * np.sin(m / 300.0)
    noise = rng.normal(0, 6, m.size)
    return pd.DataFrame(
        {
            "time": T0 + pd.to_timedelta(m, unit="min"),
            "isf_glucose": isf,
            "cgm_noise": noise,
            "glucose": np.clip(isf + noise, *LIMITS),
        }
    )


def online(trace: pd.DataFrame, faults, start: int = 0):
    raw = (trace["isf_glucose"] + trace["cgm_noise"]).to_numpy(dtype=float)
    out = [faults.apply(k, float(raw[k])) for k in range(start, raw.size)]
    obs = np.clip(np.array([o[0] for o in out]), *LIMITS)
    return obs, np.array([o[1] for o in out]), np.array([o[2] for o in out])


def assert_same(obs, label, fid, ref: pd.DataFrame):
    want = ref["cgm_observed"].to_numpy(dtype=float)
    assert np.array_equal(np.isnan(obs), np.isnan(want))
    ok = ~np.isnan(want)
    np.testing.assert_allclose(obs[ok], want[ok], rtol=0, atol=1e-9)
    assert label.tolist() == ref["fault_type"].tolist()
    assert fid.tolist() == ref["fault_id"].tolist()


# ---- equivalence with paper 1 ----
@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize("name", list(CONFIGS))
def test_schedule_reproduces_apply_run(name, seed):
    config = CONFIGS[name]
    trace = clean_trace()
    ref, instances = apply_run(trace, config, seed)
    schedule = FaultSchedule.draw(seed, DAYS, config)
    obs, label, fid = online(trace, schedule)
    assert_same(obs, label, fid, ref)
    assert (label != "none").any(), "the configuration injected nothing"
    ours = pd.DataFrame(schedule.instances)
    assert len(ours) == len(instances)
    for col in ("start", "end"):
        theirs = (instances[col] - T0) / pd.Timedelta("1min")
        np.testing.assert_allclose(ours[col], theirs, rtol=0, atol=1e-6)
    assert ours["fault_type"].tolist() == instances["fault_type"].tolist()
    assert ours["fault_id"].tolist() == instances["fault_id"].tolist()
    assert ("size" in ours) == ("size" in instances)  # stale and dropout have no size
    if "size" in ours:
        np.testing.assert_allclose(ours["size"], instances["size"], rtol=0, atol=1e-12)


def test_schedule_reproduces_apply_faults_per_run():
    config = CONFIGS["mixed_dense"]
    runs = []
    for patient, seed, mix in [("sim:adult#001", 3, "default"), ("sim:child#004", 9, "high")]:
        runs.append(clean_trace(days=3, seed=seed).assign(patient_id=patient, seed=seed, mix=mix))
    ref, _ = apply_faults(pd.concat(runs, ignore_index=True), config, seed=42)
    for run in runs:
        key = (run["patient_id"].iloc[0], run["seed"].iloc[0], run["mix"].iloc[0])
        mask = (ref["patient_id"] == key[0]) & (ref["seed"] == key[1])
        obs, label, fid = online(run, FaultSchedule.for_run(42, key, 3, config))
        assert_same(obs, label, fid, ref[mask].reset_index(drop=True))


def test_schedule_is_causal():
    """Changing readings after sample j changes no output up to j."""
    config, j = CONFIGS["mixed_dense"], 700
    trace = clean_trace()
    changed = trace.copy()
    changed.loc[j + 1 :, "isf_glucose"] += 37.0
    a = online(trace, FaultSchedule.draw(5, DAYS, config))
    b = online(changed, FaultSchedule.draw(5, DAYS, config))
    for x, y in zip(a, b, strict=True):
        assert np.array_equal(x[: j + 1], y[: j + 1], equal_nan=x.dtype.kind == "f")


def test_schedule_order_and_late_start():
    config = CONFIGS["mixed_dense"]
    trace = clean_trace()
    full = online(trace, FaultSchedule.draw(5, DAYS, config))
    schedule = FaultSchedule.draw(5, DAYS, config)
    stale_inside = {k for first, last, _ in schedule.stale for k in range(first + 1, last + 1)}
    start = next(k for k in range(600, 2000) if k not in stale_inside)
    late = online(trace, schedule, start=start)  # a branch that switches the schedule on late
    for x, y in zip(full, late, strict=True):
        assert np.array_equal(x[start:], y, equal_nan=x.dtype.kind == "f")
    with pytest.raises(ValueError, match="increasing"):
        schedule.apply(start, 100.0)
    if stale_inside:
        fresh = FaultSchedule.draw(5, DAYS, config)
        with pytest.raises(ValueError, match="stale"):
            fresh.apply(min(stale_inside), 100.0)


def test_online_in_the_loop_equals_post_hoc():
    """The loop's own fault hook, on the loop's own clean readings, gives paper 1's result.
    Nobody is shown the faulted signal, so the run is also identical to the fault-free one."""
    config = FaultConfig(pisa_fraction=1.0, severity=1.0)
    loop_config = LoopConfig(
        "adult#001",
        controller="open_loop",
        person=Person(bolus_rule="paper1", bolus_increment=None),
        routing=Routing(controller_sees="clean", person_sees="clean"),
        cgm_integer=False,
    )
    schedule = FaultSchedule.draw(11, 2, config)
    r = make_loop(loop_config, 90002, 2, faults=schedule).run().result()
    clean = make_loop(loop_config, 90002, 2).run().result()
    ref, _ = apply_run(r.samples[["time", "isf_glucose", "cgm_noise"]], config, 11)
    s = r.samples
    assert_same(
        s["cgm_observed"].to_numpy(dtype=float),
        s["fault_type"].to_numpy(),
        s["fault_id"].to_numpy(),
        ref,
    )
    assert (s["fault_type"] != "none").any()
    assert np.array_equal(r.truth, clean.truth)
    assert np.array_equal(s["cgm_clean"], clean.samples["cgm_clean"])


# ---- single events ----
ONSET = 600.0  # minutes from the run start (sample 120)


def run_single(fault: SingleFault, days: int = 1):
    trace = clean_trace(days=days)
    raw = (trace["isf_glucose"] + trace["cgm_noise"]).to_numpy(dtype=float)
    out = [fault.apply(k, float(raw[k])) for k in range(raw.size)]
    m = np.arange(raw.size) * 5.0
    obs = np.array([o[0] for o in out])
    return m, raw, obs, np.array([o[1] for o in out])


@pytest.mark.parametrize("kind", SINGLE_KINDS)
def test_single_event_changes_nothing_before_onset(kind):
    m, raw, obs, label = run_single(SingleFault.draw(kind, ONSET, seed=3))
    before = m < ONSET
    assert np.array_equal(obs[before], raw[before])
    assert (label[before] == "none").all()
    assert (obs[~before] != raw[~before]).any()  # (NaN != x for a dropout)


def test_single_pisa_is_paper1_pisa():
    fault = SingleFault.at_quantile("pisa", ONSET, 0.95)
    depth, duration = fault.params["depth"], fault.params["duration"]
    assert depth == pytest.approx(48.0) and duration == pytest.approx(70.0)  # paper 1's 95th pct
    m, raw, obs, label = run_single(fault)
    pisa = pisa_model.from_size(ONSET, depth, duration)
    np.testing.assert_allclose(raw - obs, pisa.attenuation(m), rtol=0, atol=1e-9)
    assert np.array_equal(label == "pisa", pisa.active(m))
    assert (raw - obs).max() == pytest.approx(48.0, abs=0.5)  # sampled every 5 min


def test_single_dropout_stale_spike_windows():
    m, raw, obs, label = run_single(SingleFault("dropout", ONSET, {"duration": 22.0}))
    gap = (m >= ONSET) & (m < ONSET + 22)
    assert np.isnan(obs[gap]).all() and np.array_equal(obs[~gap], raw[~gap])
    assert gap.sum() == 5
    m, raw, obs, label = run_single(SingleFault("stale", ONSET, {"duration": 30.0}))
    seg = np.flatnonzero((m >= ONSET) & (m < ONSET + 30))
    assert (obs[seg] == raw[seg[0]]).all() and label[seg[0]] == "none"
    assert (label[seg[1:]] == "stale").all() and obs[seg[-1] + 1] == raw[seg[-1] + 1]
    spike = SingleFault.draw("spike", ONSET, seed=4)
    m, raw, obs, label = run_single(spike)
    hit = (m >= ONSET) & (m < ONSET + spike.params["duration"])
    assert (obs[hit] != raw[hit]).all() and np.array_equal(obs[~hit], raw[~hit])
    assert np.array_equal(run_single(spike)[2], obs)  # the noise is fixed by the event


@pytest.mark.parametrize("kind", ["bias_pos", "bias_neg"])
def test_bias_shape(kind):
    sign = 1 if kind == "bias_pos" else -1
    fault = SingleFault.at_quantile(kind, ONSET, 0.5)
    size = fault.params["size_mg_dl"]
    assert size == sign * 50.0
    rise, plateau, fall = (BIAS_SHAPE[k] for k in ("rise_min", "plateau_min", "fall_min"))
    assert fault.bias_mg_dl(-5) == 0 and fault.bias_mg_dl(0) == 0
    assert fault.bias_mg_dl(rise / 2) == pytest.approx(size / 2)
    assert fault.bias_mg_dl(rise) == fault.bias_mg_dl(rise + plateau - 1) == size
    assert fault.bias_mg_dl(rise + plateau + fall / 2) == pytest.approx(size / 2)
    assert fault.bias_mg_dl(rise + plateau + fall) == 0
    m, raw, obs, label = run_single(fault)
    t = m - ONSET
    expect = np.array([fault.bias_mg_dl(x) for x in t])
    np.testing.assert_allclose(obs - raw, expect, rtol=0, atol=1e-12)
    assert np.array_equal(label == kind, expect != 0)
    assert fault.end == ONSET + rise + plateau + fall


def test_single_drift_is_multiplicative_and_capped():
    fault = SingleFault("drift", ONSET, {"size": -0.2, "age_days": 9.9, "sensor_days": 10.0})
    m, raw, obs, label = run_single(fault)
    t = m - ONSET
    age = np.minimum(9.9 + t / 1440, 10.0)
    after = t >= 0
    np.testing.assert_allclose(obs[after] / raw[after], 1 - 0.2 * age[after] / 10, rtol=1e-12)
    assert (label[after] == "drift").all()


def test_single_event_draws_follow_paper1_distributions():
    draws = [SingleFault.draw("pisa", 0.0, seed=s).params for s in range(400)]
    depth = np.median([d["depth"] for d in draws])
    duration = np.median([d["duration"] for d in draws])
    assert 21 < depth < 27 and 42 < duration < 49  # medians 24 mg/dL and 45 min
    median = SingleFault.at_quantile("pisa", 0.0, 0.5).params
    assert median["depth"] == pytest.approx(24.0) and median["duration"] == pytest.approx(45.0)
    gap = SingleFault.at_quantile("dropout", 0.0, 0.9).params["duration"]
    assert gap == pytest.approx(20.0)  # 90% of gaps under 20 min (paper 1's rate)
    assert SingleFault.draw("spike", 0.0, seed=1).params["sd"] == 20.0 * 0.5
    explicit = SingleFault.draw("pisa", 0.0, seed=1, depth=30.0, duration=60.0).params
    assert explicit["depth"] == 30.0 and explicit["duration"] == 60.0
