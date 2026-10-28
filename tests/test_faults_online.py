"""Online faults: paper 1's post hoc injection reproduced sample by sample, causally, and never
touching truth; single events and the new bias operators. CGM traces only (no outcomes)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from ghostlows.faults import pisa as pisa_model
from ghostlows.faults.suite import FaultConfig, apply_faults, apply_run

from ghostlows_aid.faults_online import (
    INTERFERENCE,
    SINGLE_KINDS,
    FaultSchedule,
    SingleFault,
    g6_gap_samples,
)
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


def example(kind: str) -> SingleFault:
    """One event of each kind: drawn where paper 1 has a distribution, else a memo level."""
    if kind == "gain":
        return SingleFault.draw("gain", ONSET, seed=3, gain=-0.4)
    if kind == "bias":
        return SingleFault.interference(ONSET, "large")
    return SingleFault.draw(kind, ONSET, seed=3)


@pytest.mark.parametrize("kind", SINGLE_KINDS)
def test_single_event_changes_nothing_before_onset(kind):
    m, raw, obs, label = run_single(example(kind))
    before = m < ONSET
    assert np.array_equal(obs[before], raw[before])
    assert (label[before] == "none").all()
    assert (obs[~before] != raw[~before]).any()  # (NaN != x for a dropout)


def test_single_pisa_is_paper1_pisa():
    fault = SingleFault.at_quantile("pisa", ONSET, 0.95)
    depth, duration = fault.params["depth"], fault.params["duration"]
    assert depth == pytest.approx(48.0) and duration == pytest.approx(70.0)  # PISA severe
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


def test_g6_gap_model():
    """Vettoretti 2023 Model 4 as PI_decisions §6 #4 wires it: median 20 min, 73% of gaps
    longer than oref0's 12-min staleness limit (3 samples or more)."""
    assert g6_gap_samples(q=0.5) == 4
    rng = np.random.default_rng(0)
    n = np.array([g6_gap_samples(rng) for _ in range(20000)])
    assert np.median(n) == 4
    assert np.mean(n >= 3) == pytest.approx(0.855**2, abs=0.01)
    assert np.mean(n >= 2) == pytest.approx(0.855, abs=0.01)
    fault = SingleFault.draw("dropout", ONSET, seed=5, gap_model="g6")
    assert fault.params["model"] == "g6" and fault.params["duration"] % 5 == 0
    assert SingleFault.at_quantile("dropout", 0.0, 0.5, gap_model="g6").params["duration"] == 20


@pytest.mark.parametrize("level", ["small", "large"])
def test_interference_shape(level):
    fault = SingleFault.interference(ONSET, level)
    p = INTERFERENCE[level]
    size, tp, hold, tau = p["size_mg_dl"], p["tp_min"], p["hold_until_min"], p["tau_min"]
    assert fault.bias_mg_dl(-5) == 0 and fault.bias_mg_dl(0) == 0
    assert fault.bias_mg_dl(tp / 2) == pytest.approx(size * 0.5 * np.exp(0.5))
    assert fault.bias_mg_dl(tp) == fault.bias_mg_dl(hold) == size
    assert fault.bias_mg_dl(hold + tau) == pytest.approx(size / np.e)
    assert fault.end == pytest.approx(ONSET + hold + tau * np.log(10))
    m, raw, obs, label = run_single(fault, days=2)
    expect = np.array([fault.bias_mg_dl(t) for t in m - ONSET])
    np.testing.assert_allclose(obs - raw, expect, rtol=0, atol=1e-12)
    assert np.array_equal(label == "bias", expect >= 0.1 * size)
    if level == "large":  # about 20% left at 8 h (PI_decisions §6 #17)
        assert fault.bias_mg_dl(480) / size == pytest.approx(np.exp(-1.5))


def test_negative_bias_is_a_sustained_false_low():
    fault = SingleFault("bias", ONSET, {**INTERFERENCE["large"], "size_mg_dl": -40.0})
    m, raw, obs, label = run_single(fault)
    after = m > ONSET
    assert (obs[after] < raw[after]).all() and (obs[m <= ONSET] == raw[m <= ONSET]).all()


def test_gain_and_drift_are_held_multiplicative():
    for kind, gain in [("gain", 0.4), ("gain", -0.2), ("drift", -0.2)]:
        fault = SingleFault(kind, ONSET, {"gain": gain})
        m, raw, obs, label = run_single(fault)
        after = m >= ONSET
        np.testing.assert_allclose(obs[after], raw[after] * (1 + gain), rtol=1e-15)
        assert (label[after] == kind).all()
    low = SingleFault.draw("drift", ONSET, seed=8, sign=-1.0)
    high = SingleFault.draw("drift", ONSET, seed=8, sign=1.0)
    assert low.params["gain"] == -high.params["gain"] < 0
    assert abs(low.params["gain"]) <= 0.2 * 0.5  # paper 1's maximum at s = 0.5
    with pytest.raises(ValueError, match="explicitly"):
        SingleFault.draw("gain", ONSET, seed=1)


def test_single_event_draws_follow_paper1_distributions():
    draws = [SingleFault.draw("pisa", 0.0, seed=s).params for s in range(400)]
    depth = np.median([d["depth"] for d in draws])
    duration = np.median([d["duration"] for d in draws])
    assert 21 < depth < 27 and 42 < duration < 49  # medians 24 mg/dL and 45 min
    median = SingleFault.at_quantile("pisa", 0.0, 0.5).params
    assert median["depth"] == pytest.approx(24.0) and median["duration"] == pytest.approx(45.0)
    gap = SingleFault.at_quantile("dropout", 0.0, 0.9).params["duration"]
    assert gap == pytest.approx(20.0)  # 90% of gaps under 20 min (paper 1's rate)
    spike = SingleFault.draw("spike", 0.0, seed=1).params
    assert spike["sd"] == 20.0 * 0.5 and 15 <= spike["duration"] <= 30  # typical: s = 0.5
    explicit = SingleFault.draw("pisa", 0.0, seed=1, depth=30.0, duration=60.0).params
    assert explicit["depth"] == 30.0 and explicit["duration"] == 60.0
    assert SingleFault.draw("pisa", 0.0, seed=1).describe()["fault_kind"] == "pisa"
