"""The closed-loop engine: paper 1 equivalence, causality, truth invariance, determinism, the
ledger seen through oref0's own IOB, and (slow) ledger consistency over 3 days."""

from __future__ import annotations

import numpy as np
import pytest
from ghostlows.sim.context import iob_fraction
from ghostlows.sim.engine import HypoTreatment, Therapy, simulate
from ghostlows.sim.scenario import DEFAULT, STANDARD, make_scenario

from ghostlows_aid.loop import ClosedLoop, LoopConfig, Person, Routing, make_loop, run_seeds
from ghostlows_aid.oref import REASON_KINDS, Action, iso, node_available, reason_kind
from ghostlows_aid.profile import ProfileSettings

from .conftest import requires_node

PATIENT, SEED = "adult#001", 90001  # a calibration seed (paper 2 calibration seeds are >= 90000)
TIMING = ("ms_iob", "ms_meal", "ms_determine_basal", "ms_bridge")


class Recorder:
    """A stand-in controller: records what it is sent and never acts (no node needed)."""

    name = "recorder"

    def __init__(self):
        self.calls = []

    def step(self, clock_ms, glucose, treatments, currenttemp):
        self.calls.append(
            {
                "clock": clock_ms,
                "glucose": [(g["date"], g["glucose"]) for g in glucose],
                "treatments": [t["created_at"] for t in treatments],
                "currenttemp": dict(currenttemp),
            }
        )
        reply = {"glucose_status": {}, "iob": {"iob": 0.0}, "meal_data": {}, "rT": {}}
        return Action(), reply


class Bump:
    """A test fault: adds `size` mg/dL to the readings whose sample index is in `ks`."""

    def __init__(self, ks, size=60.0):
        self.ks, self.size = set(ks), size

    def apply(self, k, raw):
        return (raw + self.size, "spike", 0) if k in self.ks else (raw, "none", -1)


def recorder_run(faults=None, routing=None, hours=8):
    config = LoopConfig(PATIENT, routing=routing or Routing())
    rec = Recorder()
    loop = make_loop(config, SEED, 1, faults=faults, controller=rec, profile={"stand-in": 1})
    loop.run(hours * 60)
    return loop.result(), rec


# ---- open loop reproduces paper 1 ----
# paper 1's own simulate takes ~9 s a day (pandas parameter lookups), so one mix runs by default
@pytest.mark.parametrize(
    "mix", [DEFAULT, pytest.param(STANDARD, marks=pytest.mark.slow)], ids=lambda m: m.name
)
def test_open_loop_reproduces_paper1_bit_for_bit(mix):
    patient, seed = "adult#004", 3
    seeds = run_seeds(seed, patient)
    scenario = make_scenario(mix, 1, seeds["meal"], seeds["perturb"])
    ref = simulate(patient, scenario, Therapy(), seeds["noise"], treatment=HypoTreatment())
    paper1_person = Person(bolus_rule="paper1", bolus_increment=None, wizard_iob=False)
    config = LoopConfig(patient, controller="open_loop", person=paper1_person, cgm_integer=False)
    r = ClosedLoop(config, scenario, seeds["noise"], seeds["person"]).run().result()
    # paper 1 stores the 1-min truth as float32 and the 5-min samples as float64
    assert np.array_equal(r.truth.astype(np.float32), ref.truth["truth_glucose"].to_numpy())
    assert np.array_equal(r.samples["truth_glucose"], ref.samples["truth_glucose"])
    assert np.array_equal(r.samples["cgm_observed"], ref.samples["glucose"])
    assert np.array_equal(r.samples["cgm_noise"], ref.samples["cgm_noise"])
    for kind in ("bolus", "rescue", "meal_eaten", "carbs_announced"):
        ours = r.events.loc[r.events["kind"] == kind, ["minute", "value"]].to_numpy(float)
        theirs = ref.events.loc[ref.events["kind"] == kind, ["minute", "value"]].to_numpy(float)
        assert np.array_equal(ours, theirs), kind


# ---- causality: the controller sees only readings up to now ----
def test_controller_sees_only_past_readings():
    r, rec = recorder_run()
    rows = r.samples.set_index("minute")
    start_ms = rec.calls[0]["clock"]
    assert len(rec.calls) == len(r.samples)  # one call per 5-min sample
    for call in rec.calls:
        clock = call["clock"]
        dates = [d for d, _ in call["glucose"]]
        assert dates == sorted(dates, reverse=True) and dates[0] == clock  # newest first, fresh
        assert all(t <= iso(clock) for t in call["treatments"])
        for date, value in call["glucose"]:
            assert value == rows.loc[(date - start_ms) // 60000, "cgm_controller"]


def test_a_future_reading_changes_nothing_before_it():
    k = 50
    _, clean = recorder_run(hours=6)
    _, bumped = recorder_run(faults=Bump([k]), hours=6)
    assert clean.calls[:k] == bumped.calls[:k]
    assert clean.calls[k] != bumped.calls[k]
    assert bumped.calls[k]["glucose"][0][1] == clean.calls[k]["glucose"][0][1] + 60


# ---- faults never touch truth ----
def test_faults_reach_only_cgm_observed():
    ks = range(20, 140)
    clean_r, clean_rec = recorder_run()
    routing = Routing(controller_sees="clean", person_sees="clean")
    r, rec = recorder_run(faults=Bump(ks), routing=routing)
    assert np.array_equal(r.truth, clean_r.truth)
    assert np.array_equal(r.insulin, clean_r.insulin)
    for col in ("truth_glucose", "isf_glucose", "cgm_noise", "cgm_clean", "cgm_controller"):
        assert np.array_equal(r.samples[col], clean_r.samples[col]), col
    assert rec.calls == clean_rec.calls
    hit = r.samples["k"].isin(list(ks))
    assert (r.samples.loc[hit, "cgm_observed"] != r.samples.loc[hit, "cgm_clean"]).all()
    assert (r.samples.loc[~hit, "cgm_observed"] == r.samples.loc[~hit, "cgm_clean"]).all()


@requires_node
def test_faults_routed_to_nobody_leave_oref0_and_truth_unchanged():
    def run(faults):
        routing = Routing(controller_sees="clean", person_sees="clean")
        loop = make_loop(LoopConfig(PATIENT, routing=routing), SEED, 1, faults=faults)
        return loop.run(6 * 60).result()

    a, b = run(None), run(Bump(range(10, 70), size=-80.0))
    assert np.array_equal(a.truth, b.truth)
    assert a.treatments == b.treatments
    assert a.samples["rT_reason"].tolist() == b.samples["rT_reason"].tolist()


# ---- determinism ----
@requires_node
def test_closed_loop_is_deterministic():
    def run():
        return make_loop(LoopConfig(PATIENT), SEED, 1).run(6 * 60).result()

    a, b = run(), run()
    assert np.array_equal(a.truth, b.truth)
    assert a.treatments == b.treatments
    pd_a, pd_b = a.samples.drop(columns=list(TIMING)), b.samples.drop(columns=list(TIMING))
    assert pd_a.equals(pd_b)


# ---- the ledger as oref0 sees it ----
@pytest.fixture(scope="module", params=["oref0_tb", "oref1_smb"])
def half_day(request):
    """12 h from midnight (breakfast included) with the rule-1800 ISF (SMBs happen)."""
    if not node_available():
        pytest.skip("node or node_modules/oref0 absent")
    config = LoopConfig(
        PATIENT, controller=request.param, profile=ProfileSettings(isf_rule="rule1800")
    )
    return make_loop(config, SEED, 1).run(12 * 60)


@requires_node
def test_oref0_iob_matches_delivered_insulin(half_day):
    """oref0's net IOB from the ledger equals the IOB of what the pump actually delivered, net of
    the profile basal (paper 1's copy of the oref0 exponential curve). The remaining difference is
    oref0's own 0.05 U discretisation of temps (lib/iob/history.js:551-565)."""
    r = half_day.result()
    net = r.insulin - half_day.profile["current_basal"] / 60
    diffs = []
    for m, iob in r.samples[["minute", "iob"]].dropna().itertuples(index=False):
        m = int(m)
        estimate = float(np.dot(net[:m], iob_fraction(m - np.arange(m))))
        diffs.append(iob - estimate)
    diffs = np.abs(diffs)
    assert len(diffs) > 100
    assert diffs.max() < 0.25
    assert diffs.mean() < 0.08


@requires_node
def test_ledger_consistent_and_every_dose_delivered(half_day):
    r = half_day.result()
    kinds = r.samples["rT_reason"].map(reason_kind)
    assert set(kinds) <= set(REASON_KINDS)
    assert not kinds.isin(["temp_mismatch", "stale"]).any()
    assert r.controller_errors == 0
    # every minute delivered the ledger's temp (or the scheduled basal) plus its boluses and SMBs
    led = half_day.state.ledger
    assert led.temps(), "oref0 set no temp in 12 h"
    expected = np.full(r.insulin.size, half_day.basal)
    for start, rate, duration in led.temps():
        expected[start : start + duration] = rate / 60
    for minute, units in led.boluses():
        expected[minute] += units
    assert np.allclose(r.insulin, expected, rtol=0, atol=1e-12)
    if half_day.config.controller == "oref1_smb":
        smb = r.events[r.events["kind"] == "smb"]
        assert len(smb) > 0
        assert all(round(u * 10, 9) == int(round(u * 10)) for u in smb["value"])  # 0.1 U steps


# ---- slow: 3 days ----
@pytest.mark.slow
@requires_node
@pytest.mark.parametrize("controller", ["oref0_tb", "oref1_smb"])
@pytest.mark.parametrize("isf_rule", ["quest", "rule1800"])
def test_three_days_without_temp_mismatch_or_stale(controller, isf_rule):
    config = LoopConfig(PATIENT, controller=controller, profile=ProfileSettings(isf_rule=isf_rule))
    r = make_loop(config, SEED, 3).run().result()
    kinds = r.samples["rT_reason"].map(reason_kind).value_counts()
    assert kinds.get("temp_mismatch", 0) == 0 and kinds.get("stale", 0) == 0, kinds.to_dict()
    assert r.controller_errors == 0
    assert len(r.samples) == 3 * 288
