"""Branching counterfactual episodes: the snapshot holds the whole state, a branched clean arm is
the continuous run, branches are reproducible, arms share their history and their random
streams. Faulted arms are compared on CGM streams, meals and controller actions only."""

from __future__ import annotations

import numpy as np
import pytest
from ghostlows.sim.scenario import DEFAULT

from ghostlows_aid import episodes
from ghostlows_aid.detectors import make_detector
from ghostlows_aid.episodes import (
    ARMS,
    FAULTS,
    SLOTS,
    arms_for,
    event_seed,
    onset,
    pre_onset_equal,
    run_unit,
)
from ghostlows_aid.faults_online import SingleFault
from ghostlows_aid.loop import CONTROLLERS, CORRECTING, LoopConfig, Routing, make_loop
from ghostlows_aid.oref import node_available

PATIENT, SEED = "adult#001", 90001
T0, END = 300, 600  # 05:00 to 10:00: breakfast, its bolus and oref0's response fall inside
TIMING = ["ms_iob", "ms_meal", "ms_determine_basal", "ms_bridge"]


def needs_node(controller):
    if controller != "open_loop" and not node_available():
        pytest.skip("node or node_modules/oref0 absent")


def samples(r):
    return r.samples.drop(columns=[c for c in TIMING if c in r.samples])


def assert_same_run(a, b):
    assert np.array_equal(a.truth, b.truth)
    for name in ("insulin", "insulin_basal", "insulin_bolus", "insulin_smb"):
        assert np.array_equal(getattr(a, name), getattr(b, name)), name
    assert samples(a).equals(samples(b))
    assert a.events.equals(b.events)
    assert a.treatments == b.treatments


def actions(r, t0):
    """What the controller and the person did after t0: CGM streams, oref0's outputs, the
    ledger. No glucose outcome."""
    s = samples(r)
    s = s[s["minute"] >= t0]
    cols = [c for c in ("cgm_observed", "cgm_controller", "cgm_sent", "rT_reason", "rT_rate",
                        "rT_duration", "rT_units", "carbs_rescue") if c in s]  # fmt: skip
    return s[cols].reset_index(drop=True), r.treatments


def trunk_with_snapshot(controller, t0=T0, end=END):
    loop = make_loop(LoopConfig(PATIENT, controller=controller), SEED, 1, mix=DEFAULT)
    loop.run(t0)
    snap = loop.snapshot()
    loop.run(end)
    return loop, snap, loop.result()


# ---- (a) a branched clean arm is the continuous run ----
@pytest.mark.parametrize("controller", CONTROLLERS)
def test_branched_clean_arm_equals_continuous_run(controller):
    needs_node(controller)
    loop, snap, trunk = trunk_with_snapshot(controller)
    loop.branch(snap).run(END)
    assert_same_run(loop.result(), trunk)
    # into a fresh ClosedLoop as well: nothing outside the snapshot carries state
    fresh = make_loop(LoopConfig(PATIENT, controller=controller), SEED, 1, mix=DEFAULT)
    fresh.branch(snap).run(END)
    assert_same_run(fresh.result(), trunk)


def test_integrator_is_rebuilt_bit_for_bit():
    """The scipy integrator is not copied but rebuilt from (y, t), mid-meal."""
    loop, snap, trunk = trunk_with_snapshot("open_loop", t0=457, end=700)
    assert snap.patient.t == 457 and "_odesolver" not in snap.patient.attrs
    continuous = loop.patient.state.copy()
    loop.branch(snap).run(700)
    assert np.array_equal(loop.patient.state, continuous)


# ---- (b) restoring twice gives identical arms ----
@pytest.mark.parametrize("controller", ["oref0_tb", "oref1_smb"])
def test_restoring_twice_gives_identical_arms(controller):
    needs_node(controller)
    loop, snap, _ = trunk_with_snapshot(controller)
    runs = []
    for _ in range(2):  # a variant arm without a fault: everything compared
        loop.branch(
            snap, routing=Routing(gate=True), person=CORRECTING, gate=make_detector("baysal")
        )
        runs.append(loop.run(END).result())
    assert_same_run(*runs)
    faulted = []
    for _ in range(2):  # a faulted arm: CGM streams and actions only
        fault = SingleFault.draw("pisa", T0 + 30, event_seed(SEED, PATIENT, 1))
        loop.branch(snap, faults=fault)
        faulted.append(actions(loop.run(END).result(), T0))
    assert faulted[0][0].equals(faulted[1][0]) and faulted[0][1] == faulted[1][1]


# ---- (c) every arm is identical before t0 ----
def test_arms_identical_before_t0():
    needs_node("oref0_tb")
    loop, snap, trunk = trunk_with_snapshot("oref0_tb")
    variants = [
        dict(faults=SingleFault("dropout", T0, {"duration": 120.0})),
        dict(faults=SingleFault.interference(T0, "large"), routing=Routing(person_sees="clean")),
        dict(routing=Routing(gate=True), gate=make_detector("matched_filter")),
        dict(person=CORRECTING),
    ]
    for kw in variants:
        r = loop.branch(snap, **kw).run(T0 + 60).result()
        assert pre_onset_equal(r, trunk, T0)


# ---- (d) after t0 the arms share the sensor-noise stream and the meal log ----
def test_arms_share_noise_and_meals_after_t0():
    needs_node("oref0_tb")
    loop, snap, trunk = trunk_with_snapshot("oref0_tb")
    meal_kinds = ["meal_eaten", "carbs_announced"]

    def meals(r):
        ev = r.events
        return ev[ev["kind"].isin(meal_kinds) & (ev["minute"] >= T0)].reset_index(drop=True)

    fault = SingleFault.at_quantile("pisa", T0 + 20, 0.95)
    for kw in (dict(faults=fault), dict(faults=SingleFault.interference(T0, "large"))):
        r = loop.branch(snap, **kw).run(END).result()
        after = r.samples["minute"] >= T0
        assert np.array_equal(
            r.samples.loc[after, "cgm_noise"], trunk.samples.loc[after, "cgm_noise"]
        )
        assert meals(r).equals(meals(trunk))
        # cgm_clean is isf + noise of the arm's own body: equal until what the body receives
        # (insulin, rescue carbohydrate) first differs
        diff = np.flatnonzero(r.insulin != trunk.insulin)
        rescue = r.samples["carbs_rescue"].to_numpy() != trunk.samples["carbs_rescue"].to_numpy()
        firsts = [END] + ([diff[0]] if diff.size else [])
        firsts += list(r.samples.loc[rescue, "minute"])[:1]
        first = min(firsts)
        upto = (r.samples["minute"] >= T0) & (r.samples["minute"] <= first)
        assert np.array_equal(
            r.samples.loc[upto, "cgm_clean"], trunk.samples.loc[upto, "cgm_clean"]
        )


# ---- the design: slots, stratified onsets, the arm catalogue ----
def test_onsets_are_stratified_on_5_min_boundaries():
    n = 20
    for slot, (lo, hi) in SLOTS.items():
        t0s = [onset(slot, 101 + k, PATIENT, k, n) for k in range(n)]
        width = (hi - lo) / n
        for k, t0 in enumerate(t0s):
            assert t0 % 5 == 0 and lo + k * width <= t0 < lo + (k + 1) * width
        assert onset(slot, 101, PATIENT, 0, n) == t0s[0]  # no controller in the key
    assert episodes.TRUNK_MIN >= SLOTS["night"][1] - 5 + episodes.WINDOW_MIN


def test_arm_catalogue_matches_the_design():
    assert [f.number for f in FAULTS] == list(range(1, 18))
    assert len({a.name for a in ARMS}) == len(ARMS)
    assert len(arms_for("oref0_tb")) == 34
    assert len(arms_for("oref1_smb")) == len(arms_for("open_loop")) == 17
    for controller in CONTROLLERS:
        assert len([a for a in arms_for(controller) if a.tier == 1]) == 7
    groups = [a.group for a in arms_for("oref0_tb")]
    assert {g: groups.count(g) for g in set(groups)} == {
        "fault": 17, "pathway": 4, "gate": 8, "detector2": 2, "guard": 1, "correcting": 2
    }  # fmt: skip


def test_run_unit_shares_onsets_and_events_across_controllers(monkeypatch):
    """A short design (monkeypatched slots) through run_unit, open loop, faults switched off."""
    monkeypatch.setattr(episodes, "SLOTS", {"night": (300, 330), "dose": (200, 230),
                                             "any": (100, 130)})  # fmt: skip
    monkeypatch.setattr(episodes, "TRUNK_MIN", 600)
    monkeypatch.setattr(episodes, "WINDOW_MIN", 120)
    chosen = [episodes.ARM_BY_NAME[n] for n in ("pisa_typical", "frozen_severe")]
    unit, runs = run_unit("open_loop", PATIENT, SEED, 0, 1, arms=chosen, with_faults=False)
    runs = list(runs)
    assert [r.arm.name for r in runs] == [
        "clean@night",
        "clean@dose",
        "clean@any",
        *[a.name for a in chosen],
    ]
    for r in runs[3:]:
        assert r.fault is None and pre_onset_equal(r.result, unit.trunk, r.t0)
        # no fault: the branch is the trunk's continuation
        end = r.t0 + 120
        assert np.array_equal(r.result.truth[:end], unit.trunk.truth[:end])
    trace = episodes.arm_trace(unit, runs[3], window=120)
    assert trace["truth_1min"].size == 120 and len(trace["samples"]) == 24
    assert set(trace["onset_state"]) == {"iob", "cob"}
