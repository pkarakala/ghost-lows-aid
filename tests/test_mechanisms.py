"""Mechanism checks and the harness-validity suite (PI_decisions rev 2 §11, N12): oref0's shipped
guards act in this harness, and every routing hook changes only what it should. Assertions are on
controller inputs, reasons and actions only; no glucose outcome is computed."""

from __future__ import annotations

import numpy as np
import pytest
from ghostlows.sim.engine import PATIENTS
from ghostlows.sim.scenario import DEFAULT

from ghostlows_aid.detectors import make_detector
from ghostlows_aid.faults_online import SingleFault
from ghostlows_aid.ledger import round_basal
from ghostlows_aid.loop import CORRECTING, LoopConfig, Routing, make_loop
from ghostlows_aid.oref import node_available, reason_kind
from ghostlows_aid.profile import ProfileSettings, build_profile, therapy_numbers

PATIENT, SEED = "adult#001", 90001
SETTINGS = ProfileSettings(isf_rule="rule1800")
NIGHT, DAWN, MORNING = 120, 360, 420  # 02:00, 06:00, 07:00 on day 1

pytestmark = pytest.mark.skipif(not node_available(), reason="node or node_modules/oref0 absent")


@pytest.fixture(scope="module")
def trunk():
    """One oref0_tb trunk with snapshots at 02:00, 06:00 and 07:00."""
    loop = make_loop(LoopConfig(PATIENT, profile=SETTINGS), SEED, 1, mix=DEFAULT)
    snaps = {}
    for t0 in (NIGHT, DAWN, MORNING):
        loop.run(t0)
        snaps[t0] = loop.snapshot()
    return loop, snaps


def arm(trunk, t0, hours, **kw):
    loop, snaps = trunk
    loop.branch(snaps[t0], **kw).run(t0 + int(hours * 60))
    s = loop.result().samples
    return loop, s[s["minute"] >= t0].reset_index(drop=True)


def kinds(s):
    return s["rT_reason"].map(reason_kind)


# ---- the stale-data stop and its ablation ----
def test_stale_data_stop_acts_during_a_long_gap(trunk):
    _, s = arm(trunk, DAWN, 3, faults=SingleFault("dropout", DAWN, {"duration": 120.0}))
    t = s["minute"] - DAWN
    k = kinds(s)
    # the last reading is at t0 - 5: oref0 stops once it is more than 12 min old
    assert (k[(t >= 10) & (t < 120)] == "stale").all()
    assert not (k[(t < 10) | (t >= 120)] == "stale").any()
    basal = round_basal(trunk[0].profile["current_basal"])
    gap = (t >= 10) & (t < 120)
    assert not (s.loc[gap, "temp_set"].fillna(0) > basal).any()  # no high temp while blind


def test_restamping_bypasses_the_stale_stop(trunk):
    fault = SingleFault("dropout", DAWN, {"duration": 120.0})
    _, s = arm(trunk, DAWN, 3, faults=fault, routing=Routing(restamp_missing=True))
    t = s["minute"] - DAWN
    assert not (kinds(s) == "stale").any()
    assert s.loc[t < 120, "restamped"].all() and not s.loc[t >= 120, "restamped"].any()
    assert (kinds(s)[(t >= 50) & (t < 120)] == "flat").all()  # the re-sent value is flat


# ---- the flat-CGM check and the fakecgm ablation ----
def test_flat_check_acts_on_a_long_freeze_unless_device_is_fakecgm(trunk):
    fault = SingleFault("stale", DAWN, {"duration": 120.0})
    _, s = arm(trunk, DAWN, 3, faults=fault)
    t = s["minute"] - DAWN
    assert (kinds(s)[(t >= 50) & (t < 120)] == "flat").all()
    loop, s = arm(trunk, DAWN, 3, faults=fault, routing=Routing(fakecgm=True))
    assert not (kinds(s) == "flat").any()
    assert {e["device"] for e in loop.state.hist_entries if e["date"] >= _ms(loop, DAWN)} == {
        "fakecgm"
    }
    assert {e["device"] for e in loop.state.hist_entries if e["date"] < _ms(loop, DAWN)} == {
        "ghostlows-aid-cgm"
    }


def _ms(loop, minute):
    return loop.start_ms + minute * 60000


# ---- low-glucose suspend ----
def test_zero_temp_follows_a_predicted_low(trunk):
    fault = SingleFault("pisa", NIGHT, {"depth": 48.0, "duration": 70.0})
    _, s = arm(trunk, NIGHT, 2, faults=fault)
    during = s["minute"] < NIGHT + 70
    suspend = during & (kinds(s) == "low_suspend")
    assert (s.loc[suspend, "rT_rate"] == 0).any()
    # a suspend either sets a zero temp or keeps the zero temp already running
    # (lib/basal-set-temp.js: "no temp required")
    assert (s.loc[suspend, "rT_rate"].fillna(0) == 0).all()


# ---- max_iob ----
def test_max_iob_is_set_for_every_patient():
    for name in PATIENTS:
        n = therapy_numbers(name, SETTINGS)
        assert n["max_iob"] > 0  # oref0's default of 0 would block every high temp
    profile = build_profile(PATIENT, SETTINGS)
    assert profile["max_iob"] == pytest.approx(therapy_numbers(PATIENT, SETTINGS)["max_iob"])


def test_max_iob_binds():
    """A deliberately low max_iob (a test setting) under a false high: oref0 caps the temp."""
    settings = ProfileSettings(isf_rule="rule1800", max_iob=1.0)
    loop = make_loop(LoopConfig(PATIENT, profile=settings), SEED, 1, mix=DEFAULT)
    loop.run(MORNING)
    loop.branch(loop.snapshot(), faults=SingleFault("gain", MORNING, {"gain": 0.4}))
    s = loop.run(MORNING + 180).result().samples
    s = s[s["minute"] >= MORNING]
    capped = s["rT_reason"].str.contains("max_iob", na=False)
    assert capped.any()
    over = s["rT_reason"].str.contains("> max_iob", na=False)
    basal = round_basal(loop.profile["current_basal"])
    assert (s.loc[over, "rT_rate"].dropna() <= basal).all()  # IOB above max_iob: no high temp


# ---- SMB ----
def test_smb_disabled_by_max_delta_on_a_jump():
    loop = make_loop(LoopConfig(PATIENT, controller="oref1_smb", profile=SETTINGS), SEED, 1,
                     mix=DEFAULT)  # fmt: skip
    loop.run(MORNING)
    jump = SingleFault(
        "bias",
        MORNING,
        {"size_mg_dl": 40.0, "tp_min": 5.0, "hold_until_min": 60.0, "tau_min": 60.0},
    )
    loop.branch(loop.snapshot(), faults=jump)
    s = loop.run(MORNING + 60).result().samples
    s = s[s["minute"] >= MORNING]
    blocked = s["rT_reason"].str.contains("SMB disabled", na=False)
    assert blocked.any()
    assert (s.loc[blocked, "rT_units"].fillna(0) == 0).all()


# ---- the detector gate ----
@pytest.mark.parametrize("policy", ["hold", "drop"])
def test_gate_policies(trunk, policy):
    fault = SingleFault("pisa", NIGHT, {"depth": 48.0, "duration": 70.0})
    routing = Routing(gate=True, gate_policy=policy)
    _, s = arm(trunk, NIGHT, 3, faults=fault, routing=routing, gate=make_detector("baysal"))
    flagged = s["gate_flag"].astype(bool)
    assert flagged.sum() >= 3
    passed = s.loc[~flagged]
    assert passed["cgm_sent"].equals(passed["cgm_controller"])
    if policy == "drop":
        assert s.loc[flagged, "cgm_sent"].isna().all()
    else:
        for i in np.flatnonzero(flagged):  # the last unflagged reading, re-sent as fresh
            last = s.loc[: i - 1].loc[~flagged.loc[: i - 1], "cgm_controller"].iloc[-1]
            assert s.loc[i, "cgm_sent"] == last
    assert s["cgm_person"].equals(s["cgm_observed"])  # the person sees the raw displayed CGM


def test_gate_detector_only_observes_without_the_routing_flag(trunk):
    fault = SingleFault("pisa", NIGHT, {"depth": 48.0, "duration": 70.0})
    _, s = arm(trunk, NIGHT, 3, faults=fault, gate=make_detector("baysal"))
    assert s["gate_flag"].any()
    assert s["cgm_sent"].equals(s["cgm_controller"])


# ---- pathway routing ----
def test_pathway_routing(trunk):
    fault = SingleFault("pisa", NIGHT, {"depth": 48.0, "duration": 70.0})
    _, ctrl = arm(trunk, NIGHT, 2, faults=fault, routing=Routing(person_sees="clean"))
    assert ctrl["cgm_controller"].equals(ctrl["cgm_observed"])
    assert ctrl["cgm_person"].equals(ctrl["cgm_clean"])
    assert (ctrl["cgm_observed"] != ctrl["cgm_clean"]).any()
    _, pers = arm(trunk, NIGHT, 2, faults=fault, routing=Routing(controller_sees="clean"))
    assert pers["cgm_controller"].equals(pers["cgm_clean"])
    assert pers["cgm_person"].equals(pers["cgm_observed"])


# ---- the correcting person ----
def test_correcting_person_follows_its_rule(trunk):
    high = SingleFault(
        "bias", MORNING, {"size_mg_dl": 150.0, "tp_min": 30.0, "hold_until_min": 480.0,
                          "tau_min": 60.0}
    )  # fmt: skip
    loop, s = arm(trunk, MORNING, 10, faults=high, person=CORRECTING)
    r = loop.result()
    ev = r.events[r.events["kind"] == "correction"]
    assert len(ev) >= 1
    rows = s.set_index("minute")
    boluses = r.events[r.events["kind"].isin(["bolus", "correction"])]["minute"].to_numpy()
    for m in ev["minute"]:
        assert rows.loc[m, "cgm_person"] >= 250 and 7 <= (m // 60) % 24 < 23
        earlier = boluses[boluses < m]
        assert earlier.size == 0 or m - earlier.max() >= 180
    ledger = [
        t for t in r.treatments if t["eventType"] == "Correction Bolus" and not t.get("isSMB")
    ]
    assert len(ledger) == len(ev)
    arm(trunk, MORNING, 10, faults=high)  # the default person never corrects
    assert not (loop.result().events["kind"] == "correction").any()


# ---- open loop ----
def test_open_loop_runs_no_controller_and_shows_the_fault_to_the_person():
    loop = make_loop(LoopConfig(PATIENT, controller="open_loop"), SEED, 1, mix=DEFAULT)
    loop.run(NIGHT)
    fault = SingleFault("pisa", NIGHT, {"depth": 48.0, "duration": 70.0})
    r = loop.branch(loop.snapshot(), faults=fault).run(NIGHT + 120).result()
    assert "rT_reason" not in r.samples
    assert np.array_equal(r.insulin_basal, np.full(r.insulin_basal.size, loop.basal))
    s = r.samples[r.samples["minute"] >= NIGHT]
    assert s["cgm_person"].equals(s["cgm_observed"]) and (s["cgm_observed"] != s["cgm_clean"]).any()
