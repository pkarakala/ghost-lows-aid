"""Each outcome of memo §8 on small hand-made traces."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from ghostlows_aid.analysis import outcomes as oc
from ghostlows_aid.analysis.schema import validate

from .helpers import BASAL, make_trace, window_values

PRE = oc.PRE_MIN

# real rT.reason strings from oref0 0.7.1 (node/parity_fixtures.js output)
STALE = (
    "If current system time Mon Jan 05 2026 12:00:00 GMT+0000 (Coordinated Universal Time) is "
    "correct, then BG data is too old. The last BG data was read 20m ago at Mon Jan 05 2026 "
    "11:40:00 GMT+0000 (Coordinated Universal Time). Replacing high temp basal of 1.5 with "
    "neutral temp of 0.9"
)
FLAT = (
    "CGM data is unchanged (150+0) for 5m w/ 0.4 mg/dL ~15m change & -0.6 mg/dL ~45m change. "
    "Replacing high temp basal of 1.5 with neutral temp of 0.9"
)
LGS = (
    "COB: 30, Dev: -5, BGI: -2, ISF: 40, CR: 10, minPredBG 54, minGuardBG 34, IOBpredBG 39, "
    "COBpredBG 107, UAMpredBG 39; minGuardBG 34<75"
)
MAXDELTA = (
    "COB: 30, Dev: 133, BGI: -2, ISF: 40, CR: 10, minPredBG 160, minGuardBG 130, IOBpredBG 100, "
    "COBpredBG 160, UAMpredBG 401; maxDelta 30 > 20% of BG 120: SMB disabled; Eventual BG 401 >= "
    "120, max_iob 2.5, temp 1.5<2.2U/hr. "
)
EVENTUAL_ZERO = (
    "COB: 0, Dev: -30, BGI: 0, ISF: 40, CR: 10, minPredBG 999, minGuardBG 999, IOBpredBG 75; "
    "Eventual BG 45 < 110, setting 60m zero temp. "
)
MISMATCH = "Warning: currenttemp rate 1.5 != lastTemp rate 0 from pumphistory; canceling temp"
MISMATCH2 = (
    "Warning: currenttemp running but lastTemp from pumphistory ended 12m ago; canceling temp"
)


def test_minutes_in_ranges_and_risk_weights():
    g = window_values(720, [(0, 10, 50), (10, 30, 60), (30, 60, 200), (60, 100, 300)])
    r = oc.glycaemia(g)
    assert r["minutes_below_54"] == 10 and r["minutes_below_70"] == 30
    assert r["minutes_above_180"] == 70 and r["minutes_above_250"] == 40
    assert r["risk_weighted_min"] == pytest.approx(3.0 * 10 + 2.4 * 20 + 1.6 * 40 + 0.8 * 30)
    assert r["peak_glucose"] == 300 and r["nadir_glucose"] == 50
    assert r["mean_glucose"] == pytest.approx(g.mean())


def test_range_boundaries_follow_the_consensus():
    # 54 is low (54-69), 70 is in range, 180 in range, 250 is high (181-250), not very high
    r = oc.glycaemia(np.array([54.0, 70.0, 180.0, 250.0, 53.999, 250.001]))
    assert r["minutes_below_54"] == 1 and r["minutes_below_70"] == 2
    assert r["minutes_above_180"] == 2 and r["minutes_above_250"] == 1
    assert r["risk_weighted_min"] == pytest.approx(3.0 + 2.4 + 1.6 + 0.8)


def test_hypo_events_need_15_minutes_and_an_onset_in_the_window():
    n = PRE + 720
    g = window_values(
        n,
        [
            (PRE - 10, PRE + 10, 65),  # began before t0: not counted
            (PRE + 100, PRE + 115, 65),  # 15 min: a level 1 event
            (PRE + 200, PRE + 214, 60),  # 14 min: no event
            (PRE + 300, PRE + 340, 50),  # level 1 and level 2
        ],
    )
    assert oc.hypo_events(g, PRE) == {"hypo_l1_events": 2.0, "hypo_l2_events": 1.0}


def test_insulin_by_route_and_temp_basal_deviation():
    n = PRE + 720
    basal = window_values(n, [(PRE + 60, PRE + 120, 0.0)], fill=BASAL)  # a 60-min zero temp
    smb = window_values(n, [(PRE + 5, PRE + 6, 0.3), (PRE + 10, PRE + 11, 0.2)], fill=0.0)
    bolus = window_values(n, [(PRE + 30, PRE + 31, 4.0), (PRE - 30, PRE - 29, 9.0)], fill=0.0)
    rescue = window_values(n, [(PRE + 400, PRE + 401, 15.0), (PRE + 415, PRE + 416, 15.0)], 0.0)
    corr = window_values(n, [(PRE + 500, PRE + 501, 1.5)], fill=0.0)
    t = make_trace(basal_u=basal, smb_u=smb, bolus_u=bolus, rescue_g=rescue, correction_u=corr)
    r = oc.insulin(t)
    assert r["insulin_tb_dev_u"] == pytest.approx(-60 * BASAL)
    assert r["insulin_smb_u"] == pytest.approx(0.5)
    assert r["insulin_bolus_u"] == pytest.approx(4.0)  # the bolus before t0 is not in the window
    assert r["insulin_correction_u"] == 1.5 and r["correction_boluses"] == 1
    assert r["insulin_total_u"] == pytest.approx(BASAL * 660 + 0.5 + 4.0 + 1.5)
    assert r["rescue_carbs_g"] == 30 and r["rescue_count"] == 2


def test_during_and_after_the_fault_with_the_clean_arm_over_the_same_minutes():
    n = PRE + 720
    clean = make_trace(truth=window_values(n, [(PRE + 100, PRE + 200, 200)]))
    arm = make_trace(
        "pisa_typical",
        fault_end_min=29.5,  # minutes 0-29 are during the fault
        basal_u=window_values(n, [(PRE, PRE + 30, 0.0)], fill=BASAL),
        truth=window_values(n, [(PRE + 20, PRE + 40, 190), (PRE + 100, PRE + 300, 220)]),
    )
    r = oc.fault_windows(arm, clean)
    assert r["insulin_during_u"] == 0 and r["ref_insulin_during_u"] == pytest.approx(30 * BASAL)
    assert r["insulin_after_u"] == pytest.approx(690 * BASAL) == r["ref_insulin_after_u"]
    assert r["minutes_above_180_during"] == 10 and r["ref_minutes_above_180_during"] == 0
    assert r["minutes_above_180_after"] == 10 + 200 and r["ref_minutes_above_180_after"] == 100
    held = make_trace("drift_low", fault_end_min=math.inf)
    r = oc.fault_windows(held, make_trace(slot="any"))
    assert r["insulin_after_u"] == 0 and r["insulin_during_u"] == pytest.approx(720 * BASAL)
    assert all(math.isnan(v) for v in oc.fault_windows(clean, None).values())


def test_time_to_reconverge():
    n = PRE + 720
    clean = make_trace()
    arm = make_trace(
        "pisa_typical", fault_end_min=30, truth=window_values(n, [(PRE, PRE + 90, 140)])
    )
    assert oc.reconvergence(arm, clean) == {"reconverge_min": 60.0, "reconverged": True}
    # differences below 5 mg/dL count as converged; a re-divergence inside 60 min restarts it
    g = window_values(
        n, [(PRE, PRE + 90, 140), (PRE + 90, PRE + 120, 124), (PRE + 120, PRE + 125, 126)]
    )
    arm = make_trace("pisa_typical", fault_end_min=30, truth=g)
    assert oc.reconvergence(arm, clean)["reconverge_min"] == 125 - 30
    # already converged when the fault ends: 0; never: censored
    assert (
        oc.reconvergence(make_trace("pisa_typical", fault_end_min=30), clean)["reconverge_min"] == 0
    )
    never = make_trace("pisa_typical", fault_end_min=30, truth=window_values(n, [(PRE, n, 150)]))
    r = oc.reconvergence(never, clean)
    assert math.isnan(r["reconverge_min"]) and r["reconverged"] is False
    late = make_trace("pisa_typical", fault_end_min=700)  # fewer than 60 min left: censored
    assert oc.reconvergence(late, clean)["reconverged"] is False
    held = make_trace("drift_low", fault_end_min=math.inf, slot="any")
    assert oc.reconvergence(held, make_trace(slot="any"))["reconverged"] is pd.NA
    assert oc.reconvergence(clean, None)["reconverged"] is pd.NA


def test_reason_strings_map_to_the_shipped_rules():
    counts = oc.classify_reasons(
        [STALE, FLAT, LGS, MAXDELTA, EVENTUAL_ZERO, MISMATCH, MISMATCH2, "", "no CGM data yet"]
    )
    assert counts == {
        "actions_stale": 1,
        "actions_flat": 1,
        "actions_temp_mismatch": 2,
        "actions_lgs_zero_temp": 1,  # "Eventual BG < target, zero temp" is not the LGS branch
        "actions_smb_maxdelta": 1,
    }
    # the header "minGuardBG 130, IOBpredBG" is not a suspend, mmol/L output is
    assert oc.classify_reasons(["... minGuardBG 3.6<3.9"])["actions_lgs_zero_temp"] == 1


def test_actions_count_cycles_in_the_window_only():
    n5 = (PRE + 720) // 5
    reasons = [""] * n5
    reasons[PRE // 5 - 1] = FLAT  # before t0
    reasons[PRE // 5] = FLAT
    reasons[PRE // 5 + 3] = LGS
    t = make_trace(reasons=reasons)
    r = oc.actions(t)
    assert r["actions_flat"] == 1 and r["actions_lgs_zero_temp"] == 1
    assert all(math.isnan(v) for v in oc.actions(make_trace(controller="open_loop")).values())


def test_crossings_of_the_controller_input_during_the_fault():
    n5 = (PRE + 720) // 5
    x = window_values(n5, [(PRE // 5 + 2, PRE // 5 + 3, 65), (PRE // 5 + 20, PRE // 5 + 21, 190)])
    t = make_trace("pisa_typical", fault_end_min=45, cgm_controller=x)
    assert oc.crossings(t) == {"ctrl_below_70_during": True, "ctrl_above_180_during": False}
    longer = make_trace("pisa_typical", fault_end_min=120, cgm_controller=x)
    assert oc.crossings(longer)["ctrl_above_180_during"] is True
    assert oc.crossings(make_trace())["ctrl_below_70_during"] is pd.NA


def test_onset_state_from_the_clean_arm():
    n, n5 = PRE + 720, (PRE + 720) // 5
    g = np.full(n, 120.0)
    g[PRE - 15] = 135.0
    cgm = np.full(n5, 110.0)
    cgm[PRE // 5 - 3] = 140.0
    iob = np.full(n5, 1.0)
    iob[PRE // 5] = 2.5
    t = make_trace(truth=g, cgm_display=cgm, iob_u=iob, onset_minute=24 * 60 + 23 * 60 + 30)
    s = oc.onset_state(t)
    assert s["onset_truth_mg_dl"] == 120 and s["onset_truth_trend"] == pytest.approx(-1.0)
    assert s["onset_cgm_mg_dl"] == 110 and s["onset_cgm_trend"] == pytest.approx(-2.0)
    assert s["onset_iob_u"] == 2.5 and s["onset_hour"] == pytest.approx(23.5)
    assert s["onset_hours_since_bolus"] == pytest.approx(4.5)


def test_realised_fault_parameters_and_artifact_area():
    t = make_trace("pisa_typical", fault={"depth": 24.0, "duration": 45.0, "tau": 10.0})
    f = oc.fault_columns(t)
    assert f["fault_depth_mg_dl"] == 24 and f["fault_duration_min"] == 45
    assert f["fault_area_mg_dl_h"] == pytest.approx(18.0) and f["fault_tau_min"] == 10
    assert math.isnan(oc.fault_columns(make_trace())["fault_area_mg_dl_h"])


def test_invalid_below_the_validity_floor():
    n = PRE + 720
    assert oc.is_invalid(make_trace(truth=window_values(n, [(PRE + 300, PRE + 301, 9.9)])))
    assert not oc.is_invalid(make_trace())
    assert oc.is_invalid(make_trace(valid=False))


def test_slot_outcomes_rows_validate_and_share_the_onset_state():
    n = PRE + 720
    clean = make_trace()
    arm = make_trace(
        "pisa_typical", fault_end_min=60.0, truth=window_values(n, [(PRE, PRE + 60, 60)])
    )
    rows = oc.slot_outcomes([clean, arm], alarm_layer=False)
    rows = validate(rows)
    c, a = rows.iloc[0], rows.iloc[1]
    assert a["minutes_below_70"] == 60 and c["minutes_below_70"] == 0
    assert a["onset_truth_mg_dl"] == c["onset_truth_mg_dl"]
    assert a["ref_insulin_during_u"] == pytest.approx(60 * BASAL)
    assert pd.isna(c["reconverged"]) and bool(a["reconverged"])


def test_slot_outcomes_refuse_arms_with_another_history():
    n = PRE + 720
    other = make_trace("pisa_typical", truth=window_values(n, [(0, 10, 100)]))
    with pytest.raises(ValueError, match="history"):
        oc.slot_outcomes([make_trace(), other], alarm_layer=False)
    with pytest.raises(ValueError, match="clean"):
        oc.slot_outcomes([other], alarm_layer=False)
    with pytest.raises(ValueError, match="another unit"):
        oc.slot_outcomes([make_trace(), make_trace("pisa_typical", seed=102)], alarm_layer=False)


def test_trace_lengths_are_checked():
    with pytest.raises(ValueError, match="truth"):
        make_trace(truth=np.zeros(10))
