"""Secondary analyses recover constructed effects (memo §10)."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from ghostlows_aid.analysis import effects as ef
from ghostlows_aid.analysis import secondary as se
from ghostlows_aid.analysis.cross_layer import RATES
from ghostlows_aid.analysis.schema import GATES, SECOND_DETECTOR
from ghostlows_aid.analysis.synthetic import PATIENTS, Design, synthetic_outcomes

H = se.HARM
EXACT = dict(patient_sd=0.0, heterogeneity_sd=0.0, seed_sd=0.0)
THETA = {
    "pisa_typical": 10.0,
    "pisa_severe": 30.0,
    "dropout_short": 1.0,
    "dropout_long": 6.0,
    "dropout_severe": 5.0,
    "frozen_typical": 2.0,
    "frozen_severe": 20.0,
    "spike_typical": 3.0,
    "drift_low": 4.0,
    "drift_high": 5.0,
    "frozen_severe+fakecgm": 50.0,
    "pisa_typical+controller_only": 7.0,
    "pisa_typical+person_only": 2.0,
}


def prepared(n_boot=200, **kw):
    base = dict(n_seeds=2, seed=2)
    base.update(kw)
    kept, _ = ef.apply_invalid_rule(synthetic_outcomes(Design(**base)))
    return kept, ef.paired_effects(kept), ef.boot_for(kept, n_boot)


@pytest.fixture(scope="module")
def exact():
    eff = {("oref0_tb", a, H): v for a, v in THETA.items()}
    eff |= {("oref1_smb", a, H): 3 * v for a, v in THETA.items()}
    return prepared(effects=eff, zero_noise_arms=True)


def row(frame, **sel):
    return ef.lookup(frame, **sel)


def test_ladders_and_guards(exact):
    _, e, boot = exact
    lad = se.ladders(e, boot)
    r = row(
        lad,
        controller="oref0_tb",
        ladder="compression",
        arm="pisa_severe",
        outcome=H,
        versus="pisa_typical",
    )
    assert r["theta"] == pytest.approx(20.0) and r["ci_lo"] == pytest.approx(20.0)
    r = row(
        lad,
        controller="oref0_tb",
        ladder="dropout",
        arm="dropout_severe",
        outcome=H,
        versus="dropout_long",
    )
    assert r["theta"] == pytest.approx(-1.0)
    g = se.guards(e, boot)
    assert row(g, controller="oref0_tb", arm="flat_cgm_check", outcome=H)["theta"] == pytest.approx(
        30.0
    )
    assert row(g, controller="oref0_tb", arm="stale_data_stop", outcome=H)[
        "theta"
    ] == pytest.approx(15.0)
    assert row(g, controller="oref1_smb", arm="stale_data_stop", outcome=H)[
        "theta"
    ] == pytest.approx(45.0)


def test_pathway_decomposition_and_wake_curve(exact):
    _, e, boot = exact
    cells, wake = se.pathway(e, boot)
    get = lambda comp: row(cells, arm="pisa_typical", component=comp, outcome=H)["theta"]  # noqa: E731
    assert (get("total"), get("controller_only"), get("person_only")) == pytest.approx((10, 7, 2))
    assert get("interaction") == pytest.approx(1.0)
    for q in se.WAKE_Q:
        assert row(wake, q=q, outcome=H)["theta"] == pytest.approx(q * 10 + (1 - q) * 7)


def test_controller_contrast_and_person_year(exact):
    _, e, boot = exact
    contrast, _ = se.controllers(e, boot)
    assert row(contrast, arm="pisa_typical", outcome=H)["theta"] == pytest.approx(20.0)
    py = se.person_year(e, boot)
    central = 0.86 * 10 + 2 * 1 + 1 * 2 + 1 * 3 + 1 * 4 + 1 * 5
    low = 0.86 * 0.03 / 0.321 * 10 + 0.13 * 6 + 2 + 3 + 4 + 5
    r = row(py, controller="oref0_tb", scenario="central", outcome=H)
    assert r["per_person_year"] == pytest.approx(365 * central)
    assert row(py, controller="oref0_tb", scenario="low", outcome=H)[
        "per_person_year"
    ] == pytest.approx(365 * low)
    assert row(py, controller="oref1_smb", scenario="central", outcome=H)[
        "per_person_year"
    ] == pytest.approx(3 * 365 * central)


def test_pairing_efficiency(exact):
    kept, _, _ = exact
    pe = se.pairing_efficiency(kept)
    r = row(pe, controller="oref0_tb", arm="pisa_typical", outcome=H)
    assert r["median_ratio"] == 0.0 and math.isinf(r["implied_runs"])  # exact pairing
    noisy, _, _ = prepared(
        controllers=("oref0_tb",), arms=("pisa_typical",), n_seeds=6, heterogeneity_sd=0.0
    )
    r = row(se.pairing_efficiency(noisy), arm="pisa_typical", outcome=H)
    assert 0 < r["median_ratio"] < 1 and r["implied_runs"] > 1 and r["n_patients"] == 30


def test_dose_response_slope_and_crossing_split():
    def dose(f):
        if f["arm"].iloc[0] != "pisa_typical":
            return {}
        return {H: 2.0 * f["fault_area_mg_dl_h"].to_numpy(float)}

    _, e, boot = prepared(
        controllers=("oref0_tb",),
        arms=("pisa_typical",),
        n_seeds=4,
        effect_fn=dose,
        zero_noise_arms=True,
    )
    slopes, split = se.dose_response(e, boot)
    assert row(slopes, outcome=H)["slope_per_mg_dl_h"] == pytest.approx(2.0)
    assert set(split["crossed_70"]) == {True, False}


def test_state_cells_and_estimability():
    def low_start(f):
        if f["arm"].iloc[0] != "pisa_typical":
            return {}
        return {H: np.where(f["onset_truth_mg_dl"].to_numpy(float) < 100, 50.0, 0.0)}

    _, e, boot = prepared(
        controllers=("oref0_tb",),
        arms=("pisa_typical",),
        n_seeds=8,
        effect_fn=low_start,
        zero_noise_arms=True,
    )
    cells = se.state_cells(e, boot)
    g = cells[(cells["grouping"] == "glucose_x_trend") & cells["estimable"]]
    assert len(g) and (g.loc[g["cell"].str.startswith("<100"), "theta"] == 50.0).all()
    assert (g.loc[~g["cell"].str.startswith("<100"), "theta"] == 0.0).all()
    assert (cells.loc[~cells["estimable"], "theta"].isna()).all()
    assert (cells.loc[~cells["estimable"], "n_patients"] < 5).all()
    _, e4, boot4 = prepared(
        controllers=("oref0_tb",), arms=("pisa_typical",), patients=PATIENTS[::8]
    )
    few = se.state_cells(e4, boot4)
    assert len(few) and not few["estimable"].any() and few["theta"].isna().all()


def test_state_model_recovers_coefficients_and_held_out_fit():
    def linear(f):
        if f["arm"].iloc[0] != "pisa_typical":
            return {}
        x = f[["onset_truth_mg_dl", "onset_iob_u"]].to_numpy(float)
        return {H: 0.5 * x[:, 0] + 3.0 * x[:, 1]}

    _, e, boot = prepared(
        controllers=("oref0_tb",),
        arms=("pisa_typical",),
        n_seeds=6,
        effect_fn=linear,
        zero_noise_arms=True,
    )
    coefs, held = se.state_model(e, boot)
    c = coefs.set_index("term")["coef"]
    assert c["onset_truth_mg_dl"] == pytest.approx(0.5) and c["onset_iob_u"] == pytest.approx(3.0)
    assert c[["onset_cgm_trend", "onset_cob_g", "onset_hours_since_bolus"]].abs().max() < 1e-8
    assert held.iloc[0]["r2_out_of_fold"] > 0.99 and held.iloc[0]["n_folds"] == 5


GATE_THETA = {"pisa_typical": 100.0}
GATE_THETA |= {"pisa_typical+baysal_aggressive": 40.0, "clean+baysal_aggressive": 30.0}
GATE_THETA |= {"pisa_typical+baysal_trial": 50.0, "clean+baysal_trial": 5.0}
GATE_THETA |= {"pisa_typical+baysal_nominal": 70.0, "pisa_typical+baysal_cautious": 90.0}
GATE_THETA |= {"pisa_typical+matched_filter": 95.0}
FA = {"trend": 1.0, "baysal_aggressive": 0.0, "baysal_trial": 0.1, "baysal_nominal": 0.3}
FA |= {"baysal_cautious": 0.6, "matched_filter": 0.2}


def gate_case(sensitivity):
    eff = {("oref0_tb", a, H): v for a, v in GATE_THETA.items()}
    eff |= {("open_loop", "pisa_typical", f"{p}_false_alarms"): v for p, v in FA.items()}
    gates = (*GATES, SECOND_DETECTOR)
    arms = ("pisa_typical", *(f"{b}+{g}" for b in ("pisa_typical", "clean") for g in gates))
    sens = {p: (0.8, 0.9) for p in FA} | sensitivity
    kept, e, boot = prepared(
        controllers=("oref0_tb", "open_loop"), arms=arms, sensitivity=sens, effects=eff, **EXACT
    )
    return se.gate(e, kept, boot)


def test_gate_net_effect_and_optima():
    curves, optima = gate_case({"baysal_aggressive": (0.8, 0.5)})  # aggressive misses lows
    c = curves.set_index("option")
    w = RATES["pisa_typical"]
    assert c.loc["aggressive", "net_per_night"] == pytest.approx(w * -60 + (1 - w) * 30)
    assert c.loc["trial", "net_per_night"] == pytest.approx(w * -50 + (1 - w) * 5)
    assert c.loc["cautious", "net_per_night"] == pytest.approx(w * -10)
    assert c.loc["no_gate", "net_per_night"] == 0.0
    assert c.loc["matched_filter", "delta_pisa"] == pytest.approx(-5.0)
    assert c.loc["no_gate", "false_alarms_per_night"] == pytest.approx(w * 1.4 + (1 - w) * 0.4)
    assert c.loc["aggressive", "l2_sensitivity"] == pytest.approx(0.5)
    o = optima.iloc[0]
    assert o["harm_optimal"] == "aggressive" and o["alarm_optimal"] == "trial"
    assert o["alarm_rule_feasible"] and o["alarm_rule_level"] == 2
    assert o["coincide_share"] == 0.0
    _, same = gate_case({})  # every setting as sensitive as the ungated alarm
    assert same.iloc[0]["alarm_optimal"] == "aggressive" and same.iloc[0]["coincide_share"] == 1.0


def test_mechanism_shares_are_proportions(exact):
    _, e, boot = exact
    m = se.mechanism(e, boot)
    s = m[m["outcome"].str.startswith("share_")]
    assert len(s) and s["theta"].between(0, 1).all()
    assert not (m["controller"].eq("open_loop") & m["outcome"].str.startswith("share_")).any()


def test_run_all_on_the_full_design():
    kept, e, boot = prepared(n_boot=100)
    out = se.run(e, kept, boot)
    assert all(isinstance(v, pd.DataFrame) for v in out.values())
    assert all(len(v) for k, v in out.items()), [k for k, v in out.items() if not len(v)]
