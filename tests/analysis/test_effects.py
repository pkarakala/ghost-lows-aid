"""Paired effects recover known effects; the invalid rule; the bootstrap is paper 1's, jointly."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from ghostlows.bench.stats import age_group, bootstrap_ci

from ghostlows_aid.analysis import effects as ef
from ghostlows_aid.analysis import schema
from ghostlows_aid.analysis.synthetic import PATIENTS, Design, synthetic_outcomes

ARMS = ("pisa_typical", "interference_large", "dropout_short")
THETA = {
    ("oref0_tb", "pisa_typical", "minutes_above_180"): 80.0,
    ("oref0_tb", "pisa_typical", "risk_weighted_min"): 40.0,
    ("oref0_tb", "interference_large", "minutes_below_70"): 20.0,
    ("oref0_tb", "interference_large", "insulin_during_u"): 3.0,
    ("oref0_tb", "dropout_short", "reconverge_min"): 35.0,
}


def small(**kw) -> Design:
    base = dict(controllers=("oref0_tb",), arms=ARMS, n_seeds=4, effects=THETA, seed=3)
    base.update(kw)
    return Design(**base)


def summary_of(table: pd.DataFrame, n_boot: int = 500) -> pd.DataFrame:
    kept, _ = ef.apply_invalid_rule(table)
    return ef.summarise(ef.paired_effects(kept), ef.boot_for(kept, n_boot))


def test_patients_are_paper1_ids_10_per_age_group():
    assert len(PATIENTS) == 30
    assert pd.Series([age_group(p) for p in PATIENTS]).value_counts().to_dict() == {
        "adolescent": 10,
        "adult": 10,
        "child": 10,
    }
    assert PATIENTS[0] == "sim:adolescent#001"


def test_known_effects_are_recovered_exactly_without_noise():
    s = summary_of(synthetic_outcomes(small(heterogeneity_sd=0.0, seed_sd=0.0)))
    for (c, arm, outcome), theta in THETA.items():
        row = ef.lookup(s, controller=c, arm=arm, outcome=outcome)
        assert row["theta"] == pytest.approx(theta, abs=1e-9)
        assert row["ci_lo"] == pytest.approx(theta, abs=1e-9)
        assert row["ci_hi"] == pytest.approx(theta, abs=1e-9)
    # an outcome without an effect is exactly zero
    row = ef.lookup(s, controller="oref0_tb", arm="pisa_typical", outcome="minutes_below_54")
    assert row["theta"] == pytest.approx(0.0, abs=1e-9)


def test_known_effects_are_recovered_with_heterogeneity_and_seed_noise():
    s = summary_of(synthetic_outcomes(small(n_seeds=6, heterogeneity_sd=0.15, seed_sd=0.3)))
    for (c, arm, outcome), theta in THETA.items():
        row = ef.lookup(s, controller=c, arm=arm, outcome=outcome)
        assert row["ci_lo"] <= theta <= row["ci_hi"], (arm, outcome)
        assert row["n_patients"] == 30 and row["n_events"] == 180
        assert row["ci_hi"] - row["ci_lo"] > 0


def test_clean_only_input_gives_exactly_zero_effects():
    table = synthetic_outcomes(small(arms=()))
    assert set(table["arm"]) == {schema.CLEAN}
    s = summary_of(table)
    diff = [o for o in schema.OUTCOMES if schema.COLUMNS[o].effect == "diff"]
    sel = s[s["outcome"].isin(diff)]
    assert len(sel) and (sel["theta"] == 0).all()
    assert (sel["ci_lo"] == 0).all() and (sel["ci_hi"] == 0).all()


def test_arms_identical_to_clean_give_exactly_zero_effects():
    s = summary_of(synthetic_outcomes(small(effects={}, zero_noise_arms=True)))
    diff = [o for o in schema.OUTCOMES if schema.COLUMNS[o].effect in ("diff", "ref")]
    sel = s[s["outcome"].isin(diff) & s["theta"].notna()]
    assert set(sel["arm"]) >= set(ARMS)
    assert (sel["theta"] == 0).all() and (sel["ci_hi"] == 0).all()


def test_invalid_rule_drops_whole_slots_with_counts():
    p, seed = "sim:child#003", 102
    d = small(invalid={("oref0_tb", p, seed, "night", "pisa_typical"): True})
    table = synthetic_outcomes(d)
    kept, counts = ef.apply_invalid_rule(table)
    unit = (table["patient_id"] == p) & (table["seed"] == seed)
    gone = table[unit & (table["slot"] == "night")]
    assert len(gone) == 2  # clean and pisa_typical of that night slot
    assert len(kept) == len(table) - 2
    left = kept[(kept["patient_id"] == p) & (kept["seed"] == seed)]
    assert set(left["slot"]) == {"dose", "any"}  # the other slots of the trunk stay
    row = counts[(counts["arm"] == "pisa_typical")].iloc[0]
    assert row["invalid_arms"] == 1 and row["dropped"] == 1
    row = counts[(counts["arm"] == "clean") & (counts["slot"] == "night")].iloc[0]
    assert row["invalid_arms"] == 0 and row["dropped"] == 1  # one unit dropped
    # the patient effect is then the mean over the patient's 3 valid seeds
    s = ef.summarise(ef.paired_effects(kept), ef.boot_for(kept, 200))
    row = ef.lookup(s, controller="oref0_tb", arm="pisa_typical", outcome="minutes_above_180")
    assert row["n_events"] == 30 * 4 - 1


def test_bootstrap_ci_is_paper1s_exactly():
    table = synthetic_outcomes(small())
    kept, _ = ef.apply_invalid_rule(table)
    e = ef.paired_effects(kept)
    boot = ef.boot_for(kept, 1000)
    s = ef.summarise(e, boot)
    for arm, outcome in (("pisa_typical", "minutes_above_180"), ("dropout_short", "rescue_count")):
        sel = e[e["arm"] == arm]
        pm = sel.groupby("patient_id")[outcome].mean().reset_index()
        ref = bootstrap_ci(pm, lambda df, o=outcome: df[o].mean(), n_boot=1000, seed=0)
        row = ef.lookup(s, controller="oref0_tb", arm=arm, outcome=outcome)
        assert (row["theta"], row["ci_lo"], row["ci_hi"]) == pytest.approx(ref, rel=1e-12)


def test_replicates_share_one_index_across_quantities():
    boot = ef.Bootstrap(PATIENTS, n_boot=300)
    rng = np.random.default_rng(0)
    x = rng.normal(size=(30, 1))
    _, a = boot.means(x)
    _, b = boot.means(np.c_[2 * x, x + 1.0])
    assert np.allclose(b[:, 0], 2 * a[:, 0], rtol=1e-14) and np.allclose(b[:, 1], a[:, 0] + 1.0)
    # stratified: every replicate draws 10 patients of each age group
    groups = np.array([age_group(p) for p in boot.patients])
    for g in ("adolescent", "adult", "child"):
        assert (boot.counts[:, groups == g].sum(1) == 10).all()
    # the same seed gives the same draws; the first 300 of 1,000 are these 300
    assert np.array_equal(ef.Bootstrap(PATIENTS, n_boot=1000).counts[:300], boot.counts)


def test_monte_carlo_error_terms():
    # two patients with known seed SDs: within SE = mean(sd / sqrt(n))
    events = pd.DataFrame(
        {
            "controller": "oref0_tb",
            "slot": "night",
            "arm": "pisa_typical",
            "patient_id": ["sim:adult#001"] * 4 + ["sim:child#001"] * 4,
            "y": [1.0, 3.0, 1.0, 3.0, 10.0, 10.0, 14.0, 14.0],
        }
    )
    row = ef.estimate(events, ["y"], boot=ef.Bootstrap(events["patient_id"], 50)).iloc[0]
    sd1, sd2 = np.std([1, 3, 1, 3], ddof=1), np.std([10, 10, 14, 14], ddof=1)
    assert row["within_se"] == pytest.approx((sd1 + sd2) / 2 / 2)
    assert row["between_sd"] == pytest.approx(np.std([2.0, 12.0], ddof=1))
    assert row["mc_se"] == pytest.approx(np.sqrt((sd1**2 / 4 + sd2**2 / 4) / 2 / 2))
    assert row["theta"] == pytest.approx(7.0) and not row["mc_flag"]


def test_contrasts_between_arms_and_controllers():
    d = small(
        controllers=("oref0_tb", "oref1_smb"),
        effects={
            ("oref0_tb", "pisa_typical", "risk_weighted_min"): 10.0,
            ("oref1_smb", "pisa_typical", "risk_weighted_min"): 25.0,
            ("oref0_tb", "pisa_typical+baysal_trial", "risk_weighted_min"): 4.0,
        },
        arms=("pisa_typical", "pisa_typical+baysal_trial"),
        heterogeneity_sd=0.0,
        seed_sd=0.0,
    )
    table = synthetic_outcomes(d)
    e = ef.paired_effects(table)
    c = ef.arm_contrast(e, "pisa_typical+baysal_trial", "pisa_typical", ["risk_weighted_min"])
    assert len(c) == 30 * 4 and np.allclose(c["risk_weighted_min"], -6.0)
    k = ef.controller_contrast(e, "oref1_smb", "oref0_tb", ["risk_weighted_min"])
    k = k[k["arm"] == "pisa_typical"]
    assert len(k) == 30 * 4 and np.allclose(k["risk_weighted_min"], 15.0)
