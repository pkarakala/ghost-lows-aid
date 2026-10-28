"""H1-H2: bootstrap p-values, Holm's step-down and the three-way rule (memo §9)."""

from __future__ import annotations

import numpy as np
import pytest

from ghostlows_aid.analysis import effects as ef
from ghostlows_aid.analysis import hypotheses as hy
from ghostlows_aid.analysis.synthetic import PATIENTS, Design, synthetic_outcomes

N_BOOT = 1000


def run(effects: dict, effect_fn=None, **kw):
    d = Design(
        controllers=("oref0_tb",),
        arms=("pisa_typical", "interference_large"),
        n_seeds=2,
        effects=effects,
        effect_fn=effect_fn,
        seed=5,
        **kw,
    )
    kept, _ = ef.apply_invalid_rule(synthetic_outcomes(d))
    return hy.run(ef.paired_effects(kept), ef.boot_for(kept, N_BOOT))


def test_bootstrap_p_is_the_memo_formula():
    assert hy.bootstrap_p(np.ones(99)) == pytest.approx(2 * 1 / 100)
    assert hy.bootstrap_p(-np.ones(99)) == pytest.approx(2 * 1 / 100)
    d = np.r_[-np.ones(30), np.ones(69)]
    assert hy.bootstrap_p(d) == pytest.approx(2 * 31 / 100)
    assert hy.bootstrap_p(np.zeros(9)) == 1.0  # all at 0: both shares are 1
    assert np.isnan(hy.bootstrap_p(np.full(5, np.nan)))


@pytest.mark.parametrize(
    ("p", "rejected"),
    [
        ({"H1": 0.01, "H2": 0.04}, {"H1": True, "H2": True}),
        ({"H1": 0.03, "H2": 0.04}, {"H1": False, "H2": False}),  # 0.03 > 0.025: stop
        ({"H1": 0.02, "H2": 0.06}, {"H1": True, "H2": False}),
        ({"H1": 0.04, "H2": 0.001}, {"H1": True, "H2": True}),  # H2 first, H1 at 0.05
        ({"H1": 0.025, "H2": 0.05}, {"H1": True, "H2": True}),  # at the level: rejected
        ({"H1": np.nan, "H2": 0.01}, {"H1": False, "H2": True}),
    ],
)
def test_holm_step_down(p, rejected):
    h = hy.holm(p)
    assert h["rejected"].to_dict() == rejected
    first = min(p, key=lambda k: (not np.isfinite(p[k]), p[k]))
    assert h.loc[first, "holm_step"] == 1 and h.loc[first, "holm_alpha"] == pytest.approx(0.025)


@pytest.mark.parametrize(
    ("theta", "rejected", "ci90", "expected"),
    [
        (100.0, True, (80.0, 120.0), hy.HARM),
        (5.0, False, (-20.0, 30.0), hy.NEGLIGIBLE),
        (30.0, True, (20.0, 40.0), hy.NEGLIGIBLE),  # significant, but within the SESOI
        (60.0, True, (40.0, 80.0), hy.INCONCLUSIVE),  # below the SESOI, not equivalent
        (100.0, False, (-10.0, 200.0), hy.INCONCLUSIVE),  # large but not rejected
        (-100.0, True, (-120.0, -80.0), hy.OPPOSITE),  # significant, beyond -SESOI
        (
            -50.0,
            True,
            (-80.0, -20.0),
            hy.INCONCLUSIVE,
        ),  # significant, inside -SESOI, not equivalent
        (np.nan, False, (np.nan, np.nan), hy.NOT_ESTIMABLE),
    ],
)
def test_three_way_rule(theta, rejected, ci90, expected):
    assert hy.decide(theta, rejected, ci90, 72.0) == expected


def test_harm_and_negligible_recovered_with_holm_levels():
    table, co = run(
        {
            ("oref0_tb", "pisa_typical", "minutes_above_180"): 150.0,
            ("oref0_tb", "pisa_typical", "minutes_below_70"): -5.0,
        }
    )
    t = table.set_index("id")
    assert t.loc["H1", "decision"] == hy.HARM
    assert t.loc["H1", "theta"] == pytest.approx(150.0, abs=10)
    assert t.loc["H1", "holm_step"] == 1 and t.loc["H1", "ci_level"] == pytest.approx(0.975)
    assert t.loc["H1", "p_value"] == pytest.approx(2 / (N_BOOT + 1))
    assert t.loc["H2", "decision"] == hy.NEGLIGIBLE and not t.loc["H2", "rejected"]
    assert t.loc["H2", "ci_level"] == pytest.approx(0.95)
    assert t.loc["H1", "ci90_lo"] > t.loc["H1", "ci_lo"]  # the 90% CI is inside the 97.5% CI
    assert set(co["outcome"]) == {"minutes_below_54", "minutes_below_70"}
    pisa = co[(co["id"] == "H1") & (co["outcome"] == "minutes_below_70")].iloc[0]
    assert pisa["theta"] == pytest.approx(-5.0, abs=1) and pisa["level"] == 0.95
    assert {"within_se", "between_sd", "mc_flag", "n_patients"} <= set(table.columns)


def test_both_harmful_both_rejected():
    table, _ = run(
        {
            ("oref0_tb", "pisa_typical", "minutes_above_180"): 150.0,
            ("oref0_tb", "interference_large", "minutes_below_70"): 40.0,
        }
    )
    assert (table["decision"] == hy.HARM).all() and table["rejected"].all()
    assert sorted(table["ci_level"]) == pytest.approx([0.95, 0.975])


def test_inconclusive_when_patients_disagree():
    sign = {p: (1.0 if i % 2 else -1.0) for i, p in enumerate(PATIENTS)}

    def spread(f):  # patient effects +190 or -110: theta 40, a wide CI
        if f["arm"].iloc[0] != "pisa_typical":
            return {}
        return {"minutes_above_180": 40.0 + 150.0 * f["patient_id"].map(sign).to_numpy()}

    table, _ = run({}, effect_fn=spread, zero_noise_arms=True)
    h1 = table.set_index("id").loc["H1"]
    assert h1["theta"] == pytest.approx(40.0)
    assert not h1["rejected"] and h1["decision"] == hy.INCONCLUSIVE
