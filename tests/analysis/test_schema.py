"""The arm registry matches memo §6, and the validator catches broken tables."""

from __future__ import annotations

import pandas as pd
import pytest

from ghostlows_aid.analysis import schema
from ghostlows_aid.analysis.schema import ARMS, CLEAN, SchemaError, arms_for, validate


def test_arm_counts_match_the_memo():
    # memo §6: 17 fault arms per controller, 17 variant arms for oref0_tb, tier 1 is 7 arms
    for c in schema.CONTROLLERS:
        faults = [a for a in arms_for(c) if a != CLEAN and not ARMS[a].variant]
        assert len(faults) == 17
        assert len([a for a in faults if ARMS[a].tier == 1]) == 7
    variants = [a for a in arms_for("oref0_tb") if ARMS[a].variant]
    assert len(variants) == 17
    assert not [a for a in arms_for("open_loop") if ARMS[a].variant]
    assert ARMS["clean+baysal_trial"].slot == "night" and ARMS["clean+correcting"].slot == "dose"
    assert ARMS["pisa_typical+matched_filter"].tier == 3 and ARMS["frozen_severe+fakecgm"].tier == 2


def _table() -> pd.DataFrame:
    rows = []
    for arm in (CLEAN, "pisa_typical", "pisa_typical+baysal_trial"):
        r = schema.empty_row(arm)
        r.update(controller="oref0_tb", patient_id="sim:child#004", seed=101, slot="night")
        r["data_source"] = "synthetic"
        rows.append(r)
    return schema.coerce(pd.DataFrame(rows))


def test_a_good_table_validates_in_schema_order():
    t = _table()
    out = validate(t.iloc[:, ::-1])
    assert list(out.columns) == list(schema.COLUMNS)


@pytest.mark.parametrize(
    "breaks, message",
    [
        (lambda t: t.drop(columns="risk_weighted_min"), "missing columns"),
        (lambda t: t.assign(extra=1.0), "unknown columns"),
        (lambda t: pd.concat([t, t.iloc[[1]]]), "duplicate keys"),
        (lambda t: t.assign(slot="any"), "outside its slot"),
        (lambda t: t.assign(controller="open_loop"), "only runs for"),
        (lambda t: t.assign(patient_id="adult#001"), "paper 1's form"),
        (lambda t: t[t["arm"] != CLEAN], "exactly one clean arm"),
        (lambda t: t.assign(fault_kind="pisa"), "fault_kind must be"),
        (lambda t: t.assign(data_source="real"), "data_source"),
        (lambda t: t.assign(seed=t["seed"].astype(float)), "seed: dtype"),
    ],
)
def test_the_validator_names_each_problem(breaks, message):
    with pytest.raises(SchemaError, match=message):
        validate(breaks(_table()))


def test_every_ref_outcome_has_its_twin():
    for name, c in schema.COLUMNS.items():
        if c.effect == "ref":
            assert schema.ref_of(name) in schema.COLUMNS
    assert {"minutes_above_180", "risk_weighted_min", "trend_false_alarms"} <= set(schema.OUTCOMES)
