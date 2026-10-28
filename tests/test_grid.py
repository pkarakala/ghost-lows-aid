"""The grid run script: the analysis contract is filled from LoopResults, finished units are skipped
on a restart, and the 10% trace sample is fixed. Open loop, faults switched off: no faulted arm
runs and every arm is the clean continuation, so its outcomes equal the clean arm's."""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pytest
from ghostlows.sim.engine import PATIENTS

from ghostlows_aid import episodes

REPO = Path(__file__).resolve().parents[1]
ANALYSIS_SRC = os.environ.get("GHOSTLOWS_AID_ANALYSIS_SRC") or str(
    REPO.parent / "ghost-lows-aid-analysis" / "src"
)
PATIENT, SEED = "adult#001", 90001


def grid():
    spec = importlib.util.spec_from_file_location("run_grid", REPO / "experiments/run_grid.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_trace_sample_is_a_fixed_tenth():
    g = grid()
    picks = [g.keeps_traces(p, s) for p in PATIENTS for s in range(101, 121)]
    assert 0.06 < np.mean(picks) < 0.14
    assert picks == [g.keeps_traces(p, s) for p in PATIENTS for s in range(101, 121)]


def test_contract_filled_from_loop_results():
    g = grid()
    try:
        outcomes, schema = g.load_analysis(ANALYSIS_SRC)
    except SystemExit:
        pytest.skip("the analysis package is not available")
    chosen = [episodes.ARM_BY_NAME[n] for n in ("pisa_typical", "dropout_short")]
    unit, runs = episodes.run_unit("open_loop", PATIENT, SEED, 0, 20, arms=chosen,
                                   with_faults=False)  # fmt: skip
    by_slot = {}
    for run in runs:
        by_slot.setdefault(run.arm.slot, []).append(g.to_trace(unit, run, outcomes, schema))
    tables = {slot: outcomes.slot_outcomes(ts) for slot, ts in by_slot.items() if len(ts) > 1}
    assert set(tables) == {"night", "any"}
    for t in tables.values():
        schema.validate(t)
        clean, arm = t[t["arm"] == "clean"].iloc[0], t[t["arm"] != "clean"].iloc[0]
        for col in ("minutes_below_70", "minutes_above_180", "mean_glucose", "insulin_total_u"):
            assert clean[col] == arm[col]


def test_resume_skips_finished_units(tmp_path, monkeypatch):
    g = grid()
    monkeypatch.setattr(g, "keeps_traces", lambda p, s: True)
    job = (tmp_path, "open_loop", PATIENT, SEED, 0, 20, 0, False, None, "rule1800")
    assert g.run_job(job) == "done"
    record = json.loads(next((tmp_path / "units").glob("*.json")).read_text())
    assert record["arm_s"] == {} and record["trunk_s"] > record["warmup_s"] > 0
    assert len(list((tmp_path / "traces").glob("*.pkl.gz"))) == 1
    assert g.run_job(job) == "skip"
    assert not list(tmp_path.rglob("*.tmp"))
