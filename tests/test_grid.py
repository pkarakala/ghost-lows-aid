"""The grid run script: the analysis contract is filled from LoopResults, finished units are skipped
on a restart, and the 10% trace sample is fixed. Open loop, faults switched off: no faulted arm
runs and every arm is the clean continuation, so its outcomes equal the clean arm's."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from ghostlows.sim.engine import PATIENTS

from ghostlows_aid import episodes

REPO = Path(__file__).resolve().parents[1]
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
    from ghostlows_aid.analysis import outcomes, schema

    chosen = [episodes.ARM_BY_NAME[n] for n in ("pisa_typical", "dropout_short")]
    unit, runs = episodes.run_unit("open_loop", PATIENT, SEED, 0, 20, arms=chosen,
                                   with_faults=False)  # fmt: skip
    by_slot = {}
    for run in runs:
        by_slot.setdefault(run.arm.slot, []).append(outcomes.trace_from_run(unit, run))
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
    job = (tmp_path, "open_loop", PATIENT, SEED, 0, 20, 0, False, "rule1800", 1.0, {})
    assert g.run_job(job) == "done"
    record = json.loads(next((tmp_path / "units").glob("*.json")).read_text())
    assert record["arm_s"] == {} and record["trunk_s"] > record["warmup_s"] > 0
    assert len(list((tmp_path / "traces").glob("*.pkl.gz"))) == 1
    assert g.run_job(job) == "skip"
    assert not list(tmp_path.rglob("*.tmp"))


def test_faults_need_a_frozen_file_naming_the_tag_and_a_tagged_head(tmp_path):
    g = grid()
    assert not g.freeze_ok(None)[0]
    f = tmp_path / "DECISIONS.md"
    f.write_text("no tag here")
    assert "does not record" in g.freeze_ok(str(f))[1]
    f.write_text("D13: protocol-v1 tagged")
    ok, why = g.freeze_ok(str(f))
    tagged = not g._git("rev-parse", "--verify", "--quiet", "refs/tags/protocol-v1").returncode
    assert ok == tagged and (ok or "protocol-v1" in why)


def test_markers_scale_map_and_trace_sample(tmp_path):
    g = grid()
    smap = tmp_path / "map.json"
    smap.write_text(json.dumps({"isf_scale": {PATIENT: 1.25}}))
    base = ["run", "--controllers", "open_loop", "--max-tier", "0", "--seeds", str(SEED),
            "--isf-scale-map", str(smap), "--workers", "1",
            "--out", str(tmp_path / "g")]  # fmt: skip
    with pytest.raises(SystemExit, match="missing from"):
        g.main([*base, "--patients", PATIENT, "child#001"])
    assert not (tmp_path / "g" / "units").exists()  # stopped before any unit
    g.main([*base, "--patients", PATIENT])
    rec = json.loads(next((tmp_path / "g" / "units").glob("*.json")).read_text())
    assert rec["isf_scale"] == 1.25 and len(rec["harness_commit"]) == 40
    assert rec["versions"]["oref0"] == "0.7.1" and rec["versions"]["ghostlows"] == "0.4.0"
    assert len(rec["isf_scale_map_sha256"]) == 64
    sample = pd.read_csv(tmp_path / "g" / "trace_sample.csv")
    assert list(sample.columns) == ["patient", "seed"]
    g.main([*base, "--patients", PATIENT])  # the same command skips the finished unit
    with pytest.raises(SystemExit, match="differ from this command"):
        g.main([*base, "--patients", PATIENT, "--max-tier", "1"])
