"""SIMULATION.md §10.4 item 5: on a shortened additivity pair (a non-study seed, one day), the
faulted run's displayed CGM equals paper 1's apply_faults applied post hoc to the run's own clean
readings, and the clean run equals the faulted run before the first fault sample. CGM streams,
truth and insulin are compared for equality only; no outcome is computed."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest
from ghostlows.faults.suite import apply_faults

from ghostlows_aid.oref import node_available

ROOT = Path(__file__).resolve().parents[1]
PATIENT, SEED, DAYS = "adult#001", 12345, 1  # outside the study's seeds

pytestmark = pytest.mark.skipif(not node_available(), reason="node or node_modules/oref0 absent")


def test_additivity_pair_equals_post_hoc_injection():
    spec = importlib.util.spec_from_file_location("additivity", ROOT / "experiments/additivity.py")
    add = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(add)
    runs = add.run_pair(PATIENT, 1.0, True, seed=SEED, days=DAYS)
    clean, faulted = runs["clean"], runs["faulted"]
    s = faulted.samples
    trace = s[["time", "isf_glucose", "cgm_noise"]].assign(
        glucose=s["cgm_clean"], patient_id=f"sim:{PATIENT}", seed=SEED, mix=add.MIX
    )
    ref, _ = apply_faults(trace, add.CONFIG, seed=SEED)
    want = ref["cgm_observed"].to_numpy(float)
    want = np.array([math.floor(v + 0.5) if v == v else v for v in want])  # the display rounds
    got = s["cgm_observed"].to_numpy(float)
    assert np.array_equal(np.isnan(got), np.isnan(want))
    assert np.array_equal(got[~np.isnan(got)], want[~np.isnan(want)])
    assert ref["fault_type"].tolist() == s["fault_type"].tolist()
    # the first reading the suite changes (an unlabelled drift changes readings from the start)
    hit = np.flatnonzero(~(s["cgm_observed"] == s["cgm_clean"]).to_numpy())
    assert hit.size > 0 and (s["fault_type"] != "none").any()
    first = int(s["minute"].iloc[hit[0]])
    assert np.array_equal(clean.truth[: first + 1], faulted.truth[: first + 1])
    assert np.array_equal(clean.insulin[: first + 1], faulted.insulin[: first + 1])
    before = slice(0, hit[0])
    for col in ("cgm_observed", "cgm_sent", "rT_reason", "rT_rate"):
        a, b = clean.samples[col].iloc[before], s[col].iloc[before]
        assert a.equals(b), col
