"""D8 additivity (secondary.additivity), the S1 and S2 options of experiments/analyse.py and the
per-day bookkeeping of experiments/additivity.py, on synthetic inputs only."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from ghostlows_aid.analysis import effects as ef
from ghostlows_aid.analysis import secondary
from ghostlows_aid.analysis.synthetic import Design, synthetic_outcomes

ROOT = Path(__file__).resolve().parents[2]
OUTCOMES = list(secondary.ADDITIVITY_OUTCOMES)
PATIENTS = [f"sim:{g}#{i:03d}" for g in ("adolescent", "adult", "child") for i in (1, 2, 3)]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def synthetic_case(extra: float = 0.0, seed: int = 0):
    """Per-event arm effects, and records whose faulted days are clean + the additive sum
    (+ `extra` on risk-weighted minutes)."""
    rng = np.random.default_rng(seed)
    arms = [*secondary.ADDITIVITY_ARMS.values(), *secondary.DRIFT_ARMS.values()]
    eff_rows, records = [], []
    for p in PATIENTS:
        unit = {a: rng.uniform(1, 20, len(OUTCOMES)) for a in arms}
        for a in arms:
            for s in (101, 102):
                eff_rows.append({"controller": "oref0_tb", "patient_id": p, "seed": s, "arm": a,
                                 **dict(zip(OUTCOMES, unit[a], strict=True))})  # fmt: skip
        sign = rng.choice([-1.0, 1.0])
        events, clean, faulted = [], [], []
        for d in range(2, 8):
            counts = {k: int(rng.integers(0, 4)) for k in secondary.ADDITIVITY_ARMS}
            events.append({"day": d, **counts, "drift_sign": sign, "drift_size": 0.1 * sign})
            base = rng.uniform(0, 300, len(OUTCOMES))
            add = sum(counts[k] * unit[a] for k, a in secondary.ADDITIVITY_ARMS.items())
            add = add + 2 * unit[secondary.DRIFT_ARMS[sign]]
            add = add + np.array([extra] + [0.0] * (len(OUTCOMES) - 1))
            clean.append({"day": d, **dict(zip(OUTCOMES, base, strict=True))})
            faulted.append({"day": d, **dict(zip(OUTCOMES, base + add, strict=True))})
        records.append({"patient": p, "synthetic": True, "events": events,
                        "clean": {"days": clean}, "faulted": {"days": faulted}})  # fmt: skip
    return pd.DataFrame(eff_rows), records


@pytest.mark.parametrize("extra", [0.0, 50.0])
def test_additivity_recovers_an_additive_suite_and_flags_an_excess(extra):
    effects, records = synthetic_case(extra)
    daily, events = secondary.read_additivity(records)
    out = secondary.additivity(daily, events, effects, ef.Bootstrap(PATIENTS, n_boot=200))
    assert len(out) == 3 * len(OUTCOMES)
    d = out[out["quantity"] == "difference"].set_index("outcome")
    assert d["n_patients"].eq(len(PATIENTS)).all()
    rest = [o for o in OUTCOMES if o != secondary.HARM]
    assert np.allclose(d.loc[rest, ["theta", "ci_lo", "ci_hi"]], 0, atol=1e-9)
    assert d.loc[secondary.HARM, "theta"] == pytest.approx(extra)
    if extra:  # an exact sum has a zero-width CI at 0 up to rounding, so only the excess is flagged
        assert not d.loc[secondary.HARM, "additive"]
    assert out.loc[out["quantity"] != "difference", "additive"].isna().all()


def test_analyse_additivity_s1_and_s2(tmp_path):
    analyse = load("analyse")
    table = synthetic_outcomes(Design(controllers=("oref0_tb",), max_tier=1, n_seeds=2))
    # S1: the literal-scale grid in its own folder, read with --input and its own --out
    grid = tmp_path / "grid_s1"
    grid.mkdir()
    scales = json.loads((ROOT / "results/calibration/selection.json").read_text())["isf_scale"]
    table["isf_scale"] = table["patient_id"].str[4:].map(scales)  # the literal reading (S1)
    table.to_parquet(grid / "units.parquet")
    _, records = synthetic_case()
    folder = tmp_path / "additivity"
    folder.mkdir()
    for r in records:
        (folder / f"{r['patient'][4:].replace('#', '')}.json").write_text(json.dumps(r))
    smap = tmp_path / "map.json"
    smap.write_text(json.dumps({"outside_acceptance": ["adult#002", "child#003"]}))
    out = tmp_path / "s1"
    args = ["--input", str(grid), "--out", str(out), "--no-figures", "--n-boot", "50",
            "--n-boot-confirmatory", "50", "--additivity", str(folder)]  # fmt: skip
    assert analyse.main(args) == 0
    add = pd.read_csv(out / "additivity.csv")
    assert len(add) == 3 * len(OUTCOMES) and add["n_patients"].max() == len(PATIENTS)
    assert "D8 additivity" in (out / "summary.md").read_text()
    # S2: the patients outside acceptance are dropped before every estimate, and it says so
    out2 = tmp_path / "s2"
    args2 = [*args[:3], str(out2), *args[4:], "--exclude-outside-acceptance", str(smap)]
    assert analyse.main(args2) == 0
    assert "S2: 2 patients outside acceptance excluded" in (out2 / "summary.md").read_text()
    manifest = json.loads((out2 / "manifest.json").read_text())
    assert manifest["excluded_outside_acceptance"] == ["sim:adult#002", "sim:child#003"]
    add2 = pd.read_csv(out2 / "additivity.csv")
    assert add2["n_patients"].max() == len(PATIENTS) - 2


def test_additivity_script_bookkeeping():
    add = load("additivity")
    sch = add.schedule("adult#001")  # a CGM-level schedule only; no loop runs
    rows = add.events_per_day(sch)
    assert [r["day"] for r in rows] == list(range(2, 8))
    for kind in add.EVENT_TYPES:
        n = sum(1 for i in sch.instances if i["fault_type"] == kind and 1440 <= i["start"] < 10080)
        assert sum(r[kind] for r in rows) == n
    assert len({r["drift_sign"] for r in rows}) == 1  # one 10-day sensor covers the run
    m = 7 * 1440
    truth = np.full(m, 120.0)
    truth[2 * 1440 : 2 * 1440 + 30] = 50.0  # 30 min below 54 on day 3
    events = pd.DataFrame({"minute": [3000], "kind": ["rescue"], "value": [15.0]})
    fake = SimpleNamespace(truth=truth, events=events, insulin_basal=np.full(m, 0.01),
                           insulin_smb=np.zeros(m), insulin_bolus=np.zeros(m))  # fmt: skip
    days = add.daily(fake)
    assert [d["minutes_below_54"] for d in days] == [0, 30, 0, 0, 0, 0]
    assert days[1]["rescue_g"] == 15.0 and days[0]["insulin_basal_u"] == pytest.approx(14.4)
