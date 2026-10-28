"""D8's additivity check: is the harm of a realistic mixed fault suite the sum of its single
events' harms?

    uv run python experiments/additivity.py run --frozen <DECISIONS.md> [--workers 4]
        [--patients ...] [--out results/additivity]
    uv run python experiments/additivity.py run --clean-only ...   # no fault: allowed any time

Per patient: two continuous 7-day `oref0_tb` closed loops on seed 200, rule1800 with the D17 ISF
map (results/calibration/isf_scale_map.json), paper 1's default mix. One is clean; the other
carries paper 1's realistic mixed suite, `FaultSchedule.for_run(200, ("sim:<patient>", 200,
"default"), 7, FaultConfig(pisa_fraction=0.86, severity=0.5))` with paper 1's default rates. The
seed fixes meals, sensor noise and the person's draws, so the pair differs only by the faults.
Faulted runs need --frozen (rule 11), as in run_grid.py.

Day 1 is warm-up (oref0 starts without history) and is not stored; days 2-7 are stored per run and
per day: minutes on truth below 54, below 70, above 180 and above 250 and risk-weighted minutes
(the analysis package's `glycaemia`, GRI weights), insulin delivered (basal or temp, SMB,
boluses) and rescue carbohydrate; and, from the schedule, the events starting that day (PISA,
dropout, stale, spike) and the drift (sign, end-of-life size and the gain at midday; a drift is
continuous, one per 10-day sensor). One JSON per patient, written atomically; a restart skips
finished patients. The analysis is `ghostlows_aid.analysis.secondary.additivity`.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import time
from pathlib import Path

import numpy as np
from ghostlows.faults.suite import FaultConfig
from ghostlows.sim.engine import PATIENTS
from ghostlows.sim.scenario import DEFAULT

OUT = Path("results/additivity")
SEED, DAYS, MIX = 200, 7, "default"  # D8 / D18
SCALE_MAP = Path("results/calibration/isf_scale_map.json")
CONFIG = FaultConfig(pisa_fraction=0.86, severity=0.5)  # paper 1's realistic mixed suite
EVENT_TYPES = ("pisa", "dropout", "stale", "spike")
GLYC = ("minutes_below_54", "minutes_below_70", "minutes_above_180", "minutes_above_250",
        "risk_weighted_min")  # fmt: skip


def schedule(patient: str):
    from ghostlows_aid.faults_online import FaultSchedule

    return FaultSchedule.for_run(SEED, (f"sim:{patient}", SEED, MIX), DAYS, CONFIG)


def daily(result) -> list[dict]:
    """Per day 2..7 of a LoopResult: glycaemia on the 1-min truth, insulin and rescue."""
    from ghostlows_aid.analysis.outcomes import glycaemia

    ev = result.events
    rows = []
    for d in range(2, DAYS + 1):
        lo, hi = (d - 1) * 1440, d * 1440
        g = glycaemia(result.truth[lo:hi])
        row = {"day": d, **{k: g[k] for k in GLYC}}
        row["insulin_basal_u"] = float(result.insulin_basal[lo:hi].sum())
        row["insulin_smb_u"] = float(result.insulin_smb[lo:hi].sum())
        row["insulin_bolus_u"] = float(result.insulin_bolus[lo:hi].sum())
        rescue = ev[(ev["kind"] == "rescue") & (ev["minute"] >= lo) & (ev["minute"] < hi)]
        row["rescue_g"] = float(rescue["value"].sum())
        rows.append(row)
    return rows


def events_per_day(sch) -> list[dict]:
    """Events starting on each day 2..7, and the drift of that day."""
    rows = []
    for d in range(2, DAYS + 1):
        lo, hi = (d - 1) * 1440, d * 1440
        row = {"day": d}
        for kind in EVENT_TYPES:
            row[kind] = sum(
                1 for i in sch.instances if i["fault_type"] == kind and lo <= i["start"] < hi
            )
        mid = (lo + hi) / 2
        drift = [x for x in sch.drift if x[0] <= mid < x[1]]
        size = drift[0][2] if drift else 0.0
        row["drift_sign"] = float(np.sign(size))
        row["drift_size"] = float(size)  # end-of-life multiplicative size, signed
        row["drift_gain_midday"] = (
            float(size * (mid - drift[0][0]) / (CONFIG.sensor_days * 1440)) if drift else 0.0
        )
        rows.append(row)
    return rows


def run_one(job) -> str:
    from ghostlows_aid.loop import LoopConfig, make_loop
    from ghostlows_aid.profile import ProfileSettings

    out, patient, scale, faulted = job
    path = out / f"{patient.replace('#', '')}.json"
    if path.exists():
        return "skip"
    t = time.perf_counter()
    settings = ProfileSettings(isf_rule="rule1800", isf_scale=scale)
    config = LoopConfig(patient, controller="oref0_tb", profile=settings)
    sch = schedule(patient)
    record = {"patient": patient, "seed": SEED, "isf_scale": scale, "config": "pisa_fraction=0.86,"
              " severity=0.5", "events": events_per_day(sch)}  # fmt: skip
    runs = {"clean": None, "faulted": sch} if faulted else {"clean": None}
    for name, faults in runs.items():
        r = make_loop(config, SEED, DAYS, mix=DEFAULT, faults=faults).run().result()
        record[name] = {"valid": bool(r.valid), "errors": int(r.controller_errors),
                        "days": daily(r)}  # fmt: skip
    record["seconds"] = time.perf_counter() - t
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(record))
    os.replace(tmp, path)
    return "done"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run"])
    ap.add_argument("--patients", nargs="+", default=PATIENTS)
    ap.add_argument("--frozen", help="the DECISIONS.md that records the freeze")
    ap.add_argument("--clean-only", action="store_true", help="clean runs only (no fault)")
    ap.add_argument("--scale-map", type=Path, default=SCALE_MAP)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args(argv)
    faulted = not a.clean_only
    if faulted and not (a.frozen and Path(a.frozen).is_file()):
        raise SystemExit("faulted runs need --frozen <the DECISIONS.md that records the freeze>")
    scales = json.loads(a.scale_map.read_text())["isf_scale"]
    out = a.out / ("clean_only" if a.clean_only else "")
    out.mkdir(parents=True, exist_ok=True)
    jobs = [(out, p, float(scales.get(p, 1.0)), faulted) for p in a.patients]
    if a.workers == 1:
        done = [run_one(j) for j in jobs]
    else:
        with mp.get_context("fork").Pool(a.workers) as pool:
            done = list(pool.imap_unordered(run_one, jobs))
    print(f"{done.count('done')} patients run, {done.count('skip')} skipped, in {out}")


if __name__ == "__main__":
    main()
