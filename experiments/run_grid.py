"""The main grid (PI_decisions §7): one unit per (controller, patient, seed), written as soon as it
finishes; a restart skips finished units.

    uv run python experiments/run_grid.py run --n 20 [--max-tier 3] [--controllers ...]
        [--patients ...] [--seeds ...] [--workers 4] [--out results/grid]   # timing mode
    uv run python experiments/run_grid.py run ... --faults --frozen DECISIONS.md  # after the freeze
    uv run python experiments/run_grid.py cost [--out results/grid] [--budget 24]

Timing mode (the default) runs every arm with its fault switched off (variants kept) and records
the cost of the warm-up, the trunk and each arm; it computes no outcome. With --faults, allowed
only once the protocol is frozen (AGENTS.md rule 11), each slot's arms become ArmTraces, the
analysis package's contract (ghostlows_aid.analysis.outcomes), and its slot_outcomes() turns them
into outcome rows: the harness computes no outcome itself. Raw 5-min traces are kept for a fixed
random 10% of the (patient, seed) pairs, for every controller (episodes.arm_trace).
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from ghostlows.sim.engine import PATIENTS

from ghostlows_aid.episodes import TRUNK_MIN, arm_trace, arms_for, run_unit
from ghostlows_aid.loop import CONTROLLERS
from ghostlows_aid.profile import ProfileSettings

OUT = Path("results/grid")
MAIN_SEED0 = 101  # main seeds 101, 102, ... (PI_decisions §7)
TRACE_SEED, TRACE_FRACTION = 20261028, 0.10  # fixed before any run
# PI_decisions §4: the ISF rule chosen by the pre-specified selection on calibration seeds
# (experiments/calibrate.py; SPIKE doc §7; DECISIONS.md D16). The per-patient ISF scale map is
# the failing-patients reading of the scale step (DECISIONS.md D17).
ISF_RULE = "rule1800"
ISF_SCALE_MAP = "results/calibration/isf_scale_map.json"


def keeps_traces(patient: str, seed: int) -> bool:
    """A fixed random 10% of (patient, seed) pairs, independent of N and of the patient list."""
    return np.random.default_rng([TRACE_SEED, PATIENTS.index(patient), seed]).random() < 0.10


def write_atomic(path: Path, write) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    write(tmp)
    os.replace(tmp, path)


def run_job(job) -> str:
    out, controller, patient, seed, k, n, max_tier, faults, isf_rule, isf_scale = job
    name = f"{controller}_{patient.replace('#', '')}_{seed}"
    marker = out / "units" / f"{name}.json"
    if marker.exists():
        return "skip"
    if faults:
        from ghostlows_aid.analysis import outcomes, schema
    arms = [a for a in arms_for(controller) if a.tier <= max_tier]
    t = time.perf_counter()
    settings = ProfileSettings(isf_rule=isf_rule, isf_scale=isf_scale)
    unit, runs = run_unit(controller, patient, seed, k, n, settings=settings, arms=arms,
                          with_faults=faults)  # fmt: skip
    record = {"controller": controller, "patient": patient, "seed": seed, "k": k, "n": n,
              "faults": faults, "isf_rule": isf_rule, "isf_scale": isf_scale,
              "max_tier": max_tier, "onsets": unit.onsets,
              "warmup_s": unit.warmup_seconds, "trunk_s": unit.trunk_seconds,
              "arm_s": {}}  # fmt: skip
    s = unit.trunk.samples
    record["trunk_stages_ms"] = {
        c: {"median": float(s[c].median()), "p95": float(s[c].quantile(0.95))}
        for c in ("ms_iob", "ms_meal", "ms_determine_basal", "ms_bridge") if c in s
    }  # fmt: skip
    by_slot, kept = defaultdict(list), []
    for run in runs:
        if run.arm.group != "clean":
            record["arm_s"][run.arm.name] = run.seconds
        if faults:
            by_slot[run.arm.slot].append(outcomes.trace_from_run(unit, run))
        if keeps_traces(patient, seed):
            kept.append(arm_trace(unit, run))
    if faults:
        rows = schema.validate(pd.concat([outcomes.slot_outcomes(ts) for ts in by_slot.values()]))
        write_atomic(out / "outcomes" / f"{name}.parquet", lambda p: rows.to_parquet(p))
    if kept:
        write_atomic(out / "traces" / f"{name}.pkl.gz", lambda p: pd.to_pickle(kept, p, "gzip"))
    record["total_s"] = time.perf_counter() - t
    write_atomic(marker, lambda p: p.write_text(json.dumps(record)))  # last: marks it finished
    return "done"


def run(a) -> None:
    if a.faults and not (a.frozen and Path(a.frozen).is_file()):
        raise SystemExit("--faults needs --frozen <the DECISIONS.md that records the freeze>")
    seeds = a.seeds or [MAIN_SEED0 + k for k in range(a.n)]
    scales = json.loads(Path(a.isf_scale_map).read_text())["isf_scale"] if a.isf_scale_map else {}
    jobs = [
        (a.out, c, p, s, seeds.index(s) if not a.seeds else a.k, a.n, a.max_tier, a.faults,
         a.isf_rule, scales.get(p, 1.0))
        for s in seeds for p in a.patients for c in a.controllers
    ]  # fmt: skip
    t = time.perf_counter()
    if a.workers == 1:
        done = [run_job(j) for j in jobs]
    else:
        with mp.get_context("fork").Pool(a.workers) as pool:  # the parent starts no bridge
            done = list(pool.imap_unordered(run_job, jobs))
    wall = time.perf_counter() - t
    log = {"jobs": len(jobs), "ran": done.count("done"), "workers": a.workers, "wall_s": wall}
    with (a.out / "run_log.jsonl").open("a") as f:
        f.write(json.dumps(log) + "\n")
    print(json.dumps(log))


def unit_costs(out: Path) -> pd.DataFrame:
    rows = [json.loads(p.read_text()) for p in sorted((out / "units").glob("*.json"))]
    return pd.DataFrame(rows)


# PI_decisions §7: N = 20, else 15, else 10; then drop tier 3, then the tier-2 arms of oref1_smb,
# then those of open_loop; tier 1 always runs at N >= 10; patients are never reduced.
STEPS = [
    (20, frozenset()),
    (15, frozenset()),
    (10, frozenset()),
    (10, frozenset({(c, 3) for c in CONTROLLERS})),
    (10, frozenset({(c, 3) for c in CONTROLLERS} | {("oref1_smb", 2)})),
    (10, frozenset({(c, 3) for c in CONTROLLERS} | {("oref1_smb", 2), ("open_loop", 2)})),
]


def projection(df: pd.DataFrame, budget_h: float, patients: int = 30) -> tuple[list, dict]:
    """CPU-hours of each step of the scaling rule from the measured median trunk cost and the
    median cost of each arm, per controller."""
    trunk = df.groupby("controller")["trunk_s"].median().to_dict()
    arm_cost = {}
    for c, d in df.groupby("controller"):
        per_arm = pd.DataFrame(list(d["arm_s"])).median().to_dict()
        arm_cost[c] = {a.name: per_arm.get(a.name, np.nan) for a in arms_for(c)}
    rows = []
    for n, dropped in STEPS:
        total = 0.0
        for c in CONTROLLERS:
            arms = [a for a in arms_for(c) if (c, a.tier) not in dropped]
            total += patients * n * (trunk[c] + sum(arm_cost[c][a.name] for a in arms))
        label = "full design" if not dropped else "drop " + ", ".join(
            sorted(f"{c} T{t}" for c, t in dropped))  # fmt: skip
        cpu_h = total / 3600
        rows.append({"N": n, "arms": label, "cpu_h": cpu_h, "fits": cpu_h <= budget_h})
    chosen = next((r for r in rows if r["fits"]), None)
    medians = {c: float(np.nanmedian(list(v.values()))) for c, v in arm_cost.items()}
    return rows, {"trunk_s": trunk, "arm_s_median": medians, "chosen": chosen}


def cost(a) -> None:
    df = unit_costs(a.out)
    rows, info = projection(df, a.budget)
    days = TRUNK_MIN / 1440
    per = df.groupby("controller").agg(
        units=("trunk_s", "size"), warmup_s=("warmup_s", "median"), trunk_s=("trunk_s", "median"),
        total_s=("total_s", "median"))  # fmt: skip
    per["s_per_patient_day"] = per["trunk_s"] / days
    per["arm_s"] = pd.Series(info["arm_s_median"])
    print(per.round(2).to_string())
    print(pd.DataFrame(rows).round(2).to_string(index=False))
    print(json.dumps(info["chosen"]))
    (a.out / "cost.json").write_text(
        json.dumps(
            {"per_controller": per.round(3).to_dict("index"), "steps": rows, **info},
            indent=2,
            default=float,
        )
    )


def main(argv=None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "cost"])
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--seeds", type=int, nargs="+", help="explicit seeds (timing samples)")
    ap.add_argument("--k", type=int, default=0, help="onset stratum for explicit seeds")
    ap.add_argument("--patients", nargs="+", default=PATIENTS)
    ap.add_argument("--controllers", nargs="+", default=list(CONTROLLERS))
    ap.add_argument("--max-tier", type=int, default=3, help="0: the trunk only")
    ap.add_argument("--isf-rule", default=ISF_RULE)
    ap.add_argument("--isf-scale-map", default=ISF_SCALE_MAP, help="per-patient ISF scales (JSON)")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--faults", action="store_true")
    ap.add_argument("--frozen")
    ap.add_argument("--budget", type=float, default=24.0)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args(argv)
    (run if a.cmd == "run" else cost)(a)


if __name__ == "__main__":
    main()
