"""TECHNICAL_DESIGN §2.1 levers: are oref0's outputs identical with shorter histories, or with
lib/meal skipped when no carbs were entered in the last 6 h? And what does each save?

    uv run python experiments/levers.py [--patients adult#001 adolescent#003 child#005]
        [--seed 90002] [--days 2] [--every 3] [--out results/levers]

Clean closed-loop runs (default mix, calibration seed) record every controller call; each recorded
state is then sent again with the full 24-h histories and with each lever, and the replies (glucose
status, IOB, meal data and the whole rT) are compared exactly. Gapped states are made by deleting a
block of readings from recorded states (no fault code). A lever may be adopted only if every reply
is identical.
"""

from __future__ import annotations

import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np
from ghostlows.sim.scenario import DEFAULT

from ghostlows_aid.loop import LoopConfig, make_loop
from ghostlows_aid.oref import shared_bridge
from ghostlows_aid.profile import ProfileSettings

TRIMS = [(6, 12), (8, 12), (8, 14), (12, 12), (12, 16)]  # (glucose h, treatment h)
SIX_H = 6 * 3600 * 1000


def record(patient: str, controller: str, seed: int, days: int, isf_rule: str) -> tuple[list, dict]:
    config = LoopConfig(patient, controller=controller, profile=ProfileSettings(isf_rule=isf_rule))
    loop = make_loop(config, seed, days, mix=DEFAULT)
    states, real = [], loop.controller.step

    def spy(clock, glucose, treatments, currenttemp):
        states.append((clock, *copy.deepcopy((glucose, treatments, currenttemp))))
        return real(clock, glucose, treatments, currenttemp)

    loop.controller.step = spy
    loop.run()
    return states, loop.controller


def request(ctrl, clock, glucose, treatments, currenttemp, skip_meal=False) -> dict:
    return {
        "op": "step",
        "clock": int(clock),
        "glucose": glucose,
        "treatments": treatments,
        "carbs": [],
        "profile": ctrl.profile,
        "currenttemp": currenttemp,
        "autosens": ctrl.autosens,
        "microBolusAllowed": ctrl.microbolus,
        "skip_meal": skip_meal,
    }


def trim(clock, glucose, treatments, g_h, t_h):
    from datetime import datetime

    def ms(t):
        return datetime.fromisoformat(t["created_at"].replace("Z", "+00:00")).timestamp() * 1000

    g = [e for e in glucose if e["date"] >= clock - g_h * 3600 * 1000]
    t = [x for x in treatments if ms(x) >= clock - t_h * 3600 * 1000]
    return g, t


def carbs_recent(clock, treatments) -> bool:
    from datetime import datetime

    for x in treatments:
        if x.get("carbs", 0) >= 1:
            t = datetime.fromisoformat(x["created_at"].replace("Z", "+00:00")).timestamp() * 1000
            if t >= clock - SIX_H:
                return True
    return False


def outputs(reply: dict) -> str:
    keep = {k: reply[k] for k in ("glucose_status", "iob", "meal_data", "rT")}
    return json.dumps(keep, sort_keys=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=["adult#001", "adolescent#003", "child#005"])
    ap.add_argument("--controllers", nargs="+", default=["oref0_tb", "oref1_smb"])
    ap.add_argument("--seed", type=int, default=90002)  # a calibration seed
    ap.add_argument("--isf", default="rule1800")
    ap.add_argument("--days", type=int, default=2)
    ap.add_argument("--every", type=int, default=3)
    ap.add_argument("--gapped", type=int, default=150)
    ap.add_argument("--out", default="results/levers")
    a = ap.parse_args(argv)
    bridge = shared_bridge()
    rng = np.random.default_rng(0)
    variants = ["full", *[f"trim_g{g}_t{t}" for g, t in TRIMS], "skip_meal"]
    stats = {c: {v: {"n": 0, "same": 0, "ms": []} for v in variants} for c in a.controllers}
    stage_ms = {c: {"glucose_status": [], "iob": [], "meal": [], "determine_basal": [],
                    "bridge_total": [], "round_trip": []} for c in a.controllers}  # fmt: skip
    for controller in a.controllers:
        for patient in a.patients:
            states, ctrl = record(patient, controller, a.seed, a.days, a.isf)
            picked = states[:: a.every]
            gapped = [states[i] for i in rng.choice(len(states), a.gapped, replace=False)]
            for kind, batch in (("clean", picked), ("gapped", gapped)):
                for clock, glucose, treatments, currenttemp in batch:
                    if kind == "gapped":  # delete a 20-120 min block of readings in the last 6 h
                        start = clock - int(rng.uniform(0, 6 * 60)) * 60000
                        length = int(rng.uniform(20, 120)) * 60000
                        glucose = [e for e in glucose if not start - length < e["date"] <= start]
                        if not glucose:
                            continue
                    t = time.perf_counter()
                    ref = bridge.request(request(ctrl, clock, glucose, treatments, currenttemp))
                    rt = (time.perf_counter() - t) * 1000
                    base = outputs(ref)
                    s = stats[controller]["full"]
                    s["n"] += 1
                    s["same"] += 1
                    s["ms"].append(ref["total_ms"])
                    for stage in ("glucose_status", "iob", "meal", "determine_basal"):
                        stage_ms[controller][stage].append(ref["timing_ms"][stage])
                    stage_ms[controller]["bridge_total"].append(ref["total_ms"])
                    stage_ms[controller]["round_trip"].append(rt)
                    for g_h, t_h in TRIMS:
                        g, tr = trim(clock, glucose, treatments, g_h, t_h)
                        rep = bridge.request(request(ctrl, clock, g, tr, currenttemp))
                        s = stats[controller][f"trim_g{g_h}_t{t_h}"]
                        s["n"] += 1
                        s["same"] += outputs(rep) == base
                        s["ms"].append(rep["total_ms"])
                    if not carbs_recent(clock, treatments):
                        rep = bridge.request(
                            request(ctrl, clock, glucose, treatments, currenttemp, skip_meal=True)
                        )
                        s = stats[controller]["skip_meal"]
                        s["n"] += 1
                        same_rt = json.dumps(rep["rT"], sort_keys=True) == json.dumps(
                            ref["rT"], sort_keys=True
                        )
                        s["same"] += same_rt  # meal_data differs by construction; rT must not
                        s["ms"].append(rep["total_ms"])
            print(f"{controller} {patient}: {len(states)} states", flush=True)
    summary = {"seed": a.seed, "isf_rule": a.isf, "days": a.days, "patients": a.patients}
    for controller in a.controllers:
        summary[controller] = {
            v: {
                "states": d["n"],
                "identical": d["same"],
                "median_ms": float(np.median(d["ms"])) if d["ms"] else None,
                "p95_ms": float(np.percentile(d["ms"], 95)) if d["ms"] else None,
            }
            for v, d in stats[controller].items()
        }
        summary[controller]["stages_full"] = {
            k: {"median_ms": float(np.median(v)), "p95_ms": float(np.percentile(v, 95))}
            for k, v in stage_ms[controller].items()
        }
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "levers.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
