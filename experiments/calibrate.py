"""Settings calibration (PI_decisions §4, D-settings): clean closed loops on calibration seeds.

    uv run python experiments/calibrate.py run [--seeds 90001 90002 90003] [--workers 4]
        [--rules quest rule1800] [--modes oref0_tb oref1_smb] [--scales 1.0]
    uv run python experiments/calibrate.py table

Every (patient, ISF rule, ISF scale, mode, seed) is one 3-day run of paper 1's default mix with no
fault, scored on the 1-min plasma truth over all 3 days. Each run writes one JSON file under
results/calibration/units/ as soon as it finishes; a restart skips finished runs. `table` pools
each patient's runs over the seeds (equal lengths, so the mean of the per-run percentages) and
applies the pre-specified selection rule on `oref0_tb` at scale 1.0. No fault code runs here.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
from ghostlows.sim.engine import PATIENTS
from ghostlows.sim.scenario import DEFAULT

OUT = Path("results/calibration")
SEEDS = (90001, 90002, 90003)  # paper 2's calibration seeds (PI_decisions §7: >= 90000)
STAGES = ("ms_iob", "ms_meal", "ms_determine_basal", "ms_bridge")
# PI_decisions §4 acceptance criteria [Battelino 2019], on truth
TIR_MIN, TBR70_MAX, TBR54_MAX = 70.0, 4.0, 1.0


def unit_path(out: Path, patient, rule, scale, mode, seed) -> Path:
    return out / "units" / f"{patient.replace('#', '')}_{rule}_x{scale:g}_{mode}_{seed}.json"


def run_one(job) -> str:
    from ghostlows_aid.loop import LoopConfig, make_loop
    from ghostlows_aid.metrics import summarise
    from ghostlows_aid.oref import REASON_KINDS, reason_kind
    from ghostlows_aid.profile import ProfileSettings

    out, patient, rule, scale, mode, seed, days = job
    path = unit_path(out, patient, rule, scale, mode, seed)
    if path.exists():
        return "skip"
    row = {"patient": patient, "rule": rule, "scale": scale, "mode": mode, "seed": seed}
    t = time.perf_counter()
    try:
        settings = ProfileSettings(isf_rule=rule, isf_scale=scale)
        loop = make_loop(LoopConfig(patient, controller=mode, profile=settings), seed, days,
                         mix=DEFAULT)  # fmt: skip
        r = loop.run().result()
        s = r.samples
        kinds = s["rT_reason"].map(reason_kind).value_counts()
        row.update(summarise(r))
        row["reasons"] = {k: int(kinds.get(k, 0)) for k in REASON_KINDS}
        row["max_iob_reasons"] = int(s["rT_reason"].str.contains("max_iob", na=False).sum())
        row["smb_u_per_day"] = float(r.insulin_smb.sum() / days)
        row["numbers"] = {k: float(v) for k, v in loop.numbers.items()}
        row["timing"] = {
            c: {"median": float(s[c].median()), "p95": float(s[c].quantile(0.95)),
                "sum_s": float(s[c].sum() / 1000)}
            for c in STAGES if c in s
        }  # fmt: skip
        row["calls"] = int(s["ms_bridge"].notna().sum())
        row["invalid"] = (not r.valid) or r.controller_errors > 0
    except Exception as e:  # an invalid run, recorded and counted
        row.update(invalid=True, error=repr(e))
    row["seconds"] = time.perf_counter() - t
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(row))
    os.replace(tmp, path)  # atomic: a killed run leaves no partial unit
    return "done"


def run(a) -> None:
    chosen = json.loads(Path(a.scale_map).read_text())["isf_scale"] if a.scale_map else {}
    jobs = [
        (a.out, p, rule, scale, mode, seed, a.days)
        for seed in a.seeds
        for p in a.patients
        for rule in a.rules
        for scale in ([chosen[p]] if chosen else a.scales)
        for mode in a.modes
    ]
    t = time.perf_counter()
    with mp.get_context("fork").Pool(a.workers) as pool:  # the parent never starts a bridge
        done = list(pool.imap_unordered(run_one, jobs))
    wall = time.perf_counter() - t
    (a.out / "run_log.jsonl").open("a").write(
        json.dumps({"jobs": len(jobs), "ran": done.count("done"), "workers": a.workers,
                    "wall_s": wall}) + "\n"
    )  # fmt: skip
    print(f"{done.count('done')} runs, {done.count('skip')} skipped, {wall:.0f} s wall")


def load(out: Path) -> pd.DataFrame:
    rows = [json.loads(p.read_text()) for p in sorted((out / "units").glob("*.json"))]
    df = pd.DataFrame(rows)
    for k in ("stale", "temp_mismatch", "flat", "low_suspend"):
        df[k] = [r.get(k, 0) if isinstance(r, dict) else 0 for r in df["reasons"]]
    df["isf"] = [n["isf"] if isinstance(n, dict) else np.nan for n in df["numbers"]]
    df["max_iob"] = [n["max_iob"] if isinstance(n, dict) else np.nan for n in df["numbers"]]
    return df


METRICS = ["tir_pct", "tbr70_pct", "tbr54_pct", "tar180_pct", "tar250_pct", "mean",
           "insulin_u_per_day", "smb_u_per_day", "rescue_g_per_day"]  # fmt: skip


def per_patient(df: pd.DataFrame) -> pd.DataFrame:
    g = df.groupby(["rule", "scale", "mode", "patient"])
    t = g[METRICS].mean()
    t["runs"] = g.size()
    t["invalid"] = g["invalid"].sum().astype(int)
    for k in ("stale", "temp_mismatch", "max_iob_reasons"):
        t[k] = g[k].sum().astype(int)
    t["isf"], t["max_iob"] = g["isf"].first(), g["max_iob"].first()
    t["passes"] = (
        (t["tir_pct"] >= TIR_MIN)
        & (t["tbr70_pct"] <= TBR70_MAX)
        & (t["tbr54_pct"] <= TBR54_MAX)
        & (t["stale"] == 0)
        & (t["temp_mismatch"] == 0)
        & (t["invalid"] == 0)
    )
    return t.reset_index()


def cohort(t: pd.DataFrame) -> pd.DataFrame:
    def rng(x):
        return f"{x.median():.1f} [{x.min():.1f}, {x.max():.1f}]"

    rows = []
    for (rule, scale, mode), d in t.groupby(["rule", "scale", "mode"]):
        rows.append({
            "rule": rule, "scale": scale, "mode": mode, "patients": len(d),
            "pass": int(d["passes"].sum()),
            "TIR median [range]": rng(d["tir_pct"]), "TBR<70 median [range]": rng(d["tbr70_pct"]),
            "TBR<54 median [range]": rng(d["tbr54_pct"]), "TAR>180 median": round(
                d["tar180_pct"].median(), 1), "TAR>250 median": round(d["tar250_pct"].median(), 1),
            "mean BG median": round(d["mean"].median(), 1),
            "insulin U/day median": round(d["insulin_u_per_day"].median(), 1),
            "rescue g/day median": round(d["rescue_g_per_day"].median(), 1),
            "invalid runs": int(d["invalid"].sum()), "stale+mismatch": int(
                (d["stale"] + d["temp_mismatch"]).sum()),
            "max_iob reasons": int(d["max_iob_reasons"].sum()),
        })  # fmt: skip
    return pd.DataFrame(rows)


def select(t: pd.DataFrame) -> dict:
    """PI_decisions §4: on oref0_tb at scale 1, the rule under which more patients pass; ties go to
    the larger ISF (less aggressive), read as the rule with the larger ISF for most patients."""
    d = t[(t["mode"] == "oref0_tb") & (t["scale"] == 1.0)].set_index(["rule", "patient"])
    n = {rule: int(d.loc[rule, "passes"].sum()) for rule in ("quest", "rule1800")}
    larger = (d.loc["rule1800", "isf"] > d.loc["quest", "isf"]).mean()
    larger_rule = "rule1800" if larger > 0.5 else "quest"
    chosen = max(n, key=lambda r: (n[r], r == larger_rule))
    return {"passing": n, "share_rule1800_isf_larger": float(larger), "chosen": chosen,
            "failing_under_chosen": int(30 - n[chosen]), "scale_step": bool(30 - n[chosen] > 5),
            "failing_patients": sorted(d.loc[chosen].index[~d.loc[chosen, "passes"]])}  # fmt: skip


def md(df: pd.DataFrame) -> str:
    """A Markdown table (tabulate is not installed)."""
    head = "| " + " | ".join(map(str, df.columns)) + " |\n|" + "---|" * len(df.columns)
    body = ["| " + " | ".join(str(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join([head, *body])


def scale_step(t: pd.DataFrame, sel: dict) -> pd.DataFrame:
    """PI_decisions §4: per patient, the smallest ISF scale in {0.8, 1.0, 1.25, 1.5} that meets the
    criteria on oref0_tb, read literally (every patient, not only the failing ones); a patient
    meeting them at no scale keeps 1.0 and is listed. Returns the runs at the chosen scales."""
    d = t[(t["rule"] == sel["chosen"]) & (t["mode"] == "oref0_tb")]
    sel["scales_evaluated"] = sorted(float(x) for x in d["scale"].unique())
    scales, none = {}, []
    for p, g in d.groupby("patient"):
        ok = g[g["passes"]].sort_values("scale")
        scales[p] = float(ok["scale"].iloc[0]) if len(ok) else 1.0
        if not len(ok):
            none.append(p)
    sel["isf_scale"] = scales
    sel["no_passing_scale"] = none
    sel["changed_from_1"] = sorted(p for p, x in scales.items() if x != 1.0)
    sel["passing_after_scale"] = int(len(scales) - len(none))
    m = t[t["rule"] == sel["chosen"]]
    final = m[[x == scales.get(p) for p, x in zip(m["patient"], m["scale"], strict=True)]].copy()
    final["scale"] = -1.0  # marks "per patient" in the cohort summary
    return final


def table(a) -> None:
    df = load(a.out)
    t = per_patient(df)
    t.to_csv(a.out / "per_patient.csv", index=False)
    c = cohort(t)
    sel = select(t)
    final = scale_step(t, sel) if sel["scale_step"] else None
    (a.out / "selection.json").write_text(json.dumps(sel, indent=2))
    with (a.out / "summary.md").open("w") as f:
        f.write(md(c) + "\n\n")
        if final is not None:
            f.write("\n**Chosen settings (per-patient ISF scale)**\n\n" + md(cohort(final)) + "\n")
            print(cohort(final).to_string(index=False))
        cols = ["patient", "tir_pct", "tbr70_pct", "tbr54_pct", "tar180_pct", "tar250_pct", "mean",
                "insulin_u_per_day", "rescue_g_per_day", "invalid", "passes"]  # fmt: skip
        for (rule, scale, mode), d in t.groupby(["rule", "scale", "mode"]):
            f.write(f"\n**{rule} x{scale:g}, {mode}**\n\n")
            f.write(md(d[cols].round(2)) + "\n")
    print(c.to_string(index=False))
    print(json.dumps(sel))
    tm = df.dropna(subset=["seconds"])
    print("seconds per 3-day run by mode (4 workers):",
          tm.groupby("mode")["seconds"].median().round(1).to_dict())  # fmt: skip


def main(argv=None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "table"])
    ap.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS))
    ap.add_argument("--patients", nargs="+", default=PATIENTS)
    ap.add_argument("--rules", nargs="+", default=["quest", "rule1800"])
    ap.add_argument("--modes", nargs="+", default=["oref0_tb", "oref1_smb"])
    ap.add_argument("--scales", type=float, nargs="+", default=[1.0])
    ap.add_argument("--scale-map", help="selection.json: run each patient at its chosen scale")
    ap.add_argument("--days", type=int, default=3)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args(argv)
    (run if a.cmd == "run" else table)(a)


if __name__ == "__main__":
    main()
