"""Task 6a: adult#001, 3 days, paper 1's standard scenario, oref0_tb, a calibration seed.

    uv run python experiments/acceptance_adult001.py [--seed 90001] [--isf quest]
        [--controller oref0_tb] [--out results/acceptance]

Acceptance: truth TIR 70-180 >= 70%, TBR < 70 <= 4%, no temp-mismatch or stale-data reasons.
Clean runs only (no faults).
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from ghostlows_aid.loop import LoopConfig, make_loop
from ghostlows_aid.metrics import summarise
from ghostlows_aid.oref import reason_kind
from ghostlows_aid.profile import ProfileSettings


def main(argv=None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=90001)  # a calibration seed (>= 90000)
    ap.add_argument("--isf", default="quest")
    ap.add_argument("--controller", default="oref0_tb")
    ap.add_argument("--days", type=int, default=3)
    ap.add_argument("--out", default="results/acceptance")
    a = ap.parse_args(argv)
    cfg = LoopConfig("adult#001", controller=a.controller, profile=ProfileSettings(isf_rule=a.isf))
    t0 = time.time()
    loop = make_loop(cfg, a.seed, a.days).run()
    r = loop.result()
    kinds = r.samples["rT_reason"].map(reason_kind).value_counts().to_dict()
    out = {
        **summarise(r),
        "reason_kinds": kinds,
        "seconds": round(time.time() - t0, 1),
        "seed": a.seed,
        "isf_rule": a.isf,
        "controller": a.controller,
    }
    ok = (
        out["tir_pct"] >= 70
        and out["tbr70_pct"] <= 4
        and kinds.get("temp_mismatch", 0) == 0
        and kinds.get("stale", 0) == 0
    )
    out["accepted"] = bool(ok)
    path = Path(a.out)
    path.mkdir(parents=True, exist_ok=True)
    r.samples.to_parquet(path / f"samples_{a.controller}_{a.isf}_s{a.seed}.parquet")
    (path / f"summary_{a.controller}_{a.isf}_s{a.seed}.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
