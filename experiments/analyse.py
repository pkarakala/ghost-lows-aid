"""Run the whole analysis of paper 2 on an outcome table: CSV tables, a Markdown summary, figures.

    python experiments/analyse.py --synthetic                         # generated synthetic data
    python experiments/analyse.py --input results/grid --protocol-tag protocol-v1

Input: parquet files (or directories of them) with the per-arm outcome table of
`ghostlows_aid.analysis.schema`. Faulted arms of simulation output are analysed only with
`--protocol-tag protocol-v1` when that tag exists in this repository and the code being run
contains it (memo §15: no faulted-versus-clean outcome before the freeze). A clean-only
simulation table is allowed and gets the descriptive tables only. `--synthetic` generates the
default synthetic design, or, with `--input`, accepts only rows marked synthetic.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from ghostlows_aid.analysis import cross_layer, figures, hypotheses, schema, secondary
from ghostlows_aid.analysis import effects as ef
from ghostlows_aid.analysis.synthetic import Design, synthetic_outcomes

REPO = Path(__file__).resolve().parents[1]
FREEZE_TAG = "protocol-v1"  # memo §15
REFUSED = 2  # exit status when the freeze guard refuses


class Refused(SystemExit):
    def __init__(self, message: str):
        print(f"analyse.py: refused: {message}", file=sys.stderr)
        super().__init__(REFUSED)


def _git(repo: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True)


def protocol_tag_ok(tag: str | None, repo: Path = REPO) -> tuple[bool, str]:
    """(ok, reason): the tag is protocol-v1, exists in the repository and is an ancestor of (or
    equal to) the commit being run."""
    if tag != FREEZE_TAG:
        return False, f"faulted arms of simulation output need --protocol-tag {FREEZE_TAG}"
    if _git(repo, "rev-parse", "--verify", "--quiet", f"refs/tags/{tag}^{{commit}}").returncode:
        return False, f"tag {tag} does not exist in {repo}"
    if _git(repo, "merge-base", "--is-ancestor", f"refs/tags/{tag}", "HEAD").returncode:
        return False, f"the code being run (HEAD) does not contain tag {tag}"
    return True, ""


def read_tables(paths: list[str]) -> pd.DataFrame:
    files = []
    for p in map(Path, paths):
        files += sorted(p.rglob("*.parquet")) if p.is_dir() else [p]
    if not files:
        raise SystemExit("analyse.py: no parquet files in the input")
    return pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)


def _commit() -> str:
    r = _git(REPO, "rev-parse", "HEAD")
    return r.stdout.strip() if r.returncode == 0 else ""


def _md_table(df: pd.DataFrame, cols: list[str], digits: int = 3) -> str:
    d = df[[c for c in cols if c in df]].copy()
    for c in d:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else f"{v:.{digits}g}")
    lines = ["| " + " | ".join(d.columns) + " |", "|" + "---|" * len(d.columns)]
    lines += ["| " + " | ".join(str(v) for v in row) + " |" for row in d.itertuples(index=False)]
    return "\n".join(lines)


def describe(table: pd.DataFrame, counts: pd.DataFrame) -> str:
    units = table[list(schema.UNIT_KEYS)].drop_duplicates()
    seeds = table.groupby(["controller", "patient_id"])["seed"].nunique()
    lines = [
        f"- rows {len(table)}, units (controller, patient, seed, slot) {len(units)}",
        f"- controllers {', '.join(sorted(table['controller'].unique()))}",
        f"- patients {table['patient_id'].nunique()}, seeds per patient {seeds.min()} to "
        f"{seeds.max()}, arms {table['arm'].nunique()}",
        f"- data source {', '.join(sorted(table['data_source'].unique()))}",
        f"- invalid rule: {int(counts['dropped'].sum())} rows dropped "
        f"({int(counts['invalid_arms'].sum())} invalid arms)",
    ]
    return "\n".join(lines)


def clean_levels(table: pd.DataFrame, boot: ef.Bootstrap) -> pd.DataFrame:
    """Descriptive: the clean arm's outcome levels per controller and slot (no contrast)."""
    clean = table[table["arm"] == schema.CLEAN]
    cols = [c for c, col in schema.COLUMNS.items() if col.effect == "diff"]
    return ef.estimate(clean, cols, ef.BY_ARM, boot)


STANDARD_MIX_NOTE = (
    "The standard scenario mix was not run: D8 allowed it only if the full design at N = 20 fitted"
    " the budget with its cost included, and it did not (D18; BENCHMARK §8.4)."
)
SCALE_MAP = REPO / "results/calibration/isf_scale_map.json"  # D17


def outside_acceptance(path: Path) -> set[str]:
    """The D17 map's patients outside the acceptance criteria, in the table's `sim:` form."""
    names = json.loads(path.read_text())["outside_acceptance"]
    return {n if n.startswith("sim:") else f"sim:{n}" for n in names}


def read_additivity_records(folder: Path, excluded: list[str]) -> list[dict]:
    """experiments/additivity.py's per-patient records, without the excluded patients."""
    records = [json.loads(p.read_text()) for p in sorted(folder.glob("*.json"))]
    drop = set(excluded)
    return [r for r in records if f"sim:{r['patient']}".replace("sim:sim:", "sim:") not in drop]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--input", nargs="*", default=[], help="parquet files or directories")
    ap.add_argument("--synthetic", action="store_true", help="synthetic data only")
    ap.add_argument("--synthetic-seeds", type=int, default=6, help="seeds of the generated design")
    ap.add_argument("--protocol-tag", default=None, help=f"must be {FREEZE_TAG} for faulted arms")
    ap.add_argument("--out", default=None, help="output directory")
    ap.add_argument("--n-boot", type=int, default=ef.N_BOOT)
    ap.add_argument("--n-boot-confirmatory", type=int, default=ef.N_BOOT_CONFIRMATORY)
    ap.add_argument("--boot-seed", type=int, default=ef.BOOT_SEED)
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--primary-input", nargs="*", default=[], help="S1: the primary tables")
    ap.add_argument("--additivity", default=None, help="additivity.py output folder (D8)")
    ap.add_argument(
        "--exclude-outside-acceptance",
        nargs="?",
        const=str(SCALE_MAP),
        default=None,
        help="S2 (D18): drop the map's outside_acceptance patients before every estimate",
    )
    args = ap.parse_args(argv)
    t_start = time.perf_counter()

    if args.input:
        table = read_tables(args.input)
    elif args.synthetic:
        table = synthetic_outcomes(Design(n_seeds=args.synthetic_seeds))
    else:
        ap.error("give --input, or --synthetic to generate synthetic data")
    table = schema.validate(schema.coerce(table))
    excluded: list[str] = []
    if args.exclude_outside_acceptance:
        excluded = sorted(outside_acceptance(Path(args.exclude_outside_acceptance)))
        table = table[~table["patient_id"].isin(excluded)].reset_index(drop=True)
    sources = set(table["data_source"].unique())
    if args.synthetic and sources != {"synthetic"}:
        raise Refused("--synthetic given but the input holds rows not marked synthetic")
    faulted = (table["arm"] != schema.CLEAN) & (table["data_source"] != "synthetic")
    if faulted.any():
        ok, why = protocol_tag_ok(args.protocol_tag)
        if not ok:
            raise Refused(why)
    synthetic = sources == {"synthetic"}
    clean_only = bool((table["arm"] == schema.CLEAN).all())
    default_out = "results/synthetic" if synthetic else "results/analysis"
    out = Path(args.out) if args.out else REPO / default_out
    fig_dir = out / "figures" if (args.out or synthetic) else REPO / "results/figures"
    out.mkdir(parents=True, exist_ok=True)

    kept, counts = ef.apply_invalid_rule(table)
    boot = ef.boot_for(kept, args.n_boot, args.boot_seed)
    results: dict[str, pd.DataFrame] = {"invalid_counts": counts}
    results["clean_levels"] = clean_levels(kept, boot)
    effects = ef.paired_effects(kept)
    if not clean_only:
        boot_c = ef.boot_for(kept, args.n_boot_confirmatory, args.boot_seed)
        summ = ef.summarise(effects[effects["arm"] != schema.CLEAN], boot)
        share = secondary.reconverged_share(effects, boot)  # BENCHMARK §7, beside reconverge_min
        summ = summ.merge(share, on=list(ef.BY_ARM), how="left")
        not_rc = summ["outcome"] != "reconverge_min"
        summ.loc[not_rc, ["reconverged_share", "reconverged_lo", "reconverged_hi"]] = np.nan
        results["effects_summary"] = summ
        results["hypotheses"], results["hypotheses_co_reported"] = hypotheses.run(effects, boot_c)
        results["c1_per_fault"], results["c1_summary"] = cross_layer.run(effects, boot_c)
        results.update(secondary.run(effects, kept, boot))
        if args.primary_input:  # S1 (BENCHMARK §8.1): tau-b against the primary ranking
            prim = schema.validate(schema.coerce(read_tables(args.primary_input)))
            if excluded:
                prim = prim[~prim["patient_id"].isin(excluded)]
            prim_kept, _ = ef.apply_invalid_rule(prim)
            tau = cross_layer.ranking_tau_b(effects, ef.paired_effects(prim_kept), boot_c)
            results["s1_tau_b"] = pd.DataFrame([tau])
        if args.additivity:
            records = read_additivity_records(Path(args.additivity), excluded)
            if any(not r.get("synthetic") and "faulted" in r for r in records):
                ok, why = protocol_tag_ok(args.protocol_tag)
                if not ok:
                    raise Refused(why)
            daily, events = secondary.read_additivity(records)
            results["additivity"] = secondary.additivity(daily, events, effects, boot)
    for name, df in results.items():
        df.to_csv(out / f"{name}.csv", index=False)
    paths = (
        [] if (args.no_figures or clean_only) else figures.all_figures(results, effects, fig_dir)
    )

    banner = "SYNTHETIC DATA: a test of the analysis code, not a result." if synthetic else ""
    md = ["# Paper 2 analysis summary", ""]
    if banner:
        md += [f"**{banner}**", ""]
    md += ["## Data", "", describe(table, counts), ""]
    md += [STANDARD_MIX_NOTE, ""]
    if excluded:
        md += [
            f"**S2: {len(excluded)} patients outside acceptance excluded before every estimate**"
            f" ({', '.join(excluded)}; {args.exclude_outside_acceptance}).",
            "",
        ]
    if clean_only:
        md += ["Clean arms only: descriptive tables (clean_levels.csv), no contrasts.", ""]
    else:
        h = results["hypotheses"]
        md += ["## Confirmatory family H1-H2 (oref0_tb; Holm, memo §9)", ""]
        md += [
            _md_table(
                h,
                ["id", "arm", "outcome", "sesoi", "theta", "ci_level", "ci_lo", "ci_hi"]
                + ["ci90_lo", "ci90_hi", "p_value", "rejected", "decision", "mc_flag"],
            ),
            "",
            "Co-reported outcomes (95% CI, uncorrected): hypotheses_co_reported.csv.",
            "",
            "## C1, the pre-specified primary comparison (estimation only)",
            "",
            _md_table(
                results["c1_summary"],
                ["analysis", "tau_b", "tau_b_lo", "tau_b_hi", "share_alarm", "share_insulin"]
                + ["share_diff", "share_diff_lo", "share_diff_hi", "discordant"],
            ),
            "",
            "Per-fault per-day harm: c1_per_fault.csv.",
            "",
            "## Secondary and exploratory (labelled, uncorrected, 95% CI)",
            "",
            "Gate optima (in-sample, descriptive; neither is a recommended setting):",
            "",
            _md_table(
                results["gate_optima"],
                ["harm_optimal", "alarm_optimal", "alarm_rule_feasible", "coincide_share"],
            ),
            "",
        ]
    if "additivity" in results:
        md += ["## D8 additivity (secondary, 95% CI)", ""]
        md += [
            _md_table(
                results["additivity"],
                ["outcome", "quantity", "theta", "ci_lo", "ci_hi", "n_patients", "additive"],
            ),
            "",
        ]
    md += ["## Files", ""]
    md += [f"- {name}.csv" for name in results]
    md += [f"- figures/{p.name}" for p in paths]
    manifest = {
        "commit": _commit(),
        "argv": sys.argv[1:] if argv is None else argv,
        "inputs": args.input,
        "data_source": sorted(sources),
        "n_boot": args.n_boot,
        "n_boot_confirmatory": args.n_boot_confirmatory,
        "boot_seed": args.boot_seed,
        "protocol_tag": args.protocol_tag,
        "excluded_outside_acceptance": excluded,
        "additivity": args.additivity,
        "seconds": round(time.perf_counter() - t_start, 1),
    }
    md += ["", "## Run", "", "```json", json.dumps(manifest, indent=1), "```", ""]
    (out / "summary.md").write_text("\n".join(md))
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1))
    print(f"analyse.py: {len(results)} tables and {len(paths)} figure files in {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
