"""Paired effects, the invalid rule, and patient-level and population estimates with CIs.

Estimand (memo §2): tau_a(p, r) = Y_a(p, r) - Y_clean(p, r) within one (controller, patient, seed,
slot); the patient effect is the mean over valid seeds and theta the mean over patients, equal
weight per patient. CIs (memo §10): percentile intervals from paper 1's patient bootstrap,
stratified by age group (`ghostlows.bench.stats.bootstrap_ci`), 10,000 replicates for H1, H2 and
C1 and 2,000 otherwise, with a fixed seed.

Joint resampling: `Bootstrap` records the patient indices that paper 1's `bootstrap_ci` draws (it
calls it with a statistic that stores each resample), so every quantity computed with one
`Bootstrap` uses the same index per replicate. Estimates for different arms, layers, faults and
controllers are therefore resampled jointly, as shares, rank correlations and contrasts need, and
percentile CIs are taken exactly as `bootstrap_ci` takes them.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from functools import lru_cache

import numpy as np
import pandas as pd
from ghostlows.bench.stats import bootstrap_ci

from .schema import CLEAN, COLUMNS, KEYS, OUTCOMES, UNIT_KEYS, ref_of

N_BOOT = 2000  # memo §10: secondary and exploratory estimates
N_BOOT_CONFIRMATORY = 10_000  # memo §10: H1, H2 and C1
BOOT_SEED = 0  # ASSUMPTION: memo §10 ("seed fixed"); paper 1's default seed
LEVEL = 0.95
MC_HALF = 0.5  # memo §10: flag when the mean within-patient SE exceeds half the between-patient SD
BY_ARM = ("controller", "slot", "arm")


# ---------------------------------------------------------------- the invalid rule (memo §10)


def apply_invalid_rule(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """D11 (BENCHMARK §2.6): an invalid arm drops every arm of its (patient, seed, slot), for all
    three controllers, so the controllers keep the same events.

    Returns (kept rows, counts per controller, slot and arm: rows, invalid arms (the triggering
    arms of that controller), dropped rows, and the dropped rows whose own controller had no
    invalid arm in that slot, i.e. dropped because another controller triggered)."""
    keys = ["patient_id", "seed", "slot"]
    bad = table.groupby(keys)["invalid"].transform("any").to_numpy(bool)
    own = table.groupby(list(UNIT_KEYS))["invalid"].transform("any").to_numpy(bool)
    counts = (
        table.assign(
            _invalid=table["invalid"].astype(int),
            _dropped=bad.astype(int),
            _other=(bad & ~own).astype(int),
        )
        .groupby(["controller", "slot", "arm"], sort=True)
        .agg(
            rows=("arm", "size"),
            invalid_arms=("_invalid", "sum"),
            dropped=("_dropped", "sum"),
            dropped_by_other_controller=("_other", "sum"),
        )
        .reset_index()
    )
    return table[~bad].reset_index(drop=True), counts


# ---------------------------------------------------------------- paired effects

CARRY = tuple(c for c, col in COLUMNS.items() if col.effect == "none" and not c.startswith("ref_"))


def paired_effects(table: pd.DataFrame) -> pd.DataFrame:
    """One row per arm row (the clean arm included, where every `diff` effect is 0): keys, design,
    fault parameters, onset state and flags, and the paired effect of every outcome under the
    outcome's own name (schema `Column.effect`)."""
    t = table.reset_index(drop=True)
    clean = t[t["arm"] == CLEAN].set_index(list(UNIT_KEYS))
    c = clean.reindex(pd.MultiIndex.from_frame(t[list(UNIT_KEYS)]))
    out = {k: t[k].to_numpy() for k in CARRY}
    for o in OUTCOMES:
        kind = COLUMNS[o].effect
        y = t[o].to_numpy(float)
        if kind == "diff":
            out[o] = y - c[o].to_numpy(float)
        elif kind == "ref":
            out[o] = y - t[ref_of(o)].to_numpy(float)
        else:  # paired
            out[o] = y
    return pd.DataFrame(out).astype({k: t[k].dtype for k in CARRY})


def arm_contrast(
    effects: pd.DataFrame, a: str, b: str, value_cols: Sequence[str], name: str | None = None
) -> pd.DataFrame:
    """Event-level tau_a - tau_b for arms a and b of the same (controller, patient, seed, slot):
    the arms share t0 and the snapshot, so the contrast is paired at the event level (memo §6)."""
    cols = list(value_cols)
    ea = effects[effects["arm"] == a].set_index(list(UNIT_KEYS))[cols]
    eb = effects[effects["arm"] == b].set_index(list(UNIT_KEYS))[cols]
    ea, eb = ea.align(eb, join="inner")
    out = (ea - eb).reset_index()
    return out.assign(contrast=name or f"{a} - {b}")


def controller_contrast(
    effects: pd.DataFrame, c1: str, c2: str, value_cols: Sequence[str]
) -> pd.DataFrame:
    """Event-level tau_c1 - tau_c2 for the same (patient, seed, slot, arm): two paired effects
    from different trunks (memo §10, controllers)."""
    keys = ["patient_id", "seed", "slot", "arm"]
    cols = list(value_cols)
    e1 = effects[effects["controller"] == c1].set_index(keys)[cols]
    e2 = effects[effects["controller"] == c2].set_index(keys)[cols]
    e1, e2 = e1.align(e2, join="inner")
    return (e1 - e2).reset_index().assign(contrast=f"{c1} - {c2}")


# ---------------------------------------------------------------- the patient bootstrap


@lru_cache(maxsize=8)
def _counts(patients: tuple[str, ...], n_boot: int, seed: int) -> np.ndarray:
    frame = pd.DataFrame({"patient_id": list(patients), "_row": np.arange(len(patients))})
    drawn: list[np.ndarray] = []

    def record(resample: pd.DataFrame) -> float:
        drawn.append(resample["_row"].to_numpy())
        return 0.0

    bootstrap_ci(frame, record, n_boot=n_boot, seed=seed)
    idx = np.vstack(drawn[1:])  # the first call is the estimate on the full sample
    counts = np.zeros((n_boot, len(patients)))
    np.add.at(counts, (np.repeat(np.arange(n_boot), idx.shape[1]), idx.ravel()), 1.0)
    counts.setflags(write=False)
    return counts


class Bootstrap:
    """Paper 1's stratified patient bootstrap as a (replicates x patients) count matrix.

    `counts[b, p]` is how often patient p is drawn in replicate b of `bootstrap_ci(n_boot, seed)`
    on the patients in sorted order. Means over resampled patients skip patients without a value
    (NaN), as a mean over the patients that have one."""

    def __init__(self, patients, n_boot: int = N_BOOT, seed: int = BOOT_SEED):
        self.patients = tuple(sorted(set(patients)))
        self.n_boot, self.seed = int(n_boot), int(seed)
        self.counts = _counts(self.patients, self.n_boot, self.seed)

    def align(self, per_patient: pd.DataFrame) -> np.ndarray:
        """(patients x quantities) values of a frame indexed by patient_id (missing: NaN)."""
        return per_patient.reindex(list(self.patients)).to_numpy(float)

    def means(self, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(estimate, replicates) of the mean over patients of each column of `values`."""
        v = np.asarray(values, float).reshape(len(self.patients), -1)
        ok = np.isfinite(v)
        filled = np.where(ok, v, 0.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            est = filled.sum(0) / ok.sum(0)
            reps = (self.counts @ filled) / (self.counts @ ok)
        return est, reps

    def sums(self, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(estimate, replicates) of sums over patients (missing values count 0)."""
        v = np.nan_to_num(np.asarray(values, float).reshape(len(self.patients), -1))
        return v.sum(0), self.counts @ v


def percentile_ci(draws: np.ndarray, level: float = LEVEL) -> tuple[np.ndarray, np.ndarray]:
    """Percentile interval per column, as paper 1's bootstrap_ci: non-finite draws are dropped
    and numpy's linear quantiles taken; NaN where no draw is finite."""
    d = np.asarray(draws, float)
    d = np.where(np.isfinite(d), d, np.nan)
    a = (1 - level) / 2
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        lo, hi = np.nanquantile(d, [a, 1 - a], axis=0)
    return lo, hi


def boot_for(table: pd.DataFrame, n_boot: int = N_BOOT, seed: int = BOOT_SEED) -> Bootstrap:
    """The bootstrap over every patient of the table (one index set for all analyses)."""
    return Bootstrap(table["patient_id"].astype(str).unique(), n_boot, seed)


# ---------------------------------------------------------------- estimates


def patient_stats(
    events: pd.DataFrame, value_cols: Sequence[str], by: Sequence[str] = BY_ARM
) -> pd.DataFrame:
    """Per (group, outcome, patient): mean over the patient's events (seeds), count and SD."""
    long = events.melt(
        id_vars=[*by, "patient_id"], value_vars=list(value_cols), var_name="outcome"
    ).dropna(subset=["value"])
    return (
        long.groupby([*by, "outcome", "patient_id"], sort=True, observed=True)["value"]
        .agg(["mean", "count", "std"])
        .reset_index()
    )


def estimate(
    events: pd.DataFrame,
    value_cols: Sequence[str],
    by: Sequence[str] = BY_ARM,
    boot: Bootstrap | None = None,
    level: float = LEVEL,
    keep_draws: bool = False,
):
    """theta (mean over patients of the patient means), its percentile CI and Monte Carlo error,
    per group of `by` and outcome. With `keep_draws`, also returns the (replicates x rows) draws
    in the row order of the table.

    Monte Carlo error (memo §10): `within_se` is the mean over patients of SD(seeds)/sqrt(n),
    `between_sd` the SD of the patient means, `mc_se` = sqrt(mean SE^2 / patients), the seed noise
    left in theta; `mc_flag` marks within_se > half of between_sd.
    """
    by = list(by)
    boot = boot or Bootstrap(events["patient_id"].unique())
    stats = patient_stats(events, value_cols, by)
    if stats.empty:
        cols = [*by, "outcome", "theta", "ci_lo", "ci_hi", "n_patients", "n_events"]
        empty = pd.DataFrame(columns=cols)
        return (empty, np.empty((boot.n_boot, 0))) if keep_draws else empty
    keys = [*by, "outcome"]
    wide = stats.pivot(index=keys, columns="patient_id", values="mean")
    wide = wide.reindex(columns=list(boot.patients))
    est, reps = boot.means(wide.to_numpy(float).T)
    lo, hi = percentile_ci(reps, level)
    se = stats.assign(se=stats["std"] / np.sqrt(stats["count"]))
    se.loc[stats["count"] < 2, "se"] = np.nan
    se["se2"] = se["se"] ** 2
    g = se.groupby(keys, sort=True, observed=True)
    mc = pd.DataFrame(
        {
            "n_patients": g["mean"].count(),
            "n_events": g["count"].sum(),
            "within_se": g["se"].mean(),
            "between_sd": g["mean"].std(ddof=1),
            "mc_se": np.sqrt(g["se2"].mean() / g["mean"].count()),
        }
    ).reindex(wide.index)
    out = wide.index.to_frame(index=False)
    out["theta"], out["ci_lo"], out["ci_hi"] = est, lo, hi
    out["level"] = level
    for k in mc:
        out[k] = mc[k].to_numpy()
    out["mc_flag"] = out["within_se"] > MC_HALF * out["between_sd"]
    out["n_events"] = out["n_events"].astype(int)
    out["n_patients"] = out["n_patients"].astype(int)
    return (out, reps) if keep_draws else out


def summarise(
    effects: pd.DataFrame, boot: Bootstrap, outcomes: Sequence[str] = OUTCOMES
) -> pd.DataFrame:
    """theta with CI and Monte Carlo error for every controller, slot, arm and outcome."""
    return estimate(effects, list(outcomes), BY_ARM, boot)


def lookup(summary: pd.DataFrame, **sel) -> pd.Series:
    """The one row of a summary that matches `sel` (raises if there is not exactly one)."""
    m = np.ones(len(summary), bool)
    for k, v in sel.items():
        m &= (summary[k] == v).to_numpy()
    if m.sum() != 1:
        raise KeyError(f"{int(m.sum())} rows match {sel}")
    return summary[m].iloc[0]


__all__ = [
    "BOOT_SEED",
    "KEYS",
    "N_BOOT",
    "N_BOOT_CONFIRMATORY",
    "Bootstrap",
    "apply_invalid_rule",
    "arm_contrast",
    "boot_for",
    "controller_contrast",
    "estimate",
    "paired_effects",
    "patient_stats",
    "percentile_ci",
    "summarise",
]
