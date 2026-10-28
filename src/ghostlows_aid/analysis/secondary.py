"""Secondary and exploratory analyses of memo §10, each labelled and uncorrected, with 95% CIs
from the secondary patient bootstrap (n_boot = 2,000; one resampling index shared by every
analysis in this module, so that the gate's two optima can be compared replicate by replicate).

Every function takes the paired effects (`effects.paired_effects` after the invalid rule) and the
bootstrap; `pairing_efficiency` and the alarm side of `gate` read the outcome table itself (levels,
not differences). Contrasts between arms of one slot are formed event by event (`combine`): the
arms share t0 and the snapshot (memo §6), so they are paired at the event level.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from ghostlows.bench.stats import age_group, choose_threshold, patient_folds

from . import effects as ef
from .cross_layer import RATES
from .schema import (
    ACTIONS,
    ARMS,
    CLEAN,
    GATE_SETTINGS,
    GATES,
    PRIMARY_CONTROLLER,
    SECOND_DETECTOR,
    UNIT_KEYS,
)

HARM = "risk_weighted_min"  # memo §8: the one harm scale
# ASSUMPTION: the "primary outcomes" of N1 are those of H1, H2 and C1.
PRIMARY_OUTCOMES = ("minutes_above_180", "minutes_below_70", HARM, "trend_false_alarms")
DOSE_OUTCOMES = ("minutes_above_180", HARM)
LADDER_OUTCOMES = (HARM, "minutes_below_70", "minutes_above_180")
STATE_OUTCOMES = (HARM,)
PATHWAY_OUTCOMES = (HARM, "minutes_above_180", "minutes_below_70", "minutes_below_54")
GUARD_OUTCOMES = (HARM, "minutes_below_70", "minutes_above_180", "insulin_total_u")
CONTROLLER_OUTCOMES = (HARM, "minutes_below_54", "minutes_below_70", "minutes_above_180")
INSULIN_PARTS = (
    "insulin_total_u",
    "insulin_smb_u",
    "insulin_tb_dev_u",
    "insulin_bolus_u",
    "insulin_correction_u",
)
MECHANISM_OUTCOMES = (
    "insulin_during_u",
    "insulin_after_u",
    "minutes_above_180_during",
    "minutes_above_180_after",
    "rescue_carbs_g",
    "rescue_count",
    *(f"actions_{a}" for a in ACTIONS),
)
CROSSINGS = ("ctrl_below_70_during", "ctrl_above_180_during")
PY_OUTCOMES = (HARM, "minutes_below_54", "minutes_below_70", "minutes_above_180")
FAULT_ARMS = tuple(a for a, x in ARMS.items() if a != CLEAN and not x.variant)
CLOSED_LOOP = ("oref0_tb", "oref1_smb")

# memo §10 N8: typical against severe for every fault; drift typical, +-20% and +-40% gain;
# interference small against large. Each rung is contrasted with the one before, event by event.
LADDERS = {
    "compression": ("pisa_typical", "pisa_severe"),
    "dropout": ("dropout_short", "dropout_long", "dropout_severe"),
    "frozen": ("frozen_typical", "frozen_severe"),
    "spikes": ("spike_typical", "spike_severe"),
    "drift_low": ("drift_low", "drift_low_severe", "gain_low"),
    "drift_high": ("drift_high", "drift_high_severe", "gain_high"),
    "interference": ("interference_small", "interference_large"),
}

# memo §10 Q4
MIN_CELL_PATIENTS = 5
GLUCOSE_CUTS = (100.0, 140.0)  # < 100, 100-140, > 140 mg/dL (true glucose at onset)
TREND_CUT = 1.0  # CGM trend falling < -1, rising > +1 mg/dL/min, else flat
NIGHT_HOURS = (23.0, 7.0)  # ASSUMPTION: night = paper 1's PISA window, 23:00 to 07:00
STATE_COVARIATES = (
    "onset_truth_mg_dl",
    "onset_cgm_trend",
    "onset_iob_u",
    "onset_cob_g",
    "onset_hours_since_bolus",
)
N_FOLDS, FOLD_SEED = 5, 0  # paper 1's patient_folds, seed 0

PATHWAY_BASES = ("pisa_typical", "interference_large")
WAKE_Q = (0.0, 0.25, 0.5, 0.75, 1.0)  # memo §10 Q3

NO_GATE = "no_gate"
GATE_OPTIONS = (NO_GATE, *GATE_SETTINGS)  # memo §10 N10: the five candidate settings
GATE_DELTA = 0.0  # paper 1's pre-registered D5 margin (ghostlows.bench.pipeline, delta = 0)

FALSE_HIGH_ARMS = (
    "drift_high",
    "drift_high_severe",
    "gain_high",
    "interference_small",
    "interference_large",
)
SPIKE_ARMS = ("spike_typical", "spike_severe")

_PAPER1 = {  # per-day rates of the faults without a published range (paper 1, as in C1)
    "frozen_typical": RATES["frozen_typical"],
    "spike_typical": RATES["spike_typical"],
    "drift_low": RATES["drift_low"],
    "drift_high": RATES["drift_high"],
}
# memo §10 Q7: paper 1's rates (central) and the lowest and highest published ones. Interference
# and gain errors are per exposure (no population rate) and are left out.
# ASSUMPTION: the published gap rates are the G6 trial's long gaps (arm 4), which replace paper
# 1's short gaps in the low and high scenarios (as in C1's sensitivity analysis).
SCENARIOS = {
    "central": {
        "pisa_typical": RATES["pisa_nights_central"],
        "dropout_short": RATES["dropout_short"],
        **_PAPER1,
    },
    "low": {
        "pisa_typical": RATES["pisa_nights_low"],
        "dropout_long": RATES["dropout_long_low"],
        **_PAPER1,
    },
    "high": {
        "pisa_typical": RATES["pisa_nights_high"],
        "dropout_long": RATES["dropout_long_high"],
        **_PAPER1,
    },
}
DAYS_PER_YEAR = 365


# ---------------------------------------------------------------- helpers


def combine(
    effects: pd.DataFrame,
    weights: Mapping[str, float],
    cols: Sequence[str],
    name: str,
    controller: str = PRIMARY_CONTROLLER,
) -> pd.DataFrame:
    """Event-level sum of weight x paired effect over arms of one slot, on the units that have
    every arm (a zero weight only restricts to common units). Columns: unit keys, `arm` = name,
    and `cols`."""
    cols = list(cols)
    e = effects[effects["controller"] == controller]
    parts = [e[e["arm"] == a].set_index(list(UNIT_KEYS))[cols] * w for a, w in weights.items()]
    if not parts or any(p.empty for p in parts):
        return pd.DataFrame(columns=[*UNIT_KEYS, "arm", *cols])
    idx = parts[0].index
    for p in parts[1:]:
        idx = idx.intersection(p.index)
    out = parts[0].loc[idx]
    for p in parts[1:]:
        out = out + p.loc[idx]
    return out.reset_index().assign(arm=name)


def _estimate(frames, cols, boot, by=ef.BY_ARM, keep_draws=False):
    frames = [f for f in frames if len(f)]
    if not frames:
        empty = pd.DataFrame(columns=[*by, "outcome", "theta", "ci_lo", "ci_hi", "n_patients"])
        return (empty, np.empty((boot.n_boot, 0))) if keep_draws else empty
    events = pd.concat(frames, ignore_index=True)
    return ef.estimate(events, list(cols), list(by), boot, keep_draws=keep_draws)


def _ci(draws) -> tuple[float, float]:
    lo, hi = ef.percentile_ci(np.asarray(draws, float).reshape(len(draws), -1))
    return float(lo[0]), float(hi[0])


def _present(effects: pd.DataFrame, cols: Sequence[str]) -> list[str]:
    return [c for c in cols if c in effects]


# ---------------------------------------------------------------- N1 pairing efficiency


def pairing_efficiency(table: pd.DataFrame, outcomes=PRIMARY_OUTCOMES) -> pd.DataFrame:
    """memo §10 N1: per patient, Var_r[Y_a(r) - Y_clean(r)] against Var_r[Y_a(r) - Y_clean(r + 1
    mod N)] over the patient's valid seeds in order; the median ratio (paired / cross-seed) over
    patients and the implied number of independent runs per paired run (1 / median ratio)."""
    cols = _present(table, outcomes)
    keys = list(UNIT_KEYS)
    clean = table[table["arm"] == CLEAN][[*keys, *cols]]
    arms = table[table["arm"] != CLEAN][[*keys, "arm", *cols]]
    m = arms.merge(clean, on=keys, suffixes=("", "_clean"))
    m = m.sort_values(["controller", "arm", "patient_id", "seed"]).reset_index(drop=True)
    g = m.groupby(["controller", "arm", "patient_id"], sort=False)
    pos, size = g.cumcount().to_numpy(), g["seed"].transform("size").to_numpy()
    row = np.arange(len(m))
    nxt = np.where(pos == size - 1, row - pos, row + 1)  # seed r + 1 mod N within the patient
    parts = {}
    for o in cols:
        y, yc = m[o].to_numpy(float), m[f"{o}_clean"].to_numpy(float)
        parts[f"{o}|paired"] = y - yc
        parts[f"{o}|cross"] = y - yc[nxt]
    d = pd.concat([m[["controller", "slot", "arm", "patient_id"]], pd.DataFrame(parts)], axis=1)
    var = d.groupby(["controller", "slot", "arm", "patient_id"], sort=True).var(ddof=1)
    rows = []
    for o in cols:
        vp, vc = var[f"{o}|paired"], var[f"{o}|cross"]
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = (vp / vc).where(~((vp == 0) & (vc == 0)))
        r = pd.DataFrame({"ratio": ratio, "vp": vp, "vc": vc}).dropna(subset=["ratio"])
        for (c, s, a), x in r.groupby(level=[0, 1, 2], sort=True):
            med = float(np.median(x["ratio"]))
            rows.append(
                {
                    "controller": c,
                    "slot": s,
                    "arm": a,
                    "outcome": o,
                    "median_ratio": med,
                    "implied_runs": (1 / med) if med > 0 else math.inf,
                    "median_var_paired": float(np.median(x["vp"])),
                    "median_var_cross": float(np.median(x["vc"])),
                    "n_patients": len(x),
                }
            )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- N8 dose-response and ladders


def dose_response(
    effects: pd.DataFrame, boot: ef.Bootstrap, outcomes=DOSE_OUTCOMES, arm: str = "pisa_typical"
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """memo §10 N8: (slopes, effects by threshold crossing). The slope is the per-patient
    least-squares slope of the event effect on the realised artifact area (depth x duration,
    mg/dL h), averaged over patients (patients with fewer than two distinct areas left out);
    the split is by whether the controller's input fell below 70 mg/dL during the fault."""
    cols = _present(effects, outcomes)
    e = effects[effects["arm"] == arm]
    rows = []
    for c, sub in e.groupby("controller", sort=True):
        x = sub["fault_area_mg_dl_h"].to_numpy(float)
        slopes = {}
        for o in cols:
            y = sub[o].to_numpy(float)
            ok = np.isfinite(x) & np.isfinite(y)
            f = pd.DataFrame({"p": sub["patient_id"].to_numpy()[ok], "x": x[ok], "y": y[ok]})
            gm = f.groupby("p")[["x", "y"]].transform("mean")
            f["sxy"] = (f["x"] - gm["x"]) * (f["y"] - gm["y"])
            f["sxx"] = (f["x"] - gm["x"]) ** 2
            s = f.groupby("p")[["sxy", "sxx"]].sum()
            slopes[o] = (s["sxy"] / s["sxx"]).where(s["sxx"] > 1e-12)
        wide = pd.DataFrame(slopes)
        if wide.empty:
            continue
        est, reps = boot.means(boot.align(wide))
        lo, hi = ef.percentile_ci(reps)
        for k, o in enumerate(cols):
            rows.append(
                {
                    "controller": c,
                    "arm": arm,
                    "outcome": o,
                    "slope_per_mg_dl_h": est[k],
                    "ci_lo": lo[k],
                    "ci_hi": hi[k],
                    "n_patients": int(wide[o].notna().sum()),
                }
            )
    crossed = e[e["controller"].isin(CLOSED_LOOP)].dropna(subset=["ctrl_below_70_during"])
    crossed = crossed.assign(crossed_70=crossed["ctrl_below_70_during"].astype(bool))
    split = _estimate([crossed], cols, boot, by=(*ef.BY_ARM, "crossed_70"))
    return pd.DataFrame(rows), split


def ladders(effects: pd.DataFrame, boot: ef.Bootstrap, outcomes=LADDER_OUTCOMES) -> pd.DataFrame:
    """memo §10 N8: each rung's effect (against clean) and each rung against the rung before,
    on shared events."""
    cols = _present(effects, outcomes)
    frames = []
    for c in sorted(effects["controller"].unique()):
        for name, rungs in LADDERS.items():
            for k, arm in enumerate(rungs):
                one = combine(effects, {arm: 1.0}, cols, arm, c)
                frames.append(one.assign(ladder=name, rung=k, versus=CLEAN))
                if k:
                    prev = rungs[k - 1]
                    step = combine(effects, {arm: 1.0, prev: -1.0}, cols, arm, c)
                    frames.append(step.assign(ladder=name, rung=k, versus=prev))
    return _estimate(
        frames, cols, boot, by=("controller", "slot", "ladder", "rung", "arm", "versus")
    )


# ---------------------------------------------------------------- Q4 state dependence


def _cells(e: pd.DataFrame) -> pd.DataFrame:
    """The event's cell in each pre-specified grouping of memo §10 Q4 (long format)."""
    g = e["onset_truth_mg_dl"].to_numpy(float)
    t = e["onset_cgm_trend"].to_numpy(float)
    lo, hi = GLUCOSE_CUTS
    band = np.select([g < lo, g <= hi, g > hi], ["<100", "100-140", ">140"], "")
    trend = np.select(
        [t < -TREND_CUT, t > TREND_CUT, np.isfinite(t)], ["falling", "rising", "flat"], ""
    )
    both = np.where((band != "") & (trend != ""), np.char.add(np.char.add(band, " & "), trend), "")
    h = e["onset_hour"].to_numpy(float)
    night = np.where(
        np.isfinite(h), np.where((h >= NIGHT_HOURS[0]) | (h < NIGHT_HOURS[1]), "night", "day"), ""
    )
    iob = e["onset_iob_u"].astype(float)
    groups = e.groupby(["controller", "slot"], sort=False)["onset_iob_u"]
    q1 = groups.transform(lambda s: s.astype(float).quantile(1 / 3)).to_numpy(float)
    q2 = groups.transform(lambda s: s.astype(float).quantile(2 / 3)).to_numpy(float)
    i = iob.to_numpy(float)
    tert = np.select([~np.isfinite(i), i <= q1, i <= q2], ["", "T1", "T2"], "T3")
    out = []
    for kind, cell in (("glucose_x_trend", both), ("iob_tertile", tert), ("night_day", night)):
        out.append(e.assign(grouping=kind, cell=cell)[cell != ""])
    return pd.concat(out, ignore_index=True)


def state_cells(
    effects: pd.DataFrame,
    boot: ef.Bootstrap,
    outcomes=STATE_OUTCOMES,
    controllers=(PRIMARY_CONTROLLER,),
    arms=FAULT_ARMS,
) -> pd.DataFrame:
    """memo §10 Q4 (exploratory): effects in the onset-state cells (true glucose x CGM trend; IOB
    tertile within controller and slot; night or day), patient-bootstrap CIs; a cell with fewer
    than 5 patients is not estimable (estimate and CI left empty)."""
    cols = _present(effects, outcomes)
    e = effects[effects["controller"].isin(controllers) & effects["arm"].isin(arms)]
    if e.empty:
        return _estimate([], cols, boot, by=(*ef.BY_ARM, "grouping", "cell"))
    est = _estimate([_cells(e)], cols, boot, by=(*ef.BY_ARM, "grouping", "cell"))
    est["estimable"] = est["n_patients"] >= MIN_CELL_PATIENTS
    est.loc[~est["estimable"], ["theta", "ci_lo", "ci_hi"]] = np.nan
    return est


def _fe_parts(x: np.ndarray, y: np.ndarray, p: np.ndarray, patients: Sequence[str]):
    """Per patient X~'X~ and X~'y~ after demeaning within patient (the fixed-effects model)."""
    k = x.shape[1]
    S, s = np.zeros((len(patients), k, k)), np.zeros((len(patients), k))
    where = {q: i for i, q in enumerate(patients)}
    for q in np.unique(p):
        m = p == q
        xd, yd = x[m] - x[m].mean(0), y[m] - y[m].mean()
        S[where[q]], s[where[q]] = xd.T @ xd, xd.T @ yd
    return S, s


def state_model(
    effects: pd.DataFrame,
    boot: ef.Bootstrap,
    outcomes=STATE_OUTCOMES,
    controllers=(PRIMARY_CONTROLLER,),
    arms=FAULT_ARMS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """memo §10 Q4: (coefficients, held-out check). A linear model of the event effect on the
    onset state with patient fixed effects, coefficients with patient-bootstrap CIs (a resampled
    patient's events enter as often as the patient is drawn); covariates without within-patient
    variation are left out. Held-out: the model without patient effects, plus age group, fitted
    in paper 1's five leave-patients-out folds; R^2 on held-out events against the training mean."""
    cols = _present(effects, outcomes)
    coefs, held = [], []
    e = effects[effects["controller"].isin(controllers) & effects["arm"].isin(arms)]
    for (c, arm), sub in e.groupby(["controller", "arm"], sort=True):
        for o in cols:
            f = sub[["patient_id", *STATE_COVARIATES, o]].astype({o: float})
            f = f.dropna()
            if f.empty:
                continue
            x = f[list(STATE_COVARIATES)].to_numpy(float)
            y, p = f[o].to_numpy(float), f["patient_id"].to_numpy(str)
            S, s = _fe_parts(x, y, p, boot.patients)
            keep = np.diagonal(S.sum(0)) > 1e-9 * max(1.0, float(np.abs(S).max()))
            idx = np.flatnonzero(keep)
            beta = np.full(len(STATE_COVARIATES), np.nan)
            reps = np.full((boot.n_boot, len(STATE_COVARIATES)), np.nan)
            if idx.size:
                Sk, sk = S[:, idx][:, :, idx], s[:, idx]
                beta[idx] = np.linalg.lstsq(Sk.sum(0), sk.sum(0), rcond=None)[0]
                A = np.einsum("bp,pij->bij", boot.counts, Sk)
                v = boot.counts @ sk
                reps[:, idx] = np.einsum("bij,bj->bi", np.linalg.pinv(A), v)
            lo, hi = ef.percentile_ci(reps)
            for k, name in enumerate(STATE_COVARIATES):
                coefs.append(
                    {
                        "controller": c,
                        "arm": arm,
                        "outcome": o,
                        "term": name,
                        "coef": beta[k],
                        "ci_lo": lo[k],
                        "ci_hi": hi[k],
                        "n_events": len(f),
                        "n_patients": len(np.unique(p)),
                    }
                )
            held.append({"controller": c, "arm": arm, "outcome": o, **_held_out(x, y, p)})
    return pd.DataFrame(coefs), pd.DataFrame(held)


def _held_out(x: np.ndarray, y: np.ndarray, p: np.ndarray) -> dict:
    groups = sorted({age_group(q) for q in p})
    dummies = [np.array([age_group(q) == g for q in p], float) for g in groups[1:]]
    design = np.column_stack([np.ones(len(y)), x, *dummies])
    folds = patient_folds(np.unique(p), N_FOLDS, FOLD_SEED)
    fold = np.array([folds[q] for q in p])
    sse = sst = 0.0
    per_fold = []
    for k in range(N_FOLDS):
        test, train = fold == k, fold != k
        if not test.any() or train.sum() <= design.shape[1]:
            continue
        b = np.linalg.lstsq(design[train], y[train], rcond=None)[0]
        e1 = float(((y[test] - design[test] @ b) ** 2).sum())
        e0 = float(((y[test] - y[train].mean()) ** 2).sum())
        sse, sst = sse + e1, sst + e0
        per_fold.append(1 - e1 / e0 if e0 > 0 else np.nan)
    r2 = 1 - sse / sst if sst > 0 else np.nan
    return {
        "r2_out_of_fold": r2,
        "r2_fold_min": np.nanmin(per_fold) if per_fold else np.nan,
        "r2_fold_max": np.nanmax(per_fold) if per_fold else np.nan,
        "n_folds": len(per_fold),
        "n_events": len(y),
    }


# ---------------------------------------------------------------- Q3 pathway


def pathway(
    effects: pd.DataFrame,
    boot: ef.Bootstrap,
    outcomes=PATHWAY_OUTCOMES,
    controller: str = PRIMARY_CONTROLLER,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """memo §10 Q3: (2x2 decomposition, wake curve). Total = both, controller-only, person-only
    and interaction = total - controller - person, on the events that have all three arms; the
    wake curve q x tau_both + (1 - q) x tau_controller-only for PISA typical."""
    cols = _present(effects, outcomes)
    frames, wake = [], []
    for base in PATHWAY_BASES:
        ctrl, person = f"{base}+controller_only", f"{base}+person_only"
        parts = {
            "total": {base: 1.0, ctrl: 0.0, person: 0.0},
            "controller_only": {base: 0.0, ctrl: 1.0, person: 0.0},
            "person_only": {base: 0.0, ctrl: 0.0, person: 1.0},
            "interaction": {base: 1.0, ctrl: -1.0, person: -1.0},
        }
        for comp, w in parts.items():
            frames.append(combine(effects, w, cols, base, controller).assign(component=comp))
        if base == "pisa_typical":
            for q in WAKE_Q:
                f = combine(effects, {base: q, ctrl: 1 - q}, cols, base, controller)
                wake.append(f.assign(q=q))
    by = (*ef.BY_ARM, "component")
    return _estimate(frames, cols, boot, by=by), _estimate(wake, cols, boot, by=(*ef.BY_ARM, "q"))


# ---------------------------------------------------------------- N8 mechanism


def mechanism(
    effects: pd.DataFrame, boot: ef.Bootstrap, outcomes=MECHANISM_OUTCOMES
) -> pd.DataFrame:
    """memo §10 N8: insulin and minutes above 180 during against after the fault, rescue
    carbohydrate and controller-action counts (paired effects), with the share of events whose
    controller input crossed 70 or 180 mg/dL during the fault (descriptive, memo §8)."""
    cols = _present(effects, outcomes)
    e = effects[effects["arm"].isin(FAULT_ARMS)]
    shares = e[e["controller"].isin(CLOSED_LOOP)]
    share_cols = {f"share_{c}": shares[c].astype("Float64").astype(float) for c in CROSSINGS}
    shares = shares.assign(**share_cols)
    a = _estimate([e], cols, boot)
    b = _estimate([shares], list(share_cols), boot)
    return pd.concat([a, b], ignore_index=True)


# ---------------------------------------------------------------- Q5 gate, N10 optima


def _alarm_choice(fa, sens2, sens1, n2) -> tuple[int, bool, int]:
    """Paper 1's operating-point rule over the five settings (index order = GATE_OPTIONS; ties go
    to the earlier one): fewest false alarms per night with level 2 sensitivity at least the
    ungated alarm's minus delta; level 1 when the nights hold no level 2 low (paper 1, A7)."""
    if not (np.all(np.isfinite(fa))):
        return -1, False, 0
    level, sens = (2, sens2) if n2 > 0 else (1, sens1)
    if not np.all(np.isfinite(sens)):
        return -1, False, level
    best, info = choose_threshold(
        np.arange(len(fa), dtype=float),
        lambda t: float(sens[int(t)]),
        lambda t: float(fa[int(t)]),
        float(sens[0]),
        GATE_DELTA,
    )
    return int(best), bool(info["feasible"]), level


def gate(
    effects: pd.DataFrame,
    table: pd.DataFrame,
    boot: ef.Bootstrap,
    controller: str = PRIMARY_CONTROLLER,
    alarm_controller: str = "open_loop",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """memo §10 Q5 and N10: (curves, optima).

    Harm (risk-weighted minutes, `controller`): with-gate minus without-gate on a compression
    night, (PISA + gate) - PISA, and on a night without, (clean + gate) - clean; net per night =
    w x the first + (1 - w) x the second, w = 0.86 (RATES); negative is a benefit; no gate is 0.
    Alarms (`alarm_controller` night-slot arms, levels): false alarms per night and level 2
    sensitivity of paper 1's ungated trend alarm (no gate) and its gated versions, the nights
    weighted w and 1 - w. The harm-optimal setting has the lowest net effect; the alarm-optimal
    one is paper 1's rule over the same five. Both are found in every bootstrap replicate; the
    second detector (matched filter) is reported alongside, outside the choice.
    """
    w = RATES["pisa_typical"]
    frames = []
    for g in (*GATES, SECOND_DETECTOR):
        dp = combine(
            effects, {f"pisa_typical+{g}": 1.0, "pisa_typical": -1.0}, [HARM], g, controller
        )
        dc = combine(effects, {f"clean+{g}": 1.0}, [HARM], g, controller)
        frames += [dp.assign(quantity="delta_pisa"), dc.assign(quantity="delta_clean")]
    harm, hd = _estimate(
        frames, [HARM], boot, by=("controller", "arm", "quantity"), keep_draws=True
    )
    look = {(r.arm, r.quantity): k for k, r in enumerate(harm.itertuples())}

    def h(g, q):
        k = look.get((g, q))
        if k is None:
            return np.nan, np.full(boot.n_boot, np.nan), 0
        return float(harm["theta"].iloc[k]), hd[:, k], int(harm["n_patients"].iloc[k])

    prefixes = {NO_GATE: "trend", **dict(zip(GATE_SETTINGS, GATES, strict=True))}
    prefixes[SECOND_DETECTOR] = SECOND_DETECTOR
    fields = ("false_alarms", "l2_caught", "l2_missed", "l1_caught", "l1_missed")
    cols = [f"{p}_{f}" for p in prefixes.values() for f in fields]
    raw = table[
        (table["controller"] == alarm_controller)
        & (table["slot"] == "night")
        & (table["arm"].isin(["pisa_typical", CLEAN]))
    ]
    alarm, ad = _estimate([raw], _present(raw, cols), boot, by=("arm",), keep_draws=True)
    alook = {(r.arm, r.outcome): k for k, r in enumerate(alarm.itertuples())}

    def night(col):  # (estimate, replicates) of w x PISA night + (1 - w) x clean night
        kp, kc = alook.get(("pisa_typical", col)), alook.get((CLEAN, col))
        if kp is None or kc is None:
            return np.nan, np.full(boot.n_boot, np.nan)
        est = w * alarm["theta"].iloc[kp] + (1 - w) * alarm["theta"].iloc[kc]
        return float(est), w * ad[:, kp] + (1 - w) * ad[:, kc]

    def alarms(option):
        p = prefixes[option]
        fa = night(f"{p}_false_alarms")
        out = {"fa": fa}
        for lv in (1, 2):
            c, m = night(f"{p}_l{lv}_caught"), night(f"{p}_l{lv}_missed")
            with np.errstate(invalid="ignore", divide="ignore"):
                out[f"s{lv}"] = (c[0] / (c[0] + m[0]), c[1] / (c[1] + m[1]))
            out[f"n{lv}"] = (c[0] + m[0], c[1] + m[1])
        return out

    rows, net_est, net_reps = [], [], []
    al = {o: alarms(o) for o in (*GATE_OPTIONS, SECOND_DETECTOR)}
    for option in (*GATE_OPTIONS, SECOND_DETECTOR):
        rec = {"option": option, "prefix": prefixes[option]}
        if option == NO_GATE:
            zero = np.zeros(boot.n_boot)
            dp = dc = (0.0, zero, 0)
        else:
            g = SECOND_DETECTOR if option == SECOND_DETECTOR else f"baysal_{option}"
            dp, dc = h(g, "delta_pisa"), h(g, "delta_clean")
        net = (w * dp[0] + (1 - w) * dc[0], w * dp[1] + (1 - w) * dc[1])
        if option != SECOND_DETECTOR:
            net_est.append(net[0])
            net_reps.append(net[1])
        for name, (est, reps) in (
            ("delta_pisa", dp[:2]),
            ("delta_clean", dc[:2]),
            ("net_per_night", net),
            ("false_alarms_per_night", al[option]["fa"]),
            ("l2_sensitivity", al[option]["s2"]),
            ("l1_sensitivity", al[option]["s1"]),
        ):
            rec[name] = est
            rec[f"{name}_lo"], rec[f"{name}_hi"] = _ci(reps)
        rec["n_patients_harm"] = min(dp[2], dc[2]) if option != NO_GATE else np.nan
        rows.append(rec)
    curves = pd.DataFrame(rows)

    net_est, net_reps = np.array(net_est), np.column_stack(net_reps)
    harm_opt = int(np.argmin(net_est)) if np.all(np.isfinite(net_est)) else -1
    ok = np.all(np.isfinite(net_reps), axis=1)
    harm_b = np.where(ok, np.argmin(np.where(np.isfinite(net_reps), net_reps, np.inf), axis=1), -1)

    def choice(b=None):
        pick = (lambda v: v[0]) if b is None else (lambda v: v[1][b])
        fa = np.array([pick(al[o]["fa"]) for o in GATE_OPTIONS])
        s2 = np.array([pick(al[o]["s2"]) for o in GATE_OPTIONS])
        s1 = np.array([pick(al[o]["s1"]) for o in GATE_OPTIONS])
        n2 = pick(al[NO_GATE]["n2"])  # the nights' level 2 lows (the same for every alarm)
        return _alarm_choice(fa, s2, s1, n2)

    alarm_opt, feasible, level = choice()
    alarm_b = np.array([choice(b)[0] for b in range(boot.n_boot)])
    both = (harm_b >= 0) & (alarm_b >= 0)
    optima = pd.DataFrame(
        [
            {
                "harm_optimal": GATE_OPTIONS[harm_opt] if harm_opt >= 0 else "",
                "alarm_optimal": GATE_OPTIONS[alarm_opt] if alarm_opt >= 0 else "",
                "alarm_rule_feasible": feasible,
                "alarm_rule_level": level,
                "coincide_share": float(np.mean(harm_b[both] == alarm_b[both]))
                if both.any()
                else np.nan,
                "n_replicates": int(both.sum()),
                "w_compression_night": w,
                "delta": GATE_DELTA,
                "harm_controller": controller,
                "alarm_controller": alarm_controller,
            }
        ]
    )
    return curves, optima


# ---------------------------------------------------------------- Q5 guards, Q6 controllers


def guards(effects: pd.DataFrame, boot: ef.Bootstrap, outcomes=GUARD_OUTCOMES) -> pd.DataFrame:
    """memo §10 Q5: flat-CGM check = (frozen severe, fakecgm) - (frozen severe) on oref0_tb;
    stale-data stop = (frozen severe) - (dropout severe) on the same t0, for both closed loops."""
    cols = _present(effects, outcomes)
    frames = [
        combine(
            effects,
            {"frozen_severe+fakecgm": 1.0, "frozen_severe": -1.0},
            cols,
            "flat_cgm_check",
            PRIMARY_CONTROLLER,
        )
    ]
    for c in CLOSED_LOOP:
        w = {"frozen_severe": 1.0, "dropout_severe": -1.0}
        frames.append(combine(effects, w, cols, "stale_data_stop", c))
    return _estimate(frames, cols, boot)


def controllers(
    effects: pd.DataFrame, boot: ef.Bootstrap, outcomes=CONTROLLER_OUTCOMES
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """memo §10 Q6: (oref1_smb - oref0_tb for every fault, paired by patient, seed and slot;
    the SMB and temp-basal parts of oref1_smb's excess insulin for false highs and spikes, with
    their contrast against oref0_tb)."""
    cols = _present(effects, [*outcomes, *INSULIN_PARTS])
    e = effects[effects["arm"].isin(FAULT_ARMS)]
    cc = ef.controller_contrast(e, "oref1_smb", PRIMARY_CONTROLLER, cols)
    contrast = _estimate([cc], cols, boot, by=("contrast", "slot", "arm"))
    parts = _present(effects, INSULIN_PARTS)
    sel = e["arm"].isin((*FALSE_HIGH_ARMS, *SPIKE_ARMS)) & (e["controller"] == "oref1_smb")
    split = _estimate([e[sel]], parts, boot)
    return contrast, split


# ---------------------------------------------------------------- Q7 person-year


def person_year(effects: pd.DataFrame, boot: ef.Bootstrap, outcomes=PY_OUTCOMES) -> pd.DataFrame:
    """memo §10 Q7: 365 x sum over faults of rate x theta per scenario (SCENARIOS), with the
    patient-bootstrap uncertainty of theta; NaN when an arm is missing. Assumes additivity."""
    cols = _present(effects, outcomes)
    arms = sorted({a for s in SCENARIOS.values() for a in s})
    rows = []
    for c in sorted(effects["controller"].unique()):
        e = effects[(effects["controller"] == c) & effects["arm"].isin(arms)]
        est, draws = _estimate([e], cols, boot, keep_draws=True)
        look = {(r.arm, r.outcome): k for k, r in enumerate(est.itertuples())}
        for name, rates in SCENARIOS.items():
            for o in cols:
                total, reps = 0.0, np.zeros(boot.n_boot)
                for a, r in rates.items():
                    k = look.get((a, o))
                    t, d = (
                        (est["theta"].iloc[k], draws[:, k]) if k is not None else (np.nan, np.nan)
                    )
                    total, reps = total + r * t, reps + r * d
                lo, hi = _ci(DAYS_PER_YEAR * np.broadcast_to(reps, (boot.n_boot,)))
                rows.append(
                    {
                        "controller": c,
                        "scenario": name,
                        "outcome": o,
                        "per_person_year": DAYS_PER_YEAR * total,
                        "ci_lo": lo,
                        "ci_hi": hi,
                        "rates_per_day": ", ".join(f"{a} {r:g}" for a, r in rates.items()),
                    }
                )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- all


def run(effects: pd.DataFrame, table: pd.DataFrame, boot: ef.Bootstrap) -> dict[str, pd.DataFrame]:
    """Every secondary analysis, by output name (each table labelled secondary or exploratory
    by the caller)."""
    out = {"pairing_efficiency": pairing_efficiency(table)}
    out["dose_response"], out["dose_by_crossing"] = dose_response(effects, boot)
    out["ladders"] = ladders(effects, boot)
    out["state_cells"] = state_cells(effects, boot)
    out["state_model"], out["state_held_out"] = state_model(effects, boot)
    out["pathway"], out["wake_curve"] = pathway(effects, boot)
    out["mechanism"] = mechanism(effects, boot)
    out["gate_curves"], out["gate_optima"] = gate(effects, table, boot)
    out["guards"] = guards(effects, boot)
    out["controllers"], out["smb_tb_split"] = controllers(effects, boot)
    out["person_year"] = person_year(effects, boot)
    return out


__all__ = [
    "FAULT_ARMS",
    "GATE_OPTIONS",
    "LADDERS",
    "SCENARIOS",
    "combine",
    "controllers",
    "dose_response",
    "gate",
    "guards",
    "ladders",
    "mechanism",
    "pairing_efficiency",
    "pathway",
    "person_year",
    "run",
    "state_cells",
    "state_model",
]
