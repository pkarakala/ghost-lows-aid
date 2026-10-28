"""Synthetic outcome tables with known effects, for the tests and the `--synthetic` dry run.

Every outcome follows an additive model per (controller, patient, seed, slot, arm):

    Y_clean = base + u_p + e_r            (patient level u_p, seed noise e_r)
    Y_arm   = Y_clean + theta + b_p + eps  (+ effect_fn(rows), if given)

with theta the known per-arm effect (`Design.effects`), b_p patient heterogeneity of the effect
and eps the arm's own seed noise. So the population estimand of memo §2 has expectation theta, and
with heterogeneity and noise set to 0 the estimate equals theta exactly. Fault-window outcomes
(`ref_` twins) and the time to re-converge follow the same model inside the arm's row. The
alarm-layer caught counts are the alarm's sensitivity times the true lows, so sensitivities are
known exactly. Patient IDs are paper 1's (`sim:adolescent#001` ... `sim:child#010`).

Synthetic tables carry data_source = "synthetic"; they are never simulation output.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from ghostlows.faults import pisa as pisa_model
from ghostlows.sim.engine import PATIENTS as PAPER1_PATIENTS

from . import schema
from .schema import ALARM_PREFIXES, ARMS, CLEAN, COLUMNS, CONTROLLERS, arms_for

PATIENTS = tuple(f"sim:{p}" for p in PAPER1_PATIENTS)  # 10 per age group

# (base, scale) per outcome: the clean level and the size of patient and seed variation
_LEVELS = {
    "minutes_below_54": (3.0, 3.0),
    "minutes_below_70": (15.0, 10.0),
    "minutes_above_180": (150.0, 60.0),
    "minutes_above_250": (30.0, 20.0),
    "risk_weighted_min": (160.0, 60.0),
    "mean_glucose": (140.0, 15.0),
    "peak_glucose": (230.0, 30.0),
    "nadir_glucose": (80.0, 10.0),
    "hypo_l1_events": (1.0, 0.1),
    "hypo_l2_events": (0.5, 0.05),
    "insulin_total_u": (18.0, 5.0),
    "insulin_bolus_u": (8.0, 3.0),
    "insulin_correction_u": (0.0, 0.5),
    "insulin_smb_u": (0.0, 0.5),
    "insulin_tb_dev_u": (0.0, 1.0),
    "rescue_carbs_g": (10.0, 8.0),
    "rescue_count": (0.7, 0.5),
    "correction_boluses": (0.0, 0.3),
    "insulin_during_u": (2.0, 1.0),
    "insulin_after_u": (16.0, 4.0),
    "minutes_above_180_during": (10.0, 8.0),
    "minutes_above_180_after": (140.0, 50.0),
    "reconverge_min": (60.0, 30.0),
    **{f"actions_{a}": (2.0, 1.0) for a in schema.ACTIONS},
    **{f"{p}_alarms": (1.0, 0.3) for p in ALARM_PREFIXES},
    **{f"{p}_false_alarms": (0.4, 0.2) for p in ALARM_PREFIXES},
}
# alarm sensitivities (level 1, level 2): caught = sensitivity x true lows
DEFAULT_SENSITIVITY = {p: (0.8, 0.9) for p in ALARM_PREFIXES}
_HOUR = {"night": (23.0, 31.0), "dose": (7.0, 15.0), "any": (0.0, 24.0)}  # memo §6 slots


@dataclass
class Design:
    """What to generate. `effects` maps (controller, arm, outcome) to theta; "*" as the
    controller applies to every controller. Noise sizes are fractions of each outcome's scale."""

    n_seeds: int = 6
    controllers: tuple[str, ...] = CONTROLLERS
    patients: tuple[str, ...] = PATIENTS
    max_tier: int = 3
    arms: tuple[str, ...] | None = None  # restrict to these arms (clean always kept)
    effects: Mapping[tuple[str, str, str], float] = field(default_factory=dict)
    patient_sd: float = 0.3
    heterogeneity_sd: float = 0.1
    seed_sd: float = 0.2
    sensitivity: Mapping[str, tuple[float, float]] = field(
        default_factory=lambda: dict(DEFAULT_SENSITIVITY)
    )
    invalid: Mapping[tuple, bool] = field(default_factory=dict)  # full keys -> True
    effect_fn: Callable[[pd.DataFrame], Mapping[str, np.ndarray]] | None = None
    zero_noise_arms: bool = False  # arms identical to their clean arm (before theta)
    first_seed: int = 101
    seed: int = 0


def _theta(d: Design, controller: str, arm: str, outcome: str, default: float = 0.0) -> float:
    for key in ((controller, arm, outcome), ("*", arm, outcome)):
        if key in d.effects:
            return float(d.effects[key])
    return default


def _units(d: Design) -> pd.DataFrame:
    rows = []
    for c in d.controllers:
        for p in d.patients:
            for r in range(d.n_seeds):
                for slot in schema.SLOTS:
                    rows.append((c, p, d.first_seed + r, slot))
    return pd.DataFrame(rows, columns=list(schema.UNIT_KEYS))


def _onset(units: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    n = len(units)
    truth = rng.normal(130, 35, n).clip(45, 350)
    trend = rng.normal(0, 1.0, n)
    lo = units["slot"].map(lambda s: _HOUR[s][0]).to_numpy()
    hi = units["slot"].map(lambda s: _HOUR[s][1]).to_numpy()
    hour = rng.uniform(lo, hi) % 24
    return pd.DataFrame(
        {
            "onset_truth_mg_dl": truth,
            "onset_truth_trend": trend,
            "onset_cgm_mg_dl": truth + rng.normal(0, 10, n),
            "onset_cgm_trend": trend + rng.normal(0, 0.5, n),
            "onset_iob_u": np.abs(rng.normal(1.0, 1.0, n)),
            "onset_cob_g": rng.exponential(15.0, n),
            "onset_hours_since_bolus": rng.uniform(0.5, 10.0, n),
            "onset_hour": hour,
            "onset_minute": (1440 + np.round(hour * 12) * 5).astype(np.int64),
        },
        index=units.index,
    )


def _fault_params(arm: str, n: int, rng: np.random.Generator) -> dict:
    """Realised parameters of the base fault, from paper 1's PISA model where it applies."""
    a = ARMS[arm]
    nan = np.full(n, np.nan)
    out = {
        "fault_end_min": nan.copy(),
        "fault_depth_mg_dl": nan.copy(),
        "fault_duration_min": nan.copy(),
        "fault_size": nan.copy(),
        "fault_tau_min": nan.copy(),
    }
    if a.fault == "none":
        return out
    if a.fault == "pisa":
        if a.level == "severe":
            depth = np.full(n, pisa_model.DEPTH.p95)
            dur = np.full(n, pisa_model.DURATION.p95)
        else:
            depth, dur = pisa_model.sample_sizes(rng, n)
        out.update(
            fault_depth_mg_dl=depth,
            fault_duration_min=dur,
            fault_tau_min=np.full(n, pisa_model.TAU_MIN),
            fault_end_min=dur + pisa_model.TAU_MIN * math.log(10.0),  # 1% tail after pressure
        )
    elif a.fault in ("dropout", "frozen", "spikes"):
        dur = {"severe": 120.0}.get(a.level, 20.0) + rng.uniform(-5, 5, n)
        out.update(fault_duration_min=dur, fault_end_min=dur)
    else:  # drift, gain, interference: held for the window
        size = {"drift": 0.1, "gain": 0.4, "interference": 40.0}[a.fault]
        out.update(
            fault_size=np.full(n, size * (1 if a.direction == "high" else -1)),
            fault_end_min=np.full(n, math.inf),
        )
    return out


def _defaults(arm: str, n: int) -> dict:
    """Column arrays of n empty rows of `arm` (schema.empty_row)."""
    out = {}
    for k, v in schema.empty_row(arm).items():
        if COLUMNS[k].dtype == "boolean":
            out[k] = pd.array([pd.NA] * n, dtype="boolean")
        else:
            out[k] = np.full(n, v, dtype=object if isinstance(v, str) else None)
    return out


def synthetic_outcomes(d: Design | None = None) -> pd.DataFrame:
    """A validated outcome table drawn from `d` (data_source "synthetic")."""
    d = d or Design()
    rng = np.random.default_rng(d.seed)
    units = _units(d)
    units = pd.concat([units, _onset(units, rng)], axis=1)
    outcome_cols = [c for c in schema.OUTCOMES if c in _LEVELS]
    ref_cols = [c for c in outcome_cols if COLUMNS[c].effect == "ref"]
    # clean level of each unit: base + patient + seed noise
    pid = {p: i for i, p in enumerate(d.patients)}
    p_idx = units["patient_id"].map(pid).to_numpy()
    clean = {}
    for o in outcome_cols:
        base, scale = _LEVELS[o]
        u = rng.normal(0, d.patient_sd * scale, len(d.patients))
        clean[o] = base + u[p_idx] + rng.normal(0, d.seed_sd * scale, len(units))
    frames = []
    for c in d.controllers:
        names = arms_for(c, d.max_tier)
        if d.arms is not None:
            names = [a for a in names if a == CLEAN or a in d.arms]
        for arm in names:
            a = ARMS[arm]
            sel = (units["controller"] == c).to_numpy()
            if a.slot:
                sel = sel & (units["slot"] == a.slot).to_numpy()
            u = units[sel].reset_index(drop=True)
            n = len(u)
            f = _defaults(arm, n)
            for k in (*schema.UNIT_KEYS, *schema.ONSET, "onset_minute"):
                f[k] = u[k].to_numpy()
            f["data_source"] = np.full(n, "synthetic", dtype=object)
            f.update(_fault_params(arm, n, rng))
            f["fault_area_mg_dl_h"] = f["fault_depth_mg_dl"] * f["fault_duration_min"] / 60.0
            ends_in_window = bool(np.isfinite(f["fault_end_min"]).all())
            has_fault = a.fault != "none"
            pi = u["patient_id"].map(pid).to_numpy()
            extra = d.effect_fn(pd.DataFrame(f)) if (d.effect_fn and arm != CLEAN) else {}
            for o in outcome_cols:
                base, scale = _LEVELS[o]
                y0 = clean[o][sel]
                if o in ref_cols:
                    f[schema.ref_of(o)] = y0 if has_fault else np.full(n, np.nan)
                if arm == CLEAN and COLUMNS[o].effect != "diff":
                    continue  # no fault window, nothing to re-converge from
                if COLUMNS[o].effect in ("ref", "paired") and not has_fault:
                    continue
                if o == "reconverge_min" and not ends_in_window:
                    continue
                if arm == CLEAN or d.zero_noise_arms:
                    het = noise = 0.0
                else:
                    het = rng.normal(0, d.heterogeneity_sd * scale, len(d.patients))[pi]
                    noise = rng.normal(0, d.seed_sd * scale, n)
                # a paired outcome is its own effect: theta defaults to the base level
                paired = COLUMNS[o].effect == "paired"
                theta = 0.0 if arm == CLEAN else _theta(d, c, arm, o, base if paired else 0.0)
                start = 0.0 if paired else y0
                f[o] = start + theta + het + noise + np.asarray(extra.get(o, 0.0))
            if has_fault and ends_in_window:
                f["reconverged"] = pd.array([True] * n, dtype="boolean")
            for prefix, (s1, s2) in d.sensitivity.items():
                f[f"{prefix}_l1_caught"] = s1 * f["hypo_l1_events"]
                f[f"{prefix}_l1_missed"] = (1 - s1) * f["hypo_l1_events"]
                f[f"{prefix}_l2_caught"] = s2 * f["hypo_l2_events"]
                f[f"{prefix}_l2_missed"] = (1 - s2) * f["hypo_l2_events"]
            if c == "open_loop":
                for act in schema.ACTIONS:
                    f[f"actions_{act}"] = np.full(n, np.nan)
            elif has_fault and arm != CLEAN:
                f["ctrl_below_70_during"] = pd.array(rng.random(n) < 0.3, dtype="boolean")
                f["ctrl_above_180_during"] = pd.array(rng.random(n) < 0.3, dtype="boolean")
            f = pd.DataFrame(f)
            frames.append(f)
    table = pd.concat(frames, ignore_index=True)
    if d.invalid:
        keys = list(zip(*(table[k] for k in schema.KEYS), strict=True))
        table["invalid"] = [bool(d.invalid.get(k, False)) for k in keys]
    return schema.validate(schema.coerce(table))
