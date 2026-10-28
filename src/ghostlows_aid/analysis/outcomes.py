"""One outcome row per arm from its traces (memo §8), and the input contract the harness fills.

The harness hands over one `ArmTrace` per arm of a slot: the 1-min plasma truth, insulin by route
and rescue carbohydrate, the 5-min displayed CGM, what the controller was sent, IOB, COB and
oref0's reason per loop cycle, from `pre_min` minutes before the onset t0 to t0 + W. Everything is
scored on [t0, t0 + W); the history before t0 is identical across the arms of a slot by
construction (memo §2), is checked to be, and only warms up the onset trend and the alarm layer.

`slot_outcomes(traces)` turns the arms of one (controller, patient, seed, slot) into rows of the
outcome table (schema.py): the clean arm is needed for the fault-window `ref_` columns, the time
to re-converge and the onset state, which every arm of the slot shares.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from ghostlows.bench.events import LEVEL_THRESHOLDS, detect_events
from ghostlows.metrics.consensus import TAR_THRESHOLDS, TBR_THRESHOLDS
from ghostlows.sim.engine import DEFAULT_START, SAMPLE_MIN, VALIDITY_FLOOR

from . import schema
from .schema import ACTIONS, ARMS, CLEAN, WINDOW_MIN

# SOURCE: Klonoff 2023 (GRI = 3.0 x VLow<54 + 2.4 x Low 54-69 + 1.6 x VHigh>250 + 0.8 x High
# 181-250, here applied to minutes, memo §8) CLINICAL-REVIEW: CR-pending (memo §14: GRI weights
# as one harm scale)
GRI_WEIGHTS = {"below_54": 3.0, "54_69": 2.4, "above_250": 1.6, "181_250": 0.8}
RECONVERGE_MG_DL = 5.0  # ASSUMPTION: memo §8 (|truth - clean truth| < 5 mg/dL ...)
RECONVERGE_HOLD_MIN = 60  # ASSUMPTION: memo §8 (... for 60 min)
# ASSUMPTION: memo §4 (min_bg = max_bg = 100 mg/dL, so oref0's suspend threshold is
# min_bg - 0.5 (min_bg - 40) = 70 mg/dL, determine-basal.js threshold)
SUSPEND_THRESHOLD = 70.0
HIGH_INPUT = TAR_THRESHOLDS[0]  # SOURCE: Battelino 2019 (180 mg/dL), memo §8 crossings
TREND_MIN = 15  # ASSUMPTION: memo §8 ("15-min trend"): (value at t0 - value 15 min before) / 15
# ASSUMPTION: history before t0 for the onset trend and the alarm layer (the trend fit, the Baysal
# Kalman estimate and the 30-min refractory period warm up on it); identical across arms
PRE_MIN = 120
PERSON_BOLUS_EVENTS = ("bolus", "correction")  # SMBs are the controller's, not the person's

# oref0 0.7.1 rT.reason texts of the shipped rules (memo §8). A cycle counts once per rule.
REASON_PATTERNS = {
    # SOURCE: oref0 0.7.1 lib/determine-basal/determine-basal.js:188 (minAgo > 12 or < -5)
    "stale": re.compile(r"BG data is too old"),
    # SOURCE: determine-basal.js:177-194 (delta 0, |avgdeltas| < 1, device not fakecgm)
    "flat": re.compile(r"CGM data is unchanged"),
    # SOURCE: determine-basal.js:381 and :388 (pump temp disagrees with pump history)
    "temp_mismatch": re.compile(r"^Warning: currenttemp .*canceling temp"),
    # SOURCE: determine-basal.js:912-920 (bg or minGuardBG below threshold: 30-120 min zero temp)
    "lgs_zero_temp": re.compile(r"minGuardBG -?\d+(?:\.\d+)?<-?\d+(?:\.\d+)?"),
    # SOURCE: determine-basal.js:877-879 (maxDelta > 20% of BG); written whatever the SMB setting
    "smb_maxdelta": re.compile(r"maxDelta -?[\d.]+ > [\d.]+% of BG -?[\d.]+: SMB disabled"),
}
assert tuple(REASON_PATTERNS) == ACTIONS


@dataclass
class ArmTrace:
    """One arm of one slot of one trunk, as the harness hands it to the analysis.

    1-min series cover minutes [t0 - pre_min, t0 + window_min) of the trunk clock (00:00 on day 1
    = minute 0); index i is minute t0 - pre_min + i. 5-min series hold the samples at minutes
    t0 - pre_min + 5 j. NaN marks a missing reading. Units: mg/dL, U per minute, g per minute.
    """

    controller: str
    patient_id: str  # paper 1's form, e.g. "sim:adult#001"
    seed: int
    slot: str  # night, dose or any
    arm: str  # schema.ARMS key; "clean" = the trunk's own continuation
    onset_minute: int  # t0, a multiple of 5
    truth: np.ndarray  # 1-min plasma glucose
    basal_u: np.ndarray  # 1-min basal or temp basal delivered
    bolus_u: np.ndarray  # 1-min meal boluses of the person (wizard)
    smb_u: np.ndarray  # 1-min SMBs (zeros without SMB)
    rescue_g: np.ndarray  # 1-min rescue carbohydrate eaten
    cgm_display: np.ndarray  # 5-min reading on the person's display
    cgm_controller: np.ndarray  # 5-min reading sent to the controller (all NaN in open loop)
    iob_u: np.ndarray  # 5-min IOB at each loop cycle (NaN where unknown)
    cob_g: np.ndarray  # 5-min COB at each loop cycle (NaN where unknown)
    reasons: Sequence[str]  # 5-min oref0 rT.reason ("" where none; all "" in open loop)
    scheduled_basal_u_per_min: float
    last_bolus_minute: float  # minute of the person's last bolus before t0 (NaN: none)
    fault_end_min: float = (
        math.nan
    )  # minutes after t0 the fault stops changing readings (inf: held)
    fault: Mapping[str, float] = field(default_factory=dict)  # realised parameters
    correction_u: np.ndarray | None = None  # 1-min correction boluses (None: none given)
    valid: bool = True  # the harness's own validity flag for the run up to t0 + W
    pre_min: int = PRE_MIN
    window_min: int = WINDOW_MIN
    clock_offset_min: int = 0  # minute of day at the trunk start
    data_source: str = "simulation"

    def __post_init__(self):
        n1 = self.pre_min + self.window_min
        if self.pre_min < TREND_MIN or self.pre_min % SAMPLE_MIN or self.window_min % SAMPLE_MIN:
            raise ValueError("pre_min >= 15 and both spans must be multiples of 5 min")
        if self.correction_u is None:
            self.correction_u = np.zeros(n1)
        for name in ("truth", "basal_u", "bolus_u", "smb_u", "rescue_g", "correction_u"):
            v = np.asarray(getattr(self, name), dtype=float)
            if v.shape != (n1,):
                raise ValueError(f"{name}: {v.shape} != ({n1},)")
            setattr(self, name, v)
        n5 = n1 // SAMPLE_MIN
        for name in ("cgm_display", "cgm_controller", "iob_u", "cob_g"):
            v = np.asarray(getattr(self, name), dtype=float)
            if v.shape != (n5,):
                raise ValueError(f"{name}: {v.shape} != ({n5},)")
            setattr(self, name, v)
        if len(self.reasons) != n5:
            raise ValueError(f"reasons: {len(self.reasons)} != {n5}")
        if self.arm not in ARMS:
            raise ValueError(f"unknown arm {self.arm!r}")

    @property
    def w(self) -> slice:
        """The scored window [t0, t0 + W) on the 1-min grid."""
        return slice(self.pre_min, self.pre_min + self.window_min)

    @property
    def j0(self) -> int:
        """Index of the sample at t0 on the 5-min grid."""
        return self.pre_min // SAMPLE_MIN

    def fault_minutes(self) -> float:
        """Whole minutes of the window during the fault, ceil(fault end), capped at W; NaN with
        no fault."""
        d = float(self.fault_end_min)
        if math.isnan(d):
            return math.nan
        if d >= self.window_min:  # held for the window (drift, gain, interference), or longer
            return float(self.window_min)
        return float(max(math.ceil(d - 1e-9), 0))


# ---------------------------------------------------------------- per-arm outcomes


def glycaemia(truth: np.ndarray) -> dict:
    """Minutes in ranges, risk-weighted minutes, mean, peak and nadir on 1-min truth."""
    g = np.asarray(truth, float)
    lo2, lo1 = TBR_THRESHOLDS[1], TBR_THRESHOLDS[0]  # 54, 70 (strictly below)
    hi1, hi2 = TAR_THRESHOLDS  # 180, 250 (strictly above)
    vlow, low = np.sum(g < lo2), np.sum((g >= lo2) & (g < lo1))
    high, vhigh = np.sum((g > hi1) & (g <= hi2)), np.sum(g > hi2)
    w = GRI_WEIGHTS
    return {
        "minutes_below_54": float(vlow),
        "minutes_below_70": float(vlow + low),
        "minutes_above_180": float(high + vhigh),
        "minutes_above_250": float(vhigh),
        "risk_weighted_min": float(
            w["below_54"] * vlow + w["54_69"] * low + w["above_250"] * vhigh + w["181_250"] * high
        ),
        "mean_glucose": float(np.mean(g)),
        "peak_glucose": float(np.max(g)),
        "nadir_glucose": float(np.min(g)),
    }


def hypo_events(truth: np.ndarray, start: int) -> dict:
    """Level 1 and 2 events (paper 1's detect_events, Battelino 2023) on the truth with its
    history, counted when their onset is at or after index `start` (t0)."""
    out = {}
    for level, name in ((1, "hypo_l1_events"), (2, "hypo_l2_events")):
        ev = detect_events(truth, LEVEL_THRESHOLDS[level])
        out[name] = float(sum(s >= start for s, _ in ev))
    return out


def insulin(t: ArmTrace) -> dict:
    """Insulin by route, rescue carbohydrate and corrections over the window."""
    w = t.w
    basal, bolus, smb, corr = t.basal_u[w], t.bolus_u[w], t.smb_u[w], t.correction_u[w]
    return {
        "insulin_total_u": float(basal.sum() + bolus.sum() + smb.sum() + corr.sum()),
        "insulin_bolus_u": float(bolus.sum()),
        "insulin_correction_u": float(corr.sum()),
        "insulin_smb_u": float(smb.sum()),
        "insulin_tb_dev_u": float(basal.sum() - t.scheduled_basal_u_per_min * basal.size),
        "rescue_carbs_g": float(t.rescue_g[w].sum()),
        "rescue_count": float(np.count_nonzero(t.rescue_g[w] > 0)),
        "correction_boluses": float(np.count_nonzero(corr > 0)),
    }


def _total_insulin(t: ArmTrace) -> np.ndarray:
    return (t.basal_u + t.bolus_u + t.smb_u + t.correction_u)[t.w]


def fault_windows(t: ArmTrace, clean: ArmTrace | None) -> dict:
    """Insulin and minutes above 180 during [t0, fault end) and after it, for the arm and (ref_)
    for the clean arm over the same minutes. NaN without a fault or a clean arm."""
    keys = ("insulin_during_u", "insulin_after_u")
    keys += ("minutes_above_180_during", "minutes_above_180_after")
    out = {k: math.nan for k in keys} | {schema.ref_of(k): math.nan for k in keys}
    d = t.fault_minutes()
    if math.isnan(d):
        return out
    d = int(d)
    for prefix, trace in (("", t), ("ref_", clean)):
        if trace is None:
            continue
        ins, high = _total_insulin(trace), trace.truth[trace.w] > HIGH_INPUT
        out[f"{prefix}insulin_during_u"] = float(ins[:d].sum())
        out[f"{prefix}minutes_above_180_during"] = float(high[:d].sum())
        if d < t.window_min:  # D9: a fault active at t0 + W is censored; nothing after it
            out[f"{prefix}insulin_after_u"] = float(ins[d:].sum())
            out[f"{prefix}minutes_above_180_after"] = float(high[d:].sum())
    return out


def reconvergence(t: ArmTrace, clean: ArmTrace | None) -> dict:
    """Minutes from the fault end to the first 60-min run with |truth - clean truth| < 5 mg/dL.

    NaN and `reconverged` False when no such run starts in the window (censored); both NA when
    there is no fault, no clean arm, or the fault lasts the whole window."""
    na = {"reconverge_min": math.nan, "reconverged": pd.NA}
    d = t.fault_minutes()
    if clean is None or math.isnan(d) or d >= t.window_min:
        return na
    ok = np.abs(t.truth[t.w] - clean.truth[clean.w]) < RECONVERGE_MG_DL
    hold = RECONVERGE_HOLD_MIN
    c = np.r_[0, np.cumsum(ok)]
    starts = np.flatnonzero(c[hold:] - c[:-hold] == hold)  # runs [s, s + hold) all within 5
    starts = starts[starts >= d]
    if starts.size == 0:
        return {"reconverge_min": math.nan, "reconverged": False}
    return {"reconverge_min": float(starts[0] - d), "reconverged": True}


def classify_reasons(reasons: Sequence[str]) -> dict:
    """Number of cycles whose reason shows each shipped rule (REASON_PATTERNS)."""
    return {
        f"actions_{k}": float(sum(bool(p.search(r or "")) for r in reasons))
        for k, p in REASON_PATTERNS.items()
    }


def actions(t: ArmTrace) -> dict:
    """Controller-action counts over the loop cycles in the window (NaN in open loop)."""
    if t.controller == "open_loop":
        return {f"actions_{k}": math.nan for k in REASON_PATTERNS}
    j0 = t.j0
    return classify_reasons(list(t.reasons)[j0 : j0 + t.window_min // SAMPLE_MIN])


def crossings(t: ArmTrace) -> dict:
    """Whether a reading sent to the controller during the fault fell below the suspend threshold
    (70) or rose above 180 mg/dL (memo §8). NA in open loop and without a fault."""
    na = {"ctrl_below_70_during": pd.NA, "ctrl_above_180_during": pd.NA}
    d = t.fault_minutes()
    if t.controller == "open_loop" or math.isnan(d):
        return na
    n = max(1, math.ceil(d / SAMPLE_MIN))  # samples at t0, t0 + 5, ... before the fault end
    x = t.cgm_controller[t.j0 : t.j0 + n]
    x = x[~np.isnan(x)]
    return {
        "ctrl_below_70_during": bool((x < SUSPEND_THRESHOLD).any()),
        "ctrl_above_180_during": bool((x > HIGH_INPUT).any()),
    }


def onset_state(clean: ArmTrace) -> dict:
    """The state at t0 from the slot's clean arm (memo §8): true glucose and its 15-min trend,
    the fault-free displayed CGM and its trend, IOB, COB, hours since the person's last bolus and
    the hour of day."""
    i, j, k = clean.pre_min, clean.j0, TREND_MIN // SAMPLE_MIN
    g, cgm = clean.truth, clean.cgm_display
    since = (clean.onset_minute - clean.last_bolus_minute) / 60.0
    return {
        "onset_truth_mg_dl": float(g[i]),
        "onset_truth_trend": float((g[i] - g[i - TREND_MIN]) / TREND_MIN),
        "onset_cgm_mg_dl": float(cgm[j]),
        "onset_cgm_trend": float((cgm[j] - cgm[j - k]) / TREND_MIN),
        "onset_iob_u": float(clean.iob_u[j]),
        "onset_cob_g": float(clean.cob_g[j]),
        "onset_hours_since_bolus": float(since),
        "onset_hour": ((clean.clock_offset_min + clean.onset_minute) % 1440) / 60.0,
    }


# realised fault parameters as the harness names them (faults_online.SingleFault.params)
_PARAM_KEYS = {
    "fault_depth_mg_dl": ("depth",),
    "fault_duration_min": ("duration",),
    "fault_size": ("size", "size_mg_dl", "sd", "gain", "peak_mg_dl"),
    "fault_tau_min": ("tau", "tau_min"),
}


def fault_columns(t: ArmTrace) -> dict:
    """Realised fault parameters; the PISA area is depth x duration in mg/dL x h (memo §10)."""
    p = dict(t.fault)
    out = {"onset_minute": int(t.onset_minute), "fault_end_min": float(t.fault_end_min)}
    for col, names in _PARAM_KEYS.items():
        out[col] = next((float(p[n]) for n in names if n in p), math.nan)
    out["fault_area_mg_dl_h"] = out["fault_depth_mg_dl"] * out["fault_duration_min"] / 60.0
    out["fault_params_json"] = json.dumps(p, sort_keys=True, default=float)
    return out


def is_invalid(t: ArmTrace) -> bool:
    """Outside model validity: simulated glucose below paper 1's 10 mg/dL floor, or the harness
    flagged the run (memo §10)."""
    return bool((not t.valid) or np.nanmin(t.truth) < VALIDITY_FLOOR)


def arm_outcomes(t: ArmTrace, clean: ArmTrace | None) -> dict:
    """The outcome row of one arm (alarm-layer columns excluded; see slot_outcomes)."""
    row = schema.empty_row(t.arm)
    row.update(
        controller=t.controller,
        patient_id=t.patient_id,
        seed=int(t.seed),
        slot=t.slot,
        data_source=t.data_source,
        invalid=is_invalid(t),
    )
    row.update(fault_columns(t))
    row.update(glycaemia(t.truth[t.w]))
    row.update(hypo_events(t.truth[: t.pre_min + t.window_min], t.pre_min))
    row.update(insulin(t))
    row.update(fault_windows(t, clean))
    row.update(reconvergence(t, clean))
    row.update(actions(t))
    row.update(crossings(t))
    return row


def _same_history(a: ArmTrace, b: ArmTrace) -> bool:
    n, k = a.pre_min, a.j0
    return (
        a.pre_min == b.pre_min
        and a.window_min == b.window_min
        and a.onset_minute == b.onset_minute
        and np.array_equal(a.truth[:n], b.truth[:n], equal_nan=True)
        and np.array_equal(a.cgm_display[:k], b.cgm_display[:k], equal_nan=True)
    )


def slot_outcomes(traces: Sequence[ArmTrace], alarm_layer: bool = True) -> pd.DataFrame:
    """Outcome rows for the arms of one (controller, patient, seed, slot); one must be clean.

    Checks that every arm shares the clean arm's identity, onset and history before t0 (pairing,
    memo §2). The onset state comes from the clean arm. With `alarm_layer`, paper 1's alarms are
    scored on every arm (alarm_layer.py).
    """
    clean = [t for t in traces if t.arm == CLEAN]
    if len(clean) != 1:
        raise ValueError("a slot needs exactly one clean arm")
    clean = clean[0]
    ident = ("controller", "patient_id", "seed", "slot", "data_source")
    for t in traces:
        if any(getattr(t, k) != getattr(clean, k) for k in ident):
            raise ValueError(f"arm {t.arm} belongs to another unit than the clean arm")
        if not _same_history(t, clean):
            raise ValueError(f"arm {t.arm}: history before t0 differs from the clean arm's")
    state = onset_state(clean)
    rows = []
    for t in traces:
        row = arm_outcomes(t, None if t.arm == CLEAN else clean)
        row.update(state)
        rows.append(row)
    if alarm_layer:
        from .alarm_layer import alarm_columns

        for row, cols in zip(rows, alarm_columns(traces), strict=True):
            row.update(cols)
    return schema.coerce(pd.DataFrame(rows))


# ---------------------------------------------------------------- adapter for the harness


def trace_from_loop(
    result,
    *,
    controller: str,
    patient_id: str,
    seed: int,
    slot: str,
    arm: str,
    onset_minute: int,
    scheduled_basal_u_per_min: float,
    fault: Mapping | None = None,
    fault_end_min: float | None = None,
    pre_min: int = PRE_MIN,
    window_min: int = WINDOW_MIN,
    clock_offset_min: int = 0,
    data_source: str = "simulation",
) -> ArmTrace:
    """An ArmTrace from a ghostlows_aid.loop.LoopResult that spans [t0 - pre_min, t0 + W).

    Reads the 1-min truth, `insulin_basal`, `insulin_bolus` (the person's meal and correction
    boluses), `insulin_smb`, the `correction` and `rescue` events, and the 5-min samples
    (`cgm_person`; in closed loop `cgm_sent`, the reading the controller received after any
    gate, and `iob`, `cob`, `rT_reason`). `fault` is SingleFault.describe(); the fault's end
    relative to t0 defaults to its `end - onset` (inf for a held fault).
    """
    lo, hi = onset_minute - pre_min, onset_minute + window_min
    if lo < 0 or hi > len(result.truth):
        raise ValueError(f"the result covers {len(result.truth)} min, need [{lo}, {hi})")
    ev = result.events
    in_span = (ev["minute"] >= lo) & (ev["minute"] < hi)

    def per_minute(kind: str) -> np.ndarray:
        e = ev[in_span & (ev["kind"] == kind)]
        out = np.zeros(hi - lo)
        np.add.at(out, e["minute"].to_numpy(int) - lo, e["value"].to_numpy(float))
        return out

    corr = per_minute("correction")
    bolus = np.asarray(result.insulin_bolus[lo:hi], float) - corr
    bolus[np.abs(bolus) < 1e-12] = 0.0
    s = result.samples.set_index("minute").reindex(np.arange(lo, hi, SAMPLE_MIN))

    def col(name: str) -> np.ndarray:
        return s[name].to_numpy(float) if name in s else np.full(len(s), np.nan)

    reasons = s["rT_reason"] if "rT_reason" in s else pd.Series("", index=s.index)
    person = ev[(ev["kind"].isin(PERSON_BOLUS_EVENTS)) & (ev["minute"] < onset_minute)]
    person = person[person["value"] > 0]
    # a trunk (the clean arm) runs on past t0 + W: its validity is judged up to t0 + W only
    valid = bool(result.valid)
    if len(result.truth) > hi:
        valid = bool(np.nanmin(np.asarray(result.truth[:hi], float)) >= VALIDITY_FLOOR)
    fault = dict(fault or {})
    if fault_end_min is None:
        end, start = fault.get("end", math.nan), fault.get("onset", onset_minute)
        fault_end_min = float(end) - float(start) if fault else math.nan
    return ArmTrace(
        controller=controller,
        patient_id=patient_id,
        seed=seed,
        slot=slot,
        arm=arm,
        onset_minute=onset_minute,
        truth=np.asarray(result.truth[lo:hi], float),
        basal_u=np.asarray(result.insulin_basal[lo:hi], float),
        bolus_u=bolus,
        smb_u=np.asarray(result.insulin_smb[lo:hi], float),
        rescue_g=per_minute("rescue"),
        correction_u=corr,
        cgm_display=col("cgm_person"),
        cgm_controller=col("cgm_sent"),
        iob_u=col("iob"),
        cob_g=col("cob"),
        reasons=[r if isinstance(r, str) else "" for r in reasons],
        scheduled_basal_u_per_min=float(scheduled_basal_u_per_min),
        last_bolus_minute=float(person["minute"].max()) if len(person) else math.nan,
        fault_end_min=fault_end_min,
        fault=fault,
        valid=valid,
        pre_min=pre_min,
        window_min=window_min,
        clock_offset_min=clock_offset_min,
        data_source=data_source,
    )


def trace_from_run(unit, run, *, pre_min: int = PRE_MIN, window_min: int = WINDOW_MIN) -> ArmTrace:
    """An ArmTrace from the harness's episodes.Unit and one of its episodes.ArmRun.

    The harness's `clean@<slot>` arm becomes the table's `clean`; its patient `adult#001` becomes
    `sim:adult#001`; the clock starts at paper 1's DEFAULT_START, as the loop's does.
    """
    patient = unit.patient if unit.patient.startswith("sim:") else f"sim:{unit.patient}"
    return trace_from_loop(
        run.result,
        controller=unit.controller,
        patient_id=patient,
        seed=int(unit.seed),
        slot=run.arm.slot,
        arm=schema.table_arm(run.arm.name),
        onset_minute=int(run.t0),
        scheduled_basal_u_per_min=float(unit.scheduled_basal_u_per_min),
        fault=run.fault,
        pre_min=pre_min,
        window_min=window_min,
        clock_offset_min=DEFAULT_START.hour * 60 + DEFAULT_START.minute,
    )
