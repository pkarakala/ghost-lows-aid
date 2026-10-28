"""The per-arm outcome table: arms, keys, columns, units and dtypes in one place, and a validator.

One row per (controller, patient_id, seed, slot, arm) (memo §8, "stored per unit"). Arms are named
in `ARMS` (memo §6), with the harness's catalogue names (`ghostlows_aid.episodes.ARMS`). The
fault-free arm of every slot is `clean` (the harness's `clean@<slot>`, see `table_arm`): the
trunk's own continuation from that slot's onset t0. A variant arm is named `<base>+<variant>`, for
example `pisa_typical+controller_only` or `clean+baysal_trial`.

How a column's paired effect is formed (`Column.effect`, used by `effects.paired_effects`):

- `diff`: the arm's row minus the clean row of the same (controller, patient_id, seed, slot).
- `ref`: the column minus its `ref_` twin in the same row, which holds the clean arm's value over
  this arm's fault window (insulin and hyperglycaemia during and after the fault).
- `paired`: already a paired quantity (time to re-converge with the clean arm).
- `none`: keys, design, realised fault parameters, onset state and descriptive flags.

Counts are stored as float64 so that "not applicable" can be NaN (controller actions in open loop).
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- design (memo §4, §6)

CONTROLLERS = ("oref0_tb", "oref1_smb", "open_loop")  # memo §4
PRIMARY_CONTROLLER = "oref0_tb"
SLOTS = ("night", "dose", "any")  # memo §6
CLEAN = "clean"
DATA_SOURCES = ("synthetic", "simulation")
WINDOW_MIN = 720  # ASSUMPTION: memo §2 (W = 12 h)
PATIENT_ID = re.compile(r"^sim:(adolescent|adult|child)#\d{3}$")  # paper 1's simulated IDs


@dataclass(frozen=True)
class Arm:
    """One arm of memo §6. `base` is the arm a variant modifies (its own name for a fault arm)."""

    name: str
    slot: str  # "" for the clean arm, which exists in every slot
    fault: str  # none, pisa, dropout, frozen, spikes, drift, gain, interference
    level: str  # none, typical, typical_g6, severe, stress, small, large
    tier: int
    controllers: tuple[str, ...] = CONTROLLERS
    direction: str = ""  # low or high for drift, gain and interference
    base: str = ""
    variant: str = ""  # controller_only, person_only, baysal_<set>, matched_filter, fakecgm, ...


_FAULTS = (
    Arm("pisa_typical", "night", "pisa", "typical", 1),
    Arm("pisa_severe", "night", "pisa", "severe", 2),
    Arm("dropout_short", "any", "dropout", "typical", 1),
    Arm("dropout_long", "any", "dropout", "typical_g6", 2),
    Arm("dropout_severe", "any", "dropout", "severe", 2),
    Arm("frozen_typical", "any", "frozen", "typical", 1),
    Arm("frozen_severe", "any", "frozen", "severe", 2),
    Arm("spike_typical", "any", "spikes", "typical", 1),
    Arm("spike_severe", "any", "spikes", "severe", 2),
    Arm("drift_low", "any", "drift", "typical", 1, direction="low"),
    Arm("drift_low_severe", "any", "drift", "severe", 2, direction="low"),
    Arm("drift_high", "any", "drift", "typical", 1, direction="high"),
    Arm("drift_high_severe", "any", "drift", "severe", 2, direction="high"),
    Arm("gain_low", "any", "gain", "stress", 2, direction="low"),
    Arm("gain_high", "any", "gain", "stress", 2, direction="high"),
    Arm("interference_small", "dose", "interference", "small", 2, direction="high"),
    Arm("interference_large", "dose", "interference", "large", 1, direction="high"),
)
FAULT_ARMS = tuple(a.name for a in _FAULTS)
GATE_SETTINGS = ("aggressive", "trial", "nominal", "cautious")  # paper 1's Baysal 2014 sets
GATES = tuple(f"baysal_{s}" for s in GATE_SETTINGS)  # the harness's detector names
SECOND_DETECTOR = "matched_filter"  # paper 1's causal matched filter (Idi 2023 adaptation)


def _variant(base: str, variant: str, tier: int) -> Arm:
    b = next((a for a in _FAULTS if a.name == base), None)
    slot = {"clean_night": "night", "clean_dose": "dose"}.get(base, b.slot if b else "")
    name = base.split("_")[0] if base.startswith("clean") else base
    return Arm(
        f"{name}+{variant}",
        slot,
        b.fault if b else "none",
        b.level if b else "none",
        tier,
        controllers=(PRIMARY_CONTROLLER,),
        direction=b.direction if b else "",
        base=name,
        variant=variant,
    )


_VARIANTS = (
    # pathway (Q3, tier 2)
    _variant("pisa_typical", "controller_only", 2),
    _variant("pisa_typical", "person_only", 2),
    _variant("interference_large", "controller_only", 2),
    _variant("interference_large", "person_only", 2),
    # detector gate (Q5, tier 2) on the night-slot clean continuation and on PISA typical
    *(_variant("clean_night", g, 2) for g in GATES),
    *(_variant("pisa_typical", g, 2) for g in GATES),
    # second detector (Q5, tier 3)
    _variant("clean_night", SECOND_DETECTOR, 3),
    _variant("pisa_typical", SECOND_DETECTOR, 3),
    # guard ablation (Q5, tier 2): flat-CGM check bypassed by the device name
    _variant("frozen_severe", "fakecgm", 2),
    # correcting person (N7, tier 3)
    _variant("clean_dose", "correcting", 3),
    _variant("interference_large", "correcting", 3),
)
ARMS: dict[str, Arm] = {
    CLEAN: Arm(CLEAN, "", "none", "none", 1),
    **{a.name: a for a in _FAULTS},
    **{a.name: a for a in _VARIANTS},
}


def table_arm(name: str) -> str:
    """The outcome table's arm name for a harness arm name (`clean@night` -> `clean`)."""
    return CLEAN if name.startswith(f"{CLEAN}@") else name


def arms_for(controller: str, max_tier: int = 3) -> list[str]:
    """Arm names that run for `controller` up to `max_tier` (clean first)."""
    return [n for n, a in ARMS.items() if controller in a.controllers and a.tier <= max_tier]


# ---------------------------------------------------------------- columns


@dataclass(frozen=True)
class Column:
    name: str
    dtype: str  # str, int64, float64, bool, boolean (nullable)
    unit: str
    effect: str  # diff, ref, paired, none
    group: str
    doc: str


def _c(name, dtype, unit, effect, group, doc) -> Column:
    return Column(name, dtype, unit, effect, group, doc)


KEYS = ("controller", "patient_id", "seed", "slot", "arm")
UNIT_KEYS = ("controller", "patient_id", "seed", "slot")  # one trunk branch point

# Paper 1's alarm layer (alarm_layer.py maps each prefix to a ghostlows alarm).
ALARM_PREFIXES = ("trend", "threshold", *GATES, SECOND_DETECTOR)
ALARM_FIELDS = {
    "alarms": "alarms fired in the window (floor included)",
    "false_alarms": "alarms in the window matched to no true level 1 low (paper 1 §4)",
    "l1_caught": "true level 1 lows with onset in the window that the alarm caught",
    "l1_missed": "true level 1 lows with onset in the window that it missed",
    "l2_caught": "true level 2 lows with onset in the window that the alarm caught",
    "l2_missed": "true level 2 lows with onset in the window that it missed",
}
ACTIONS = ("stale", "flat", "temp_mismatch", "lgs_zero_temp", "smb_maxdelta")

_COLUMNS = [
    # keys
    _c("controller", "str", "", "none", "key", "controller (memo §4)"),
    _c("patient_id", "str", "", "none", "key", "virtual patient, paper 1's form sim:adult#001"),
    _c("seed", "int64", "", "none", "key", "main seed (memo §7)"),
    _c("slot", "str", "", "none", "key", "onset slot: night, dose or any (memo §6)"),
    _c("arm", "str", "", "none", "key", "arm name (ARMS); clean = the fault-free continuation"),
    _c("data_source", "str", "", "none", "key", "synthetic or simulation"),
    # design (from ARMS, stored so the table is self-describing)
    _c("fault_kind", "str", "", "none", "design", "fault type (ARMS)"),
    _c("fault_level", "str", "", "none", "design", "fault level (ARMS)"),
    _c("variant", "str", "", "none", "design", "variant of the base arm, '' for none"),
    # realised fault parameters
    _c("onset_minute", "int64", "min", "none", "fault", "t0, minutes from the trunk start"),
    _c(
        "fault_end_min", "float64", "min", "none", "fault", "fault end after t0; inf held; NaN none"
    ),
    _c("fault_depth_mg_dl", "float64", "mg/dL", "none", "fault", "PISA depth (peak attenuation)"),
    _c("fault_duration_min", "float64", "min", "none", "fault", "PISA duration or gap/hold/burst"),
    _c("fault_size", "float64", "", "none", "fault", "drift/gain fraction, bias mg/dL or spike SD"),
    _c("fault_tau_min", "float64", "min", "none", "fault", "PISA time constant"),
    _c("fault_area_mg_dl_h", "float64", "mg/dL*h", "none", "fault", "PISA depth x duration"),
    _c("fault_params_json", "str", "", "none", "fault", "every realised parameter, as given"),
    # validity (memo §10)
    _c("invalid", "bool", "", "none", "validity", "simulated glucose < 10 mg/dL (paper 1 floor)"),
    # glycaemia on 1-min plasma truth over [t0, t0 + W) (memo §8)
    _c("minutes_below_54", "float64", "min", "diff", "glycaemia", "minutes < 54 mg/dL"),
    _c("minutes_below_70", "float64", "min", "diff", "glycaemia", "minutes < 70 mg/dL"),
    _c("minutes_above_180", "float64", "min", "diff", "glycaemia", "minutes > 180 mg/dL"),
    _c("minutes_above_250", "float64", "min", "diff", "glycaemia", "minutes > 250 mg/dL"),
    _c("risk_weighted_min", "float64", "min", "diff", "glycaemia", "GRI-weighted minutes"),
    _c("mean_glucose", "float64", "mg/dL", "diff", "glycaemia", "mean truth"),
    _c("peak_glucose", "float64", "mg/dL", "diff", "glycaemia", "maximum truth"),
    _c("nadir_glucose", "float64", "mg/dL", "diff", "glycaemia", "minimum truth"),
    _c("hypo_l1_events", "float64", "count", "diff", "glycaemia", "level 1 events, onset in W"),
    _c("hypo_l2_events", "float64", "count", "diff", "glycaemia", "level 2 events, onset in W"),
    # insulin and carbohydrate over the window
    _c("insulin_total_u", "float64", "U", "diff", "insulin", "all insulin delivered"),
    _c("insulin_bolus_u", "float64", "U", "diff", "insulin", "the person's meal boluses"),
    _c("insulin_correction_u", "float64", "U", "diff", "insulin", "the person's corrections"),
    _c("insulin_smb_u", "float64", "U", "diff", "insulin", "SMBs"),
    _c("insulin_tb_dev_u", "float64", "U", "diff", "insulin", "basal delivered - scheduled"),
    _c("rescue_carbs_g", "float64", "g", "diff", "insulin", "rescue carbohydrate eaten"),
    _c("rescue_count", "float64", "count", "diff", "insulin", "rescue treatments"),
    _c("correction_boluses", "float64", "count", "diff", "insulin", "correction boluses given"),
    # during and after the fault, with the clean arm over the same minutes (ref_ twins)
    _c("insulin_during_u", "float64", "U", "ref", "window", "insulin in [t0, fault end)"),
    _c("ref_insulin_during_u", "float64", "U", "none", "window", "clean arm, same minutes"),
    _c("insulin_after_u", "float64", "U", "ref", "window", "insulin in [fault end, t0 + W)"),
    _c("ref_insulin_after_u", "float64", "U", "none", "window", "clean arm, same minutes"),
    _c("minutes_above_180_during", "float64", "min", "ref", "window", "> 180 during the fault"),
    _c("ref_minutes_above_180_during", "float64", "min", "none", "window", "clean, same minutes"),
    _c("minutes_above_180_after", "float64", "min", "ref", "window", "> 180 after the fault"),
    _c("ref_minutes_above_180_after", "float64", "min", "none", "window", "clean, same minutes"),
    _c("reconverge_min", "float64", "min", "paired", "window", "fault end to |truth - clean| < 5"),
    _c("reconverged", "boolean", "", "none", "window", "re-converged within the window"),
    # controller actions: loop cycles in the window whose reason shows the rule
    *(
        _c(f"actions_{a}", "float64", "count", "diff", "actions", f"cycles with rule {a}")
        for a in ACTIONS
    ),
    # threshold crossings of the controller's input during the fault (descriptive)
    _c("ctrl_below_70_during", "boolean", "", "none", "crossing", "input < 70 during the fault"),
    _c("ctrl_above_180_during", "boolean", "", "none", "crossing", "input > 180 during the fault"),
    # onset state at t0, from the slot's clean arm (identical for every arm of a slot)
    _c("onset_truth_mg_dl", "float64", "mg/dL", "none", "onset", "true glucose at t0"),
    _c("onset_truth_trend", "float64", "mg/dL/min", "none", "onset", "true 15-min trend"),
    _c("onset_cgm_mg_dl", "float64", "mg/dL", "none", "onset", "fault-free displayed CGM at t0"),
    _c("onset_cgm_trend", "float64", "mg/dL/min", "none", "onset", "its 15-min trend"),
    _c("onset_iob_u", "float64", "U", "none", "onset", "IOB at t0"),
    _c("onset_cob_g", "float64", "g", "none", "onset", "COB at t0"),
    _c("onset_hours_since_bolus", "float64", "h", "none", "onset", "since the last bolus"),
    _c("onset_hour", "float64", "h", "none", "onset", "hour of day at t0"),
    # paper 1's alarm layer on the arm's displayed CGM against its truth (memo §8)
    *(
        _c(f"{p}_{f}", "float64", "count", "diff", "alarm", f"{p}: {doc}")
        for p in ALARM_PREFIXES
        for f, doc in ALARM_FIELDS.items()
    ),
]
COLUMNS: dict[str, Column] = {c.name: c for c in _COLUMNS}
OUTCOMES = tuple(c.name for c in _COLUMNS if c.effect != "none")
ONSET = tuple(c.name for c in _COLUMNS if c.group == "onset")
FAULT_PARAMS = tuple(c.name for c in _COLUMNS if c.group == "fault")


def ref_of(name: str) -> str:
    """The `ref_` twin of a `ref`-effect column."""
    return f"ref_{name}"


def describe() -> pd.DataFrame:
    """The column dictionary as a table (for the protocol documents)."""
    return pd.DataFrame([vars(c) for c in _COLUMNS])


# ---------------------------------------------------------------- building and validating


class SchemaError(ValueError):
    pass


def coerce(table: pd.DataFrame) -> pd.DataFrame:
    """Cast the schema's columns to their dtypes and put them in schema order."""
    out = {}
    for c in _COLUMNS:
        s = table[c.name]
        if c.dtype == "str":
            out[c.name] = s.astype("str")
        elif c.dtype == "boolean":
            out[c.name] = s.astype("boolean")
        else:
            out[c.name] = s.astype(c.dtype)
    return pd.DataFrame(out, index=table.index)


def empty_row(arm: str) -> dict:
    """Defaults for one row: NaN outcomes, empty strings, the arm's design fields."""
    a = ARMS[arm]
    row = {}
    for c in _COLUMNS:
        row[c.name] = {"str": "", "int64": 0, "bool": False, "boolean": pd.NA}.get(c.dtype, np.nan)
    row.update(arm=arm, fault_kind=a.fault, fault_level=a.level, variant=a.variant)
    row["fault_params_json"] = json.dumps({})
    return row


def _dtype_ok(s: pd.Series, dtype: str) -> bool:
    t = pd.api.types
    return {
        "str": t.is_string_dtype(s) or (s.dtype == object and s.map(type).eq(str).all()),
        "int64": t.is_integer_dtype(s) and not t.is_bool_dtype(s),
        "float64": t.is_float_dtype(s),
        "bool": s.dtype == bool,
        "boolean": t.is_bool_dtype(s),
    }[dtype]


def validate(table: pd.DataFrame) -> pd.DataFrame:
    """Check the outcome table against the schema and the design; returns it in schema order.

    Raises SchemaError listing every problem: missing or unknown columns, wrong dtypes, duplicate
    keys, unknown controllers, slots, arms or data sources, arms in the wrong slot or for the
    wrong controller, design columns that disagree with ARMS, and units without a clean arm.
    """
    problems = []
    missing = [c for c in COLUMNS if c not in table]
    unknown = [c for c in table if c not in COLUMNS]
    if missing:
        problems.append(f"missing columns: {missing}")
    if unknown:
        problems.append(f"unknown columns: {unknown}")
    if missing:
        raise SchemaError("; ".join(problems))
    for c in _COLUMNS:
        if not _dtype_ok(table[c.name], c.dtype):
            problems.append(f"{c.name}: dtype {table[c.name].dtype}, expected {c.dtype}")
    if table.duplicated(list(KEYS)).any():
        problems.append(f"{int(table.duplicated(list(KEYS)).sum())} duplicate keys")
    for col, allowed in (
        ("controller", CONTROLLERS),
        ("slot", SLOTS),
        ("arm", tuple(ARMS)),
        ("data_source", DATA_SOURCES),
    ):
        bad = sorted(set(table[col].astype(str)) - set(allowed))
        if bad:
            problems.append(f"{col}: unknown values {bad}")
    bad_ids = sorted({p for p in table["patient_id"].astype(str) if not PATIENT_ID.match(p)})
    if bad_ids:
        problems.append(f"patient_id not in paper 1's form: {bad_ids[:5]}")
    known = table[table["arm"].isin(list(ARMS))]
    for arm, g in known.groupby("arm"):
        a = ARMS[arm]
        if a.slot and (g["slot"] != a.slot).any():
            problems.append(f"arm {arm} outside its slot {a.slot}")
        if not g["controller"].isin(a.controllers).all():
            problems.append(f"arm {arm} only runs for {a.controllers}")
        for col, want in (
            ("fault_kind", a.fault),
            ("fault_level", a.level),
            ("variant", a.variant),
        ):
            if (g[col] != want).any():
                problems.append(f"arm {arm}: {col} must be {want!r}")
    n_clean = (table["arm"] == CLEAN).groupby([table[k] for k in UNIT_KEYS]).sum()
    if (n_clean != 1).any():
        problems.append(
            f"{int((n_clean != 1).sum())} (controller, patient, seed, slot) units "
            "without exactly one clean arm"
        )
    if problems:
        raise SchemaError("; ".join(problems))
    return table[list(COLUMNS)]
