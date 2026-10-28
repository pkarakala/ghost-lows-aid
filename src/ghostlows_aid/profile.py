"""The oref0 profile of a simglucose virtual patient (TECHNICAL_DESIGN §6).

Python only decides the pump settings and preferences; oref0's own lib/profile (through the bridge,
wired as bin/oref0-get-profile.js) turns them into the profile, so every preference we do not set
keeps its lib/profile/index.js default.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from ghostlows.sim.engine import patient_params, quest
from ghostlows.sim.scenario import MEAL_AMOUNT_MU_G, MEAL_PROB

from .oref import Bridge, shared_bridge

ISF_RULES = ("quest", "rule1800", "rule1700")
CR_RULES = ("quest", "rule500")
MODES = ("oref0_tb", "oref1_smb")


@dataclass(frozen=True)
class ProfileSettings:
    """Settings that turn a virtual patient into an oref0 profile. All configurable."""

    # ASSUMPTION: Quest correction factor by default; 1800/TDI and 1700/TDI are the rule-of-thumb
    # alternatives (TECHNICAL_DESIGN §6 risk: Quest CF is far below 1800/TDI for some patients).
    # Chosen on calibration seeds and fixed by amendment (TECHNICAL_DESIGN §7.2).
    isf_rule: str = "quest"
    cr_rule: str = "quest"  # ASSUMPTION: Quest carb ratio, or 500/TDI
    basal_factor: float = 1.0  # paper 1's Therapy.basal_factor (D17); 1 = steady-state basal
    bolus_scale: float = 1.0  # carb ratio and ISF are multiplied by this (paper 1's D17 knob)
    # ASSUMPTION: 5 h insulin action, rapid-acting exponential curve (peak 75 min), matching paper
    # 1's IOB assumption (TECHNICAL_DESIGN §6); lib/iob/total.js:59-63 forces >= 5 h anyway.
    dia_h: float = 5.0
    curve: str = "rapid-acting"
    # ASSUMPTION: common OpenAPS target 100 mg/dL, min_bg = max_bg (TECHNICAL_DESIGN §6); oref0
    # forces max_bg = min_bg anyway (lib/profile/targets.js:31). Suspend threshold = 70 mg/dL.
    target: float = 100.0
    max_basal_multiple: float = 4.0  # ASSUMPTION: max_basal = 4 x scheduled basal (§6)
    # ASSUMPTION: max_iob = average meal bolus + 3 x max daily basal (OpenAPS documentation
    # guidance, TECHNICAL_DESIGN §6; a SOURCE key is still to be added). None = that rule.
    max_iob: float | None = None
    mode: str = "oref0_tb"  # oref0_tb: temp basals only; oref1_smb: SMB + UAM
    # ASSUMPTION: classic oref0 (AMA) has UAM off; oref1 turns on SMB with COB and UAM
    # (TECHNICAL_DESIGN §5.3 arms 1 and 2). lib/profile's own default is enableUAM true.
    enable_uam: bool | None = None  # None = by mode

    def __post_init__(self):
        if self.isf_rule not in ISF_RULES:
            raise ValueError(f"isf_rule {self.isf_rule!r} not in {ISF_RULES}")
        if self.cr_rule not in CR_RULES:
            raise ValueError(f"cr_rule {self.cr_rule!r} not in {CR_RULES}")
        if self.mode not in MODES:
            raise ValueError(f"mode {self.mode!r} not in {MODES}")


def scheduled_basal_u_per_min(name: str, basal_factor: float = 1.0) -> float:
    """Paper 1's basal exactly: basal_factor * u2ss * BW / 6000 U/min (ghostlows engine.py:133)."""
    p = patient_params(name)
    return basal_factor * float(p["u2ss"]) * float(p["BW"]) / 6000


def quest_tdi(name: str) -> float:
    """Total daily insulin. SOURCE: simglucose (Quest.csv)."""
    from ghostlows.sim.engine import _tables

    _, q, _ = _tables()
    return float(q.loc[q["Name"] == name, "TDI"].iloc[0])


def therapy_numbers(name: str, s: ProfileSettings) -> dict:
    """Basal (U/h), ISF (mg/dL/U), carb ratio (g/U) and max_iob for one patient."""
    basal = scheduled_basal_u_per_min(name, s.basal_factor) * 60
    cr_quest, cf_quest = quest(name)
    tdi = quest_tdi(name)
    isf = {"quest": cf_quest, "rule1800": 1800 / tdi, "rule1700": 1700 / tdi}[s.isf_rule]
    cr = {"quest": cr_quest, "rule500": 500 / tdi}[s.cr_rule]
    isf, cr = isf * s.bolus_scale, cr * s.bolus_scale
    if s.max_iob is None:
        # ASSUMPTION: "average meal bolus" = the mean of the main meals' mean sizes (slots eaten
        # with probability 0.95 in paper 1's scenario), divided by the carb ratio.
        main = [mu for mu, p in zip(MEAL_AMOUNT_MU_G, MEAL_PROB, strict=True) if p > 0.5]
        max_iob = float(np.mean(main)) / cr + 3 * basal
    else:
        max_iob = s.max_iob
    return {"basal_u_per_h": basal, "isf": isf, "carb_ratio": cr, "tdi": tdi, "max_iob": max_iob}


def pump_inputs(name: str, s: ProfileSettings) -> dict:
    """The inputs bin/oref0-get-profile.js assembles: preferences, pump settings, schedules."""
    n = therapy_numbers(name, s)
    smb = s.mode == "oref1_smb"
    uam = smb if s.enable_uam is None else s.enable_uam
    preferences = {
        "max_iob": n["max_iob"],
        "curve": s.curve,
        "enableUAM": uam,
        # ASSUMPTION: the oref1 arm turns on enableSMB_with_COB (TECHNICAL_DESIGN §5.3); the other
        # SMB triggers and maxSMBBasalMinutes 30, maxUAMSMBBasalMinutes 30, SMBInterval 3,
        # bolus_increment 0.1 keep the lib/profile/index.js:40-56 defaults.
        "enableSMB_with_COB": smb,
        "skip_neutral_temps": False,  # TECHNICAL_DESIGN §1.1: avoids the getMinutes() branch
    }
    start = {"i": 0, "start": "00:00:00", "offset": 0}
    return {
        "preferences": preferences,
        "settings": {
            "insulin_action_curve": s.dia_h,
            "maxBasal": s.max_basal_multiple * n["basal_u_per_h"],
        },
        "targets": {
            "units": "mg/dL",
            "user_preferred_units": "mg/dL",
            "targets": [{**start, "low": s.target, "high": s.target}],
        },
        "basals": [{"i": 0, "start": "00:00:00", "minutes": 0, "rate": n["basal_u_per_h"]}],
        "isf": {
            "units": "mg/dL",
            "user_preferred_units": "mg/dL",
            "sensitivities": [{**start, "sensitivity": n["isf"], "x": 0, "endOffset": 1440}],
        },
        "carbratio": {
            "units": "grams",
            "schedule": [{**start, "ratio": n["carb_ratio"], "x": 0, "r": 0}],
        },
    }


def build_profile(name: str, s: ProfileSettings, bridge: Bridge | None = None) -> dict:
    """The oref0 profile, generated by oref0's lib/profile (single-entry schedules, so its
    wall-clock schedule lookups cannot change the result)."""
    bridge = bridge or shared_bridge()
    return bridge.request({"op": "profile", "clock": 0, **pump_inputs(name, s)})["profile"]
