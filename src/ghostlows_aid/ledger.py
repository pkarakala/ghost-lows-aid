"""The pump and its delivery ledger, written as Nightscout-style treatments (TECHNICAL_DESIGN §3).

oref0 reads the ledger as pump history: "Temp Basal" records with rate, duration and created_at
(lib/iob/history.js:326-337) and "Meal Bolus"/"Correction Bolus" records with insulin
(history.js:295-303); carbs ride on the Meal Bolus record (lib/meal/history.js:85-95). A temp
that is superseded or cancelled early gets its actual duration, so oref0's lastTemp agrees with the
pump's currenttemp (determine-basal.js:369-391 cancels temps when they disagree).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .oref import NO_TEMP, iso

ENTERED_BY = "ghostlows-aid"


def js_round(x: float) -> float:
    """JavaScript Math.round for the non-negative values used here (halves round up)."""
    return float(math.floor(x + 0.5))


def round_basal(rate: float) -> float:
    """Temp rate the pump runs. SOURCE: oref0 lib/round-basal.js:11-39 (profile.model unset:
    0.05 U/h steps below 10 U/h, 0.1 U/h above); idempotent on oref0's own output."""
    if rate < 10:
        return js_round(rate * 20) / 20
    return js_round(rate * 10) / 10


def round_bolus(units: float, increment: float | None) -> float:
    """Wizard bolus on the pump's grid, to the nearest increment (None: no rounding)."""
    if increment is None:
        return units
    scale = round(1 / increment)
    return js_round(units * scale) / scale


@dataclass
class Temp:
    rate: float  # U/h
    start: int  # minute
    duration: int  # planned minutes
    record: dict  # its ledger record, whose duration is cut if it ends early


class Ledger:
    """Everything the pump delivered on top of its schedule, in time order."""

    def __init__(self, start_ms: int):
        self.start_ms = int(start_ms)
        self.records: list[tuple[int, dict]] = []  # (minute, Nightscout treatment)
        self.temp: Temp | None = None

    def ms(self, minute: int) -> int:
        return self.start_ms + int(minute) * 60000

    def _add(self, minute: int, record: dict) -> dict:
        record = {**record, "created_at": iso(self.ms(minute)), "enteredBy": ENTERED_BY}
        self.records.append((int(minute), record))
        return record

    # -- temp basals --
    def active_temp(self, minute: int) -> Temp | None:
        t = self.temp
        if t is not None and minute < t.start + t.duration:
            return t
        return None

    def currenttemp(self, minute: int) -> dict:
        """What the pump reports: rate and whole minutes left (Medtronic read_temp_basal)."""
        t = self.active_temp(minute)
        if t is None:
            return dict(NO_TEMP)
        return {"duration": t.start + t.duration - minute, "rate": t.rate, "temp": "absolute"}

    def _end_running(self, minute: int) -> None:
        t = self.active_temp(minute)
        if t is not None:
            actual = minute - t.start
            t.record["duration"] = actual  # the actual duration of a temp cut short
            if actual == 0:  # superseded in the minute it started: it never ran
                self.records = [(m, r) for m, r in self.records if r is not t.record]
        self.temp = None

    def set_temp(self, minute: int, rate: float, duration: int) -> float:
        """Start a temp (duration 0 cancels the running one); returns the rate the pump runs."""
        self._end_running(minute)
        if duration <= 0:
            return 0.0
        rate = round_basal(max(rate, 0.0))
        rec = self._add(
            minute, {"eventType": "Temp Basal", "rate": rate, "duration": int(duration)}
        )
        self.temp = Temp(rate, int(minute), int(duration), rec)
        return rate

    # -- boluses and carbs --
    def bolus(self, minute: int, units: float, carbs: float = 0.0, smb: bool = False) -> None:
        """A bolus (or a carb entry with insulin 0). SMBs are Correction Boluses (TECHNICAL_DESIGN
        §3); meal boluses carry the announced carbs."""
        if carbs > 0:
            self._add(minute, {"eventType": "Meal Bolus", "insulin": units, "carbs": carbs})
        elif units > 0:
            rec = {"eventType": "Correction Bolus", "insulin": units}
            if smb:
                rec["isSMB"] = True
            self._add(minute, rec)

    # -- what oref0 is sent --
    def treatments(self, minute: int, hours: float) -> list[dict]:
        """Records from the last `hours` up to `minute`, newest first: lib/iob/history.js:274-280
        drops any record newer than the one before it."""
        lo = minute - hours * 60
        return [r for m, r in reversed(self.records) if lo <= m <= minute]

    def boluses(self) -> list[tuple[int, float]]:
        return [(m, r["insulin"]) for m, r in self.records if r.get("insulin")]

    def temps(self) -> list[tuple[int, float, float]]:
        """(start minute, rate U/h, actual or planned duration) of every temp."""
        return [
            (m, r["rate"], r["duration"]) for m, r in self.records if r["eventType"] == "Temp Basal"
        ]
