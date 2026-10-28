"""The pump ledger: Nightscout treatments as oref0 reads them, temps with their actual duration."""

from __future__ import annotations

import json
import shutil
import subprocess

from ghostlows_aid.ledger import Ledger, round_basal, round_bolus
from ghostlows_aid.oref import NO_TEMP, NODE_DIR, iso

from .conftest import requires_node

START = 1767571200000  # 2026-01-05T00:00:00Z


def test_temp_is_rounded_and_stamped():
    led = Ledger(START)
    assert led.set_temp(10, 1.2735, 30) == 1.25
    ((minute, rec),) = led.records
    assert minute == 10
    assert rec == {
        "eventType": "Temp Basal",
        "rate": 1.25,
        "duration": 30,
        "created_at": iso(START + 10 * 60000),
        "enteredBy": "ghostlows-aid",
    }
    assert led.set_temp(20, -0.4, 30) == 0.0  # oref0 never asks for less than zero; clamp anyway


def test_superseded_temp_keeps_its_actual_duration():
    led = Ledger(START)
    led.set_temp(0, 2.0, 30)
    led.set_temp(10, 0.5, 30)
    assert led.temps() == [(0, 2.0, 10), (10, 0.5, 30)]
    assert led.currenttemp(15) == {"duration": 25, "rate": 0.5, "temp": "absolute"}
    assert led.active_temp(39).rate == 0.5
    assert led.active_temp(40) is None
    assert led.currenttemp(40) == NO_TEMP


def test_cancel_records_actual_duration_and_no_new_record():
    led = Ledger(START)
    led.set_temp(5, 3.0, 30)
    led.set_temp(20, 0, 0)  # duration 0: cancel (lib/basal-set-temp.js)
    assert led.temps() == [(5, 3.0, 15)]
    assert led.currenttemp(20) == NO_TEMP
    led.set_temp(25, 0, 0)  # cancelling with nothing running changes nothing
    assert led.temps() == [(5, 3.0, 15)]


def test_expired_temp_is_not_shortened():
    led = Ledger(START)
    led.set_temp(0, 2.0, 30)
    led.set_temp(45, 1.0, 30)
    assert led.temps() == [(0, 2.0, 30), (45, 1.0, 30)]


def test_replacing_a_temp_in_its_first_minute_drops_it():
    led = Ledger(START)
    led.set_temp(10, 2.0, 30)
    led.set_temp(10, 0.0, 30)
    assert led.temps() == [(10, 0.0, 30)]


def test_boluses_as_nightscout_treatments():
    led = Ledger(START)
    led.bolus(0, 4.5, carbs=45.0)
    led.bolus(5, 0.3, smb=True)
    led.bolus(10, 1.0)
    led.bolus(15, 0.0)  # nothing to record
    led.bolus(20, 0.0, carbs=10.0)  # carbs entered with no insulin (wizard floored at 0)
    events = [
        (m, r["eventType"], r.get("insulin"), r.get("carbs"), r.get("isSMB"))
        for m, r in led.records
    ]
    assert events == [
        (0, "Meal Bolus", 4.5, 45.0, None),
        (5, "Correction Bolus", 0.3, None, True),
        (10, "Correction Bolus", 1.0, None, None),
        (20, "Meal Bolus", 0.0, 10.0, None),
    ]
    assert led.boluses() == [(0, 4.5), (5, 0.3), (10, 1.0)]


def test_treatments_newest_first_within_window():
    led = Ledger(START)
    led.set_temp(0, 1.0, 30)
    led.bolus(60, 2.0, carbs=20.0)
    led.set_temp(120, 0.0, 30)
    led.bolus(120, 0.2, smb=True)  # same minute: kept by history.js:274-280 (equal times pass)
    out = led.treatments(180, hours=2)
    assert [r["created_at"] for r in out] == [iso(START + m * 60000) for m in (120, 120, 60)]
    assert [r["eventType"] for r in out] == ["Correction Bolus", "Temp Basal", "Meal Bolus"]
    times = [r["created_at"] for r in led.treatments(180, hours=24)]
    assert times == sorted(times, reverse=True)
    assert led.treatments(100, hours=24)[0]["created_at"] == iso(START + 60 * 60000)  # no future


def test_round_bolus():
    assert round_bolus(1.04, 0.1) == 1.0
    assert round_bolus(1.05, 0.1) == 1.1
    assert round_bolus(0.049, 0.05) == 0.05
    assert round_bolus(1.2345, None) == 1.2345


RATES = [
    0.0,
    0.012,
    0.024,
    0.025,
    0.026,
    0.49,
    0.975,
    0.999,
    1.0,
    1.0249,
    1.025,
    1.2735,
    2.675,
    3.801,
    5.0694,
    9.97,
    9.975,
    10.0,
    10.04,
    10.05,
    12.345,
]


@requires_node
def test_round_basal_matches_oref0():
    script = (
        "const r=require('oref0/lib/round-basal');"
        f"console.log(JSON.stringify({json.dumps(RATES)}.map(x=>r(x))));"
    )
    out = subprocess.run(
        [shutil.which("node"), "-e", script],
        cwd=NODE_DIR,
        capture_output=True,
        text=True,
        check=True,
    )
    assert json.loads(out.stdout) == [round_basal(x) for x in RATES]
