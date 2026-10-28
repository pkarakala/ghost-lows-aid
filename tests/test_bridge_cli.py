"""Bridge parity with the shipped oref0 CLIs (task 2).

bin/oref0-calculate-iob.js, bin/oref0-meal.js and bin/oref0-determine-basal.js run on the same
inputs, written as the JSON files an OpenAPS rig uses, with their wall clock pinned to the simulated
clock by node/preload_clock.js. SMB is off (the CLI disables it without an autotune directory,
bin/oref0-determine-basal.js:105-114), so microBolusAllowed is false in the bridge too. The full
IOB array, meal data, glucose status and rT must be identical.
"""

from __future__ import annotations

import json
import os
import subprocess

import pytest

from ghostlows_aid.oref import NODE_DIR, iso
from ghostlows_aid.profile import ProfileSettings, pump_inputs

from .conftest import CLOCK, requires_node, synthetic_state

pytestmark = requires_node

BIN = NODE_DIR / "node_modules" / "oref0" / "bin"
CASES = {
    "meal 95 min ago": dict(meal=(95, 4.5, 50)),
    "meal 30 min ago": dict(meal=(30, 3.0, 40)),
    "no carbs in 6 h": dict(meal=None),
    "CGM gap": dict(gap_minutes=(5, 10, 15, 20)),
    "under 36 readings": dict(hours=2.5),
}


def cli(args: list[str], cwd, clock_iso: str) -> subprocess.CompletedProcess:
    env = {**os.environ, "TZ": "UTC", "GHOSTLOWS_AID_FAKE_NOW": clock_iso}
    cmd = ["node", "--require", str(NODE_DIR / "preload_clock.js"), *args]
    return subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True, timeout=60)


@pytest.mark.parametrize("isf_rule", ["quest", "rule1800"])
@pytest.mark.parametrize("case", list(CASES))
def test_bridge_equals_shipped_clis(bridge, tmp_path, case, isf_rule):
    inputs = pump_inputs("adult#001", ProfileSettings(isf_rule=isf_rule))
    profile = bridge.request({"op": "profile", "clock": CLOCK, **inputs})["profile"]
    state = synthetic_state(**CASES[case])
    clock_iso = iso(CLOCK)
    files = {
        "pumphistory.json": state["treatments"],
        "profile.json": profile,
        "clock.json": clock_iso,
        "autosens.json": {"ratio": 1.0},
        "glucose.json": state["glucose"],
        "basalprofile.json": profile["basalprofile"],
        "carbhistory.json": [],
        "currenttemp.json": state["currenttemp"],
    }
    for name, data in files.items():
        (tmp_path / name).write_text(json.dumps(data))
    iob = cli(
        [str(BIN / "oref0-calculate-iob.js"), "pumphistory.json", "profile.json", "clock.json"]
        + ["autosens.json"],
        tmp_path,
        clock_iso,
    )
    meal = cli(
        [str(BIN / "oref0-meal.js"), "pumphistory.json", "profile.json", "clock.json"]
        + ["glucose.json", "basalprofile.json", "carbhistory.json"],
        tmp_path,
        clock_iso,
    )
    (tmp_path / "iob.json").write_text(iob.stdout)
    (tmp_path / "meal.json").write_text(meal.stdout)
    db = cli(
        [str(BIN / "oref0-determine-basal.js"), "iob.json", "currenttemp.json", "glucose.json"]
        + ["profile.json", "--auto-sens", "autosens.json", "--meal", "meal.json"]
        + ["--currentTime", clock_iso],
        tmp_path,
        clock_iso,
    )
    assert iob.returncode == meal.returncode == db.returncode == 0
    status_lines = [ln for ln in db.stderr.splitlines() if ln.startswith('{"delta"')]

    reply = bridge.request(
        {
            "op": "step",
            "clock": CLOCK,
            "profile": profile,
            "autosens": {"ratio": 1.0},
            "microBolusAllowed": False,
            "carbs": [],
            "iob_array": True,
            **state,
        }
    )
    assert reply["iob_array"] == json.loads(iob.stdout)
    assert reply["meal_data"] == json.loads(meal.stdout)
    assert reply["rT"] == json.loads(db.stdout)
    assert [reply["glucose_status"]] == [json.loads(ln) for ln in status_lines]
