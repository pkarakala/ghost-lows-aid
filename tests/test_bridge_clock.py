"""Wall-clock independence of the bridge (task 1).

The same requests run in two bridges whose whole process wall clock (Date.now and new Date()) is
moved by node/preload_clock.js, once to years before and once to years after the simulated time.
Replies must be identical. A negative control shows that the one wall-clock read the bridge
neutralises (lib/meal/total.js:22-31 -> cob.js:26 -> lib/iob/history.js:175, moment's clock) is
real: left on the wall clock (`wallclock: true`, a test-only switch), lib/meal's output changes.
"""

from __future__ import annotations

from datetime import datetime

import pytest

from ghostlows_aid.oref import NODE_DIR, Bridge
from ghostlows_aid.profile import ProfileSettings, pump_inputs

from .conftest import CLOCK, requires_node, synthetic_state

pytestmark = requires_node

WALL_CLOCKS = ("2019-06-01T03:17:00Z", "2034-11-20T21:41:00Z")


@pytest.fixture(scope="module")
def shifted():
    bridges = [
        Bridge(preload=NODE_DIR / "preload_clock.js", env={"GHOSTLOWS_AID_FAKE_NOW": w})
        for w in WALL_CLOCKS
    ]
    for b in bridges:
        b.start()
    yield bridges
    for b in bridges:
        b.close()


def steps(profile: dict) -> list[dict]:
    out = []
    for smb, gap in ((False, ()), (True, ()), (False, (10, 15, 20))):
        out.append(
            {
                "op": "step",
                "clock": CLOCK,
                "profile": profile,
                "autosens": {"ratio": 1.0},
                "microBolusAllowed": smb,
                "carbs": [],
                "iob_array": True,
                **synthetic_state(gap_minutes=gap),
            }
        )
    return out


def without_timing(reply: dict) -> dict:
    return {k: v for k, v in reply.items() if k not in ("timing_ms", "total_ms", "id")}


def test_preload_moves_the_wall_clock(shifted):
    for b, wall in zip(shifted, WALL_CLOCKS, strict=True):
        reported = b.request({"op": "ping", "clock": CLOCK})["wall_ms"]
        assert reported == int(datetime.fromisoformat(wall).timestamp() * 1000)


def test_same_requests_same_replies_at_two_wall_clocks(shifted):
    inputs = {"op": "profile", "clock": CLOCK, **pump_inputs("adult#001", ProfileSettings())}
    profiles = [b.request(inputs)["profile"] for b in shifted]
    assert profiles[0] == profiles[1]
    for req in steps(profiles[0]):
        a, b = (without_timing(br.request(req)) for br in shifted)
        assert a == b
        assert a["iob"]["iob"] != 0 and a["meal_data"]["carbs"] == 50  # a non-trivial state


def test_negative_control_meal_reads_the_wall_clock_without_the_pin(shifted):
    inputs = {"op": "profile", "clock": CLOCK, **pump_inputs("adult#001", ProfileSettings())}
    profile = shifted[0].request(inputs)["profile"]
    # No carbs in the last 6 h: lib/meal's deviation pass is then cob.js's first call, the one
    # whose IOB history is read without a clock.
    meal = {**steps(profile)[0], **synthetic_state(meal=None), "op": "meal"}
    unpinned = [b.request({**meal, "wallclock": True})["meal_data"] for b in shifted]
    assert unpinned[0] != unpinned[1]  # the read is real
    pinned = [b.request(meal)["meal_data"] for b in shifted]
    assert pinned[0] == pinned[1]  # and the simulated clock removes it
