"""Bridge parity with oref0's own determine-basal fixtures (task 2).

node/parity_fixtures.js records every determine_basal call made by oref0's
tests/determine-basal.test.js, adds derived fixtures for SMB, stale, flat and noisy CGM, and calls
the real function directly. Each fixture is replayed through the bridge's determine_basal op, which
must return the identical rT.
"""

from __future__ import annotations

import json
import os
import subprocess

import pytest

from ghostlows_aid.oref import NODE_DIR

from .conftest import requires_node

pytestmark = requires_node


@pytest.fixture(scope="module")
def fixtures():
    out = subprocess.run(
        ["node", "parity_fixtures.js"],
        cwd=NODE_DIR,
        env={**os.environ, "TZ": "UTC"},
        capture_output=True,
        text=True,
        check=True,
        timeout=60,
    )
    return json.loads(out.stdout)


def test_fixture_set_covers_the_branches(fixtures):
    fx = fixtures["fixtures"]
    assert len(fx) >= 10
    reasons = [f["rT"].get("reason", "") for f in fx]
    rts = [f["rT"] for f in fx]
    assert any("minGuardBG" in r for r in reasons), "low suspend"
    assert any(rt.get("rate", 0) > 2.0 and rt.get("duration") == 30 for rt in rts), "high temp"
    assert any(rt.get("units", 0) > 0 for rt in rts), "SMB"
    assert any("BG data is too old" in r for r in reasons), "stale data"
    assert any("CGM data is unchanged" in r for r in reasons), "flat CGM"
    assert sum(1 for f in fx if not f["name"].startswith("derived")) >= 40


def test_bridge_reproduces_direct_determine_basal(bridge, fixtures):
    mismatches = []
    for f in fixtures["fixtures"]:
        reply = bridge.request(
            {"op": "determine_basal", "currentTime": fixtures["currentTime"], **f["args"]}
        )
        if reply["rT"] != f["rT"]:
            mismatches.append(f["name"])
    assert mismatches == []
