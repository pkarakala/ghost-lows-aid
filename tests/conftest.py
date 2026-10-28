from __future__ import annotations

import pytest

from ghostlows_aid.oref import Bridge, glucose_entry, iso, node_available

requires_node = pytest.mark.skipif(not node_available(), reason="node or node_modules/oref0 absent")

CLOCK = 1767614400000  # 2026-01-05T12:00:00Z, a simulated time in paper 1's calendar


@pytest.fixture(scope="session")
def bridge():
    if not node_available():
        pytest.skip("node or node_modules/oref0 absent")
    with Bridge() as b:
        yield b


def synthetic_state(
    clock: int = CLOCK, hours: float = 8.0, gap_minutes=(), meal=(95, 4.5, 50)
) -> dict:
    """A plausible loop state: CGM every 5 min (whole mg/dL, with an optional gap), and a pump
    history with temps of varying rates, a meal bolus with carbs and two SMB-sized boluses.
    Newest first, as oref0 wants both."""
    import math

    glucose = []
    n = int(hours * 12)
    for i in range(n):
        minutes_ago = 5 * i
        if minutes_ago in gap_minutes:
            continue
        value = round(130 + 45 * math.sin(minutes_ago / 70.0) - 0.08 * minutes_ago)
        glucose.append(glucose_entry(clock - minutes_ago * 60000, value))
    treatments = []
    rates = [0.0, 0.4, 1.6, 2.3, 0.9, 0.0, 1.25, 3.1, 0.6, 1.8]
    for j, minutes_ago in enumerate(range(5, int(hours * 60), 15)):
        treatments.append(
            {
                "eventType": "Temp Basal",
                "rate": rates[j % len(rates)],
                "duration": 15 if j else 30,  # the running temp has time left
                "created_at": iso(clock - minutes_ago * 60000),
                "enteredBy": "ghostlows-aid",
            }
        )
    boluses = [(40, 0.3, 0), (20, 0.2, 0)] + ([meal] if meal else [])
    for minutes_ago, insulin, carbs in boluses:
        rec = {
            "eventType": "Meal Bolus" if carbs else "Correction Bolus",
            "insulin": insulin,
            "created_at": iso(clock - minutes_ago * 60000),
            "enteredBy": "ghostlows-aid",
        }
        if carbs:
            rec["carbs"] = carbs
        treatments.append(rec)
    treatments.sort(key=lambda r: r["created_at"], reverse=True)
    currenttemp = {"duration": 25, "rate": treatments_rate(treatments), "temp": "absolute"}
    return {"glucose": glucose, "treatments": treatments, "currenttemp": currenttemp}


def treatments_rate(treatments: list[dict]) -> float:
    return next(t["rate"] for t in treatments if t["eventType"] == "Temp Basal")
