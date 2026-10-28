"""Small hand-made traces for the outcome tests."""

from __future__ import annotations

import numpy as np

from ghostlows_aid.analysis.outcomes import PRE_MIN, ArmTrace

BASAL = 0.02  # U/min


def make_trace(arm: str = "clean", pre: int = PRE_MIN, window: int = 720, **kw) -> ArmTrace:
    """A flat, quiet arm (truth and CGM 120 mg/dL, scheduled basal only) with overrides."""
    n1, n5 = pre + window, (pre + window) // 5
    slot = kw.pop("slot", "night" if arm.startswith(("clean", "pisa")) else "any")
    base = dict(
        controller="oref0_tb",
        patient_id="sim:adult#001",
        seed=101,
        slot=slot,
        arm=arm,
        onset_minute=24 * 60 + 23 * 60,
        truth=np.full(n1, 120.0),
        basal_u=np.full(n1, BASAL),
        bolus_u=np.zeros(n1),
        smb_u=np.zeros(n1),
        rescue_g=np.zeros(n1),
        cgm_display=np.full(n5, 120.0),
        cgm_controller=np.full(n5, 120.0),
        iob_u=np.full(n5, 0.5),
        cob_g=np.full(n5, 0.0),
        reasons=[""] * n5,
        scheduled_basal_u_per_min=BASAL,
        last_bolus_minute=24 * 60 + 19 * 60,
        pre_min=pre,
        window_min=window,
    )
    base.update(kw)
    return ArmTrace(**base)


def window_values(n: int, segments: list[tuple[int, int, float]], fill: float = 120.0):
    """An array of n values, `fill` except value v on [a, b) for each (a, b, v)."""
    x = np.full(n, fill)
    for a, b, v in segments:
        x[a:b] = v
    return x
