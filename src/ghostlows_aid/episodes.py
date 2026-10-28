"""Branching counterfactual episodes (PI_decisions.md rev 2 §2, §6, §7).

A trunk is one fault-free closed-loop run of a (controller, patient, seed) from 00:00 on day 1 for
67 h; the first 24 h are warm-up. At one onset t0 per slot (night, dose, any) the whole simulated
state is snapshotted (ClosedLoop.snapshot: body and integrator, sensor-noise generator, pump
ledger, every history the controller and the person read, the detector and fault objects, the
person's rescue timer). Each arm restores the snapshot of its slot, sets its fault, routing,
person and gate, and runs for W = 12 h. The fault-free arm of a slot is the trunk's own
continuation over [t0, t0 + W]; a branched arm with no fault equals it exactly
(tests/test_branching.py).

Random numbers, and why the arms stay in lockstep:
- meals, carb-count errors, missed boluses and skipped meals, sensitivity drift and exercise are
  drawn before the run by paper 1's `make_scenario` (meal and perturbation seeds of the seed and
  patient): identical in every arm;
- sensor noise: simglucose's CGMNoise is one sequential stream with one draw per 5-min sample,
  drawn whether or not the reading is faulted, missing, gated or shown; its state is in the
  snapshot, so sample j of every arm gets the same draw;
- the person's night-response draw is keyed by (person seed, 7, sample index): the same number
  for the same sample in every arm;
- onsets: (seed, patient, 12, slot); fault shapes: (seed, patient, 13, fault number), so every
  controller and every variant arm of a fault meets the same realised event; spike noise:
  (event seed, sample index);
- the controller (oref0 through the stateless bridge) and the detectors have no randomness.
No random draw can diverge between arms. What diverges is only what the shared numbers are used
for: a rescue taken or not, a wizard dose computed from a different reading.
"""

from __future__ import annotations

import time
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field

import numpy as np
from ghostlows.sim.engine import PATIENTS
from ghostlows.sim.scenario import DEFAULT, Mix

from .detectors import make_detector
from .faults_online import SingleFault
from .loop import CONTROLLERS, CORRECTING, LoopConfig, LoopResult, Person, Routing, make_loop
from .profile import ProfileSettings

# ---- the design (PI_decisions rev 2 §2, §6) ----
WARMUP_MIN = 24 * 60
WINDOW_MIN = 12 * 60  # W  ASSUMPTION: PI_decisions §2
TRUNK_MIN = 67 * 60  # 24-h warm-up + the latest night onset (07:00 day 3) + W
SCENARIO_DAYS = 3
# minutes from 00:00 on day 1, [start, end)
SLOTS = {"night": (2820, 3300), "dose": (1860, 2340), "any": (1440, 2880)}
SLOT_ORDER = ("night", "dose", "any")
ONSET_STREAM, FAULT_STREAM = 12, 13


def onset(slot: str, seed: int, patient: str, k: int, n: int) -> int:
    """t0 of `slot` for seed k (0-based) of the n seeds of a patient: uniform over the 5-min
    boundaries in the k-th of n equal sub-intervals of the slot (stratified onsets)."""
    lo, hi = SLOTS[slot]
    a, b = lo + (hi - lo) * k / n, lo + (hi - lo) * (k + 1) / n
    grid = np.arange(int(np.ceil(a / 5)) * 5, b, 5)
    if grid.size == 0:
        raise ValueError(f"slot {slot} has no 5-min boundary in sub-interval {k} of {n}")
    stream = [seed, PATIENTS.index(patient), ONSET_STREAM, SLOT_ORDER.index(slot)]
    return int(np.random.default_rng(np.random.SeedSequence(stream)).choice(grid))


def event_seed(seed: int, patient: str, number: int) -> int:
    """The seed of fault number `number`'s event in a (seed, patient), shared by controllers."""
    ss = np.random.SeedSequence([seed, PATIENTS.index(patient), FAULT_STREAM, number])
    return int(ss.generate_state(1)[0])


# ---- faults and arms (PI_decisions rev 2 §6) ----
@dataclass(frozen=True)
class FaultSpec:
    number: int  # the row of PI_decisions §6's fault table
    name: str
    slot: str
    tier: int
    make: Callable[[float, int], SingleFault] = field(repr=False, compare=False)


def _f(kind, **kw):
    return lambda t0, seed: SingleFault.draw(kind, t0, seed, **kw)


def _interference(level):
    return lambda t0, seed: SingleFault.interference(t0, level)


FAULTS = [
    FaultSpec(1, "pisa_typical", "night", 1, _f("pisa")),
    FaultSpec(2, "pisa_severe", "night", 2, _f("pisa", depth=48.0, duration=70.0)),
    FaultSpec(3, "dropout_short", "any", 1, _f("dropout")),
    FaultSpec(4, "dropout_long", "any", 2, _f("dropout", gap_model="g6")),
    FaultSpec(5, "dropout_severe", "any", 2, _f("dropout", duration=120.0)),
    FaultSpec(6, "frozen_typical", "any", 1, _f("stale")),
    FaultSpec(7, "frozen_severe", "any", 2, _f("stale", duration=120.0)),
    FaultSpec(8, "spike_typical", "any", 1, _f("spike")),
    FaultSpec(9, "spike_severe", "any", 2, _f("spike", duration=30.0, sd=20.0)),
    FaultSpec(10, "drift_low", "any", 1, _f("drift", sign=-1.0)),
    FaultSpec(11, "drift_low_severe", "any", 2, _f("drift", gain=-0.2)),
    FaultSpec(12, "drift_high", "any", 1, _f("drift", sign=1.0)),
    FaultSpec(13, "drift_high_severe", "any", 2, _f("drift", gain=0.2)),
    FaultSpec(14, "gain_low", "any", 2, _f("gain", gain=-0.4)),
    FaultSpec(15, "gain_high", "any", 2, _f("gain", gain=0.4)),
    FaultSpec(16, "interference_small", "dose", 2, _interference("small")),
    FaultSpec(17, "interference_large", "dose", 1, _interference("large")),
]  # fmt: skip
FAULT_BY_NAME = {f.name: f for f in FAULTS}


@dataclass(frozen=True)
class Arm:
    name: str
    slot: str
    tier: int
    fault: str | None  # a FaultSpec name; None = no fault (a variant on the clean continuation)
    controllers: tuple[str, ...]
    group: str = "fault"  # fault, pathway, gate, detector2, guard, correcting
    routing: Routing | None = None
    person: Person | None = None
    gate: str | None = None  # detectors.make_detector spec


def _arms() -> list[Arm]:
    arms = [Arm(f.name, f.slot, f.tier, f.name, CONTROLLERS) for f in FAULTS]
    tb = ("oref0_tb",)
    for fault in ("pisa_typical", "interference_large"):
        spec = FAULT_BY_NAME[fault]
        for name, routing in (
            ("controller_only", Routing(person_sees="clean")),
            ("person_only", Routing(controller_sees="clean")),
        ):
            arms.append(Arm(f"{fault}+{name}", spec.slot, 2, fault, tb, "pathway", routing))
    gated = Routing(gate=True, gate_policy="hold")
    for base in (None, "pisa_typical"):
        for params in ("aggressive", "trial", "nominal", "cautious"):
            name = f"{base or 'clean'}+baysal_{params}"
            arms.append(Arm(name, "night", 2, base, tb, "gate", gated, gate=f"baysal_{params}"))
    for base in (None, "pisa_typical"):
        name = f"{base or 'clean'}+matched_filter"
        arms.append(Arm(name, "night", 3, base, tb, "detector2", gated, gate="matched_filter"))
    arms.append(Arm("frozen_severe+fakecgm", "any", 2, "frozen_severe", tb, "guard",
                    Routing(fakecgm=True)))  # fmt: skip
    for base in (None, "interference_large"):
        name = f"{base or 'clean'}+correcting"
        arms.append(Arm(name, "dose", 3, base, tb, "correcting", person=CORRECTING))
    return arms


ARMS = _arms()
ARM_BY_NAME = {a.name: a for a in ARMS}


def clean_arm(slot: str) -> Arm:
    """The fault-free arm of a slot: the trunk's own continuation (no extra simulation)."""
    return Arm(f"clean@{slot}", slot, 1, None, CONTROLLERS, "clean")


def arms_for(controller: str, tiers=(1, 2, 3)) -> list[Arm]:
    return [a for a in ARMS if controller in a.controllers and a.tier in tiers]


# ---- running one unit ----
@dataclass
class ArmRun:
    arm: Arm
    t0: int
    fault: dict | None  # SingleFault.describe(): the realised event
    result: LoopResult  # from minute 0 to t0 + W (the trunk for a clean arm: to its end)
    seconds: float  # wall time of the branch (0 for a clean arm)


@dataclass
class Unit:
    controller: str
    patient: str
    seed: int
    k: int  # index of the seed among the patient's n seeds (stratified onsets)
    n: int
    onsets: dict  # slot -> t0
    trunk: LoopResult
    trunk_seconds: float
    scheduled_basal_u_per_min: float
    profile: dict | None


def run_unit(
    controller: str,
    patient: str,
    seed: int,
    k: int,
    n: int,
    settings: ProfileSettings | None = None,
    arms: list[Arm] | None = None,
    mix: Mix = DEFAULT,
    with_faults: bool = True,
    config_overrides: dict | None = None,
) -> tuple[Unit, Iterator[ArmRun]]:
    """Run the trunk of a (controller, patient, seed) and return it with an iterator over its arms
    (the clean arm of each slot first, then the branched arms, run lazily, one at a time).
    `with_faults=False` runs every arm without its fault (variants kept): for timing only."""
    settings = settings or ProfileSettings()
    config = LoopConfig(
        patient, controller=controller, profile=settings, **(config_overrides or {})
    )
    t_start = time.perf_counter()
    loop = make_loop(config, seed, SCENARIO_DAYS, mix=mix)
    onsets = {slot: onset(slot, seed, patient, k, n) for slot in SLOT_ORDER}
    snapshots = {}
    for slot, t0 in sorted(onsets.items(), key=lambda kv: kv[1]):
        loop.run(t0)
        snapshots[slot] = loop.snapshot()
    loop.run(TRUNK_MIN)
    trunk = loop.result()
    unit = Unit(controller, patient, seed, k, n, onsets, trunk, time.perf_counter() - t_start,
                loop.basal, loop.profile)  # fmt: skip
    chosen = arms_for(controller) if arms is None else arms

    def arm_runs() -> Iterator[ArmRun]:
        for slot in SLOT_ORDER:
            yield ArmRun(clean_arm(slot), onsets[slot], None, trunk, 0.0)
        for arm in chosen:
            t0 = onsets[arm.slot]
            t = time.perf_counter()
            fault = None
            if arm.fault is not None and with_faults:
                spec = FAULT_BY_NAME[arm.fault]
                fault = spec.make(t0, event_seed(seed, patient, spec.number))
            gate = make_detector(arm.gate) if arm.gate else None
            loop.branch(snapshots[arm.slot], fault, arm.routing, arm.person, gate)
            loop.run(t0 + WINDOW_MIN)
            described = fault.describe() if fault is not None else None
            yield ArmRun(arm, t0, described, loop.result(), time.perf_counter() - t)

    return unit, arm_runs()


# ---- what the analysis reads ----
SAMPLE_COLUMNS = (
    "minute", "truth_glucose", "cgm_clean", "cgm_observed", "cgm_person", "cgm_controller",
    "cgm_sent", "restamped", "gate_flag", "fault_type", "rT_reason", "rT_rate", "rT_duration",
    "rT_units", "temp_set", "iob", "cob", "eventualBG", "insulinReq", "insulin_basal",
    "insulin_bolus", "insulin_smb", "carbs_meal", "carbs_rescue",
)  # fmt: skip


def arm_trace(unit: Unit, run: ArmRun, window: int = WINDOW_MIN) -> dict:
    """One arm's traces over [t0, t0 + W]: everything PI_decisions §8 says to log (truth at 1 min;
    displayed CGM and controller input at 5 min; insulin split into basal or temp, boluses and SMB
    at 1 min; events; oref0's reason per cycle; IOB and COB at onset from the trunk; the realised
    fault). Raw traces only: outcomes are the analysis package's (src/ghostlows_aid/analysis/)."""
    r, t0 = run.result, run.t0
    lo, hi = t0, t0 + window
    s = r.samples
    s = s.loc[(s["minute"] >= lo) & (s["minute"] < hi), [c for c in SAMPLE_COLUMNS if c in s]]
    ev = r.events
    onset_row = unit.trunk.samples.loc[unit.trunk.samples["minute"] == t0]
    onset_state = {
        "iob": float(onset_row["iob"].iloc[0]) if "iob" in onset_row else float("nan"),
        "cob": float(onset_row["cob"].iloc[0]) if "cob" in onset_row else float("nan"),
    }
    return {
        "controller": unit.controller,
        "patient": unit.patient,
        "seed": unit.seed,
        "arm": run.arm.name,
        "group": run.arm.group,
        "slot": run.arm.slot,
        "tier": run.arm.tier,
        "t0": t0,
        "window_min": window,
        "fault": run.fault,
        "truth_1min": r.truth[lo:hi],
        "insulin_basal_1min": r.insulin_basal[lo:hi],
        "insulin_bolus_1min": r.insulin_bolus[lo:hi],
        "insulin_smb_1min": r.insulin_smb[lo:hi],
        "scheduled_basal_u_per_min": unit.scheduled_basal_u_per_min,
        "samples": s.reset_index(drop=True),
        "events": ev[(ev["minute"] >= lo) & (ev["minute"] < hi)].reset_index(drop=True),
        "onset_state": onset_state,
        "seconds": run.seconds,
    }


def pre_onset_equal(a: LoopResult, b: LoopResult, t0: int) -> bool:
    """Bit-identical histories before t0: truth, insulin, events and every logged sample value
    (timing columns aside; a column only one run ever filled in counts as a difference)."""
    if not (
        np.array_equal(a.truth[:t0], b.truth[:t0])
        and np.array_equal(a.insulin[:t0], b.insulin[:t0])
    ):
        return False
    ea, eb = (list(r.events[r.events["minute"] < t0].itertuples(index=False)) for r in (a, b))
    if ea != eb:
        return False
    sa, sb = (r.samples[r.samples["minute"] < t0].reset_index(drop=True) for r in (a, b))
    filled = [{c for c in x if not c.startswith("ms_") and x[c].notna().any()} for x in (sa, sb)]
    if filled[0] != filled[1]:
        return False
    return all(sa[c].astype(object).equals(sb[c].astype(object)) for c in filled[0])
