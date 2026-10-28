"""The closed-loop engine: paper 1's minute-step simglucose physiology, insulin from a controller.

Reused from paper 1 (ghostlows 0.4.0) by import: the T1DPatient set-up and parameter tables, the
scenario (meals, recording errors, sensitivity drift, exercise), the Breton sensor noise on a 5-min
grid, the device range and validity floor, HypoTreatment (15-15 rescue), the seeds of a
patient-run, and the IOB curve. Re-implemented, because insulin now comes from the controller: the
minute loop (it mirrors ghostlows.sim.engine.simulate:150-185 and, in open loop with paper 1's
bolus rule and display, reproduces it bit for bit), the pump, the person's wizard and the routing.

Order within a minute (TECHNICAL_DESIGN §5.2): truth is recorded; on a 5-min sample the CGM is read
(noise, then the fault hook, then the device range), the person may treat a low, the controller runs
on what it is routed; then meals are eaten and bolused; then the pump delivers the minute's insulin.
The controller runs every 5 simulated minutes, also when its reading is missing.
"""

from __future__ import annotations

import bisect
import copy
import math
import types
from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd
from ghostlows.sim.context import iob as python_iob
from ghostlows.sim.engine import (
    DEFAULT_START,
    DEVICE_LIMITS,
    PATIENTS,
    SAMPLE_MIN,
    VALIDITY_FLOOR,
    HypoTreatment,
    Therapy,
    patient_params,
    quest,
    sensor_params,
)
from ghostlows.sim.run import task_seeds
from ghostlows.sim.scenario import STANDARD, Mix, Scenario, make_scenario

from .ledger import Ledger, round_bolus
from .oref import DEVICE, OrefController, OrefError, glucose_entry
from .profile import ProfileSettings, build_profile, scheduled_basal_u_per_min, therapy_numbers

CONTROLLERS = ("oref0_tb", "oref1_smb", "open_loop")
SIGNALS = ("faulted", "clean")
PERSON_STREAM = 7  # SeedSequence branch for the person's own draws (paper 1 uses 0-4 for faults)


# ---- settings ----
@dataclass(frozen=True)
class Person:
    """The virtual person: meal entry and boluses, rescue carbs (paper 1's behaviour otherwise)."""

    # "wizard": announced/CR + max(0, (cgm - target)/ISF) - max(0, IOB), floored at 0, with the
    # profile's CR, ISF and target (TECHNICAL_DESIGN §5.3). "paper1": paper 1's open-loop dose
    # (ghostlows engine.py:178-184: announced/(CR k) + correction above 150 to 140, no IOB).
    bolus_rule: str = "wizard"
    # ASSUMPTION: the wizard subtracts the pump's IOB (oref0's last iob in closed loop; bolus IOB
    # with the same exponential curve in open loop) from the whole bolus, as TECHNICAL_DESIGN §5.3
    # writes it; needs a CLINICAL-REVIEW item (many pumps subtract it from the correction only).
    wizard_iob: bool = True
    # ASSUMPTION: wizard boluses are rounded to the nearest 0.1 U (pump bolus increment; oref0's
    # own default bolus_increment is 0.1, lib/profile/index.js:50). None = no rounding (paper 1).
    bolus_increment: float | None = 0.1
    therapy: Therapy = field(default_factory=Therapy)  # paper 1's settings for "paper1"
    treatment: HypoTreatment | None = field(default_factory=HypoTreatment)  # None: no rescues
    # ASSUMPTION (TECHNICAL_DESIGN §0.2, CR17): rescues are triggered on what the person sees
    # ("observed": their displayed CGM) or, in a sensitivity arm, on true plasma glucose ("truth").
    rescue_signal: str = "observed"
    # ASSUMPTION: probability that an urgent-low alarm at night wakes the person, drawn per 5-min
    # sample from a time-keyed stream. 1.0 reproduces paper 1. Needs a CLINICAL-REVIEW item.
    night_response: float = 1.0
    # ASSUMPTION: the wizard and the rescue rule use the newest displayed reading up to 15 min old;
    # with none, the wizard gives no correction and no rescue is triggered.
    reading_max_age_min: int = 15
    # Correcting-person variant (PI_decisions rev 2 §5, tier 3), off by default: a correction bolus
    # of max(0, (CGM - target)/ISF - IOB) when the displayed CGM is at or above `correct_above` in
    # waking hours and the person gave no bolus (meal or correction; SMBs are the pump's) in the
    # last 3 h. ASSUMPTION: 250 mg/dL, level 2 hyperglycaemia (Battelino 2019); needs a
    # CLINICAL-REVIEW item.
    correct_above: float | None = None
    correct_hours: tuple[int, int] = (7, 23)
    correct_gap_min: int = 180

    def __post_init__(self):
        if self.bolus_rule not in ("wizard", "paper1"):
            raise ValueError(f"bolus_rule {self.bolus_rule!r}")
        if self.rescue_signal not in ("observed", "truth"):
            raise ValueError(f"rescue_signal {self.rescue_signal!r}")


CORRECTING = Person(correct_above=250.0)


@dataclass(frozen=True)
class Routing:
    """Who sees which signal, and the deliberate input manipulations (all off by default)."""

    controller_sees: str = "faulted"  # "clean": the fault reaches only the person
    person_sees: str = "faulted"  # "clean": the fault reaches only the controller
    # Gate a detector (detectors.py) between the displayed CGM and the controller; without `gate`
    # an attached detector only observes (its flags are logged). The person always sees the raw
    # displayed CGM.
    gate: bool = False
    # flagged reading: "hold" re-sends the last unflagged value with a fresh timestamp (oref0's
    # flat-CGM check can still act); "drop" sends nothing (oref0 then sees the data age).
    gate_policy: str = "hold"
    # Guard ablations. These manipulate oref0's INPUTS, never its code (rule 6):
    fakecgm: bool = False  # report device "fakecgm": skips the flat-CGM check (db.js:177-185)
    restamp_missing: bool = False  # re-send the last reading stamped "now" when one is missing:
    # bypasses the stale-data stop (determine-basal.js:187, minAgo > 12)

    def __post_init__(self):
        for v in (self.controller_sees, self.person_sees):
            if v not in SIGNALS:
                raise ValueError(f"routing signal {v!r} not in {SIGNALS}")
        if self.gate_policy not in ("hold", "drop"):
            raise ValueError(f"gate_policy {self.gate_policy!r}")


@dataclass(frozen=True)
class LoopConfig:
    patient: str
    controller: str = "oref0_tb"
    profile: ProfileSettings = field(default_factory=ProfileSettings)
    person: Person = field(default_factory=Person)
    routing: Routing = field(default_factory=Routing)
    start: pd.Timestamp = DEFAULT_START  # read as UTC; the bridge runs with TZ=UTC
    limits: tuple[float, float] = DEVICE_LIMITS
    # ASSUMPTION: the CGM reports whole mg/dL (as Dexcom/Nightscout sgv do); False keeps paper 1's
    # unrounded display. It matters for oref0's flat-CGM check (delta == 0).
    cgm_integer: bool = True
    glucose_hours: float = 24.0  # history sent to oref0 per call
    treatment_hours: float = 24.0
    skip_meal: bool = False  # TECHNICAL_DESIGN §2.1 lever (off; see SPIKE doc)
    fast_params: bool = True  # plain-attribute patient parameters (bitwise equal, 14x faster)

    def __post_init__(self):
        if self.controller not in CONTROLLERS:
            raise ValueError(f"controller {self.controller!r} not in {CONTROLLERS}")
        if self.patient not in PATIENTS:
            raise ValueError(f"unknown patient {self.patient!r}")


# ---- physiology helpers ----
def fast_params(series: pd.Series) -> types.SimpleNamespace:
    """The same parameter values with plain attribute access. T1DPatient.model only reads
    attributes, and pandas Series attribute lookup dominated its cost."""
    return types.SimpleNamespace(**series.to_dict())


def new_patient(name: str):
    from simglucose.patient.t1dpatient import T1DPatient

    base = patient_params(name)
    init = np.asarray(base.iloc[2:15].to_numpy(), dtype=float)  # as ghostlows engine.py:129
    return T1DPatient(base.copy(), init_state=init), base


@dataclass
class PatientState:
    """Everything a T1DPatient carries between steps. The scipy `ode` integrator itself is not
    stored: it is rebuilt from (y, t), which continues bit for bit (dopri5 restarts its step-size
    control on every integrate call; tests/test_branching.py checks the equality)."""

    y: np.ndarray
    t: float
    attrs: dict  # _params, _last_Qsto, _last_foodtaken, _last_action, is_eating, planned_meal, ...

    @classmethod
    def capture(cls, patient) -> PatientState:
        attrs = {k: copy.deepcopy(v) for k, v in patient.__dict__.items() if k != "_odesolver"}
        return cls(np.array(patient.state, dtype=float), float(patient.t), attrs)

    def rebuild(self, patient) -> None:
        from scipy.integrate import ode

        patient.__dict__.update(copy.deepcopy(self.attrs))
        solver = ode(patient.model).set_integrator("dopri5")
        solver.set_initial_value(self.y.copy(), self.t)
        patient._odesolver = solver


def display(value: float, limits: tuple[float, float], integer: bool) -> float:
    """What the device shows: NaN stays missing, values are clipped to the reportable range."""
    if value != value:  # NaN
        return value
    v = min(max(value, limits[0]), limits[1])
    return math.floor(v + 0.5) if integer else v


# ---- state ----
@dataclass
class LoopState:
    minute: int  # next minute to simulate
    patient: PatientState
    params_key: tuple
    noise: object  # simglucose CGMNoise (plain Python state)
    plasma: np.ndarray  # 1-min truth (x3/Vg), filled up to `minute`
    isf: np.ndarray  # 1-min interstitial glucose (x12/Vg)
    insulin: np.ndarray  # 1-min insulin delivered (U): basal_1min + bolus_1min + smb_1min
    basal_1min: np.ndarray  # scheduled basal or the running temp
    bolus_1min: np.ndarray  # the person's meal and correction boluses
    smb_1min: np.ndarray  # oref1 microboluses
    rows: list  # one dict per 5-min sample
    events: list  # (minute, kind, value)
    ledger: Ledger
    hist_minutes: list  # readings sent to the controller (time order)
    hist_entries: list
    person_minutes: list  # readings the person saw (for the wizard and rescues)
    person_values: list
    last_treatment: float
    last_iob: float | None  # the controller's latest IOB (the wizard's "pump IOB")
    faults: object | None
    gate: object | None  # a detector with flag(minute, value) (detectors.py)
    last_unflagged: tuple | None  # (minute, value) of the last reading the gate passed
    controller_errors: int = 0
    last_person_bolus: float = -np.inf  # minute of the person's last bolus (correcting variant)


class ClosedLoop:
    """One patient-run. `run(until)` advances it; `snapshot()`/`restore()` copy its whole state."""

    def __init__(
        self,
        config: LoopConfig,
        scenario: Scenario,
        noise_seed: int,
        person_seed: int,
        faults=None,
        gate=None,
        controller: OrefController | None = None,
        profile: dict | None = None,
    ):
        from simglucose.sensor.noise_gen import CGMNoise

        self.config = config
        self.scenario = scenario
        self.person_seed = int(person_seed)
        self.patient, self.base = new_patient(config.patient)
        self.vg = float(self.base["Vg"])
        self.basal = scheduled_basal_u_per_min(config.patient, config.profile.basal_factor)
        self.start_ms = int(pd.Timestamp(config.start).tz_localize("UTC").value // 10**6)
        self.minutes = scenario.days * 1440
        self.meals = {
            m: [row for _, row in g.iterrows()] for m, g in scenario.meals.groupby("minute")
        }
        self.numbers = therapy_numbers(config.patient, config.profile)
        self.cr_quest, self.cf_quest = quest(config.patient)
        if config.controller == "open_loop":
            self.controller = None
            self.profile = profile
        else:
            mode = "oref1_smb" if config.controller == "oref1_smb" else "oref0_tb"
            settings = replace(config.profile, mode=mode)
            self.profile = profile or build_profile(config.patient, settings)
            self.controller = controller or OrefController(
                self.profile,
                microbolus=mode == "oref1_smb",
                glucose_hours=config.glucose_hours,
                treatment_hours=config.treatment_hours,
                skip_meal=config.skip_meal,
                name=config.controller,
            )
        if config.fast_params:
            self.patient._params = fast_params(self.base)
        n = self.minutes
        self.state = LoopState(
            minute=0,
            patient=None,
            params_key=(1.0, 1.0),
            noise=CGMNoise(sensor_params(), seed=noise_seed),
            plasma=np.full(n, np.nan),
            isf=np.full(n, np.nan),
            insulin=np.zeros(n),
            basal_1min=np.zeros(n),
            bolus_1min=np.zeros(n),
            smb_1min=np.zeros(n),
            rows=[],
            events=[],
            ledger=Ledger(self.start_ms),
            hist_minutes=[],
            hist_entries=[],
            person_minutes=[],
            person_values=[],
            last_treatment=-np.inf,
            last_iob=None,
            faults=faults,
            gate=gate,
            last_unflagged=None,
        )
        self.base_config = config
        self.routing = config.routing
        self.bridge_ms: list[float] = []  # wall time per controller call (diagnostics only)

    # -- snapshot / restore --
    def snapshot(self) -> LoopState:
        s = copy.copy(self.state)
        s.patient = PatientState.capture(self.patient)
        return copy.deepcopy(s)

    def restore(self, state: LoopState) -> None:
        s = copy.deepcopy(state)
        s.patient.rebuild(self.patient)
        s.patient = None
        self.state = s

    def set_faults(self, faults) -> None:
        self.state.faults = faults

    def set_routing(self, routing: Routing) -> None:
        self.routing = routing

    def set_person(self, person: Person) -> None:
        self.config = replace(self.config, person=person)

    def branch(self, snapshot: LoopState, faults=None, routing=None, person=None, gate=None):
        """Restore `snapshot` and set one arm's fault, routing, person and gate (each back to the
        trunk's setting when not given). A gate's detector is first fed every reading the
        controller was shown before the snapshot, so it starts in the state it would have had
        watching from the start; it acts only from here on."""
        self.restore(snapshot)
        self.state.faults = faults
        self.routing = routing or self.base_config.routing
        self.config = replace(self.base_config, person=person or self.base_config.person)
        if gate is not None:
            last = None
            for row in self.state.rows:
                m, v = row["minute"], row["cgm_controller"]
                if not gate.flag(m, v) and v == v:
                    last = (m, v)
            self.state.last_unflagged = last
        self.state.gate = gate
        return self

    # -- the minute loop --
    def run(self, until: int | None = None) -> ClosedLoop:
        until = self.minutes if until is None else min(int(until), self.minutes)
        while self.state.minute < until:
            self._minute(self.state.minute)
            self.state.minute += 1
        return self

    def _set_params(self, m: int) -> None:
        sens, uptake = self.scenario.sensitivity[m], self.scenario.uptake[m]
        if (sens, uptake) != self.state.params_key:
            # as ghostlows engine.py:151-157
            params = self.base.copy()
            params["Vmx"] = float(self.base["Vmx"]) * sens
            params["kp3"] = float(self.base["kp3"]) * sens
            params["Vm0"] = float(self.base["Vm0"]) * uptake
            self.patient._params = fast_params(params) if self.config.fast_params else params
            self.state.params_key = (sens, uptake)

    def _minute(self, m: int) -> None:
        from simglucose.patient.t1dpatient import Action

        s = self.state
        self._set_params(m)
        s.plasma[m] = self.patient.state[3] / self.vg
        s.isf[m] = self.patient.state[12] / self.vg
        cho = 0.0
        row = None
        if m % SAMPLE_MIN == 0:
            row = self._sample(m)
            cho += row["carbs_rescue"]
        else:
            row = s.rows[-1] if s.rows else None
        bolus = (
            row.pop("_correction_now", 0.0) if (row is not None and m % SAMPLE_MIN == 0) else 0.0
        )
        for meal in self.meals.get(m, ()):
            if meal["eaten"]:
                cho += meal["carbs_true"]
                s.events.append((m, "meal_eaten", float(meal["carbs_true"])))
                if row is not None:
                    row["carbs_meal"] += float(meal["carbs_true"])
            if meal["bolus_given"]:
                dose = self._wizard(m, float(meal["carbs_announced"]))
                bolus += dose
                s.ledger.bolus(m, dose, carbs=float(meal["carbs_announced"]))
                s.last_person_bolus = m
                s.events.append((m, "carbs_announced", float(meal["carbs_announced"])))
                s.events.append((m, "bolus", dose))
                if row is not None:
                    row["insulin_bolus"] += dose
        temp = s.ledger.active_temp(m)
        basal = self.basal if temp is None else temp.rate / 60
        smb = row.pop("_smb_now", 0.0) if (row is not None and m % SAMPLE_MIN == 0) else 0.0
        insulin = basal + bolus + smb
        s.insulin[m] = insulin
        s.basal_1min[m], s.bolus_1min[m], s.smb_1min[m] = basal, bolus, smb
        if row is not None:
            row["insulin_basal"] += basal
        self.patient.step(Action(CHO=cho, insulin=insulin))

    # -- the 5-min sample --
    def _sample(self, m: int) -> dict:
        s, cfg, routing = self.state, self.config, self.routing
        k = m // SAMPLE_MIN
        noise = next(s.noise)
        raw = s.isf[m] + noise  # as ghostlows engine.py:163 (before the device range)
        label, fid = "none", -1
        faulted_raw = raw
        if s.faults is not None:
            faulted_raw, label, fid = s.faults.apply(k, raw)
        clean = display(raw, cfg.limits, cfg.cgm_integer)
        faulted = display(faulted_raw, cfg.limits, cfg.cgm_integer)
        to_controller = faulted if routing.controller_sees == "faulted" else clean
        to_person = faulted if routing.person_sees == "faulted" else clean
        row = {
            "k": k,
            "minute": m,
            "truth_glucose": s.plasma[m],
            "isf_glucose": s.isf[m],
            "cgm_noise": noise,
            "cgm_clean": clean,
            "cgm_observed": faulted,
            "fault_type": label,
            "fault_id": fid,
            "cgm_controller": to_controller,
            "cgm_person": to_person,
            "gate_flag": False,
            "insulin_basal": 0.0,
            "insulin_bolus": 0.0,
            "insulin_smb": 0.0,
            "carbs_meal": 0.0,
            "carbs_rescue": 0.0,
        }
        # the person's display and the rescue rule (paper 1's HypoTreatment)
        if to_person == to_person:
            s.person_minutes.append(m)
            s.person_values.append(to_person)
        treat = self.config.person.treatment
        if treat is not None:
            signal = s.plasma[m] if self.config.person.rescue_signal == "truth" else self._seen(m)
            if signal == signal and treat.treats(
                signal, (self._day_offset() + m) % 1440, m - s.last_treatment
            ):
                if self._responds(k, m, treat):
                    row["carbs_rescue"] = treat.grams
                    s.last_treatment = m
                    s.events.append((m, "rescue", treat.grams))
        # the controller
        if self.controller is not None:
            self._control(m, k, to_controller, row)
        if self.config.person.correct_above is not None:
            self._correct(m, row)
        s.rows.append(row)
        return row

    def _correct(self, m: int, row: dict) -> None:
        """The correcting-person variant: a correction bolus on a high displayed reading."""
        person, s = self.config.person, self.state
        cgm = self._seen(m)
        hour = ((self._day_offset() + m) % 1440) // 60
        lo, hi = person.correct_hours
        if not (cgm == cgm and cgm >= person.correct_above and lo <= hour < hi):
            return
        if m - s.last_person_bolus < person.correct_gap_min:
            return
        dose = (cgm - self.config.profile.target) / self.numbers["isf"]
        dose = round_bolus(max(0.0, dose - max(0.0, self._pump_iob(m))), person.bolus_increment)
        if dose <= 0:
            return
        s.ledger.bolus(m, dose)
        s.last_person_bolus = m
        s.events.append((m, "correction", dose))
        row["insulin_bolus"] += dose
        row["_correction_now"] = dose

    def _day_offset(self) -> int:
        start = pd.Timestamp(self.config.start)
        return start.hour * 60 + start.minute

    def _seen(self, m: int) -> float:
        """The person's newest displayed reading, if recent enough."""
        s = self.state
        if s.person_minutes and m - s.person_minutes[-1] <= self.config.person.reading_max_age_min:
            return s.person_values[-1]
        return float("nan")

    def _responds(self, k: int, m: int, treat: HypoTreatment) -> bool:
        p = self.config.person.night_response
        if p >= 1.0:
            return True
        night = treat.threshold((self._day_offset() + m) % 1440, m - self.state.last_treatment)
        if night != treat.night_threshold:
            return True  # awake: the daytime rule applies
        u = np.random.default_rng(np.random.SeedSequence([self.person_seed, PERSON_STREAM, k]))
        return bool(u.random() < p)

    def _control(self, m: int, k: int, value: float, row: dict) -> None:
        s, routing = self.state, self.routing
        clock = self.start_ms + m * 60000
        device = "fakecgm" if routing.fakecgm else DEVICE
        sent = value
        if s.gate is not None:
            flagged = s.gate.flag(m, value)
            row["gate_flag"] = flagged
            if routing.gate and flagged:
                if routing.gate_policy == "drop" or s.last_unflagged is None:
                    sent = float("nan")
                else:
                    sent = s.last_unflagged[1]
            elif value == value:
                s.last_unflagged = (m, value)
        restamped = False
        if sent != sent and routing.restamp_missing and s.hist_entries:
            sent, restamped = s.hist_entries[-1]["glucose"], True
        if sent == sent:
            s.hist_minutes.append(m)
            s.hist_entries.append(glucose_entry(clock, sent, device))
        row["cgm_sent"] = sent
        row["restamped"] = restamped
        lo = bisect.bisect_left(s.hist_minutes, m - self.config.glucose_hours * 60)
        glucose = s.hist_entries[lo:][::-1]
        if not glucose:
            row["rT_reason"] = "no CGM data yet"
            return
        treatments = s.ledger.treatments(m, self.config.treatment_hours)
        currenttemp = s.ledger.currenttemp(m)
        try:
            action, reply = self.controller.step(clock, glucose, treatments, currenttemp)
        except OrefError as e:
            s.controller_errors += 1
            row["rT_reason"] = f"oref0 error: {e}"
            return
        self.bridge_ms.append(reply.get("total_ms", float("nan")))
        self._log_reply(row, reply)
        if action.rate is not None and action.duration is not None:
            row["temp_set"] = s.ledger.set_temp(m, action.rate, action.duration)
        if action.smb > 0:
            s.ledger.bolus(m, action.smb, smb=True)
            s.events.append((m, "smb", action.smb))
            row["insulin_smb"] = action.smb
            row["_smb_now"] = action.smb
        s.last_iob = reply["iob"]["iob"]

    @staticmethod
    def _log_reply(row: dict, reply: dict) -> None:
        gs, iob, meal, rT = reply["glucose_status"], reply["iob"], reply["meal_data"], reply["rT"]
        t = reply.get("timing_ms", {})
        row.update(
            {
                "gs_glucose": gs.get("glucose"),
                "gs_delta": gs.get("delta"),
                "gs_short_avgdelta": gs.get("short_avgdelta"),
                "gs_long_avgdelta": gs.get("long_avgdelta"),
                "gs_date": gs.get("date"),
                "iob": iob.get("iob"),
                "basaliob": iob.get("basaliob"),
                "bolusiob": iob.get("bolusiob"),
                "activity": iob.get("activity"),
                "cob": meal.get("mealCOB"),
                "meal_carbs": meal.get("carbs"),
                "eventualBG": rT.get("eventualBG"),
                "insulinReq": rT.get("insulinReq"),
                "rT_rate": rT.get("rate"),
                "rT_duration": rT.get("duration"),
                "rT_units": rT.get("units"),
                "rT_reason": rT.get("reason"),
                "ms_iob": t.get("iob"),
                "ms_meal": t.get("meal"),
                "ms_determine_basal": t.get("determine_basal"),
                "ms_bridge": reply.get("total_ms"),
            }
        )

    # -- the person's bolus --
    def _wizard(self, m: int, announced: float) -> float:
        person = self.config.person
        if person.bolus_rule == "paper1":
            # ghostlows engine.py:179-182, with the person's displayed reading
            t, cgm = person.therapy, self._seen(m)
            dose = announced / (self.cr_quest * t.bolus_scale)
            if cgm == cgm and cgm > t.correction_above:
                dose += (cgm - t.target) / (self.cf_quest * t.bolus_scale)
            return round_bolus(dose, person.bolus_increment)
        n, cgm = self.numbers, self._seen(m)
        target = self.config.profile.target
        dose = announced / n["carb_ratio"]
        if cgm == cgm:
            dose += max(0.0, (cgm - target) / n["isf"])
        if person.wizard_iob:
            dose -= max(0.0, self._pump_iob(m))
        return round_bolus(max(0.0, dose), person.bolus_increment)

    def _pump_iob(self, m: int) -> float:
        """oref0's latest IOB in closed loop; in open loop the bolus IOB with oref0's curve
        (paper 1's exponential curve, dia 5 h, peak 75 min: the same formula as
        lib/iob/calculate.js:83-143)."""
        if self.controller is not None:
            return self.state.last_iob or 0.0
        boluses = self.state.ledger.boluses()
        if not boluses:
            return 0.0
        t, u = zip(*boluses, strict=True)
        return float(python_iob(t, u, [m], dia_h=self.config.profile.dia_h)[0])

    # -- results --
    def result(self, until: int | None = None) -> LoopResult:
        s = self.state
        n = s.minute if until is None else until
        samples = pd.DataFrame(s.rows)
        if len(samples):
            samples["time"] = pd.Timestamp(self.config.start) + pd.to_timedelta(
                samples["minute"], unit="min"
            )
            samples = samples.drop(columns=[c for c in samples if c.startswith("_")])
        truth = s.plasma[:n]
        valid = bool(np.nanmin(truth) >= VALIDITY_FLOOR) if n else True
        return LoopResult(
            samples=samples,
            truth=truth.copy(),
            insulin=s.insulin[:n].copy(),
            insulin_basal=s.basal_1min[:n].copy(),
            insulin_bolus=s.bolus_1min[:n].copy(),
            insulin_smb=s.smb_1min[:n].copy(),
            events=pd.DataFrame(s.events, columns=["minute", "kind", "value"]),
            treatments=[dict(r) for _, r in s.ledger.records],
            valid=valid,
            nadir=float(np.nanmin(truth)) if n else float("nan"),
            controller_errors=s.controller_errors,
        )


@dataclass
class LoopResult:
    samples: pd.DataFrame  # one row per 5-min sample
    truth: np.ndarray  # 1-min plasma glucose
    insulin: np.ndarray  # 1-min insulin delivered (U)
    insulin_basal: np.ndarray  # of which: scheduled basal or temp
    insulin_bolus: np.ndarray  # the person's boluses
    insulin_smb: np.ndarray  # microboluses
    events: pd.DataFrame
    treatments: list
    valid: bool  # False if plasma fell below VALIDITY_FLOOR (paper 1)
    nadir: float
    controller_errors: int


# ---- building a run ----
def run_seeds(seed: int, patient: str) -> dict:
    """Seeds of one (seed, patient): paper 1's meal, perturbation and noise seeds (so a seed has
    the same meals and sensor noise as in paper 1), and the person's own stream."""
    meal_seed, perturb_seed, noise_seed = task_seeds(seed, patient)
    person = int(
        np.random.SeedSequence([seed, PATIENTS.index(patient), PERSON_STREAM]).generate_state(1)[0]
    )
    return {"meal": meal_seed, "perturb": perturb_seed, "noise": noise_seed, "person": person}


def make_loop(
    config: LoopConfig,
    seed: int,
    days: int,
    mix: Mix = STANDARD,
    faults=None,
    gate=None,
    profile: dict | None = None,
    controller: OrefController | None = None,
) -> ClosedLoop:
    seeds = run_seeds(seed, config.patient)
    scenario = make_scenario(mix, days, seeds["meal"], seeds["perturb"])
    return ClosedLoop(
        config,
        scenario,
        seeds["noise"],
        seeds["person"],
        faults=faults,
        gate=gate,
        profile=profile,
        controller=controller,
    )
