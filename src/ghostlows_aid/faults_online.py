"""Online fault injection: paper 1's fault suite as a pre-drawn schedule applied sample by sample.

`FaultSchedule.draw(seed, days, config)` draws every fault instance exactly as
ghostlows.faults.suite.apply_run does (same streams, same order of draws; glucose is never
looked at); `apply(k, raw_k)` then corrupts reading k using only readings up to k. Applied to a
clean trace it reproduces apply_run / apply_faults (tests/test_faults_online.py).

No paper 1 operator needs a future sample: drift, PISA and dropout are functions of time, spike
noise is pre-drawn per sample, and stale holds the first reading of its segment (suite.py:190-193).
The realisation does depend on the planned run length (`days`: the number of candidate faults,
their start times and the complete nights), so a schedule is drawn for the whole planned run.

`SingleFault` is one event of one type at onset t0, for branching episodes (PI_decisions.md rev 2
§6): its shape is drawn from paper 1's own distributions, set at a quantile of them (a severity
level), or set explicitly. Beyond paper 1's operators it has a held multiplicative gain (drift
held for the window, the +/-20% and +/-40% levels), the G6 long-gap model, and a sustained
additive bias with a rise, a hold and an exponential offset: positive for drug interference (a
false high, which pushes oref0 to over-deliver), negative for a sustained false low.

Faults only ever produce `cgm_observed`; truth is never an input (rule 4).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from ghostlows.faults import pisa as pisa_model
from ghostlows.faults.suite import (
    DROPOUT_RATE_PER_MIN,
    LABEL_PRIORITY,
    FaultConfig,
    _candidates,
    draw_pisas,
    run_seed,
    stream,
)

SAMPLE_MIN = 5
RANK = {kind: i for i, kind in enumerate(LABEL_PRIORITY)}
PAPER1_KINDS = ("pisa", "dropout", "stale", "spike", "drift")
SINGLE_KINDS = PAPER1_KINDS + ("gain", "bias")
SINGLE_STREAM = 11  # SeedSequence branch of single events (paper 1 uses 0-4, the person 7)
NONE = ("none", -1)


def _mark(label: tuple, kind: str, fid: int) -> tuple:
    """suite.py:119-124: a sample keeps the highest-priority kind, and within a kind the first
    instance that marked it."""
    if label[0] == "none" or RANK[kind] < RANK[label[0]]:
        return (kind, fid)
    return label


@dataclass
class FaultSchedule:
    """All fault instances of one run, drawn from (seed, days, config) only."""

    config: FaultConfig
    seed: int
    days: float
    drift: list = field(default_factory=list)  # (start min, end min, size, id)
    pisas: list = field(default_factory=list)  # ghostlows.faults.pisa.Pisa, id = position
    spikes: dict = field(default_factory=dict)  # sample index -> [(mg/dL added, id), ...]
    stale: list = field(default_factory=list)  # (first sample, last sample, id)
    dropout: list = field(default_factory=list)  # (start min, end min, id)
    instances: list = field(default_factory=list)  # as apply_run's, times in minutes
    _held: dict = field(default_factory=dict)
    _last_k: int = -1

    @classmethod
    def draw(cls, seed: int, days: float, config: FaultConfig) -> FaultSchedule:
        """`days` is the planned run length (apply_run takes it from the trace, suite.py:108)."""
        sch = cls(config, int(seed), float(days))
        n = int(round(days * 1440 / SAMPLE_MIN))
        minutes = np.arange(n) * float(SAMPLE_MIN)
        s, on = config.severity, config.enabled
        if "drift" in on and s > 0:  # suite.py:126-144
            rng = stream(seed, "drift")
            life = config.sensor_days * 1440
            for k, start in enumerate(np.arange(0, days * 1440, life)):
                size = rng.uniform(*config.drift_range) * rng.choice([-1.0, 1.0]) * s
                sch.drift.append((start, start + life, size, k))
                sch.instances.append(
                    {
                        "fault_type": "drift",
                        "fault_id": k,
                        "start": start,
                        "end": start + life,
                        "size": size,
                    }
                )
        if "pisa" in on and config.pisa_fraction > 0:  # suite.py:147-161
            for k, p in enumerate(draw_pisas(stream(seed, "pisa"), days, config)):
                sch.pisas.append(p)
                sch.instances.append(
                    {
                        "fault_type": "pisa",
                        "fault_id": k,
                        "start": p.onset,
                        "end": p.onset + p.duration,
                        "size": p.depth,
                        "duration": p.duration,
                        "tau": p.tau,
                    }
                )
        if "spike" in on and s > 0:  # suite.py:163-183
            rng = stream(seed, "spike")
            spikes = _candidates(
                rng, days, config.spike_per_day, lambda r: r.uniform(*config.spike_minutes)
            )
            for k, (a, b, u) in enumerate(spikes):
                if u >= s:
                    continue
                idx = np.flatnonzero((minutes >= a) & (minutes < b))
                add = stream(seed, "spike", k).standard_normal(idx.size) * config.spike_sd * s
                for i, v in zip(idx, add, strict=True):
                    sch.spikes.setdefault(int(i), []).append((float(v), k))
                sch.instances.append(
                    {
                        "fault_type": "spike",
                        "fault_id": k,
                        "start": a,
                        "end": b,
                        "size": config.spike_sd * s,
                    }
                )
        if "stale" in on and s > 0:  # suite.py:185-195
            rng = stream(seed, "stale")
            stale = _candidates(
                rng, days, config.stale_per_day, lambda r: r.uniform(*config.stale_minutes)
            )
            for k, (a, b, u) in enumerate(stale):
                seg = np.flatnonzero((minutes >= a) & (minutes < b))
                if u < s and seg.size > 1:
                    sch.stale.append((int(seg[0]), int(seg[-1]), k))
                    sch.instances.append(
                        {"fault_type": "stale", "fault_id": k, "start": a, "end": b}
                    )
        if "dropout" in on and s > 0:  # suite.py:197-208
            rng = stream(seed, "dropout")
            gaps = _candidates(
                rng, days, config.dropout_per_day, lambda r: r.exponential(1 / DROPOUT_RATE_PER_MIN)
            )
            for k, (a, b, u) in enumerate(gaps):
                if u >= s:
                    continue
                sch.dropout.append((a, b, k))
                sch.instances.append({"fault_type": "dropout", "fault_id": k, "start": a, "end": b})
        return sch

    @classmethod
    def for_run(cls, suite_seed: int, key: tuple, days: float, config: FaultConfig):
        """The schedule apply_faults gives the run with identity `key` = (patient_id, seed, mix)
        (suite.py:212-253)."""
        return cls.draw(run_seed(suite_seed, key), days, config)

    def apply(self, k: int, raw: float) -> tuple[float, str, int]:
        """Reading k (sample index from the run start), given its clean value before the device
        range. Returns (observed value or NaN, fault_type, fault_id). Samples must come in
        increasing order; a run may start late (a branch), except inside a stale segment."""
        if k <= self._last_k:
            raise ValueError(
                f"samples must be applied in increasing order: {k} after {self._last_k}"
            )
        self._last_k = k
        m = float(k * SAMPLE_MIN)
        cfg = self.config
        label = NONE
        bias = 1.0
        for start, end, size, fid in self.drift:
            if start <= m < end:
                bias = 1 + size * (m - start) / (cfg.sensor_days * 1440)
                if abs(bias - 1) >= cfg.drift_label_above:
                    label = _mark(label, "drift", fid)
        obs = raw * bias
        minute = np.array([m])
        for fid, p in enumerate(self.pisas):
            obs = obs - float(p.attenuation(minute)[0])
            if p.active(minute)[0]:
                label = _mark(label, "pisa", fid)
        for add, fid in self.spikes.get(k, ()):
            obs = obs + add
            label = _mark(label, "spike", fid)
        for first, last, fid in self.stale:
            if k == first:
                self._held[fid] = obs
            elif first < k <= last:
                if fid not in self._held:
                    raise ValueError(f"stale segment {fid} began before the first applied sample")
                obs = self._held[fid]
                label = _mark(label, "stale", fid)
        for start, end, fid in self.dropout:
            if start <= m < end:
                obs = float("nan")
                label = _mark(label, "dropout", fid)
        return float(obs), label[0], label[1]


# ---- one event at onset t0 ----
# Levels of PI_decisions.md rev 2 §6 (arm numbers in brackets). Gaps of the G6 long-gap model
# [4], in 5-min samples: 1 + continuations, each continuing with probability 0.855 up to the
# third sample and 0.727 after (median 20 min; 73% longer than oref0's 12-min staleness limit).
# SOURCE: Vettoretti 2023 (Model 4 transition probabilities)
# ASSUMPTION: the transition wiring (PI_decisions §6 #4: "to be checked against the paper")
G6_CONTINUE = (0.855, 0.727)
G6_SWITCH_SAMPLE = 3
# Interference, additive [16, 17]: rise g(u) = (u/tp) exp(1 - u/tp) to the peak at tp, hold at
# the peak until `hold_until_min`, then exponential decay with time constant `tau_min`.
# SOURCE: FDA K193371 (repeated 1000 mg vitamin C on Libre 2: +18.4 and +19.7 mg/dL after doses 2
# and 3); Maahs 2015 (1 g acetaminophen on Dexcom G4: largest mean difference 61 mg/dL at 120 min)
# ASSUMPTION: the hold and the decay time constants (no washout data, PI_decisions §6)
INTERFERENCE = {
    "small": {"size_mg_dl": 20.0, "tp_min": 150.0, "hold_until_min": 630.0, "tau_min": 180.0,
              "doses_min": [0.0, 240.0, 480.0]},
    "large": {"size_mg_dl": 60.0, "tp_min": 120.0, "hold_until_min": 120.0, "tau_min": 240.0,
              "doses_min": [0.0]},
}  # fmt: skip
# ASSUMPTION: an additive bias "lasts" until it decays to 10% of its peak, as a PISA does
# (ghostlows.faults.pisa.END_FRACTION)
BIAS_END_FRACTION = pisa_model.END_FRACTION


def g6_gap_samples(rng: np.random.Generator | None = None, q: float | None = None) -> int:
    """Length of one G6 long gap in 5-min samples: drawn, or its quantile q."""

    def p_next(n):  # probability that a gap of n samples continues to n + 1
        return G6_CONTINUE[0] if n < G6_SWITCH_SAMPLE else G6_CONTINUE[1]

    n = 1
    if q is None:
        while rng.random() < p_next(n):
            n += 1
        return n
    survive = 1.0  # P(length > n)
    while True:
        survive *= p_next(n)
        if 1 - survive >= q:
            return n
        n += 1


@dataclass
class SingleFault:
    """One fault event. `onset` in minutes from the run start; `params` set its shape:

    - pisa: depth (mg/dL), duration (min), tau (paper 1's Facchinetti shape, subtractive);
    - dropout: duration (min), readings missing; stale: duration, the first reading repeated;
    - spike: duration, sd (mg/dL), noise_seed (N(0, sd) added per reading);
    - drift and gain: gain (signed fraction), the CGM multiplied by 1 + gain from onset on;
    - bias: size_mg_dl (signed), tp_min, hold_until_min, tau_min (additive; interference).

    Apply in sample order; readings before the onset pass unchanged."""

    kind: str
    onset: float
    params: dict
    fault_id: int = 0
    label_above: float = 0.05  # drift is labelled from |gain| >= 5%, as paper 1
    _held: float | None = None

    def __post_init__(self):
        if self.kind not in SINGLE_KINDS:
            raise ValueError(f"fault kind {self.kind!r} not in {SINGLE_KINDS}")
        if self.kind == "pisa":
            p = self.params
            tau = p.get("tau", pisa_model.TAU_MIN)
            self._pisa = pisa_model.from_size(self.onset, p["depth"], p["duration"], tau)

    # -- shapes --
    @classmethod
    def draw(
        cls,
        kind: str,
        onset: float,
        seed: int,
        config: FaultConfig | None = None,
        gap_model: str = "paper1",
        sign: float | None = None,
        **overrides,
    ) -> SingleFault:
        """Shape drawn from paper 1's distributions (suite.py, pisa.py) at `config.severity`, or
        for a dropout from the G6 model (`gap_model="g6"`); `sign` fixes a drift's direction;
        anything in `overrides` is set explicitly instead."""
        rng = np.random.default_rng(
            np.random.SeedSequence([seed, SINGLE_STREAM, SINGLE_KINDS.index(kind)])
        )
        return cls._make(
            kind, onset, config or FaultConfig(), rng, None, gap_model, sign, overrides
        )

    @classmethod
    def at_quantile(
        cls,
        kind: str,
        onset: float,
        q: float,
        config: FaultConfig | None = None,
        gap_model: str = "paper1",
        sign: float | None = None,
        **overrides,
    ) -> SingleFault:
        """Shape at quantile q of each of the size distributions: a severity level (q = 0.5 the
        median event, q = 0.95 the 95th percentiles, as PISA severe)."""
        if not 0 < q < 1:
            raise ValueError("q must be in (0, 1)")
        cfg = config or FaultConfig()
        return cls._make(kind, onset, cfg, None, q, gap_model, sign, overrides)

    @classmethod
    def interference(cls, onset: float, level: str, **overrides) -> SingleFault:
        """Drug interference at level "small" or "large" (PI_decisions §6 #16-17)."""
        return cls("bias", float(onset), {**INTERFERENCE[level], **overrides})

    @classmethod
    def _make(cls, kind, onset, cfg, rng, q, gap_model, sign, overrides) -> SingleFault:
        def uniform(lo_hi):
            lo, hi = lo_hi
            return float(rng.uniform(lo, hi)) if q is None else lo + q * (hi - lo)

        s = cfg.severity
        if kind == "pisa":
            if q is None:
                depth, duration = pisa_model.sample_sizes(rng, 1, cfg.pisa_tau)
                depth, duration = float(depth[0]), float(duration[0])
            else:
                depth, duration = pisa_model.DEPTH.quantile(q), pisa_model.DURATION.quantile(q)
            params = {"depth": depth, "duration": duration, "tau": cfg.pisa_tau}
        elif kind == "dropout":
            if gap_model == "g6":
                params = {"duration": float(SAMPLE_MIN * g6_gap_samples(rng, q)), "model": "g6"}
            elif gap_model == "paper1":
                rate = DROPOUT_RATE_PER_MIN
                gap = rng.exponential(1 / rate) if q is None else -np.log(1 - q) / rate
                params = {"duration": float(gap), "model": "paper1"}
            else:
                raise ValueError(f"gap_model {gap_model!r} not in ('paper1', 'g6')")
        elif kind == "stale":
            params = {"duration": uniform(cfg.stale_minutes)}
        elif kind == "spike":
            seed = int(rng.integers(2**31)) if q is None else 0
            params = {"duration": uniform(cfg.spike_minutes), "sd": cfg.spike_sd * s,
                      "noise_seed": seed}  # fmt: skip
        elif kind == "drift":
            # ASSUMPTION (PI_decisions §6 #10, #12): paper 1's drift at a uniformly drawn sensor
            # age, held for the window: end-of-life size U(drift_range) x severity (suite.py:132),
            # reached linearly over the sensor's life; a quantile event is at the end of life.
            if sign is None:
                sign = float(rng.choice([-1.0, 1.0])) if q is None else 1.0
            age = float(rng.uniform(0, cfg.sensor_days)) if q is None else cfg.sensor_days
            gain = sign * uniform(cfg.drift_range) * s * age / cfg.sensor_days
            params = {"gain": gain, "sensor_age_days": age}
        else:  # gain and bias: explicit sizes only (stress levels and biases have no distribution)
            needed = "gain" if kind == "gain" else "size_mg_dl"
            if needed not in overrides:
                raise ValueError(f"{kind!r} has no size distribution: give {needed} explicitly")
            params = {}
        params.update(overrides)
        if kind == "bias":
            params = {**INTERFERENCE["large"], "doses_min": [0.0], **params}
        return cls(kind, float(onset), params, label_above=cfg.drift_label_above)

    @property
    def end(self) -> float:
        """Minute at which the fault stops changing readings (PISA: below 1% of its depth; a bias:
        down to 10% of its peak; inf for a gain held for the window)."""
        p = self.params
        if self.kind == "pisa":
            return self.onset + self._pisa.length
        if self.kind in ("dropout", "stale", "spike"):
            return self.onset + p["duration"]
        if self.kind == "bias":
            return self.onset + p["hold_until_min"] + p["tau_min"] * np.log(1 / BIAS_END_FRACTION)
        return float("inf")

    def describe(self) -> dict:
        """The realised fault, for the per-arm log."""
        return {"fault_kind": self.kind, "onset": self.onset, "end": self.end, **self.params}

    def bias_mg_dl(self, t: float) -> float:
        """The additive bias t minutes after onset (kind "bias")."""
        p = self.params
        size, tp, hold = p["size_mg_dl"], p["tp_min"], p["hold_until_min"]
        if t <= 0:
            return 0.0
        if t < tp:
            u = t / tp
            return size * u * float(np.exp(1 - u))
        if t <= hold:
            return size
        return size * float(np.exp(-(t - hold) / p["tau_min"]))

    # -- the operator --
    def apply(self, k: int, raw: float) -> tuple[float, str, int]:
        m = float(k * SAMPLE_MIN)
        t = m - self.onset
        if t < 0:
            return float(raw), "none", -1
        p, kind, fid = self.params, self.kind, self.fault_id
        if kind == "pisa":
            minute = np.array([m])
            obs = raw - float(self._pisa.attenuation(minute)[0])
            if self._pisa.active(minute)[0]:
                return obs, "pisa", fid
            return obs, "none", -1
        if kind in ("drift", "gain"):
            gain = p["gain"]
            label = (kind, fid) if abs(gain) >= self.label_above else NONE
            return raw * (1 + gain), *label
        if kind == "bias":
            b = self.bias_mg_dl(t)
            active = abs(b) >= BIAS_END_FRACTION * abs(p["size_mg_dl"])
            return raw + b, *((kind, fid) if active else NONE)
        if t >= p["duration"]:
            return float(raw), "none", -1
        if kind == "dropout":
            return float("nan"), "dropout", fid
        if kind == "stale":
            if self._held is None:
                self._held = float(raw)  # the first reading of the segment is still real
                return float(raw), "none", -1
            return self._held, "stale", fid
        # spike: N(0, sd) added per reading (paper 1's form), from the event's own stream
        rng = np.random.default_rng(np.random.SeedSequence([p["noise_seed"], k]))
        return raw + float(rng.standard_normal()) * p["sd"], "spike", fid
