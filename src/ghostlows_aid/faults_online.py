"""Online fault injection: paper 1's fault suite as a pre-drawn schedule applied sample by sample.

`FaultSchedule.draw(seed, days, config)` draws every fault instance exactly as
ghostlows.faults.suite.apply_run does (same streams, same order of draws; glucose is never
looked at); `apply(k, raw_k)` then corrupts reading k using only readings up to k. Applied to a
clean trace it reproduces apply_run / apply_faults (tests/test_faults_online.py).

No paper 1 operator needs a future sample: drift, PISA and dropout are functions of time, spike
noise is pre-drawn per sample, and stale holds the first reading of its segment (suite.py:190-193).
The realisation does depend on the planned run length (`days`: the number of candidate faults,
their start times and the complete nights), so a schedule is drawn for the whole planned run.

`SingleFault` is one event of one type at onset t0, for branching episodes: its shape is drawn
from paper 1's own distributions, set at a quantile of them (a severity level), or set explicitly.
It adds two fault types paper 1 lacks: a sustained positive bias (false high, for example drug
interference), which pushes oref0 to over-deliver, and a sustained negative bias.

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
BIAS_KINDS = ("bias_pos", "bias_neg")
SINGLE_KINDS = PAPER1_KINDS + BIAS_KINDS
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
# ASSUMPTION: placeholder shape and sizes of the sustained biases until the protocol fixes them
# from the literature (drug interference such as acetaminophen reads falsely high for hours):
# additive, linear rise over 30 min, plateau of 6 h, linear return over 60 min; size uniform on
# 20-80 mg/dL, so quantile 0.5 is 50 mg/dL.
BIAS_SHAPE = {"rise_min": 30.0, "plateau_min": 360.0, "fall_min": 60.0}
BIAS_SIZE_MG_DL = (20.0, 80.0)


@dataclass
class SingleFault:
    """One fault event. `onset` in minutes from the run start; `params` set its shape (see
    `draw`). Apply in sample order; samples before the onset pass unchanged."""

    kind: str
    onset: float
    params: dict
    fault_id: int = 0
    label_above: float = 0.05  # drift is labelled from |bias| >= 5%, as paper 1
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
        cls, kind: str, onset: float, seed: int, config: FaultConfig | None = None, **overrides
    ) -> SingleFault:
        """Shape drawn from paper 1's distributions (suite.py, pisa.py) at `config.severity`;
        anything in `overrides` is set explicitly instead."""
        rng = np.random.default_rng(
            np.random.SeedSequence([seed, SINGLE_STREAM, SINGLE_KINDS.index(kind)])
        )
        return cls._make(kind, onset, config or FaultConfig(), rng, None, overrides)

    @classmethod
    def at_quantile(
        cls, kind: str, onset: float, q: float, config: FaultConfig | None = None, **overrides
    ) -> SingleFault:
        """Shape at quantile q of each of paper 1's size distributions: a severity level
        (q = 0.5 the median event, q = 0.95 a large one)."""
        if not 0 < q < 1:
            raise ValueError("q must be in (0, 1)")
        return cls._make(kind, onset, config or FaultConfig(), None, q, overrides)

    @classmethod
    def _make(cls, kind, onset, cfg, rng, q, overrides) -> SingleFault:
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
            rate = DROPOUT_RATE_PER_MIN
            gap = rng.exponential(1 / rate) if q is None else -np.log(1 - q) / rate
            params = {"duration": float(gap)}
        elif kind == "stale":
            params = {"duration": uniform(cfg.stale_minutes)}
        elif kind == "spike":
            seed = int(rng.integers(2**31)) if q is None else 0
            params = {
                "duration": uniform(cfg.spike_minutes),
                "sd": cfg.spike_sd * s,
                "noise_seed": seed,
            }
        elif kind == "drift":
            # ASSUMPTION: paper 1's end-of-life bias (suite.py:132), reached linearly over the
            # sensor's life; the event is the sensor at age `age_days` (uniform when drawn, end of
            # life at a quantile), its bias capped at the end-of-life value within the window.
            sign = float(rng.choice([-1.0, 1.0])) if q is None else 1.0
            age = float(rng.uniform(0, cfg.sensor_days)) if q is None else cfg.sensor_days
            params = {
                "size": sign * uniform(cfg.drift_range) * s,
                "age_days": age,
                "sensor_days": cfg.sensor_days,
            }
        else:
            sign = 1.0 if kind == "bias_pos" else -1.0
            params = {**BIAS_SHAPE, "size_mg_dl": sign * uniform(BIAS_SIZE_MG_DL)}
        params.update(overrides)
        return cls(kind, float(onset), params, label_above=cfg.drift_label_above)

    @property
    def end(self) -> float:
        """First minute at which the fault no longer changes a reading (inf for drift)."""
        p = self.params
        if self.kind == "pisa":
            return self.onset + self._pisa.length
        if self.kind in ("dropout", "stale", "spike"):
            return self.onset + p["duration"]
        if self.kind == "drift":
            return float("inf")
        return self.onset + p["rise_min"] + p["plateau_min"] + p["fall_min"]

    def bias_mg_dl(self, t: float) -> float:
        """The additive bias t minutes after onset (bias kinds)."""
        p = self.params
        rise, plateau, fall = p["rise_min"], p["plateau_min"], p["fall_min"]
        if t < 0 or t >= rise + plateau + fall:
            return 0.0
        if t < rise:
            return p["size_mg_dl"] * t / rise
        if t < rise + plateau:
            return p["size_mg_dl"]
        return p["size_mg_dl"] * (1 - (t - rise - plateau) / fall)

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
        if kind == "drift":
            age = min(p["age_days"] + t / 1440, p["sensor_days"])
            bias = p["size"] * age / p["sensor_days"]
            label = ("drift", fid) if abs(bias) >= self.label_above else NONE
            return raw * (1 + bias), *label
        if kind in BIAS_KINDS:
            b = self.bias_mg_dl(t)
            return (raw + b, kind, fid) if b != 0.0 else (float(raw), "none", -1)
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
