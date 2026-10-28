"""Paper 1's alarm layer on each arm's displayed CGM, scored against that arm's truth (memo §8).

Everything is paper 1's own code, imported from ghostlows 0.4.0: the `trend` and `threshold`
alarms with the ungated 55 mg/dL urgent-low floor and the 30-min refractory period
(`ghostlows.alarms`), the gated trend alarm for each Baysal 2014 parameter set and the causal
matched filter at their default threshold 0.5 (the published rule), the true low events
(`ghostlows.bench.events.true_events`, Battelino 2023), and the matching and false-alarm rules
(`ghostlows.bench.scoring.fire_and_score`, paper 1 BENCHMARK.md §1-4).

The alarms run on the whole trace, history included, so the trend fit, the Baysal Kalman
estimate and the refractory period carry over from the trunk exactly as in a continuous run.
Counted are alarms fired in [t0, t0 + W) and true lows whose onset is in that window.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from ghostlows.alarms import GatedAlarm, ThresholdAlarm, TrendAlarm, view
from ghostlows.alarms.learned import BaysalDetector, MatchedFilterDetector
from ghostlows.bench.events import true_events
from ghostlows.bench.scoring import DEFAULT_WINDOWS, RUN_KEYS, Windows, fire_and_score
from ghostlows.sim.engine import DEFAULT_START, SAMPLE_MIN

from .schema import ALARM_FIELDS, ALARM_PREFIXES, GATE_SETTINGS, SECOND_DETECTOR

UNGATED = ("trend", "threshold")
GATED = tuple(p for p in ALARM_PREFIXES if p not in UNGATED)


def make_alarms() -> dict:
    """Paper 1's alarms by column prefix (fresh instances: gates cache per run)."""
    alarms = {"trend": TrendAlarm(), "threshold": ThresholdAlarm()}
    for s in GATE_SETTINGS:
        alarms[f"baysal_{s}"] = GatedAlarm(BaysalDetector(s))
    alarms[SECOND_DETECTOR] = GatedAlarm(MatchedFilterDetector())
    return alarms


def gated_by_default(controller: str, slot: str) -> bool:
    """memo §8: the gated alarms are reported on the open_loop night-slot arms."""
    return controller == "open_loop" and slot == "night"


def _frames(traces) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Paper 1-style samples and 1-min truth, one run per arm (run key `mix` = arm name)."""
    obs, truth = [], []
    for t in traces:
        first = DEFAULT_START + pd.Timedelta(minutes=t.onset_minute - t.pre_min)
        keys = {"patient_id": t.patient_id, "seed": int(t.seed), "mix": t.arm}
        n1 = t.pre_min + t.window_min
        minutes5 = np.arange(0, n1, SAMPLE_MIN)
        obs.append(
            pd.DataFrame(
                {
                    **keys,
                    "time": (first + pd.to_timedelta(minutes5, unit="min")).astype(
                        "datetime64[ns]"
                    ),
                    "cgm_observed": t.cgm_display,
                    "fault_type": "none",
                }
            )
        )
        truth.append(
            pd.DataFrame(
                {
                    **keys,
                    "time": (first + pd.to_timedelta(np.arange(n1), unit="min")).astype(
                        "datetime64[ns]"
                    ),
                    "truth_glucose": t.truth,
                }
            )
        )
    return pd.concat(obs, ignore_index=True), pd.concat(truth, ignore_index=True)


def score_traces(
    traces: Sequence, prefixes: Sequence[str], windows: Windows = DEFAULT_WINDOWS
) -> dict[tuple[str, str], dict]:
    """{(arm, prefix): counts} for the arms of one unit: alarms and false alarms fired in the
    window, and caught and missed true lows of each level with onset in the window."""
    observed, truth = _frames(traces)
    events = true_events(truth)
    faults = pd.DataFrame()
    starts = {t.arm: DEFAULT_START + pd.Timedelta(minutes=t.onset_minute) for t in traces}
    ends = {t.arm: starts[t.arm] + pd.Timedelta(minutes=t.window_min) for t in traces}
    alarms = make_alarms()
    out = {}
    for prefix in prefixes:
        alarm = alarms[prefix]
        frame = view(observed, "cgm_observed")
        _, classified, scored = fire_and_score(
            alarm, frame, alarm.default_threshold, events, faults, observed, windows
        )
        for t in traces:
            a = classified[classified["mix"] == t.arm] if len(classified) else classified
            if len(a):
                at = pd.DatetimeIndex(a["time"])
                a = a[(at >= starts[t.arm]) & (at < ends[t.arm])]
            e = scored[scored["mix"] == t.arm] if len(scored) else scored
            if len(e):
                e = e[pd.DatetimeIndex(e["onset"]) >= starts[t.arm]]
            counts = {
                "alarms": float(len(a)),
                "false_alarms": float((~a["matched"].astype(bool)).sum()) if len(a) else 0.0,
            }
            for level in (1, 2):
                lev = e[e["level"] == level] if len(e) else e
                caught = int(lev["caught"].astype(bool).sum()) if len(lev) else 0
                counts[f"l{level}_caught"] = float(caught)
                counts[f"l{level}_missed"] = float(len(lev) - caught)
            out[(t.arm, prefix)] = counts
    return out


def alarm_columns(traces: Sequence, gated: bool | None = None) -> list[dict]:
    """Alarm-layer columns for each trace (schema `<prefix>_<field>`), in the order given.

    `gated` (default: open_loop night slot only, memo §8) adds the gated trend alarms; their
    columns are NaN otherwise.
    """
    if not traces:
        return []
    t0 = traces[0]
    if gated is None:
        gated = gated_by_default(t0.controller, t0.slot)
    prefixes = list(UNGATED) + (list(GATED) if gated else [])
    counts = score_traces(traces, prefixes)
    rows = []
    for t in traces:
        row = {f"{p}_{f}": np.nan for p in ALARM_PREFIXES for f in ALARM_FIELDS}
        for p in prefixes:
            row.update({f"{p}_{f}": v for f, v in counts[(t.arm, p)].items()})
        rows.append(row)
    return rows


__all__ = ["RUN_KEYS", "alarm_columns", "gated_by_default", "make_alarms", "score_traces"]
