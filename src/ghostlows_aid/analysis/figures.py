"""Figures from the analysis tables: the forest plot of per-event effects, C1's two layers, the
PISA dose-response, the pathway 2x2 and the gate curves.

Each title says what is plotted, never what it shows. Colours: slots 1-3 of a categorical
palette validated for colour-vision deficiency on all pairs (dataviz reference palette; worst CVD
Delta E 9.2), one per controller, always with a marker shape and a legend; single-series figures
use slot 1. Every figure is written as PNG and SVG; the numbers behind it are in the CSV tables.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .schema import ARMS, CLEAN, COLUMNS, CONTROLLERS, PRIMARY_CONTROLLER  # noqa: E402

COLOURS = {"oref0_tb": "#2a78d6", "oref1_smb": "#eb6834", "open_loop": "#1baf7a"}
MARKERS = {"oref0_tb": "o", "oref1_smb": "s", "open_loop": "D"}
ONE = COLOURS[PRIMARY_CONTROLLER]
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE = "#e1e0d9", "#c3c2b7", "#ffffff"
FAULT_ARMS = tuple(a for a, x in ARMS.items() if a != CLEAN and not x.variant)
LABELS = {
    "risk_weighted_min": "risk-weighted minutes",
    "minutes_above_180": "minutes above 180 mg/dL",
    "minutes_below_70": "minutes below 70 mg/dL",
    "minutes_below_54": "minutes below 54 mg/dL",
    "trend_false_alarms": "false alarms (trend alarm)",
    "threshold_false_alarms": "false alarms (threshold alarm)",
}

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelcolor": INK2,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK2,
        "ytick.labelcolor": INK2,
        "axes.edgecolor": AXIS,
        "axes.titlecolor": INK,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "svg.fonttype": "none",
    }
)


def _label(outcome: str) -> str:
    return LABELS.get(outcome, outcome.replace("_", " "))


def _unit(outcome: str) -> str:
    col = COLUMNS.get(outcome)
    return f" ({col.unit})" if col is not None and col.unit else ""


def _axes(ax, grid: str = "x") -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.8)
    ax.grid(axis=grid, color=GRID, linewidth=0.6, linestyle="-")
    ax.set_axisbelow(True)
    ax.tick_params(length=3, width=0.6)


def _zero(ax, vertical: bool = True) -> None:
    (ax.axvline if vertical else ax.axhline)(0, color=AXIS, linewidth=1.0, zorder=1)


def _dots(ax, y, x, lo, hi, colour=ONE, marker="o", label=None, horizontal=True) -> None:
    x, lo, hi = (np.asarray(v, float) for v in (x, lo, hi))
    err = np.vstack([x - lo, hi - x])
    kw = dict(
        fmt=marker,
        color=colour,
        ecolor=colour,
        elinewidth=1.4,
        capsize=0,
        markersize=5.5,
        markeredgecolor=SURFACE,
        markeredgewidth=1.0,
        label=label,
        zorder=3,
    )
    if horizontal:
        ax.errorbar(x, y, xerr=err, **kw)
    else:
        ax.errorbar(y, x, yerr=err, **kw)


def _save(fig, out_dir: Path, name: str) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = [out_dir / f"{name}.png", out_dir / f"{name}.svg"]
    fig.savefig(paths[0], dpi=200, bbox_inches="tight")
    fig.savefig(paths[1], bbox_inches="tight")
    plt.close(fig)
    return paths


def forest(summary: pd.DataFrame, outcome: str, out_dir: Path, arms=FAULT_ARMS) -> list[Path]:
    """theta and 95% CI per fault arm (rows) and controller (series)."""
    s = summary[(summary["outcome"] == outcome) & summary["arm"].isin(arms)]
    if s.empty:
        return []
    rows = [a for a in arms if a in set(s["arm"])]
    present = [c for c in CONTROLLERS if c in set(s["controller"])]
    fig, ax = plt.subplots(figsize=(6.4, 0.32 * len(rows) + 1.3))
    step = 0.22
    for k, c in enumerate(present):
        sc = s[s["controller"] == c].set_index("arm").reindex(rows)
        y = np.arange(len(rows)) + (k - (len(present) - 1) / 2) * step
        _dots(ax, y, sc["theta"], sc["ci_lo"], sc["ci_hi"], COLOURS[c], MARKERS[c], c)
    _zero(ax)
    ax.set_yticks(np.arange(len(rows)), [a.replace("_", " ") for a in rows])
    ax.invert_yaxis()
    ax.set_xlabel(f"paired effect per event, {_label(outcome)}{_unit(outcome)}")
    ax.set_title(f"Per-event effect of each fault on {_label(outcome)}, 95% CI", loc="left")
    _axes(ax)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.55, 0.0),
        ncol=len(present),
        frameon=False,
        fontsize=8,
        labelcolor=INK2,
    )
    return _save(fig, out_dir, f"forest_{outcome}")


C1_FAULTS = ("compression", "dropout", "frozen", "spikes", "drift")


def c1_layers(per_fault: pd.DataFrame, out_dir: Path, analysis: str = "main") -> list[Path]:
    """Per-day harm with 95% CI of the five fault types in each layer, with each layer's rank."""
    p = per_fault[per_fault["analysis"] == analysis]
    layers = [x for x in ("alarm", "insulin") if x in set(p["layer"])]
    if not layers:
        return []
    fig, axes = plt.subplots(1, len(layers), figsize=(7.2, 2.8), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, layer in zip(axes, layers, strict=True):
        q = p[p["layer"] == layer].set_index("fault").reindex(list(C1_FAULTS))
        y = np.arange(len(C1_FAULTS))
        _dots(ax, y, q["per_day"], q["ci_lo"], q["ci_hi"])
        for yi, (v, r) in enumerate(zip(q["ci_hi"], q["rank"], strict=True)):
            if np.isfinite(v) and np.isfinite(r):
                ax.annotate(
                    f"rank {int(r)}",
                    (v, yi),
                    xytext=(6, 0),
                    textcoords="offset points",
                    va="center",
                    fontsize=7.5,
                    color=INK2,
                )
        _zero(ax)
        outcome = q["outcome"].dropna().iloc[0] if q["outcome"].notna().any() else ""
        controller = q["controller"].dropna().iloc[0] if q["controller"].notna().any() else ""
        ax.set_xlabel(f"excess {_label(outcome)} per day")
        ax.set_title(f"{layer} layer ({controller})", loc="left")
        ax.margins(x=0.25)
        _axes(ax)
    axes[0].set_yticks(np.arange(len(C1_FAULTS)), list(C1_FAULTS))
    axes[0].invert_yaxis()
    fig.suptitle(
        f"C1: per-day harm of five fault types at their event rates ({analysis})",
        x=0.02,
        ha="left",
        fontsize=10,
        color=INK,
    )
    fig.tight_layout()
    return _save(fig, out_dir, f"c1_layers_{analysis}")


def dose_response(
    effects: pd.DataFrame,
    slopes: pd.DataFrame,
    out_dir: Path,
    outcome: str = "minutes_above_180",
    controller: str = PRIMARY_CONTROLLER,
    arm: str = "pisa_typical",
) -> list[Path]:
    """Event effects of PISA typical against the realised artifact area, with the mean of the
    per-patient least-squares slopes drawn through the centroid."""
    e = effects[(effects["controller"] == controller) & (effects["arm"] == arm)]
    e = e[["fault_area_mg_dl_h", outcome]].astype(float).dropna()
    if e.empty:
        return []
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.scatter(
        e["fault_area_mg_dl_h"],
        e[outcome],
        s=14,
        color=ONE,
        alpha=0.45,
        linewidths=0,
        label="events",
        zorder=2,
    )
    sl = slopes[(slopes["controller"] == controller) & (slopes["outcome"] == outcome)]
    if len(sl) and np.isfinite(sl["slope_per_mg_dl_h"].iloc[0]):
        b = float(sl["slope_per_mg_dl_h"].iloc[0])
        x0, y0 = e["fault_area_mg_dl_h"].mean(), e[outcome].mean()
        xs = np.linspace(e["fault_area_mg_dl_h"].min(), e["fault_area_mg_dl_h"].max(), 50)
        ax.plot(xs, y0 + b * (xs - x0), color=INK, linewidth=1.5, label="mean per-patient slope")
    _zero(ax, vertical=False)
    ax.set_xlabel("artifact area, depth x duration (mg/dL h)")
    ax.set_ylabel(f"paired effect, {_label(outcome)}")
    ax.set_title(f"PISA typical: event effect against artifact area ({controller})", loc="left")
    _axes(ax, grid="both")
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2)
    return _save(fig, out_dir, f"dose_response_{outcome}")


PATHWAY_ROWS = (
    ("reference", "fault-free"),
    ("controller_only", "controller sees the fault"),
    ("person_only", "person sees the fault"),
    ("total", "both see the fault"),
    ("interaction", "interaction"),
)


def pathway(table: pd.DataFrame, out_dir: Path, outcome: str = "risk_weighted_min") -> list[Path]:
    """The 2x2 pathway cells (paired effects against the fault-free arm) and the interaction."""
    t = table[table["outcome"] == outcome]
    bases = [b for b in ("pisa_typical", "interference_large") if b in set(t["arm"])]
    if not bases:
        return []
    fig, axes = plt.subplots(1, len(bases), figsize=(7.2, 2.6), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, base in zip(axes, bases, strict=True):
        q = t[t["arm"] == base].set_index("component")
        y = np.arange(len(PATHWAY_ROWS))
        vals = []
        for k, _ in PATHWAY_ROWS:
            if k == "reference":
                vals.append((0.0, 0.0, 0.0))
            elif k in q.index:
                vals.append(tuple(q.loc[k, ["theta", "ci_lo", "ci_hi"]]))
            else:
                vals.append((np.nan, np.nan, np.nan))
        v = np.array(vals, float)
        _dots(ax, y[1:], v[1:, 0], v[1:, 1], v[1:, 2])
        ax.plot([0], [0], "o", markersize=5.5, markerfacecolor=SURFACE, markeredgecolor=ONE)
        ax.axhline(3.5, color=GRID, linewidth=0.6)
        _zero(ax)
        ax.set_title(base.replace("_", " "), loc="left")
        ax.set_xlabel(f"paired effect, {_label(outcome)}")
        ax.margins(x=0.15)
        _axes(ax)
    axes[0].set_yticks(np.arange(len(PATHWAY_ROWS)), [label for _, label in PATHWAY_ROWS])
    axes[0].invert_yaxis()
    fig.suptitle(
        f"Pathway 2x2 on {_label(outcome)} ({PRIMARY_CONTROLLER}), 95% CI",
        x=0.02,
        ha="left",
        fontsize=10,
        color=INK,
    )
    fig.tight_layout()
    return _save(fig, out_dir, f"pathway_{outcome}")


GATE_PANELS = (
    ("net_per_night", "net effect per night\n(risk-weighted minutes)"),
    ("false_alarms_per_night", "false alarms per night\n(open loop)"),
    ("l2_sensitivity", "level 2 sensitivity\n(open loop)"),
)


def gate_curves(curves: pd.DataFrame, out_dir: Path) -> list[Path]:
    """Net harm per night, false alarms per night and level 2 sensitivity of each gate setting
    (no gate first), with the second detector beside them; 95% CI."""
    if curves.empty:
        return []
    c = curves.set_index("option")
    order = [o for o in ("no_gate", "aggressive", "trial", "nominal", "cautious") if o in c.index]
    extra = [o for o in c.index if o not in order]
    fig, axes = plt.subplots(1, len(GATE_PANELS), figsize=(8.4, 2.9))
    for ax, (col, label) in zip(axes, GATE_PANELS, strict=True):
        x = np.arange(len(order))
        q = c.loc[order]
        ax.plot(x, q[col], color=ONE, linewidth=1.0, zorder=2)
        _dots(ax, x, q[col], q[f"{col}_lo"], q[f"{col}_hi"], horizontal=False)
        xe = len(order) + 0.6 + np.arange(len(extra))
        if extra:
            r = c.loc[extra]
            _dots(ax, xe, r[col], r[f"{col}_lo"], r[f"{col}_hi"], marker="D", horizontal=False)
        if col == "net_per_night":
            _zero(ax, vertical=False)
        ax.set_xticks(
            [*x, *xe], [o.replace("_", " ") for o in (*order, *extra)], rotation=35, ha="right"
        )
        ax.set_title(label, loc="left", fontsize=9)
        _axes(ax, grid="y")
    fig.suptitle(
        "Detector gate settings (Baysal sets; matched filter beside), 95% CI",
        x=0.02,
        ha="left",
        fontsize=10,
        color=INK,
    )
    fig.tight_layout()
    return _save(fig, out_dir, "gate_curves")


def all_figures(
    results: dict[str, pd.DataFrame],
    effects: pd.DataFrame,
    out_dir: Path,
    forest_outcomes: Sequence[str] = ("risk_weighted_min", "minutes_above_180", "minutes_below_70"),
) -> list[Path]:
    """Every figure whose table is in `results` (keys as written by experiments/analyse.py)."""
    out_dir = Path(out_dir)
    paths: list[Path] = []
    if "effects_summary" in results:
        for o in forest_outcomes:
            paths += forest(results["effects_summary"], o, out_dir)
    if "c1_per_fault" in results:
        paths += c1_layers(results["c1_per_fault"], out_dir)
    if "dose_response" in results:
        paths += dose_response(effects, results["dose_response"], out_dir)
    if "pathway" in results:
        paths += pathway(results["pathway"], out_dir)
    if "gate_curves" in results:
        paths += gate_curves(results["gate_curves"], out_dir)
    return paths


__all__ = [
    "COLOURS",
    "all_figures",
    "c1_layers",
    "dose_response",
    "forest",
    "gate_curves",
    "pathway",
]
