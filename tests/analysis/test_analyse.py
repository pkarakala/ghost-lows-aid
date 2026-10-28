"""experiments/analyse.py: the synthetic run end to end, and the freeze guard."""

from __future__ import annotations

import importlib.util
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from ghostlows_aid.analysis import schema
from ghostlows_aid.analysis.synthetic import Design, synthetic_outcomes

SCRIPT = Path(__file__).resolve().parents[2] / "experiments" / "analyse.py"
spec = importlib.util.spec_from_file_location("analyse", SCRIPT)
analyse = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analyse)


def as_simulation(tmp_path: Path, clean_only: bool = False) -> Path:
    t = synthetic_outcomes(Design(controllers=("oref0_tb",), arms=("pisa_typical",), n_seeds=1))
    if clean_only:
        t = t[t["arm"] == schema.CLEAN]
    t = t.assign(data_source="simulation")
    path = tmp_path / "table.parquet"
    t.to_parquet(path)
    return path


def test_synthetic_run_writes_tables_summary_and_figures(tmp_path):
    out = tmp_path / "out"
    args = ["--synthetic", "--synthetic-seeds", "2", "--out", str(out)]
    assert analyse.main([*args, "--n-boot", "100", "--n-boot-confirmatory", "200"]) == 0
    summary = (out / "summary.md").read_text()
    assert "SYNTHETIC DATA" in summary and "H1" in summary
    h = pd.read_csv(out / "hypotheses.csv")
    assert list(h["id"]) == ["H1", "H2"] and h["decision"].notna().all()
    for name in ("c1_summary", "c1_per_fault", "effects_summary", "gate_optima", "person_year"):
        assert len(pd.read_csv(out / f"{name}.csv")), name
    figs = {p.name for p in (out / "figures").iterdir()}
    for stem in ("forest_risk_weighted_min", "c1_layers_main", "pathway_risk_weighted_min"):
        assert {f"{stem}.png", f"{stem}.svg"} <= figs
    assert {"gate_curves.png", "dose_response_minutes_above_180.svg"} <= figs


def test_faulted_simulation_output_is_refused(tmp_path):
    path = as_simulation(tmp_path)
    out = tmp_path / "out"
    for extra in ([], ["--protocol-tag", "protocol-v2"]):
        with pytest.raises(SystemExit) as stop:
            analyse.main(["--input", str(path), "--out", str(out), *extra])
        assert stop.value.code == analyse.REFUSED
    assert not out.exists()  # nothing computed or written
    with pytest.raises(SystemExit) as stop:
        analyse.main(["--input", str(path), "--synthetic", "--out", str(out)])
    assert stop.value.code == analyse.REFUSED


def test_clean_only_simulation_output_is_descriptive(tmp_path):
    out = tmp_path / "out"
    path = as_simulation(tmp_path, clean_only=True)
    assert analyse.main(["--input", str(path), "--out", str(out), "--n-boot", "50"]) == 0
    assert (out / "clean_levels.csv").exists() and not (out / "hypotheses.csv").exists()


def git(repo: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(repo), *args], check=True, capture_output=True)


def test_protocol_tag_must_exist_and_be_contained(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    git(repo, "init", "-q", "-b", "main")
    git(
        repo,
        "-c",
        "user.name=t",
        "-c",
        "user.email=t@t",
        "commit",
        "-q",
        "--allow-empty",
        "-m",
        "a",
    )
    assert not analyse.protocol_tag_ok("protocol-v1", repo)[0]  # no tag yet
    git(repo, "tag", "protocol-v1")
    assert analyse.protocol_tag_ok("protocol-v1", repo) == (True, "")
    assert not analyse.protocol_tag_ok("protocol-v2", repo)[0]
    assert not analyse.protocol_tag_ok(None, repo)[0]
    git(repo, "checkout", "-q", "-b", "side")
    git(
        repo,
        "-c",
        "user.name=t",
        "-c",
        "user.email=t@t",
        "commit",
        "-q",
        "--allow-empty",
        "-m",
        "b",
    )
    git(repo, "tag", "-f", "protocol-v1")
    git(repo, "checkout", "-q", "main")
    ok, why = analyse.protocol_tag_ok("protocol-v1", repo)
    assert not ok and "does not contain" in why
