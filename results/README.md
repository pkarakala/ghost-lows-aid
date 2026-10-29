# Results (protocol-v1)

Computed from the frozen protocol (git tag `protocol-v1`, commit 69c6c13) with the pre-specified analysis. Simulation results only: 30 virtual patients (UVA/Padova 2008 model), OpenAPS oref0 0.7.1, 10 seeds. Read with the rules and limitations in [docs/reference/BENCHMARK.md](../docs/reference/BENCHMARK.md) §7 and [RESEARCH.md](../docs/reference/RESEARCH.md) §6.

| Test | Result | Reading |
|---|---|---|
| H1: typical compression low at night, minutes above 180 mg/dL in 12 h (oref0 temp basal) | +5.4 min (97.5% CI 2.4 to 9.0); SESOI 72 | negligible |
| H2: large drug-interference false high (+60 mg/dL), minutes below 70 mg/dL in 12 h | +50.1 min (95% CI 41.8 to 58.6); SESOI 15 | meaningful harm |
| C1: agreement of fault rankings, false alarms vs insulin harm | Kendall tau-b 0 (95% CI -0.6 to 0.4) | rankings do not agree |

Both H1 and H2 carry the Monte Carlo flag (within-patient error large relative to between-patient spread), a stated limitation.

Sensitivity: S1 (every patient at its smallest passing ISF scale) H1 +4.9 min, negligible; H2 +57.3 min, meaningful harm; insulin ranking tau-b vs primary 0.8. S2 (25 patients within the glucose targets) changes no reading. Additivity (7-day mixed-fault runs): no departure detected, but too imprecise to establish it, so person-year estimates are conditional.

- `tables/analysis/`: primary analysis (summary.md, hypotheses.csv, c1_*.csv, all secondary tables).
- `tables/analysis_s1/`, `tables/analysis_s2/`, `tables/analysis_additivity/`: sensitivity analyses and the additivity check.
- `figures/`: forest plots, C1 layers, dose-response, pathway, gate curves (PNG and SVG; each drawn from a CSV in tables/analysis).
- `calibration/`, `timing/`: the fault-free settings calibration and timing behind DECISIONS.md D16-D18.

Reproduce: `uv run python experiments/run_grid.py run --n 10 --faults --frozen docs/reference/DECISIONS.md --out results/grid`, then `uv run python experiments/analyse.py --input results/grid/outcomes --protocol-tag protocol-v1 --out results/analysis` (SIMULATION.md §10, BENCHMARK.md §6).
