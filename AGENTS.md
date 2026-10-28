# AGENTS.md

The project rules for any AI agent, whatever the tool; `CLAUDE.md`, `GEMINI.md` and `.github/copilot-instructions.md` only point here. People start at the [README](README.md).

## Two layers

- People: [README.md](README.md) and [docs/SUMMARY.md](docs/SUMMARY.md). Short, plain English, no jargon left unexplained; the README stays under 350 words and the summary under one page.
- Agents and reviewers: this file and the specifications in [docs/reference/](docs/reference/). Complete and exact: every number cited or marked ASSUMPTION, every rule testable.
- A fact lives in one layer's source document; the other layer links to it rather than restating it differently.

## What this is

Paper 2 of Ghost Lows. It runs the Ghost Lows fault suite, injected online, inside unmodified open-source automated insulin delivery (oref0 first) in simglucose, and measures the harm each fault causes on true glucose against the identical run without the fault. Ghost Lows ([repository](https://github.com/pkarakala/ghost-lows), concept DOI 10.5281/zenodo.23573707) is paper 1 and a pinned dependency.

## Map

The layout fills in as the protocol and code land.

- [docs/reference/RESEARCH.md](docs/reference/RESEARCH.md): questions, hypotheses, limitations.
- [docs/reference/SIMULATION.md](docs/reference/SIMULATION.md): cohort, controllers, faults, run design.
- [docs/reference/BENCHMARK.md](docs/reference/BENCHMARK.md): outcomes and statistics.
- [docs/reference/DECISIONS.md](docs/reference/DECISIONS.md): decisions and dated amendments.
- [docs/CLINICAL_REVIEW.md](docs/CLINICAL_REVIEW.md): clinical claims awaiting sign-off.
- [docs/CITATIONS.md](docs/CITATIONS.md) and [docs/references.bib](docs/references.bib): every source key.
- [docs/development/ROADMAP.md](docs/development/ROADMAP.md): milestones.
- [docs/development/SPIKE_2026-10-28.md](docs/development/SPIKE_2026-10-28.md): the harness
  hand-over: what is built, how it was checked, calibration, timing, deviations, open problems.
- [node/bridge.js](node/bridge.js): oref0 0.7.1 as a stateless JSON-lines worker on the simulated
  clock (oref0 itself is installed from node/package-lock.json, never edited).
- [src/ghostlows_aid/](src/ghostlows_aid/): the harness. `oref.py` bridge client and controller
  step; `profile.py` pump settings and oref0 profile; `ledger.py` pump history; `loop.py` minute
  engine, person, routing hooks, snapshots; `faults_online.py` paper 1's faults online and single
  events; `detectors.py` online Baysal and matched filter; `episodes.py` trunk, onsets, arm
  catalogue, `run_unit`; `metrics.py` clean-run summaries for calibration only; `analysis/`
  outcomes, effects and statistics (see
  [docs/development/ANALYSIS_NOTES.md](docs/development/ANALYSIS_NOTES.md)).
- [experiments/](experiments/): `acceptance_adult001.py`; `calibrate.py` settings selection;
  `analyse.py` the analysis of outcome tables;
  `additivity.py` D8's additivity runs (7-day clean and mixed-suite pairs, seed 200);
  `levers.py` history and lib/meal levers; `run_grid.py` the grid (timing mode, resumable units,
  the analysis package's ArmTrace contract). Outputs go to results/<name>/ (gitignored).
- [tests/](tests/): fast by default; `-m slow` runs the multi-day loops.

## Commands

- `cd node && npm ci --ignore-scripts`: oref0 and its pinned dependencies, once per checkout.
- `uv run pytest -q`: the fast tests; `uv run pytest -q -m slow`: the multi-day loops.
- `uv run ruff check . && uv run ruff format --check .`: both must pass before a commit.
- `uv run python experiments/calibrate.py run`, then `... table`: settings calibration on the
  calibration seeds (results/calibration/).
- `uv run python experiments/levers.py`: replays recorded oref0 requests with each §2.1 lever.
- `uv run python experiments/run_grid.py run --isf-rule rule1800 --seeds 90001 --patients ...`,
  then `... cost`: timing mode and the §7 projection. `--faults --frozen <DECISIONS.md>` runs
  the faulted arms and is for after the freeze only.
- `uv run python experiments/additivity.py run --frozen <DECISIONS.md>`: D8's paired 7-day
  runs (results/additivity/); `--clean-only` runs the clean half, allowed before the freeze.
- `uv run python experiments/analyse.py --input results/grid ... --additivity results/additivity`
  adds the D8 table; `--exclude-outside-acceptance` is S2; S1 is `--input results/grid_s1
  --out results/analysis_s1`.

## Hard rules

1. Literature values carry `# SOURCE: <key>` from [CITATIONS.md](docs/CITATIONS.md); never invent one.
2. Values we chose carry `# ASSUMPTION: <why>` and are configurable.
3. Code depending on a clinical claim carries `# CLINICAL-REVIEW: CRn`.
4. Faults never touch truth: the controller sees only the faulted CGM; every outcome is scored on true glucose.
5. Paired runs: a faulted run and its clean run share patient, meals, sensor noise and rescue draws; only the fault differs.
6. Controllers run unmodified, at a pinned release. Harness settings that change their behaviour (profile, device name, timing) are documented and tested.
7. Paper 1 is untouched from here: import `ghostlows` at a pinned version; never change the ghost-lows repository from this one.
8. No leakage: anything fitted is fitted on calibration seeds or held-out patients; CIs across patients, not seeds.
9. No data in git, never patient data.
10. Do not overclaim: effect sizes with CIs; limitations from RESEARCH.md; nothing about commercial systems.
11. Amend before computing: a dated amendment in DECISIONS first; `protocol-v1` never moves.

## Writing

- Authors: Pranav Reddy and Cohen Ferrari, co-authors, University of California, Santa Barbara; an independent project, not funded or endorsed by UCSB. `pkarakala` is a username, not a name.
- Cohen Ferrari's work is a literature fact-check; never describe it as clinical review. The team has no clinician.
- Decisions an agent makes under delegated authority are labelled as such in DECISIONS.md until the project lead reviews them.
- Commits and pull requests carry no AI attribution: no `Co-Authored-By` trailer, no "generated with" footer, no session link. The disclosure in the README covers AI use.
- Short, plain sentences. Cite literature claims; give effect sizes with CIs (rule 10).

## Propose a change

- One topic per pull request against `main`.
- Never rewrite git history or move a tag.
- When a rule and a request conflict, stop and ask the project lead.
