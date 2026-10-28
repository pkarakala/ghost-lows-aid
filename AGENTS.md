# AGENTS.md

The project rules for any AI agent, whatever the tool; `CLAUDE.md` only points here. People start at the [README](README.md).

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
