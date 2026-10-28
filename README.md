# Ghost Lows AID

What CGM sensor faults do to the insulin an open-source closed loop delivers.

A continuous glucose monitor (CGM) can read falsely low, for example when someone lies on the sensor at night (a compression low). [Ghost Lows](https://github.com/pkarakala/ghost-lows) (concept DOI [10.5281/zenodo.23573707](https://doi.org/10.5281/zenodo.23573707)) measured how such faults set off false low alarms. An automated insulin delivery system acts on the same readings: during a false low it holds back insulin, which can cause a high later, and during a false high it can give too much.

This study asks how much harm each fault type causes when unmodified open-source dosing software, OpenAPS oref0, is in charge. A simulated person (30 virtual patients of the UVA/Padova simulator) lives a day with the software. At a chosen moment the run splits in two: one copy gets a sensor fault, the other does not. Everything else stays the same, so the difference in true glucose is the fault's harm. The study also asks whether the faults that cause the most false alarms are the ones that cause the most insulin harm.

**Status:** protocol written; no results yet. It is frozen by the git tag `protocol-v1` before any faulted run is simulated, so no choice can follow the results.

Research software, not a medical device. Do not use it to make treatment decisions.

## Where to start

- People: [docs/SUMMARY.md](docs/SUMMARY.md), one page.
- Reviewers: [docs/reference/](docs/reference/), the exact protocol and every decision.
- AI agents: [AGENTS.md](AGENTS.md).

## Authors and disclosure

Independent project by UCSB students; not funded or endorsed by UCSB. Co-authors: Pranav Reddy and Cohen Ferrari. Research question and study direction: Pranav Reddy. The protocol, software and documentation are being developed with AI coding agents (Claude, Anthropic) under Pranav Reddy's direction. Design decisions made by the agents under authority he delegated are logged in docs/reference/DECISIONS.md and await his review. Literature fact-check: Cohen Ferrari (not started). Clinical review: none yet.

## Licence

MIT. Cite this repository with [CITATION.cff](CITATION.cff), and cite Ghost Lows for the fault suite.
