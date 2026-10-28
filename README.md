# Ghost Lows AID

What CGM sensor faults do to the insulin an open-source closed loop delivers.

[Ghost Lows](https://github.com/pkarakala/ghost-lows) (concept DOI [10.5281/zenodo.23573707](https://doi.org/10.5281/zenodo.23573707)) measured how sensor faults, mainly compression lows, make real-time CGM low alarms fire falsely. A closed loop acts on the same faulted readings: during a false low it withholds insulin, and the insulin it withholds can cause a high later. This study runs the Ghost Lows fault suite inside unmodified open-source automated insulin delivery (OpenAPS oref0 first) in the same virtual cohort. Each fault is injected online, so the controller sees only the faulted reading, and its harm is measured on true glucose against the identical run without the fault.

**Status:** protocol in preparation; no results yet. The protocol is tagged `protocol-v1` before any faulted-versus-clean outcome is computed.

Research software, not a medical device. Do not use it to make treatment decisions.

## Where to start

- People: this README; the protocol in [docs/reference/](docs/reference/) as it lands.
- AI agents: [AGENTS.md](AGENTS.md).

## Authors and disclosure

Independent project by UCSB students; not funded or endorsed by UCSB. Co-authors: Pranav Reddy and Cohen Ferrari. Research question and study direction: Pranav Reddy. The protocol, software and documentation are being developed with AI coding agents (Claude, Anthropic) under Pranav Reddy's direction. Design decisions made by the agents under authority he delegated are logged in docs/reference/DECISIONS.md and await his review. Literature fact-check: Cohen Ferrari (not started). Clinical review: none yet.

## Licence

MIT. Cite this repository with [CITATION.cff](CITATION.cff), and cite Ghost Lows for the fault suite.
