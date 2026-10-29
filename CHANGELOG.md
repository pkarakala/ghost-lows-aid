# Changelog

## v0.1.0 (2026-10-29)

First archived release: the frozen protocol and its results.

- Protocol frozen by the git tag `protocol-v1` (commit 69c6c13, 2026-10-28) before any faulted arm was simulated: docs/reference/ (RESEARCH, SIMULATION, BENCHMARK, DECISIONS D1-D19), the closed-loop harness running unmodified OpenAPS oref0 0.7.1, and the analysis code.
- Results of the pre-specified analysis (results/README.md): H1, typical compression low at night, +5.4 min above 180 mg/dL per event, negligible; H2, large drug-interference false high, +50.1 min below 70 mg/dL per event, meaningful harm; C1, false-alarm and insulin-harm rankings of faults do not agree. Sensitivity analyses S1, S2 and the additivity check.
- Simulation results only; no clinician has reviewed the protocol.
