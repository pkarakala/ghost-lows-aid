# Analysis notes (agent layer)

Hand-over record for the analysis code of paper 2 (`src/ghostlows_aid/analysis/`,
`tests/analysis/`, `experiments/analyse.py`). The spec is the study lead's memo
(`ghost-lows-aid-plan/research/PI_decisions.md`, revision 2, §2, §6, §8-§10). Where the memo is
silent, the decision taken is listed in §3 with its reason. No faulted-versus-clean outcome of a
real simulation has been computed: every test and run uses synthetic data or clean-only runs.

## 1. Milestones

| Step | State |
|---|---|
| schema.py (arms, columns, validator), outcomes.py (ArmTrace, per-arm outcomes), alarm_layer.py | done |
| synthetic generator, effects.py (invalid rule, paired effects, bootstrap, Monte Carlo error) | |
| hypotheses.py (H1, H2, three-way rule, Holm) | |
| cross_layer.py (C1 and its sensitivity analyses) | |
| secondary.py | |
| figures.py, experiments/analyse.py | |

## 2. The input contract (ArmTrace) and the outcome table

`outcomes.ArmTrace` is what the harness hands over per arm; `outcomes.trace_from_loop` builds one
from a `LoopResult` of today's loop. `outcomes.slot_outcomes(traces)` turns the arms of one
(controller, patient, seed, slot) into rows of the table defined in `schema.py`.

- 1-min series over [t0 - pre_min, t0 + W): truth, basal delivered, the person's meal boluses,
  correction boluses, SMBs, rescue carbohydrate. 5-min series: the person's displayed CGM, the
  reading sent to the controller, IOB, COB and rT.reason per loop cycle.
- Scalars: scheduled basal (U/min), minute of the person's last bolus before t0, the fault's end
  after t0 (inf when held for the window, NaN for no fault), the realised fault parameters (dict),
  the harness's validity flag, `data_source`.
- pre_min defaults to 120 min (at least 15 for the onset trend). The history before t0 must be
  identical across the arms of a slot; `slot_outcomes` checks it and refuses otherwise.

Fields the harness does not give today: the scheduled basal in LoopResult; a `correction` event
kind for the correcting person; IOB and COB per sample in open loop (onset state); the fault end
and realised parameters per arm (SingleFault.end - onset, SingleFault.params); the arm, slot and
onset of each branch.

## 3. Decisions beyond the memo

Outcomes
- O1 Minutes are counts of 1-min truth samples in [t0, t0 + W); thresholds strictly below 54 and
  70 and strictly above 180 and 250 (paper 1's consensus module). GRI bands: < 54, 54 to < 70,
  > 180 to 250, > 250.
- O2 Hypoglycaemic events use paper 1's `detect_events` on the truth with its history; an event
  counts when its onset lies in the window, so one already running at t0 is not counted.
- O3 "During the fault" is minutes [t0, t0 + ceil(fault end)); after is the rest of the window.
  The clean arm's insulin and minutes above 180 over the same minutes are stored in the arm's row
  (`ref_` columns), so their paired effects compare like with like.
- O4 Time to re-converge is minutes from ceil(fault end) to the start of the first 60-min run
  with |truth - clean truth| < 5 mg/dL; censored (NaN, `reconverged` False) when no such run
  starts inside the window; not applicable (NA) for faults held for the window and arms without a
  fault.
- O5 "Zero temp for a predicted low" is oref0's low-glucose-suspend branch (bg or minGuardBG below
  threshold, determine-basal.js:911-920, reason "minGuardBG x<y"). Zero temps that correct an
  eventual BG below target ("Eventual BG x < y, setting Nm zero temp") are not counted. The
  maxDelta text is written whatever the SMB setting, so it is counted for oref0_tb too.
- O6 Threshold crossings use the reading actually sent to the controller (after a gate) at the
  samples in [t0, fault end); NA without a fault and in open loop.
- O7 Onset state comes from the slot's clean arm (the reading at t0 may already be faulted):
  15-min trends are (value at t0 - value 15 min before) / 15; hours since the person's last meal
  or correction bolus (SMBs excluded); hour of day from the trunk clock.
- O8 Invalid: any truth (history or window) below paper 1's 10 mg/dL floor, or the harness flag.

Alarm layer
- A1 Alarms run on the whole trace with its history (trend fit, Baysal Kalman estimate, refractory
  carry over from the trunk); counted are alarms fired in [t0, t0 + W) and lows with onset there.
  An alarm late in the window before a low that starts after the window counts as false: the same
  edge for every arm.
- A2 Gated alarms (four Baysal sets, matched filter) at their default threshold 0.5, the published
  rule, on open_loop night-slot arms by default (memo §8); trend and threshold on every arm.

## 4. Points in the memo to raise

(filled in as the work goes)
