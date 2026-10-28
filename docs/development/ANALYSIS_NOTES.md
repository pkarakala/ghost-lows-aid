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
| synthetic generator, effects.py (invalid rule, paired effects, bootstrap, Monte Carlo error) | done |
| arm names and traces fitted to the harness catalogue (episodes.ARMS, Unit, ArmRun) | done |
| hypotheses.py (H1, H2, three-way rule, Holm) | done |
| cross_layer.py (C1 and its sensitivity analyses) | done |
| secondary.py | |
| figures.py, experiments/analyse.py | |

## 2. The input contract (ArmTrace) and the outcome table

`outcomes.ArmTrace` is what the analysis reads per arm. `outcomes.trace_from_run(unit, run)`
builds one from the harness's `episodes.Unit` and `episodes.ArmRun` (via `trace_from_loop` on the
run's `LoopResult`, which spans minute 0 to t0 + W); `outcomes.slot_outcomes(traces)` turns the
arms of one (controller, patient, seed, slot) into rows of the table defined in `schema.py`.
Checked on clean-only harness units (open_loop and oref0_tb, adult#001, seed 101).

- 1-min series over [t0 - pre_min, t0 + W): truth, basal or temp delivered (`insulin_basal`),
  the person's meal boluses (`insulin_bolus` minus `correction` events), correction boluses,
  SMBs (`insulin_smb`), rescue carbohydrate (`rescue` events). 5-min series: the person's
  displayed CGM (`cgm_person`), the reading sent to the controller after any gate (`cgm_sent`),
  IOB, COB and rT.reason per loop cycle.
- Scalars: scheduled basal (U/min, `Unit.scheduled_basal_u_per_min`), minute of the person's last
  meal or correction bolus before t0, the fault's end after t0 (`describe()["end"] - onset`; inf
  when held, NaN for no fault), the realised fault (`SingleFault.describe()`), the validity flag,
  `data_source`.
- Names: the harness's `clean@<slot>` is the table's `clean` (`schema.table_arm`); patients get
  paper 1's `sim:` prefix; variant arms keep the harness's `<base>+<variant>` names.
- pre_min defaults to 120 min (at least 15 for the onset trend). The history before t0 must be
  identical across the arms of a slot; `slot_outcomes` checks it and refuses otherwise.

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

Effects and bootstrap
- E1 Paired effects: `diff` columns are arm minus the clean row of the same unit; `ref` columns
  (during and after the fault) are the arm minus the clean arm over the same minutes; time to
  re-converge is already paired.
- E2 Invalid rule: a (controller, patient, seed, slot) with any invalid arm is dropped for every
  arm of that slot; counts by controller, slot and arm.
- E3 One joint bootstrap: paper 1's `bootstrap_ci` resampling is recorded once per (patients,
  n_boot, seed) as a replicate x patient count matrix and reused for every quantity, so all
  estimates of an analysis share each replicate's index. In a replicate, patients without a value
  are skipped (mean over the drawn patients that have one). Percentile CIs drop non-finite
  replicates, as bootstrap_ci does; equality with bootstrap_ci is a test.
- E4 Monte Carlo error: within-patient SE = mean over patients (with 2 or more seeds) of
  SD/sqrt(n); between-patient SD = SD of the patient means; flagged when the former exceeds half
  the latter. Also mc_se = sqrt(mean SE^2 / patients).
- E5 Bootstrap seed 0 everywhere; n_boot 10,000 for H1, H2 and C1 (one index shared by all
  three), 2,000 for the rest (another shared index).

Hypotheses
- H-a Holm rejects at p <= the step's level; equal p are ordered by id; a NaN p goes last and is
  never rejected. The CI level follows the step (97.5% first, 95% second), rejected or not.
- H-b TOST uses the 90% percentile CI of the same 10,000 replicates, unadjusted (memo §9).
- H-c Co-reported outcomes use the same replicates, 95% CI, labelled uncorrected.

C1
- C-a Per-day harm of a fault = sum over its arms of rate x theta (drift: 1 x theta_low + 1 x
  theta_high), per replicate; shares and tau-b are computed per replicate from these. With every
  patient present for every fault this equals the ratio of patient-level sums (memo §10); with a
  patient missing for one fault (all seeds invalid) it stays a mean over the patients present.
- C-b Paper 1's s = 0.5 rates: FaultConfig's full-severity rates thinned by s: dropout 2.0,
  frozen 1.0, spikes 1.0 per day (rates paper 1 itself marks as assumptions).
- C-c The long-gap sensitivity analysis runs at both ends of the G6 range (0.13 and 0.24 per day)
  with arm 4 (dropout_long) in place of arm 3.
- C-d Estimability counts the replicates with a positive total among all B (a NaN total is not
  positive). Discordance is NA unless both shares are estimable.
- C-e Kendall's tau-b over the five per-day harms (NaN if one layer is all ties); harms equal to
  1e-9 rank as ties.

## 4. Points in the memo to raise

- M1 Two compression rates are called paper 1's: C1 and the gate use 0.86 per night [Reddy 2026],
  the person-year harm (Q7) uses 0.321 per day [Lobo 2026] as paper 1's central rate. The code
  follows each section as written; the memo should say why they differ or pick one.
- M2 The three-way rule has no category for a meaningful benefit (a significant estimate at or
  below -SESOI): it reads "inconclusive". Implausible for H1 and H2, but the label would mislead.
- M3 C1 sensitivity analysis 3 names a range (0.13 to 0.24 gaps per day), not a rate: run at both
  ends.
- M4 Shares "from patient-level sums" are not comparable across faults when a patient is missing
  for one fault; the code uses means over the patients present (identical when complete).
