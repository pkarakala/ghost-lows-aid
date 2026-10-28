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
| hypotheses.py (H1, H2, four-way rule, Holm) | done |
| cross_layer.py (C1 and its sensitivity analyses) | done |
| secondary.py | done |
| figures.py, experiments/analyse.py | done |
| BENCHMARK §9 requirements (D11 invalid rule, censored after-fault columns, correcting person, per-fault Q7 terms, S1 tau-b, standard-mix note, re-converged share, in-sample gate label, additivity R) | done |

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
  The clean arm's LoopResult is the whole trunk, which runs past t0 + W; its validity is judged
  on the truth up to t0 + W (as the harness's grid runner does), so a later dip does not drop the
  slot.

Alarm layer
- A1 Alarms run on each arm's trace from t0 - 120 min (`pre_min`), not from the trunk's start:
  the trend fit, the Baysal Kalman estimate and the refractory period warm up on those 120 min,
  which are identical across the arms of an event (a limitation copied into every results text,
  BENCHMARK §7); counted are alarms fired in [t0, t0 + W) and lows with onset there.
  An alarm late in the window before a low that starts after the window counts as false: the same
  edge for every arm.
- A2 Gated alarms (four Baysal sets, matched filter) at their default threshold 0.5, the published
  rule, on open_loop night-slot arms by default (memo §8); trend and threshold on every arm.

Effects and bootstrap
- E1 Paired effects: `diff` columns are arm minus the clean row of the same unit; `ref` columns
  (during and after the fault) are the arm minus the clean arm over the same minutes; time to
  re-converge is already paired.
- E2 Invalid rule (D11): a (patient, seed, slot) with any invalid arm, in any controller, is
  dropped for every arm of every controller, so the controllers keep the same events; counts by
  controller, slot and arm, with the rows dropped because another controller triggered.
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

Secondary (each labelled, uncorrected, 95% CIs from the secondary bootstrap)
- S1 Pairing efficiency (N1): "primary outcomes" = those of H1, H2 and C1 (minutes above 180,
  below 70, risk-weighted minutes, trend false alarms), every controller and arm; seeds in
  increasing order for r + 1 mod N; ratio = paired / cross-seed variance, implied independent
  runs per paired run = 1 / median ratio; a patient with both variances 0 is left out.
- S2 Dose-response (N8) for every controller, on minutes above 180 and risk-weighted minutes;
  patients with fewer than two distinct artifact areas are left out of the mean slope. The
  crossing split is for the closed loops (no controller input in open loop).
- S3 Ladders: compression typical-severe; dropout short-long-severe; frozen and spikes
  typical-severe; drift low typical, -20% (drift severe), -40% (gain low), and the same high;
  interference small-large. Each rung against clean and against the rung before, event by event,
  every controller, on risk-weighted minutes and minutes below 70 and above 180.
- S4 State dependence (Q4) on `oref0_tb`'s 17 fault arms and risk-weighted minutes. Bands < 100,
  100 to 140 inclusive, > 140 mg/dL of true glucose at onset; trend of the displayed CGM; IOB
  tertiles within (controller, slot) over the units; night = onset hour in [23, 7). A cell
  with fewer than 5 patients keeps its counts and loses its estimate. Fixed-effects model:
  within-patient demeaning, covariates without within-patient variation left out; in a
  replicate a patient's events count as often as the patient is drawn. Held-out check:
  intercept, the five covariates and age-group dummies; R^2 pooled over the five folds (sum of
  squared errors over sum of squares about each training-fold mean), with the fold range.
- S5 Pathway: the four components on the events that have all three arms (so they add up
  exactly); the wake curve on the events that have both arms.
- S6 Mechanism on the fault arms of every controller; the threshold-crossing shares (patient
  mean of the event indicator) for the closed loops.
- S7 Gate: net effect per night from the two thetas (0.86 and 0.14 from `RATES`); no gate is 0;
  ties go to the earlier option (no gate first). Alarm rule: paper 1's `choose_threshold` over the
  five options in that order, reference = the ungated trend alarm's level 2 sensitivity on the
  same nights, delta = 0 (paper 1's D5 margin); sensitivity pools caught and missed lows of the
  two kinds of night with the same weights; level 1 when the nights hold no level 2 low (paper
  1's A7). The matched filter is reported beside the curves, outside the five-way choice.
  Coincidence share over the replicates in which both optima are defined.
- S8 Guards: the flat-CGM check on `oref0_tb` (the only controller with the variant); the
  stale-data stop for both closed loops.
- S9 Controllers: oref1_smb - oref0_tb for every fault arm on risk-weighted minutes, minutes
  below 54, below 70 and above 180 and the insulin parts; the SMB / temp-basal split for
  oref1_smb's false highs (drift high typical and severe, gain high, interference small and
  large) and spikes.
- S10 Person-year harm per controller on risk-weighted minutes and minutes below 54, below 70 and
  above 180. Central: compression on 0.86 of nights (paper 1's rate, as C1 and the gate use),
  short gaps 2.0, frozen 1.0, spikes 1.0, drift 1 + 1 windows. Low and high: compression on
  0.86 x 0.03 / 0.321 = 0.080 and 0.86 x 0.33 / 0.321 = 0.884 of nights, long gaps (arm 4) 0.13
  and 0.24 in place of short gaps, the rest at paper 1's rates. Each fault's term is reported
  beside the total (D11). The additivity check (D8) runs continuous 7-day pairs
  (experiments/additivity.py) and is analysed by `secondary.additivity` (BENCHMARK §8.5).

Figures and entry point
- F1 Colours: slots 1-3 of the dataviz reference palette, one per controller, validated on all
  pairs against the white figure surface (worst CVD Delta E 9.2, normal vision 24.0); the aqua
  slot's 2.8:1 contrast is relieved by marker shapes, a legend and the CSV behind every figure.
  Single-series figures use slot 1. PNG at 200 dpi and SVG.
- F2 Figures: forest plots of every fault arm by controller for risk-weighted minutes and minutes
  above 180 and below 70; C1's main analysis, both layers with each fault's rank; the PISA typical
  dose-response (oref0_tb, minutes above 180) with the mean per-patient slope drawn through the
  centroid; the pathway 2x2 with its interaction; the gate curves (net effect per night, false
  alarms per night, level 2 sensitivity) with the matched filter beside the four Baysal sets.
- P1 Freeze guard: a row with a faulted arm and data_source other than synthetic needs
  `--protocol-tag protocol-v1`; the tag must resolve to a commit of this repository that HEAD
  contains (`git merge-base --is-ancestor`). Otherwise the run exits with status 2 before
  anything is computed or written. A clean-only simulation table gets the descriptive tables
  only (invalid counts, clean-arm levels). `--synthetic` alone generates the default design;
  with `--input` it requires every row to be synthetic.
- P2 Outputs: synthetic runs to results/synthetic/ (tables, summary.md, manifest.json,
  figures/); simulation runs to results/analysis/ with figures in results/figures/; `--out`
  overrides both (figures in <out>/figures). A full synthetic run (13,860 rows, 10,000 and
  2,000 replicates) takes about 17 s.
- P3 In this worktree run with `PYTHONPATH=src`: the venv's editable install points at the
  harness worktree until the branches are merged.

## 4. Points in the memo to raise

- M1 Resolved in the protocol: Q7, C1 and the gate all use compression on 0.86 of nights; the
  published per-day rates enter only the low and high scenarios, mapped through 0.321 per day
  (BENCHMARK §4.15).
- M2 Resolved in the protocol: the four-way rule adds "opposite direction" (BENCHMARK §4.3),
  as `hypotheses.py` implements.
- M3 C1 sensitivity analysis 3 names a range (0.13 to 0.24 gaps per day), not a rate: run at both
  ends.
- M4 Shares "from patient-level sums" are not comparable across faults when a patient is missing
  for one fault; the code uses means over the patients present (identical when complete).

- M5 N10 weights false alarms 0.86 / 0.14 over the two kinds of night but does not say how the
  level 2 sensitivity combines them; the code pools caught and missed lows with the same
  weights.
- M6 Q7 gives published ranges only for compression and long gaps; using the long-gap arm in
  place of the short gaps in the low and high scenarios is an assumption (as in C1).
