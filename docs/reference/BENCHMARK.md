# Benchmark specification

The analysis of paper 2, frozen by the git tag `protocol-v1`: the unit of analysis, the estimands, every outcome, the paired effects, the statistics, the alarm layer, the sensitivity analyses, the output files and the rules for reading results. It is written for whoever implements, runs, reviews or replicates the analysis. Every formula here is implemented in `src/ghostlows_aid/analysis/` (entry point `experiments/analyse.py`; hand-over record `docs/development/ANALYSIS_NOTES.md`), and each quantity names the function that computes it. Where the code does not yet do what this file says, the section says so under **Requirement on the code**; this file is the specification and the code must be brought to it before the tag.

Why each rule was chosen is in [DECISIONS.md](DECISIONS.md) (D1-D19); the questions, hypotheses and claims are in [RESEARCH.md](RESEARCH.md); the simulator, controllers, person, faults and arms are in [SIMULATION.md](SIMULATION.md). Literature keys resolve in [CITATIONS.md](../CITATIONS.md). Clinical dependencies are marked **[clinical review CRn]** and listed in [CLINICAL_REVIEW.md](../CLINICAL_REVIEW.md). Values we chose are marked **[assumption]** with the reason. "Our arithmetic" marks a number we derived from cited numbers. "SPIKE §n" is a section of the harness record, `docs/development/SPIKE_2026-10-28.md`.

**Status on 2026-10-28.** No faulted arm has run; no faulted-versus-clean outcome exists. The analysis code has run only on synthetic data with known effects and on clean-only harness units. `protocol-v1` is not yet tagged. D1-D19 were made by an AI agent (Claude, Anthropic) acting as study lead under authority the project lead, Pranav Reddy, delegated on 2026-10-28; the choices this file adds are marked *Added here* and were made by an AI agent under the same delegation. All are pending his review. No clinician has reviewed anything here.

## Contents

1. Unit of analysis and estimands (ADEMP)
2. Outcomes
3. Paired estimands
4. Statistics
5. The alarm layer
6. Outputs
7. Rules for reading results
8. Sensitivity analyses
9. Requirements on the code before the tag

---

## 1. Unit of analysis and estimands (ADEMP)

The study is reported in ADEMP form [Morris 2019].

**Aims.** The primary question and Q2-Q7 ([RESEARCH.md](RESEARCH.md) §2): the per-event harm, on true glucose, of single CGM fault events when unmodified OpenAPS oref0 0.7.1 acts on the faulty reading, with its timing, pathway, state dependence and mechanism.

**Data-generating mechanism** ([SIMULATION.md](SIMULATION.md); D1-D8, D16-D19).
- 30 virtual patients of simglucose 0.2.11 (10 children, 10 adolescents, 10 adults) [simglucose; Kovatchev 2009] through paper 1's minute-step engine [Reddy 2026], paper 1's default scenario mix.
- 3 controllers c: `oref0_tb` (primary), `oref1_smb`, `open_loop` (D3). ISF rule `rule1800` with the D17 scale map (`results/calibration/isf_scale_map.json`).
- Main seeds r = 101, ..., 110 (N = 10; D18), onset strata k = 0, ..., 9.
- A **trunk** is one continuous fault-free 67-h closed-loop run per (c, p, r). In each of 3 onset slots s (night, dose, any) an onset t0 is drawn; the full state is snapshotted at t0 and each **arm** a restarts from it for W = 720 min (12 h) **[assumption; clinical review CR7]**. The fault-free arm `clean` is the trunk's own continuation over [t0, t0 + W).
- An **event** (or **unit**) is one (c, p, r, s): one branch point with its clean arm and its fault and variant arms. Up to 30 × 10 = 300 events per (controller, arm).

**Arm names.** The analysis uses the harness's catalogue names (`schema.ARMS`; D6, D7). Fault arms, all three controllers:

| # | Arm | Slot | Tier | # | Arm | Slot | Tier |
|---|---|---|---|---|---|---|---|
| 1 | `pisa_typical` | night | 1 | 10 | `drift_low` | any | 1 |
| 2 | `pisa_severe` | night | 2 | 11 | `drift_low_severe` | any | 2 |
| 3 | `dropout_short` | any | 1 | 12 | `drift_high` | any | 1 |
| 4 | `dropout_long` | any | 2 | 13 | `drift_high_severe` | any | 2 |
| 5 | `dropout_severe` | any | 2 | 14 | `gain_low` | any | 2 |
| 6 | `frozen_typical` | any | 1 | 15 | `gain_high` | any | 2 |
| 7 | `frozen_severe` | any | 2 | 16 | `interference_small` | dose | 2 |
| 8 | `spike_typical` | any | 1 | 17 | `interference_large` | dose | 1 |
| 9 | `spike_severe` | any | 2 | | | | |

Variant arms, `oref0_tb` only, named `<base>+<variant>`: `pisa_typical+controller_only`, `pisa_typical+person_only`, `interference_large+controller_only`, `interference_large+person_only` (pathway); `clean+baysal_<set>` and `pisa_typical+baysal_<set>` for set in {aggressive, trial, nominal, cautious} (gate, night slot); `clean+matched_filter` and `pisa_typical+matched_filter` (second detector, night slot); `frozen_severe+fakecgm` (guard); `clean+correcting` and `interference_large+correcting` (correcting person, dose slot). The table's `clean` is the harness's `clean@<slot>` (`schema.table_arm`). N = 10 runs every tier (D18).

**Estimands** (D2).
- Individual effect: τ_a(c, p, r, s) = Y_a(c, p, r, s) − Y_clean(c, p, r, s), for an outcome Y of §2 over the window.
- Patient effect: τ̄_a(c, p) = mean over the valid seeds r of τ_a(c, p, r, s) (one slot per arm, so s is implied by a).
- Population estimand: θ_a(c) = (1/P_a) Σ_p τ̄_a(c, p), equal weight per patient, over the P_a patients with at least one valid seed for a (P_a = 30 when no patient is lost; §2.6).
- The patient is the unit of analysis. Seeds only reduce within-patient noise. CIs run across patients, never across seeds or events ([AGENTS.md](../../AGENTS.md) rule 8).
- The estimands are per event: the expected change in an outcome over 12 h caused by one fault event at an onset drawn as in the slot, averaged over this cohort's patients. They describe this simulator's cohort, not a population of people.

**Methods.** §3 and §4.

**Performance measures.** This is not a method-comparison simulation, so there is no bias or coverage target to score. The performance measures are each estimate's Monte Carlo error (§4.2), its CI width and the pairing efficiency (§4.6).

**The estimates, by role.**

| Role | What | Section |
|---|---|---|
| Confirmatory (H1, H2; Holm) | θ(`pisa_typical`, `oref0_tb`) on minutes above 180; θ(`interference_large`, `oref0_tb`) on minutes below 70 | §4.3 |
| Pre-specified primary comparison, estimation only (C1) | per-day harm of paper 1's five fault types in the alarm and insulin layers; Kendall's tau-b; compression's shares | §4.4 |
| Secondary, labelled, uncorrected | every θ_a(c); pairing efficiency; dose-response and ladders; pathway; mechanism; gate; guards; correcting person; controllers; person-year harm; additivity | §4.5-§4.16, §8.5 |
| Exploratory, labelled | state cells and the state model (Q4) | §4.8 |
| Sensitivity | S1, S2, C1's sensitivity analyses; standard mix (not run) | §8 |

## 2. Outcomes

### 2.1 What the analysis reads

One `outcomes.ArmTrace` per arm, built from the harness's `episodes.Unit` and `episodes.ArmRun` by `outcomes.trace_from_run`. Each covers [t0 − 120 min, t0 + W) (`PRE_MIN = 120` **[assumption: at least the 15 min the onset trends need; long enough for the alarms' trend fit and refractory period to carry over]**):
- 1-min series: true (plasma) glucose, the truth; basal or temp-basal insulin delivered; the person's meal boluses; correction boluses; SMBs; rescue carbohydrate.
- 5-min series: the CGM the person sees (the **displayed CGM**), the reading sent to the controller after any gate (the **controller input**), IOB, COB and oref0's `rT.reason` per loop cycle.
- Scalars: the scheduled basal (U/min), the minute of the person's last meal or correction bolus before t0, the fault's end, the realised fault parameters, the validity flag, `data_source` (synthetic or simulation).

`outcomes.slot_outcomes` turns the arms of one event into rows of the outcome table (`schema.COLUMNS`, one row per (controller, patient_id, seed, slot, arm)). It refuses an event whose arms do not share an identical history before t0. Patients carry paper 1's form `sim:<group>#<nnn>`.

### 2.2 Windows

All windows are half-open, on the 1-min grid, with minute 0 = t0. Glucose outcomes count 1-min truth samples.

- **Window:** [t0, t0 + W), W = 720 min.
- **Fault end** (D9): the fault ends at the last reading it changes. Compression and interference, whose offsets decay, end when the offset falls to 10% of its peak, as paper 1 defines the realised compression duration (SPIKE §4). `fault_end_min` is minutes after t0; `inf` when the fault is held for the window (drift #10-#13, gain #14-#15); NaN for an arm without a fault.
- **During the fault:** [t0, t0 + D), D = min(W, ⌈fault_end_min⌉) whole minutes (`ArmTrace.fault_minutes`).
- **After the fault:** [t0 + D, t0 + W).
- A fault still active at t0 + W is censored there: D = W, and every after-the-fault quantity is not observed (NaN) (D9). **Requirement on the code:** `outcomes.fault_windows` writes 0 for the after-the-fault columns when D = W; it must write NaN for the arm and its `ref_` twin.
- **Onset state** is read at t0 from the slot's clean arm (§2.5), because the reading at t0 may already be faulted.

### 2.3 Outcome definitions

Every outcome is on the truth, or on insulin and carbohydrate delivered, never on the CGM, except the alarm layer (§5), which scores alarms fired on the displayed CGM against true lows. "Effect" says how the paired effect is formed (§3): `diff` = arm minus the clean row of the same event; `ref` = arm minus the clean arm over the arm's own during or after minutes; `paired` = already paired.

**Glycaemia** (`outcomes.glycaemia`, on the truth over the window; thresholds strictly below 54 and 70, strictly above 180 and 250 mg/dL, paper 1's consensus module) [Battelino 2019; Battelino 2023] **[clinical review CR11]**.

| Column | Definition | Unit | Effect |
|---|---|---|---|
| `minutes_below_54` | #{minutes with g < 54} | min | diff |
| `minutes_below_70` | #{g < 70} | min | diff |
| `minutes_above_180` | #{g > 180} | min | diff |
| `minutes_above_250` | #{g > 250} | min | diff |
| `risk_weighted_min` | 3.0 × #{g < 54} + 2.4 × #{54 ≤ g < 70} + 1.6 × #{g > 250} + 0.8 × #{180 < g ≤ 250}: the glycaemia risk index (GRI) weights applied to minutes [Klonoff 2023] **[clinical review CR6]** | min (weighted) | diff |
| `mean_glucose`, `peak_glucose`, `nadir_glucose` | mean, maximum, minimum of g | mg/dL | diff |
| `hypo_l1_events`, `hypo_l2_events` | level 1 (< 70) and level 2 (< 54) hypoglycaemic events: at least 15 consecutive minutes below the threshold, ending after at least 15 consecutive minutes at or above it [Battelino 2023], found by paper 1's `detect_events` on the truth from t0 − 120 min to t0 + W; an event counts when its onset lies in the window (`outcomes.hypo_events`), so one already running at t0 does not count | count | diff |

Risk-weighted minutes are the one harm scale for ranking faults (C1, gate, Q7, state analysis). They add across events, so per-event effects can be turned into per-day (C1) and per-year (Q7) harm (D9). They are not the GRI itself, which is a capped percentage-based index.

**Insulin and carbohydrate** over the window (`outcomes.insulin`).

| Column | Definition | Unit | Effect |
|---|---|---|---|
| `insulin_total_u` | basal or temp + meal boluses + SMBs + corrections delivered | U | diff |
| `insulin_bolus_u` | the person's meal boluses | U | diff |
| `insulin_correction_u` | the person's correction boluses (correcting-person arms) | U | diff |
| `insulin_smb_u` | SMBs (`oref1_smb`; 0 otherwise) | U | diff |
| `insulin_tb_dev_u` | basal or temp delivered − scheduled basal × 720 min: the temp-basal deviation | U | diff |
| `rescue_carbs_g` | rescue carbohydrate eaten | g | diff |
| `rescue_count` | rescue treatments (minutes with rescue carbohydrate) | count | diff |
| `correction_boluses` | correction boluses given | count | diff |

Identity used as a test: the paired effect of `insulin_total_u` equals the sum of the paired effects of `insulin_tb_dev_u`, `insulin_smb_u`, `insulin_bolus_u` and `insulin_correction_u`, because the scheduled basal is the same in both arms of an event.

**During and after the fault** (`outcomes.fault_windows`; total insulin and truth above 180 mg/dL). Each column has a `ref_` twin holding the clean arm's value over the same minutes, so the effect compares like with like.

| Column | Definition | Unit | Effect |
|---|---|---|---|
| `insulin_during_u` | total insulin in [t0, t0 + D) | U | ref |
| `insulin_after_u` | total insulin in [t0 + D, t0 + W); NaN when D = W | U | ref |
| `minutes_above_180_during` | #{g > 180} in [t0, t0 + D) | min | ref |
| `minutes_above_180_after` | #{g > 180} in [t0 + D, t0 + W); NaN when D = W | min | ref |

**Time to re-converge** (`outcomes.reconvergence`). `reconverge_min` = s* − D, where s* is the first minute s ≥ D such that |g_a(m) − g_clean(m)| < 5 mg/dL for every m in [s, s + 60) **[assumption]**; the 60-min run must lie inside the window. Unit min; effect `paired` (it is a comparison with the clean arm already). `reconverged` is True when such a run exists. Not reached by t0 + W: censored, `reconverge_min` NaN and `reconverged` False. Not applicable (both NA) for an arm without a fault and for a fault with D = W.

**Controller actions** (`outcomes.actions`, `outcomes.classify_reasons`). For each shipped rule, the number of the 144 loop cycles in the window whose `rT.reason` matches it (a cycle counts once per rule). NaN in open loop. Effect `diff`.

| Column | Rule | Pattern (oref0 0.7.1 `lib/determine-basal/determine-basal.js`) |
|---|---|---|
| `actions_stale` | stale-data stop (reading older than 12 min) | "BG data is too old" (line 188) |
| `actions_flat` | flat-CGM check | "CGM data is unchanged" (lines 177-194) |
| `actions_temp_mismatch` | temp-mismatch cancel | "Warning: currenttemp ... canceling temp" (lines 381, 388) |
| `actions_lgs_zero_temp` | zero temp for a predicted low: bg or minGuardBG below the threshold | "minGuardBG x<y" (lines 912-920) |
| `actions_smb_maxdelta` | SMB disabled by maxDelta (> 20% of BG) | "maxDelta ... SMB disabled" (lines 877-879) |

*Added here* (ANALYSIS_NOTES O5): zero temps that correct an eventual BG below target are not counted as `lgs_zero_temp`; the maxDelta text is written whatever the SMB setting, so it is counted for `oref0_tb` too.

**Threshold crossings** (`outcomes.crossings`; descriptive, D9). `ctrl_below_70_during` and `ctrl_above_180_during`: whether any controller input at the 5-min samples t0, t0 + 5, ..., before t0 + D (at least the sample at t0; missing readings skipped) was below 70 mg/dL (oref0's suspend threshold) or above 180 mg/dL. NA in open loop and without a fault. Not differenced.

**Realised fault parameters** (`outcomes.fault_columns`). `fault_depth_mg_dl`, `fault_duration_min`, `fault_size`, `fault_tau_min`, `fault_end_min`, and `fault_area_mg_dl_h` = depth (mg/dL) × duration (min) / 60, the PISA artifact area used by the dose-response (§4.7).

**Alarm layer** (`alarm_layer.alarm_columns`): `<alarm>_<field>` for alarm in {trend, threshold, baysal_aggressive, baysal_trial, baysal_nominal, baysal_cautious, matched_filter} and field in {alarms, false_alarms, l1_caught, l1_missed, l2_caught, l2_missed}; counts, effect `diff`. Defined in §5.

### 2.4 Validity flag

`invalid` (`outcomes.is_invalid`) is True when any truth sample of the trace (history included) is below 10 mg/dL, paper 1's validity floor, outside the model's validity, or when the harness flagged the run. The clean arm's validity is judged on its truth up to t0 + W only, so a later dip of the trunk does not drop the slot (ANALYSIS_NOTES O8).

### 2.5 Onset state (Q4)

From the slot's clean arm (`outcomes.onset_state`), identical for every arm of the event:

| Column | Definition | Unit |
|---|---|---|
| `onset_truth_mg_dl` | g(t0) | mg/dL |
| `onset_truth_trend` | (g(t0) − g(t0 − 15)) / 15 | mg/dL/min |
| `onset_cgm_mg_dl` | the fault-free displayed CGM at t0 | mg/dL |
| `onset_cgm_trend` | (CGM(t0) − CGM(t0 − 15)) / 15 | mg/dL/min |
| `onset_iob_u`, `onset_cob_g` | IOB and COB of the trunk's loop cycle at t0 (SPIKE §5) | U, g |
| `onset_hours_since_bolus` | (t0 − the person's last meal or correction bolus) / 60; SMBs excluded | h |
| `onset_hour` | hour of day at t0 | h |

### 2.6 Invalid arms

The rule (D11): an invalid arm drops the whole (patient, seed, slot) for every arm of that slot, and "every arm" includes every controller's arms, so all controllers keep the same events. A patient left with no valid seed for an arm is left out of that arm's estimate (P_a < 30), and the count is reported. Counts are reported by controller, slot and arm in `invalid_counts.csv` (rows, invalid arms, dropped rows).

**Requirement on the code:** `effects.apply_invalid_rule` groups by (controller, patient, seed, slot), so it drops a slot within one controller only. It must group by (patient, seed, slot) across the three controllers, as D11 says, and report both the triggering controller and the rows dropped per controller. With no invalid arm the two rules agree.

## 3. Paired estimands

**Event-level effect** (`effects.paired_effects`), for every outcome column o, row (c, p, r, s, a):
- `diff`: τ = Y_a − Y_clean, with Y_clean the clean row of the same (c, p, r, s);
- `ref`: τ = Y_a − Y_a^ref, where Y_a^ref (`ref_<o>`) is the clean arm over the arm's own during or after minutes;
- `paired`: τ = Y_a (time to re-converge).
The clean row's own `diff` effects are 0.

**Patient effect** (`effects.patient_stats`): τ̄_a(c, p) = (1/n_p) Σ_r τ_a(c, p, r), over the n_p valid seeds with a finite value (NaN dropped); also n_p and the SD over seeds.

**Population estimate** (`effects.estimate`; `effects.summarise` for every controller, slot, arm and outcome): θ̂_a(c) = mean over patients with a value of τ̄_a(c, p).

**Contrasts between arms of one slot** (`secondary.combine`, `effects.arm_contrast`). Arms sharing a slot share t0 and the snapshot, so a contrast is formed event by event and then averaged as above: for weights w_a over arms of one slot, the event value is Σ_a w_a τ_a(c, p, r, s), on the events that have every arm named (a zero weight only restricts to those events).

**Controller contrast** (`effects.controller_contrast`): for the same (p, r, s, a), Δ = τ_a(oref1_smb, p, r, s) − τ_a(oref0_tb, p, r, s). The trunks differ by controller, so this is a difference of two paired effects, not a within-trunk contrast.

## 4. Statistics

### 4.1 The stratified patient bootstrap

Paper 1's `ghostlows.bench.stats.bootstrap_ci` [Reddy 2026], used through one joint resampling index (`effects.Bootstrap`):
1. Patients are put in sorted order and grouped by age group (child, adolescent, adult; 10 each).
2. With `rng = numpy.random.default_rng(seed)`, seed = 0 **[assumption: paper 1's default]**, each replicate b = 1, ..., B draws, within each age group in turn, as many patients as the group holds, uniformly with replacement. `Bootstrap` records the drawn patients as a count matrix n_{b,p} (how often patient p is drawn in replicate b).
3. A replicate of a mean over patients is θ*_b = Σ_p n_{b,p} v_p 1[v_p finite] / Σ_p n_{b,p} 1[v_p finite] (`Bootstrap.means`): patients without a value are skipped.
4. A percentile CI at level L takes numpy's linear-interpolation quantiles (1 − L)/2 and (1 + L)/2 of the finite replicates (`effects.percentile_ci`), exactly as `bootstrap_ci` does; equality with `bootstrap_ci` is a test.

B = 10,000 for H1, H2 and C1 (one index shared by all three and by C1's sensitivity analyses); B = 2,000 for everything else (one other index shared by every secondary and exploratory analysis, so the gate's two optima can be compared replicate by replicate) (D11; ANALYSIS_NOTES E3, E5). Shares, rank correlations, sums over faults and contrasts are computed per replicate from the replicate means on the same index. Stratifying keeps 10 patients per age group in every replicate, as in paper 1.

### 4.2 Monte Carlo error

For every estimate (`effects.estimate`), with SE_p = SD_p / √n_p over the seeds of patient p (patients with n_p ≥ 2):
- within-patient SE = mean over patients of SE_p (`within_se`);
- between-patient SD = SD (ddof 1) of the patient means τ̄(p) (`between_sd`);
- `mc_se` = √(mean_p SE_p² / P), the seed noise left in θ̂;
- `mc_flag` = within_se > 0.5 × between_sd. A flagged primary estimate (H1, H2, C1's ten per-fault per-day harms) is reported as a limitation (D11).

### 4.3 The confirmatory family H1-H2

`hypotheses.run` on `oref0_tb`, B = 10,000. Definitions in `hypotheses.HYPOTHESES`:

| ID | Arm | Outcome | SESOI S_k | Co-reported |
|---|---|---|---|---|
| H1 | `pisa_typical` (#1) | `minutes_above_180` | 72 min per event: the 5-percentage-point change in time in range (about 1 h 12 min per day) that the consensus calls clinically meaningful [Battelino 2019] **[assumption: per-event reading; clinical review CR1]** | `minutes_below_54`, `minutes_below_70` |
| H2 | `interference_large` (#17) | `minutes_below_70` | 15 min per event, about 1% of a day **[assumption: between the consensus targets of under 1% below 54 and under 4% below 70 mg/dL [Battelino 2019]; clinical review CR1]** | `minutes_below_54` |

**Bootstrap p-value** (`hypotheses.bootstrap_p`). From the B' finite replicates θ*_b of θ_k: a_k = (1 + #{θ*_b ≤ 0}) / (B' + 1), b_k = (1 + #{θ*_b ≥ 0}) / (B' + 1), p_k = min(1, 2 min(a_k, b_k)).

**Holm** across H1 and H2 at familywise α = 0.05 [Holm 1979] (`hypotheses.holm`). Order by p (ties by ID, so H1 first; a NaN p goes last and never rejects). Step 1 is rejected when p_(1) ≤ 0.025; step 2 only if step 1 was rejected and p_(2) ≤ 0.05. With B' = 10,000 a p-value is 2j/10,001 for an integer j and never equals 0.025 or 0.05 exactly (our arithmetic), so "≤" and "<" give the same decisions unless replicates are lost.

**CIs.** The estimate's CI is at the Holm-adjusted level of its step: 97.5% at step 1, 95% at step 2, whether or not rejected. The 90% percentile CI of the same replicates, unadjusted, serves the equivalence reading.

**The four-way rule** (`hypotheses.decide`, applied in this order):
1. *Not estimable* if θ̂_k or either 90% CI bound is not finite.
2. *Meaningful harm* if Holm rejects θ_k = 0, θ̂_k > 0 and θ̂_k ≥ S_k.
3. *Opposite direction* if Holm rejects, θ̂_k < 0 and θ̂_k ≤ −S_k: reported as an unexpected direction, never as harm or benefit of the fault.
4. *Negligible* if the 90% CI lies strictly inside (−S_k, +S_k): two one-sided tests, each at α = 0.05 [Lakens 2017], not adjusted for multiplicity.
5. *Inconclusive* otherwise.

The readings exclude each other: harm needs θ̂ ≥ S, negligible needs the 90% upper bound, and so θ̂, below S; opposite needs θ̂ ≤ −S, negligible needs θ̂ > −S. The Holm test, not a plain 95% CI, decides harm.

**Co-reported outcomes** use the same 10,000 replicates, 95% CIs, labelled "co-reported, 95% CI, uncorrected" (`hypotheses_co_reported.csv`).

### 4.4 C1, the cross-layer comparison (pre-specified primary comparison; estimation only)

`cross_layer.run` with `cross_layer.specs()`, on the B = 10,000 index.

**Faults, arms and rates** (`cross_layer.RATES`; D10). Per-day harm of fault f in layer L:

H_L(f) = Σ_{a ∈ f} rate_a × θ_a(c_L, o_L),

| Fault f | Arms a and rate_a (events per person-day) | Source |
|---|---|---|
| compression | `pisa_typical` × 0.86 | compression on 86% of nights, one per night, the fraction at which paper 1's PISA model reproduces the published 0.321 PISAs per day with a sensor minimum below 85 mg/dL [Reddy 2026; Lobo 2026] |
| dropout | `dropout_short` × 2.0 | paper 1's 4 per day at s = 1, thinned by s = 0.5 [Reddy 2026]; no source for the rate itself (paper 1 marks it **[assumption]**) |
| frozen | `frozen_typical` × 1.0 | paper 1's 2 per day × 0.5, as above |
| spikes | `spike_typical` × 1.0 | paper 1's 2 per day × 0.5, as above |
| drift | `drift_low` × 1.0 + `drift_high` × 1.0 | drift as a continuous state: two 12-h windows per day, half in each direction, as paper 1 draws the sign with probability 1/2 [Reddy 2026] (our arithmetic) |

**Layers** (`cross_layer.LAYERS`):
- alarm layer: c = `open_loop`, o = `trend_false_alarms` (excess false alarms of paper 1's trend alarm, §5), so H is in false alarms per day;
- insulin layer: c = `oref0_tb`, o = `risk_weighted_min`, so H is in risk-weighted minutes per day.

**Ranks.** Within an analysis and layer, rank 1 = the largest H (pandas average ranks for ties).

**Kendall's tau-b** [Kendall 1945] between the two layers' five per-day harms (`cross_layer.kendall_tau_b`): with x_f = H_alarm(f), y_f = H_insulin(f), rounded to 9 decimals so that values equal to 1e-9 rank as ties **[assumption: rounding, not effects]**,

τ_b = Σ_{i<j} sgn(x_i − x_j) sgn(y_i − y_j) / √(n_x n_y),

where n_x and n_y are the numbers of the 10 pairs not tied in x and in y. NaN if either layer is all ties or a value is missing. Its CI is the 95% percentile CI of τ_b computed per replicate from the replicate per-day harms.

**Compression's share** of layer L: share_L = H_L(compression) / Σ_f H_L(f), with replicates from the replicate harms. A share is **estimable** only when Σ_f H*_L(f) > 0 in at least 97.5% of all B replicates (a NaN total counts as not positive); otherwise it is reported as not estimable (D10; `ESTIMABLE_SHARE`). A share can lie outside [0, 1] when some fault's per-day harm is negative; it is reported as computed.

**Discordance.** When both shares are estimable: difference = share_alarm − share_insulin, with its 95% percentile CI from the per-replicate differences. *Discordant* if that CI excludes 0; NA when either share is not estimable.

**Patients.** Per replicate, θ*_a is a mean over the drawn patients that have a value for a (§4.1), identical to the ratio of patient-level sums when every patient is present for every fault (D11; ANALYSIS_NOTES C-a). The per-fault table reports `n_patients` = the minimum P_a over the fault's arms.

C1 has no test and no decision rule beyond "discordant". Its sensitivity analyses are in §8.3.

### 4.5 All per-arm effects

`effects.summarise` gives θ̂_a(c) with its 95% CI (B = 2,000) and Monte Carlo error for every controller, slot, non-clean arm and outcome (`effects_summary.csv`). This is the table behind the forest plots. Secondary, uncorrected.

### 4.6 Pairing efficiency (N1)

`secondary.pairing_efficiency`, on the outcome table after the invalid rule, for every controller, non-clean arm and the primary outcomes `minutes_above_180`, `minutes_below_70`, `risk_weighted_min` and `trend_false_alarms` (those of H1, H2 and C1) **[assumption: D11's "primary outcome" read as the outcomes of H1, H2 and C1]**. For patient p, with its n valid seeds in increasing order r_1, ..., r_n and i + 1 taken mod n:
- d_i = Y_a(p, r_i) − Y_clean(p, r_i) (paired);
- x_i = Y_a(p, r_i) − Y_clean(p, r_{i+1}) (cross-seed, standing in for independent runs);
- V_p = Var(d), V'_p = Var(x) (ddof 1); ratio_p = V_p / V'_p, a patient with V_p = V'_p = 0 left out.
Reported: the median ratio over patients, and the implied number of independent runs per paired run = 1 / median ratio (∞ when the median is 0), with the medians of V_p and V'_p. No CI. The cross-seed difference also contains the between-seed variation of the clean arm, which is what independent runs would add.

### 4.7 Dose-response and level ladders (N8)

**Dose-response** (`secondary.dose_response`), for `pisa_typical`, every controller, outcomes `minutes_above_180` and `risk_weighted_min`. Per patient, with x = `fault_area_mg_dl_h` and y = the event effect over the patient's valid events:

b_p = Σ (x − x̄_p)(y − ȳ_p) / Σ (x − x̄_p)²,

defined when Σ (x − x̄_p)² > 1e-12 (patients with fewer than two distinct areas are left out). Reported: the mean of b_p over patients with a value, in minutes per mg/dL·h, with its 95% bootstrap CI (B = 2,000) and the patient count.

**Crossing split** (same function; closed loops only): θ for `pisa_typical` events split by `ctrl_below_70_during` (did the controller's input fall below the 70 mg/dL suspend threshold during the fault), each group a mean over the patients with events in it.

**Ladders** (`secondary.ladders`; `LADDERS`), every controller, outcomes `risk_weighted_min`, `minutes_below_70`, `minutes_above_180`. Each rung's θ against clean, and each rung against the rung before it on shared events (`combine` with weights +1, −1):

| Ladder | Rungs |
|---|---|
| compression | `pisa_typical` → `pisa_severe` |
| dropout | `dropout_short` → `dropout_long` → `dropout_severe` |
| frozen | `frozen_typical` → `frozen_severe` |
| spikes | `spike_typical` → `spike_severe` |
| drift low | `drift_low` → `drift_low_severe` (−20%) → `gain_low` (−40%) |
| drift high | `drift_high` → `drift_high_severe` (+20%) → `gain_high` (+40%) |
| interference | `interference_small` → `interference_large` |

*Added here* (ANALYSIS_NOTES S3): the dropout ladder includes the G6 long-gap arm between the short and the severe gap.

### 4.8 State dependence (Q4; exploratory)

On `oref0_tb`'s 17 fault arms and `risk_weighted_min`, per (slot, arm); *Added here:* `oref1_smb` and `open_loop` are not analysed for Q4, because Q4 asks about the primary controller and every further controller would multiply the exploratory cells (`state_cells(controllers=("oref0_tb",))`, the default) **[clinical review CR18]**.

**Cells** (`secondary.state_cells`):
- glucose × trend (9 cells): true glucose at onset < 100, 100-140 (both ends included), > 140 mg/dL × displayed fault-free CGM trend falling (< −1 mg/dL/min), flat (−1 to +1), rising (> +1);
- IOB tertile: T1 ≤ the 1/3 quantile, T2 ≤ the 2/3 quantile, T3 above, the quantiles taken over the events of the same controller and slot (*Added here* **[assumption: tertiles of the analysed events give three cells of equal size; IOB has no consensus bands]**);
- night or day: night when the onset hour is in [23, 24) or [0, 7) **[assumption: paper 1's PISA window]**.
Each cell: θ with 95% CI (B = 2,000), counts of events and patients. A cell with fewer than 5 patients keeps its counts and is reported as not estimable (estimate and CI empty).

**Fixed-effects model** (`secondary.state_model`). For each arm, on the events with all values present: y_e = α_p + β′x_e + ε_e, with x = (`onset_truth_mg_dl`, `onset_cgm_trend`, `onset_iob_u`, `onset_cob_g`, `onset_hours_since_bolus`). With X̃_p and ỹ_p demeaned within patient, S_p = X̃_p′X̃_p and s_p = X̃_p′ỹ_p:

β̂ = (Σ_p S_p)⁻¹ Σ_p s_p (least squares), and per replicate β*_b = (Σ_p n_{b,p} S_p)⁺ Σ_p n_{b,p} s_p (pseudo-inverse),

so a resampled patient's events count as often as the patient is drawn. A covariate whose diagonal of Σ_p S_p is at most 1e-9 × max(1, max |S|) (no within-patient variation) is left out. Reported: coefficients with 95% percentile CIs, events and patients.

**Held-out check** (same function). The model without patient effects: design = intercept, the five covariates and age-group indicators (the first group alphabetically, adolescent, as reference). Paper 1's five leave-patients-out folds (`ghostlows.bench.stats.patient_folds`, stratified by age group, seed 0). For each fold k with test events and more training events than columns, ordinary least squares on the training folds; SSE_k = Σ_test (y − ŷ)², SST_k = Σ_test (y − ȳ_train)². Reported: out-of-fold R² = 1 − Σ_k SSE_k / Σ_k SST_k, the fold minimum and maximum of 1 − SSE_k / SST_k, and the folds used. R² ≤ 0 means the onset state predicts held-out patients' effects no better than the training mean.

### 4.9 Pathway 2 × 2 and the wake probability (Q3)

`secondary.pathway`, `oref0_tb`, for `pisa_typical` and `interference_large`, outcomes `risk_weighted_min`, `minutes_above_180`, `minutes_below_70`, `minutes_below_54`. The four cells of one event: fault-free (clean), controller-only (the controller sees the faulted CGM, the person the fault-free CGM), person-only (the reverse), both (the default fault arm). On the events that have all three faulted arms, so the components add up exactly:
- total = τ_both;
- controller-only = τ_ctrl;
- person-only = τ_person;
- interaction = τ_both − τ_ctrl − τ_person.
Each with θ and 95% CI (B = 2,000).

**Wake probability (derived; PISA typical only).** In the default arms the person wakes to every night alarm (night response probability 1.0; D5 **[assumption; clinical review CR2]**). Suppose instead that for each event the person independently wakes to a false night low with probability q and otherwise sleeps through it, while responses to true lows are unchanged **[assumption: one wake decision per event]**. If the person wakes, the event follows the both arm; if not, the person acts as if they saw the fault-free CGM, which is the controller-only arm. The expected effect of the event is then

E[τ | q] = q τ_both + (1 − q) τ_ctrl,

and because patient means and θ are linear, θ(q) = q θ_both + (1 − q) θ_ctrl. It is computed event by event on the events that have both arms (`combine` with weights q and 1 − q) for q = 0, 0.25, 0.5, 0.75 and 1, with 95% CIs. It is an approximation: a person who sleeps through the false low but would have woken for a true low later in the same night is not modelled separately.

### 4.10 Mechanism (N8)

`secondary.mechanism`, the 17 fault arms of every controller: θ and 95% CI of `insulin_during_u`, `insulin_after_u`, `minutes_above_180_during`, `minutes_above_180_after` (all `ref` effects, §3), `rescue_carbs_g`, `rescue_count` and the five controller-action counts. For the closed loops also the share of events whose controller input crossed below 70 or above 180 mg/dL during the fault (`share_ctrl_below_70_during`, `share_ctrl_above_180_during`): per patient the mean of the event indicator, then the mean over patients with its CI (descriptive).

### 4.11 Gate (Q5) and its least-harm against fewest-false-alarms setting (N10)

`secondary.gate`; harm on `oref0_tb`, alarms on `open_loop`; B = 2,000; w = 0.86, the fraction of nights with compression (`RATES["pisa_typical"]`) [Reddy 2026].

**Harm side**, on `risk_weighted_min`, for each detector setting g in {baysal_aggressive, baysal_trial, baysal_nominal, baysal_cautious, matched_filter}, event by event:
- on a compression night: Δ_P(g) = τ(`pisa_typical+g`) − τ(`pisa_typical`);
- on a night without compression: Δ_C(g) = τ(`clean+g`);
- net effect per night: N(g) = w θ(Δ_P(g)) + (1 − w) θ(Δ_C(g)), with replicates formed the same way. Negative is a benefit. No gate has N = 0 by definition.

**Alarm side**, on the `open_loop` night-slot arms `pisa_typical` and `clean`, levels (not differences), for the ungated trend alarm (no gate) and paper 1's gated trend alarm with each setting (§5). For a column X, the per-night value is w θ_level(X, pisa_typical) + (1 − w) θ_level(X, clean), with θ_level the mean over patients of the patient-mean level:
- false alarms per night from `<alarm>_false_alarms`;
- level 2 sensitivity = C2 / (C2 + M2), with C2 and M2 the per-night weighted caught and missed level 2 lows (`<alarm>_l2_caught`, `<alarm>_l2_missed`), so caught and missed lows of the two kinds of night are pooled with the same weights (D11); level 1 sensitivity likewise.

**Optima** among the five options {no gate, aggressive, trial, nominal, cautious} (in that order):
- *harm-optimal*: the lowest N; ties to the earlier option;
- *alarm-optimal*: paper 1's operating-point rule (`ghostlows.bench.stats.choose_threshold`, δ = 0, its pre-registered margin [Reddy 2026]): among options whose level 2 sensitivity is at least that of no gate (minus 1e-12), the fewest false alarms per night, ties to the higher sensitivity, then the earlier option. No gate is the reference, so a feasible option always exists. If the drawn nights hold no level 2 low (weighted C2 + M2 = 0 for no gate), level 1 sensitivity is used instead (paper 1's amendment A7).
- Both optima are found on the estimate and in every replicate. Reported: both curves (N, false alarms per night, level 2 and level 1 sensitivity, each with 95% CI), the two optima, whether the alarm rule was feasible and at which level, and the **coincidence share**: the fraction of replicates, among those in which both optima are defined, in which the two optima are the same setting.
- The matched filter is reported beside the curves, outside the five-way choice (it has one setting).
- *Added here:* both optima are **in-sample descriptive quantities**. Paper 1's `choose_threshold` is meant for training data only; here it is applied to the same events it is scored on, and so is the harm-optimal choice. Neither optimum is a recommended or validated setting, and neither is reported as one.

### 4.12 Guards (Q5)

`secondary.guards`, event by event, outcomes `risk_weighted_min`, `minutes_below_70`, `minutes_above_180`, `insulin_total_u`:
- flat-CGM check (`oref0_tb` only): τ(`frozen_severe+fakecgm`) − τ(`frozen_severe`): the harm when the check cannot act, minus the harm when it can. The check also fires on 3.6-5.6% of clean loop cycles (SPIKE §3.6), so this contrast measures the check's whole effect in the window, not only its effect on the freeze (D7).
- stale-data stop (`oref0_tb` and, *Added here*, `oref1_smb`, since both arms exist for both closed loops and the stop is shared code): τ(`frozen_severe`) − τ(`dropout_severe`) on the same t0: the harm of acting on frozen values with fresh timestamps instead of stopping after 12 min.

### 4.13 Correcting person (N7; tier 3)

*Added here* in D11: the correcting person's interference effect = τ(`interference_large+correcting`) − τ(`clean+correcting`), event by event on the dose slot (the excess harm of the false high when the person corrects at a displayed CGM ≥ 250 mg/dL), set beside the default person's θ(`interference_large`) and θ(`clean+correcting`) alone, on `risk_weighted_min`, `minutes_below_54`, `minutes_below_70`, `minutes_above_180`, `insulin_correction_u` and `correction_boluses`; 95% CIs, B = 2,000 **[clinical review CR4]**.

*Added here:* beside the contrast, the share of events with at least one correction bolus is reported for each correcting-person arm (per patient the mean of the indicator `correction_boluses` > 0, then the mean over patients with its 95% CI). In a harness check (one patient, calibration seed, controller actions only) interference large and a +40% gain never brought the displayed CGM to 250 mg/dL within 12 h under oref0, so the rule did not act (SPIKE §6; CR4); the contrast may therefore be 0 by construction. If the share is 0 for every patient, the contrast is reported as "the correcting-person rule did not act", never as an effect of correcting. Whether the 250 mg/dL trigger stays is decided before the tag and recorded in DECISIONS.md; without such a decision it stays as D5 sets it.

**Requirement on the code:** not implemented at df50a9f. Add `secondary.correcting_person(effects, boot)` built on `combine`, with the correction share above, called from `secondary.run`, written to `correcting_person.csv`.

### 4.14 Controllers (Q6)

`secondary.controllers`:
- for each of the 17 fault arms, θ of Δ = τ_a(`oref1_smb`) − τ_a(`oref0_tb`) (§3) on `risk_weighted_min`, `minutes_below_54`, `minutes_below_70`, `minutes_above_180` and the insulin parts `insulin_total_u`, `insulin_smb_u`, `insulin_tb_dev_u`, `insulin_bolus_u`, `insulin_correction_u`;
- the SMB and temp-basal split: for `oref1_smb` on the false-high arms (`drift_high`, `drift_high_severe`, `gain_high`, `interference_small`, `interference_large`) and the spike arms (`spike_typical`, `spike_severe`), θ of the paired effect of each insulin part. By the identity in §2.3, the parts sum to the excess total insulin.

### 4.15 Expected excess harm per person-year (Q7)

`secondary.person_year`, for every controller, on `risk_weighted_min` and `minutes_below_54` (D11's two outcomes) and also `minutes_below_70` and `minutes_above_180`:

PY(c, o, scenario) = 365 × Σ_a rate_a × θ_a(c, o),

with replicates formed the same way and 95% percentile CIs (B = 2,000). Scenarios (`secondary.SCENARIOS`; rates per person-day):

| Arm | Central | Low | High |
|---|---|---|---|
| `pisa_typical` (fraction of nights) | 0.86 | 0.86 × 0.03 / 0.321 = 0.080 | 0.86 × 0.33 / 0.321 = 0.884 |
| `dropout_short` | 2.0 | not used | not used |
| `dropout_long` | not used | 0.13 | 0.24 |
| `frozen_typical` | 1.0 | 1.0 | 1.0 |
| `spike_typical` | 1.0 | 1.0 | 1.0 |
| `drift_low`, `drift_high` | 1.0 each | 1.0 each | 1.0 each |

Sources: central = paper 1's rates, as in C1 [Reddy 2026]; compression low and high = the lowest and highest published detected-PISA rates, 0.03 per person-day [LoopEval 2026] and about 0.33 per day [Idi 2025], mapped to a fraction of nights by paper 1's calibration (0.86 nights reproduce 0.321 per day [Lobo 2026]) **[assumption: the published detection definitions are comparable with Lobo 2026's; clinical review CR13]** (our arithmetic); long gaps 0.13 to 0.24 per day [Vettoretti 2023], in place of paper 1's short gaps, not added to them (D11). Frozen values, spikes and drift keep paper 1's rates in every scenario, since no published rates exist. Interference and gain errors are per exposure, with no population rate, and are left out. For `oref0_tb` on `risk_weighted_min`, the central value equals 365 × C1's insulin-layer total; this identity is a test (up to the different B).

The estimate assumes the harm of a day's events is the sum of their single-event harms; the additivity check (§8.5) tests that. The θ uncertainty comes from the patient bootstrap; the rate uncertainty is shown by the scenarios, not by a CI.

**Requirement on the code:** D11 says each fault's term is reported with the total. `person_year` reports only the total and the rates. Add one row per (controller, scenario, outcome, arm) with 365 × rate_a × θ_a and its CI.

### 4.16 Labels

Everything outside H1 and H2 is reported with 95% CIs and labelled: C1 "pre-specified primary comparison, estimation only"; §4.5-§4.15 "secondary, uncorrected", except §4.8 "exploratory"; §8 "sensitivity analysis". No secondary or exploratory estimate is called significant.

## 5. The alarm layer

Paper 1's alarms, truth and scoring, imported unchanged from `ghostlows` v0.4.0 [Reddy 2026] (paper 1's BENCHMARK.md §1-4), run on each arm's displayed CGM and scored against that arm's own truth (`alarm_layer.score_traces`, `alarm_layer.alarm_columns`).

**Alarms** (`alarm_layer.make_alarms`), each at its default threshold:
- `trend`: linear extrapolation over the last 15 to 20 min predicts below 55 mg/dL within 20 min (a stand-in for Dexcom Urgent Low Soon, whose model is not public [Garg 2022]);
- `threshold`: fires when the CGM is below 70 mg/dL;
- both with paper 1's ungated urgent-low floor (fires on any reading ≤ 55 mg/dL) **[clinical review CR20]** and its 30-min refractory period **[assumption, paper 1]**; floor alarms are included in the counts;
- gated trend alarms: paper 1's `GatedAlarm` with the causal Baysal detector [Baysal 2014] at each of its four parameter sets (`baysal_<set>`), and with the causal matched filter (a causal adaptation of [Idi 2023]) (`matched_filter`), each at its default threshold 0.5, the published rule.

**Which arms.** `trend` and `threshold` on every arm of every controller. The gated alarms on the `open_loop` night-slot arms only (`alarm_layer.gated_by_default`; D9); their columns are NaN elsewhere. The alarm layer reads the displayed CGM, which is what the person sees: in a controller-only pathway arm it is the fault-free CGM, so that arm's alarm effects are 0 by construction; in the gate arms the person sees the raw CGM, so their alarm columns are those of the raw displayed CGM.

**Inputs.** For each arm of an event: the 5-min displayed CGM and the 1-min truth over [t0 − 120 min, t0 + W), on paper 1's clock (`alarm_layer._frames`). The alarms run over this whole trace, so the trend fit, the Baysal Kalman estimate and the refractory period carry over from the 120 min before t0, identically in every arm of the event.

**True lows** (`ghostlows.bench.events.true_events`; paper 1 BENCHMARK.md §3): on the 1-min truth, at least 15 consecutive minutes below 70 (level 1) or below 54 mg/dL (level 2), ending after at least 15 consecutive minutes at or above the threshold [Battelino 2023]. Level 2 events lie inside level 1 events; each level is scored separately; false-alarm matching uses level 1 events.

**Matching** (`ghostlows.bench.scoring.fire_and_score` with `DEFAULT_WINDOWS`; paper 1 BENCHMARK.md §4, its D12):
- an alarm matches the first true low with the alarm time in [onset − 30 min, event end];
- a low is **caught** if a matched alarm falls in [onset − 30 min, onset + 15 min] **[assumption on windows, paper 1]**; matched alarms after onset + 15 min are late, not false;
- an alarm is **false** when it matches no true level 1 low.

**What is counted per arm** (one row each, effect `diff`):
- `<alarm>_alarms`: alarms fired at times in [t0, t0 + W);
- `<alarm>_false_alarms`: those of them that are false;
- `<alarm>_l{1,2}_caught`, `<alarm>_l{1,2}_missed`: true lows of each level with onset in [t0, t0 + W), caught or not.

**Excess false alarms and missed lows per event** are the paired effects: τ = count in the arm − count in the clean arm of the same event (§3). The alarm-layer per-day harm of C1 is θ(`trend_false_alarms`, `open_loop`) × rate (§4.4).

*Added here* (ANALYSIS_NOTES A1, A2):
- False alarms are not split by cause, as paper 1 does (compression, other fault, near-event, noise): the paired difference already attributes the excess to the fault, because the clean arm carries every other cause.
- Edge rule: an alarm fired late in the window before a low that starts after t0 + W counts as false. The edge is the same for every arm of an event, so it cancels in expectation in the paired difference.
- Paper 1's sensitivity analyses of the matching windows (early 20/30/45 min, late 15/30 min) are not run.

## 6. Outputs

`python experiments/analyse.py --input <outcome tables> --protocol-tag protocol-v1` (`experiments/analyse.py`). A table with a faulted arm whose `data_source` is not `synthetic` is analysed only with `--protocol-tag protocol-v1`, when that tag exists in the repository and HEAD contains it; otherwise the run exits with status 2 before computing or writing anything (ANALYSIS_NOTES P1). A clean-only simulation table gets the descriptive tables only (`invalid_counts`, `clean_levels`). `--synthetic` writes to `results/synthetic/`, with "SYNTHETIC DATA: a test of the analysis code, not a result." at the top of `summary.md`.

**Tables**, CSV in `results/analysis/` (`--out` overrides):

| File | Content | Section |
|---|---|---|
| `invalid_counts.csv` | rows, invalid arms and dropped rows by controller, slot and arm | §2.6 |
| `clean_levels.csv` | the clean arm's outcome levels (no contrast) per controller and slot, θ with CI; face validity **[clinical review CR17]** | §1 |
| `effects_summary.csv` | θ, 95% CI and Monte Carlo error for every controller, slot, arm and outcome | §4.5 |
| `hypotheses.csv` | H1, H2: θ, p, Holm step and level, rejection, Holm-level CI, 90% CI, decision, Monte Carlo error | §4.3 |
| `hypotheses_co_reported.csv` | co-reported outcomes, 95% CI, uncorrected | §4.3 |
| `c1_per_fault.csv` | per analysis, layer and fault: arms, rates, per-day harm, CI, rank, patients | §4.4, §8.3 |
| `c1_summary.csv` | per analysis: tau-b with CI, totals, share estimability, shares with CIs, share difference with CI, discordant | §4.4, §8.3 |
| `pairing_efficiency.csv` | median variance ratio and implied runs | §4.6 |
| `dose_response.csv`, `dose_by_crossing.csv` | mean per-patient slopes; effects by suspend-threshold crossing | §4.7 |
| `ladders.csv` | rungs against clean and against the rung before | §4.7 |
| `state_cells.csv`, `state_model.csv`, `state_held_out.csv` | Q4 cells, fixed-effects coefficients, held-out R² | §4.8 |
| `pathway.csv`, `wake_curve.csv` | the 2 × 2 components; θ(q) | §4.9 |
| `mechanism.csv` | during and after the fault, rescue, actions, crossing shares | §4.10 |
| `gate_curves.csv`, `gate_optima.csv` | per option: Δ_P, Δ_C, net per night, false alarms per night, sensitivities; the two optima and the coincidence share | §4.11 |
| `guards.csv` | the two guard contrasts | §4.12 |
| `controllers.csv`, `smb_tb_split.csv` | `oref1_smb` − `oref0_tb`; insulin parts | §4.14 |
| `person_year.csv` | per controller, scenario and outcome | §4.15 |
| `summary.md`, `manifest.json` | the headline tables, file list, commit, arguments, n_boot, bootstrap seed, protocol tag, run time | |

Written since df50a9f: `additivity.csv` (with `--additivity results/additivity`, §8.5); the S1 and S2 tables, each a full analysis written to its own folder (`results/analysis_s1/`, `results/analysis_s2/`, §8.1, §8.2) *Added here*. Required and not yet written (§9): `correcting_person.csv` (§4.13); the per-fault rows of `person_year.csv` (§4.15); the S1-against-primary tau-b (§8.1); Δ as a share of the observed difference in `additivity.csv` (§8.5).

**Figures**, PNG (200 dpi) and SVG, in `results/figures/` (`figures.all_figures`), each with the CSV behind it:

| File stem | Content |
|---|---|
| `forest_risk_weighted_min`, `forest_minutes_above_180`, `forest_minutes_below_70` | θ with 95% CI for every fault arm, one marker per controller (`figures.forest`) |
| `c1_layers_main` | C1's main analysis: per-day harm in both layers, each fault's rank (`figures.c1_layers`) |
| `dose_response_minutes_above_180` | `pisa_typical` event effects on `oref0_tb` against artifact area, with the mean per-patient slope through the centroid (`figures.dose_response`) |
| `pathway_risk_weighted_min` | the pathway 2 × 2 with its interaction (`figures.pathway`) |
| `gate_curves` | net effect per night, false alarms per night and level 2 sensitivity for the five options, the matched filter beside them (`figures.gate_curves`) |

Colours: three palette slots, one per controller, checked for colour-vision deficiency; marker shapes and a legend carry the same information (ANALYSIS_NOTES F1).

**Stored per unit** by the harness, not by the analysis: one outcome row per arm, and 5-min traces (truth, displayed CGM, controller input, insulin, carbohydrate) for a fixed 10% of (patient, seed) pairs drawn with seed 20261028 for every controller (D9, D19), for figures of example events and audits.

## 7. Rules for reading results

**What may be said.**
- Every effect is "per event, in this simulator's cohort of 30 virtual patients, under unmodified oref0 0.7.1 (or its SMB mode) with the settings of D16-D17". CIs describe variability across this cohort, not a population of people.
- H1 and H2 are read only by §4.3's four-way rule, with the exact category words: *meaningful harm*, *negligible*, *opposite direction*, *inconclusive* (or *not estimable*). A meaningful-harm reading gives its CI at the Holm-adjusted level. *Negligible* means smaller than the per-event SESOI, not zero. *Inconclusive* says neither harm nor negligibility was shown and gives the CI. *Opposite direction* is an unexpected direction, never "the fault helps". The SESOIs are per-event readings of per-day consensus numbers and are not yet clinically reviewed **[clinical review CR1]**; every H1 and H2 text says so. A well-measured null is a result.
- C1 is estimation only. "Discordant" is used only when §4.4's rule holds; otherwise the shares and their difference are reported without that word, and "concordant" is never claimed from a CI that includes 0. Tau-b over five faults takes few distinct values, so its CI is coarse; the per-fault per-day harms are reported beside it.
- Secondary and exploratory estimates are labelled, carry 95% CIs and are uncorrected; none is called significant. Q4 results are exploratory and are read with the held-out R²; a cell reported as not estimable is not interpreted.
- Interference (#16, #17) and gain errors (#14, #15) are per exposure: there is no population rate, and their harm is never put into per-day or per-year totals. Whenever #17 is reported, the text states that the modern Dexcom G6 barely responds to 1 g acetaminophen (+3.1 mg/dL [Calhoun 2018]). The ±40% gains are stress levels, not measured sizes.
- Q7 is conditional on the rate scenario and on additivity: it is reported as a range across scenarios, each with its θ CI, never as one number; it may be called additive only under §8.5's rule.
- The gate's two optima are compared with the coincidence share, never by eye. Both are in-sample descriptive quantities (§4.11), never recommended settings. The matched filter uses the same compression model that generates the injected artifacts, so its results carry a home advantage ([RESEARCH.md](RESEARCH.md) §6).
- Mechanism and controller-action counts explain effects; they are not separate findings of harm. A null with no threshold crossing (§4.10) is reported as "the fault did not reach the controller's decision threshold", not as "the controller is robust".
- Monte Carlo: every flagged primary estimate (§4.2) is reported as a limitation. Invalid counts (§2.6) are reported with every result.
- Claims follow the exact sentences and credits in [RESEARCH.md](RESEARCH.md) §4 (D15). Nothing is said about the real-world safety of any commercial system, and no result is called "first" unless the novelty gate agrees.

**What may not be said.** That a fault "causes" harm in people; that a result transfers to commercial systems, to Trio or Loop, or to other settings or targets; that a CI covers a population; any significance for a secondary estimate; any per-day or per-year harm for interference or gain errors; that pre-registration makes the study novel.

**Limitations to copy.** Every results text copies the limitations of [RESEARCH.md](RESEARCH.md) §6 in full. It also copies these analysis limitations, *Added here*:
- The alarm layer starts each arm's alarms 120 min before t0, not from the start of the trunk, so detectors whose state builds over hours (the Baysal Kalman estimate) start from that history; the start is identical across the arms of an event.
- Time to re-converge is summarised over the events that re-converged inside the window; events that did not are censored, and their share is reported beside it.
- The SESOIs read per-day consensus numbers as per-event thresholds, which makes them conservative but is not validated **[clinical review CR1]**.
- The standard scenario mix was not run (§8.4).

## 8. Sensitivity analyses

Each is labelled "sensitivity analysis" and reported beside the primary analysis whatever it shows. None changes a primary reading.

### 8.1 S1: the literal reading of the ISF scale step (D17)

- **Data.** The tier-1 arms (#1, #3, #6, #8, #10, #12, #17) of `oref0_tb` at the main seeds 101-110 and all 30 patients, rerun with every patient at its smallest passing scale (19 at 0.8, 7 at 1.0, 3 at 1.25, 1 at 1.5; `results/calibration/selection.json`). Run after the main grid, outside the 24 CPU-h (D18).
- **Analysis.** H1 and H2 exactly as §4.3 (B = 10,000, Holm, four-way rule) and C1's insulin layer (per-day harms, CIs and ranks of the five faults, §4.4). The alarm layer and so tau-b and the shares are not computed, because S1 does not rerun `open_loop`. *Added here:* Kendall's tau-b between S1's and the primary analysis's insulin-layer rankings, with its 95% CI on the shared 10,000-replicate index (the patients are the same), as a summary of whether the ranking moves.
- **How it is run** (SPIKE §14): `python experiments/analyse.py --input results/grid_s1 --out results/analysis_s1 --protocol-tag protocol-v1`. Every analysis runs on what the S1 table contains; those whose arms or controllers are absent are not computed.
- **Requirement on the code:** the tau-b against the primary insulin-layer ranking is not computed at df50a9f. Add it, reading the primary run's replicates or recomputing them on the shared index, with a test.

### 8.2 S2: the 25 patients that meet every acceptance criterion (D17)

- **Data.** The main grid, restricted to the 25 patients not flagged `outside_acceptance` in `results/calibration/isf_scale_map.json` (leaving out adolescent#002, adolescent#007, child#002, child#003 and child#008).
- **Analysis.** H1, H2 and C1 (both layers, tau-b, shares) exactly as §4.3-§4.4, on the 25 patients. The bootstrap stays stratified by age group, with 7 children, 8 adolescents and 10 adults per replicate; the population estimate is the equal-weight mean over the 25 patients.
- **How it is run** (met at df50a9f, SPIKE §14): `python experiments/analyse.py --input results/grid --exclude-outside-acceptance --out results/analysis_s2 --protocol-tag protocol-v1`. The filter reads the five patients from the scale map's `outside_acceptance` list, not from the command line, drops them before every estimate, and records this in `summary.md` and `manifest.json`.

### 8.3 C1's sensitivity analyses (D10)

Computed by `cross_layer.specs()` on the same 10,000-replicate index as C1 and written to the same tables, with `analysis` =:
- `every_night`: compression on every night (rate 1.0), paper 1's headline setting;
- `threshold_alarm`: paper 1's threshold alarm (`threshold_false_alarms`) as the alarm layer;
- `long_gaps_low` and `long_gaps_high`: the G6 long-gap arm (#4, `dropout_long`) at 0.13 and at 0.24 gaps per day [Vettoretti 2023] in place of paper 1's short gaps (#3), whose rate has no source. *Added here* in D10: both ends of the range are run.

Each reports the per-day harms, ranks, tau-b, shares and the discordance reading of §4.4.

### 8.4 The standard scenario mix: not run

D1 named paper 1's standard mix as a sensitivity analysis if compute allowed. D8 allowed it only if the full design at N = 20 fitted with its cost included; it did not (D18). Every results text states that the standard-mix analysis was not run and why.

### 8.5 The additivity check (D8)

**Data** (`experiments/additivity.py`). For each of the 30 patients, two continuous 7-day `oref0_tb` runs on seed 200, `rule1800` with the D17 scale map, paper 1's default mix: one clean, one with paper 1's realistic mixed suite (compression on 86% of nights, the other faults at s = 0.5 with paper 1's rates [Reddy 2026]). The seed fixes meals, sensor noise and the person's draws, so the pair differs only by the faults. Day 1 is warm-up and is not stored. For days d = 2, ..., 7 the script stores per run the minutes below 54, below 70, above 180 and above 250 mg/dL and the risk-weighted minutes (`outcomes.glycaemia`), insulin and rescue, and from the schedule the events starting each day (pisa, dropout, stale, spike) and that day's drift (sign, size, gain at midday). One JSON per patient in `results/additivity/`.

**Analysis**, for o in {`risk_weighted_min`, `minutes_below_54`, `minutes_below_70`, `minutes_above_180`}, with D = 6 days:
- observed: O_p(o) = (1/D) Σ_{d=2}^{7} [Y_faulted(p, d, o) − Y_clean(p, d, o)];
- predicted: Pr_p(o) = (1/D) Σ_f n_f(p) × τ̄_{a(f)}(oref0_tb, p, o), where n_f(p) is the run's count over days 2-7 and τ̄ is the patient effect from the main grid (§3, mean over p's valid seeds 101-110):

| Schedule type f | Arm a(f) | n_f(p) |
|---|---|---|
| pisa | `pisa_typical` (#1) | PISAs starting on days 2-7 |
| dropout | `dropout_short` (#3) | gaps starting on days 2-7 |
| stale | `frozen_typical` (#6) | frozen episodes starting on days 2-7 |
| spike | `spike_typical` (#8) | spike bursts starting on days 2-7 |
| drift | `drift_low` (#10) if the drawn sign is negative, `drift_high` (#12) if positive | 2 × the number of days 2-7 with a drift present (two 12-h windows per day, as C1 counts drift) *Added here* |

- difference: Δ_p(o) = O_p(o) − Pr_p(o).
- Reported: the means over patients of O_p, Pr_p and Δ_p, each with the 95% percentile CI of the stratified patient bootstrap (§4.1, B = 2,000, seed 0); the same replicate index for all three, so Δ's CI accounts for the pairing of observed and predicted within a patient. Units: the outcome's minutes per day.
- A patient whose additivity run is invalid (truth below 10 mg/dL), or who lacks τ̄ for an arm whose events it has, is left out, and the count is reported.
- *Added here:* also reported, the relative departure R(o) = mean Δ(o) / mean O(o), with its 95% percentile CI from the same replicates; R is not estimable when the 95% CI of mean O(o) includes 0, and is then reported as such.
- **Reading.** Q7 on outcome o is called additive only if the 95% CI of mean Δ(o) includes 0 and, *Added here*, the 95% CI of R(o) lies inside (−0.5, +0.5) **[assumption: a departure of half the observed harm would change the reading of Q7; one seed per patient cannot exclude smaller ones]**. If the CI of Δ includes 0 but that of R does not lie inside (−0.5, +0.5), or R is not estimable, the text says the check was too imprecise to detect a departure from additivity of that size, gives both CIs, and does not call Q7 additive; the text then says "the additivity check did not detect a departure from additivity" and gives the CI, never "additivity was shown". If the CI excludes 0, Q7 on o is reported with the observed mean difference and the statement that the per-event effects do not add to the observed harm of a mixed-fault week in this cohort.
- **Limitations of the check.** One seed per patient; the mixed run's drift grows linearly over the sensor's life while arms #10 and #12 hold a gain drawn at a uniform sensor age; faults in the mixed run can overlap and fall in any state, while the main grid's events are single and drawn by slot.

**Code** (met at df50a9f, SPIKE §14): `secondary.additivity` with `secondary.read_additivity`, called by `python experiments/analyse.py --additivity results/additivity`, written to `additivity.csv`; `tests/analysis/test_additivity.py` covers a suite that is exactly additive (Δ = 0) and one that is not. The run definition is in [SIMULATION.md §10.4](SIMULATION.md#104-the-additivity-check-d8-secondary).

**Requirement on the code:** `secondary.additivity` sets `additive` from the CI of Δ alone and does not compute R(o). Add R(o) with its CI and the second condition of the reading, with a test.

## 9. Requirements on the code before the tag

The analysis code must match this file before `protocol-v1` is tagged (D13). Audit against harness commit df50a9f (2026-10-28).

Met at df50a9f: `secondary.additivity` and `analyse.py --additivity` (§8.5, except item 9 below); S2 through `analyse.py --exclude-outside-acceptance` (§8.2); S1 through `analyse.py --input results/grid_s1` (§8.1, except item 5 below).

Still to do:

1. `effects.apply_invalid_rule` groups within a controller (`effects.py:39-52`, grouping by the unit keys); D11 drops the (patient, seed, slot) across all three controllers (§2.6).
2. `outcomes.fault_windows` writes 0 instead of NaN for the after-the-fault columns when the fault lasts the whole window, D = W (`outcomes.py:200-217`) (§2.2).
3. `secondary.correcting_person` is missing, with the correction share (§4.13).
4. `secondary.person_year` does not report each fault's term (§4.15).
5. The S1-against-primary tau-b is not computed (§8.1).
6. `experiments/analyse.py` does not state in `summary.md` that the standard mix was not run (§8.4).
7. The share of events that re-converged within the window is not reported beside `reconverge_min` (§7); add it per controller and arm to `effects_summary.csv` or `mechanism.csv`.
8. The alarm-optimal and harm-optimal rows of `gate_optima.csv` and the gate section of `summary.md` carry the label "in-sample, descriptive" (§4.11).
9. `secondary.additivity` computes R(o) and applies both conditions of the reading (§8.5).

Each change needs a test, and the synthetic-data run of `experiments/analyse.py` must still pass.
