# Simulation specification

The simulator, cohort, controllers, bridge, virtual person, faults, run design, arms, seeds, grid, storage and reproducibility of paper 2, for reviewers and for agents changing the harness. This file is the source of truth for each rule it states; [DECISIONS.md](DECISIONS.md) records why (D1-D19), [RESEARCH.md](RESEARCH.md) holds the questions and [BENCHMARK.md](BENCHMARK.md) the outcomes and statistics. Every rule here is checked by a named test, or is listed in §15 as a requirement on the harness before the freeze.

Conventions. Literature keys are in [CITATIONS.md](../CITATIONS.md) and cited as [Key]. **[assumption]** marks a value we chose, with its reason. **[clinical review CRn]** marks a clinical dependency waiting for a clinician in [CLINICAL_REVIEW.md](../CLINICAL_REVIEW.md); no clinician has reviewed any item, and Cohen Ferrari's literature fact-check is not clinical review. "SPIKE §n" is a section of the harness record [SPIKE_2026-10-28.md](../development/SPIKE_2026-10-28.md). "Our arithmetic" marks a value derived from cited or measured numbers. Code references are to the harness repository (`src/ghostlows_aid/`, `node/`, `experiments/`, `tests/`) at the commit that `protocol-v1` tags. oref0 file:line references are to oref0 0.7.1 as installed in `node/node_modules/oref0`. The design decisions behind this file (D1-D19) were made on 2026-10-28 by an AI agent acting as study lead under authority delegated by the project lead, Pranav Reddy, and are pending his review.

Status on 2026-10-28: the harness and the analysis code are built and tested; the fault-free calibration and the timing are done (D16-D18); no faulted arm has run and `protocol-v1` is not yet tagged.

## 1. Simulator and cohort (D1)

| Item | Value | Source |
|---|---|---|
| Paper 1 | Ghost Lows v0.4.0, imported as the Python package `ghostlows` 0.4.0; version DOI 10.5281/zenodo.23720382, concept DOI 10.5281/zenodo.23573707. Never changed from this repository. | [Reddy 2026] |
| Simulator | simglucose 0.2.11, which implements the 2008 UVA/Padova type 1 diabetes model. | [simglucose; Dalla Man 2007; Kovatchev 2009] |
| Cohort | simglucose's 30 virtual patients: 10 children, 10 adolescents, 10 adults (`ghostlows.sim.engine.PATIENTS`). No patient is ever dropped (D8). | [simglucose; Reddy 2026] |
| Engine | Paper 1's minute-step loop on `T1DPatient` (ghostlows `engine.py:150-185`), re-implemented as `loop.ClosedLoop` with paper 1's pieces imported: patient tables, Quest.csv, scenario, sensor noise, device range, validity floor, `HypoTreatment`, seeds, IOB curve (SPIKE §3.3). | [Reddy 2026] |
| Truth | Plasma glucose, recorded every 1 min. Every outcome is scored on it (AGENTS.md rule 4). | [Reddy 2026] |
| Sensor | Paper 1's 5-min CGM: interstitial glucose plus simglucose's Johnson-distributed noise. | [Breton 2008; Reddy 2026] |
| Device range | 40-400 mg/dL; faulted readings are clipped to it. | [Reddy 2026] |
| Displayed CGM | Rounded to whole mg/dL **[assumption: Dexcom and Nightscout report integer `sgv` values]**; `LoopConfig.cgm_integer=False` restores paper 1's unrounded display (SPIKE §12 item 4). | D3 |
| Validity floor | A run whose simulated glucose falls below 10 mg/dL is outside model validity (`invalid` in the outcome table; rule in BENCHMARK.md). | [Reddy 2026] |
| Scheduled basal | Paper 1's calibrated basal, basal_factor × u2ss × BW / 6000 U/min (`profile.scheduled_basal_u_per_min`), the same for every controller. | [Reddy 2026] (its D17) |
| Clock | Simulated time starts at 2026-01-05 00:00 (paper 1's `DEFAULT_START`, a Monday; calendar dates carry no meaning), read as UTC. | [Reddy 2026] |

Applicability. The 2008 model lacks counter-regulation and intraday variability [Dalla Man 2014; Visentin 2018]. The FDA's acceptance of the UVA/Padova simulator concerns the commercial T1DMS versions, not simglucose [Kovatchev 2009]. Untreated lows in this model last until the insulin wears off, which is why the person treats lows (§5) [Reddy 2026].

### 1.1 Scenario

Paper 1's default mix, unchanged for comparability (`ghostlows` `make_scenario`, mix `DEFAULT`): random meals; carbohydrate miscounting with SD 20%; a missed bolus with probability 0.05 per meal; a skipped meal after its bolus with probability 0.02 per meal; a day-level insulin-sensitivity multiplier drawn from U(0.7, 1.3); 2-4 sessions of unannounced exercise per week **[assumption in paper 1, kept unchanged for comparability]** [Reddy 2026]. A trunk draws its scenario for 3 days (`episodes.SCENARIO_DAYS = 3`). Paper 1's standard mix was to be a sensitivity analysis only if the full design fit at N = 20; it did not, so it does not run and is reported as not run (D18).

Note for implementers: `loop.make_loop` defaults to `mix=STANDARD`; every grid, calibration and additivity call passes `DEFAULT` explicitly (`episodes.run_unit` does).

## 2. Controllers (D3)

| Controller | What runs | Role |
|---|---|---|
| `oref0_tb` | OpenAPS oref0 0.7.1 (npm, MIT) as shipped, temp basals only (`microBolusAllowed` false; SMB and UAM off) [oref0] | Primary |
| `oref1_smb` | The same code with `microBolusAllowed` true and `enableSMB_with_COB` and `enableUAM` on [oref0] | Secondary (Q6) |
| `open_loop` | Paper 1's scheduled basal with the person's meal boluses and no automatic adjustment | Alarm layer (Q2) and the person-only pathway |

All three share the scheduled basal (§1), the person (§5) and the fault hook (§6). oref0 runs unmodified; every harness setting that changes its behaviour is listed below and tested (AGENTS.md rule 6). `loop.CONTROLLERS` holds the three names.

### 2.1 Profile settings

The profile is built by oref0's own `lib/profile` (the bridge's `profile` op, wired as `bin/oref0-get-profile.js:217-244`) from the pump settings and preferences below (`profile.ProfileSettings`, `profile.pump_inputs`, `profile.build_profile`). Every preference not listed keeps its `lib/profile/index.js` default.

| Setting | Value | Marker |
|---|---|---|
| Scheduled basal | Paper 1's basal as U/h (adult#001: 1.26736 U/h). oref0's profile stores it rounded to 0.001 U/h (1.267), while the pump delivers paper 1's exact basal; the difference is at most 0.0005 U/h (D19; SPIKE §12 item 3). | [Reddy 2026]; D19 **[assumption]** |
| ISF | `rule1800`: 1800 / TDI, TDI from Quest.csv, times the patient's ISF scale (§2.3). The harness also offers `quest` (Quest CF) and `rule1700`. | D4, D16, D17 **[assumption; clinical review CR5]** |
| Carb ratio (CR) | Quest CR (Quest.csv). The harness also offers 500 / TDI, not used. | D4 **[assumption; clinical review CR5]** |
| DIA and curve | 5 h, rapid-acting curve (peak 75 min), as paper 1's IOB | **[assumption: comparability with paper 1]** |
| Targets | min_bg = max_bg = 100 mg/dL, so oref0's suspend threshold, min_bg − 0.5 × (min_bg − 40), is 70 mg/dL [OpenAPS docs] | **[assumption: a common OpenAPS default]** |
| max_basal | 4 × the scheduled basal. oref0 then applies min(max_basal, 3 × max daily basal, 4 × current basal); for adult#001 the 3 × cap, 3.801 U/h, binds (SPIKE §3.1). | **[assumption]** |
| max_iob | 65 g / CR + 3 h × the scheduled basal rate (U), the OpenAPS guidance "no more than 3 times your highest basal rate PLUS your typical meal bolus" [OpenAPS docs], with the typical meal bolus read as the mean main-meal bolus; 65 g is the mean of simglucose's three main meals (45, 70 and 80 g) [simglucose]. Set explicitly because the default 0 (`lib/profile/index.js:15`) blocks every high temp. Values per patient in §2.4: 3.774 U (child#009) to 17.923 U (adolescent#008). | **[assumption]** |
| autosens | Off (ratio 1), since autosens would also read the faulted CGM | **[assumption]** |
| skip_neutral_temps | false (avoids the `getMinutes()` branch) | **[assumption]** |
| SMB settings (`oref1_smb` only) | `enableSMB_with_COB` true, `enableUAM` true; the other SMB triggers, maxSMBBasalMinutes 30 min, maxUAMSMBBasalMinutes 30 min, SMBInterval 3 min and bolus_increment 0.1 U at `lib/profile` defaults | **[assumption]** |
| Noise field | Not sent, as from CGM sources that report none, so oref0's noise rules (noise 2 raises the target, noise 3 or more stops looping [oref0 PR 833]) stay inactive. A stated limitation (RESEARCH.md). | D3 |
| Device name | `ghostlows-aid-cgm` (`oref.DEVICE`), never `fakecgm` except in the guard-ablation arm (§8.2), because `determine-basal.js:177-185` skips the flat-CGM check for `fakecgm`. *Departure from D3, recorded here:* D3 asks for "a real CGM name"; the harness sends a project name instead. Behaviour is identical: in oref0 0.7.1 the only branch on the device name is the `fakecgm` test at `lib/determine-basal/determine-basal.js:179`, and `lib/glucose-get-last.js:32` only compares each reading's device with the newest reading's, which a constant name always passes (checked in the installed `node_modules/oref0` on 2026-10-28). A project name also avoids implying that a particular commercial sensor was simulated. | D3 **[assumption: a non-`fakecgm` name is equivalent to a real one in oref0 0.7.1]** |

Tests: `tests/test_mechanisms.py` checks that max_iob is above 0 for all 30 patients and binds when set low (SPIKE §6); `tests/test_bridge_parity.py` and `tests/test_bridge_cli.py` check the profile path (SPIKE §2.2).

### 2.2 Settings rule and its calibration procedure (D4)

Pre-specified before any faulted arm, and run only on calibration seeds without faults (AGENTS.md rule 8). Script: `experiments/calibrate.py run`, then `experiments/calibrate.py table`; one JSON per run in `results/calibration/units/`.

1. Runs: the 30 patients × calibration seeds 90001, 90002 and 90003 × 3 days of the default mix, no fault. `quest` and `rule1800` at scale 1.0 under `oref0_tb` and `oref1_smb`; then `rule1800` under `oref0_tb` at scales 0.8, 1.25 and 1.5 for the scale step.
2. Scoring: on the 1-min plasma truth over all 3 days, the first included, since oref0 starts without history in every arm alike (D19). A patient's three runs are pooled as the mean of the per-run percentages (the runs have equal length) **[assumption: three runs steady the counts at small cost]**.
3. Acceptance criteria, per patient: time in range 70-180 mg/dL ≥ 70%; time below 70 mg/dL ≤ 4%; time below 54 mg/dL ≤ 1% [Battelino 2019]; no loop cycle with a stale-data or temp-mismatch reason; no invalid run **[clinical review CR5]**.
4. Selection: the global ISF rule under which more patients meet every criterion under `oref0_tb` at scale 1.0. Ties go to the rule with the larger median ISF over the 30 patients (less aggressive).
5. Scale step (D4 fallback, as read by D17): if more than 5 of 30 patients fail under the chosen rule, each failing patient takes the smallest scale in {0.8, 1.0, 1.25, 1.5} at which it meets every criterion on the calibration seeds. Passing patients keep 1.0. A failing patient that no scale satisfies keeps 1.0, stays in the cohort and is flagged `outside_acceptance` **[assumption: mimics a clinician leaving settings that meet targets and tuning those that do not; clinical review CR5]**.
6. No per-patient tuning on main seeds, ever.

### 2.3 Calibration result (D16, D17)

630 runs, all valid; no cycle of any run gave a stale-data or temp-mismatch reason; 1,940 s + 1,230 s wall on 4 workers (SPIKE §7).

| ISF rule | Scale | Mode | Patients meeting every criterion | TIR % median [range] | TBR < 70 % | TBR < 54 % | max_iob cycles |
|---|---|---|---|---|---|---|---|
| `quest` | 1 | `oref0_tb` | 13 | 88.5 [36.2, 99.0] | 4.8 [0.0, 13.3] | 0.7 [0.0, 3.7] | 2,131 |
| `quest` | 1 | `oref1_smb` | 13 | 89.2 [35.7, 99.0] | 4.6 [0.0, 11.3] | 0.7 [0.0, 4.4] | 3,526 |
| `rule1800` | 1 | `oref0_tb` | 21 | 89.8 [34.4, 99.0] | 2.3 [0.0, 7.7] | 0.2 [0.0, 2.4] | 0 |
| `rule1800` | 1 | `oref1_smb` | 21 | 90.1 [34.8, 99.1] | 2.4 [0.0, 7.1] | 0.2 [0.0, 2.6] | 252 |
| `rule1800` | D17 map | `oref0_tb` | 25 | 90.6 [34.4, 99.0] | 2.3 [0.0, 7.7] | 0.1 [0.0, 2.4] | 0 |
| `rule1800` | literal map (S1) | `oref0_tb` | 25 | 90.9 [34.4, 99.1] | 2.6 [0.0, 7.7] | 0.1 [0.0, 2.4] | 4 |

Medians and ranges are over the 30 patients' pooled values; "max_iob cycles" counts cycles whose reason names max_iob, over all 90 patient-days (SPIKE §7; D16; D17). The D17 row is from D17; the other rows are SPIKE §7's table.

- **Rule: `rule1800`** (21 of 30 against 13; the same patients under `oref1_smb`). It is also the larger ISF for all 30 patients, so the tie-break would have agreed (`results/calibration/selection.json`: `share_rule1800_isf_larger` = 1.0).
- **Scale step applies:** 9 patients fail under `rule1800` (adolescent#002, adolescent#007, adult#002, adult#009, child#001, child#002, child#003, child#004, child#008).
- **Scale map (D17), `results/calibration/isf_scale_map.json`, the grid script's default:** adult#002, child#001 and child#004 take 1.25; adult#009 takes 1.5; the other 26 keep 1.0. The five failing patients that no scale satisfies (adolescent#002, adolescent#007, child#002, child#003, child#008) keep 1.0 and are flagged `outside_acceptance`. Children 001-004 and adults 002 and 009 fail on time below range (4.1-7.7% below 70 mg/dL with 40-75 g/day of rescue carbohydrate); adolescents 002 and 007, child#003 and child#008 fail on time above range, which an ISF scale cannot fix because meal boluses use the CR (SPIKE §7).
- 25 of 30 then meet every criterion.

### 2.4 Therapy numbers per patient

From `results/calibration/therapy_numbers.csv` (`profile.therapy_numbers` with `rule1800` and the D17 map). ISF is 1800 / TDI × the ISF scale (so adult#002's 38.881 mg/dL/U is 31.105 × 1.25). max_iob is 65 g / CR + 3 h × basal. Per-patient calibration outcomes at scale 1.0 are in SPIKE §7.

| Patient | Basal (U/h) | ISF (mg/dL/U) | CR (g/U) | TDI (U/day) | max_iob (U) | ISF scale | outside_acceptance |
|---|---|---|---|---|---|---|---|
| adolescent#001 | 0.8361 | 49.001 | 12.0 | 36.734 | 7.925 | 1.0 | no |
| adolescent#002 | 0.9180 | 29.018 | 5.0 | 62.031 | 15.754 | 1.0 | **yes** |
| adolescent#003 | 0.6478 | 74.249 | 23.0 | 24.243 | 4.769 | 1.0 | no |
| adolescent#004 | 0.8736 | 51.060 | 14.0 | 35.253 | 7.264 | 1.0 | no |
| adolescent#005 | 0.7224 | 52.934 | 12.0 | 34.005 | 7.584 | 1.0 | no |
| adolescent#006 | 0.8754 | 36.304 | 7.0 | 49.581 | 11.912 | 1.0 | no |
| adolescent#007 | 0.7766 | 41.248 | 8.0 | 43.638 | 10.455 | 1.0 | **yes** |
| adolescent#008 | 0.5578 | 28.397 | 4.0 | 63.387 | 17.923 | 1.0 | no |
| adolescent#009 | 0.6064 | 74.757 | 21.0 | 24.078 | 4.915 | 1.0 | no |
| adolescent#010 | 0.7870 | 54.260 | 14.0 | 33.174 | 7.004 | 1.0 | no |
| adult#001 | 1.2674 | 35.702 | 10.0 | 50.417 | 10.302 | 1.0 | no |
| adult#002 | 1.3695 | 38.881 | 8.0 | 57.869 | 12.234 | **1.25** | no |
| adult#003 | 1.4253 | 31.898 | 9.0 | 56.430 | 11.498 | 1.0 | no |
| adult#004 | 0.8878 | 53.242 | 16.0 | 33.808 | 6.726 | 1.0 | no |
| adult#005 | 1.1798 | 26.348 | 5.0 | 68.316 | 16.539 | 1.0 | no |
| adult#006 | 1.7245 | 29.321 | 10.0 | 61.389 | 11.674 | 1.0 | no |
| adult#007 | 1.3715 | 42.850 | 22.0 | 42.007 | 7.069 | 1.0 | no |
| adult#008 | 1.1414 | 42.077 | 13.0 | 42.779 | 8.424 | 1.0 | no |
| adult#009 | 1.1338 | 40.172 | 5.0 | 67.211 | 16.401 | **1.5** | no |
| adult#010 | 1.0187 | 27.929 | 5.0 | 64.449 | 16.056 | 1.0 | no |
| child#001 | 0.3947 | 128.771 | 25.0 | 17.473 | 3.784 | **1.25** | no |
| child#002 | 0.3951 | 99.022 | 23.0 | 18.178 | 4.011 | 1.0 | **yes** |
| child#003 | 0.2888 | 112.346 | 22.0 | 16.022 | 3.821 | 1.0 | **yes** |
| child#004 | 0.4923 | 113.550 | 25.0 | 19.815 | 4.077 | **1.25** | no |
| child#005 | 0.5151 | 43.973 | 7.0 | 40.935 | 10.831 | 1.0 | no |
| child#006 | 0.4041 | 89.003 | 19.0 | 20.224 | 4.633 | 1.0 | no |
| child#007 | 0.4672 | 49.706 | 8.0 | 36.213 | 9.527 | 1.0 | no |
| child#008 | 0.3400 | 83.743 | 15.0 | 21.494 | 5.353 | 1.0 | **yes** |
| child#009 | 0.3914 | 103.483 | 25.0 | 17.394 | 3.774 | 1.0 | no |
| child#010 | 0.3975 | 87.160 | 18.0 | 20.652 | 4.804 | 1.0 | no |

### 2.5 Settings sensitivity analyses (D17)

- **S1, the literal reading of the scale step:** every patient at its smallest passing scale (19 at 0.8, 7 at 1.0, 3 at 1.25, 1 at 1.5; the five unsatisfiable patients at 1.0), map `results/calibration/selection.json` (its `isf_scale` field). Run for the tier-1 arms of `oref0_tb` at the main seeds and N = 10 (§10.3). It reports H1, H2 and C1's insulin-layer ranking beside the primary analysis whatever they show (BENCHMARK.md).
- **S2:** H1, H2 and C1 on the 25 patients that meet every criterion (the primary grid's rows without the five `outside_acceptance` patients; no extra run).

## 3. The bridge (D3; SPIKE §2)

oref0 is called unmodified through its library by `node/bridge.js`, one long-lived, stateless JSON-lines worker per Python process; the client is `oref.Bridge` and `oref.OrefController.step`. Python sends the whole state every call; oref0 holds no state between calls.

| Item | Rule | Where |
|---|---|---|
| Stages per cycle | `glucose-get-last` (glucose status), `lib/iob`, `lib/meal`, `determine-basal`, each wired as its CLI and given its own fresh copy of the inputs, because oref0 mutates them (`lib/iob/history.js:269-270`, `glucose-get-last.js:48`, `cob.js:143-144`). IOB and meal results pass through a JSON round trip, as through files in the CLI. | `bridge.js` op `step` |
| CLI-only checks kept | Carb ratio below 3 (`bin/oref0-meal.js:88-92`); fewer than 36 readings (`bin/oref0-meal.js:131-135`) | `bridge.js` |
| SMB | The library is called directly, so `microBolusAllowed` is honoured; the CLI turns SMB off without an `autotune` directory (`bin/oref0-determine-basal.js:105-114`) | `bridge.js` |
| Cadence | Once every 5 simulated minutes, also when the reading is missing; oref0 then sees the last reading's age, and its stale-data stop (`determine-basal.js:187`, more than 12 min) acts as on a rig | `loop.ClosedLoop` |
| Glucose input | The readings the controller was shown over the last 24 h (`LoopConfig.glucose_hours = 24`), newest first, each with `date`, `dateString` (needed by `cob.js:46-50`), `sgv` and `device` (`oref.glucose_entry`) | `oref.py`, `loop.py` |
| Treatments input | The pump ledger over the last 24 h (`treatment_hours = 24`), newest first (`history.js:274-280` drops a record newer than the one before it) | `ledger.py` |
| currenttemp and lastTemp | `ledger.currenttemp(minute)` reports the running temp as the pump would; a temp superseded or cancelled early gets its actual duration, so oref0's lastTemp agrees with the pump and the temp-mismatch cancel (`determine-basal.js:369-391`) never fires in clean runs (SPIKE §3.5, §7) | `ledger.py` |
| Meal input | The ledger as pump history and an empty carb history; `meal/history.js:85` reads the carbs of `Meal Bolus` records. Giving the ledger as both would count every meal twice. | SPIKE §3.2 |
| Ledger records | Nightscout treatments: `Temp Basal` (rate U/h, duration min, `created_at`), `Meal Bolus` (insulin, carbs), `Correction Bolus` (insulin; SMBs carry `isSMB`). A temp replaced in its first minute is dropped. | `ledger.py` |
| Pump rounding | Temps as oref0's `lib/round-basal.js` (0.05 U/h below 10 U/h, `profile.model` unset; checked against node on 21 rates); wizard and correction boluses to 0.1 U **[assumption: pump bolus increment, oref0's own default]**; SMBs arrive in 0.1 U steps from determine-basal. determine-basal rounds the basal to 0.05 U/h (`:142`), so a neutral temp for adult#001 runs 1.25 U/h, as on a pump. | `ledger.py`; SPIKE §3.6 |
| Clock | `currentTime` and `inputs.clock` = the simulated clock; `moment.now` is set to the simulated clock per request (moment's documented clock hook), which neutralises the wall-clock read reached through `cob.js:26` via `lib/meal/total.js:22-31`. Every reachable wall-clock read is listed in SPIKE §2.1. | `bridge.js` |
| Timezone | The bridge runs with `TZ=UTC` (set by the Python client), since oref0 uses local-time `getHours()` | `oref.Bridge.start` |
| Logging | stderr output is dropped (returned with `debug: true`); `rT.reason` is kept (`oref.reason_kind` classifies it for the harness's own checks; the analysis counts controller actions with `analysis.outcomes.classify_reasons`, BENCHMARK.md) | `bridge.js`, `oref.py` |
| Install | oref0 0.7.1 from npm through the committed `node/package-lock.json` with install scripts off (`npm ci --ignore-scripts`), not vendored (D19; SPIKE §12 item 1) | `node/` |

Tests (all pass, SPIKE §2.2): `tests/test_bridge_parity.py` (all 52 determine-basal fixtures: oref0's 42 test calls and 10 for branches they never reach, identical `rT`); `tests/test_bridge_cli.py` (the three CLIs on 5 states × 2 ISF rules, identical IOB array, meal data, glucose status and `rT`, SMB off); `tests/test_bridge_clock.py` (identical replies with the wall clock moved to 2019 and to 2034; a negative control differs); oref0's own mocha suite, 84 passing.

## 4. Order within a simulated minute (SPIKE §3.3)

1. Truth is recorded.
2. On a 5-min sample: the CGM is read (raw = interstitial glucose + noise; the fault hook on the raw value; the device display: clip to 40-400 mg/dL, round to whole mg/dL); the person may treat a low (§5); the controller runs.
3. Meals are eaten and bolused.
4. The pump delivers the minute's basal or temp, boluses and SMBs.

Tests (`tests/test_loop.py`, SPIKE §3.5): the open loop with paper 1's bolus rule and display reproduces `ghostlows.sim.engine.simulate` bit for bit (1-min truth, 5-min truth, CGM, noise, boluses, rescues, meals); every controller call gets only readings up to its clock; a fault on sample 50 leaves calls 0-49 identical; with both readers routed to the clean signal a fault leaves truth, insulin and every oref0 treatment and reason identical; two runs are identical; oref0's net IOB equals the IOB of the insulin actually delivered within 0.25 U (max) and 0.08 U (mean).

## 5. The virtual person (D5, D19)

One scripted person for every controller and arm (`loop.Person`), acting on the displayed CGM, which is what people act on **[clinical review CR3]**.

- **Meal bolus (wizard)** at meal start: bolus = announced carbs / CR + max(0, (CGM − target) / ISF) − max(0, IOB), floored at 0 and rounded to 0.1 U, with the profile's CR, ISF and target (100 mg/dL) **[assumption; clinical review CR2]**. IOB is subtracted from the whole bolus, as the plan's technical design writes it (SPIKE §3.4); how real bolus calculators apply IOB is not yet sourced and is a question for CR2. IOB is oref0's latest in closed loop and the ledger's bolus IOB in open loop. The CGM is the newest displayed reading up to 15 min old; with none, the correction term is 0.
- **Meal errors:** paper 1's carbohydrate miscounting, missed boluses and skipped meals (§1.1).
- **Low treatment** (paper 1's rule, its amendment A6, `HypoTreatment`): 15 g of fast carbohydrate when the displayed CGM is below 70 mg/dL by day, or at or below 55 mg/dL at night (23:00-07:00 simulated UTC clock); recheck after 15 min; repeat while below 70 mg/dL [CDC 2024; Garg 2022] **[assumption; clinical review CR2]**. The two triggers are paper 1's threshold alarm (70 mg/dL) and its ungated urgent-low floor (55 mg/dL). The person acts on the newest displayed reading up to 15 min old (`observed`, the default; `truth` exists for tests only). Rescue carbohydrate is eaten but not entered, so oref0 does not count it as COB.
- **Night response probability:** 1.0, paper 1's assumption **[assumption; clinical review CR2]**, drawn from the person stream keyed by sample index, so arms draw the same number for a sample. Lower wake probabilities are derived without new runs (BENCHMARK.md).
- **No corrections outside meals** (paper 1's person).
- **Correcting person** (tier 3 variant only, `episodes.CORRECTING`): a correction bolus of max(0, (CGM − 100) / ISF − max(0, IOB)), rounded to 0.1 U, when the displayed CGM is at or above 250 mg/dL between 07:00 and 23:00 and the person gave no bolus in the last 3 h; SMBs do not count as boluses **[assumption; clinical review CR4]**. 250 mg/dL is the consensus level 2 hyperglycaemia threshold [Battelino 2019]; correction boluses are the harm route named in the 2024 Libre 3 false-high recall [FDA 2024 Libre 3 false-high recall]. IOB is floored at 0 as in the wizard (D19).
- **Person's random numbers:** `SeedSequence([seed, patient index, 7])`, keyed by sample index (SPIKE §3.3).

Tests: `tests/test_mechanisms.py` (correcting person: every correction at a displayed CGM ≥ 250 mg/dL, 07:00-23:00, ≥ 3 h after the last bolus; none with the default person; pathway routing).

## 6. Faults (D1, D6)

### 6.1 How paper 1's operators run online

Faults are injected online: each 5-min reading is corrupted as it is produced, on the raw value (interstitial glucose + noise, before the device range), before the controller or the person reads it, with paper 1's operators and distributions unchanged wherever they exist (`faults_online.py`). Faults never touch true glucose (AGENTS.md rule 4).

- `FaultSchedule.draw(seed, days, config)` draws every instance of paper 1's mixed suite exactly as paper 1's `apply_run` does (same per-type streams, draw order, thinning by severity and label priority), from (seed, planned run length, `FaultConfig`) only. `FaultSchedule.for_run(suite_seed, key, days, config)` is the schedule paper 1's `apply_faults` gives the run with identity `key` = (patient_id, seed, mix), through paper 1's `run_seed`. `apply(k, raw_k)` returns (value or NaN, fault_type, fault_id).
- No operator needs future samples: drift, compression (PISA) and dropout are functions of time; spike noise is pre-drawn per sample from the instance's stream; a frozen segment holds its first reading, a past value. A schedule depends on the planned run length, so it is drawn once for the whole planned run and never extended. It may be switched on late (a branch) except inside a frozen segment, which raises.
- `SingleFault`: one event of a type at onset t0, used by every branched arm. `SingleFault.draw(kind, onset, seed, ...)` takes the shape from paper 1's distributions from its own stream; `SingleFault.at_quantile` sets each size at a quantile; any parameter can be set explicitly; `SingleFault.interference(t0, level)` builds the two interference levels; `describe()` returns the realised parameters, which reach the outcome table (`fault_params_json`).
- Display after the fault: clip to 40-400 mg/dL, round to whole mg/dL. Gains are multiplicative and held from t0 to the end of the window.
- A fault ends at the last reading it changes; compression and interference, whose offsets decay, end when the offset falls to 10% of its peak, as paper 1 defines the realised compression duration (`END_FRACTION`) **[assumption, as paper 1]**.

### 6.2 Online equals post hoc (D1)

On the same schedule, online injection must equal paper 1's post-hoc injection. `tests/test_faults_online.py` (47 tests): each fault type alone, the mixed suite at moderate severity and a dense mixed suite (all types, severity 1, compression every night, 3-day sensors), at seeds 0, 1, 7 and 90001, on a 7-day trace crossing both device limits: `cgm_observed`, `fault_type` and `fault_id` equal `apply_run` sample by sample within 1e-9 (measured: bitwise equal on 55,702 of 55,702 readings) and the instance tables agree; `apply_faults` over two runs is reproduced run by run with `for_run`; changing readings after sample j changes no output up to j; a late start equals the tail of a full run; the loop's own hook on its own clean readings equals `apply_run` post hoc, and with nobody shown the faulted signal truth equals the fault-free run. Single events: nothing changes before onset; PISA equals paper 1's `from_size` attenuation; the shapes and level checks listed in §6.3 (SPIKE §4).

### 6.3 The 17 fault arms

One fault event per arm, starting at the arm's onset t0 (`episodes.FAULTS`; the fault's shape stream is keyed by (seed, patient, 13, fault number), §7.4). Slots are defined in §7.2. Tier 1 always runs; tiers 2 and 3 follow the scaling rule (§10.2).

| # | Arm (`FaultSpec.name`) | Slot | Level | Tier | Generator and parameters | Sources and markers |
|---|---|---|---|---|---|---|
| 1 | `pisa_typical` | night | typical | 1 | `SingleFault.draw("pisa")`: depth and duration from paper 1's two-piece log-normals, median 24 mg/dL and 45 min, 5th-95th percentiles 10-48 mg/dL and 30-70 min; first-order time constant τ = 10 min (`ghostlows.faults.pisa.TAU_MIN`), amplitude A = depth / (1 − e^(−D/τ)), realised duration D + τ ln 10 ([paper 1 SIMULATION.md §4.1 at v0.4.0](https://github.com/pkarakala/ghost-lows/blob/v0.4.0/docs/reference/SIMULATION.md#41-compression-pisa-facchinetti-model), where τ is an **[assumption]**); spikes and drift (#8-#13) follow [paper 1 SIMULATION.md §4.2 at v0.4.0](https://github.com/pkarakala/ghost-lows/blob/v0.4.0/docs/reference/SIMULATION.md#42-other-faults-each-switchable-with-a-severity-knob) | [Facchinetti 2016; Reddy 2026] |
| 2 | `pisa_severe` | night | severe | 2 | PISA with depth 48 mg/dL and duration 70 min (the 95th percentiles; `at_quantile(q = 0.95)` gives the same) | [Facchinetti 2016] **[clinical review CR8]** |
| 3 | `dropout_short` | any | typical | 1 | One gap with paper 1's length model at s = 0.5: exponential, 90% under 20 min; measured median 6.8 min, 93% under 20 min over 300 draws (SPIKE §13) | [Facchinetti 2016; Reddy 2026] |
| 4 | `dropout_long` | any | typical (G6) | 2 | `gap_model="g6"`: 1 + C readings of 5 min; each continuation with probability 0.855 up to the third reading and 0.727 after; median 20 min; 73% last at least 15 min and so exceed oref0's 12-min stale limit (measured over 300 draws: median 20 min, P(≥ 15 min) = 0.72) | [Vettoretti 2023] (Model 4) **[assumption: transition wiring]** |
| 5 | `dropout_severe` | any | severe | 2 | One 120-min gap; readings are missing, so their timestamps are old | [Vettoretti 2023] (longest G6 trial gaps about 95 min in adults, 120 min in children) **[clinical review CR8]** |
| 6 | `frozen_typical` | any | typical | 1 | The last value repeated with fresh timestamps for U(10, 30) min | [Reddy 2026] |
| 7 | `frozen_severe` | any | severe | 2 | Frozen for 120 min, long enough for oref0's flat-CGM check to act while glucose moves (about 45 min then) | [openFDA recall 96684] **[assumption; clinical review CR8]** |
| 8 | `spike_typical` | any | typical | 1 | Paper 1's burst at s = 0.5: U(15, 30) min of extra Gaussian noise, SD 10 mg/dL (20 × s) | [Reddy 2026] |
| 9 | `spike_severe` | any | severe | 2 | 30 min with SD 20 mg/dL (paper 1 at s = 1) | [Reddy 2026] **[clinical review CR8]** |
| 10 | `drift_low` | any | typical | 1 | Paper 1's drift at s = 0.5: end-of-life size U(10%, 20%) × 0.5, reached linearly over a 10-day sensor, at a sensor age drawn uniformly; negative; held for the window as a constant gain | [Reddy 2026] **[assumption: held gain]** |
| 11 | `drift_low_severe` | any | severe | 2 | Held gain −20%: the iCGM 20% agreement bound and paper 1's maximum drift | [21 CFR 862.1355; Reddy 2026] **[clinical review CR8]** |
| 12 | `drift_high` | any | typical | 1 | As #10, positive | [Reddy 2026] |
| 13 | `drift_high_severe` | any | severe | 2 | Held gain +20% | [21 CFR 862.1355] **[clinical review CR8]** |
| 14 | `gain_low` | any | stress | 2 | Held gain −40%: the iCGM 40% bound; the 2025 Libre 3/3 Plus false-low recall published no magnitude | [21 CFR 862.1355; FDA 2025 Libre 3 false-low alert] **[assumption; clinical review CR8]** |
| 15 | `gain_high` | any | stress | 2 | Held gain +40%; the 2024 Libre 3 false-high recall published no magnitude | [21 CFR 862.1355; FDA 2024 Libre 3 false-high recall] **[assumption; clinical review CR8]** |
| 16 | `interference_small` | dose | small | 2 | Additive, peak +20 mg/dL, doses at t0, t0 + 4 h and t0 + 8 h; shape in §6.4 | [FDA K193371] **[assumption: hold and decay; clinical review CR9]** |
| 17 | `interference_large` | dose | large | 1 | Additive, peak +60 mg/dL, one dose; shape in §6.4 | [Maahs 2015; Healio 2015; Matsuyama 2025; MAUDE 2026] **[assumption: decay; clinical review CR9]** |

Arms that share a slot share t0 and the snapshot, so the level contrasts (typical against severe) and dropout against frozen are paired at the event level.

Not modelled in protocol-v1 (D6): a single-dose +10 mg/dL false high, intravenous acetaminophen, a daily hydroxyurea pulse, sustained gain loss of 10-30%, additive offsets of −10 to −30 mg/dL, early-wear offsets, a dying-sensor terminal event, freezes of 60 and 180 min, and paper 1's sweep of the compression fraction. Early-wear error has no established direction and the signal before a sensor fails is not quantified [Staal 2018; FDA K213919]. An amendment can add any of them before they are computed.

### 6.4 The interference model (D6; from research note R4's positive-bias model)

The offset added to the raw reading is A × g(u), where u is the time since t0 and g is a unit-peak curve that rises, holds and decays:

- g(u) = (u / t_p) × exp(1 − u / t_p) for 0 ≤ u ≤ t_p;
- g(u) = 1 for t_p < u ≤ t_p + H;
- g(u) = exp(−(u − t_p − H) / τ) after.

| Level | A | t_p | H | τ | Doses | Check (SPIKE §4) |
|---|---|---|---|---|---|---|
| small (#16) | +20 mg/dL | 2.5 h | 8 h (held to 10.5 h) | 3 h | t0, t0 + 4 h, t0 + 8 h; repeated doses extend the hold | peak, hold and decay |
| large (#17) | +60 mg/dL | 2 h | 0 h | 4 h | one, at t0 | 22% of the peak left at 8 h |

- Additive because the vitamin C effect is independent of background glucose [FDA K193371] and acetaminophen raised the CGM while plasma glucose stayed flat [Basu 2016].
- Small: repeated 1000 mg vitamin C on Libre 2, +9.3 mg/dL 2-3 h after one dose and 18.4 and 19.7 mg/dL after a second and third dose 4 h apart [FDA K193371]; hold and decay **[assumption: no washout data; clinical review CR9]**.
- Large: 1 g acetaminophen on the Dexcom G4, largest mean difference 61 mg/dL at 120 min [Maahs 2015; numbers from a secondary report, Healio 2015]; consistent with intravenous acetaminophen, +55 to +114 mg/dL [Matsuyama 2025], and a hydroxyurea report of about +63 mg/dL [MAUDE 2026]; the size manufacturers warn can cause over-delivery [Medtronic Guardian 4 UG; FDA K203774]. Decay **[assumption; clinical review CR9]**.
- The Dexcom G6 barely responds to 1 g acetaminophen (+3.1 mg/dL [Calhoun 2018]); this is stated wherever the large level is reported. Nothing here is a statement about any commercial system's real-world safety.
- The event ends when the offset falls to 10% of its peak (§6.1). Interference starts at t0 in the dose slot (07:00-15:00 day 2), so every dose of the small level falls in waking hours **[assumption: no source on when people take these drugs]**.

## 7. Run design (D2)

### 7.1 Trunk and warm-up

One continuous fault-free closed-loop run per unit, a unit being (controller c, patient p, seed r): 67 h from 00:00 on day 1 (`episodes.TRUNK_MIN = 4020` min), on the default mix with its scenario drawn for 3 days. The first 24 h (`WARMUP_MIN = 1440`) are warm-up: oref0 gets insulin history over its 5-h DIA and more than the 36 readings `lib/meal` requires (`bin/oref0-meal.js:131-135`). 67 h = the latest night onset (07:00 day 3) + 12 h.

### 7.2 Slots and stratified onsets

| Slot | Minutes from 00:00 day 1 (`episodes.SLOTS`) | Clock | Arms |
|---|---|---|---|
| night | [2820, 3300) | 23:00 day 2 to 07:00 day 3 (paper 1's compression window) | compression, gate arms, second detector |
| dose | [1860, 2340) | 07:00 to 15:00 day 2 | interference, correcting person |
| any | [1440, 2880) | 00:00 to 24:00 day 2 | dropout, frozen, spikes, drift, gains, guard ablation |

`episodes.onset(slot, seed, patient, k, n)`: seed k (0-based) of the N main seeds draws t0 uniformly among the 5-min boundaries in the k-th of N equal sub-intervals of the slot, from the stream `SeedSequence([seed, patient index, 12, slot index])` **[assumption: spreads each patient's events evenly over the slot for Q4 without changing the marginal distribution]**. At N = 10 the sub-intervals are 48 min (night, dose) and 144 min (any) wide (our arithmetic). No controller or arm is in the key, so the three controllers and every arm meet the same onset. The latest possible t0 is minute 3295, and 3295 + 720 = 4015 < 4020, so every window ends inside the trunk (our arithmetic).

### 7.3 Snapshot, restore and window

- `episodes.run_unit` runs the trunk to each slot's t0 in time order and takes `ClosedLoop.snapshot()` there: a deep copy of the whole loop state, namely the body (state vector, time, every `T1DPatient` attribute), the sensor-noise generator, the 1-min arrays, the logged samples and events, the pump ledger and the running temp, the controller's and the person's reading histories, the rescue timer and the person's last bolus, oref0's last IOB, and the fault and detector objects. The scipy dopri5 integrator is rebuilt from (y, t), which continues bit for bit; oref0 holds no state.
- `ClosedLoop.branch(snapshot, faults, routing, person, gate)` restores and sets one arm, which then runs for W = 12 h (`WINDOW_MIN = 720`) **[assumption: covers a compression low (median 45 min [Facchinetti 2016]), the controller's response, the rebound and the next meal; clinical review CR7]**. Outcomes use the half-open window [t0, t0 + 720 min).
- The fault-free arm of a slot (`clean@<slot>`, the table's `clean`) is the trunk's own continuation over the window, at no extra cost. Arms that share a slot share t0 and the snapshot.
- A detector attached to an arm is first fed every reading the controller was shown before t0: it has watched from the start and acts only from t0, so the history before t0 stays the trunk's (SPIKE §5).

Tests (`tests/test_branching.py`, SPIKE §5): (a) for each controller a branched clean arm equals the continuous run bit for bit (truth, 1-min insulin split, every logged sample, events, ledger), restored into the same and into a fresh `ClosedLoop`, including from a snapshot taken mid-meal; (b) restoring twice gives identical arms; (c) dropout, interference with the person routed clean, matched-filter gate and correcting-person arms are bit-identical to the trunk before t0; (d) after t0 arms share the sensor-noise draws and the meal log, and `cgm_clean` is equal until the insulin or rescue carbohydrate the body receives first differs; onsets are stratified on 5-min boundaries; the catalogue has 34/17/17 arms and 7 tier-1 arms per controller. The analysis refuses a slot whose arms' histories before t0 differ (`outcomes.slot_outcomes`).

### 7.4 Common random numbers

After t0 all arms of a unit share every random draw. Which draws can diverge between arms: none. What diverges is only what the shared numbers are used for: whether a rescue is taken, what a wizard dose comes to, what oref0 decides.

| Draw | Source | Key | Same in every arm because |
|---|---|---|---|
| Meals, carbohydrate errors, missed boluses, skipped meals, sensitivity, exercise | paper 1's `make_scenario` | paper 1's meal and perturbation seeds (`run_seeds(seed, patient)`) | drawn before the run |
| Sensor noise | simglucose CGM noise | paper 1's noise seed | one draw per 5-min sample, whether the reading is faulted, missing, gated or shown; state in the snapshot |
| Night response | the person | (person seed, 7, sample index) | keyed by sample, not by draw count |
| Onsets | `episodes.onset` | (seed, patient, 12, slot) | no controller or arm in the key |
| Fault shape | `SingleFault.draw` via `episodes.event_seed` | (seed, patient, 13, fault number) | shared by controllers and variant arms |
| Spike noise | `SingleFault` | (event seed, sample index) | keyed by sample |

Across controllers the trunks differ (each controller doses differently), so a contrast between controllers compares two paired effects on the same patient, seed, scenario, noise, onset and fault shape (BENCHMARK.md).

## 8. Arms per controller (D6, D7)

`episodes.ARMS`, filtered by `episodes.arms_for(controller)`. Variant names are `<base>+<variant>`; a base of `clean` means the slot's fault-free continuation.

### 8.1 Fault arms (all three controllers)

The 17 arms of §6.3: 7 at tier 1 (#1, #3, #6, #8, #10, #12, #17) and 10 at tier 2.

### 8.2 Variant arms (`oref0_tb` only; 17)

| Group (`Arm.group`) | Arms | Slot | Tier | Rule | Markers |
|---|---|---|---|---|---|
| pathway (Q3) | `pisa_typical+controller_only`, `pisa_typical+person_only`, `interference_large+controller_only`, `interference_large+person_only` | night, dose | 2 | Controller-only: the controller sees the faulted CGM, the person the fault-free CGM (`Routing(person_sees="clean")`); person-only: the reverse. With the default arm (both) and the fault-free arm, all four cells of the 2 × 2. | D7 |
| gate (Q5) | `clean+baysal_<set>`, `pisa_typical+baysal_<set>` for the sets aggressive (g′_in −1.9, g′_out −2.8 mg/dL/min, ratio 1.2, tolerance 0.1), trial (−2.0, −3.1, 1.5, 0.3; paper 1's pre-specified set), nominal (−2.5, −2.8, 1.3, 0.3), cautious (−3.1, −2.8, 1.3, 0.3): 8 arms | night | 2 | Paper 1's causal Baysal detector (`detectors.BaysalOnline`) between the displayed CGM and the controller; policy hold (`Routing(gate=True, gate_policy="hold")`): a flagged reading is replaced by the last unflagged reading, sent with a fresh timestamp so that oref0's flat-CGM check can still act. The person sees the raw displayed CGM. | [Baysal 2014; Reddy 2026] **[assumption; clinical review CR10]** |
| detector2 (Q5) | `clean+matched_filter`, `pisa_typical+matched_filter` | night | 3 | Paper 1's causal matched filter on the last five readings (`detectors.MatchedFilterOnline`), a causal adaptation of [Idi 2023], flag above 0.5, same policy | [Idi 2023; Reddy 2026] |
| guard (Q5) | `frozen_severe+fakecgm` | any | 2 | Flat-CGM check ablated through the input: every glucose entry sent to oref0 from t0 on, history included, carries device `fakecgm` (`Routing(fakecgm=True)`), because `glucose-get-last` keeps only the newest reading's device (`lib/glucose-get-last.js:32`). The stale-data stop is ablated without an extra arm: `dropout_severe` (old timestamps, the stop acts after 12 min) against `frozen_severe` (fresh timestamps) on the same t0. | D7 |
| correcting (N7) | `clean+correcting`, `interference_large+correcting` | dose | 3 | The correcting person of §5 | **[clinical review CR4]** |

Each detector equals paper 1's batch flags bit for bit for all four Baysal sets and the matched filter, and a warmed detector equals one that watched from the start (`tests/test_detectors.py`, SPIKE §5). The routing hook `restamp_missing` (re-sends the last reading stamped "now" when one is missing) exists for the harness-validity suite only and is in no arm.

### 8.3 Counts per (patient, seed)

| Controller | Tier 1 | Tier 2 | Tier 3 | Branched arms | Clean arms (free) |
|---|---|---|---|---|---|
| `oref0_tb` | 7 | 23 (10 fault, 4 pathway, 8 gate, 1 guard) | 4 (2 matched filter, 2 correcting) | 34 | 3 |
| `oref1_smb` | 7 | 10 | 0 | 17 | 3 |
| `open_loop` | 7 | 10 | 0 | 17 | 3 |

68 branched arms of 12 h and 3 trunks of 67 h per (patient, seed): 1,017 simulated hours (D8; our arithmetic).

## 9. Seeds (D4, D8)

| Use | Seeds | Rule |
|---|---|---|
| Calibration (settings, timing, harness tests) | 90001, 90002, 90003 (all ≥ 90000) | No fault in calibration; timing ran with faults off |
| Main grid | 101-110 (N = 10, D18); seed 100 + k + 1 has onset stratum k = 0, ..., 9 | `run_grid.MAIN_SEED0 = 101` |
| Additivity check | 200 | **[assumption: outside the main and calibration seeds; fixed so the check is reproducible]** |
| Trace sample | 20261028 | §11 |
| Paper 1 (not used here) | benchmark seeds 1-5, calibration seed 9001 | disjoint from all of the above |

A seed sets paper 1's meal, perturbation and noise seeds for the patient (`task_seeds`), so meals and sensor noise equal paper 1's for the same (patient, seed).

## 10. Grid, cost and the scaling rule (D8, D18)

### 10.1 Measured cost (SPIKE §8; `results/timing/cost.json`)

Measured on 4 CPUs (15 GB RAM) with `experiments/run_grid.py` in timing mode: every arm with its fault off, variants kept, `rule1800` at scale 1, calibration seed 90001, onset stratum k = 0 of N = 20; 6 patients (adolescent, adult and child #001 and #006) × 3 controllers on 4 workers.

| Controller | 24-h warm-up (s) | Trunk, 67 h (s) | s per patient-day | 12-h arm (s) | Arms per unit | Unit (s) |
|---|---|---|---|---|---|---|
| `oref0_tb` | 5.36 | 20.0 | 7.16 | 3.65 | 34 | 147.2 |
| `oref1_smb` | 5.78 | 22.2 | 7.94 | 4.12 | 17 | 92.2 |
| `open_loop` | 0.64 | 2.13 | 0.76 | 0.41 | 17 | 9.4 |

Bridge call: 19.0 ms median (32.8 ms p95) for `oref0_tb` and 21.7 ms (38.6 ms) for `oref1_smb` on 4 workers; IOB about 10 ms and meal about 8 ms of it. 4-worker speedup 3.50× (8 trunks: 153.7 s on 1 worker, 43.9 s on 4). The ODE alone runs at 0.56 s per simulated day. Not in the projection: fault-specific costs (fewer readings in a dropout, a stale-data exit) and the analysis step, 0.03-0.04 s per arm, about 0.2 CPU-h at N = 10 (SPIKE §8, §12 item 10).

### 10.2 Scaling rule and its result

Rule, fixed before the freeze: a budget of 24 CPU-hours, projected from the measured costs. N = 20 if the full design fits; else N = 15; else N = 10. If N = 10 does not fit, drop, in this order, until it fits: (1) tier 3; (2) the tier-2 arms of `oref1_smb`; (3) the tier-2 arms of `open_loop`. Tier 1 always runs at N ≥ 10. Patients are never reduced. If the design left after the three drops still does not fit at N = 10, it runs at N = 10 anyway and the overrun is recorded.

| N | Arms | CPU-h | Wall h on 4 workers | Fits 24 CPU-h |
|---|---|---|---|---|
| 20 | all tiers | 41.4 | 10.4 | no |
| 15 | all tiers | 31.1 | 7.8 | no |
| 10 | all tiers | 20.7 | 5.2 | **yes** |
| 10 | tier 3 dropped | 19.4 | 4.8 | yes |
| 10 | and `oref1_smb` tier 2 | 15.9 | 4.0 | yes |
| 10 | and `open_loop` tier 2 | 15.6 | 3.9 | yes |

**Result (D18): N = 10 with every tier.** 30 patients × seeds 101-110 × 3 controllers = 900 units, 20.7 CPU-h (14% under the budget), about 5.2 h on 4 workers. The grid simulates 30 × 10 × 1,017 = 305,100 h (our arithmetic). `oref0_tb` is 59% of the cost. The standard-mix sensitivity analysis does not run (§1.1).

### 10.3 Commands

All commands run from the harness repository root after `cd node && npm ci --ignore-scripts`, and only after `protocol-v1` is tagged and recorded in DECISIONS.md (D13).

Main grid (defaults: `--isf-rule rule1800`, `--isf-scale-map results/calibration/isf_scale_map.json`, `--max-tier 3`, all three controllers, all 30 patients, `--workers 4`, `--out results/grid`; `--n` defaults to 20, so `--n 10` is required):

```
uv run python experiments/run_grid.py run --n 10 --faults --frozen docs/reference/DECISIONS.md
uv run python experiments/run_grid.py cost --out results/grid
```

S1 (D17; tier-1 arms of `oref0_tb`, the literal scale map, its own output folder; about 3.8 CPU-h, D18):

```
uv run python experiments/run_grid.py run --n 10 --faults --frozen docs/reference/DECISIONS.md \
    --isf-scale-map results/calibration/selection.json --max-tier 1 --controllers oref0_tb \
    --out results/grid_s1
```

S1 uses the same seeds and onset strata as the main grid, so its units pair with the main grid's `oref0_tb` units patient by patient and seed by seed. Its 19 patients whose scale differs from the D17 map (0.8 against 1.0) have different trunks; the 11 whose scale is the same (the 4 rescaled patients, the 5 unsatisfiable patients, and adolescent#005 and adult#006 at 1.0) must reproduce the main grid's `oref0_tb` clean and tier-1 rows bit for bit, which is a check, not a result.

The additivity check (§10.4) and S1 run after the main grid and outside the 24 CPU-h, which serves only to choose N (D18).

### 10.4 The additivity check (D8; secondary)

Question: does the per-day harm of paper 1's realistic mixed suite in a continuous run equal the sum of per-event effects × event counts? The per-person-year estimate (Q7) and C1 assume it.

This section defines the runs; what is stored, how they are analysed and how the result is read are specified once, in [BENCHMARK.md §8.5](BENCHMARK.md#85-the-additivity-check-d8). `experiments/additivity.py` exists and is tested (SPIKE §14); it has not run.

1. Units: the 30 patients × suite and loop seed 200 × `oref0_tb`, `rule1800` with the D17 scale map, the default person, mix `DEFAULT` passed explicitly to `loop.make_loop(config, 200, 7, mix=DEFAULT, faults=...)`. Seed 200 **[assumption: outside the main seeds 101-110 and the calibration seeds 90001-90003, so the check reuses no draw of the grid or of the calibration; D18]**.
2. Two continuous 7-day runs per patient (10,080 min each, scenario drawn for 7 days), no snapshot or branch: (a) clean, `faults=None`; (b) faulted, with `faults = FaultSchedule.for_run(200, ("sim:<patient>", 200, "default"), 7, FaultConfig(pisa_fraction=0.86, severity=0.5))`, where `FaultConfig` is paper 1's (`ghostlows.faults.suite.FaultConfig`) with its default rates and sizes, thinned and scaled by s = 0.5 [Reddy 2026]. The schedule is drawn once for the planned 7 days.
3. Both runs share the seed, so they share the meal scenario, the meal errors, the sensor noise and the person's draws (§7.4); only the fault schedule differs.
4. Day 1 is warm-up, as in the trunks (§7.1), and is not stored; days 2-7 are stored, one JSON file per patient in `results/additivity/`, written atomically, so a restart skips finished patients (§11).
5. A test checks, on a shortened run, that the faulted run's `cgm_observed` equals paper 1's `apply_faults` applied post hoc to the run's own clean readings, and that the clean run equals the faulted run before the first fault sample (§15 item 1).

Cost: 30 patients × 2 runs × 7 days at 7.16 s per patient-day is about 0.8 CPU-h (our arithmetic; D18 estimates about 0.7 CPU-h).

## 11. What is stored, and how a run resumes (D8, D9, D19; SPIKE §10)

Per unit (controller, patient, seed), in the output folder:

| File | Content | When |
|---|---|---|
| `outcomes/<unit>.parquet` | One row per arm and slot, the table of `analysis.schema` (validated by `schema.validate`): keys (controller, patient_id, seed, slot, arm, data_source); design (fault kind, level, variant); realised fault parameters (onset minute, fault end, depth, duration, size, time constant, area, `fault_params_json`); validity; glycaemia on truth over [t0, t0 + W); insulin and carbohydrate; during- and after-fault quantities with the clean arm's twins; re-convergence; controller actions (cycles whose `rT.reason` shows each of stale, flat, temp_mismatch, lgs_zero_temp, smb_maxdelta; `analysis.outcomes.classify_reasons`); threshold crossings of the controller's input (below 70, above 180 mg/dL during the fault); the onset state (true glucose and 15-min trend, fault-free displayed CGM and trend, IOB, COB, hours since the last bolus, hour of day, all from the slot's clean arm); paper 1's alarm layer on the arm's displayed CGM. Definitions in BENCHMARK.md. | `--faults` only |
| `traces/<unit>.pkl.gz` | `episodes.arm_trace` of every arm, cut to [t0, t0 + W]: 1-min truth; 1-min insulin split into basal or temp, the person's boluses and SMBs; the 5-min displayed CGM, the person's and the controller's input, what was sent, the gate flag, oref0's outputs and reason; events; IOB and COB at onset; the realised fault | for the 10% sample only |
| `units/<unit>.json` | Cost record (warm-up, trunk and per-arm seconds, bridge stage times), onsets, ISF rule and scale, `faults`, `max_tier`, k and n | written last; its presence marks the unit finished |
| `run_log.jsonl` | One line per invocation: jobs, units run, workers, wall time | appended |

- Outcomes are computed by the analysis package (`analysis.outcomes.trace_from_run`, `slot_outcomes`), not by the harness; `run_grid.run_job` only passes the traces on.
- **The 10% trace sample** (D19) **[assumption: 10% gives about 30 example pairs for figures and audits at a small fraction of the storage of keeping every trace]**: (patient, seed) is kept when `numpy.random.default_rng([20261028, patient index, seed]).random() < 0.10` (`run_grid.keeps_traces`), for every controller. Fixed before any run and independent of N and of the patient list. At N = 10 about 30 of the 300 pairs are expected (our arithmetic); the drawn list is recorded with the grid's results.
- **Resume:** every file is written to a temporary name and renamed (`run_grid.write_atomic`), and the unit marker is written last, so an interrupted run leaves no partial file under a final name. A restarted command skips every unit whose marker exists, so a container restart loses at most the units in progress (4 with 4 workers) **[engineering requirement]**. Re-running a unit gives identical outputs (determinism, §12).
- Invalid arms are flagged in the table; the rule that drops a whole (patient, seed, slot) for every arm of every controller is in BENCHMARK.md.

## 12. Reproducibility and caching

- Pins: oref0 0.7.1 from npm via `node/package-lock.json` (moment 2.31.0, moment-timezone 0.5.x as resolved there), node v22.22.0, simglucose 0.2.11, `ghostlows` 0.4.0, and the Python environment in `uv.lock` (SPIKE header).
- Determinism: two runs give identical logs, timing columns aside (`tests/test_loop.py`); restoring twice gives identical arms (`tests/test_branching.py`); replies do not depend on the wall clock (`tests/test_bridge_clock.py`). Every random draw is keyed as in §7.4.
- Speed: patient parameters are held as a `SimpleNamespace` instead of a pandas Series (`LoopConfig.fast_params`), bitwise identical and 14× faster (SPIKE §3.3). The §2.1 levers of the technical design (trimmed histories, skipping `lib/meal`) were checked and not adopted, because neither gives identical outputs in every state of both modes (SPIKE §9); `OrefController.skip_meal` stays off.
- Caching: calibration runs (`results/calibration/units/`), grid units (`results/grid/units/`) and S1 units (`results/grid_s1/units/`) are cached by their markers. Outputs are gitignored; no data go into git (AGENTS.md rule 9).
- Tests: `uv run pytest -q` (209 passed, 6 deselected at the calibration commit: 125 harness, 84 analysis) and `uv run pytest -q -m slow` (6 passed); `ruff check` and `ruff format --check` pass (SPIKE §11).

## 13. Harness-validity suite (D12; `tests/test_mechanisms.py`, SPIKE §6)

Built and passing (adult#001, calibration seed 90001, default mix, `rule1800`; each check asserts on controller inputs, reasons and actions only, never on a faulted-versus-clean glucose difference):

- the stale-data stop acts on every cycle from 10 min after the last reading to the end of a 120-min gap (22 cycles), none outside, with no high temp while blind; the `restamp_missing` ablation removes it;
- the flat-CGM check acts on every cycle from 50 min into a 120-min freeze with the real device name, and never with `fakecgm`;
- a predicted low (PISA 48 mg/dL, 70 min, at 02:00) sets or keeps a zero temp;
- max_iob is above 0 for all 30 patients, and binds when set to 1 U under a +40% gain;
- SMB is blocked by maxDelta on a +40 mg/dL jump (`oref1_smb`);
- the gate's hold and drop policies, a gate attached but off, the pathway routing, the correcting person and the open loop behave as specified.

In clean runs the flat-CGM branch fires on 3.6-5.6% of calls, because whole-mg/dL readings and smooth noise give delta 0; the first call of every run is one of them (SPIKE §3.6).

## 14. What protocol-v1 does not include

Trio's oref, Loop and AndroidAPS; unannounced meals; a target sweep; autosens on; a noise field; per-patient tuning on main seeds; paper 1's standard mix; harm measured on real data (D3, D14). Nothing here or in the results is a statement about commercial systems.

## 15. Requirements on the harness before the freeze

Each must be met, with its test, on the commit that `protocol-v1` tags. None needs a faulted-versus-clean outcome.

1. `experiments/additivity.py` meets §10.4 items 1-4 (met at harness commit df50a9f, SPIKE §14) and has the test of §10.4 item 5 (not yet written: `tests/analysis/test_additivity.py` checks the schedule's bookkeeping and the analysis only). **Met at 2b981a3:** `tests/test_additivity_run.py` (seed 12345, one day: the displayed CGM equals `apply_faults` post hoc, and the clean run equals the faulted run before the first changed reading).
2. `run_grid.py run --faults` refuses to start unless the `--frozen` file records the `protocol-v1` tag (today it checks only that the file exists), and each unit marker records the harness commit, the `ghostlows`, simglucose and oref0 versions, and the path and SHA-256 of the ISF scale map used. **Met at 2b981a3:** `run_grid.freeze_ok` and `provenance`; `tests/test_grid.py`.
3. On a restart, `run_grid.py` stops with an error, rather than skipping, when an existing unit marker's `faults`, `isf_rule`, `isf_scale`, `max_tier`, `k` or `n` differ from the command's, so a timing or S1 marker can never stand in for a main-grid unit. **Met at 2b981a3:** `run_grid.check_markers`; `tests/test_grid.py`.
4. The list of (patient, seed) pairs in the 10% trace sample at N = 10 is written to the grid's output folder before the first unit runs. **Met at 2b981a3:** `trace_sample.csv` in the output folder; `tests/test_grid.py`.
5. `run_grid.py run` stops with an error, before any unit runs, when a patient in the command's patient list is missing from the ISF scale map it reads. Today `run_grid.py:102-105` uses `scales.get(p, 1.0)`, which silently gives such a patient scale 1.0. **Met at 2b981a3** (also in `additivity.py`); `tests/test_grid.py`.
6. **Source checks still open.** Five frozen values rest on a source whose check is still open in [CITATIONS.md](../CITATIONS.md) or [CLINICAL_REVIEW.md](../CLINICAL_REVIEW.md): the transition wiring of arm #4's G6 gap model [Vettoretti 2023] (§6.3); the +60 mg/dL size of arm #17, read through a secondary report [Maahs 2015] (§6.3, §6.4); the high compression rate of the Q7 scenarios [Idi 2025] (BENCHMARK.md §4.15); the wording of the OpenAPS max_iob guidance [OpenAPS docs] (§2.1); and the GRI weights and bands [Klonoff 2023] (BENCHMARK.md §2.3). Each check is done, and its outcome recorded in [CITATIONS.md](../CITATIONS.md), before `protocol-v1` is tagged; any correction it calls for is made before the tag. A check not done by the tag freezes the value as written, and the unchecked source is reported as a limitation in every results text. No value is changed after the tag on the strength of such a check except by a dated amendment made before the affected result is computed (AGENTS.md rule 11).
