# Paper 2 analysis summary

## Data

- rows 23100, units (controller, patient, seed, slot) 2700
- controllers open_loop, oref0_tb, oref1_smb
- patients 30, seeds per patient 10 to 10, arms 35
- data source simulation
- invalid rule: 99 rows dropped (11 invalid arms)

The standard scenario mix was not run: D8 allowed it only if the full design at N = 20 fitted the budget with its cost included, and it did not (D18; BENCHMARK §8.4).

## Confirmatory family H1-H2 (oref0_tb; Holm, memo §9)

| id | arm | outcome | sesoi | theta | ci_level | ci_lo | ci_hi | ci90_lo | ci90_hi | p_value | rejected | decision | mc_flag |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| H1 | pisa_typical | minutes_above_180 | 72 | 5.43 | 0.975 | 2.4 | 9.03 | 3.05 | 8.01 | 0.0002 | True | negligible | True |
| H2 | interference_large | minutes_below_70 | 15 | 50.1 | 0.95 | 41.8 | 58.6 | 43.1 | 57.3 | 0.0002 | True | meaningful harm | True |

Co-reported outcomes (95% CI, uncorrected): hypotheses_co_reported.csv.

## C1, the pre-specified primary comparison (estimation only)

| analysis | tau_b | tau_b_lo | tau_b_hi | share_alarm | share_insulin | share_diff | share_diff_lo | share_diff_hi | discordant |
|---|---|---|---|---|---|---|---|---|---|
| main | 0 | -0.6 | 0.4 | 0.706 |  |  |  |  | <NA> |
| every_night | 0 | -0.6 | 0.4 | 0.736 |  |  |  |  | <NA> |
| threshold_alarm | 0 | -0.6 | 0.6 | 0.479 |  |  |  |  | <NA> |
| long_gaps_low | 0 | -0.8 | 0.4 | 0.65 |  |  |  |  | <NA> |
| long_gaps_high | 0 | -0.8 | 0.4 | 0.652 |  |  |  |  | <NA> |

Per-fault per-day harm: c1_per_fault.csv.

## Secondary and exploratory (labelled, uncorrected, 95% CI)

Gate optima (in-sample, descriptive; neither is a recommended setting):

| harm_optimal | alarm_optimal | alarm_rule_feasible | coincide_share |
|---|---|---|---|
| no_gate | trial | True | 0.003 |

## D8 additivity (secondary, 95% CI)

| outcome | quantity | theta | ci_lo | ci_hi | n_patients | additive |
|---|---|---|---|---|---|---|
| risk_weighted_min | observed | -1.49 | -13.6 | 9.94 | 30 | <NA> |
| minutes_below_54 | observed | -0.494 | -1.38 | 0.317 | 30 | <NA> |
| minutes_below_70 | observed | -3.98 | -8.03 | -0.422 | 30 | <NA> |
| minutes_above_180 | observed | 9.16 | 1.3 | 18 | 30 | <NA> |
| risk_weighted_min | predicted | 1.84 | -12.7 | 15.4 | 30 | <NA> |
| minutes_below_54 | predicted | -0.322 | -2.21 | 1.13 | 30 | <NA> |
| minutes_below_70 | predicted | -4.61 | -10.7 | 1.5 | 30 | <NA> |
| minutes_above_180 | predicted | 14.7 | 8.26 | 21.4 | 30 | <NA> |
| risk_weighted_min | difference | -3.33 | -22.1 | 17.6 | 30 | False |
| minutes_below_54 | difference | -0.172 | -1.31 | 1.15 | 30 | False |
| minutes_below_70 | difference | 0.629 | -6.42 | 7.32 | 30 | False |
| minutes_above_180 | difference | -5.56 | -13.4 | 3.48 | 30 | False |
| risk_weighted_min | relative_departure |  |  |  | 30 | <NA> |
| minutes_below_54 | relative_departure |  |  |  | 30 | <NA> |
| minutes_below_70 | relative_departure | -0.158 | -4.91 | 2.06 | 30 | <NA> |
| minutes_above_180 | relative_departure | -0.607 | -5.64 | 0.32 | 30 | <NA> |

## Files

- invalid_counts.csv
- clean_levels.csv
- effects_summary.csv
- hypotheses.csv
- hypotheses_co_reported.csv
- c1_per_fault.csv
- c1_summary.csv
- pairing_efficiency.csv
- dose_response.csv
- dose_by_crossing.csv
- ladders.csv
- state_cells.csv
- state_model.csv
- state_held_out.csv
- pathway.csv
- wake_curve.csv
- mechanism.csv
- gate_curves.csv
- gate_optima.csv
- correcting_person.csv
- guards.csv
- controllers.csv
- smb_tb_split.csv
- person_year.csv
- additivity.csv

## Run

```json
{
 "commit": "69c6c1364afc1033b0d0edd9cb84fb99e678a010",
 "argv": [
  "--input",
  "results/grid/outcomes",
  "--protocol-tag",
  "protocol-v1",
  "--additivity",
  "results/additivity",
  "--no-figures",
  "--out",
  "results/analysis_additivity"
 ],
 "inputs": [
  "results/grid/outcomes"
 ],
 "data_source": [
  "simulation"
 ],
 "n_boot": 2000,
 "n_boot_confirmatory": 10000,
 "boot_seed": 0,
 "protocol_tag": "protocol-v1",
 "excluded_outside_acceptance": [],
 "additivity": "results/additivity",
 "seconds": 19.6
}
```
