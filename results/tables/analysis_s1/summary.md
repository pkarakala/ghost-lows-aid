# Paper 2 analysis summary

## Data

- rows 3000, units (controller, patient, seed, slot) 900
- controllers oref0_tb
- patients 30, seeds per patient 10 to 10, arms 8
- data source simulation
- invalid rule: 2 rows dropped (1 invalid arms)

The standard scenario mix was not run: D8 allowed it only if the full design at N = 20 fitted the budget with its cost included, and it did not (D18; BENCHMARK §8.4).

## Confirmatory family H1-H2 (oref0_tb; Holm, memo §9)

| id | arm | outcome | sesoi | theta | ci_level | ci_lo | ci_hi | ci90_lo | ci90_hi | p_value | rejected | decision | mc_flag |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| H1 | pisa_typical | minutes_above_180 | 72 | 4.94 | 0.975 | 1.87 | 8.61 | 2.59 | 7.6 | 0.0002 | True | negligible | True |
| H2 | interference_large | minutes_below_70 | 15 | 57.3 | 0.95 | 48.7 | 66 | 50.1 | 64.6 | 0.0002 | True | meaningful harm | True |

Co-reported outcomes (95% CI, uncorrected): hypotheses_co_reported.csv.

## C1, the pre-specified primary comparison (estimation only)

| analysis | tau_b | tau_b_lo | tau_b_hi | share_alarm | share_insulin | share_diff | share_diff_lo | share_diff_hi | discordant |
|---|---|---|---|---|---|---|---|---|---|
| main |  |  |  |  |  |  |  |  | <NA> |
| every_night |  |  |  |  |  |  |  |  | <NA> |
| threshold_alarm |  |  |  |  |  |  |  |  | <NA> |
| long_gaps_low |  |  |  |  |  |  |  |  | <NA> |
| long_gaps_high |  |  |  |  |  |  |  |  | <NA> |

Per-fault per-day harm: c1_per_fault.csv.

## Secondary and exploratory (labelled, uncorrected, 95% CI)

Gate optima (in-sample, descriptive; neither is a recommended setting):

| harm_optimal | alarm_optimal | alarm_rule_feasible | coincide_share |
|---|---|---|---|
|  |  | False |  |

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
- s1_tau_b.csv

## Run

```json
{
 "commit": "69c6c1364afc1033b0d0edd9cb84fb99e678a010",
 "argv": [
  "--input",
  "results/grid_s1/outcomes",
  "--primary-input",
  "results/grid/outcomes",
  "--protocol-tag",
  "protocol-v1",
  "--no-figures",
  "--out",
  "results/analysis_s1"
 ],
 "inputs": [
  "results/grid_s1/outcomes"
 ],
 "data_source": [
  "simulation"
 ],
 "n_boot": 2000,
 "n_boot_confirmatory": 10000,
 "boot_seed": 0,
 "protocol_tag": "protocol-v1",
 "excluded_outside_acceptance": [],
 "additivity": null,
 "seconds": 13.9
}
```
