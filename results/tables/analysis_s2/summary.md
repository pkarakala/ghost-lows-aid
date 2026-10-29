# Paper 2 analysis summary

## Data

- rows 19250, units (controller, patient, seed, slot) 2250
- controllers open_loop, oref0_tb, oref1_smb
- patients 25, seeds per patient 10 to 10, arms 35
- data source simulation
- invalid rule: 13 rows dropped (2 invalid arms)

The standard scenario mix was not run: D8 allowed it only if the full design at N = 20 fitted the budget with its cost included, and it did not (D18; BENCHMARK §8.4).

**S2: 5 patients outside acceptance excluded before every estimate** (sim:adolescent#002, sim:adolescent#007, sim:child#002, sim:child#003, sim:child#008; /home/user/ghost-lows-aid/results/calibration/isf_scale_map.json).

## Confirmatory family H1-H2 (oref0_tb; Holm, memo §9)

| id | arm | outcome | sesoi | theta | ci_level | ci_lo | ci_hi | ci90_lo | ci90_hi | p_value | rejected | decision | mc_flag |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| H1 | pisa_typical | minutes_above_180 | 72 | 4.15 | 0.975 | 1.52 | 7.24 | 2.15 | 6.34 | 0.0002 | True | negligible | True |
| H2 | interference_large | minutes_below_70 | 15 | 51.9 | 0.95 | 42.7 | 61.3 | 44.1 | 59.8 | 0.0002 | True | meaningful harm | True |

Co-reported outcomes (95% CI, uncorrected): hypotheses_co_reported.csv.

## C1, the pre-specified primary comparison (estimation only)

| analysis | tau_b | tau_b_lo | tau_b_hi | share_alarm | share_insulin | share_diff | share_diff_lo | share_diff_hi | discordant |
|---|---|---|---|---|---|---|---|---|---|
| main | 0.2 | -0.6 | 0.6 | 0.718 |  |  |  |  | <NA> |
| every_night | 0.2 | -0.6 | 0.6 | 0.747 | 0.0589 | 0.688 | 0.208 | 1.55 | True |
| threshold_alarm | 0.2 | -0.6 | 0.8 | 0.436 |  |  |  |  | <NA> |
| long_gaps_low | 0.2 | -0.8 | 0.4 | 0.641 | 0.0508 | 0.59 | 0.223 | 1.18 | True |
| long_gaps_high | 0.2 | -0.8 | 0.4 | 0.645 | 0.05 | 0.595 | 0.237 | 1.2 | True |

Per-fault per-day harm: c1_per_fault.csv.

## Secondary and exploratory (labelled, uncorrected, 95% CI)

Gate optima (in-sample, descriptive; neither is a recommended setting):

| harm_optimal | alarm_optimal | alarm_rule_feasible | coincide_share |
|---|---|---|---|
| no_gate | trial | True | 0.001 |

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

## Run

```json
{
 "commit": "69c6c1364afc1033b0d0edd9cb84fb99e678a010",
 "argv": [
  "--input",
  "results/grid/outcomes",
  "--protocol-tag",
  "protocol-v1",
  "--exclude-outside-acceptance",
  "--no-figures",
  "--out",
  "results/analysis_s2"
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
 "excluded_outside_acceptance": [
  "sim:adolescent#002",
  "sim:adolescent#007",
  "sim:child#002",
  "sim:child#003",
  "sim:child#008"
 ],
 "additivity": null,
 "seconds": 20.8
}
```
