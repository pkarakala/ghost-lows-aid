# Summary

**Question.** When a glucose sensor gives a wrong reading, how much harm follows if open-source insulin-dosing software acts on it? And are the sensor faults that cause the most false alarms the same ones that cause the most insulin harm?

**Why it matters.** A continuous glucose monitor (CGM) is a small sensor that reads glucose every 5 minutes. Automated insulin delivery (AID) systems use those readings to raise or lower insulin on their own. A wrong reading can therefore change the insulin a person gets. The best-known wrong reading is a *compression low*: lying on the sensor at night squeezes the tissue around it, and the reading drops although blood glucose has not. Our first study, [Ghost Lows](https://github.com/pkarakala/ghost-lows), showed how such faults set off false low alarms. This study follows the same faults into the insulin.

**Design in plain words.**
- *Virtual patients.* 30 simulated people (10 children, 10 adolescents, 10 adults) from the UVA/Padova type 1 diabetes simulator, eating the same meals as in Ghost Lows. Their true glucose is always known.
- *Dosing software.* OpenAPS oref0 0.7.1, run exactly as released, in two modes: temporary basal rates only (the main one) and with small automatic boluses (SMB). A third arm with no automation, just the person's scheduled basal insulin and meal boluses, is the reference.
- *Settings.* Each patient's insulin settings follow one fixed rule chosen before any fault was simulated. 25 of the 30 patients then meet the international targets for glucose control (time in range, time low); the other 5 stay in and are checked separately.
- *Faults.* 17 kinds and sizes of sensor error: compression lows, gaps in the data, frozen values, spikes, slow drift up or down, gain errors, and false highs from drugs such as acetaminophen in some sensors.
- *Paired runs.* Each simulated day runs once without faults. At a chosen moment it splits: one copy gets one fault, the other continues clean. The 12 hours after the split are compared on true glucose. Because meals, sensor noise and the person are identical, the difference is the fault's effect.
- *Size.* 30 patients × 10 random seeds × 3 controllers, with every fault arm.

**Main tests.** These two tests are fixed in advance:
1. **H1:** does a typical compression low at night lead to meaningful extra time above 180 mg/dL, because insulin was held back? "Meaningful" is 72 extra minutes or more in the 12 hours.
2. **H2:** does a large drug-induced false high lead to meaningful extra time below 70 mg/dL, because too much insulin was given? "Meaningful" is 15 extra minutes or more.

Each test ends in one of four readings: meaningful harm, negligible (the effect is shown to be smaller than meaningful), opposite direction, or inconclusive.

**Main comparison (C1).** Faults are ranked twice, once by the false alarms they cause per day and once by the insulin harm they cause per day. The study reports how far the two rankings agree.

**What each outcome would mean.**
- *Harm found:* the fault matters for AID users, and these simulations say how large the harm is and through which route (the software, the person, or both).
- *Negligible:* the software absorbs the fault in this simulator.
- *The rankings disagree:* work aimed at reducing false alarms is not the same as work aimed at reducing insulin harm.

Every result is a simulation result, not a clinical finding.

**Main limitations.**
- All patients are virtual. The 2008 simulator model has known gaps, for example in how it handles exercise and illness.
- One open-source algorithm in two modes. Commercial systems are not studied, and nothing here describes their safety.
- Fault sizes and rates come from the literature where it exists and are our assumptions where it does not.
- No clinician has reviewed the protocol yet. The open clinical questions are listed in [CLINICAL_REVIEW.md](CLINICAL_REVIEW.md).

**Where the exact specification lives.** [docs/reference/](reference/): the questions and hypotheses are in RESEARCH.md, the simulation in SIMULATION.md, the analysis in BENCHMARK.md, and every decision and its reason in DECISIONS.md.
