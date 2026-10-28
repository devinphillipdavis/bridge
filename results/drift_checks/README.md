# Drift checks: results

Pre-registration: [`PREREGISTRATION.md`](PREREGISTRATION.md), sha256 `655472339132…ec989`, with its
OpenTimestamps proof. It was frozen before both runs. Code: `bridge/drift_check.py` and `osiris_governance.telemetry`, pinned by hash
in the pre-registration and recorded again in each run's output.

## E1: simulation. All five predictions pass at both base rates

See `simulation/oc.json` and `simulation/verdict.json`: 26 cells, each 1000 trials × 400 replicas.

| prediction | alarm rates (both base rates) | verdict |
|---|---|---|
| P1: under the null, every predictive check alarms in ≤ 0.064 of trials | 0.042–0.057 | PASS |
| P2: the matched check alarms in ≥ 0.80 of trials at the largest linear, step and walk sizes | 0.918–1.000 | PASS |
| P3: under 4 SD jitter, the time-order checks stay ≤ 0.064 … | 0.042–0.056 | PASS |
| … while dispersion alarms in ≥ 0.80 | 1.000 | PASS |
| P4: under the same jitter, the squared-slope check alarms in > 0.064 | 0.595–0.633 | PASS |
| P5: the fixed-limit CUSUM (h = 27) alarms in < 0.50 under an 8 SD linear drift | 0.000 | PASS |

What this says about the claims from the paste:
- **C1 (the checks detect drift):** supported in simulation.
- **C2 (time structure means drift, not noise misspecification):** holds only for checks that ignore scale. The squared slope (P4) mistakes job-to-job jitter for drift. So does a CUSUM on model-standardised residuals with h = 5, which alarms in 0.754–1.000 of trials under 2–4 SD jitter against 0.002–0.003 under the null (exploratory).
- **C3 (h = 27 flags drift):** not supported at this design (K = 20, allowance 0.5). h = 27 never alarmed on an 8 SD linear drift. It fired only on the 8 SD step (0.999–1.000) and the largest random walk (0.259–0.322).

## E2: ibm_fez hardware, 2026-10-27 14:47–16:22Z

- **Runs:** 20 of 20 jobs finished, using 100 QPU-seconds. Qubits 4, 22, 54 and 82; each job ran `prep0`, `prep1` and `x100` with 4000 shots.
- **Calibration:** one calibration hash (`d54e808bae58`) throughout. The session check therefore didn't apply, and there were 48 tests.
- **Ledger:** `hardware/ledger.jsonl` verifies. Every job's `submit-intent` row precedes its submission. The rows are unsigned.
- **Control:** the same 20 jobs on Aer with ibm_fez's static noise model flagged **0 of 48**.
- **Hardware:** **18 of 48** flagged, with S = 4000 replicas and Holm's correction at family-wise error 0.05. The smallest possible p is 1/4001 ≈ 0.00025.

| hypothesis | result |
|---|---|
| HD: dispersion flagged on at least one series (variation between jobs beyond the binomial model) | **TRUE**, on 7 of 12 series |
| HT: a time-order check flagged on at least one series | **TRUE**, on 6 of 12 series |

| series | error rate | variance ÷ binomial | slope (pts/h) | p: disp / lag-1 / trend / CUSUM | flagged after Holm |
|---|---|---|---|---|---|
| `q4_prep0` | 0.116% | 4.2 | -0.151 | 0.00025 / 0.02 / 0.005 / 0.0045 | disp |
| `q4_prep1` | 0.684% | 1.4 | -0.178 | 0.13 / 0.1 / 0.022 / 0.019 | — |
| `q4_x100` | 12.411% | 167.4 | +6.343 | 0.00025 / 0.00025 / 0.057 / 0.0065 | disp, lag-1 |
| `q22_prep0` | 0.126% | 0.9 | +0.009 | 0.55 / 0.099 / 0.75 / 0.31 | — |
| `q22_prep1` | 0.671% | 1.3 | -0.049 | 0.16 / 0.17 / 0.55 / 0.084 | — |
| `q22_x100` | 5.958% | 211.1 | +7.809 | 0.00025 / 0.0005 / 0.0035 / 0.0012 | disp, lag-1, CUSUM |
| `q54_prep0` | 0.338% | 2.1 | -0.118 | 0.0075 / 0.092 / 0.08 / 0.0067 | — |
| `q54_prep1` | 0.946% | 1.2 | -0.113 | 0.24 / 0.46 / 0.19 / 0.61 | — |
| `q54_x100` | 36.914% | 140.3 | +9.530 | 0.00025 / 0.00025 / 0.035 / 0.002 | disp, lag-1 |
| `q82_prep0` | 0.365% | 2.9 | +0.283 | 0.00025 / 0.00025 / 0.00025 / 0.00025 | disp, lag-1, trend, CUSUM |
| `q82_prep1` | 1.413% | 7.1 | +0.867 | 0.00025 / 0.0042 / 0.00025 / 0.00025 | disp, trend, CUSUM |
| `q82_x100` | 11.925% | 192.2 | +8.719 | 0.00025 / 0.00025 / 0.012 / 0.0015 | disp, lag-1, CUSUM |

"slope" is the least-squares change in error rate per hour, in percentage points.

**What the flags show** [OBSERVED]:
- **x100 drifts on all four qubits.** The variance between jobs is 140–211 times what shot noise allows, and lag-1 autocorrelation is flagged on every qubit.
- **The four qubits move together.** Errors are high in jobs 1–2, roughly flat through job 14, then rise together over jobs 15–20:

  | qubit | median of jobs 3–13 | job 15 | job 20 |
  |---|---|---|---|
  | 54 | 1264 | 1356 | 2388 |
  | 4 | 332 | 444 | 1236 |
  | 82 | 290 | 501 | 1276 |
  | 22 | 113 | 234 | 857 |

- **Qubit 82's readout drifted.** `prep1` errors stepped from about 30 to about 70 per 4000 at job 7, and `prep0` errors rose from 4–10 to about 20.
- **Qubit 4's `prep0` varied but did not trend.** It was high in jobs 1–5 (up to 18), then 1–5 for the rest of the run.
- **Qubits 22 and 54 read out cleanly.** Their readout is consistent with shot noise. Qubit 54's `prep0` dispersion has p = 0.0075, which does not survive Holm's correction.
- **No recalibration was published.** Everything above happened while ibm_fez kept the same calibration.

**What this means for operational directive 5.** Its premise holds on this backend. Pooling shots across jobs here would understate the x100 standard error by a factor of about √140–√211, roughly 12–15.

**Possible cause** [HYPOTHESIS]. Because all four qubits moved together, the cause is likely shared, such as a common drive, electronics or thermal change. Two-level-system defects act on single qubits, so they fit this pattern poorly. E3, pre-registered separately, tests whether the change is the X gate's rotation error or the circuit's contrast.

## E3: x100 amplitude drift, ibm_fez, 2026-10-27 16:38–17:33Z

Pre-registration: [`AMPLITUDE_PREREGISTRATION.md`](AMPLITUDE_PREREGISTRATION.md), sha256 `261586b7aa41…3da`, with its OpenTimestamps proof. It was frozen before any E3 job. Code: `bridge/amplitude_drift.py`, pinned by hash in the pre-registration and in the ledger's first row.

- **Runs:** 12 of 12 jobs finished, using 96 QPU-seconds, on E2's backend and qubits (4, 22, 54, 82). Each job ran `prep0`, `prep1`, eight `ampN` circuits (√X, then N X gates) and `x100`, with 2000 shots each.
- **Timing:** jobs were submitted 300 s apart, but IBM queued jobs 3–8, which all ran between 17:11:40 and 17:13:35Z. The time-order checks use execution midpoints, as pre-registered.
- **Calibration:** one hash (`d54e808bae58`) throughout, the same as E2's. The session check therefore didn't apply, and there were 32 tests.
- **Ledger:** `amplitude/ledger.jsonl` verifies. The rows are unsigned.
- **Software control** (`amplitude/control.json`, before dispatch): the constant control flagged nothing, so dispatch was allowed. The ε step flagged HA and HAt only; the contrast step flagged HC and HCt only.
- **Hardware:** **6 of 32** flagged, with S = 1000 replicas and Holm's correction at family-wise error 0.05. The smallest possible p is 1/1001 ≈ 0.001.

| hypothesis | result |
|---|---|
| HA: ε_k dispersion flagged on at least one qubit (the rotation error varies between jobs) | **TRUE**, on all 4 qubits |
| HAt: an ε_k time-order check flagged on at least one qubit | **TRUE**, trend on qubits 4 and 82 |
| HC: a_k dispersion flagged on at least one qubit (the contrast varies between jobs) | **FALSE** |
| HCt: an a_k time-order check flagged on at least one qubit | **FALSE** |
| prediction: HA on at least one of qubits 4, 54 and 82 | **met**, on all three |

**Outcome: HA and not HC.** By the pre-registered table, the job-to-job variation in `x100` during E3 is a change in the X gate's coherent rotation error, not in contrast. This supports the hypothesis. [OBSERVED]

| qubit | pooled ε (rad / deg) | ε_k jobs 1–2 | median ε_k jobs 3–12 | a_k range | p ε: disp / lag-1 / trend / CUSUM | p a: disp / lag-1 / trend / CUSUM | flagged after Holm |
|---|---|---|---|---|---|---|---|
| 4 | +0.0080 / +0.46° | +0.0128, +0.0127 | +0.0075 | 0.990–1.000 | 0.001 / 0.005 / 0.001 / 0.022 | 0.011 / 0.67 / 0.63 / 0.22 | ε disp, ε trend |
| 22 | +0.0052 / +0.30° | +0.0099, +0.0098 | +0.0045 | 0.990–0.995 | 0.001 / 0.005 / 0.002 / 0.041 | 0.59 / 0.67 / 0.32 / 0.38 | ε disp |
| 54 | +0.0123 / +0.70° | +0.0098, +0.0193 | +0.0128 | 0.980–0.990 | 0.001 / 0.72 / 0.26 / 0.69 | 0.13 / 0.64 / 0.08 / 0.28 | ε disp |
| 82 | +0.0028 / +0.16° | +0.0118, +0.0127 | +0.0015 | 0.960–0.980 | 0.001 / 0.019 / 0.001 / 0.018 | 0.002 / 0.11 / 0.087 / 0.007 | ε disp, ε trend |

| qubit | pooled c | pooled a | pooled γ | pooled δ (rad) | x100 errors, jobs 1–12 (of 2000) | Spearman x100 vs sin²(50ε_k) |
|---|---|---|---|---|---|---|
| 4 | 0.498 | 0.993 | 0.00000 | -0.040 | 745, 697, 348, 368, 369, 336, 336, 301, 269, 200, 269, 339 | +0.93 |
| 22 | 0.498 | 0.993 | 0.00000 | -0.058 | 492, 458, 165, 183, 179, 162, 145, 147, 129, 86, 182, 240 | +0.94 |
| 54 | 0.496 | 0.991 | 0.00011 | -0.048 | 1203, 1356, 872, 893, 918, 835, 843, 857, 720, 699, 722, 794 | +0.59 |
| 82 | 0.499 | 0.972 | 0.00000 | -0.087 | 782, 840, 147, 166, 164, 175, 170, 205, 144, 190, 192, 199 | +0.55 |

**Near miss on HC.** Qubit 82's contrast dispersion has p = 0.0020. Holm's procedure stopped at rank 7 (qubit 22's ε trend, p = 0.0020 against 0.05/26 = 0.0019), so it is not flagged and HC stays FALSE. Qubit 82's a_k fell from 0.98 to 0.96–0.97 over the run, in line with its readout drift in E2.

**What moved** [OBSERVED]:
- **The X gate over-rotates on all four qubits** (ε > 0).
- **ε dropped between job 2 (16:43Z) and job 3 (17:11Z)** on qubits 4, 22 and 82, by 0.005, 0.005 and 0.011 rad. Their `x100` errors halved or more (qubit 82: 840 → 147). On qubits 4, 22 and 54 the contrast stayed within 0.01.
- **Qubit 54** read ε = 0.0098 in job 1, 0.0193 in job 2 and about 0.013 afterwards. Its job-1 estimate disagrees with its `x100` count (see the post hoc check below).

**Exploratory** (pre-registered, not judged):
- **Link to `x100`:** the Spearman correlation between `x100` and sin²(50ε_k) is positive on all four qubits, strongly on qubits 4 and 22.
- **δ and γ:** δ is −0.04 to −0.09 rad on all four (√X under-rotates). γ sits at its lower bound 0 on three qubits and is 0.0001 on qubit 54, so no contrast loss per X gate is resolved.
- **Comparison with E2,** as `x100` error rates (E2 used 4000 shots, E3 2000):

  | qubit | E2 jobs 1–2 (14:47–15:00Z) | E2 median jobs 3–13 (15:05–15:52Z) | E2 jobs 19–20 (16:11–16:21Z) | E3 jobs 1–2 (16:38–16:43Z) | E3 median jobs 3–12 (17:11–17:33Z) |
  |---|---|---|---|---|---|
  | 4 | 16.6% | 8.3% | 28.9% | 36.0% | 16.8% |
  | 22 | 4.0% | 2.8% | 19.2% | 23.8% | 8.2% |
  | 54 | 43.7% | 31.6% | 56.7% | 64.0% | 42.0% |
  | 82 | 12.1% | 7.2% | 29.3% | 40.5% | 8.6% |

  In clock order, `x100` was high at 14:47–15:00Z, lower from 15:05Z, rising from about 15:57Z and still high at 16:38–16:43Z. By 17:11Z it had fallen again, though on all four qubits it stayed above E2's mid-run plateau. So both runs began on a high stretch of one continuous timeline. Nothing here points to an effect of a run's start.

**Post hoc check** (not pre-registered) [OBSERVED]. Each job's ε_k and a_k, with the pooled c and γ, predict `x100` as P(1) = c − (a/2)·e^(−100γ)·cos(100ε). The observed counts exceed that prediction:
- In jobs 3–12, the median excess is 62, 58, 133 and 135 per 2000 on qubits 4, 22, 54 and 82.
- On qubit 82 the model predicts 36 errors and about 170 are seen.
- On qubit 54's job 1 the excess is 758, because the per-job ε there is far below what its `x100` count implies.

So the rotation error accounts for most of the *change* in `x100`, but not for all of its level. `x100` carries another 3–7 percentage points of error that the `ampN` circuits do not pick up.

**Possible cause** [HYPOTHESIS]. Within E3, the X pulse amplitude changed on several qubits at once, with no change in the calibration hash. `dnalang.backends.ibm.calibration_hash` covers reported gate errors and durations but not pulse amplitudes. So either a slow drift in drive amplitude or an unpublished amplitude re-tune would look like this. Telling them apart needs a denser time series or the backend's pulse amplitudes.

**What this means for operational directive 5.** For long single-qubit gate sequences, a between-job model needs a rotation error for each job. A term for readout or contrast alone would not capture the variation seen here.

## E4: pulse-amplitude test, ibm_fez, submitted 2026-10-27 22:45–23:58Z, run 2026-10-28 00:16–00:37Z

Pre-registration: [`E4_PREREGISTRATION.md`](E4_PREREGISTRATION.md), sha256 `7e20ee97bdff…ae07`, with its OpenTimestamps proof. It was frozen and pushed (`7be9176`) before any E4 job. Code: `bridge/pulse_amplitude.py` (`431602d8df17…8955`), pinned by hash in the pre-registration and in the ledger's first row.

- **Runs:** 30 of 30 jobs finished, using 150 QPU-seconds. Qubits 4, 22, 54 and 82 (E3's) plus 144, 24, 14 and 71, chosen by `choose_extras`. Each job ran 12 circuits with 1000 shots each. The Ramsey delay was 200 ns.
- **Timing:** jobs were submitted 150 s apart, but IBM queued all 30. They ran between 00:16:48 and 00:36:33Z, 4–209 s apart (median 18 s). The analysis uses execution times, as pre-registered. All 29 intervals were 300 s or shorter, so H3 was not underpowered.
- **Ledger:** `pulse_amplitude/ledger.jsonl` verifies (122 rows). A `calibration-snapshot` row precedes every `submit-intent`. The rows are unsigned.
- **Software control** (`pulse_amplitude/control.json`, before dispatch): the constant control flagged nothing (H4 p = 0.51), so dispatch was allowed.
- **Hardware:** S = 2000 replicas, with Holm's correction at family-wise error 0.05. The smallest possible p is 1/2001 ≈ 0.0005.

| hypothesis | result |
|---|---|
| H1: ε_k dispersion flagged on at least one qubit | **TRUE**, on 7 of 8 qubits. Qubit 71 is not flagged (p = 0.0085 against a Holm bar of 0.0071). |
| HC (reported): a_k dispersion flagged | **TRUE**, on qubits 4 and 22 |
| H2: β = Δδ/Δε, pooled over the 7 flagged qubits | **0.438**, 95% interval [0.373, 0.504]: the **gain** region |
| H2F: frequency does not drive ε | **fails on qubits 24 and 82**; passes on 4, 14, 22, 54 and 144 |
| H3: jumps beyond a smooth drift | **TRUE**, on qubit 82 only (max\|z\| = 4.29, p = 0.0005), between jobs 18 and 19 |
| H4: the changes are shared across qubits | **TRUE**: the mean pairwise correlation of ε increments is 0.116, p = 0.0045. The coincidence test did not apply, because H3 flagged one qubit. |

**Pre-registered reading** (`pulse_amplitude.reading`): "step changes in a shared gain, or a joint x and sx re-tune; H2F falsifies amplitude on qubits [24, 82]". So the amplitude hypothesis is supported with its mechanism open, except on qubits 24 and 82, where the frequency test falsifies it. [OBSERVED]

| qubit | pooled ε (rad) | pooled δ (rad) | ε_k range | ε_k SD ÷ measurement SD | a_k range | p: ε disp / a disp | H2 β (weight share) | max\|Δf\| (kHz) | Spearman Δf vs ε_k (p) | `x100` errors, median of 1000 |
|---|---|---|---|---|---|---|---|---|---|---|
| 4 | +0.0059 | +0.0016 | +0.0034 to +0.0093 | 2.4 | 0.887–0.993 | 0.0005 / 0.0005 | 1.46 (0.00) | 118 | +0.26 (0.16) | 183 |
| 22 | +0.0021 | +0.0001 | −0.0001 to +0.0040 | 1.9 | 0.978–0.999 | 0.0005 / 0.0010 | 0.31 (0.05) | 51 | +0.36 (0.050) | 97 |
| 54 | +0.0117 | +0.0042 | +0.0098 to +0.0132 | 1.9 | 0.977–0.993 | 0.0005 / 0.22 | 0.73 (0.06) | 50 | −0.14 (0.49) | 524 |
| 82 | +0.0017 | −0.0017 | −0.0022 to +0.0056 | 4.8 | 0.975–0.993 | 0.0005 / 0.40 | 0.41 (0.73) | 223 | +0.64 (0.0005) | 119 |
| 144 | +0.0049 | +0.0021 | +0.0034 to +0.0062 | 1.5 | 0.983–0.997 | 0.0010 / 0.60 | 0.17 (0.03) | 56 | +0.11 (0.55) | 190 |
| 24 | +0.0021 | +0.0039 | +0.0002 to +0.0038 | 1.9 | 0.980–0.997 | 0.0005 / 0.095 | 0.58 (0.06) | 56 | −0.47 (0.006) | 122 |
| 14 | +0.0056 | +0.0024 | +0.0033 to +0.0071 | 1.9 | 0.989–0.998 | 0.0005 / 0.90 | 0.58 (0.07) | 72 | −0.04 (0.85) | 230 |
| 71 | +0.0064 | −0.0003 | +0.0046 to +0.0076 | 1.3 | 0.984–0.997 | 0.0085 / 0.58 | — | — | — | 291 |

Pooled γ sits at its lower bound 0 on all eight qubits. Measurement SDs of ε_k are 0.0005–0.0006 rad. Ramsey Δf is relative to each qubit's run mean.

**What moved** [OBSERVED]:
- **ε varies between jobs on seven of eight qubits,** by 1.5–4.8 times its measurement SD, and the X gate over-rotates on all eight (pooled ε +0.0017 to +0.0117 rad).
- **Qubit 82 stepped once,** between job 18 (ending 00:33:31Z) and job 19 (starting 00:33:57Z). Medians of jobs 1–18 against jobs 19–30:
  - ε: −0.0006 → +0.0048 rad, so \|ε\| grew;
  - δ: −0.0026 → −0.0006 rad;
  - Ramsey frequency: −133 → +207 kHz, a step of about 340 kHz;
  - `x100` errors: 85 → 235 per 1000.
- **Qubit 4's readout of \|1⟩ degraded from job 6** (about 00:19:40Z). The \|1⟩ error e1 = 1 − c − a/2 rose from 0.003–0.010 in jobs 1–5 to 0.018–0.108, while e0 stayed below 0.007. This is HC's flag on qubit 4. Qubit 22's contrast flag reflects a small spread (0.978–0.999).
- **√X:** with E4's √X ladders, pooled δ is within ±0.005 rad on all eight qubits. E3 reported δ of −0.04 to −0.09 rad from a design whose δ had a Fisher SD of 0.012 rad and absorbed the static misfit (groundwork G2, G3). E4's values supersede it as a measure of δ.

**Secondary checks** (pre-registered, reported, not judged):
- **\|ε\| across the flagged jump:** it grew on qubit 82 (0.0004 → 0.0049 rad), against the re-tune prediction that a re-tune moves ε toward zero.
- **Calibration snapshots dated inside qubit 82's jump interval:** none. But the snapshots were taken at submission (22:45–23:58Z), before any job ran, so they cannot show an entry dated in the execution window (see Deviations). Across the 30 snapshots, no recorded field on any of the eight qubits changed. The backend-wide `last_update_date` moved twice, to 22:50:54 and 23:19:32Z.
- **E4a** (`pulse_amplitude/retrospective/`, 0 QPU-seconds): the dated properties query works and returned different property sets over E2 and E3. For qubits 4, 22, 54 and 82 no recorded field changed or was dated inside either change interval (`dated_inside` is empty for both). Only the backend-wide `last_update_date` moved, from 16:17:30 to 16:50:32Z, which falls inside E3's change interval. In E2 it moved from 14:41:17 to 14:52:51Z between jobs 5 and 6, outside E2's interval. E2's `x100` rise over jobs 15–20 came under one unchanged set of reported properties.
- **Residual pattern (G2 repeated):** the mean residual z per circuit is within ±1.6 on every qubit, with the same signs on all eight (`amp16`, `amp32` and `sx17` low; `amp48` and `sx49` high). The SD across jobs is 0.4–1.2 on the gate circuits, near shot noise. Qubit 54, whose pattern was not static in E3, is static here.

**Post hoc checks** (not pre-registered) [OBSERVED]:
- **Qubit 82 carries 73% of H2's pooled weight.** Pooled the same way without it, β = 0.508; without qubits 82 and 24, β = 0.488. Per-qubit β on the five qubits that pass H2F is 0.17–0.73 (qubit 4's 1.46 has weight ≈ 0).
- **H4 without qubit 82:** mean correlation 0.114, p = 0.018. This uses the same function and seed with a fresh generator, which gives p = 0.007 for all eight qubits against the recorded 0.0045.
- **What H2F's failures detect.** On qubit 82, ε and frequency step at the same interval, and the rank correlation registers that coincidence.
  - Under the pre-registration's square-pulse derivation (32 ns X assumed), a 340 kHz shift changes ε by 5 mrad only if qubit 82 was already about 1 MHz off resonance [DERIVED]. Ramsey here measures only changes in frequency, so that cannot be checked.
  - Qubit 82's β (about 0.4) is near gain's 0.5, not detuning's 2.
  - On qubit 24 the frequency moved by at most 56 kHz, too little to change ε through detuning unless the qubit was more than about 1 MHz off resonance [DERIVED]. A slow trend shared by both series would give the same correlation, and the pre-registered test does not remove one.

**Interpretation** [HYPOTHESIS]:
- **On the pre-registered reading,** E4 supports the amplitude hypothesis with its mechanism open, and falsifies it on qubits 24 and 82.
- **With the post hoc checks,** the five qubits that pass H2F show correlated ε changes, a ratio near gain's 0.5 (0.49–0.51 pooled; E4-sim's count-level runs read a true 0.5 as about 0.45) and no jumps. That is the pattern pre-registered for a shared analog gain drift. It is post hoc and not judged here; testing it needs a new pre-registration.
- **Qubit 82's step is a separate event:** a simultaneous change in one qubit's frequency and rotation error, with \|ε\| growing and no reported property change in snapshots that predate it. A two-level-system defect near qubit 82 or a recalibration of that qubit alone would both fit; E4 cannot tell them apart [INCONCLUSIVE].
- **E4 did not see a shared drop like E3's** (0.005–0.011 rad). Outside qubit 82, the ε_k SD is 0.0007–0.0015 rad.

**What this means for operational directive 5.** As in E3, a between-job model for single-qubit gate sequences needs a rotation error for each job. Qubit 4 shows that readout can also change within 20 minutes, so a term for each job's readout is needed too.

## Deviations and disclosures

- **E2 interim look.** At about 15:40Z, 6 jobs into E2, I took a read-only, descriptive look at the finished jobs. It wrote nothing and changed nothing in the plan or the analysis. It is not in the pre-registration.
- **No E3 interim look.** During E3 I read only job status. I first looked at E3 counts at 17:36Z, after the run and its analysis had finished.
- **E3 queueing.** IBM queued jobs 3–8, so they ran within two minutes of each other rather than 300 s apart.
- **E3 post hoc check.** The `x100` prediction check above was added after the results were in. It is labelled as post hoc and did not change any hypothesis.
- **No E4 interim look.** During E4b I read only job status: the submission log and the ledger's row kinds. The user asked twice what E4 showed while it ran; I did not look. I first read E4 results at about 00:38Z, after the run and its analysis had finished.
- **E4 queueing.** IBM queued all 30 jobs, so they ran 4–209 s apart within 20 minutes, rather than 150 s apart over 75 minutes. The pre-registration's Limits anticipated queueing; H3 needed 15 intervals of 300 s or less and had 29.
- **E4 calibration snapshots predate execution.** The snapshots were taken before each submission, as pre-registered, but the jobs ran 39–91 minutes later. So they do not cover the execution window, and the secondary check of dates inside qubit 82's jump interval is uninformative. A design for a queued backend would need snapshots at execution time.
- **E4 field name.** `abs_eps_before_after` in `analysis.json` holds the *signed* ε values before and after the jump. `abs_eps_decreased` compares their absolute values, correctly.
- **E4 post hoc checks.** The weight shares, β and H4 without qubits 82 and 24, the step medians on qubit 82, qubit 4's readout errors and the detuning arithmetic were added after the results were in. They are labelled as post hoc and do not change the pre-registered reading.
- **Unsigned ledger.** The ledger rows in all three hardware runs are hash-chained but unsigned, which deviates from directive 4. They need a signature from an enrolled node key.
  - [`ATTESTATION_E2-E4.json`](ATTESTATION_E2-E4.json), written on 2026-10-28 after the runs by `python -m bridge.attest build`, binds each run's pre-registration hash, ledger (SHA-256, rows, final chain hash), telemetry head and every file in its folder.
  - It is unsigned until the holder of an enrolled node key runs `python -m bridge.attest sign`, which writes `ATTESTATION_E2-E4.json.sig`. `python -m bridge.attest verify` checks the signature and recomputes every hash.
  - A signature on it is retrospective. It shows that the key holder vouched for these files when signing, not that any row was signed while its job ran or that no row is missing.
- **Base rates.** E1's base rates (0.5% and 2%) are above E2's readout error rates (0.1–0.4%). E1's operating characteristics therefore apply only approximately to the readout series.
- **QPU budget.** E2 used 100 s, E3 96 s and E4 150 s, so 346 s of the open plan's 600 s per 28 days. The account reported 398 s remaining at E4b's start (6 s fewer than 600 − 196), so about 248 s should remain.
- **Agent rules.**
  - E2 used IBM Quantum, lasted about 95 minutes, and its thresholds were chosen by the agent. The user asked for the run and approved letting it continue past the node horizon.
  - E3 was requested in the user's own words. With its software control, it ran about 70 minutes, past the node horizon, as flagged in its pre-registration.
  - E4's values, the agent's build choices, its duration (about 110 minutes from submission to analysis, past the node horizon) and its QPU budget were approved by the user before the freeze. The user started E4a and E4b ("approve the agent choices, freeze it and run E4").
  - The thresholds still need the human's review: α and the 0.80 and 0.50 in E2's predictions, and the replica count and control step sizes in E3's.
