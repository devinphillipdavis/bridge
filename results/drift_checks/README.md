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

## Deviations and disclosures

- **Interim look.** At about 15:40Z, 6 jobs into E2, I took a read-only, descriptive look at the finished jobs. It wrote nothing and changed nothing in the plan or the analysis. It is not in the pre-registration.
- **Unsigned ledger.** The ledger rows are hash-chained but unsigned, which deviates from directive 4. They need a signature from an enrolled node key.
- **Base rates.** E1's base rates (0.5% and 2%) are above E2's readout error rates (0.1–0.4%). E1's operating characteristics therefore apply only approximately to the readout series.
- **Agent rules.**
  - The run used IBM Quantum, lasted about 95 minutes, and its thresholds were chosen by the agent.
  - The user asked for the run and approved letting it continue past the node horizon.
  - The thresholds (α, and 0.80 and 0.50 in the predictions) still need the human's review.
