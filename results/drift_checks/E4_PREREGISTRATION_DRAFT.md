# E4, pulse-amplitude test: pre-registration DRAFT

Status: **draft, not frozen.** Nothing on hardware has run.
- **Values:** the proposed values (K, cadence, shots, S, α, the β bands, the Ramsey limit and the analysis seed) were approved by the human on 2026-10-27 ("approve the proposed values and build it").
- **Built:** the code and tests are written, and E4-sim has run (results below).
- **Details for review:** implementation details chosen by the agent during the build are listed under "Agent choices to review".
- **Freezing means:** recording the SHA-256 of the final file (renamed `E4_PREREGISTRATION.md`), stamping it with OpenTimestamps, and pinning the code by hash.

## Why

E3 found that, on ibm_fez qubits 4, 22, 54 and 82, the job-to-job variation in `x100` is a change in the X gate's coherent rotation error ε, not in contrast (HA TRUE, HC FALSE) [OBSERVED]. ε dropped between 16:43 and 17:11Z on three qubits, while the calibration hash stayed the same.

E3's README states the hypothesis tested here [HYPOTHESIS]:
- The X pulse amplitude changed, through either a slow drift in drive amplitude or an unpublished amplitude re-tune.
- `dnalang.backends.ibm.calibration_hash` covers reported gate errors and durations, not pulse amplitudes, so it would not register either one.

### Groundwork on E3's committed data

This was done on 2026-10-27 while designing E4. It is exploratory and post hoc. The scripts were scratch files and are not part of E4's code.

| | finding | consequence for E4 |
|---|---|---|
| G1 | Jobs 3–8 ran within 2 minutes (17:11:40–17:13:35Z). Their ε_k varies with SD 0.00015–0.00037 rad, which matches the Fisher SD of E3's design (0.00025). In jobs 8–12, about 5 minutes apart, ε_k moved by up to 0.0026 rad (qubit 22). | ε is stable over 2 minutes to within measurement noise and wanders over 5 minutes. A cadence under 5 minutes is needed to see the shape of a change. |
| G2 | E3's model does not fit. The per-job Pearson χ² has a median of 33–54 against about 8 expected, with the same pattern on all four qubits: `amp0` and `amp100` read high, `amp32`–`amp64` low. On qubits 4, 22 and 82 the pattern is **static**: residuals vary across jobs with SD 0.6–1.5, close to shot noise. On qubit 54 it does not (SD up to 4.7). Neither a leakage-like term nor an offset between the X and √X axes removes it. | E4 tests **changes** in each parameter from job to job, and not absolute values. A static misfit biases every job's fit alike, so it cancels in the changes. This was checked in the tests with E3's residuals injected (see Estimation). Absolute ε and δ are reported but not judged. |
| G3 | E3's design cannot estimate the √X error δ per job (Fisher SD 0.012 rad). Adding √X ladders brings it to about 0.0006 rad at 1000 shots. | E4 adds √X ladders, which makes the ratio test below possible. |
| G4 | The installed `qiskit-ibm-runtime` (0.50.0) exposes `backend.properties()`: per-gate error and length, and per-qubit frequency, T1, T2 and readout, each with a `date`. It does not expose pulse amplitudes. Its documentation says `properties(datetime=…)` may raise `NotImplementedError`. | E4 cannot read amplitudes; it tests what the amplitude hypothesis predicts. It records a calibration snapshot before each job, and treats the retrospective query (E4a) as best-effort. |

## Question

Are the job-to-job changes in apparent X rotation error changes in drive amplitude? If they are, are they a slow drift of a shared analog gain, or discrete re-tunes?

## What each explanation predicts

### How the predicted ratios are derived [DERIVED]

**Gain.** A relative drive-gain change g scales every pulse on the line:
- X rotates by π(1+g), so Δε = πg.
- √X rotates by (π/2)(1+g), so Δδ = (π/2)g.

Hence **β = Δδ/Δε = 1/2**.

**X-only re-tune.** A re-tune that changes only the X amplitude gives **β = 0**.

**Detuning.** A detuning Δ, with X and √X of equal duration (the square-pulse approximation), gives:
- ε = πα²/2 and δ = πα², where α = Δ/Ω_X;
- hence **β = 2**.

Explaining Δε = 0.005 rad this way needs α ≈ 0.056, which is about 0.9 MHz for a 32 ns X pulse. The backend's reported pulse lengths are recorded in every calibration snapshot.

### Predictions

| observable | shared gain drift | unpublished re-tune of the X amplitude | detuning (alternative) |
|---|---|---|---|
| β = Δδ/Δε | ≈ 0.5 | ≈ 0 if X alone is re-tuned; any value if √X is re-tuned too | ≈ 2 |
| Ramsey frequency | flat | flat | moves with ε, by ~MHz |
| increments of ε_k | Gaussian, scaling with the time between jobs | mostly small, with a few large jumps | — |
| across qubits | correlated if the drive chain is shared | simultaneous jumps if a backend-wide cycle runs | — |
| reported properties | no dated update at the change | possibly a dated update at a jump | frequency update |
| \|ε\| across a jump | no preferred direction | decreases (a re-tune moves toward zero) | — |

The alternative already excluded by E3 is readout or contrast, which E3 found does not vary (HC FALSE). E4 re-tests it as part of H1.

## Design

### E4a: retrospective, 0 QPU-seconds

- **Query:** `backend.properties(datetime=t)` at the execution midpoint of each of the 32 jobs in E2 and E3.
- **Recorded, for qubits 4, 22, 54 and 82:** every reported qubit property, and X and √X error and length, each as [value, date]; plus `last_update_date` and the SHA-256 of the full properties.
- **If the call raises or returns nothing:** E4a records `available: false` and the exception's type (not its message, which may echo request details), then ends.
- **What the re-tune hypothesis predicts:** an entry dated inside E3's change interval (16:43:29–17:11:40Z, from the end of job 2 to the start of job 3) or E2's (15:00:16–15:05:31Z). Finding none is weak evidence, because amplitudes are not reported.
- **Command (from `bridge/`):**
  `python -m bridge.pulse_amplitude retrospective --key-file <key> --out results/drift_checks/pulse_amplitude/retrospective --runs results/drift_checks/hardware results/drift_checks/amplitude`

### E4b: hardware

- **Backend:** ibm_fez. If it is not operational at start, E4b is not run.
- **Qubits:** 4, 22, 54 and 82, plus four more. The extras are taken in order of lowest reported readout error at start (ties by index). A qubit is skipped if it is one of the first four, neighbours one of them in the coupling map, or neighbours an extra already chosen. (`choose_extras`)
  - Extra qubits cost no QPU time and show whether a change is shared across the chip.
- **Jobs:** K = **30**, submitted **150 s** apart, so about 75 minutes. Each submission goes through `dnalang.backends.ibm.submit`, ledger-first.
  - Before each `submit-intent`, a `calibration-snapshot` row records, for all eight qubits, the fields listed under E4a, from freshly refreshed properties.
- **Shots:** **1000** per circuit.
- **Ramsey delay:** 200 ns, rounded to the backend's timing granularity. The rounded value is recorded in the `design` row and `state.json` and used in the analysis.

**Circuits, 12 per job, all qubits in parallel:**

| circuit | gates | P(1) model | measures |
|---|---|---|---|
| `prep0`, `prep1` | measure; X, measure | c − a/2; c + (a/2)cos ε | readout c and contrast a |
| `ampN`, N ∈ {0, 16, 32, 48} | √X, then N X | c + (a/2) e^(−γN) sin(δ + Nε) | ε (and δ) |
| `sxM`, M ∈ {17, 33, 49} | M √X (M = 4m + 1) | c + (a/2) e^(−γM) sin(Mδ) | δ |
| `ram0`, `ram90` | √X, delay, [rz(π/2)], √X | c + r cos φ; c − r sin φ | frequency, ±2.5 MHz unambiguous |
| `x100` | 100 X | not fitted; kept for continuity with E2 and E3 | — |

The model was checked against exact rotation matrices in `tests/test_pulse_amplitude.py::test_model_matches_exact_rotations`.

**Precision per job at 1000 shots (Fisher, at E3's values) [DERIVED]:**
- ε: 0.00054 rad;
- δ: 0.00053–0.00060 rad;
- a: 0.0030;
- Ramsey frequency: about 25 kHz.

E3's change of 0.005–0.011 rad is 10–20 SD of ε. The 0.9 MHz detuning needed to mimic it is about 35 SD.

**Budget.** E2's 12,000-shot jobs were billed 5 s each, so E4b should cost about 150 s. After E2 and E3, about 404 s of the open plan's 600 s per 28 days remain, if nothing else has used the instance.

**Start and stop rules:**
- E4b is not started if the account reports fewer than 200 s remaining, or reports no figure (fail closed).
- At least 20 jobs must finish, or the run is reported as incomplete.

## Estimation (per qubit)

1. **Pooled fit:** maximum-likelihood (c, a, γ, δ, ε) over all jobs by L-BFGS-B, from a grid start:
   - c and a from `prep0`/`prep1`;
   - δ from `amp0` and the √X ladders on [−0.3, 0.3] in steps of 0.0005;
   - ε from the X ladders on [−0.1, 0.1] in steps of 0.0002.

   Readout is fitted as the errors e0 = c − a/2 and e1 = 1 − c − a/2, each in [0, 0.5], so no circuit's P(1) can leave [0, 1]. Bounds: γ ∈ [0, 0.01], δ ∈ [−0.3, 0.3], ε ∈ [−0.1, 0.1].
2. **Per job:** (c_k, a_k, δ_k, ε_k) by 5 weighted Gauss–Newton steps from the pooled fit, with γ held at the pooled value.
   - **Weights:** taken from the run's mean frequency of each circuit, smoothed by half a count, so none is unbounded.
   - **Misfit check:** with E3's residuals injected (scaled by 4p(1−p) so they vanish at P = 0 or 1), the tests recover:
     - a 0.005 rad ε step to within 0.0003 rad;
     - a gain change as β within 0.1 of 0.5, on all four sets of E3 parameters.

     The development check, near noiseless, gave β = 0.43–0.49.
3. **Ramsey:** φ_k = atan2(−(f90 − c_k), f0 − c_k) and Δf_k = wrap(φ_k − circular mean) / (2π · delay). Only changes, and their absolute values, are used, so the sign convention does not matter.
4. **Null replicas:** S = **2000** draws per qubit. Each job's counts are drawn binomially at the run's mean frequencies (no drift, the same static residual), then passed through steps 2–3. The pooled fit is held at the observed one.
   - **Seed:** 2026102801.
   - **Why S matters:** with S = 2000 the smallest p is 0.0005, below H1's smallest Holm bar (0.05/16 ≈ 0.0031).

## Confirmatory hypotheses

The hypotheses are tested in a fixed sequence. H2–H4 are evaluated only if H1 holds.

| | hypothesis | test |
|---|---|---|
| **H1** (gate; replicates E3's HA and HC) | ε_k dispersion is flagged on at least one qubit | Dispersion (`osiris_governance.telemetry.dispersion`, exact, micro-units) of ε_k and of a_k on each of the 8 qubits (16 tests), against the replicas, with Holm's correction at **α = 0.05**. If H1 is FALSE, E4 stops: "no variation in E4's window", which is not evidence of absence. HC (contrast flagged) is reported. |
| **H2** (mechanism) | β = Δδ/Δε | For each qubit whose ε dispersion is flagged: a Deming regression of δ_k on ε_k, with error-variance ratio set by the replica variances of δ_k and ε_k. Per-qubit slopes are pooled by inverse pairs-bootstrap variance (S draws, the same resampled jobs for every qubit), and the 95% interval is the percentile interval of the pooled slope. Regions: **gain** if the interval lies inside [0.25, 1.0]; **X-only** if it contains 0 and ends below 0.25; **detuning-like** if it starts above 1.25; otherwise **inconclusive**. |
| **H2F** (frequency exclusion) | detuning does not drive ε | On the same qubits: max\|Δf_k\| < **200 kHz**, and the Spearman correlation of Δf_k with ε_k is not flagged (two-sided, permutation null with S draws, Holm). If either fails, the amplitude hypothesis is **falsified** for that qubit. |
| **H3** (steps) | ε changes in jumps beyond a smooth drift | On the same qubits, fit a Gaussian random walk with drift to the ε_k series: drift μ = mean increment per second; rate q per second by moments from the mean squared increment, with measurement variances known. Standardise each increment by its own Δt; the statistic is max\|z_k\|. The null is a parametric bootstrap from the fitted walk, refitted the same way (S draws), with Holm's correction across qubits. If fewer than 15 intervals are 300 s or shorter, H3 is reported as **underpowered** and counted as no jumps. |
| **H4** (common cause) | the changes are shared across qubits | The mean pairwise correlation of the ε increments across all 8 qubits, against independent circular shifts of each qubit's increments (S draws), at α. If H3 flags two or more qubits, the largest number of flagged jumps in one interval is also tested against jumps placed at random (reported). |

**Secondary (reported, not judged):**
- whether \|ε\| is smaller after each flagged jump than before it;
- which `calibration-snapshot` dates fall inside a flagged jump's interval;
- E4a's result;
- the residual pattern per circuit (the G2 check, repeated).

## Reading the outcome

This is `pulse_amplitude.reading`, applied in order:

| condition | reading |
|---|---|
| H1 FALSE | not evaluable: no variation in rotation error in E4's window (not evidence of absence) |
| H2 detuning-like, or H2F fails on every tested qubit | **frequency, not amplitude**: falsifies the amplitude hypothesis |
| gain, H4 shared, no H3 jumps | **shared analog gain drift**: supports the amplitude hypothesis, as a slow drift |
| gain, H4 shared, H3 jumps | step changes in a shared gain, or a joint X and √X re-tune: amplitude supported, mechanism open |
| gain, not shared | gain change, not shown to be shared |
| X-only, H3 jumps, H4 shared | **X-amplitude re-tunes**: supports the unpublished re-tune |
| X-only, H3 jumps, not shared | X-amplitude jumps, not shown to be shared |
| X-only, no jumps | an X-only change that is not gain: mechanism inconclusive |
| inconclusive | amplitude not excluded; mechanism inconclusive |

If H2F fails on some tested qubits but not all, the reading names them.

## Simulation stage before dispatch (E4-sim)

`python -m bridge.pulse_amplitude simulate --out results/drift_checks/pulse_amplitude/simulation`

**Synthetic qubits.** Eight, using E3's four pooled parameter sets twice, with K = 30 jobs 150 s apart. E3's mean residuals are injected, scaled by 4p(1−p). The √X ladders reuse the amp16–48 residuals.

**Seeds:** control 2026102802, scenarios 2026102803.

| scenario | level, trials | expectation | consequence |
|---|---|---|---|
| S1: constant | counts, 1 (the software control) | H1, H3 and H4 flag nothing (every test run, ungated) | if anything is flagged, **E4b is not dispatched** |
| S1: constant | series, 200 | H3 and H4 false-alarm rates | reported |
| S2: Gaussian random walk at E3's rate (4.3 × 10⁻⁹ rad²/s) | series, 200 | H3 alarms in ≤ 0.06 of trials | reported |
| S3: S2 plus a 0.005 rad step on all qubits at a random interval | series, 200 | H3 power ≥ 0.80; H4 rate reported | if power < 0.80, K or the cadence is revised **before** the freeze |
| S3_one: the same step on one qubit | series, 200 | H3 power | reported (agent addition) |
| S4: gain drift, 0.005 rad over the run plus the S2 walk, with Δδ = Δε/2 | series, 200 | the β interval covers 0.5 in ≥ 0.90 of trials | reported |
| S4_counts: S4 through the full count-level analysis | counts, 20 | H1 rate; coverage of 0.5; regions | reported (agent addition) |
| S5: X-only step, 0.005 rad, plus the S2 walk | series, 200 | the β interval covers 0 in ≥ 0.90 of trials | reported |
| S6: 0.9 MHz detuning step (32 ns X assumed) | series, 200 | H2F fails on every qubit | reported |

At series level, the estimates get their Fisher noise, with the δ and ε errors correlated. The whole stage ran in 177 s, inside the node's 300 s horizon.

### E4-sim result

Recorded run: `results/drift_checks/pulse_amplitude/simulation/oc.json`. It took 177 s, with code sha256 `431602d8df17…8955`. The earlier run with the rejected median rate is kept as `oc_rejected_median_rate.json`.

| scenario | result | expectation | verdict |
|---|---|---|---|
| S1, counts (the control) | H1 flagged 0 of 16 tests; H3 no jumps; H4 p = 0.51 | nothing flagged | **PASS**: dispatch allowed |
| S1, series | H3 alarms 0.025; H4 alarms 0.070 | reported | — |
| S2 | H3 alarms **0.070** (H4 0.045) | ≤ 0.06 | **FAIL** (consequence: reported) |
| S3 | H3 power 1.000; H4 1.000 | ≥ 0.80 | **PASS** |
| S3_one | H3 power 0.600 (H4 0.070) | reported | — |
| S4 | covers 0.5 in 0.915; all 200 read gain | ≥ 0.90 | **PASS** |
| S4_counts | H1 in 20 of 20; all 20 read gain; median β 0.452, median interval width 0.073; covers 0.5 in 0.35 | reported | — |
| S5 | covers 0 in 0.910; 182 X-only, 18 inconclusive | ≥ 0.90 | **PASS** |
| S6 | H2F detects it in 1.000 | reported | — |

**S2 failed its threshold.**
- **Uncertainty:** at 200 trials the binomial SE is about 0.015, so 0.070 cannot tell a calibrated 0.05 from a mildly liberal test. The development seed gave 0.048 over 400 trials.
- **Consequence:** under the approved table, S2 is reported and does not gate dispatch. It was not rerun on another seed or with more trials.
- **Decision:** the human accepted H3 as it is on 2026-10-27 ("accept S2, commit and push it").

**H4's false-alarm rate** under S1 (0.070) is above α by the same margin, with the same caveat.

**S4_counts.** At count level, the β interval is narrow (about 0.07) and centred near 0.45, so it seldom contains exactly 0.5. Every trial still lands inside the gain band. The pre-registered decision uses the bands, so this bias does not change a reading; it is listed under Limits.

## Development record

- **The series-level simulation missed an estimator failure.** The first version weighted each circuit by the pooled model's probability, and allowed contrast a > 1. With E3's residuals injected unscaled, `sx33`'s model probability fell to the clip (10⁻⁹). Its weight became about 10¹⁶, which froze every job's δ at the pooled value: β read 0 for a gain change. A series-level simulation cannot see this, because it adds noise straight to the parameters.
  - **Fixes:** weights from the smoothed observed frequencies, readout fitted as errors in [0, 0.5], and simulated residuals scaled by 4p(1−p).
  - **Coverage:** count-level tests on all four parameter sets, and S4_counts.
- **The misfit correction was dropped.** It subtracted the static residual before the per-job fits, but it changed per-job estimates by less than 10⁻⁴ (checked on all four parameter sets), because a static misfit biases every job alike.
- **A median-based walk rate for H3 was tried and rejected.** It was meant to stop one jump from inflating the fitted rate. In the first recorded E4-sim run it failed S2: under a pure walk, H3 false-alarmed in 0.205 of trials, against the approved 0.06 (`simulation/oc_rejected_median_rate.json`). At E3's walk rate the walk and the measurement noise are comparable, so the median estimate often fell to zero, and the null then had no walk.
  - **Development seed (777, 400 trials, not the recorded seed):**

    | estimator | S2 false alarms | S3 power | S3_one power |
    |---|---|---|---|
    | mean | 0.048 | 0.995 | 0.540 |
    | median | 0.223 | 0.975 | 0.555 |
    | leave-one-out mean | 0.052 | 0.995 | 0.540 |

  - **Decision:** H3 went back to the mean, the simplest calibrated choice; the median added no power. E4-sim was then rerun on its recorded seed.
- **Planted bugs:** 20 of 22 were caught by `tests/test_pulse_amplitude.py`. The two survivors are equivalent under the current design:
  - **Model-based weights:** they behave like the smoothed weights now that readout is constrained. The smoothed weights are kept as a safeguard.
  - **A null without the walk:** H3's statistic standardises each increment by the walk fitted to the same series, so the null's walk size barely matters. That is what keeps H3 calibrated, and what the median version broke.

## Agent choices to review

These were chosen by the agent during the build. They are not among the approved values.
- the control and simulation seeds (2026102802, 2026102803);
- 200 series-level trials and 20 count-level trials;
- the scenarios S3_one and S4_counts;
- the 4p(1−p) scaling of the simulated residuals;
- 32 ns as the X length in S6;
- 5 Gauss–Newton steps, and half-count smoothing of the weights;
- the extra-qubit rule's exclusion of neighbours of earlier extras;
- rounding the Ramsey delay to the backend's granularity;
- H4's decision on the correlation test alone.

## Analysis rules

These are as in E3.
- **One run.** Once frozen: no reruns, and no change of seeds, sizes, α, replica counts, grids, code, backend or qubits.
- **Everything is kept.** All output goes under `results/drift_checks/pulse_amplitude/`, and the ledger must verify.
- **No interim look.** During the run, only job status is read.

## Limits

- **No direct view of amplitudes.** Pulse amplitudes are not observable, so every conclusion rests on these signatures. A re-tune that scales X and √X together reads as gain in H2; H3 and the secondary checks then carry the distinction.
- **Queueing.** IBM's queue can break the cadence, as it did for E3's jobs 3–8. H3 standardises by Δt but needs at least 15 intervals of 300 s or less.
- **Power on one qubit.** H3's power for a single 0.005 rad step on one qubit, against E3's walk, is moderate (S3_one). A shared step is found far more reliably (S3).
- **Static misfit.** Changes cancel the misfit only if it is static. E4 repeats the residual analysis. Qubit 54 already fails this assumption in E3, so its results are reported with that caveat.
- **β bias.** With E3's residuals injected, a gain change reads as β of about 0.45 (0.43–0.49 near noiseless; median 0.452 in S4_counts), not 0.5. Its interval seldom contains 0.5, but stays inside the gain band.
- **Scope.** One backend and one afternoon. A finding here is about ibm_fez during E4's window.

## Operational directives and agent rules

- **1. Invariance:** there is no tunable parameter, so this is not applicable.
- **2. Noiseless baseline:** Aer without noise gives:
  - `prep0`: 0 ones; `prep1` and `ram0`: all ones;
  - `ampN`, `sxM` and `ram90`: about half;
  - `x100`: 0 errors.

  See `tests/test_pulse_amplitude.py::test_noiseless_baseline`.
- **3. Fail closed:** ruff and the bridge and osiris-governance suites must be clean before dispatch. E4b also stops on:
  - the constant control flagging anything;
  - a non-operational backend;
  - a low or unknown budget.
- **4. Ledger:** rows are written ahead and hash-chained. E2 and E3 are still unsigned; **an enrolled node key should sign E4's rows.**
- **5. Heron variance:** K = 30 ≥ 10.
- **GEMINI.md:**
  - E4a and E4b call IBM Quantum, so the human must start them.
  - E4b (about 75 minutes) is longer than the node horizon and needs explicit approval, as does its QPU budget (about 150 s). Neither has been given yet.

## Before freezing

1. ~~The human approves the proposed values~~ (done 2026-10-27).
2. ~~Implement `bridge/pulse_amplitude.py` and `tests/test_pulse_amplitude.py`~~ (done; 47 tests).
3. ~~Run E4-sim and record its operating characteristics~~ (done; see above; **S2 failed its threshold**).
4. ~~The human decides on S2~~ (accepted 2026-10-27). The human reviews the agent choices above.
5. Freeze: rename this file `E4_PREREGISTRATION.md`, record its SHA-256, stamp it with `ots`, and pin the code hashes:
   - `bridge/pulse_amplitude.py`;
   - `bridge/drift_check.py`;
   - `osiris_governance/telemetry.py`.
6. The human approves E4b's duration and budget, then starts E4a and E4b:
   `python -m bridge.pulse_amplitude hardware --key-file <key> --out results/drift_checks/pulse_amplitude --prereg results/drift_checks/E4_PREREGISTRATION.md --after results/drift_checks/amplitude/state.json`
