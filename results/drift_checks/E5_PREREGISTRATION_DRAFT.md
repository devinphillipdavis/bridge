# E5, shared gain drift: pre-registration DRAFT

Status: **draft, not frozen.** Nothing has been sent to IBM.
- **Values:** approved by the human on 2026-10-28 ("approve the proposed values and build it"). Values marked "approved" below were the proposals; one E5-sim value was changed at the build, for the reason given under "Build".
- **Built:** `bridge/gain_drift.py` and `tests/test_gain_drift.py` (see "Build"), then changed after E5-sim v1 (see "Changes after E5-sim v1"). Results are under "E5-sim v2 results"; v1's are kept.
- **Open before the freeze:** the human's review of E5-sim v2 (R_min = 5 mrad, every check passing), of Cochran's false-heterogeneity rate and of the agent's build choices; then setting `r_min_rad`.
- **Freezing means:** recording the SHA-256 of the final file (renamed `E5_PREREGISTRATION.md`), stamping it with OpenTimestamps, pinning the code by hash, and committing it before any E5 job.

## Why

E4's pre-registered reading was "step changes in a shared gain, or a joint x and sx re-tune; H2F falsifies amplitude on qubits [24, 82]" [OBSERVED]. Its post hoc checks suggested something narrower [HYPOTHESIS]:
- Qubit 82's step, the only H3 jump, was a separate event: frequency and ε moved together on one qubit.
- On the other qubits, the changes in ε were correlated, had a ratio Δδ/Δε near a gain change's 0.5 (0.508 without qubit 82), and showed no jumps.

**Hypothesis under test (shared gain drift):** on ibm_fez, the job-to-job changes in single-qubit rotation error are mostly a slow drift in an analog drive gain shared across qubits. A gain change g scales every drive pulse, so Δε = πg and Δδ = (π/2)g [DERIVED]. E4 found this pattern post hoc. E5 tests it on fresh data (directive 8).

### Groundwork on E4's committed data

This was done on 2026-10-28 while designing E5. It is exploratory and post hoc; E5 does not reuse E4's data. The scripts were scratch files and are not part of E5's code.

| | finding | consequence for E5 |
|---|---|---|
| G1 | Between the medians of E4's first and last 10 jobs (about 15 minutes apart), ε rose by 1.41–2.03 mrad on 6 of the 7 qubits other than 82 (4: 1.90, 22: 1.78, 54: 2.03, 24: 1.87, 14: 1.56, 71: 1.41), and by 0.72 mrad on 144. Their pooled ε ranged from 2 to 12 mrad. | A gain shared by every qubit predicts **equal changes in radians**, whatever each qubit's starting ε. A per-qubit re-tune predicts changes toward zero, proportional to each qubit's ε. E5 tests equality directly (H3). |
| G2 | At E4's median spacing of 18 s, the ε increments correlate weakly across qubits (mean pairwise 0.11), but 5-job running means correlate at 0.90. On 6 of the 7 qubits, each qubit's increments vary around the cross-qubit mean by 0.6–1.4 times measurement noise (qubit 4: 2.7). The cross-qubit mean increment varies by 2.3 times its noise. | At 18 s spacing, measurement noise dominates the increments. E5 spaces jobs at least 600 s apart, so a drift at E3's rate (4.3 × 10⁻⁹ rad²/s) gives increments of about 1.6 mrad, against an increment noise of about 0.76 mrad. |
| G3 | Without qubit 82, E4's per-qubit β values are homogeneous: pooled 0.508 (SE 0.066), Cochran's Q = 8.0 on 5 df (p = 0.16). | β's homogeneity across qubits can be tested; a common gain predicts it, and independent re-tunes of X and √X do not. |
| G4 | In a noiseless check of E4's estimator, a contrast drop from 0.99 to 0.89 (as on qubit 4) shifts ε by −0.14 to −0.21 mrad under E3's static misfit, and by 0 without it. Qubit 4's ε rose by 1.6 mrad, in the opposite direction. | Readout leaks into ε only slightly. To be safe, qubits whose contrast varies are left out of E5's primary set. |
| G5 | Billing: 5 s per 12,000-shot job (E2, E4), 8 s per 22,000-shot job (E3). Queue waits had a median under 1 minute in E2 and E3 (14:47–17:33Z), but 38–91 minutes in E4 (22:45–00:37Z). Queueing packed E4's 30 jobs into 20 minutes. | E5 keeps **one job in flight**: each job is submitted only after the previous one has finished. Execution times are the time variable. Starting in UTC daytime is suggested but cannot be relied on. |
| G6 | E4's frequency test (H2F) flagged qubit 24 on a rank correlation of *levels*, although its frequency moved by at most 56 kHz. A slow trend shared by both series produces such a correlation. | E5 correlates **increments** of frequency and ε. |
| G7 | E4's shared test (H4) subtracts each qubit's mean increment, so a common linear trend is invisible to it [DERIVED]. | E5's sharing test (H2) uses uncentred increments, with a sign-flip and circular-shift null. |

## Question

Are the changes in ε a gain drift shared across qubits? If they are shared, is the gain fully common (equal changes in radians on every qubit), and does it drift smoothly or step?

## What each explanation predicts

| observable | **fully common gain drift** (hypothesis) | partly shared gain | independent per-qubit gain drift | per-qubit re-tunes of X and √X | frequency (detuning, TLS) |
|---|---|---|---|---|---|
| β = Δδ/Δε, per qubit | ≈ 0.5 on every qubit | ≈ 0.5 | ≈ 0.5 | differs between qubits and between re-tunes | ≈ 2 |
| sharing of ε increments across qubits | strong | present | none beyond chance | only if re-tunes coincide | none |
| increments equal in radians across qubits | **yes** | no | no | no; each moves toward its own zero | no |
| common component over time | smooth | — | — | steps | — |
| Ramsey frequency increments vs ε increments | uncorrelated | uncorrelated | uncorrelated | uncorrelated | correlated |
| net change against starting ε, across qubits | slope 0 | slope 0 | slope 0 | slope < 0 | — |

## Design

### E5: hardware

- **Backend:** ibm_fez. If it is not operational at start, E5 is not run.
- **Qubits (approved: 16):** E4's eight (4, 22, 54, 82, 144, 24, 14, 71), plus eight more chosen by E4's `choose_extras` rule: lowest reported readout error at start, skipping E4's qubits, their neighbours and neighbours of extras already chosen.
  - More qubits make the sharing and equality tests stronger, and parallel qubits cost no QPU time.
- **Circuits:** E4's 12, unchanged (`pulse_amplitude.pulse_circuits`): `prep0`, `prep1`, `amp0/16/32/48`, `sx17/33/49`, `ram0`, `ram90`, `x100`. Ramsey delay 200 ns, rounded to the backend's granularity.
- **Shots (approved):** 1000 per circuit, as in E4.
- **Jobs (approved):** K = **24**.
- **Schedule (approved):** **one job in flight.** Job k+1 is submitted at the later of two times: when job k finishes, and **600 s** after job k was submitted. Without queueing, jobs run about 600 s apart, so the run takes about 4 hours.
  - **Wall-clock cap:** no submission later than **8 hours** after the first.
  - **Minimum:** at least **16** jobs must finish, and at least 12 intervals between executions must be 300 s or longer. Otherwise the run is reported as incomplete.
- **Calibration records:**
  - a `calibration-snapshot` row from freshly refreshed properties before each `submit-intent`, as in E4;
  - **new:** a second snapshot as soon as each job's result arrives;
  - **new (E5a, after the run, 0 QPU-seconds):** `properties(datetime=t)` at each job's execution midpoint, through E4a's code. This fixes E4's gap, where every snapshot predated execution.
- **Resume:** if the process dies, a `resume` command continues from `state.json` and the ledger. A `submit-intent` without a matching `submitted` row is recorded as lost, never silently resubmitted.

**Budget:**
- **Expected cost:** about 5 s per job, so about **120 s**.
- **Available:** about 248 s should remain of the open plan's 600 s per 28 days. E2's 100 s should leave the window around 2026-11-24.
- **Fail closed (approved):** E5 is not started if the account reports fewer than **150 s** or no figure. Submission stops if a check before a job reports fewer than **20 s**.

**Precision per job [DERIVED, from E4's Fisher values]:**
- ε: 0.00054 rad, so a difference between two jobs has an SD of about 0.00076 rad;
- δ: 0.00053–0.00060 rad;
- Ramsey frequency: about 25 kHz.

## Estimation (per qubit)

Steps 2–4 are E4's, reused unchanged from `bridge.pulse_amplitude` (pinned by hash). Step 1 adds a joint fit of γ, approved on 2026-10-28 (see "Changes after E5-sim v1"):
1. E4's pooled maximum-likelihood fit of (e0, e1, γ, δ, ε). γ is then **refitted jointly with the per-job values**: it minimises the weighted squared residual of step 2's per-job fit (E4's inverse binomial variances as weights), by a bounded scalar search over E4's γ bounds. γ = 0 is kept if it fits at least as well. The pooled fit's γ is reported.
2. per-job (c_k, a_k, δ_k, ε_k) by 5 weighted Gauss–Newton steps, at the joint γ;
3. Ramsey frequency shift Δf_k;
4. **null replicas:** S = **2000** (approved) binomial draws per qubit at the run's mean frequencies, passed through steps 2–3 at the same joint γ. Analysis seed (approved): 2026102901.

New in E5:
5. **Increments:** Δε_{q,k} = ε_{q,k+1} − ε_{q,k}, Δδ likewise, and ΔΔf likewise. Their null variances v_{q,k} come from the same replicas.
6. **Common component:** for each interval, the precision-weighted mean of Δε_{q,k} over the primary set P (defined below), and its cumulative sum, the common series ε̄_k.

## Confirmatory tests

Every test uses S draws and α = **0.05** (approved). The hypotheses are evaluated in a fixed sequence. H1–H4 are evaluated only if the gate holds.

**Exclusions (H5 and HC, computed first on all 16 qubits):**
- **H5 (frequency):** a qubit fails if max|Δf_k| ≥ **200 kHz**, or if the Spearman correlation of ΔΔf with Δε is flagged. The correlation test is two-sided, with a permutation null and Holm's correction over 16 qubits.
- **HC (contrast):** a_k dispersion, with Holm's correction jointly with the gate's ε tests (32 tests).
- **Primary set P:** the qubits that fail neither H5 nor HC. All tests below use P; results on all 16 qubits are reported as secondary.

| | hypothesis | test |
|---|---|---|
| **Gate** | ε varies | ε_k dispersion (`osiris_governance.telemetry.dispersion`, as in E4) flagged, after the joint Holm correction, on at least **3** qubits of P. If not: "no variation in E5's window", which is not evidence of absence. |
| **H1** (gain ratio) | β = 0.5 on every qubit | Per-qubit Deming β of δ_k on ε_k, as in E4's H2, pooled over P by inverse pairs-bootstrap variance, with a 95% percentile interval. **Gain** if the interval lies inside **[0.30, 0.75]** (approved). **Homogeneous** unless Cochran's Q (inverse bootstrap variances, χ² on \|P\|−1 df) is flagged at α. |
| **H2** (shared) | the ε increments are shared | Mean pairwise *uncentred* correlation of Δε across P. Null: each qubit's increment series is independently sign-flipped (±1) and circularly shifted. The null holds if qubits are independent and none has a preferred direction. Flagged at α. |
| **H3** (fully common) | equal increments in radians | **Loadings test** (made confirmatory on 2026-10-28; see "Changes after E5-sim v1"). Each qubit's deviation from the common increment, projected on the common increment: L = Σ_{q∈P} w_q [Σ_k (Δε_{q,k} − Δε̄_k) Δε̄_k]² / Σ_k Δε̄_k², with w_q = 1/v_q. Null: the replicas' deviations from their own common increment, projected on the observed Δε̄_k, i.e. "fully common, along the observed path". **Rejected** (not fully common) if flagged at α. Loadings κ_q (each qubit's increments regressed on the leave-one-out mean) are reported. |
| **H4** (smooth) | the common component has no steps | E4's H3 jump test (Gaussian walk with drift, max\|z\| statistic, parametric-bootstrap null), applied to the common series ε̄_k, with variances from the replicas. Flagged at α. |

**Power condition for H3.** "Not rejected" counts as support only if the common series moved by at least **R_min** over the run, measured as the range of ε̄_k. R_min is the smallest common range at which E5-sim gives H3 a power of at least 0.80 against loadings spread uniformly over [0.5, 1.5]. It is set from E5-sim before the freeze. (This draft first expected about 2 mrad; that was wrong, see "Build".) Below R_min, H3 is reported as **underpowered**.

**Secondary (reported, not judged):**
- Per-qubit jump tests (E4's H3, with Holm), naming single-qubit events.
- The net change in ε against starting ε across P: least-squares slope with a bootstrap interval. Re-tunes predict a negative slope; common gain predicts 0.
- The loadings κ_q, and the correlation of increment residuals against distance on the coupling map.
- H3 as first drafted: the spread statistic T = Σ_k Σ_{q∈P} (Δε_{q,k} − Δε̄_k)² / v_q against its replica distribution.
- Calibration snapshots (before and after each job) and E5a entries dated inside any interval where the common series moves by more than 3 of its SDs.
- E4's residual pattern per circuit (the G2 check), repeated.
- All tests repeated on all 16 qubits.

## Reading the outcome

This will be `gain_drift.reading`, applied in order:

| condition | reading |
|---|---|
| gate FALSE | not evaluable: no variation in E5's window |
| H5 fails on more than half the 16 qubits | **frequency, not gain** |
| H1 interval outside [0.30, 0.75] | **not a gain change**: the region is named (X-only if the interval contains 0 and ends below 0.25; detuning-like if it starts above 1.25; otherwise inconclusive) |
| H1 gain, not homogeneous | gain-like on average, but β differs between qubits: re-tunes or mixed mechanisms |
| H1 gain and homogeneous, H2 not flagged | **independent per-qubit gain drift**: falsifies sharing |
| H1, H2 flagged, H3 rejected | **partly shared gain**: a common component plus qubit-specific changes |
| H1, H2 flagged, H3 underpowered | shared gain; whether it is fully common is untested |
| H1, H2 flagged, H3 not rejected, H4 not flagged | **fully common gain drift: supports the hypothesis** |
| H1, H2 flagged, H3 not rejected, H4 flagged | **fully common gain, changing in steps**: a shared change, not a drift |

Qubits excluded by H5 or HC are named in the reading.

## Simulation stage before dispatch (E5-sim)

`python -m bridge.gain_drift simulate --out results/drift_checks/gain_drift/simulation --scenario S<n> [--trials a:b] [--range R --path walk|linear]`, one invocation per scenario, S3 grid point and block of trials (each under 300 s), then `python -m bridge.gain_drift summarize --out results/drift_checks/gain_drift/simulation`, which merges them, sets R_min and writes `oc.json`.

**Synthetic qubits:** 16, cycling E3's four parameter sets, with E3's residuals scaled by 4p(1−p), as in E4-sim. There are K = 24 jobs, 600 s apart. **Seeds (approved):** control 2026102902, scenarios 2026102903.

**Lesson from E4:** a series-level simulation missed an estimator failure, so every scenario that bears on a confirmatory decision runs at **count level** (counts → estimator → tests).

| scenario | trials (approved) | expectation | consequence |
|---|---|---|---|
| S0: constant (the software control) | 1 | the gate, H1–H5 and HC flag nothing (all run, ungated) | if anything is flagged, **E5 is not dispatched** |
| S1: fully common gain walk at E3's rate, δ moving by half of ε | 50 | H2 power ≥ 0.80; H3 false rejections ≤ 0.10; H1 reads gain in ≥ 0.90; H4 false alarms ≤ 0.10 | if H2 power < 0.80, K or the spacing is revised **before** the freeze |
| S2: independent gain walks per qubit, at the same rate | 50 | H2 false alarms ≤ 0.10 | reported |
| S3: shared gain with loadings drawn uniformly from [0.5, 1.5] | 50 per common range, over a grid | H3 power as a function of the common range; this sets R_min | sets R_min |
| S4: per-qubit re-tunes: 1–2 steps per qubit toward zero, each with β drawn from [0, 1] | 50 | the reading is not "fully common gain" in ≥ 0.90 | reported |
| S5: S1 plus a qubit-82-like event on one qubit (a 5 mrad ε step together with a 340 kHz frequency step) | 50 | H5 excludes that qubit in ≥ 0.80; the reading stays "fully common" in ≥ 0.80 | reported |
| S6: S1 plus a qubit-4-like contrast drop on one qubit | 50 | HC excludes it in ≥ 0.80; the reading is unchanged in ≥ 0.80 | reported |
| S7: one common gain step of 2 mrad at a random interval | 50 | H4 power ≥ 0.80 | reported |

The bounds of 0.10 on false alarms allow for the binomial error of 50 trials (SE about 0.03 at 0.05). To fit inside the 300 s node horizon, E5-sim may run in several invocations, each under 300 s, recorded separately. Replica counts inside E5-sim may be reduced (approved: 500; changed at the build to 1000, see "Build"), and the reduction is reported.

## Analysis rules

As in E3 and E4:
- **One run.** Once frozen: no reruns, and no change of seeds, sizes, α, replicas, bands, code, backend or qubits.
- **Everything is kept.** All output goes under `results/drift_checks/gain_drift/`, and the ledger must verify.
- **No interim look.** During the run, only job status is read: the submission log and the ledger's row kinds.

## Limits

- **No direct view of amplitudes.** Pulse amplitudes are not observable. "Gain" means the signature Δδ = Δε/2, shared across qubits.
- **The physical source stays open.** Even a fully common gain does not identify its source (a shared reference, an amplifier, or the temperature of the control electronics). E5 tests the signature only.
- **Scope.** One backend for a few hours. A daily cycle, if there is one, would show up only in part.
- **β bias.** Under E3's injected misfit, the estimator reads a true 0.5 as about 0.45 (E4-sim). The approved band contains both values.
- **Absence of sharing.** H2 not flagged is evidence against sharing only as far as E5-sim's power at the observed drift size allows. E5-sim's S1 and S2 results are reported next to the outcome.
- **The sign-flip null** assumes that, without a shared cause, qubits have no preferred direction of drift. A common physical cause with qubit-specific sizes would count as sharing; that is intended.
- **Queueing.** One job in flight keeps the schedule, but a long queue stretches the run toward the 8-hour cap.
- **Run length on the phone.** About 4 hours or more under PRoot. The `resume` command covers a killed process; a sleeping device delays submissions but does not lose jobs.

## Operational directives and agent rules

- **1. Invariance:** there is no tunable parameter, so this is not applicable.
- **2. Noiseless baseline:** E4's circuits, unchanged; see `tests/test_pulse_amplitude.py::test_noiseless_baseline`.
- **3. Fail closed:** ruff and the bridge and osiris-governance suites must be clean before dispatch. E5 also stops on:
  - the constant control flagging anything;
  - a non-operational backend;
  - a low or unknown budget.
- **4. Ledger:** rows are written ahead and hash-chained. E2–E4 are still unsigned; **an enrolled node key should sign E5's rows.**
- **5. Heron variance:** K = 24 ≥ 10.
- **8. Fresh re-measurement:** E5 tests on new data a pattern found post hoc in E4. E4's data is used only for this groundwork.
- **GEMINI.md:**
  - E5 and E5a call IBM Quantum, so the human must start them.
  - E5 (about 4 hours, capped at 8) is longer than the node horizon and needs explicit approval, as does its QPU budget (about 120 s, about half of what remains).
  - The bands, thresholds and sizes above were approved on 2026-10-28. R_min and the agent's build choices still need review.

## Considered, not included

- **A programmed gain change as a positive control.** Heron's fractional `rx(θ)` gates could inject a known gain change (RX(π − η) and RX(π/2 − η/2)) and measure β's bias on hardware. Whether ibm_fez exposes them on the open plan, at which angles, and how they are calibrated against X and √X is unverified. Finding out needs an IBM query (0 QPU-seconds) that the human would start. It can be added before the freeze if wanted.
- **More shots per job.** At about 5 s per job, the budget favours more jobs over more shots.
- **New circuits.** E4's circuits and estimator are tested and pinned. Changing them would need the estimator re-validated.

## Values (approved 2026-10-28)

| item | approved |
|---|---|
| qubits | 16: E4's eight plus eight by `choose_extras` |
| jobs, schedule | K = 24; one in flight; at least 600 s between submissions; cap 8 h |
| minimum for analysis | 16 jobs, and 12 intervals of 300 s or longer |
| shots | 1000 per circuit |
| replicas, α | S = 2000; α = 0.05 |
| gate | ε dispersion flagged on at least 3 qubits of P |
| H1 band | [0.30, 0.75], with Cochran's Q homogeneity at α |
| H5 frequency limit | 200 kHz on max\|Δf\|, plus the increments' Spearman correlation |
| H3 power condition | R_min from E5-sim (power ≥ 0.80 against loadings in [0.5, 1.5]) |
| budget | start ≥ 150 s; stop below 20 s |
| seeds | analysis 2026102901; control 2026102902; scenarios 2026102903 |
| E5-sim | the scenarios and thresholds above; count level; 50 trials; 500 replicas (changed at the build to 1000: see "Build") |

## Build (2026-10-28)

**Tests.** `tests/test_gain_drift.py` has 53 tests, all offline, with fake IBM objects and a fake clock. They cover:
- the new statistics, checked against their formulas;
- the reading table;
- the scheduler: one job in flight, spacing, queue delay, wall-clock cap, budget stop, unknown budget, retries before and after the intent, and snapshots on both sides of each job;
- `resume`, including refusing changed code;
- `finish`'s minimums;
- E5a;
- the E5-sim scenarios and summary.

**Planted bugs:** 37, all caught. The first pass caught 30 of 36. Five tests were then added for the six survivors (Cochran's degrees of freedom; H3's weighting and p-value; the loadings null; HC exclusion; Holm correction). A 37th bug, the earlier loadings null, was also caught.

`bridge/gain_drift.py` reuses `bridge.pulse_amplitude` unchanged: its circuits, pooled fit, per-job estimator, replicas, Deming ratio test, jump test, snapshots and calibration-date matching. `bridge/pulse_amplitude.py` was not edited, so E4's pinned hash still holds.

### Changed at the build

- **E5-sim replicas: 500 → 1000.** With 500 replicas the smallest p-value is 1/501 = 0.0020. That is above Holm's first bar over the 32 dispersion tests, 0.05/32 = 0.0016, so no ε dispersion could ever be flagged and the gate could never pass in E5-sim. 1000 is the smallest round count that clears the bar (1/1001 = 0.0010). The hardware analysis keeps S = 2000.
- **R_min's expected size: about 2 mrad was wrong.** A series-level check at the build (Gaussian noise of 0.54 mrad per job, 16 qubits, loadings uniform on [0.5, 1.5]) gave H3 a power of 0.80 at a common range of about **10 mrad** when the common component is a random walk, and about **32 mrad** when it is a straight line. A straight line is the least favourable path for a given range. E5-sim measures this at count level (S3); R_min is read from the walk grid.

### Agent's choices (for review)

Analysis:
- **Increment variances.** v_q, each qubit's increment variance among the replicas, is pooled over intervals: the replica jobs are identically distributed, so one variance per qubit is the efficient estimate of the v_{q,k} named above.
- **Common series.** Weights are 1/v_q. H4's common series is the weighted mean of the ε levels, less its first value; its per-job variance comes from the same weighted mean of the replicas.
- **Random streams.** Each qubit's replicas get their own stream from the analysis seed, and the tests get one more, so results do not depend on qubit order.
- **Secondary checks** (reported, not judged):
  - *Net change:* the mean of the last quarter of jobs (6) less the mean of the first quarter, against the first quarter's mean; bootstrap over qubits.
  - *Large common moves:* intervals where the common increment lies more than 3 robust SDs (1.4826 × MAD) from the median increment.
  - *Distance:* the Spearman correlation, over pairs in P, between residual-increment correlation and coupling-map distance.
  - *All 16 qubits:* reported as the reading's category only.
- **Loadings test.** It asks H3's question of the loadings alone: each qubit's deviation from the common increment is projected on the common increment, and the null is replicas carrying the observed common increments. At the build's series-level check it reached 0.80 power at about 6 mrad (walk) and 12 mrad (line), against the spread test's 10 and 32. It was built as a secondary test; the human made it the confirmatory H3 after E5-sim v1 (see "Changes after E5-sim v1").

Schedule:
- **Waiting and retries.**
  - Poll every 30 s.
  - Once submissions stop, wait up to **2 hours** for the job in flight; `resume` or `finish` collects anything later.
  - An unknown budget figure or a failed calibration snapshot submits nothing, and the next poll tries again. The run is not ended for a network blip.
- **Failed submissions.**
  - A failure before the `submit-intent` row retries the same job after the spacing.
  - A failure after that row records the job as lost and moves on, because a job may exist at IBM without its ID.
- **Polling errors.** Errors while collecting results or taking after-snapshots are retried at the next poll. An after-snapshot taken late is still recorded, with its own timestamp.
- **Queue-delayed jobs.** Under the approved rule (the later of "previous job finished" and "600 s after the previous submission"), a job delayed in IBM's queue is followed at once by the next. Those two executions can be close together; the minimum of 12 intervals of 300 s or longer covers this.
- **Resume safeguards.** `resume` refuses to run if the code hashes differ from the run's `preregistration` row, or if the circuits differ from its `design` row.
- **Software control size.** The control runs on as many synthetic qubits as the run has (16) and at the hardware's 2000 replicas.

E5-sim scenarios:
- S1, S3, S5 and S6 share each trial's common walk, and each qubit's counts have their own stream. So S5's and S6's trial n differs from S1's trial n only on the event qubit. "The reading stays fully common" (S5) and "the reading is unchanged" (S6) are both measured as the same reading category as S1's trial n. The fraction reading "fully common" is reported as well.
- S3's paths:
  - *walk:* the common walk rescaled to the range; this sets R_min;
  - *straight line:* reported.
  - Grids: in v1, 6, 8, 10, 12 and 16 mrad (walk) and 16 and 32 mrad (line), chosen from the series-level check; in v2, 3, 4, 5, 6 and 8 mrad (walk) and 8, 12 and 16 mrad (line), chosen from the loadings test's power in v1.
- S4: each qubit gets 1 or 2 re-tunes at random intervals. Each re-tune removes 50–100% of the current ε and moves δ by β times that, with β uniform on [0, 1].
- S5: ε (+5 mrad) and frequency (+340 kHz) step together on one qubit, as drafted; δ does not move.
- S6: from a random job onward, one qubit's |1⟩ readout error rises by an amount drawn uniformly from [0.02, 0.10] in each job, as on qubit 4 in E4.
- S7: the 2 mrad common step has no background walk, as drafted.

### Found while building

- **A conservative null in the loadings test, fixed.** The first version's null added the observed common increments to the replicas. Those replicas kept their own mean as well, which roughened the common direction. Under a smooth common path with equal loadings, its mean p was 0.75 instead of 0.5. The null now projects the replicas' deviations on the observed common increments. Over 40 seeds its mean p is 0.51 (straight line), 0.53 (walk) and 0.51 (no path), with rejections at 0.05, 0.025 and 0.05. A test now holds it to this.

- **The sign-flip null and re-tunes.** Re-tunes toward zero move every qubit in the same direction when all qubits start with positive ε, as E4's did. A smoke trial of S4 then flagged H2 (p = 0.002), as the "Limits" section anticipates. The reading still separated it: H1's homogeneity failed (Cochran p ≈ 10⁻⁹³), giving "gain-like on average, but β differs". S4's rate measures how often this holds.
- **A straight-line drift and H2.** For a straight-line common drift, circular shifts change nothing, so H2's power comes from the sign flips alone. The p-value then cannot fall below 2^(1−|P|), which is 3 × 10⁻⁵ for |P| = 16 and 0.25 for |P| = 3.

## E5-sim v1 results (2026-10-28, superseded)

v1 used H3 as first drafted (the spread test) and E4's pooled γ. Its outputs were moved to `results/drift_checks/gain_drift/simulation_v1/` before v2 and are kept.

**Run.**
- Code `bridge/gain_drift.py` sha256 `5b7efb7de942680c018edf21582f0dfd863c5d110ed3574e27deb59409624a37`, the same in every output file.
- 27 invocations: S0 at 2000 replicas, and 26 blocks of 25 trials at 1000 replicas.
- Six ran in parallel, each taking 157–292 s, under the 300 s horizon but close to it; next time, fewer at once.
- 6677 s of compute in all. Merged by `summarize` into `oc.json`.

| check | expected | result | |
|---|---|---|---|
| S0: constant control flags nothing | dispatch allowed | no flags (H1 inconclusive, Cochran p 1.00; H2 p 0.40; H3 p 0.52; H4 p 0.37) | PASS |
| S1: H2 power | ≥ 0.80 | 1.00 | PASS |
| S1: H3 false rejections | ≤ 0.10 | **0.16** | **FAIL** |
| S1: H1 reads gain | ≥ 0.90 | 1.00 | PASS |
| S1: H4 false alarms | ≤ 0.10 | 0.02 | PASS |
| S2: H2 false alarms | ≤ 0.10 | 0.04 | PASS |
| S4: reading not "fully common" | ≥ 0.90 | 1.00 | PASS |
| S5: H5 excludes the event qubit | ≥ 0.80 | 0.98 | PASS |
| S5: reading unchanged from S1's trial | ≥ 0.80 | 0.94 | PASS |
| S6: HC excludes the event qubit | ≥ 0.80 | 1.00 | PASS |
| S6: reading unchanged from S1's trial | ≥ 0.80 | 0.96 | PASS |
| S7: H4 power | ≥ 0.80 | 1.00 | PASS |

**R_min = 10 mrad** (agent-computed, by the approved rule). H3's power against loadings uniform on [0.5, 1.5]:

| common range | 6 mrad | 8 | 10 | 12 | 16 | 32 |
|---|---|---|---|---|---|---|
| H3, walk path | 0.42 | 0.78 | **0.84** | 0.98 | 1.00 | — |
| H3, straight line | — | — | — | — | 0.30 | 0.46 |
| loadings test (secondary), walk | 0.94 | 1.00 | 1.00 | 1.00 | 1.00 | — |
| loadings test (secondary), straight line | — | — | — | — | 1.00 | 1.00 |

The loadings test rejected a fully common gain in 0.02 of S1's trials, and in 0.00 of S5's and S6's.

**Readings when the hypothesis is true (S1, 50 trials):**

| reading | trials |
|---|---|
| fully common gain drift | 17 |
| shared gain, untested (common range below R_min) | 15 |
| gain-like, β differs | 13 |
| partly shared | 5 |

So E5 as built would support a truly fully common gain about a third of the time. In 36% of trials it would wrongly read β heterogeneity or qubit-specific changes.

Two things drive the misreadings:
- H3's false rejections (0.16).
- **Cochran's homogeneity test rejected in 0.26 of S1's trials.** No bound was set for it, but it comes first in the reading table. It rejected in 0.32 of S2's trials, in 0.72 of S3's 16 mrad straight-line trials, and in all of the 32 mrad trials.

**Cause [DERIVED, noise-free check of E4's estimator at 10⁹ shots, exploratory]:** the per-qubit estimates respond to a common gain change in slightly different ways.
- **E3's static misfit,** which the simulation injects, gives each parameter set its own response. Over a 10 mrad straight line, the loadings are 0.997–1.030 and β is 0.429–0.479. Without the misfit, the loadings are 1.002–1.003 and β is 0.482–0.500.
- **The pooled fit absorbs the drift into the decay rate γ.** The per-job fit then holds that γ fixed, which biases each job's ε. Over a straight line, the largest error is:

  | common range | 5 mrad | 10 mrad | 16 mrad | 32 mrad |
  |---|---|---|---|---|
  | largest ε error | 3–9 µrad | 20–52 µrad | 84–201 µrad | 0.72–1.20 mrad |

  With γ held at 0, the error is 0 at every range. On hardware, E4's pooled γ was at its bound of 0 on all eight qubits, so whether this happens on ibm_fez is open.
- **When β and the loadings are measured this precisely,** these offsets of a few percent register as heterogeneity.

**Other observations:**
- S2 read "independent per-qubit gain drift" in 34 of 50 trials and "β differs" in 16.
- S4 read "not a gain change" in 11 and "β differs" in 39. H2 flagged S4 in 0.40 of trials, the co-directional re-tunes described under "Build".
- S7's gate held in 0.90 of trials. H1 read gain in only 0.54, since a single 2 mrad step gives one long lever arm.

**Approved consequence.** The only consequence set for an S1 failure is that K or the spacing is revised if H2's power is below 0.80. H2's power was 1.00, so nothing is triggered. The H3 failure is reported, as E4's S2 failure was. The agent offered four options; the human chose two of them (see the next section).

## Changes after E5-sim v1 (2026-10-28)

The human chose options 2 and 3 of the four offered ("do 2 and 3 and rerun the simulation"):
- **Option 2: H3 is the loadings test.** The spread test T (H3 as first drafted) is now a secondary. In v1, T reached 0.80 power only at 10 mrad (walk), stayed below 0.50 on both straight lines, and rejected in 0.16 of S1's trials. The loadings test reached 0.94 at 6 mrad (walk) and 1.00 on both straight lines, and rejected in 0.02 of S1's trials.
- **Option 3: γ is refitted jointly with the per-job values** (estimation step 1). The pooled fit takes a common drift for decay, and the per-job fit then biases ε by holding γ at that value.
  - Build check [DERIVED, noise-free, without the misfit; exploratory]: over a straight line, the largest per-job ε error falls from 20–52 µrad to at most 1.2 µrad at 10 mrad, and from 0.72–1.20 mrad to at most 35 µrad at 32 mrad.
  - A real decay is kept. Without drift, the tests recover γ = 0, 2 × 10⁻⁴ and 10⁻³. With a drift and no decay, the joint γ is below 10⁻⁶, where the pooled fit's is above 10⁻⁴.
  - **Not removed:** E3's static misfit, which the simulation injects. With it, the joint γ goes to 0, but each parameter set still responds to a common change in its own way. Over a 32 mrad line the ε error stays at 130–551 µrad and β at 0.43–0.48 [DERIVED, exploratory]. Option 4 (a tolerance for this) was not chosen.
- **Tests:** 58, all passing. New or changed: the joint fit's handling of drift and decay, `estimate` using the joint γ for the jobs and the replicas, and `summarize` refusing outputs from mixed code versions.
- **Planted bugs:** 43, all caught. Six are new for this change. The first pass caught 42; the survivor (replicas fitted at the pooled γ) led to a stronger `estimate` test.
- **E5-sim v2 grids:** S3 walk 3, 4, 5, 6 and 8 mrad; straight line 8, 12 and 16 mrad.

## E5-sim v2 results (2026-10-28)

**Run.**
- Code `bridge/gain_drift.py` sha256 `cff696133d8adef89853206822e2b2f3c1a6e1e2c27e723f64f1ec3a163e9b2e`, the same in every output file (`summarize` now refuses mixed versions).
- 130 invocations at 1000 replicas, plus S0 at 2000:
  - S1 in three blocks of 16–17 trials;
  - S2 in two blocks of 17 trials and four of 4;
  - everything else in blocks of 5.
- Each trial's draws depend only on its scenario and index, so the block sizes do not change the results.
- **Run time.**
  - Five 16–17-trial blocks reached the 295 s limit and wrote nothing, while up to seven processes were running. They were rerun in smaller blocks, together with blocks stopped by hand or by a session restart.
  - The 5-trial blocks, four at a time, took 30–139 s (median 42 s, about 8 s per trial; v1 took about 10 s per trial, six at a time). The phone's speed varied: some groups of four blocks took 2–3 times the median.
  - The completed blocks took 5738 s in all. The joint fit adds about 31 per-job fits per qubit, but per-trial time did not rise much.
- Outputs are in `results/drift_checks/gain_drift/simulation/`, merged by `summarize` into `oc.json`.

| check | expected | result | |
|---|---|---|---|
| S0: constant control flags nothing | dispatch allowed | no flags (H1 inconclusive, Cochran p 1.00; H2 p 0.40; H3 p 0.78; H4 p 0.37) | PASS |
| S1: H2 power | ≥ 0.80 | 1.00 | PASS |
| S1: H3 false rejections | ≤ 0.10 | 0.02 | PASS |
| S1: H1 reads gain | ≥ 0.90 | 1.00 | PASS |
| S1: H4 false alarms | ≤ 0.10 | 0.02 | PASS |
| S2: H2 false alarms | ≤ 0.10 | 0.04 | PASS |
| S4: reading not "fully common" | ≥ 0.90 | 1.00 | PASS |
| S5: H5 excludes the event qubit | ≥ 0.80 | 0.98 | PASS |
| S5: reading unchanged from S1's trial | ≥ 0.80 | 0.94 | PASS |
| S6: HC excludes the event qubit | ≥ 0.80 | 1.00 | PASS |
| S6: reading unchanged from S1's trial | ≥ 0.80 | 0.96 | PASS |
| S7: H4 power | ≥ 0.80 | 1.00 | PASS |

**R_min = 5 mrad** (agent-computed, by the approved rule: the walk grid). H3's power against loadings uniform on [0.5, 1.5]:

| common range | 3 mrad | 4 | 5 | 6 | 8 | 12 | 16 |
|---|---|---|---|---|---|---|---|
| H3 (loadings test), walk | 0.48 | 0.78 | **0.86** | 0.94 | 1.00 | — | — |
| H3, straight line | — | — | — | — | 0.54 | 1.00 | 1.00 |
| spread test (secondary), walk | 0.12 | 0.22 | 0.40 | 0.42 | 0.78 | — | — |
| spread test, straight line | — | — | — | — | 0.06 | 0.12 | 0.30 |

R_min holds for a walk, not for a straight line. A straight line's power is 0.54 at 8 mrad and 1.00 at 12 mrad. So a straight-line common drift of 5–12 mrad that is not fully common can be read as fully common more often than one time in five. The grid does not show where between 8 and 12 mrad the line reaches 0.80.

**Readings when the hypothesis is true (S1, 50 trials; v1 in brackets):**

| reading | trials |
|---|---|
| fully common gain drift | 35 (17) |
| fully common gain, changing in steps | 1 (0) |
| gain-like, β differs | 13 (13) |
| partly shared | 1 (5) |
| shared gain, untested (common range below R_min) | 0 (15) |

So E5 as now built would support a truly fully common gain in 0.70 of trials, against 0.34 in v1.

**What is left: Cochran's homogeneity test.** It rejected in 0.26 of S1's trials, as in v1, and this causes 13 of the 15 misreadings. No bound was set for it, but it comes first in the reading table. It also rejected in:
- 0.30 of S2's trials;
- 0.46 of the 12 mrad and 0.64 of the 16 mrad straight-line trials;
- 0.98 of S4's, where it should, since S4 draws β per qubit.

**What the joint γ changed.** Per-trial results moved: the spread test's p differs from v1 in 28 of S1's 50 trials, and Cochran's p by up to 0.7 decades. The rates at α did not move (spread test 0.16, Cochran 0.26 in both runs). This fits the build check: the β spread in S1 comes from E3's misfit, not from γ.

**Other observations:**
- S2 read "independent per-qubit gain drift" in 34 of 50 trials, "β differs" in 15 and "partly shared" in 1.
- S4 read "not a gain change" in 11 and "β differs" in 39. H2 flagged S4 in 0.40 of trials, the co-directional re-tunes described under "Build".
- S5 read "fully common" in 36 trials and S6 in 37, each with one "in steps" and the rest "β differs".
- S7: the gate held in 0.90 of trials and H1 read gain in 0.54. Its 2 mrad step is below R_min, so 25 trials read "untested". H3 rejected in 0.10 of S7's trials, though its step is fully common; no bound was set for this.

**Approved consequence.** Every check passes, so nothing is triggered.

## Before freezing

1. ~~The human approves the proposed values, or changes them.~~ Done 2026-10-28.
2. ~~Implement `bridge/gain_drift.py` and `tests/test_gain_drift.py`, catching planted bugs as for E4.~~ Done 2026-10-28 (see "Build").
3. ~~Run E5-sim, record its operating characteristics, and set R_min.~~ Done 2026-10-28. v1: R_min = 10 mrad, and S1's H3 check failed. After options 2 and 3, v2: R_min = 5 mrad, and every check passes (see "E5-sim v2 results").
4. The human reviews E5-sim v2, Cochran's false-heterogeneity rate and the agent's build choices, and decides whether to change anything before the freeze. A change to code or tests means rerunning the tests, the planted bugs and E5-sim. Then set `r_min_rad` in `E5Design`. It is still `None` in the built code, so H3 currently always reads as underpowered. Setting it changes `gain_drift.py`'s hash from the one E5-sim v2 ran under; the frozen file should record both, and that this one value is the only difference.
5. Freeze: rename this file `E5_PREREGISTRATION.md`, record its SHA-256, stamp it with `ots`, pin the code hashes, commit and push.
6. The human approves E5's duration and budget and starts it:
   `python -m bridge.gain_drift hardware --key-file <key> --out results/drift_checks/gain_drift --prereg results/drift_checks/E5_PREREGISTRATION.md --after results/drift_checks/pulse_amplitude/state.json`
   then, after it finishes, E5a:
   `python -m bridge.gain_drift retrospective --key-file <key> --out results/drift_checks/gain_drift/retrospective --run results/drift_checks/gain_drift`
