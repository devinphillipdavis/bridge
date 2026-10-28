# Drift checks: pre-registration

Status: **frozen 2026-10-27**, before E1 has run and before any hardware job has been submitted.

How the freeze is recorded:
- E1 writes the SHA-256 of this file into `simulation/oc.json`.
- E2 writes it as the first row of `hardware/ledger.jsonl`.
- `PREREGISTRATION.md.ots` holds an OpenTimestamps proof of the file, if the stamp succeeds.

None of these files was committed when the freeze was made. Committing is left to the human.

Runs so far were smoke and unit tests only. None used IBM hardware and none used the seeds below:
- ≤ 30 trials;
- fake backends;
- one preflight query of the account (backend list, calibration data and usage, with no job submitted).

## Question

Can posterior predictive checks on time-ordered telemetry tell calibration drift apart from other ways the model is wrong? And does a current IBM backend show variation between jobs, or drift, beyond shot noise over about 100 minutes?

## Claims under test

The claims come from an AI-generated summary pasted on 2026-10-27. It recommends four time-ordered checks:
- slow trends in residuals;
- lag-1 autocorrelation compared against replicas;
- step shifts between sessions;
- a CUSUM monitor that "flags unmitigated calibration drift" when it passes a "pre-declared threshold (h = 27.0)".

It also says that using the same posterior draws and operating conditions for observed and replicated data "ensures that any divergence in time structure is strictly due to model misspecification (e.g., unmodeled drift) rather than noise distribution shifts".

Until tested, each claim below is [UNSUPPORTED]:

| | claim |
|---|---|
| C1 | time-ordered predictive checks detect calibration drift |
| C2 | a divergence in time structure means drift, not a misspecified noise distribution |
| C3 | a CUSUM with h = 27.0 flags drift |

**How h is read.** The paste gives no units, allowance or run length for h. C3 is read as a CUSUM on model-standardised residuals with the textbook allowance k = 0.5. That reading is ours, not the paste's.

## Checks

All code is pinned by hash (see below).

**Data and null model.** A series is K error counts y_k out of N shots, in time order. The no-drift model M0 has one rate for the whole series:
- p ~ Beta(1 + Σy, 1 + Σ(N − y));
- y_k ~ Binomial(N, p) for every job, in the observed order.

**p-values.** Each check is an upper-tail predictive p-value against S replicas from M0:

  p = (1 + #{T(rep) ≥ T(obs)}) / (1 + S)

The statistics are computed exactly, by `osiris_governance.telemetry`:

| check | statistic | changes when the series is scaled? |
|---|---|---|
| `dispersion` | sample variance | yes (heterogeneity check) |
| `lag1_autocorrelation` | Σ d_t d_{t+1} / Σ d_t² | no |
| `trend_r2` | squared correlation of value with time | no |
| `cusum_peak` | max_k S_k² / Σ d², S_k the running sum of deviations | no |
| `session_eta2` | between-session SS / total SS | no |

Comparators, of the scale-dependent kind the paste recommends:
- `slope_sq`: the squared least-squares slope, as a predictive check.
- `cusum_h27` and `cusum_h5`: a tabular CUSUM on residuals standardised by the pooled binomial model, z_k = (y_k/N − p̂) / √(p̂(1 − p̂)/N), with allowance k = 0.5. It alarms when either one-sided sum passes h = 27 (the paste) or h = 5 (a common textbook value).

The alarm level is α = 0.05 for every predictive check.

## E1: simulation

### Design

- K = 20 jobs, N = 4000 shots.
- Base rates: 0.005 and 0.02.
- M = 1000 trials per cell, S = 400 replicas per trial.
- Seed 20261027. Each cell uses `numpy.random.default_rng([seed, rate index, scenario index])`.
- Times are the job indices; the sessions are the first and second halves.

Sizes are in σ_b, the binomial SD of one job's rate, √(p(1 − p)/N). Rates are clipped to [1e-9, 1 − 1e-9].

| scenario | per-job rate | sizes |
|---|---|---|
| null | p | — |
| linear | p + Δσ_b (k/(K−1) − ½) | Δ = 2, 4, 8 |
| step | p ∓ Δσ_b/2, first half / second half | Δ = 2, 4, 8 |
| walk | p + random walk with step SD τσ_b, centred | τ = 0.5, 1, 2 |
| jitter | p + independent N(0, (jσ_b)²) per job | j = 1, 2, 4 |

### Predictions

Each prediction must hold at **both** base rates. The bound is α + 2√(α(1 − α)/M) = 0.0638.

| | prediction | tests |
|---|---|---|
| P1 | under null, every predictive check (the five plus `slope_sq`) alarms in ≤ 0.0638 of trials | calibration |
| P2 | alarm rate ≥ 0.80 for: linear 8 → `trend_r2`, `cusum_peak`; step 8 → `cusum_peak`, `session_eta2`; walk 2 → `lag1_autocorrelation` | C1 |
| P3 | under jitter 4, each scale-free time-order check alarms in ≤ 0.0638, while `dispersion` alarms in ≥ 0.80 | C2, scale-free checks |
| P4 | under jitter 4, `slope_sq` alarms in > 0.0638 | C2, scale-dependent check |
| P5 | under linear 8, `cusum_h27` alarms in < 0.50 | C3 |

How the results bear on the claims:
- **C1** is supported if P2 holds.
- **C2** holds for the scale-free checks if P3 holds. It fails for scale-dependent checks if P4 holds, since jitter with no time order would then look like drift.
- **C3** is not supported at this design if P5 holds, and supported if P5 fails.

## E2: hardware

### Choices fixed at the start
- **Backend:** the least-busy operational IBM backend for this account (open plan) at start. At the 14:37Z preflight that was ibm_fez.
- **Layout:** the 4 qubits with the lowest reported readout error at start, chosen by the same rule. At preflight these were 4, 22, 54 and 82 (readout error 0.0024–0.0029).

### Circuits
| circuit | what it does | ideal outcome | counted as an error |
|---|---|---|---|
| `prep0` | measure | 0 | outcome 1 |
| `prep1` | X, then measure | 1 | outcome 0 |
| `x100` | 100 X gates, then measure | 0 | outcome 1 |

### Running the jobs
- 4000 shots per circuit; one job holds the three circuits.
- SamplerV2 with default options, optimisation level 0, trivial layout and no routing (`dnalang.backends.ibm.transpile_fixed`).
- K = 20 jobs, submitted 300 s apart. Before each submission, `dnalang.backends.ibm.submit` writes a `submit-intent` ledger row that includes the backend's calibration hash.

### Inputs to the checks
- **Times:** the midpoint of each job's execution span, as IBM reports it.
- **Sessions:** the calibration hash at submission.

### Control
The same 20 jobs on `AerSimulator.from_backend(backend)` at start, with seeds 2026102702 + j. The noise is static, so there is no drift by construction.

### Analysis
- 12 series (4 qubits × 3 circuits), each with 4 checks, or 5 when the run spans two or more calibration sessions.
- S = 4000 replicas, seed 2026102701.
- Holm's correction over every p-value in the run, at family-wise error rate 0.05.
- At least 10 jobs must finish, or the run is reported as incomplete.

### Budget
594 QPU-seconds remain in the 28-day window. Expected use is under 120 s.

### Hypotheses

| | hypothesis |
|---|---|
| HD | at least one series is flagged by `dispersion`: variation between jobs beyond the binomial model. This is the premise of operational directive 5. |
| HT | at least one series is flagged by a time-order check. |
| control | nothing is flagged. If the control flags anything, the pipeline has a false-alarm problem and the hardware flags are not interpreted. |

No direction is predicted, and both outcomes are reported. Our expectation, not a criterion: HD is more likely than HT, and for `prep1` and `x100` rather than `prep0`. The reason is that T1 fluctuations driven by two-level systems are documented on transmons (Klimov et al., PRL 121, 090502 (2018)).

How each outcome is read:
- **HD only:** jobs vary more than shot noise, but no time order was detected. Pooling shots across jobs then understates the uncertainty, which supports directive 5. It does not show drift.
- **HT:** the flagged series drifted over the run in a time-structured way.
- **Neither:** there is no evidence of anything beyond binomial noise at this scale (K = 20 over about 100 min). That is not evidence of absence.

## Analysis rules

- **One run each.** After this file is frozen:
  - no reruns;
  - no change of seeds, sizes, α, replica counts or code;
  - no swapping the backend or qubits.
- **Everything is kept.** All output is written under `results/drift_checks/` whatever the outcome, and the hardware ledger must verify.
- **Failed jobs** are recorded. The analysis uses the jobs that finished, if there are at least 10.

## Exploratory (reported, not judged)

- E1 power at sizes 2 and 4, and the `cusum_h5` rates.
- E2:
  - overdispersion factors (sample variance / binomial variance);
  - slopes per hour;
  - the readout error IBM reports against the measured `prep0` rate.

## Operational directives

- **1. Invariance:** the circuits have no tunable parameter, so ∂ℱ/∂θ is not defined. Not applicable.
- **2. Noiseless baseline:** Aer without noise gives zero errors on all three circuits (`tests/test_drift_check.py::test_noiseless_baseline_has_no_errors`). Any variation on hardware is therefore not circuit structure.
- **3. Fail closed:**
  - bridge: ruff clean, 45 tests pass;
  - osiris-governance: 271 pass, 3 skipped;
  - `git diff --check` clean in both.
- **4. Ledger:** rows are written ahead and hash-chained, but **unsigned**. An agent may not sign, so they need a signature from an enrolled node key. This is a deviation from directive 4.
- **5. Heron variance:** K = 20 ≥ 10 jobs. HD tests this directive's premise.

## Code and environment (pinned)

| file | sha256 |
|---|---|
| `bridge/drift_check.py` | `95b01edcf5ba54dce154236ba48cb386a1a073295b64f6568615b042219cbe60` |
| `osiris_governance/telemetry.py` | `f805f8bd9871f87c7a5463071a6b28b52036a8fd04e5044b7b8f5e51d378a420` |

**Repositories.** Both files carry uncommitted changes on top of these heads:
- bridge `2b4043d`
- osiris-governance `352e1cc`
- dnalang-core `d3c454e`

**Versions:** Python 3.13.5, qiskit 2.5.2, qiskit-aer 0.17.2, qiskit-ibm-runtime 0.50.0, numpy 2.5.3.

Commands (from `bridge/`):

```bash
python -m bridge.drift_check simulate --out results/drift_checks/simulation \
    --prereg results/drift_checks/PREREGISTRATION.md
python -m bridge.drift_check hardware --key-file <IBM key JSON> --out results/drift_checks/hardware \
    --prereg results/drift_checks/PREREGISTRATION.md
```
