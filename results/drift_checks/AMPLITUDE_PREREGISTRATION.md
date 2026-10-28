# E3, x100 amplitude drift: pre-registration

Status: **frozen 2026-10-27, before any E3 job**. Two things had already happened:
- one descriptive interim look at 6 of E2's 20 jobs (disclosed below);
- smoke and unit tests on model data and a noiseless simulator.

E2's own pre-registered analysis had not yet run.

How the freeze is recorded:
- the SHA-256 of this file is the first row of `amplitude/ledger.jsonl`;
- `AMPLITUDE_PREREGISTRATION.md.ots` holds an OpenTimestamps proof of the file.

The file was not committed when frozen.

## Why

E2's interim look (6 jobs, descriptive, not a test) found two things:
- The `x100` error counts varied between jobs 31–100 times more than shot noise allows, while readout varied only 0.8–3.9 times.
- On qubits 4, 54 and 82, the `x100` errors dropped between jobs 2 and 3 (qubit 54: 1740, 1758 → 1149, 1162, 1119, 1281). The backend's calibration hash did not change.

[HYPOTHESIS] The X gate's rotation error drifted.

## Question

Is the job-to-job variation in `x100` a change in the X gate's coherent rotation error ε, or a change in contrast a (readout or decoherence)?

## Circuits and model

Each job runs these circuits on every qubit in parallel, with 2000 shots each:

| circuit | gates | P(1) |
|---|---|---|
| `prep0` | measure | c − a/2 |
| `prep1` | X, measure | c + (a/2) cos ε |
| `ampN`, N ∈ {0, 4, 8, 16, 32, 48, 64, 100} | √X, then N X gates, measure | c + (a/2) e^(−γN) sin(δ + Nε) |
| `x100` | 100 X gates, measure | not fitted; kept for continuity with E2 |

Here each X rotates by π + ε and √X by π/2 + δ:
- c is the readout offset and a the contrast;
- γ is the contrast lost per X gate.

The two anchors are what make the fit possible:
- `prep0` and `prep1` pin c and a within every job.
- The √X prefix makes `ampN` move linearly with ε and shows its sign. `x100` alone gives only sin²(50ε).

Without the anchors, c, a, δ and ε cannot be told apart. On model data, a fit without `prep0` and `prep1` returned c = 0.43 and a = 1.08 for true values 0.50 and 0.95.

## Design

- **Backend and qubits:** ibm_fez, qubits 4, 22, 54 and 82, taken from E2's `state.json`, so the device and qubits are the same as E2's. If ibm_fez is not operational at start, E3 is not run.
- **Jobs:** K = 12, submitted 300 s apart through `dnalang.backends.ibm.submit`, which writes a ledger row before each submission. E3 starts after E2's last job has finished.
- **Budget:** about 9 QPU-seconds per job, so about 110 s in total. About 494 s should remain after E2's roughly 100 s.

## Estimation (per qubit)

1. **Pooled fit:** maximum-likelihood (c, a, γ, δ, ε) over all jobs, by L-BFGS-B from a grid start. The grid start takes c and a from `prep0`/`prep1`, and (δ, ε) by least squares over δ ∈ [−0.6, 0.6] in steps of 0.01 and ε ∈ [−0.08, 0.08] in steps of 0.0002. The bounds are c ∈ [0, 1], a ∈ [0, 1.2], γ ∈ [0, 0.01], δ ∈ [−1, 1] and ε ∈ [−0.1, 0.1].
2. **Per job:** (ε_k, a_k) jointly, by maximum likelihood, with c, γ and δ held at the pooled values. Holding none of the two fixed means a change in one can't be absorbed by the other.
   - ε is searched over the pooled ε ± 0.02 rad in steps of 1e-4; a over [0, 1.2] in steps of 0.005.
   - The search goes coarse grid first, then fine, and gives the same result as an exhaustive search (tested).
3. The estimates enter the exact checks in micro-units, as integers.

## Tests

For each qubit and each of the two parameters (ε_k and a_k), `osiris_governance.telemetry.drift_p_values` runs:
- dispersion;
- lag-1 autocorrelation;
- trend r²;
- CUSUM peak;
- session η², only if the jobs span two or more calibration hashes.

**Replicas.** There are S = 1000 per qubit, drawn by parametric bootstrap from the pooled fit, so they contain no drift. Each replica is refitted (pooled fit warm-started at the observed one) and re-estimated per job, exactly as the data are. The seed is 2026102703.

This is a plug-in bootstrap, not a full posterior predictive.

**Correction.** Holm's correction over every p-value in the run (32, or 40 with sessions), at family-wise error rate 0.05.

## Software control (before dispatch)

**Synthetic data.** Seed 2026102704. Four synthetic qubits:
- ε = 0.0068, 0.0038, 0.0125 and 0.0063, which spans what E2's interim `x100` rates imply;
- c = 0.5, a = 0.95, γ = 0.0005, δ = 0.02;
- 12 jobs of 2000 shots each.

The three controls and what they decide:

| control | expectation | consequence |
|---|---|---|
| constant | flags nothing | if it flags anything, **E3 is not dispatched** |
| ε step, +0.003 in the second half | HA true, HC false | reported, does not gate dispatch |
| contrast step, −0.1 in the second half | HC true, HA false | reported, does not gate dispatch |

## Hypotheses

| | hypothesis |
|---|---|
| HA | ε_k dispersion is flagged on at least one qubit: the rotation error varies between jobs |
| HAt | a time-order check on ε_k is flagged on at least one qubit |
| HC | a_k dispersion is flagged on at least one qubit: the contrast varies between jobs |
| HCt | a time-order check on a_k is flagged on at least one qubit |

**Prediction** (from the interim look): HA holds on at least one of qubits 4, 54 and 82.

| outcome | reading |
|---|---|
| HA and not HC | the `x100` variation is coherent rotation-error (pulse-amplitude) drift; supports the hypothesis |
| HC and not HA | it is contrast (readout or decoherence), not rotation error; falsifies the hypothesis |
| both | both vary |
| neither | no variation detected in E3's window, which runs roughly 1.5 h after E2's step; not evidence of absence |

## Exploratory (reported, not judged)

- **Link to `x100`:** per qubit, the Spearman correlation between the `x100` error counts and sin²(50 ε_k) across jobs. A positive correlation would mean ε accounts for `x100`.
- **ε in degrees,** for readability.
- **Comparison with E2:** E3's `x100` against E2's.
- **Signs and sizes** of the pooled δ and γ.

## Analysis rules

- **One run.** After this file is frozen:
  - no reruns;
  - no change of seeds, sizes, α, replica counts, grids or code;
  - no swapping the backend or qubits.
- **Everything is kept.** All output is written under `results/drift_checks/amplitude/` whatever the outcome, and the ledger must verify.
- **Minimum jobs.** At least 10 jobs must finish, or the run is reported as incomplete.

## Operational directives and agent rules

- **1. Invariance:** there is no tunable parameter, so this is not applicable.
- **2. Noiseless baseline:** Aer without noise gives:
  - `prep0`: 0 ones;
  - `prep1`: all ones;
  - `ampN`: about half;
  - `x100`: 0 errors.

  See `tests/test_amplitude_drift.py::test_noiseless_baseline`.
- **3. Fail closed:**
  - bridge: ruff clean, 57 tests pass;
  - osiris-governance: 271 pass, 3 skipped;
  - `git diff --check` clean.
- **4. Ledger:** rows are written ahead and hash-chained, but unsigned. They need a signature from an enrolled node key.
- **5. Heron variance:** K = 12 ≥ 10 jobs.
- **`GEMINI.md`:**
  - **Run:** the user requested it in their own words ("run x100 amplitude drift test after this finishes").
  - **Duration:** about 1 h, longer than the node horizon. The user approved E2's equally long run; this one is flagged here.
  - **Thresholds:** α, the control steps and S were chosen by the agent and need the human's review.

## Code and environment (pinned)

| file | sha256 |
|---|---|
| `bridge/amplitude_drift.py` | `19aeaaf63a330c89ce3dea16b4517ffb7bb5b0639f0d8fde1c3d6a4b556dea98` |
| `bridge/drift_check.py` | `95b01edcf5ba54dce154236ba48cb386a1a073295b64f6568615b042219cbe60` |
| `osiris_governance/telemetry.py` | `f805f8bd9871f87c7a5463071a6b28b52036a8fd04e5044b7b8f5e51d378a420` |

**Repositories:**
- bridge `6b39f0c`, with `amplitude_drift.py` and its tests uncommitted
- osiris-governance `e8a2eb3`

**Versions:** Python 3.13.5, numpy 2.5.3, scipy 1.18.1, qiskit 2.5.2, qiskit-aer 0.17.2, qiskit-ibm-runtime 0.50.0.

Command (from `bridge/`):

```bash
python -m bridge.amplitude_drift hardware --key-file <IBM key JSON> --out results/drift_checks/amplitude \
    --prereg results/drift_checks/AMPLITUDE_PREREGISTRATION.md --after results/drift_checks/hardware/state.json
```
