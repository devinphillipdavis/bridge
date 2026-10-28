"""
bridge.pulse_amplitude — is the drifting X rotation error a drive-amplitude change? (E4)
========================================================================================

Pre-registered in ``results/drift_checks/E4_PREREGISTRATION.md``; follows E3 on the same backend.
Each job runs, on every qubit in parallel:

    prep0, prep1   measure; x, measure                   P(1) = c - a/2; c + (a/2) cos(eps)
    ampN           sx, then N x gates, then measure      P(1) = c + (a/2) exp(-gamma N) sin(delta + N eps)
    sxM            M = 4m + 1 sx gates, then measure     P(1) = c + (a/2) exp(-gamma M) sin(M delta)
    ram0, ram90    sx, delay, [rz(pi/2)], sx, measure    P(1) = c + r cos(phi); c - r sin(phi)
    x100           100 x gates, then measure             (E2's circuit; continuity only, not fitted)

where x rotates by pi + eps and sx by pi/2 + delta. A shared drive gain g moves both
(d delta / d eps = 1/2), a re-tune of the x amplitude alone moves eps only (0), and a detuning
moves delta twice as much as eps (2) while shifting the Ramsey phase phi.

Per qubit: a pooled fit of (c, a, gamma, delta, eps); then each job's (c, a, delta, eps) by
weighted Gauss-Newton on frequencies that first lose the run's static residual, which E3 showed
every job shares. Replicas are binomial draws at the run's mean frequencies (no drift, same static
residual) through the same estimator. H1 (variation) gates H2 (ratio), H2F (frequency), H3 (jumps)
and H4 (shared changes).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from osiris_governance import telemetry
from osiris_governance.telemetry import TelemetryRun, quantity
from scipy.optimize import minimize
from scipy.stats import rankdata

from dnalang.ledger import Ledger
from dnalang.qc_ir import Circuit, Op

from . import drift_check as dc

X_REPS = (0, 16, 32, 48)
SX_REPS = (17, 33, 49)                  # 4m + 1 sx gates: a quarter turn plus M delta
FITTED = ("prep0", "prep1") + tuple(f"amp{n}" for n in X_REPS) + tuple(f"sx{m}" for m in SX_REPS)
RAMSEY = ("ram0", "ram90")
CIRCUIT_NAMES = FITTED + RAMSEY + ("x100",)
PER_JOB = ("c", "a", "delta", "eps")    # re-estimated in every job; gamma stays pooled
RAMSEY_DELAY_S = 200e-9                 # +-2.5 MHz unambiguous
X_LENGTH_S = 32e-9                      # simulation only; hardware records the backend's lengths
MICRO = 10 ** 6
BOUNDS = [(0.0, 0.5), (0.0, 0.5), (0.0, 0.01), (-0.3, 0.3), (-0.1, 0.1)]   # e0, e1, gamma, delta, eps
_NX = np.asarray(X_REPS, float)
_MS = np.asarray(SX_REPS, float)
_AMP = slice(2, 2 + len(X_REPS))
_SX = slice(_AMP.stop, len(FITTED))

# E3's pooled (c, a, gamma, delta, eps) and mean residual (observed minus fitted P(1)) for prep0,
# prep1, amp0, amp16, amp32 and amp48 on qubits 4, 22, 54 and 82 (results/drift_checks/amplitude,
# computed 2026-10-27). The simulations reuse the amp16-48 residuals for sx17-49 and scale every
# residual by 4 p (1 - p), so that it vanishes where P(1) is 0 or 1.
E3_THETAS = ((0.498, 0.993, 0.0, -0.040, 0.0080), (0.498, 0.993, 0.0, -0.058, 0.0052),
             (0.496, 0.991, 0.0, -0.048, 0.0123), (0.499, 0.972, 0.0, -0.087, 0.0028))
E3_MISFIT = ((-0.0003, 0.0, 0.0225, -0.0145, -0.0273, -0.0203),
             (-0.0005, 0.0005, 0.0290, -0.0089, -0.0230, -0.0357),
             (-0.0009, 0.0009, 0.0154, -0.0122, -0.0204, -0.0070),
             (0.0, -0.0001, 0.0418, -0.0163, -0.0326, -0.0298))
E3_WALK_Q = 4.3e-9                      # rad^2 / s: median random-walk rate of E3's eps, jobs 8-12


@dataclass(frozen=True)
class E4Design:
    jobs: int = 30
    spacing_s: int = 150
    shots: int = 1000
    replicas: int = 2000
    alpha: float = dc.ALPHA
    analysis_seed: int = 2026102801
    control_seed: int = 2026102802
    sim_seed: int = 2026102803
    sim_trials: int = 200
    sim_count_trials: int = 20           # S4_counts: the gain scenario through the full count-level analysis
    extra_qubits: int = 4
    ramsey_delay_s: float = RAMSEY_DELAY_S
    iterations: int = 5                  # Gauss-Newton steps per job
    gain_band: Tuple[float, float] = (0.25, 1.0)
    x_only_max: float = 0.25
    detuning_min: float = 1.25
    ramsey_limit_hz: float = 200e3
    min_jobs: int = 20
    min_budget_s: float = 200.0
    short_gap_s: float = 300.0           # H3 needs min_short_gaps intervals at most this long
    min_short_gaps: int = 15
    step_rad: float = 0.005              # S3 and S5 step; S4 drift over the run
    detuning_hz: float = 0.9e6           # S6 step


def pulse_circuits(n: int, delay_s: float = RAMSEY_DELAY_S) -> List[Circuit]:
    measure = [Op("measure", (q,), cbit=q) for q in range(n)]
    half = [Op("sx", (q,)) for q in range(n)]
    flips = [Op("x", (q,)) for q in range(n)]
    wait = [Op("delay", (q,), duration_s=delay_s) for q in range(n)]
    quarter = [Op("rz", (q,), (float(np.pi / 2),)) for q in range(n)]
    circuits = [Circuit(n, n, list(measure), name="prep0"), Circuit(n, n, flips + measure, name="prep1")]
    circuits += [Circuit(n, n, half + flips * reps + measure, name=f"amp{reps}") for reps in X_REPS]
    circuits += [Circuit(n, n, half * reps + measure, name=f"sx{reps}") for reps in SX_REPS]
    circuits += [Circuit(n, n, half + wait + half + measure, name="ram0"),
                 Circuit(n, n, half + wait + quarter + half + measure, name="ram90")]
    return circuits + [Circuit(n, n, flips * 100 + measure, name="x100")]


def ones(counts: Dict[str, int], n: int) -> List[int]:
    """Per-qubit count of outcome 1; qubit q is classical bit q, read from the right."""
    return dc.error_counts(counts, n, "prep0")


# --------------------------------------------------------------------------------------------
# Model and estimators
# --------------------------------------------------------------------------------------------

def _join(parts: Sequence[np.ndarray]) -> np.ndarray:
    lead = np.broadcast_shapes(*(p.shape[:-1] for p in parts))
    return np.concatenate([np.broadcast_to(p, lead + p.shape[-1:]) for p in parts], axis=-1)


def _columns(*values):
    return tuple(np.asarray(v, float)[..., None] for v in values)


def p_one(c, a, gamma, delta, eps) -> np.ndarray:
    """P(1) of each FITTED circuit; arguments broadcast over leading axes."""
    c, a, gamma, delta, eps = _columns(c, a, gamma, delta, eps)
    out = _join([c - a / 2, c + a / 2 * np.cos(eps),
                 c + a / 2 * np.exp(-gamma * _NX) * np.sin(delta + _NX * eps),
                 c + a / 2 * np.exp(-gamma * _MS) * np.sin(_MS * delta)])
    return np.clip(out, 1e-9, 1 - 1e-9)


def jacobian(c, a, gamma, delta, eps) -> np.ndarray:
    """d p_one / d (c, a, delta, eps), shaped (..., len(FITTED), 4)."""
    c, a, gamma, delta, eps = _columns(c, a, gamma, delta, eps)
    ex, es = np.exp(-gamma * _NX), np.exp(-gamma * _MS)
    sx, cx = np.sin(delta + _NX * eps), np.cos(delta + _NX * eps)
    ss, cs = np.sin(_MS * delta), np.cos(_MS * delta)
    zero, one = 0 * c, 0 * c + 1
    d_c = _join([one, one, one + 0 * ex, one + 0 * es])
    d_a = _join([zero - 0.5, np.cos(eps) / 2, ex * sx / 2, es * ss / 2])
    d_delta = _join([zero, zero, a / 2 * ex * cx, a / 2 * es * _MS * cs])
    d_eps = _join([zero, -a / 2 * np.sin(eps), a / 2 * ex * _NX * cx, zero * es])
    return np.stack([d_c, d_a, d_delta, d_eps], axis=-1)


def nll(theta: Sequence[float], y: np.ndarray, shots: int) -> float:
    """Binomial negative log-likelihood of counts y (jobs, len(FITTED)), up to a constant."""
    p = p_one(*theta)
    return float(-(y * np.log(p) + (shots - y) * np.log1p(-p)).sum())


def _from_errors(e0: float, e1: float) -> Tuple[float, float]:
    """(c, a) from the readout errors P(1|0) = e0 and P(0|1) = e1, so that 0 <= P(1) <= 1."""
    return (1 + e0 - e1) / 2, 1 - e0 - e1


def _grid_start(y: np.ndarray, shots: int) -> Tuple[float, ...]:
    """c and a from prep0 and prep1; delta from amp0 and the sx ladders; then eps from the x ladders."""
    f = y.sum(0) / (shots * len(y))
    c, a = (f[0] + f[1]) / 2, max(f[1] - f[0], 1e-3)
    deltas = np.arange(-0.3, 0.3 + 1e-12, 0.0005)
    turns = np.concatenate([[1.0], _MS])
    target = np.concatenate([f[_AMP][:1], f[_SX]])
    delta = deltas[np.argmin(((c + a / 2 * np.sin(deltas[:, None] * turns) - target) ** 2).sum(1))]
    epss = np.arange(-0.1, 0.1 + 1e-12, 0.0002)
    eps = epss[np.argmin(((c + a / 2 * np.sin(delta + epss[:, None] * _NX) - f[_AMP]) ** 2).sum(1))]
    return float(c), float(a), 0.0, float(delta), float(eps)


def fit_pooled(y: np.ndarray, shots: int) -> np.ndarray:
    """Maximum-likelihood (c, a, gamma, delta, eps) over every job, from the grid start; fitted as
    readout errors (e0, e1) in [0, 0.5] so that no circuit's P(1) leaves [0, 1]."""
    c, a, gamma, delta, eps = _grid_start(y, shots)
    start = (min(max(c - a / 2, 0.0), 0.5), min(max(1 - c - a / 2, 0.0), 0.5), gamma, delta, eps)

    def cost(x):
        return nll((*_from_errors(x[0], x[1]), *x[2:]), y, shots)

    x = minimize(cost, start, method="L-BFGS-B", bounds=BOUNDS).x
    return np.array([*_from_errors(x[0], x[1]), *x[2:]])


def per_job(freq: np.ndarray, theta: Sequence[float], shots: int, iterations: int) -> np.ndarray:
    """Each job's (c, a, delta, eps) by weighted Gauss-Newton from the pooled fit; gamma held.

    freq is (..., jobs, len(FITTED)). Every job is fitted with the same model, so a static misfit
    biases every job alike and cancels in the changes. Weights come from the run's mean
    frequencies, smoothed by half a count, so none is unbounded.
    """
    c, a, gamma, delta, eps = theta
    mean = freq.mean(-2, keepdims=True)
    n = freq.shape[-2] * shots
    smooth = (mean * n + 0.5) / (n + 1)
    weight = shots / (smooth * (1 - smooth))
    phi = np.broadcast_to(np.array([c, a, delta, eps], float), freq.shape[:-1] + (4,)).copy()
    for _ in range(iterations):
        args = (phi[..., 0], phi[..., 1], gamma, phi[..., 2], phi[..., 3])
        jac = jacobian(*args)
        weighted = jac * weight[..., None]
        normal = np.swapaxes(jac, -1, -2) @ weighted
        gradient = np.swapaxes(weighted, -1, -2) @ (freq - p_one(*args))[..., None]
        phi += np.linalg.solve(normal, gradient)[..., 0]
    return phi


def fisher_cov(theta: Sequence[float], shots: int) -> np.ndarray:
    """Covariance of one job's (c, a, delta, eps) from the binomial Fisher information."""
    jac = jacobian(*theta)
    p = p_one(*theta)
    return np.linalg.inv(jac.T @ (jac * (shots / (p * (1 - p)))[:, None]))


def ramsey_phase(f0: np.ndarray, f90: np.ndarray, c: np.ndarray) -> np.ndarray:
    """phi from P(ram0) = c + r cos(phi) and P(ram90) = c - r sin(phi). Only changes in phi are used."""
    return np.arctan2(-(f90 - c), f0 - c)


def frequency_shift(phase: np.ndarray, delay_s: float) -> np.ndarray:
    """Each job's frequency change from the run's circular-mean phase, in Hz."""
    reference = np.angle(np.exp(1j * phase).mean(-1, keepdims=True))
    return np.angle(np.exp(1j * (phase - reference))) / (2 * np.pi * delay_s)


def _micro(values: np.ndarray) -> List[int]:
    return [int(round(v * MICRO)) for v in values]


# --------------------------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------------------------

def deming(x: np.ndarray, y: np.ndarray, ratio: float) -> np.ndarray:
    """Deming slope of y on x over the last axis; ratio = error variance of y / error variance of x."""
    xc, yc = x - x.mean(-1, keepdims=True), y - y.mean(-1, keepdims=True)
    n = x.shape[-1] - 1
    sxx, syy, sxy = (xc * xc).sum(-1) / n, (yc * yc).sum(-1) / n, (xc * yc).sum(-1) / n
    d = syy - ratio * sxx
    with np.errstate(divide="ignore", invalid="ignore"):
        return (d + np.sqrt(d * d + 4 * ratio * sxy * sxy)) / (2 * sxy)


def ratio_region(low: float, high: float, design: E4Design) -> str:
    if design.gain_band[0] <= low and high <= design.gain_band[1]:
        return "gain"
    if low <= 0.0 <= high and high < design.x_only_max:
        return "x_only"
    if low > design.detuning_min:
        return "detuning_like"
    return "inconclusive"


def ratio_test(eps: Dict[int, np.ndarray], delta: Dict[int, np.ndarray], var_eps: Dict[int, float],
               var_delta: Dict[int, float], rng: np.random.Generator, design: E4Design) -> Dict[str, Any]:
    """H2: Deming beta = d delta / d eps per qubit, pooled by inverse bootstrap variance, with a
    95% pairs-bootstrap interval (the same resampled jobs for every qubit)."""
    qubits = sorted(eps)
    k = len(eps[qubits[0]])
    draws = rng.integers(0, k, size=(design.replicas, k))
    per, boot = {}, {}
    for q in qubits:
        ratio = var_delta[q] / var_eps[q]
        per[q] = float(deming(eps[q], delta[q], ratio))
        boot[q] = deming(eps[q][draws], delta[q][draws], ratio)
    usable = [q for q in qubits if np.isfinite(per[q]) and np.isfinite(boot[q]).sum() > 1]
    if not usable:
        return {"per_qubit": per, "pooled": None, "interval": None, "region": "inconclusive"}
    weights = {q: 1.0 / float(np.var(boot[q][np.isfinite(boot[q])], ddof=1)) for q in usable}
    total = sum(weights.values())
    pooled = sum(weights[q] * per[q] for q in usable) / total
    stacked = np.array([weights[q] * boot[q] for q in usable]).sum(0) / total
    stacked = stacked[np.isfinite(stacked)]
    low, high = (float(v) for v in np.percentile(stacked, [2.5, 97.5]))
    return {"per_qubit": per, "weights": weights, "pooled": float(pooled), "interval": [low, high],
            "draws_used": int(stacked.size), "region": ratio_region(low, high, design)}


def frequency_test(shift_hz: Dict[int, np.ndarray], eps: Dict[int, np.ndarray], rng: np.random.Generator,
                   design: E4Design) -> Dict[str, Any]:
    """H2F: per qubit, max |frequency change| below the limit and no Holm-flagged Spearman
    correlation (two-sided, permutation null) between frequency change and eps."""
    p, max_hz, rho = {}, {}, {}
    for q in sorted(shift_hz):
        rf, re = rankdata(shift_hz[q]), rankdata(eps[q])
        rf, re = rf - rf.mean(), re - re.mean()
        norm = np.sqrt((rf * rf).sum() * (re * re).sum())
        rho[q] = float((rf * re).sum() / norm) if norm else 0.0
        perms = rng.permuted(np.broadcast_to(rf, (design.replicas, rf.size)), axis=1)
        null = np.abs(perms @ re) / norm if norm else np.zeros(design.replicas)
        p[q] = float((1 + (null >= abs(rho[q]) - 1e-12).sum()) / (1 + design.replicas))
        max_hz[q] = float(np.abs(shift_hz[q]).max())
    flagged = dc.holm({str(q): v for q, v in p.items()}, design.alpha)
    passed = {q: max_hz[q] < design.ramsey_limit_hz and not flagged[str(q)] for q in p}
    return {"max_abs_hz": max_hz, "spearman": rho, "p": p, "passed": passed}


def _walk_statistic(x: np.ndarray, t: np.ndarray, v: np.ndarray):
    d, dt, vv = np.diff(x, axis=-1), np.diff(t), v[:-1] + v[1:]
    mu = d.sum(-1, keepdims=True) / dt.sum()
    q = np.maximum(0.0, ((d - mu * dt) ** 2 - vv).sum(-1, keepdims=True) / dt.sum())
    z = (d - mu * dt) / np.sqrt(q * dt + vv)
    return np.abs(z).max(-1), z, mu, q


def jump_test(x: np.ndarray, t: np.ndarray, v: np.ndarray, rng: np.random.Generator, draws: int) -> Dict[str, Any]:
    """H3 for one series: max |standardised increment| against a Gaussian random walk with drift.

    The walk (drift mu, rate q per second) is fitted by moments with the measurement variances v
    known; the null is a parametric bootstrap from that fit, refitted the same way. (A median-based
    rate was tried and rejected: under a pure walk it false-alarmed in 0.2 of trials.)
    """
    stat, z, mu, q = _walk_statistic(x, t, v)
    dt = np.diff(t)
    steps = rng.normal(mu * dt, np.sqrt(q * dt), size=(draws, dt.size))
    walks = np.concatenate([np.zeros((draws, 1)), np.cumsum(steps, axis=1)], axis=1)
    null = _walk_statistic(walks + rng.normal(0.0, np.sqrt(v), size=(draws, t.size)), t, v)[0]
    return {"statistic": float(stat), "p": float((1 + (null >= stat).sum()) / (1 + draws)),
            "interval": int(np.abs(z).argmax()), "drift_per_s": float(mu[0]), "walk_q_per_s": float(q[0])}


def shared_test(x: np.ndarray, rng: np.random.Generator, draws: int) -> Dict[str, Any]:
    """H4: mean pairwise correlation of the qubits' increments (rows of x) against independent
    circular shifts of each qubit's increments."""
    if len(x) < 2:
        return {"mean_correlation": None, "p": 1.0}
    d = np.diff(x, axis=-1)
    sd = d.std(-1, ddof=1, keepdims=True)
    z = np.divide(d - d.mean(-1, keepdims=True), sd, out=np.zeros_like(d), where=sd > 0)
    qn, n = z.shape

    def mean_corr(rows: np.ndarray) -> np.ndarray:
        return (((rows.sum(-2)) ** 2).sum(-1) - (rows ** 2).sum((-1, -2))) / (qn * (qn - 1) * (n - 1))

    stat = float(mean_corr(z))
    shifts = rng.integers(0, n, size=(draws, qn))
    index = (np.arange(n)[None, None, :] + shifts[:, :, None]) % n
    null = mean_corr(z[np.arange(qn)[None, :, None], index])
    return {"mean_correlation": stat, "p": float((1 + (null >= stat - 1e-12).sum()) / (1 + draws))}


def coincidence_test(intervals: Sequence[int], n_intervals: int, rng: np.random.Generator,
                     draws: int) -> Dict[str, Any]:
    """How many flagged jumps share one interval, against jumps placed independently at random."""
    stat = int(np.bincount(intervals).max())
    null = rng.integers(0, n_intervals, size=(draws, len(intervals)))
    null_max = np.array([np.bincount(row).max() for row in null])
    return {"most_in_one_interval": stat, "p": float((1 + (null_max >= stat).sum()) / (1 + draws))}


def reading(result: Dict[str, Any]) -> str:
    """The pre-registered reading of H1-H4."""
    if not result["H1_rotation_error_varies"]:
        return "not evaluable: no variation in rotation error in E4's window (not evidence of absence)"
    region = result["H2"]["region"]
    passed = result["H2F"]["passed"]
    falsified = sorted(q for q, ok in passed.items() if not ok)
    if region == "detuning_like" or (passed and len(falsified) == len(passed)):
        return "frequency, not amplitude: falsifies the amplitude hypothesis"
    jumps, shared = result["H3_jumps"], result["H4_shared"]
    if region == "gain":
        text = ("shared analog gain drift" if shared and not jumps else
                "step changes in a shared gain, or a joint x and sx re-tune" if shared else
                "gain change, not shown to be shared")
    elif region == "x_only":
        text = ("x-amplitude re-tunes" if jumps and shared else
                "x-amplitude jumps, not shown to be shared" if jumps else
                "an x-only change that is not gain: mechanism inconclusive")
    else:
        text = "amplitude not excluded; mechanism inconclusive"
    if falsified:
        text += f"; H2F falsifies amplitude on qubits {falsified}"
    return text


# --------------------------------------------------------------------------------------------
# Analysis
# --------------------------------------------------------------------------------------------

def estimate(y: np.ndarray, shots: int, design: E4Design, delay_s: float,
             rng: np.random.Generator) -> Dict[str, Any]:
    """One qubit: pooled fit, per-job estimates, Ramsey shifts, and null replicas."""
    fitted = y[:, :len(FITTED)]
    freq = fitted / shots
    theta = fit_pooled(fitted, shots)
    phi = per_job(freq, theta, shots, design.iterations)
    phase = ramsey_phase(y[:, len(FITTED)] / shots, y[:, len(FITTED) + 1] / shots, phi[:, 0])
    replicas = rng.binomial(shots, freq.mean(0), size=(design.replicas,) + freq.shape) / shots
    rep = per_job(replicas, theta, shots, design.iterations)
    model = p_one(*phi.T[:2], theta[2], *phi.T[2:])
    z = (freq - model) / np.sqrt(model * (1 - model) / shots)
    return {"theta": theta, "phi": phi, "shift_hz": frequency_shift(phase, delay_s),
            "replicas": rep, "var": rep.var(0, ddof=1), "residual_z": z}


def analyze(times: Sequence[float], y_by_qubit: Dict[int, np.ndarray], shots: int, design: E4Design,
            seed: int, delay_s: float, gate: bool = True,
            x100_by_qubit: Optional[Dict[int, Sequence[int]]] = None) -> Dict[str, Any]:
    """H1, then (if H1 holds, or always with gate=False) H2, H2F, H3 and H4."""
    rng = np.random.default_rng(seed)
    t = np.asarray(times, float)
    est = {q: estimate(np.asarray(y), shots, design, delay_s, rng) for q, y in y_by_qubit.items()}
    p1, qubits = {}, {}
    for q, e in est.items():
        record: Dict[str, Any] = {
            "pooled": dict(zip(("c", "a", "gamma", "delta", "eps"), map(float, e["theta"]))),
            "per_job": {name: [float(v) for v in e["phi"][:, i]] for i, name in enumerate(PER_JOB)},
            "frequency_shift_hz": [float(v) for v in e["shift_hz"]],
            "measurement_sd": {name: float(np.sqrt(e["var"][:, i].mean())) for i, name in enumerate(PER_JOB)},
            "residual_z": {"mean": dict(zip(FITTED, map(float, e["residual_z"].mean(0)))),
                           "sd": dict(zip(FITTED, map(float, e["residual_z"].std(0, ddof=1))))},
            "p": {},
        }
        for name, column in (("rotation_error", 3), ("contrast", 1)):
            observed = telemetry.dispersion(_micro(e["phi"][:, column]))
            null = [telemetry.dispersion(_micro(r)) for r in e["replicas"][:, :, column]]
            p = float(telemetry.predictive_p_value(observed, null))
            record["p"][name] = {"dispersion": p}
            p1[f"q{q}/{name}/dispersion"] = p
        if x100_by_qubit is not None:
            record["x100_errors"] = [int(v) for v in x100_by_qubit[q]]
        qubits[str(q)] = record
    rejected = dc.holm(p1, design.alpha)
    flagged = sorted(k for k, r in rejected.items() if r)
    varying = sorted(q for q in est if rejected[f"q{q}/rotation_error/dispersion"])
    result: Dict[str, Any] = {
        "qubits": qubits, "H1_tests": len(p1), "alpha": design.alpha, "replicas": design.replicas, "seed": seed,
        "flagged": flagged, "H1_rotation_error_varies": bool(varying),
        "HC_contrast_varies": any(k.endswith("/contrast/dispersion") for k in flagged),
    }
    if gate and not varying:
        result["reading"] = reading(result)
        return result
    tested = varying if varying else sorted(est)
    eps = {q: est[q]["phi"][:, 3] for q in tested}
    result["H2"] = ratio_test(eps, {q: est[q]["phi"][:, 2] for q in tested},
                              {q: float(est[q]["var"][:, 3].mean()) for q in tested},
                              {q: float(est[q]["var"][:, 2].mean()) for q in tested}, rng, design)
    result["H2F"] = frequency_test({q: est[q]["shift_hz"] for q in tested}, eps, rng, design)
    short = int((np.diff(t) <= design.short_gap_s).sum())
    if short < design.min_short_gaps:
        result["H3"] = {"underpowered": True, "short_intervals": short}
        result["H3_jumps"] = False
    else:
        jumps = {q: jump_test(eps[q], t, est[q]["var"][:, 3], rng, design.replicas) for q in tested}
        rejected3 = dc.holm({str(q): j["p"] for q, j in jumps.items()}, design.alpha)
        for q, j in jumps.items():
            j["flagged"] = rejected3[str(q)]
        result["H3"] = {"underpowered": False, "short_intervals": short, "per_qubit": jumps}
        result["H3_jumps"] = any(rejected3.values())
    allq = sorted(est)
    result["H4"] = shared_test(np.array([est[q]["phi"][:, 3] for q in allq]), rng, design.replicas)
    result["H4_shared"] = result["H4"]["p"] <= design.alpha
    if result["H3_jumps"]:
        hit = [j for j in result["H3"]["per_qubit"].values() if j["flagged"]]
        if len(hit) >= 2:
            result["H4"]["coincidence"] = coincidence_test([j["interval"] for j in hit], len(t) - 1, rng,
                                                           design.replicas)
        for q, j in result["H3"]["per_qubit"].items():
            if j["flagged"]:
                before, after = est[q]["phi"][j["interval"]:j["interval"] + 2, 3]
                j["abs_eps_before_after"] = [float(before), float(after)]
                j["abs_eps_decreased"] = bool(abs(after) < abs(before))
    result["reading"] = reading(result)
    return result


# --------------------------------------------------------------------------------------------
# Software control and simulation stage (E4-sim)
# --------------------------------------------------------------------------------------------

def sim_theta(i: int) -> Tuple[float, ...]:
    return E3_THETAS[i % len(E3_THETAS)]


def sim_misfit(i: int) -> np.ndarray:
    row = E3_MISFIT[i % len(E3_MISFIT)]
    return np.array(row + row[3:6])


def synthetic_counts(i: int, eps: np.ndarray, delta: np.ndarray, shift_hz: np.ndarray, design: E4Design,
                     rng: np.random.Generator) -> np.ndarray:
    """Counts (jobs, len(CIRCUIT_NAMES)) for sim qubit i: per-job eps and delta, E3's static residual
    for that qubit, Ramsey phase 2 pi f tau, x100 from the model."""
    c, a, gamma, _, _ = sim_theta(i)
    eps, delta, shift_hz = (np.broadcast_to(np.asarray(v, float), (design.jobs,)) for v in (eps, delta, shift_hz))
    ideal = p_one(c, a, gamma, delta, eps)
    p = np.clip(ideal + sim_misfit(i) * 4 * ideal * (1 - ideal), 1e-6, 1 - 1e-6)
    phase = 2 * np.pi * shift_hz * design.ramsey_delay_s
    ram = np.stack([c + a / 2 * np.cos(phase), c - a / 2 * np.sin(phase)], axis=-1)
    x100 = (c - a / 2 + a * np.sin(50 * eps) ** 2)[:, None]
    return rng.binomial(design.shots, np.concatenate([p, ram, x100], axis=-1))


def software_control(design: E4Design, qubits: int = 8) -> Dict[str, Any]:
    """S1 at count level: constant parameters with E3's static residuals; every test, ungated."""
    rng = np.random.default_rng(design.control_seed)
    k = design.jobs
    y = {}
    for i in range(qubits):
        _, _, _, delta, eps = sim_theta(i)
        y[i] = synthetic_counts(i, np.full(k, eps), np.full(k, delta), np.zeros(k), design, rng)
    times = [design.spacing_s * j for j in range(k)]
    result = analyze(times, y, design.shots, design, design.control_seed, design.ramsey_delay_s, gate=False)
    result["dispatch_allowed"] = not (result["flagged"] or result["H3_jumps"] or result["H4_shared"])
    return result


def _series_noise(i: int, design: E4Design) -> Tuple[np.ndarray, float]:
    c, a, gamma, delta, eps = sim_theta(i)
    cov = fisher_cov((c, a, gamma, delta, eps), design.shots)[2:, 2:]     # (delta, eps)
    return cov, 1.0 / (a * np.sqrt(design.shots)) / (2 * np.pi * design.ramsey_delay_s)


def sim_truth(scenario: str, design: E4Design, rng: np.random.Generator, qubits: int = 8):
    """Times and each sim qubit's true per-job (eps, delta, frequency shift) for one trial."""
    k = design.jobs
    t = design.spacing_s * np.arange(k, dtype=float)
    step_at = rng.integers(1, k)
    late = (np.arange(k) >= step_at).astype(float)
    out = {}
    for i in range(qubits):
        _, _, _, delta0, eps0 = sim_theta(i)
        walk = np.concatenate([[0.0], np.cumsum(rng.normal(0, np.sqrt(E3_WALK_Q * np.diff(t))))])
        eps, delta, shift = eps0 + 0 * t, delta0 + 0 * t, 0 * t
        if scenario in ("S2", "S3", "S3_one", "S5"):
            eps = eps + walk
        if scenario == "S3" or (scenario == "S3_one" and i == 0):
            eps = eps + design.step_rad * late
        if scenario == "S4":
            gain = (design.step_rad * t / t[-1] + walk) / np.pi
            eps, delta = eps0 + np.pi * gain, delta0 + np.pi / 2 * gain
        if scenario == "S5":
            eps = eps + design.step_rad * late
        if scenario == "S6":
            alpha = 2 * np.pi * design.detuning_hz * X_LENGTH_S / np.pi
            eps, delta = eps + np.pi * alpha ** 2 / 2 * late, delta + np.pi * alpha ** 2 * late
            shift = design.detuning_hz * late
        out[i] = (eps, delta, shift)
    return t, out


def sim_series(scenario: str, design: E4Design, rng: np.random.Generator, qubits: int = 8):
    """Per-qubit (eps_hat, delta_hat, shift_hat, var_eps, var_delta) for one trial of a scenario,
    with the estimates' Fisher noise (delta and eps errors correlated)."""
    t, truth = sim_truth(scenario, design, rng, qubits)
    out = {}
    for i, (eps, delta, shift) in truth.items():
        cov, sd_hz = _series_noise(i, design)
        noise = rng.multivariate_normal([0.0, 0.0], cov, size=t.size)
        out[i] = (eps + noise[:, 1], delta + noise[:, 0], shift + rng.normal(0, sd_hz, t.size), cov[1, 1], cov[0, 0])
    return t, out


def sim_trial(scenario: str, design: E4Design, rng: np.random.Generator) -> Dict[str, Any]:
    t, series = sim_series(scenario, design, rng)
    qs = sorted(series)
    eps = {q: series[q][0] for q in qs}
    row: Dict[str, Any] = {}
    if scenario in ("S1", "S2", "S3", "S3_one"):
        jumps = {q: jump_test(eps[q], t, np.full(t.size, series[q][3]), rng, design.replicas)["p"] for q in qs}
        row["H3"] = any(dc.holm({str(q): p for q, p in jumps.items()}, design.alpha).values())
        row["H4"] = shared_test(np.array([eps[q] for q in qs]), rng, design.replicas)["p"] <= design.alpha
    if scenario in ("S4", "S5"):
        h2 = ratio_test(eps, {q: series[q][1] for q in qs}, {q: series[q][3] for q in qs},
                        {q: series[q][4] for q in qs}, rng, design)
        target = 0.5 if scenario == "S4" else 0.0
        row["covers"] = h2["interval"] is not None and h2["interval"][0] <= target <= h2["interval"][1]
        row["region"] = h2["region"]
    if scenario == "S6":
        h2f = frequency_test({q: series[q][2] for q in qs}, eps, rng, design)
        row["detected"] = not any(h2f["passed"].values())
    return row


SCENARIOS = ("S1", "S2", "S3", "S3_one", "S4", "S5", "S6")   # S3_one: the step on one qubit; reported only


def simulate(out: Path, design: E4Design = E4Design()) -> Dict[str, Any]:
    """E4-sim: S1 at count level (the dispatch gate), S1-S6 at series level over sim_trials, and the
    gain scenario through the full count-level analysis over sim_count_trials (reported only)."""
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    control = software_control(design)
    rng = np.random.default_rng(design.sim_seed)
    rates: Dict[str, Dict[str, float]] = {}
    regions: Dict[str, Dict[str, int]] = {}
    for scenario in SCENARIOS:
        rows = [sim_trial(scenario, design, rng) for _ in range(design.sim_trials)]
        keys = [k for k in rows[0] if k != "region"]
        rates[scenario] = {k: float(np.mean([r[k] for r in rows])) for k in keys}
        if "region" in rows[0]:
            regions[scenario] = {name: sum(r["region"] == name for r in rows)
                                 for name in ("gain", "x_only", "detuning_like", "inconclusive")}
    counted = []
    for _ in range(design.sim_count_trials):
        t, truth = sim_truth("S4", design, rng)
        y = {i: synthetic_counts(i, *truth[i], design, rng) for i in truth}
        result = analyze(t, y, design.shots, design, int(rng.integers(2 ** 31)), design.ramsey_delay_s)
        h2 = result.get("H2") or {}
        interval = h2.get("interval")
        counted.append({"H1": result["H1_rotation_error_varies"],
                        "covers": bool(interval and interval[0] <= 0.5 <= interval[1]),
                        "region": h2.get("region", "not evaluated"), "beta": h2.get("pooled"), "interval": interval})
    rates["S4_counts"] = {k: float(np.mean([r[k] for r in counted])) for k in ("H1", "covers")}
    betas = [r["beta"] for r in counted if r["beta"] is not None]
    rates["S4_counts"]["median_beta"] = float(np.median(betas)) if betas else None
    rates["S4_counts"]["median_interval_width"] = (float(np.median([r["interval"][1] - r["interval"][0]
                                                                    for r in counted if r["interval"]]))
                                                   if betas else None)
    regions["S4_counts"] = {name: sum(r["region"] == name for r in counted)
                            for name in ("gain", "x_only", "detuning_like", "inconclusive", "not evaluated")}
    verdict = {
        "S1_dispatch_allowed": control["dispatch_allowed"],
        "S2_H3_alarm_le_0.06": rates["S2"]["H3"] <= 0.06,
        "S3_H3_power_ge_0.80": rates["S3"]["H3"] >= 0.80,
        "S4_covers_half_ge_0.90": rates["S4"]["covers"] >= 0.90,
        "S5_covers_zero_ge_0.90": rates["S5"]["covers"] >= 0.90,
    }
    summary = {"design": asdict(design), "control": control, "rates": rates, "regions": regions,
               "verdict": verdict, "seconds": round(time.monotonic() - started, 1),
               "code_sha256": dc.sha256_file(Path(__file__))}
    (out / "oc.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    return summary


# --------------------------------------------------------------------------------------------
# Hardware
# --------------------------------------------------------------------------------------------

def choose_extras(backend, base: Sequence[int], n: int) -> List[int]:
    """The n qubits with the lowest reported readout error, skipping base, and taking none that
    neighbours base or an earlier pick in the coupling map."""
    errors = {qs[0]: props.error for qs, props in backend.target["measure"].items()
              if props is not None and props.error is not None}
    neighbours: Dict[int, set] = {}
    for a, b in backend.coupling_map.get_edges():
        neighbours.setdefault(a, set()).add(b)
        neighbours.setdefault(b, set()).add(a)
    blocked = set(base).union(*(neighbours.get(q, set()) for q in base))
    picks: List[int] = []
    for q in sorted(errors, key=lambda q: (errors[q], q)):
        if len(picks) == n:
            break
        if q not in blocked:
            picks.append(q)
            blocked |= {q} | neighbours.get(q, set())
    return picks


def aligned_delay(backend, delay_s: float) -> float:
    """delay_s rounded to the backend's timing granularity."""
    target = backend.target
    step = target.dt * max(1, int(getattr(target, "granularity", 1) or 1))
    return max(1, round(delay_s / step)) * step


def snapshot(props, qubits: Sequence[int]) -> Dict[str, Any]:
    """Reported calibration for the qubits: each qubit property and the x and sx gate properties,
    as [value, date]; plus last_update_date and the SHA-256 of the full properties."""
    full = json.dumps(props.to_dict(), sort_keys=True, default=str)
    per: Dict[str, Dict[str, Any]] = {}
    for q in qubits:
        record = {name: [value, str(date)] for name, (value, date) in props.qubit_property(q).items()}
        for gate in ("x", "sx"):
            try:
                found = props.gate_property(gate, [q])
            except Exception:
                continue
            record.update({f"{gate}_{name}": [value, str(date)] for name, (value, date) in found.items()})
        per[str(q)] = record
    return {"last_update_date": str(props.last_update_date), "properties_sha256": hashlib.sha256(full.encode()).hexdigest(),
            "qubits": per}


def _moment(stamp: str) -> datetime:
    return datetime.fromisoformat(stamp.replace("Z", "+00:00"))


def dated_inside(snapshots: Sequence[Dict[str, Any]], start: str, end: str,
                 qubits: Optional[Sequence[int]] = None) -> List[Dict[str, str]]:
    """Every (qubit, field, date) in the snapshots whose date falls in [start, end], once each."""
    lo, hi = _moment(start), _moment(end)
    found = {}
    for snap in snapshots:
        for q, fields in snap["qubits"].items():
            if qubits is not None and int(q) not in qubits:
                continue
            for name, (_, date) in fields.items():
                try:
                    moment = _moment(date)
                except ValueError:
                    continue
                if lo <= moment <= hi:
                    found[(q, name, date)] = {"qubit": q, "field": name, "date": date}
    return list(found.values())


def _midpoint(window: Sequence[str]) -> datetime:
    start, end = (_moment(w) for w in window)
    return start + (end - start) / 2


def retrospective(out: Path, key_file: Path, runs: Sequence[Path]) -> Dict[str, Any]:
    """E4a: reported calibration at each earlier job's midpoint, and what was dated inside each
    run's change interval (end of job 2 to start of job 3). Stops at the first failed query."""
    from dnalang.backends import ibm
    dc._use_key_file(key_file)
    svc = ibm.service()
    record: Dict[str, Any] = {"available": True, "runs": {}, "code_sha256": dc.sha256_file(Path(__file__))}
    for run in runs:
        state = json.loads((run / "state.json").read_text())
        results = dc._ledger_results(Ledger(run / "ledger.jsonl"))
        windows = sorted(results[j]["window"] for j in state["job_ids"] if results.get(j, {}).get("status") == "DONE")
        backend = svc.backend(state["backend"])
        layout = state["layout"][:4]
        snaps = []
        for window in windows:
            try:
                props = backend.properties(datetime=_midpoint(window))
            except Exception as err:            # the type only: messages may echo request details
                props, record["error"] = None, type(err).__name__
            if props is None:
                record["available"] = False
                break
            snaps.append({"at": window, **snapshot(props, layout)})
        interval = [windows[1][1], windows[2][0]]
        record["runs"][run.name] = {"change_interval": interval, "snapshots": snaps,
                                    "dated_inside": dated_inside(snaps, *interval)}
        if not record["available"]:
            break
    out.mkdir(parents=True, exist_ok=True)
    (out / "retrospective.json").write_text(json.dumps(record, indent=1, default=str) + "\n")
    return record


def write_telemetry(jobs: List[Dict[str, Any]], layout: Sequence[int], shots: int, run_id: str,
                    source: Dict[str, Any]) -> List[bytes]:
    run = TelemetryRun(run_id, source, "MEASURED" if source["kind"] == "HARDWARE" else "SIMULATED")
    lines = []
    for job in sorted(jobs, key=lambda j: j["window"][0]):
        metrics = {"shots": quantity(shots, 0, "1")}
        for name, counts in zip(CIRCUIT_NAMES, job["counts"]):
            label = "errors" if name == "x100" else "ones"
            for physical, count in zip(layout, ones(counts, len(layout))):
                metrics[f"q{physical}_{name}_{label}"] = quantity(count, 0, "1")
        start, end = job["window"]
        lines.append(run.frame(dc._now(), {"start": start, "end": end, "samples": shots}, metrics))
    return lines


def finish(out: Path, design: E4Design = E4Design()) -> Dict[str, Any]:
    """Telemetry and analysis from the ledger's result rows. Needs no network."""
    state = json.loads((out / "state.json").read_text())
    ledger = Ledger(out / "ledger.jsonl")
    if ledger.verify() is not None:
        raise RuntimeError(f"ledger does not verify: {ledger.verify()}")
    results = dc._ledger_results(ledger)
    done = sorted((results[j] for j in state["job_ids"] if results.get(j, {}).get("status") == "DONE"),
                  key=lambda r: r["window"][0])
    summary: Dict[str, Any] = {"jobs_submitted": len(state["job_ids"]), "jobs_done": len(done),
                               "qpu_seconds": sum(r.get("qpu_seconds") or 0 for r in done)}
    if len(done) < design.min_jobs:
        summary["verdict"] = f"INCOMPLETE: {len(done)} jobs finished, fewer than {design.min_jobs}"
    else:
        layout = state["layout"]
        source = {"kind": "HARDWARE", "node": state["backend"],
                  "generator": f"bridge.pulse_amplitude {dc.sha256_file(Path(__file__))[:12]}"}
        lines = write_telemetry([{"window": r["window"], "counts": r["counts"]} for r in done],
                                layout, design.shots, f"pulse-amplitude-hw-{state['backend']}", source)
        (out / "telemetry.jsonl").write_bytes(b"".join(line + b"\n" for line in lines))
        times = [float(v) for v in telemetry.metric_series(lines, f"q{layout[0]}_prep0_ones")[0]]
        y = {q: np.array([[int(v) for v in telemetry.metric_series(lines, f"q{q}_{name}_ones")[1]]
                          for name in FITTED + RAMSEY]).T for q in layout}
        x100 = {q: [int(v) for v in telemetry.metric_series(lines, f"q{q}_x100_errors")[1]] for q in layout}
        analysis = analyze(times, y, design.shots, design, design.analysis_seed, state["ramsey_delay_s"],
                           x100_by_qubit=x100)
        snaps = [row for row in ledger if row.get("kind") == "calibration-snapshot"]
        for q, jump in (analysis.get("H3", {}).get("per_qubit") or {}).items():
            if jump.get("flagged"):
                k = jump["interval"]
                jump["dated_inside"] = dated_inside(snaps, done[k]["window"][1], done[k + 1]["window"][0], [q])
        summary["analysis"] = analysis
    (out / "analysis.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    return summary


def hardware(out: Path, key_file: Path, prereg: Path, after: Path, design: E4Design = E4Design()) -> Dict[str, Any]:
    """Software control, budget and status checks, then design.jobs ledger-first jobs, each after a
    calibration snapshot, on E3's backend and qubits plus design.extra_qubits more; then results."""
    from dnalang.backends import ibm
    previous = json.loads(after.read_text())
    out.mkdir(parents=True, exist_ok=True)
    ledger = Ledger(out / "ledger.jsonl")
    hashes = {**dc.code_hashes(), "bridge/pulse_amplitude.py": dc.sha256_file(Path(__file__))}
    ledger.append("preregistration", {"sha256": dc.sha256_file(prereg), "code_sha256": hashes})
    control = software_control(design)
    (out / "control.json").write_text(json.dumps(control, indent=1, default=float) + "\n")
    if not control["dispatch_allowed"]:
        ledger.append("aborted", {"reason": "software control flagged constant data"})
        return {"verdict": "NOT RUN: the constant control flagged", "control": control["flagged"]}

    dc._use_key_file(key_file)
    svc = ibm.service()
    backend = svc.backend(previous["backend"])
    if not backend.status().operational:
        ledger.append("aborted", {"reason": f"{backend.name} not operational"})
        return {"verdict": f"NOT RUN: {backend.name} not operational"}
    remaining = svc.usage().get("usage_remaining_seconds")
    if remaining is None or remaining < design.min_budget_s:
        ledger.append("aborted", {"reason": "QPU budget below the minimum or unknown",
                                  "usage_remaining_seconds": remaining})
        return {"verdict": f"NOT RUN: {remaining} QPU-seconds remaining, minimum {design.min_budget_s}"}
    layout = list(previous["layout"]) + choose_extras(backend, previous["layout"], design.extra_qubits)
    delay = aligned_delay(backend, design.ramsey_delay_s)
    circuits = pulse_circuits(len(layout), delay)
    ledger.append("design", {**asdict(design), "backend": backend.name, "layout": layout, "ramsey_delay_s": delay,
                             "usage_remaining_seconds": remaining,
                             "circuits": {c.name: c.sha256() for c in circuits}})
    state = {"backend": backend.name, "layout": layout, "ramsey_delay_s": delay, "job_ids": []}
    for j in range(design.jobs):
        began = time.monotonic()
        props = svc.backend(backend.name).properties(refresh=True)
        ledger.append("calibration-snapshot", {"job": j, **snapshot(props, layout)})
        job_id = ibm.submit(circuits, backend.name, layout, design.shots, ledger,
                            tags={"experiment": "pulse_amplitude", "job": j})
        state["job_ids"].append(job_id)
        (out / "state.json").write_text(json.dumps(state, indent=1) + "\n")
        print(f"submitted {j + 1}/{design.jobs}: {job_id}", flush=True)
        if j + 1 < design.jobs:
            time.sleep(max(0.0, design.spacing_s - (time.monotonic() - began)))
    deadline = time.monotonic() + 8 * 3600
    while time.monotonic() < deadline:
        dc.collect(svc, ledger, state["job_ids"])
        if len(dc._ledger_results(ledger)) >= len(state["job_ids"]):
            break
        time.sleep(60)
    return finish(out, design)


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(prog="python -m bridge.pulse_amplitude")
    sub = parser.add_subparsers(dest="command", required=True)
    sim = sub.add_parser("simulate")
    sim.add_argument("--out", type=Path, required=True)
    retro = sub.add_parser("retrospective")
    retro.add_argument("--key-file", type=Path, required=True)
    retro.add_argument("--out", type=Path, required=True)
    retro.add_argument("--runs", type=Path, nargs="+", required=True)
    hw = sub.add_parser("hardware")
    hw.add_argument("--key-file", type=Path, required=True)
    hw.add_argument("--out", type=Path, required=True)
    hw.add_argument("--prereg", type=Path, required=True)
    hw.add_argument("--after", type=Path, required=True)
    fin = sub.add_parser("finish")
    fin.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "simulate":
        summary = simulate(args.out)
        print(json.dumps({"verdict": summary["verdict"], "rates": summary["rates"], "seconds": summary["seconds"]},
                         indent=1))
    elif args.command == "retrospective":
        record = retrospective(args.out, args.key_file, args.runs)
        print(json.dumps({"available": record["available"], "error": record.get("error"),
                          "dated_inside": {k: v["dated_inside"] for k, v in record["runs"].items()}}, indent=1))
    elif args.command == "hardware":
        summary = hardware(args.out, args.key_file, args.prereg, args.after)
        print(json.dumps({k: v for k, v in summary.items() if k != "analysis"}, indent=1))
    else:
        summary = finish(args.out)
        print(json.dumps({k: v for k, v in summary.items() if k != "analysis"}, indent=1))


if __name__ == "__main__":
    main()
