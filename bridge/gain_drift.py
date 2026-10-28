"""
bridge.gain_drift — is the drifting X rotation error a drive gain shared across qubits? (E5)
===========================================================================================

Pre-registered in ``results/drift_checks/E5_PREREGISTRATION.md``; follows E4 on the same backend
and reuses E4's circuits and per-job estimator from ``bridge.pulse_amplitude`` unchanged. A drive
gain g scales every pulse, so x over-rotates by eps = pi g and sx by delta = (pi / 2) g. A gain
common to every qubit moves eps by the same number of radians on all of them, whatever each
started at; a re-tune moves each qubit's eps toward its own zero.

Jobs run one at a time: each is submitted once the previous one has a result and at least
spacing_s after the previous submission, so that a queue cannot pack them together. Each job gets
a calibration snapshot before submission and another when its result arrives.

Analysis: E4's estimates on every qubit, except that the decay gamma is fitted jointly with each
job's values (E4's pooled fit holds eps fixed across jobs, so it mistakes a drift for decay);
exclusions for frequency (H5) and contrast (HC), which leave the primary set P; the gate (eps
varies on enough of P); then on P, H1 (the ratio beta = d delta / d eps is a gain's, and the same
on every qubit), H2 (the eps increments are shared), H3 (every qubit carries the common increment
with the same loading) and H4 (the common component has no steps).
"""
from __future__ import annotations

import argparse
import calendar
import json
import time
from collections import deque
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from osiris_governance import telemetry
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, rankdata

from dnalang.ledger import Ledger

from . import drift_check as dc
from . import pulse_amplitude as pa


@dataclass(frozen=True)
class E5Design:
    jobs: int = 24
    spacing_s: float = 600.0             # at least this long from one submission to the next
    cap_s: float = 8 * 3600.0            # no submission later than this after the first
    drain_s: float = 7200.0              # once submissions stop, wait this long for the job in flight
    poll_s: float = 30.0
    shots: int = 1000
    replicas: int = 2000
    alpha: float = dc.ALPHA
    analysis_seed: int = 2026102901
    control_seed: int = 2026102902
    sim_seed: int = 2026102903
    sim_trials: int = 50
    sim_replicas: int = 1000             # Holm's first bar over 32 tests is 0.05 / 32, below 1 / 501
    extra_qubits: int = 8
    ramsey_delay_s: float = pa.RAMSEY_DELAY_S
    iterations: int = 5                  # Gauss-Newton steps per job
    gate_min: int = 3
    gain_band: Tuple[float, float] = (0.30, 0.75)
    x_only_max: float = 0.25
    detuning_min: float = 1.25
    ramsey_limit_hz: float = 200e3
    min_jobs: int = 16
    long_gap_s: float = 300.0
    min_long_gaps: int = 12
    min_budget_s: float = 150.0
    stop_budget_s: float = 20.0
    r_min_rad: Optional[float] = None    # H3's power condition; set from E5-sim before the freeze
    move_sd: float = 3.0                 # secondary: common moves beyond this many robust SDs
    walk_q: float = pa.E3_WALK_Q         # S1, S2, S5, S6: rad^2 / s
    loadings: Tuple[float, float] = (0.5, 1.5)      # S3
    event_eps_rad: float = 0.005         # S5
    event_hz: float = 340e3              # S5
    e1_rise: Tuple[float, float] = (0.02, 0.10)     # S6: per-job rise of the |1> readout error
    step_rad: float = 0.002              # S7


def code_hashes() -> Dict[str, str]:
    return {**dc.code_hashes(), "bridge/pulse_amplitude.py": dc.sha256_file(Path(pa.__file__)),
            "bridge/gain_drift.py": dc.sha256_file(Path(__file__))}


# --------------------------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------------------------

def dispersion_p(series: np.ndarray, replicas: np.ndarray) -> float:
    """E4's H1 statistic: the predictive p of the series' variance, in micro-units, among the replicas'."""
    observed = telemetry.dispersion(pa._micro(series))
    return float(telemetry.predictive_p_value(observed, [telemetry.dispersion(pa._micro(r)) for r in replicas]))


def frequency_test(shift_hz: Dict[int, np.ndarray], eps: Dict[int, np.ndarray], rng: np.random.Generator,
                   design: E5Design) -> Dict[str, Any]:
    """H5, per qubit: fails if max |frequency shift| reaches the limit, or if the Spearman correlation
    of the frequency increments with the eps increments is Holm-flagged (two-sided, permutation
    null). Increments, not levels: a slow trend in both series correlates their levels (G6)."""
    p, max_hz, rho = {}, {}, {}
    for q in sorted(shift_hz):
        rf, re = rankdata(np.diff(shift_hz[q])), rankdata(np.diff(eps[q]))
        rf, re = rf - rf.mean(), re - re.mean()
        norm = np.sqrt((rf * rf).sum() * (re * re).sum())
        rho[q] = float((rf * re).sum() / norm) if norm else 0.0
        perms = rng.permuted(np.broadcast_to(rf, (design.replicas, rf.size)), axis=1)
        null = np.abs(perms @ re) / norm if norm else np.zeros(design.replicas)
        p[q] = float((1 + (null >= abs(rho[q]) - 1e-12).sum()) / (1 + design.replicas))
        max_hz[q] = float(np.abs(shift_hz[q]).max())
    flagged = dc.holm({str(q): v for q, v in p.items()}, design.alpha)
    failed = sorted(q for q in p if max_hz[q] >= design.ramsey_limit_hz or flagged[str(q)])
    return {"max_abs_hz": max_hz, "spearman_increments": rho, "p": p, "failed": failed}


def cochran(h1: Dict[str, Any], alpha: float) -> Dict[str, Any]:
    """H1's homogeneity: Cochran's Q of the per-qubit betas about the pooled one, weighted by their
    inverse bootstrap variances, on (qubits - 1) degrees of freedom."""
    weights = h1.get("weights") or {}
    if len(weights) < 2:
        return {"cochran_q": None, "cochran_df": 0, "cochran_p": 1.0, "homogeneous": True}
    stat = sum(w * (h1["per_qubit"][q] - h1["pooled"]) ** 2 for q, w in weights.items())
    p = float(chi2.sf(stat, len(weights) - 1))
    return {"cochran_q": float(stat), "cochran_df": len(weights) - 1, "cochran_p": p, "homogeneous": bool(p > alpha)}


def sharing_test(d: np.ndarray, rng: np.random.Generator, draws: int) -> Dict[str, Any]:
    """H2: mean pairwise uncentred correlation of the qubits' increments (rows of d), against each row
    independently sign-flipped and circularly shifted. Uncentred, so that a trend common to every
    qubit counts as shared; E4's centred test could not see one (G7)."""
    qn, n = d.shape
    if qn < 2:
        return {"mean_correlation": None, "p": 1.0}
    rms = np.sqrt((d * d).mean(-1, keepdims=True))
    z = np.divide(d, rms, out=np.zeros_like(d), where=rms > 0)

    def mean_corr(rows: np.ndarray) -> np.ndarray:
        return ((rows.sum(-2) ** 2).sum(-1) - (rows ** 2).sum((-1, -2))) / (qn * (qn - 1) * n)

    stat = float(mean_corr(z))
    signs = rng.choice((-1.0, 1.0), size=(draws, qn, 1))
    shifts = rng.integers(0, n, size=(draws, qn))
    index = (np.arange(n)[None, None, :] + shifts[:, :, None]) % n
    null = mean_corr(signs * z[np.arange(qn)[None, :, None], index])
    return {"mean_correlation": stat, "p": float((1 + (null >= stat - 1e-12).sum()) / (1 + draws))}


def _weighted_mean(d: np.ndarray, w: np.ndarray) -> np.ndarray:
    """The w-weighted mean over the qubit axis (-2) of d (..., qubits, intervals), kept as an axis."""
    return (d * w[:, None]).sum(-2, keepdims=True) / w.sum()


def _spread(d: np.ndarray, w: np.ndarray) -> np.ndarray:
    return (w[:, None] * (d - _weighted_mean(d, w)) ** 2).sum((-1, -2))


def equality_test(d: np.ndarray, rep: np.ndarray, v: np.ndarray) -> Dict[str, Any]:
    """Secondary (H3 as first drafted): T = sum over intervals k and qubits q of (d_qk - dbar_k)^2 / v_q,
    dbar the precision-weighted mean increment, against T on the no-drift replicas (rep: replicas,
    qubits, intervals). Anything common to every qubit cancels in T."""
    w = 1.0 / v
    stat = float(_spread(d, w))
    null = _spread(rep, w)
    return {"statistic": stat, "null_mean": float(null.mean()),
            "p": float((1 + (null >= stat).sum()) / (1 + len(null)))}


def _loading_statistic(d: np.ndarray, w: np.ndarray) -> np.ndarray:
    mean = _weighted_mean(d, w)
    projection = ((d - mean) * mean).sum(-1)
    return (w * projection ** 2).sum(-1) / (mean ** 2).sum((-1, -2))


def loading_test(d: np.ndarray, rep: np.ndarray, v: np.ndarray) -> Dict[str, Any]:
    """H3: each qubit's deviation from the common increment, projected on the common increment,
    against the replicas' deviations projected on the same common increment: the null is 'fully
    common, along this path'. (Replicas with the path added keep their own mean as well, which
    roughens the direction and made the test conservative.)"""
    w = 1.0 / v
    stat = float(_loading_statistic(d, w))
    null = _loading_statistic(rep - _weighted_mean(rep, w) + _weighted_mean(d, w), w)
    return {"statistic": stat, "p": float((1 + (null >= stat).sum()) / (1 + len(null)))}


def loadings(d: np.ndarray, v: np.ndarray) -> List[Optional[float]]:
    """Each qubit's increments regressed, through 0, on the weighted mean of the other qubits'."""
    w = 1.0 / v
    total = (d * w[:, None]).sum(0)
    out = []
    for q in range(len(d)):
        rest = (total - w[q] * d[q]) / (w.sum() - w[q])
        norm = float(rest @ rest)
        out.append(float(d[q] @ rest) / norm if norm > 0 else None)
    return out


def common_series(levels: np.ndarray, rep_levels: np.ndarray, w: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """The w-weighted mean of the qubits' eps (levels: qubits, jobs) less its first value, and each
    job's variance of that mean among the replicas (rep_levels: replicas, qubits, jobs)."""
    x = (levels * w[:, None]).sum(0) / w.sum()
    xr = (rep_levels * w[:, None]).sum(1) / w.sum()
    return x - x[0], xr.var(0, ddof=1)


def net_change(eps: Dict[int, np.ndarray], rng: np.random.Generator, draws: int) -> Dict[str, Any]:
    """Secondary: the least-squares slope, across qubits, of eps's net change (the last quarter's mean
    less the first quarter's) on its start (the first quarter's mean), with a 95% bootstrap interval
    over qubits. Re-tunes toward zero predict a negative slope; a common gain, zero."""
    rows = np.array([eps[q] for q in sorted(eps)])
    quarter = max(1, rows.shape[1] // 4)
    start = rows[:, :quarter].mean(1)
    change = rows[:, -quarter:].mean(1) - start

    def slope(x, y):
        xc = x - x.mean(-1, keepdims=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            return (xc * (y - y.mean(-1, keepdims=True))).sum(-1) / (xc * xc).sum(-1)

    pick = rng.integers(0, len(start), size=(draws, len(start)))
    boot = slope(start[pick], change[pick])
    boot = boot[np.isfinite(boot)]
    return {"start": dict(zip(sorted(eps), map(float, start))), "change": dict(zip(sorted(eps), map(float, change))),
            "slope": float(slope(start, change)), "interval": [float(v) for v in np.percentile(boot, [2.5, 97.5])]}


def large_moves(common: np.ndarray, k_sd: float) -> List[int]:
    """Secondary: intervals where the common series moves by more than k_sd robust SDs (1.4826 MAD)
    of its own increments away from their median."""
    m = np.diff(common)
    centre = np.median(m)
    sd = 1.4826 * np.median(np.abs(m - centre))
    return [int(k) for k in np.flatnonzero(np.abs(m - centre) > k_sd * sd)] if sd > 0 else []


def _distances(edges: Sequence[Sequence[int]], qubits: Sequence[int]) -> Dict[Tuple[int, int], int]:
    """Shortest-path distances on the coupling map between every pair of the qubits."""
    graph: Dict[int, set] = {}
    for a, b in edges:
        graph.setdefault(a, set()).add(b)
        graph.setdefault(b, set()).add(a)
    out = {}
    for source in qubits:
        seen, queue = {source: 0}, deque([source])
        while queue:
            node = queue.popleft()
            for nxt in graph.get(node, ()):
                if nxt not in seen:
                    seen[nxt] = seen[node] + 1
                    queue.append(nxt)
        out.update({(source, q): seen[q] for q in qubits if q in seen and q != source})
    return out


def distance_check(residuals: Dict[int, np.ndarray], edges: Sequence[Sequence[int]]) -> Dict[str, Any]:
    """Secondary: Spearman correlation, over pairs of qubits, between the correlation of their residual
    increments (each less the common increment) and their distance on the coupling map."""
    qubits = sorted(residuals)
    dist = _distances(edges, qubits)
    pairs = [(a, b) for i, a in enumerate(qubits) for b in qubits[i + 1:] if (a, b) in dist]
    if len(pairs) < 3:
        return {"pairs": len(pairs), "spearman": None}
    corr = [float(np.corrcoef(residuals[a], residuals[b])[0, 1]) for a, b in pairs]
    rc, rd = rankdata(corr), rankdata([dist[p] for p in pairs])
    rc, rd = rc - rc.mean(), rd - rd.mean()
    norm = np.sqrt((rc * rc).sum() * (rd * rd).sum())
    return {"pairs": len(pairs), "spearman": float((rc * rd).sum() / norm) if norm else None,
            "by_pair": [{"qubits": list(p), "distance": dist[p], "correlation": c} for p, c in zip(pairs, corr)]}


READINGS = {
    "not_evaluable": "not evaluable: no variation in E5's window (not evidence of absence)",
    "frequency": "frequency, not gain",
    "not_gain": "not a gain change",
    "gain_heterogeneous": "gain-like on average, but beta differs between qubits: re-tunes or mixed mechanisms",
    "independent": "independent per-qubit gain drift: falsifies sharing",
    "partly_shared": "partly shared gain: a common component plus qubit-specific changes",
    "shared_untested": "shared gain; whether it is fully common is untested (H3 underpowered)",
    "fully_common": "fully common gain drift: supports the hypothesis",
    "fully_common_steps": "fully common gain, changing in steps: a shared change, not a drift",
}
REGIONS = {"x_only": "the ratio reads X-only", "detuning_like": "the ratio is detuning-like",
           "inconclusive": "the ratio is inconclusive"}


def category(result: Dict[str, Any]) -> str:
    """The pre-registered reading table, applied in order; returns its row's key."""
    if not result["gate"]:
        return "not_evaluable"
    if 2 * len(result["excluded"]["H5"]) > result["n_qubits"]:
        return "frequency"
    if result["H1"]["region"] != "gain":
        return "not_gain"
    if not result["H1"]["homogeneous"]:
        return "gain_heterogeneous"
    if not result["H2"]["flagged"]:
        return "independent"
    if result["H3"]["rejected"]:
        return "partly_shared"
    if result["H3"]["underpowered"]:
        return "shared_untested"
    return "fully_common_steps" if result["H4"]["flagged"] else "fully_common"


def reading(result: Dict[str, Any]) -> str:
    key = category(result)
    text = READINGS[key]
    if key == "not_gain":
        text += f": {REGIONS[result['H1']['region']]}"
    for test in ("H5", "HC"):
        if result["excluded"][test]:
            text += f"; excluded by {test}: qubits {result['excluded'][test]}"
    return text


# --------------------------------------------------------------------------------------------
# Analysis
# --------------------------------------------------------------------------------------------

def _weights(freq: np.ndarray, shots: int) -> np.ndarray:
    """pa.per_job's weights: inverse binomial variances at the run's mean frequencies, smoothed by
    half a count."""
    n = freq.shape[-2] * shots
    smooth = (freq.mean(-2, keepdims=True) * n + 0.5) / (n + 1)
    return shots / (smooth * (1 - smooth))


def fit_joint(fitted: np.ndarray, shots: int, design: E5Design) -> Tuple[np.ndarray, float]:
    """E4's pooled fit with gamma re-estimated jointly with each job's (c, a, delta, eps): the gamma
    that minimises the per-job fits' weighted squared error, with the per-job values refitted at
    every trial gamma (a profile over gamma). Holding eps fixed across jobs, E4's pooled fit takes
    a drift for decay; this does not. Returns the parameters and the pooled fit's own gamma."""
    freq = fitted / shots
    theta = pa.fit_pooled(fitted, shots)
    weight = _weights(freq, shots)

    def cost(gamma: float) -> float:
        start = np.r_[theta[:2], gamma, theta[3:]]
        phi = pa.per_job(freq, start, shots, design.iterations)
        return float((weight * (freq - pa.p_one(phi[:, 0], phi[:, 1], gamma, phi[:, 2], phi[:, 3])) ** 2).sum())

    gamma = float(minimize_scalar(cost, bounds=pa.BOUNDS[2], method="bounded", options={"xatol": 1e-9}).x)
    if cost(0.0) <= cost(gamma):         # the bounded search never lands on the bound itself
        gamma = 0.0
    return np.r_[theta[:2], gamma, theta[3:]], float(theta[2])


def estimate(y: np.ndarray, shots: int, design: E5Design, delay_s: float, rng: np.random.Generator) -> Dict[str, Any]:
    """pa.estimate with fit_joint's gamma: per-job estimates, Ramsey shifts and null replicas, which
    are binomial draws at the run's mean frequencies through the same per-job fit."""
    fitted = y[:, :len(pa.FITTED)]
    freq = fitted / shots
    theta, pooled_gamma = fit_joint(fitted, shots, design)
    phi = pa.per_job(freq, theta, shots, design.iterations)
    phase = pa.ramsey_phase(y[:, len(pa.FITTED)] / shots, y[:, len(pa.FITTED) + 1] / shots, phi[:, 0])
    replicas = rng.binomial(shots, freq.mean(0), size=(design.replicas,) + freq.shape) / shots
    rep = pa.per_job(replicas, theta, shots, design.iterations)
    model = pa.p_one(*phi.T[:2], theta[2], *phi.T[2:])
    z = (freq - model) / np.sqrt(model * (1 - model) / shots)
    return {"theta": theta, "pooled_gamma": pooled_gamma, "phi": phi, "shift_hz": pa.frequency_shift(phase, delay_s),
            "replicas": rep, "var": rep.var(0, ddof=1), "residual_z": z}


def _increments(est: Dict[int, Dict[str, Any]], qubits: Sequence[int]):
    """eps levels (qubits, jobs) and their replicas (replicas, qubits, jobs); the increments of both;
    and each qubit's increment variance among the replicas, pooled over intervals (replica jobs are
    identically distributed)."""
    levels = np.array([est[q]["phi"][:, 3] for q in qubits])
    rep_levels = np.stack([est[q]["replicas"][:, :, 3] for q in qubits], axis=1)
    d, rep = np.diff(levels, axis=-1), np.diff(rep_levels, axis=-1)
    return levels, rep_levels, d, rep, rep.var(0, ddof=1).mean(-1)


def confirmatory(t: np.ndarray, est: Dict[int, Dict[str, Any]], qubits: Sequence[int], rng: np.random.Generator,
                 design: E5Design) -> Dict[str, Any]:
    """H1-H4 on the given qubits."""
    if len(qubits) < 2:
        return {"H1": {"region": "inconclusive", "homogeneous": True, "pooled": None, "interval": None},
                "H2": {"p": 1.0, "flagged": False}, "H4": {"p": 1.0, "flagged": False},
                "H3": {"p": 1.0, "rejected": False, "common_range": 0.0, "underpowered": True}}
    eps = {q: est[q]["phi"][:, 3] for q in qubits}
    h1 = pa.ratio_test(eps, {q: est[q]["phi"][:, 2] for q in qubits},
                       {q: float(est[q]["var"][:, 3].mean()) for q in qubits},
                       {q: float(est[q]["var"][:, 2].mean()) for q in qubits}, rng, design)
    h1.update(cochran(h1, design.alpha))
    levels, rep_levels, d, rep, v = _increments(est, qubits)
    h2 = sharing_test(d, rng, design.replicas)
    h2["flagged"] = bool(h2["p"] <= design.alpha)
    h3 = loading_test(d, rep, v)
    common, common_var = common_series(levels, rep_levels, 1.0 / v)
    span = float(np.ptp(common))
    h3.update({"rejected": bool(h3["p"] <= design.alpha), "common_range": span, "r_min": design.r_min_rad,
               "underpowered": design.r_min_rad is None or span < design.r_min_rad,
               "loadings": dict(zip(qubits, loadings(d, v))), "spread_test": equality_test(d, rep, v)})
    h4 = pa.jump_test(common, t, common_var, rng, design.replicas)
    h4["flagged"] = bool(h4["p"] <= design.alpha)
    return {"H1": h1, "H2": h2, "H3": h3, "H4": h4,
            "common": {"series": common.tolist(), "variance": common_var.tolist(),
                       "increment_variance": dict(zip(qubits, map(float, v)))}}


def secondary_checks(t: np.ndarray, est: Dict[int, Dict[str, Any]], result: Dict[str, Any], rng: np.random.Generator,
                     design: E5Design, edges: Optional[Sequence[Sequence[int]]] = None) -> Dict[str, Any]:
    """Reported, not judged: per-qubit jumps, net change against start, large common moves, residual
    variance ratios (E4's G2 check), distance on the coupling map, and H1-H4 on every qubit."""
    primary = result["primary"]
    out: Dict[str, Any] = {}
    jumps = {q: pa.jump_test(est[q]["phi"][:, 3], t, est[q]["var"][:, 3], rng, design.replicas) for q in primary}
    flagged = dc.holm({str(q): j["p"] for q, j in jumps.items()}, design.alpha)
    for q, j in jumps.items():
        j["flagged"] = flagged[str(q)]
    out["per_qubit_jumps"] = jumps
    if len(primary) >= 3:
        out["net_change"] = net_change({q: est[q]["phi"][:, 3] for q in primary}, rng, design.replicas)
    if len(primary) >= 2:
        _, _, d, _, v = _increments(est, primary)
        m = _weighted_mean(d, 1.0 / v)[0]
        out["large_common_moves"] = large_moves(np.array(result["common"]["series"]), design.move_sd)
        out["residual_variance_ratio"] = {q: float((d[i] - m).var(ddof=1) / v[i]) for i, q in enumerate(primary)}
        out["common_increment_variance_ratio"] = float(m.var(ddof=1) * (1.0 / v).sum())
        if edges:
            out["distance"] = distance_check({q: d[i] - m for i, q in enumerate(primary)}, edges)
    every = confirmatory(t, est, result["analysed"], rng, design)
    out["all_qubits"] = {k: v for k, v in every.items() if k != "common"}
    out["all_qubits"]["category"] = category({**result, **out["all_qubits"]})
    return out


def analyze(times: Sequence[float], y_by_qubit: Dict[int, np.ndarray], shots: int, design: E5Design, seed: int,
            delay_s: float, gate: bool = True, x100_by_qubit: Optional[Dict[int, Sequence[int]]] = None,
            secondary: bool = True, edges: Optional[Sequence[Sequence[int]]] = None) -> Dict[str, Any]:
    """Exclusions, the gate, then (if it holds, or always with gate=False) H1-H4 on the primary set."""
    t = np.asarray(times, float)
    qubits = sorted(y_by_qubit)
    streams = np.random.SeedSequence(seed).spawn(len(qubits) + 1)    # replicas per qubit; then the tests
    est = {q: estimate(np.asarray(y_by_qubit[q]), shots, design, delay_s, np.random.default_rng(s))
           for q, s in zip(qubits, streams)}
    rng = np.random.default_rng(streams[-1])
    p_disp, records = {}, {}
    for q in qubits:
        e = est[q]
        record: Dict[str, Any] = {
            "pooled": dict(zip(("c", "a", "gamma", "delta", "eps"), map(float, e["theta"]))),
            "pooled_fit_gamma": e["pooled_gamma"],
            "per_job": {name: [float(v) for v in e["phi"][:, i]] for i, name in enumerate(pa.PER_JOB)},
            "frequency_shift_hz": [float(v) for v in e["shift_hz"]],
            "measurement_sd": {name: float(np.sqrt(e["var"][:, i].mean())) for i, name in enumerate(pa.PER_JOB)},
            "residual_z": {"mean": dict(zip(pa.FITTED, map(float, e["residual_z"].mean(0)))),
                           "sd": dict(zip(pa.FITTED, map(float, e["residual_z"].std(0, ddof=1))))},
            "p": {},
        }
        for name, column in (("rotation_error", 3), ("contrast", 1)):
            p = dispersion_p(e["phi"][:, column], e["replicas"][:, :, column])
            record["p"][name] = {"dispersion": p}
            p_disp[f"q{q}/{name}/dispersion"] = p
        if x100_by_qubit is not None:
            record["x100_errors"] = [int(v) for v in x100_by_qubit[q]]
        records[str(q)] = record
    rejected = dc.holm(p_disp, design.alpha)
    h5 = frequency_test({q: est[q]["shift_hz"] for q in qubits}, {q: est[q]["phi"][:, 3] for q in qubits}, rng, design)
    hc = [q for q in qubits if rejected[f"q{q}/contrast/dispersion"]]
    primary = [q for q in qubits if q not in h5["failed"] and q not in hc]
    varying = [q for q in primary if rejected[f"q{q}/rotation_error/dispersion"]]
    result: Dict[str, Any] = {
        "qubits": records, "n_qubits": len(qubits), "analysed": qubits, "dispersion_tests": len(p_disp),
        "alpha": design.alpha, "replicas": design.replicas, "seed": seed,
        "flagged": sorted(k for k, r in rejected.items() if r), "H5": h5,
        "excluded": {"H5": h5["failed"], "HC": hc}, "primary": primary, "varying": varying,
        "gate": len(varying) >= design.gate_min,
    }
    if gate and not result["gate"]:
        result["reading"] = reading(result)
        return result
    result.update(confirmatory(t, est, primary, rng, design))
    if secondary:
        result["secondary"] = secondary_checks(t, est, result, rng, design, edges)
    result["reading"] = reading(result)
    return result


def control_flags(result: Dict[str, Any]) -> List[str]:
    """Everything the constant control flagged: dispersions, H5 failures, H1 reading gain or
    heterogeneous, H2, H3 and H4."""
    flags = list(result["flagged"]) + [f"q{q}/H5" for q in result["excluded"]["H5"]]
    if result["H1"]["region"] == "gain":
        flags.append("H1/gain")
    if not result["H1"]["homogeneous"]:
        flags.append("H1/heterogeneous")
    return flags + [name for name, key in (("H2", "flagged"), ("H3", "rejected"), ("H4", "flagged")) if result[name][key]]


# --------------------------------------------------------------------------------------------
# Software control and simulation stage (E5-sim)
# --------------------------------------------------------------------------------------------

def synthetic_counts(i: int, eps: np.ndarray, delta: np.ndarray, shift_hz: np.ndarray, design: E5Design,
                     rng: np.random.Generator, e1_rise: Any = 0.0) -> np.ndarray:
    """pa.synthetic_counts for any number of jobs, with an optional per-job rise of the |1> readout
    error e1, which lowers the contrast a by the rise and the midpoint c by half of it."""
    c0, a0, gamma, _, _ = pa.sim_theta(i)
    jobs = np.asarray(eps).size
    eps, delta, shift_hz, e1_rise = (np.broadcast_to(np.asarray(v, float), (jobs,))
                                     for v in (eps, delta, shift_hz, e1_rise))
    c, a = c0 - e1_rise / 2, a0 - e1_rise
    ideal = pa.p_one(c, a, gamma, delta, eps)
    p = np.clip(ideal + pa.sim_misfit(i) * 4 * ideal * (1 - ideal), 1e-6, 1 - 1e-6)
    phase = 2 * np.pi * shift_hz * design.ramsey_delay_s
    ram = np.stack([c + a / 2 * np.cos(phase), c - a / 2 * np.sin(phase)], axis=-1)
    x100 = (c - a / 2 + a * np.sin(50 * eps) ** 2)[:, None]
    return rng.binomial(design.shots, np.concatenate([p, ram, x100], axis=-1))


def software_control(design: E5Design, qubits: int = 16) -> Dict[str, Any]:
    """S0 at count level: constant parameters with E3's static residuals; every test, ungated."""
    rng = np.random.default_rng(design.control_seed)
    k = design.jobs
    y = {i: synthetic_counts(i, np.full(k, pa.sim_theta(i)[4]), np.full(k, pa.sim_theta(i)[3]), np.zeros(k),
                             design, rng) for i in range(qubits)}
    result = analyze(design.spacing_s * np.arange(k), y, design.shots, design, design.control_seed,
                     design.ramsey_delay_s, gate=False, secondary=False)
    result["control_flags"] = control_flags(result)
    result["dispatch_allowed"] = not result["control_flags"]
    return result


SCENARIOS = ("S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7")
S3_WALK = (0.003, 0.004, 0.005, 0.006, 0.008)       # common ranges, rad; R_min is read from these
S3_LINEAR = (0.008, 0.012, 0.016)                   # reported: the least favourable path


def sim_truth(scenario: str, design: E5Design, seeds: Sequence[np.random.SeedSequence], qubits: int = 16,
              common_range: Optional[float] = None, path: str = "walk"):
    """Times, each sim qubit's true per-job (eps, delta, frequency shift, e1 rise), and the event
    qubit (S5, S6) for one trial. seeds: (walk, event, per-scenario) streams; one trial index gives
    S1, S3, S5 and S6 the same common walk, so S5 and S6 can be compared with S1 trial by trial."""
    walk_rng, event_rng, own_rng = (np.random.default_rng(s) for s in seeds)
    k = design.jobs
    t = design.spacing_s * np.arange(k, dtype=float)
    sd = np.sqrt(design.walk_q * np.diff(t))
    common = np.concatenate([[0.0], np.cumsum(walk_rng.normal(0, sd))])
    eps = {i: pa.sim_theta(i)[4] + np.zeros(k) for i in range(qubits)}
    delta = {i: pa.sim_theta(i)[3] + np.zeros(k) for i in range(qubits)}
    shift = {i: np.zeros(k) for i in range(qubits)}
    rise = {i: np.zeros(k) for i in range(qubits)}
    event = None
    if scenario in ("S1", "S5", "S6"):
        for i in range(qubits):
            eps[i] += common
            delta[i] += common / 2
    if scenario == "S2":
        for i in range(qubits):
            own = np.concatenate([[0.0], np.cumsum(own_rng.normal(0, sd))])
            eps[i] += own
            delta[i] += own / 2
    if scenario == "S3":
        line = common * common_range / np.ptp(common) if path == "walk" else common_range * t / t[-1]
        load = own_rng.uniform(*design.loadings, qubits)
        for i in range(qubits):
            eps[i] += load[i] * line
            delta[i] += load[i] * line / 2
    if scenario == "S4":                  # 1-2 re-tunes per qubit, each removing 50-100% of eps
        for i in range(qubits):
            for j in sorted(own_rng.choice(np.arange(1, k), int(own_rng.integers(1, 3)), replace=False)):
                change = -own_rng.uniform(0.5, 1.0) * eps[i][j]
                eps[i][j:] += change
                delta[i][j:] += own_rng.uniform(0.0, 1.0) * change
    if scenario in ("S5", "S6"):
        event = int(event_rng.integers(qubits))
        late = np.arange(k) >= event_rng.integers(1, k)
        if scenario == "S5":
            eps[event] += design.event_eps_rad * late
            shift[event] += design.event_hz * late
        else:
            rise[event] = late * event_rng.uniform(*design.e1_rise, k)
    if scenario == "S7":
        late = np.arange(k) >= event_rng.integers(1, k)
        for i in range(qubits):
            eps[i] += design.step_rad * late
            delta[i] += design.step_rad / 2 * late
    return t, {i: (eps[i], delta[i], shift[i], rise[i]) for i in range(qubits)}, event


def compact(result: Dict[str, Any], event: Optional[int] = None) -> Dict[str, Any]:
    """The parts of one analysis that the E5-sim rates and the reading need."""
    keep = {"H1": ("region", "homogeneous", "pooled", "interval", "cochran_p"),
            "H2": ("mean_correlation", "p", "flagged"), "H3": ("p", "rejected", "common_range", "underpowered"),
            "H4": ("p", "flagged")}
    row = {"gate": result["gate"], "n_qubits": result["n_qubits"], "excluded": result["excluded"],
           "flagged": result["flagged"], "event_qubit": event}
    for test, keys in keep.items():
        row[test] = {k: result[test].get(k) for k in keys}
    row["H3"]["spread_p"] = (result["H3"].get("spread_test") or {}).get("p")
    row["category"] = category(row)
    return row


def sim_trial(scenario: str, trial: int, design: E5Design, qubits: int = 16, common_range: Optional[float] = None,
              path: str = "walk") -> Dict[str, Any]:
    """One count-level trial: truth, counts, the full analysis (ungated, no secondary checks)."""
    children = np.random.SeedSequence([design.sim_seed, trial]).spawn(4 + qubits)
    t, truth, event = sim_truth(scenario, design, children[:3], qubits, common_range, path)
    y = {i: synthetic_counts(i, *truth[i][:3], design, np.random.default_rng(children[4 + i]), truth[i][3])
         for i in truth}
    seed = int(children[3].generate_state(1)[0])
    result = analyze(t, y, design.shots, design, seed, design.ramsey_delay_s, gate=False, secondary=False)
    return {"trial": trial, **compact(result, event)}


def sim_key(scenario: str, common_range: Optional[float] = None, path: str = "walk") -> str:
    return f"S3_{path}_{common_range * 1e3:g}mrad" if scenario == "S3" else scenario


def simulate(out: Path, scenario: str, trials: Tuple[int, int], design: E5Design = E5Design(),
             common_range: Optional[float] = None, path: str = "walk") -> Dict[str, Any]:
    """One E5-sim invocation: S0 (the software control, at the hardware's replica count) or trials
    [start, stop) of one scenario at sim_replicas. Several invocations keep each under the horizon."""
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    key = sim_key(scenario, common_range, path)
    record: Dict[str, Any] = {"key": key, "scenario": scenario, "common_range": common_range, "path": path,
                              "code_sha256": dc.sha256_file(Path(__file__))}
    if scenario == "S0":
        control = software_control(design)
        record["control"] = {k: control[k] for k in ("control_flags", "dispatch_allowed", "flagged", "excluded")}
        record["control"]["category"] = category(control)
        record["control"]["tests"] = {name: {k: v for k, v in control[name].items()
                                             if k in ("p", "flagged", "rejected", "region", "interval", "cochran_p",
                                                      "homogeneous", "common_range")}
                                      for name in ("H1", "H2", "H3", "H4")}
        name = "S0.json"
    else:
        sim = replace(design, replicas=design.sim_replicas)
        record["trials"] = list(trials)
        record["rows"] = [sim_trial(scenario, n, sim, common_range=common_range, path=path) for n in range(*trials)]
        name = f"{key}_{trials[0]}-{trials[1]}.json"
    record["seconds"] = round(time.monotonic() - started, 1)
    (out / name).write_text(json.dumps(record, indent=1, default=float) + "\n")
    return record


def _rate(rows: List[Dict[str, Any]], test) -> Optional[float]:
    return float(np.mean([bool(test(r)) for r in rows])) if rows else None


def r_min_from(power: Dict[float, float], target: float = 0.80) -> Optional[float]:
    """The smallest grid range whose power, and every larger range's, reaches the target."""
    ranges = sorted(power)
    for i, r in enumerate(ranges):
        if all(power[s] >= target for s in ranges[i:]):
            return r
    return None


def summarize(out: Path, design: E5Design = E5Design()) -> Dict[str, Any]:
    """Merges the E5-sim invocations under out, sets R_min from S3's walk grid, applies it to every
    trial's reading, and writes oc.json."""
    rows: Dict[str, Dict[int, Dict[str, Any]]] = {}
    seconds = 0.0
    control = None
    codes = set()
    for path in sorted(out.glob("S*.json")):
        record = json.loads(path.read_text())
        codes.add(record["code_sha256"])
        seconds += record["seconds"]
        if record["scenario"] == "S0":
            control = record["control"]
            continue
        rows.setdefault(record["key"], {}).update({r["trial"]: r for r in record["rows"]})
    if len(codes) > 1:
        raise RuntimeError(f"{out} mixes invocations from {len(codes)} code versions")
    by_key = {key: [found[n] for n in sorted(found)] for key, found in rows.items()}
    walk = {r: _rate(by_key[sim_key("S3", r)], lambda x: x["H3"]["rejected"])
            for r in S3_WALK if sim_key("S3", r) in by_key}
    r_min = r_min_from(walk) if walk else None
    for found in by_key.values():
        for row in found:
            row["H3"]["underpowered"] = r_min is None or row["H3"]["common_range"] < r_min
            row["category"] = category(row)
    rates: Dict[str, Dict[str, Any]] = {}
    for key, found in by_key.items():
        rates[key] = {"trials": len(found), "gate": _rate(found, lambda x: x["gate"]),
                      "H1_gain": _rate(found, lambda x: x["H1"]["region"] == "gain"),
                      "H1_homogeneous": _rate(found, lambda x: x["H1"]["homogeneous"]),
                      "H2": _rate(found, lambda x: x["H2"]["flagged"]),
                      "H3_rejected": _rate(found, lambda x: x["H3"]["rejected"]),
                      "H3_spread_rejected": _rate(found, lambda x: (x["H3"]["spread_p"] or 1.0) <= design.alpha),
                      "H4": _rate(found, lambda x: x["H4"]["flagged"]),
                      "median_common_range": float(np.median([x["H3"]["common_range"] for x in found])),
                      "categories": {c: sum(x["category"] == c for x in found) for c in READINGS}}
    paired = rows.get("S1", {})
    for key, test in (("S5", "H5"), ("S6", "HC")):
        if key not in rates:
            continue
        found = by_key[key]
        rates[key]["event_excluded"] = _rate(found, lambda x: x["event_qubit"] in x["excluded"][test])
        both = [x for x in found if x["trial"] in paired]
        rates[key]["reading_unchanged"] = _rate(both, lambda x: x["category"] == paired[x["trial"]]["category"])
    common = ("fully_common", "fully_common_steps")
    if "S4" in rates:
        rates["S4"]["not_fully_common"] = _rate(by_key["S4"], lambda x: x["category"] not in common)

    def at_least(key, name, bound):
        value = (rates.get(key) or {}).get(name)
        return None if value is None else bool(value >= bound)

    def at_most(key, name, bound):
        value = (rates.get(key) or {}).get(name)
        return None if value is None else bool(value <= bound)

    verdict = {
        "S0_dispatch_allowed": (control or {}).get("dispatch_allowed"),
        "S1_H2_power_ge_0.80": at_least("S1", "H2", 0.80),
        "S1_H3_false_rejection_le_0.10": at_most("S1", "H3_rejected", 0.10),
        "S1_H1_gain_ge_0.90": at_least("S1", "H1_gain", 0.90),
        "S1_H4_false_alarm_le_0.10": at_most("S1", "H4", 0.10),
        "S2_H2_false_alarm_le_0.10": at_most("S2", "H2", 0.10),
        "S4_not_fully_common_ge_0.90": at_least("S4", "not_fully_common", 0.90),
        "S5_H5_excludes_ge_0.80": at_least("S5", "event_excluded", 0.80),
        "S5_reading_unchanged_ge_0.80": at_least("S5", "reading_unchanged", 0.80),
        "S6_HC_excludes_ge_0.80": at_least("S6", "event_excluded", 0.80),
        "S6_reading_unchanged_ge_0.80": at_least("S6", "reading_unchanged", 0.80),
        "S7_H4_power_ge_0.80": at_least("S7", "H4", 0.80),
    }
    summary = {"design": asdict(replace(design, r_min_rad=r_min)), "r_min_rad": r_min,
               "S3_power": {"walk": {str(r): p for r, p in walk.items()},
                            "linear": {str(r): _rate(by_key[sim_key("S3", r, "linear")],
                                                     lambda x: x["H3"]["rejected"])
                                       for r in S3_LINEAR if sim_key("S3", r, "linear") in by_key},
                            "spread_walk": {str(r): rates[sim_key("S3", r)]["H3_spread_rejected"]
                                            for r in S3_WALK if sim_key("S3", r) in rates},
                            "spread_linear": {str(r): rates[sim_key("S3", r, "linear")]["H3_spread_rejected"]
                                              for r in S3_LINEAR if sim_key("S3", r, "linear") in rates}},
               "code_versions": sorted(codes),
               "control": control, "rates": rates, "verdict": verdict, "seconds": round(seconds, 1),
               "code_sha256": dc.sha256_file(Path(__file__))}
    (out / "oc.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    return summary


# --------------------------------------------------------------------------------------------
# Hardware
# --------------------------------------------------------------------------------------------

def _clock() -> float:
    return time.time()


def _sleep(seconds: float) -> None:
    time.sleep(seconds)


def _ts_epoch(row: Dict[str, Any]) -> float:
    return float(calendar.timegm(time.strptime(row["ts"], "%Y-%m-%dT%H:%M:%SZ")))


def _save(out: Path, state: Dict[str, Any]) -> None:
    (out / "state.json").write_text(json.dumps(state, indent=1) + "\n")


def _budget(svc) -> Optional[float]:
    try:
        return svc.usage().get("usage_remaining_seconds")
    except Exception:
        return None


def _stop(out: Path, ledger: Ledger, state: Dict[str, Any], reason: str, **extra) -> None:
    state["stopped"], state["stopped_at"] = reason, _clock()
    ledger.append("stopped", {"reason": reason, "next_job": state["next_job"], **extra})
    _save(out, state)
    print(f"submissions stopped: {reason}", flush=True)


def _after_snapshots(svc, ledger: Ledger, state: Dict[str, Any]) -> None:
    """A calibration snapshot for every job that has a result and none yet."""
    results = dc._ledger_results(ledger)
    have = {r["job"] for r in ledger if r.get("kind") == "calibration-snapshot" and r.get("when") == "after"}
    for job in state["jobs"]:
        if job["job_id"] in results and job["job"] not in have:
            props = svc.backend(state["backend"]).properties(refresh=True)
            ledger.append("calibration-snapshot", {"job": job["job"], "when": "after", "job_id": job["job_id"],
                                                   **pa.snapshot(props, state["layout"])})


def _submit_next(out: Path, svc, circuits, state: Dict[str, Any], ledger: Ledger, design: E5Design) -> None:
    """Budget check, calibration snapshot, then ledger-first submission of the next job. An unknown
    budget or a failed snapshot submits nothing (the next poll retries); a failure after the
    submit-intent row records the job as lost and moves on; one before it retries the same job."""
    from dnalang.backends import ibm
    j = state["next_job"]
    remaining = _budget(svc)
    if remaining is not None and remaining < design.stop_budget_s:
        _stop(out, ledger, state, "QPU budget below the stop level", usage_remaining_seconds=remaining)
        return
    if remaining is None:
        print("QPU budget unknown; not submitting", flush=True)
        return
    try:
        props = svc.backend(state["backend"]).properties(refresh=True)
    except Exception as err:            # the type only: messages may echo request details
        print(f"calibration snapshot failed ({type(err).__name__}); not submitting", flush=True)
        return
    ledger.append("calibration-snapshot", {"job": j, "when": "before", "usage_remaining_seconds": remaining,
                                           **pa.snapshot(props, state["layout"])})
    state["last_attempt_at"] = _clock()
    try:
        job_id = ibm.submit(circuits, state["backend"], state["layout"], design.shots, ledger,
                            tags={"experiment": "gain_drift", "job": j})
    except Exception as err:
        last = list(ledger)[-1]
        if last.get("kind") == "submit-intent" and last.get("tags", {}).get("job") == j:
            ledger.append("lost", {"job": j, "intent_hash": last["hash"], "error": type(err).__name__})
            state["lost"].append(j)
            state["next_job"] = j + 1
        print(f"submission {j + 1} failed ({type(err).__name__})", flush=True)
        _save(out, state)
        return
    state["jobs"].append({"job": j, "job_id": job_id, "submitted_at": _clock()})
    state["job_ids"].append(job_id)
    state["next_job"] = j + 1
    if state["first_submitted_at"] is None:
        state["first_submitted_at"] = state["jobs"][-1]["submitted_at"]
    _save(out, state)
    print(f"submitted {j + 1}/{design.jobs}: {job_id}", flush=True)


def run_schedule(out: Path, svc, circuits, state: Dict[str, Any], ledger: Ledger, design: E5Design) -> None:
    """One job in flight. The next job is submitted once every submitted job has a result and
    spacing_s after the last attempt; none later than cap_s after the first. Once submissions
    stop, waits up to drain_s for the job in flight. Polling errors are retried at the next poll."""
    while True:
        try:
            dc.collect(svc, ledger, state["job_ids"])
            _after_snapshots(svc, ledger, state)
        except Exception as err:
            print(f"poll failed ({type(err).__name__}); retrying", flush=True)
        results = dc._ledger_results(ledger)
        pending = [job for job in state["jobs"] if job["job_id"] not in results]
        now = _clock()
        if state["stopped"] is None:
            if state["next_job"] >= design.jobs:
                _stop(out, ledger, state, "all jobs submitted")
            elif state["first_submitted_at"] is not None and now - state["first_submitted_at"] > design.cap_s:
                _stop(out, ledger, state, "wall-clock cap")
        if state["stopped"] is not None:
            if not pending or now - state["stopped_at"] > design.drain_s:
                return
        elif not pending and (state["last_attempt_at"] is None or now - state["last_attempt_at"] >= design.spacing_s):
            _submit_next(out, svc, circuits, state, ledger, design)
        _sleep(design.poll_s)


def reconcile(ledger: Ledger, state: Dict[str, Any]) -> None:
    """After a crash, brings state in line with the ledger: a submitted job missing from state is
    added; a submit-intent without a submitted row is recorded once as lost, never resubmitted."""
    rows = list(ledger)
    submitted = {r["intent_hash"]: r for r in rows if r.get("kind") == "submitted"}
    lost = {r["job"] for r in rows if r.get("kind") == "lost"}
    known = {job["job"] for job in state["jobs"]}
    attempts = [state["last_attempt_at"]] if state["last_attempt_at"] is not None else []
    for row in rows:
        if row.get("kind") != "submit-intent":
            continue
        j = row.get("tags", {}).get("job")
        attempts.append(_ts_epoch(row))
        done = submitted.get(row["hash"])
        if done is not None and j not in known:
            state["jobs"].append({"job": j, "job_id": done["job_id"], "submitted_at": _ts_epoch(done)})
            state["job_ids"].append(done["job_id"])
            known.add(j)
        elif done is None and j not in lost:
            ledger.append("lost", {"job": j, "intent_hash": row["hash"], "error": "no submitted row"})
            lost.add(j)
    state["jobs"].sort(key=lambda job: job["job"])
    state["job_ids"] = [job["job_id"] for job in state["jobs"]]
    state["lost"] = sorted(lost)
    seen = known | lost
    state["next_job"] = max(state["next_job"], max(seen) + 1 if seen else 0)
    if state["jobs"] and state["first_submitted_at"] is None:
        state["first_submitted_at"] = state["jobs"][0]["submitted_at"]
    state["last_attempt_at"] = max(attempts) if attempts else None


def hardware(out: Path, key_file: Path, prereg: Path, after: Path, design: E5Design = E5Design()) -> Dict[str, Any]:
    """Software control, budget and status checks, then up to design.jobs ledger-first jobs, one in
    flight, on E4's backend and qubits plus design.extra_qubits more; then the analysis."""
    from dnalang.backends import ibm
    previous = json.loads(after.read_text())
    out.mkdir(parents=True, exist_ok=True)
    ledger = Ledger(out / "ledger.jsonl")
    if next(iter(ledger), None) is not None:
        raise RuntimeError(f"{out / 'ledger.jsonl'} already has rows; use resume")
    ledger.append("preregistration", {"sha256": dc.sha256_file(prereg), "code_sha256": code_hashes()})
    control = software_control(design, len(previous["layout"]) + design.extra_qubits)
    (out / "control.json").write_text(json.dumps(control, indent=1, default=float) + "\n")
    if not control["dispatch_allowed"]:
        ledger.append("aborted", {"reason": "software control flagged constant data", "flags": control["control_flags"]})
        return {"verdict": "NOT RUN: the constant control flagged", "control": control["control_flags"]}

    dc._use_key_file(key_file)
    svc = ibm.service()
    backend = svc.backend(previous["backend"])
    if not backend.status().operational:
        ledger.append("aborted", {"reason": f"{backend.name} not operational"})
        return {"verdict": f"NOT RUN: {backend.name} not operational"}
    remaining = _budget(svc)
    if remaining is None or remaining < design.min_budget_s:
        ledger.append("aborted", {"reason": "QPU budget below the minimum or unknown",
                                  "usage_remaining_seconds": remaining})
        return {"verdict": f"NOT RUN: {remaining} QPU-seconds remaining, minimum {design.min_budget_s}"}
    layout = list(previous["layout"]) + pa.choose_extras(backend, previous["layout"], design.extra_qubits)
    delay = pa.aligned_delay(backend, design.ramsey_delay_s)
    circuits = pa.pulse_circuits(len(layout), delay)
    edges = sorted([int(a), int(b)] for a, b in backend.coupling_map.get_edges())
    ledger.append("design", {**asdict(design), "backend": backend.name, "layout": layout, "ramsey_delay_s": delay,
                             "usage_remaining_seconds": remaining, "coupling_edges": edges,
                             "circuits": {c.name: c.sha256() for c in circuits}})
    state = {"backend": backend.name, "layout": layout, "ramsey_delay_s": delay, "coupling_edges": edges,
             "jobs": [], "job_ids": [], "lost": [], "next_job": 0, "first_submitted_at": None,
             "last_attempt_at": None, "stopped": None, "stopped_at": None}
    _save(out, state)
    run_schedule(out, svc, circuits, state, ledger, design)
    return finish(out, design)


def resume(out: Path, key_file: Path, design: E5Design = E5Design()) -> Dict[str, Any]:
    """Continues a run from state.json and the ledger, with the same code and circuits."""
    from dnalang.backends import ibm
    ledger = Ledger(out / "ledger.jsonl")
    problem = ledger.verify()
    if problem is not None:
        raise RuntimeError(f"ledger does not verify: {problem}")
    rows = list(ledger)
    if rows[0]["code_sha256"] != code_hashes():
        raise RuntimeError("code differs from the run's preregistration row")
    state = json.loads((out / "state.json").read_text())
    circuits = pa.pulse_circuits(len(state["layout"]), state["ramsey_delay_s"])
    designed = next(r for r in rows if r.get("kind") == "design")
    if {c.name: c.sha256() for c in circuits} != designed["circuits"]:
        raise RuntimeError("circuits differ from the run's design row")
    reconcile(ledger, state)
    _save(out, state)
    ledger.append("resumed", {"next_job": state["next_job"], "lost": state["lost"], "stopped": state["stopped"]})
    dc._use_key_file(key_file)
    run_schedule(out, ibm.service(), circuits, state, ledger, design)
    return finish(out, design)


def finish(out: Path, design: E5Design = E5Design()) -> Dict[str, Any]:
    """Telemetry and analysis from the ledger's result rows. Needs no network."""
    state = json.loads((out / "state.json").read_text())
    ledger = Ledger(out / "ledger.jsonl")
    problem = ledger.verify()
    if problem is not None:
        raise RuntimeError(f"ledger does not verify: {problem}")
    results = dc._ledger_results(ledger)
    done = sorted((results[j] for j in state["job_ids"] if results.get(j, {}).get("status") == "DONE"),
                  key=lambda r: r["window"][0])
    summary: Dict[str, Any] = {"jobs_submitted": len(state["job_ids"]), "jobs_lost": len(state["lost"]),
                               "jobs_done": len(done), "jobs_pending": [j for j in state["job_ids"] if j not in results],
                               "stopped": state["stopped"],
                               "qpu_seconds": sum(r.get("qpu_seconds") or 0 for r in done), "long_intervals": 0}
    layout = state["layout"]
    if len(done) >= 2:
        source = {"kind": "HARDWARE", "node": state["backend"],
                  "generator": f"bridge.gain_drift {dc.sha256_file(Path(__file__))[:12]}"}
        lines = pa.write_telemetry([{"window": r["window"], "counts": r["counts"]} for r in done],
                                   layout, design.shots, f"gain-drift-hw-{state['backend']}", source)
        (out / "telemetry.jsonl").write_bytes(b"".join(line + b"\n" for line in lines))
        times = [float(v) for v in telemetry.metric_series(lines, f"q{layout[0]}_prep0_ones")[0]]
        summary["long_intervals"] = int((np.diff(times) >= design.long_gap_s).sum())
    if len(done) < design.min_jobs or summary["long_intervals"] < design.min_long_gaps:
        summary["verdict"] = (f"INCOMPLETE: {len(done)} jobs finished (minimum {design.min_jobs}) and "
                              f"{summary['long_intervals']} intervals of {design.long_gap_s:g} s or longer "
                              f"(minimum {design.min_long_gaps})")
    else:
        y = {q: np.array([[int(v) for v in telemetry.metric_series(lines, f"q{q}_{name}_ones")[1]]
                          for name in pa.FITTED + pa.RAMSEY]).T for q in layout}
        x100 = {q: [int(v) for v in telemetry.metric_series(lines, f"q{q}_x100_errors")[1]] for q in layout}
        analysis = analyze(times, y, design.shots, design, design.analysis_seed, state["ramsey_delay_s"],
                           x100_by_qubit=x100, edges=state.get("coupling_edges"))
        snaps = [row for row in ledger if row.get("kind") == "calibration-snapshot"]
        windows = [r["window"] for r in done]
        secondary = analysis.get("secondary")
        if secondary:
            secondary["dated_inside_large_moves"] = {
                str(k): pa.dated_inside(snaps, windows[k][1], windows[k + 1][0])
                for k in secondary["large_common_moves"]}
            for q, jump in secondary["per_qubit_jumps"].items():
                if jump["flagged"]:
                    k = jump["interval"]
                    jump["dated_inside"] = pa.dated_inside(snaps, windows[k][1], windows[k + 1][0], [q])
        summary["analysis"] = analysis
    (out / "analysis.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    return summary


def retrospective(out: Path, key_file: Path, run: Path) -> Dict[str, Any]:
    """E5a: reported calibration at each finished job's execution midpoint; with the run's analysis,
    also what was dated inside its large common moves and flagged jumps. Stops at the first failed
    query."""
    from dnalang.backends import ibm
    dc._use_key_file(key_file)
    svc = ibm.service()
    state = json.loads((run / "state.json").read_text())
    results = dc._ledger_results(Ledger(run / "ledger.jsonl"))
    windows = sorted(results[j]["window"] for j in state["job_ids"] if results.get(j, {}).get("status") == "DONE")
    backend = svc.backend(state["backend"])
    record: Dict[str, Any] = {"available": True, "snapshots": [], "code_sha256": dc.sha256_file(Path(__file__))}
    for window in windows:
        try:
            props = backend.properties(datetime=pa._midpoint(window))
        except Exception as err:            # the type only: messages may echo request details
            props, record["error"] = None, type(err).__name__
        if props is None:
            record["available"] = False
            break
        record["snapshots"].append({"at": window, **pa.snapshot(props, state["layout"])})
    found = run / "analysis.json"
    analysis = json.loads(found.read_text()).get("analysis") if found.exists() else None
    if analysis and analysis.get("secondary"):
        secondary = analysis["secondary"]
        intervals = set(secondary["large_common_moves"])
        intervals |= {j["interval"] for j in secondary["per_qubit_jumps"].values() if j["flagged"]}
        if analysis["H4"]["flagged"]:
            intervals.add(analysis["H4"]["interval"])
        record["dated_inside"] = {str(k): pa.dated_inside(record["snapshots"], windows[k][1], windows[k + 1][0])
                                  for k in sorted(intervals)}
    out.mkdir(parents=True, exist_ok=True)
    (out / "retrospective.json").write_text(json.dumps(record, indent=1, default=str) + "\n")
    return record


def _trials(text: str) -> Tuple[int, int]:
    start, stop = (int(v) for v in text.split(":"))
    return start, stop


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(prog="python -m bridge.gain_drift")
    sub = parser.add_subparsers(dest="command", required=True)
    sim = sub.add_parser("simulate")
    sim.add_argument("--out", type=Path, required=True)
    sim.add_argument("--scenario", choices=SCENARIOS, required=True)
    sim.add_argument("--trials", type=_trials, default=(0, E5Design().sim_trials), help="start:stop")
    sim.add_argument("--range", type=float, dest="common_range", help="S3: the common range, rad")
    sim.add_argument("--path", choices=("walk", "linear"), default="walk", help="S3: the common path")
    summ = sub.add_parser("summarize")
    summ.add_argument("--out", type=Path, required=True)
    hw = sub.add_parser("hardware")
    hw.add_argument("--key-file", type=Path, required=True)
    hw.add_argument("--out", type=Path, required=True)
    hw.add_argument("--prereg", type=Path, required=True)
    hw.add_argument("--after", type=Path, required=True)
    res = sub.add_parser("resume")
    res.add_argument("--key-file", type=Path, required=True)
    res.add_argument("--out", type=Path, required=True)
    fin = sub.add_parser("finish")
    fin.add_argument("--out", type=Path, required=True)
    retro = sub.add_parser("retrospective")
    retro.add_argument("--key-file", type=Path, required=True)
    retro.add_argument("--out", type=Path, required=True)
    retro.add_argument("--run", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "simulate":
        if (args.scenario == "S3") != (args.common_range is not None):
            parser.error("--range is required for S3 and only for S3")
        record = simulate(args.out, args.scenario, args.trials, common_range=args.common_range, path=args.path)
        print(json.dumps({k: v for k, v in record.items() if k not in ("rows",)}, indent=1, default=float))
    elif args.command == "summarize":
        summary = summarize(args.out)
        print(json.dumps({k: summary[k] for k in ("r_min_rad", "S3_power", "verdict", "seconds")}, indent=1))
    elif args.command == "retrospective":
        record = retrospective(args.out, args.key_file, args.run)
        print(json.dumps({"available": record["available"], "error": record.get("error"),
                          "dated_inside": record.get("dated_inside")}, indent=1))
    else:
        summary = {"hardware": lambda: hardware(args.out, args.key_file, args.prereg, args.after),
                   "resume": lambda: resume(args.out, args.key_file),
                   "finish": lambda: finish(args.out)}[args.command]()
        print(json.dumps({k: v for k, v in summary.items() if k != "analysis"}, indent=1))


if __name__ == "__main__":
    main()
