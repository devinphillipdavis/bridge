"""
bridge.amplitude_drift — is the job-to-job x100 variation a drifting rotation error? (E3)
=========================================================================================

Pre-registered in ``results/drift_checks/AMPLITUDE_PREREGISTRATION.md``; follows E2 on the same
backend and qubits. Each job runs, on every qubit in parallel:

    prep0  measure                                   P(1) = c - a/2
    prep1  x, measure                                P(1) = c + (a/2) cos(eps)
    ampN   sx, then N x gates, then measure          P(1) = c + (a/2) exp(-gamma N) sin(delta + N eps)
    x100   100 x gates, then measure                 (E2's circuit; continuity only, not fitted)

where every x rotates by pi + eps and sx by pi/2 + delta; c and a are the readout offset and
contrast, gamma the contrast lost per x gate. prep0 and prep1 fix c and a in every job, so the
ampN series, which moves linearly with eps and shows its sign, fixes eps and delta.

Per qubit: a pooled fit of (c, a, gamma, delta, eps) over all jobs; then each job's (eps_k, a_k),
fitted jointly with c, gamma and delta held at the pooled values, so a change in one is not
absorbed by the other. Both series get the osiris_governance.telemetry drift checks. Replicas
come from a parametric bootstrap of the pooled fit, so they have no drift, and pass through the
same refit and per-job estimates. Every p-value in the run is Holm-corrected together.
"""
from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from osiris_governance import telemetry
from osiris_governance.telemetry import TelemetryRun, quantity
from scipy.optimize import minimize
from scipy.stats import spearmanr

from dnalang.ledger import Ledger
from dnalang.qc_ir import Circuit, Op

from . import drift_check as dc

REPS = (0, 4, 8, 16, 32, 48, 64, 100)
FITTED = ("prep0", "prep1") + tuple(f"amp{n}" for n in REPS)
CIRCUIT_NAMES = FITTED + ("x100",)
EPS_RANGE = 0.08            # |eps| searched by the pooled fit, rad per x gate
EPS_WINDOW = 0.02           # per-job eps searched within pooled eps +- this
EPS_STEP = 0.0001           # per-job eps grid, rad
CONTRAST_STEP = 0.005       # per-job contrast grid, over [0, 1.2]
MICRO = 10 ** 6             # estimates enter the exact checks in micro-units
PARAMETERS = ("rotation_error", "contrast")
BOUNDS = [(0.0, 1.0), (0.0, 1.2), (0.0, 0.01), (-1.0, 1.0), (-0.1, 0.1)]
_N = np.asarray(REPS, float)


@dataclass(frozen=True)
class AmpDesign:
    jobs: int = 12
    spacing_s: int = 300
    shots: int = 2000
    replicas: int = 1000
    alpha: float = dc.ALPHA
    analysis_seed: int = 2026102703
    control_seed: int = 2026102704
    control_eps_step: float = 0.003      # eps step in the second half of one control
    control_contrast_step: float = -0.1  # contrast step in the second half of another
    min_jobs: int = 10


def amplitude_circuits(n: int) -> List[Circuit]:
    measure = [Op("measure", (q,), cbit=q) for q in range(n)]
    half = [Op("sx", (q,)) for q in range(n)]
    flips = [Op("x", (q,)) for q in range(n)]
    circuits = [Circuit(n, n, list(measure), name="prep0"), Circuit(n, n, flips + measure, name="prep1")]
    circuits += [Circuit(n, n, half + flips * reps + measure, name=f"amp{reps}") for reps in REPS]
    return circuits + [Circuit(n, n, flips * 100 + measure, name="x100")]


def ones(counts: Dict[str, int], n: int) -> List[int]:
    """Per-qubit count of outcome 1; qubit q is classical bit q, read from the right."""
    return dc.error_counts(counts, n, "prep0")


# --------------------------------------------------------------------------------------------
# Model and fits
# --------------------------------------------------------------------------------------------

def p_one(c, a, gamma, delta, eps) -> np.ndarray:
    """P(1) of each FITTED circuit; a and eps may be arrays shaped (..., 1) to evaluate a grid."""
    prep0 = c - a / 2 + 0 * eps
    prep1 = c + a / 2 * np.cos(eps)
    amp = c + a / 2 * np.exp(-gamma * _N) * np.sin(delta + _N * eps)
    shape = np.broadcast(prep0, amp).shape
    out = np.concatenate([np.broadcast_to(prep0, shape[:-1] + (1,)), np.broadcast_to(prep1, shape[:-1] + (1,)),
                          np.broadcast_to(amp, shape)], axis=-1)
    return np.clip(out, 1e-9, 1 - 1e-9)


def nll(theta: Sequence[float], y: np.ndarray, shots: int) -> float:
    """Binomial negative log-likelihood of counts y (jobs, len(FITTED)), up to a constant."""
    p = p_one(*theta)
    return float(-(y * np.log(p) + (shots - y) * np.log1p(-p)).sum())


def _grid_start(y: np.ndarray, shots: int) -> Tuple[float, ...]:
    """c and a from prep0 and prep1; (delta, eps) by least squares on the pooled amp frequencies."""
    freq = y.sum(0) / (shots * y.shape[0])
    c, a = (freq[0] + freq[1]) / 2, max(freq[1] - freq[0], 1e-3)
    deltas = np.arange(-0.6, 0.6 + 1e-12, 0.01)
    epss = np.arange(-EPS_RANGE, EPS_RANGE + 1e-12, 0.0002)
    model = c + a / 2 * np.sin(deltas[:, None, None] + epss[None, :, None] * _N)
    sse = ((model - freq[2:]) ** 2).sum(-1)
    i, j = np.unravel_index(np.argmin(sse), sse.shape)
    return float(c), float(a), 0.0, float(deltas[i]), float(epss[j])


def fit_pooled(y: np.ndarray, shots: int, start: Optional[Sequence[float]] = None) -> np.ndarray:
    """Maximum-likelihood (c, a, gamma, delta, eps) over every job; grid start unless one is given."""
    x0 = _grid_start(y, shots) if start is None else tuple(start)
    return minimize(nll, x0, args=(y, shots), method="L-BFGS-B", bounds=BOUNDS).x


def _grid_best(y: np.ndarray, shots: int, c: float, gamma: float, delta: float,
               eps_grid: np.ndarray, a_grid: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Each row of y's (eps, a) maximizing the binomial likelihood over the grid's points."""
    p = p_one(c, a_grid[:, None, None], gamma, delta, eps_grid[None, :, None]).reshape(-1, len(FITTED))
    best = (y @ np.log(p).T + (shots - y) @ np.log1p(-p).T).argmax(1)
    a_index, eps_index = np.unravel_index(best, (len(a_grid), len(eps_grid)))
    return eps_grid[eps_index], a_grid[a_index]


def per_job(y: np.ndarray, shots: int, theta: Sequence[float]) -> Tuple[np.ndarray, np.ndarray]:
    """Each job's (eps, a), jointly, by maximum likelihood; c, gamma, delta pooled.

    A coarse grid (5 eps steps, 4 contrast steps) over eps +- EPS_WINDOW and contrast 0-1.2, then
    the fine grid (EPS_STEP, CONTRAST_STEP) one coarse cell either side of each job's coarse best.
    """
    c, _, gamma, delta, eps = theta
    coarse_eps, coarse_a = _grid_best(y, shots, c, gamma, delta,
                                      eps + np.arange(-EPS_WINDOW, EPS_WINDOW + 1e-12, 5 * EPS_STEP),
                                      np.arange(0.0, 1.2 + 1e-12, 4 * CONTRAST_STEP))
    eps_k, a_k = np.empty(len(y)), np.empty(len(y))
    for k in range(len(y)):
        fine = _grid_best(y[k:k + 1], shots, c, gamma, delta, coarse_eps[k] + np.arange(-5, 6) * EPS_STEP,
                          np.clip(coarse_a[k] + np.arange(-4, 5) * CONTRAST_STEP, 0.0, 1.2))
        eps_k[k], a_k[k] = fine[0][0], fine[1][0]
    return eps_k, a_k


def _micro(values: np.ndarray) -> List[int]:
    return [int(round(v * MICRO)) for v in values]


# --------------------------------------------------------------------------------------------
# Analysis
# --------------------------------------------------------------------------------------------

def analyze(times: Sequence[Any], y_by_qubit: Dict[int, np.ndarray], shots: int, design: AmpDesign,
            seed: int, sessions: Optional[Sequence[str]] = None,
            x100_by_qubit: Optional[Dict[int, Sequence[int]]] = None) -> Dict[str, Any]:
    """Drift checks on each qubit's eps_k and a_k series, Holm-corrected over the whole run."""
    rng = np.random.default_rng(seed)
    qubits, flat = {}, {}
    for q, y in y_by_qubit.items():
        theta = fit_pooled(y, shots)
        eps_k, a_k = per_job(y, shots, theta)
        expected = p_one(*theta)
        reps: Dict[str, List[List[int]]] = {name: [] for name in PARAMETERS}
        for _ in range(design.replicas):
            y_rep = rng.binomial(shots, np.broadcast_to(expected, y.shape))
            e_rep, a_rep = per_job(y_rep, shots, fit_pooled(y_rep, shots, start=theta))
            reps["rotation_error"].append(_micro(e_rep))
            reps["contrast"].append(_micro(a_rep))
        record: Dict[str, Any] = {
            "pooled": dict(zip(("c", "a", "gamma", "delta", "eps"), map(float, theta))),
            "rotation_error": [float(v) for v in eps_k], "contrast": [float(v) for v in a_k],
            "at_window_edge": bool((np.abs(eps_k - theta[4]) >= EPS_WINDOW - EPS_STEP / 2).any()),
            "p": {},
        }
        for name, series in (("rotation_error", eps_k), ("contrast", a_k)):
            p = telemetry.drift_p_values(times, _micro(series), reps[name], sessions)
            record["p"][name] = {k: float(v) for k, v in p.items()}
            flat.update({f"q{q}/{name}/{k}": float(v) for k, v in p.items()})
        if x100_by_qubit is not None:
            rho = spearmanr(np.sin(50 * eps_k) ** 2, x100_by_qubit[q]).statistic
            record["x100_vs_sin2_50eps_spearman"] = None if np.isnan(rho) else float(rho)
        qubits[str(q)] = record
    rejected = dc.holm(flat, design.alpha)
    flagged = sorted(k for k, r in rejected.items() if r)
    return {"qubits": qubits, "tests": len(flat), "alpha": design.alpha, "replicas": design.replicas,
            "seed": seed, "flagged": flagged, **hypotheses(flagged)}


def hypotheses(flagged: Sequence[str]) -> Dict[str, bool]:
    """HA/HAt/HC/HCt from flagged keys "q<n>/<parameter>/<check>": dispersion vs time-order checks."""
    def any_flag(parameter: str, time_order: bool) -> bool:
        return any(k.split("/")[1] == parameter and (k.endswith("/dispersion") != time_order) for k in flagged)

    return {"HA_rotation_error_varies": any_flag("rotation_error", False),
            "HAt_rotation_error_time_order": any_flag("rotation_error", True),
            "HC_contrast_varies": any_flag("contrast", False),
            "HCt_contrast_time_order": any_flag("contrast", True)}


def synthetic(theta: Sequence[float], design: AmpDesign, rng: np.random.Generator,
              eps_step: float = 0.0, contrast_step: float = 0.0) -> np.ndarray:
    """Counts (jobs, len(REPS)) from the model; the second half of the jobs gets the steps."""
    c, a, gamma, delta, eps = theta
    y = []
    for k in range(design.jobs):
        late = k >= design.jobs // 2
        p = p_one(c, a + (contrast_step if late else 0.0), gamma, delta, eps + (eps_step if late else 0.0))
        y.append(rng.binomial(design.shots, p))
    return np.array(y)


# (c, a, gamma, delta, eps) per control qubit; eps spans what E2's interim x100 rates suggest
CONTROL_THETAS = ((0.5, 0.95, 0.0005, 0.02, 0.0068), (0.5, 0.95, 0.0005, 0.02, 0.0038),
                  (0.5, 0.95, 0.0005, 0.02, 0.0125), (0.5, 0.95, 0.0005, 0.02, 0.0063))


def software_control(design: AmpDesign) -> Dict[str, Any]:
    """The analysis on model data: constant (must flag nothing), an eps step and a contrast step."""
    rng = np.random.default_rng(design.control_seed)
    times = [design.spacing_s * k for k in range(design.jobs)]
    out = {}
    for label, steps in (("constant", (0.0, 0.0)), ("eps_step", (design.control_eps_step, 0.0)),
                         ("contrast_step", (0.0, design.control_contrast_step))):
        y = {q: synthetic(theta, design, rng, *steps) for q, theta in enumerate(CONTROL_THETAS)}
        out[label] = analyze(times, y, design.shots, design, design.control_seed)
    out["dispatch_allowed"] = not out["constant"]["flagged"]
    return out


# --------------------------------------------------------------------------------------------
# Hardware
# --------------------------------------------------------------------------------------------

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


def finish(out: Path, design: AmpDesign = AmpDesign()) -> Dict[str, Any]:
    """Telemetry and analysis from the ledger's result rows. Needs no network."""
    state = json.loads((out / "state.json").read_text())
    ledger = Ledger(out / "ledger.jsonl")
    if ledger.verify() is not None:
        raise RuntimeError(f"ledger does not verify: {ledger.verify()}")
    results, sessions = dc._ledger_results(ledger), dc._sessions(ledger)
    done = [results[j] for j in state["job_ids"] if results.get(j, {}).get("status") == "DONE"]
    summary: Dict[str, Any] = {"jobs_submitted": len(state["job_ids"]), "jobs_done": len(done),
                               "qpu_seconds": sum(r.get("qpu_seconds") or 0 for r in done)}
    if len(done) < design.min_jobs:
        summary["verdict"] = f"INCOMPLETE: {len(done)} jobs finished, fewer than {design.min_jobs}"
    else:
        layout = state["layout"]
        source = {"kind": "HARDWARE", "node": state["backend"],
                  "generator": f"bridge.amplitude_drift {dc.sha256_file(Path(__file__))[:12]}"}
        lines = write_telemetry([{"window": r["window"], "counts": r["counts"]} for r in done],
                                layout, design.shots, f"amplitude-hw-{state['backend']}", source)
        (out / "telemetry.jsonl").write_bytes(b"".join(line + b"\n" for line in lines))
        times = telemetry.metric_series(lines, f"q{layout[0]}_amp0_ones")[0]
        y = {q: np.array([[int(v) for v in telemetry.metric_series(lines, f"q{q}_{name}_ones")[1]]
                          for name in FITTED]).T for q in layout}
        x100 = {q: [int(v) for v in telemetry.metric_series(lines, f"q{q}_x100_errors")[1]] for q in layout}
        labels = [sessions[r["job_id"]][:12] for r in sorted(done, key=lambda r: r["window"][0])]
        summary["sessions"] = sorted(set(labels))
        summary["analysis"] = analyze(times, y, design.shots, design, design.analysis_seed, labels, x100)
    (out / "analysis.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def hardware(out: Path, key_file: Path, prereg: Path, after: Path, design: AmpDesign = AmpDesign()) -> Dict[str, Any]:
    """Software control, then design.jobs ledger-first jobs on E2's backend and qubits, then results."""
    from dnalang.backends import ibm
    previous = json.loads(after.read_text())
    out.mkdir(parents=True, exist_ok=True)
    ledger = Ledger(out / "ledger.jsonl")
    hashes = {**dc.code_hashes(), "bridge/amplitude_drift.py": dc.sha256_file(Path(__file__))}
    ledger.append("preregistration", {"sha256": dc.sha256_file(prereg), "code_sha256": hashes})
    control = software_control(design)
    (out / "control.json").write_text(json.dumps(control, indent=1) + "\n")
    if not control["dispatch_allowed"]:
        ledger.append("aborted", {"reason": "software control flagged drift in constant-eps data"})
        return {"verdict": "NOT RUN: the constant-eps control flagged drift", "control": control["constant"]["flagged"]}

    dc._use_key_file(key_file)
    svc = ibm.service()
    backend = svc.backend(previous["backend"])
    if not backend.status().operational:
        ledger.append("aborted", {"reason": f"{backend.name} not operational"})
        return {"verdict": f"NOT RUN: {backend.name} not operational"}
    layout = previous["layout"]
    circuits = amplitude_circuits(len(layout))
    ledger.append("design", {**asdict(design), "backend": backend.name, "layout": layout, "reps": list(REPS),
                             "circuits": {c.name: c.sha256() for c in circuits}})
    state = {"backend": backend.name, "layout": layout, "job_ids": []}
    for j in range(design.jobs):
        began = time.monotonic()
        job_id = ibm.submit(circuits, backend.name, layout, design.shots, ledger,
                            tags={"experiment": "amplitude_drift", "job": j})
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
    parser = argparse.ArgumentParser(prog="python -m bridge.amplitude_drift")
    sub = parser.add_subparsers(dest="command", required=True)
    hw = sub.add_parser("hardware")
    hw.add_argument("--key-file", type=Path, required=True)
    hw.add_argument("--out", type=Path, required=True)
    hw.add_argument("--prereg", type=Path, required=True)
    hw.add_argument("--after", type=Path, required=True, help="E2's state.json: backend and layout")
    fin = sub.add_parser("finish")
    fin.add_argument("--out", type=Path, required=True)
    fin.add_argument("--key-file", type=Path)
    args = parser.parse_args(argv)
    if args.command == "hardware":
        print(json.dumps(hardware(args.out, args.key_file, args.prereg, args.after), indent=1, default=str))
    else:
        if args.key_file:
            from dnalang.backends import ibm
            dc._use_key_file(args.key_file)
            state = json.loads((args.out / "state.json").read_text())
            dc.collect(ibm.service(), Ledger(args.out / "ledger.jsonl"), state["job_ids"])
        print(json.dumps(finish(args.out), indent=1, default=str))


if __name__ == "__main__":
    main()
