"""
bridge.drift_check — do posterior predictive drift checks catch calibration drift?
=================================================================================

Pre-registered in ``results/drift_checks/PREREGISTRATION.md``. Both experiments use the checks
in ``osiris_governance.telemetry`` (dispersion, lag-1 autocorrelation, trend r², CUSUM peak,
session eta²): each is an upper-tail predictive p-value of one series against replicas of it
drawn from a model without drift.

No-drift model M0 for K error counts y_k out of N shots each: one rate for the whole series,
p ~ Beta(1 + Σy, 1 + Σ(N − y)), then y_k ~ Binomial(N, p) for every job, in the observed order.

E1 (simulation) measures how often each check alarms (p ≤ α) when the per-job rate is
constant, drifts linearly, steps between two sessions, wanders as a random walk, or jitters
independently from job to job. Two scale-dependent comparators run beside them: the squared
least-squares slope as a predictive check, and tabular CUSUMs on model-standardised residuals
with fixed limits h = 27 and h = 5 (allowance k = 0.5).

E2 (hardware) runs K jobs on one IBM backend, spaced in time. Each job measures readout and
gate-error counts on a fixed set of qubits. The jobs are written as an OSIRIS-TELEMETRY-V1 run
and every series is checked, with Holm's correction over all p-values. An Aer run of the same
jobs under the backend's static noise model is the control.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from osiris_governance import telemetry
from osiris_governance.telemetry import TelemetryRun, quantity

from dnalang.ledger import Ledger
from dnalang.qc_ir import Circuit, Op

ALPHA = 0.05
CHECKS = telemetry.DRIFT_CHECKS
TIME_ORDER_CHECKS = tuple(c for c in CHECKS if c != "dispersion")
PPC_CHECKS = CHECKS + ("slope_sq",)
FIXED_LIMITS = {"cusum_h27": 27.0, "cusum_h5": 5.0}
CUSUM_ALLOWANCE = 0.5

SCENARIOS: Tuple[Tuple[str, float], ...] = (
    ("null", 0), ("linear", 2), ("linear", 4), ("linear", 8), ("step", 2), ("step", 4), ("step", 8),
    ("walk", 0.5), ("walk", 1), ("walk", 2), ("jitter", 1), ("jitter", 2), ("jitter", 4))

CIRCUITS = ("prep0", "prep1", "x100")
ERROR_BIT = {"prep0": "1", "prep1": "0", "x100": "1"}       # the outcome that counts as an error
HERE = Path(__file__).resolve()


@dataclass(frozen=True)
class SimDesign:
    jobs: int = 20
    shots: int = 4000
    rates: Tuple[float, ...] = (0.005, 0.02)
    trials: int = 1000
    replicas: int = 400
    alpha: float = ALPHA
    seed: int = 20261027


@dataclass(frozen=True)
class HardwareDesign:
    jobs: int = 20
    spacing_s: int = 300
    shots: int = 4000
    qubits: int = 4
    replicas: int = 4000
    alpha: float = ALPHA
    analysis_seed: int = 2026102701
    control_seed: int = 2026102702
    min_jobs: int = 10


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def code_hashes() -> Dict[str, str]:
    return {"bridge/drift_check.py": sha256_file(HERE),
            "osiris_governance/telemetry.py": sha256_file(Path(telemetry.__file__))}


# --------------------------------------------------------------------------------------------
# Replicas and fast statistics
# --------------------------------------------------------------------------------------------

def binomial_sd(rate: float, shots: int) -> float:
    return float(np.sqrt(rate * (1 - rate) / shots))


def posterior_predictive(y: np.ndarray, shots: int, n: int, rng: np.random.Generator) -> np.ndarray:
    """Replicas under M0 for each series in y (..., K): integer counts shaped (..., n, K)."""
    total = y.sum(-1)
    p = rng.beta(1 + total, 1 + shots * y.shape[-1] - total, size=(n,) + total.shape)
    p = np.moveaxis(p, 0, -1)
    return rng.binomial(shots, p[..., None], size=p.shape + (y.shape[-1],))


def np_checks(x: np.ndarray, t: np.ndarray, sessions: np.ndarray) -> Dict[str, np.ndarray]:
    """The governance checks plus the squared slope, over the last axis of x, in float64.

    Matches osiris_governance.telemetry to rounding (tests/test_drift_check.py); E1 needs the speed.
    """
    k = x.shape[-1]
    d = x - x.mean(-1, keepdims=True)
    ss = (d * d).sum(-1)
    flat = ss <= 0
    safe = np.where(flat, 1.0, ss)
    dt = t - t.mean()
    stt = float(dt @ dt)
    sxy = d @ dt
    run = np.cumsum(d, -1)
    between = sum(m.sum() * x[..., m].mean(-1) ** 2 for m in (sessions == g for g in np.unique(sessions)))
    between = between - k * x.mean(-1) ** 2
    return {
        "dispersion": ss / (k - 1),
        "lag1_autocorrelation": np.where(flat, 0.0, (d[..., :-1] * d[..., 1:]).sum(-1) / safe),
        "trend_r2": np.where(flat, 0.0, sxy * sxy / (stt * safe)),
        "cusum_peak": np.where(flat, 0.0, (run * run).max(-1) / safe),
        "session_eta2": np.where(flat, 0.0, between / safe),
        "slope_sq": (sxy / stt) ** 2,
    }


def np_p_values(stats: np.ndarray) -> np.ndarray:
    """Upper-tail predictive p-values from (..., 1 + S) statistics, observed first; ties count as extreme."""
    return (1 + (stats[..., 1:] >= stats[..., :1]).sum(-1)) / stats.shape[-1]


def fixed_cusum_alarm(y: np.ndarray, shots: int, limit: float, allowance: float = CUSUM_ALLOWANCE) -> np.ndarray:
    """Tabular CUSUM on residuals standardised by the pooled binomial model, with a fixed limit."""
    p = y.sum(-1, keepdims=True) / (shots * y.shape[-1])
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.nan_to_num((y / shots - p) / np.sqrt(p * (1 - p) / shots))
    upper = lower = np.zeros(y.shape[:-1])
    hit = np.zeros(y.shape[:-1], dtype=bool)
    for j in range(y.shape[-1]):
        upper = np.maximum(0.0, upper + z[..., j] - allowance)
        lower = np.maximum(0.0, lower - z[..., j] - allowance)
        hit |= (upper > limit) | (lower > limit)
    return hit


def holm(pvalues: Dict[str, float], alpha: float) -> Dict[str, bool]:
    """Holm's step-down rejections at family-wise error rate alpha."""
    order = sorted(pvalues, key=lambda key: (pvalues[key], key))
    rejected, going = {}, True
    for i, key in enumerate(order):
        going = going and pvalues[key] <= alpha / (len(order) - i)
        rejected[key] = going
    return rejected


# --------------------------------------------------------------------------------------------
# E1: operating characteristics
# --------------------------------------------------------------------------------------------

def true_rates(kind: str, size: float, rate: float, design: SimDesign, rng: np.random.Generator,
               trials: int) -> np.ndarray:
    """Per-job error rates (trials, K) around rate; size is in binomial SDs of one job's rate."""
    k = np.arange(design.jobs)
    sigma = binomial_sd(rate, design.shots)
    if kind == "null":
        shift = np.zeros((trials, design.jobs))
    elif kind == "linear":
        shift = np.broadcast_to(size * sigma * (k / (design.jobs - 1) - 0.5), (trials, design.jobs))
    elif kind == "step":
        shift = np.broadcast_to(size * sigma * np.where(k < design.jobs // 2, -0.5, 0.5), (trials, design.jobs))
    elif kind == "walk":
        walk = np.cumsum(rng.normal(0, size * sigma, (trials, design.jobs)), axis=1)
        shift = walk - walk.mean(axis=1, keepdims=True)
    elif kind == "jitter":
        shift = rng.normal(0, size * sigma, (trials, design.jobs))
    else:
        raise ValueError(f"unknown scenario {kind!r}")
    return np.clip(rate + shift, 1e-9, 1 - 1e-9)


def operating_characteristics(design: SimDesign, scenarios: Sequence[Tuple[str, float]] = SCENARIOS,
                              chunk: int = 50) -> List[Dict[str, Any]]:
    t = np.arange(design.jobs, dtype=float)
    sessions = np.where(np.arange(design.jobs) < design.jobs // 2, 0, 1)
    rows = []
    for (r, rate), (index, (kind, size)) in itertools.product(enumerate(design.rates), enumerate(scenarios)):
        rng = np.random.default_rng([design.seed, r, index])
        alarms = dict.fromkeys(PPC_CHECKS + tuple(FIXED_LIMITS), 0)
        done = 0
        while done < design.trials:
            n = min(chunk, design.trials - done)
            y = rng.binomial(design.shots, true_rates(kind, size, rate, design, rng, n))
            reps = posterior_predictive(y, design.shots, design.replicas, rng)
            stats = np_checks(np.concatenate([y[:, None, :], reps], axis=1).astype(float), t, sessions)
            for name in PPC_CHECKS:
                alarms[name] += int((np_p_values(stats[name]) <= design.alpha).sum())
            for name, limit in FIXED_LIMITS.items():
                alarms[name] += int(fixed_cusum_alarm(y, design.shots, limit).sum())
            done += n
        rows.append({"rate": rate, "scenario": kind, "size_sd": size,
                     "alarm_rate": {name: count / design.trials for name, count in alarms.items()}})
    return rows


def e1_verdict(rows: List[Dict[str, Any]], design: SimDesign) -> Dict[str, Any]:
    """The pre-registered predictions P1-P5, judged from the operating-characteristics table.

    Each prediction must hold at every base rate in design.rates.
    """
    table = {(r["rate"], r["scenario"], r["size_sd"]): r["alarm_rate"] for r in rows}
    bound = design.alpha + 2 * float(np.sqrt(design.alpha * (1 - design.alpha) / design.trials))

    def judge(items, test):
        detail = {f"{rate:g}/{s}/{size:g}/{check}": table[(rate, s, size)][check]
                  for rate in design.rates for s, size, check in items}
        return {"pass": all(test(v) for v in detail.values()), "rates": detail}

    return {
        "false_alarm_bound": bound,
        "P1_calibrated_under_null": judge([("null", 0, c) for c in PPC_CHECKS], lambda v: v <= bound),
        "P2_sensitive_to_drift": judge([("linear", 8, "trend_r2"), ("linear", 8, "cusum_peak"),
                                        ("step", 8, "cusum_peak"), ("step", 8, "session_eta2"),
                                        ("walk", 2, "lag1_autocorrelation")], lambda v: v >= 0.80),
        "P3_time_order_checks_ignore_jitter": judge([("jitter", 4, c) for c in TIME_ORDER_CHECKS],
                                                    lambda v: v <= bound),
        "P3_dispersion_detects_jitter": judge([("jitter", 4, "dispersion")], lambda v: v >= 0.80),
        "P4_slope_check_confuses_jitter_with_drift": judge([("jitter", 4, "slope_sq")], lambda v: v > bound),
        "P5_h27_misses_linear_drift": judge([("linear", 8, "cusum_h27")], lambda v: v < 0.50),
    }


def simulate(out: Path, prereg: Path, design: SimDesign = SimDesign()) -> Dict[str, Any]:
    out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    frozen = {"preregistration_sha256": sha256_file(prereg), "code_sha256": code_hashes()}
    rows = operating_characteristics(design)
    record = {**frozen, "design": asdict(design), "scenarios": rows,
              "numpy": np.__version__, "seconds": round(time.time() - started, 1)}
    (out / "oc.json").write_text(json.dumps(record, indent=1) + "\n")
    verdict = e1_verdict(rows, design)
    (out / "verdict.json").write_text(json.dumps(verdict, indent=1) + "\n")
    return verdict


# --------------------------------------------------------------------------------------------
# E2: hardware
# --------------------------------------------------------------------------------------------

def probe_circuits(n: int) -> List[Circuit]:
    """prep0 (readout of |0>), prep1 (X then readout) and x100 (100 X gates, ideally |0>) on n qubits."""
    measure = [Op("measure", (q,), cbit=q) for q in range(n)]
    flips = [Op("x", (q,)) for q in range(n)]
    return [Circuit(n, n, list(measure), name="prep0"),
            Circuit(n, n, flips + measure, name="prep1"),
            Circuit(n, n, flips * 100 + measure, name="x100")]


def error_counts(counts: Dict[str, int], n: int, circuit: str) -> List[int]:
    """Per-qubit error counts; qubit q is classical bit q, the (q+1)-th character from the right."""
    bad = ERROR_BIT[circuit]
    return [sum(c for key, c in counts.items() if key.replace(" ", "")[-1 - q] == bad) for q in range(n)]


def _stamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def _now() -> str:
    return _stamp(datetime.now(timezone.utc))


def metric_name(physical: int, circuit: str) -> str:
    return f"q{physical}_{circuit}_errors"


def write_telemetry(jobs: List[Dict[str, Any]], layout: Sequence[int], shots: int, run_id: str,
                    source: Dict[str, Any], evidence_class: str) -> List[bytes]:
    """One frame per job, in execution order. A job is {"window": (start, end), "counts": [...]}."""
    run = TelemetryRun(run_id, source, evidence_class)
    lines = []
    for job in sorted(jobs, key=lambda j: j["window"][0]):
        metrics = {"shots": quantity(shots, 0, "1")}
        for circuit, counts in zip(CIRCUITS, job["counts"]):
            for physical, errors in zip(layout, error_counts(counts, len(layout), circuit)):
                metrics[metric_name(physical, circuit)] = quantity(errors, 0, "1")
        start, end = job["window"]
        lines.append(run.frame(_now(), {"start": start, "end": end, "samples": shots}, metrics))
    return lines


def analyze_run(lines: Sequence[bytes], layout: Sequence[int], shots: int, replicas: int, seed: int,
                alpha: float, sessions: Optional[Sequence[str]] = None) -> Dict[str, Any]:
    """Every check on every series of a verified run, Holm-corrected over all p-values together."""
    rng = np.random.default_rng(seed)
    series, flat = {}, {}
    for physical in layout:
        for circuit in CIRCUITS:
            name = metric_name(physical, circuit)
            times, values = telemetry.metric_series(lines, name)
            counts = [int(v) for v in values]
            reps = posterior_predictive(np.array(counts)[None], shots, replicas, rng)[0].tolist()
            p = telemetry.drift_p_values(times, counts, reps, sessions)
            rate = sum(counts) / (shots * len(counts))
            binomial_var = shots * rate * (1 - rate)
            series[name] = {
                "jobs": len(counts), "counts": counts, "rate": rate,
                "overdispersion": float(telemetry.dispersion(counts)) / binomial_var if binomial_var else None,
                "slope_per_hour": float(telemetry.trend_slope(times, counts)) * 3600 / shots,
                "p": {k: float(v) for k, v in p.items()}, "p_exact": {k: str(v) for k, v in p.items()},
            }
            flat.update({f"{name}/{k}": float(v) for k, v in p.items()})
    rejected = holm(flat, alpha)
    flagged = sorted(k for k, r in rejected.items() if r)
    return {"series": series, "tests": len(flat), "alpha": alpha, "replicas": replicas, "seed": seed,
            "flagged": flagged,
            "heterogeneity": any(k.endswith("/dispersion") for k in flagged),
            "time_structure": any(not k.endswith("/dispersion") for k in flagged)}


def _use_key_file(path: Path) -> None:
    """Loads the API key into this process's environment for dnalang.backends.ibm; never printed."""
    with open(path) as f:
        os.environ["IBM_QUANTUM_TOKEN"] = json.load(f)["apikey"]


def choose_layout(backend, n: int) -> Tuple[List[int], Dict[int, float]]:
    """The n qubits with the lowest reported readout error, in index order."""
    errors = {qs[0]: props.error for qs, props in backend.target["measure"].items()
              if props is not None and props.error is not None}
    best = sorted(errors, key=lambda q: (errors[q], q))[:n]
    return sorted(best), {q: errors[q] for q in sorted(best)}


def preflight(key_file: Path, design: HardwareDesign = HardwareDesign()) -> Dict[str, Any]:
    """Backend, layout and account usage, without submitting anything."""
    from dnalang.backends import ibm
    _use_key_file(key_file)
    svc = ibm.service()
    backend = svc.least_busy(operational=True, simulator=False)
    layout, readout = choose_layout(backend, design.qubits)
    usage = svc.usage()
    status = backend.status()
    return {"backend": backend.name, "num_qubits": backend.num_qubits, "pending_jobs": status.pending_jobs,
            "layout": layout, "readout_error": readout, "calibration_hash": ibm.calibration_hash(backend),
            "usage": {k: v for k, v in usage.items() if "usage" in k or "period" in k}}


def control_run(backend, layout: Sequence[int], design: HardwareDesign) -> List[Dict[str, Any]]:
    """The same jobs on Aer under the backend's noise model at dispatch time: no drift by construction."""
    from qiskit_aer import AerSimulator

    from dnalang.backends import ibm
    sim = AerSimulator.from_backend(backend)
    tqcs = ibm.transpile_fixed(probe_circuits(len(layout)), backend, layout)
    jobs = []
    for j in range(design.jobs):
        start = _now()
        result = sim.run(tqcs, shots=design.shots, seed_simulator=design.control_seed + j).result()
        jobs.append({"window": (start, _now()), "counts": [result.get_counts(i) for i in range(len(tqcs))]})
    return jobs


def _ledger_results(ledger: Ledger) -> Dict[str, Dict[str, Any]]:
    return {row["job_id"]: row for row in ledger if row.get("kind") == "result"}


def _sessions(ledger: Ledger) -> Dict[str, str]:
    """job_id -> calibration hash recorded in its write-ahead submit-intent row."""
    intents = {row["hash"]: row["calibration_hash"] for row in ledger if row.get("kind") == "submit-intent"}
    return {row["job_id"]: intents[row["intent_hash"]] for row in ledger if row.get("kind") == "submitted"}


def collect(svc, ledger: Ledger, job_ids: Sequence[str]) -> None:
    """Appends one result row per finished job (counts, execution window, QPU seconds), once."""
    done = _ledger_results(ledger)
    for job_id in job_ids:
        if job_id in done:
            continue
        job = svc.job(job_id)
        status = str(job.status())
        if status not in ("DONE", "ERROR", "CANCELLED"):
            continue
        row: Dict[str, Any] = {"job_id": job_id, "status": status}
        if status == "DONE":
            result = job.result()
            spans = result.metadata["execution"]["execution_spans"]
            row["counts"] = [pub.data.c.get_counts() for pub in result]
            row["counts_hash"] = hashlib.sha256(json.dumps(row["counts"], sort_keys=True).encode()).hexdigest()
            row["window"] = [_stamp(spans.start), _stamp(spans.stop)]
            try:
                row["qpu_seconds"] = job.usage()
            except Exception:
                row["qpu_seconds"] = None
        ledger.append("result", row)


def finish(out: Path, design: HardwareDesign = HardwareDesign()) -> Dict[str, Any]:
    """Telemetry and analysis from the ledger's result rows. Needs no network."""
    state = json.loads((out / "state.json").read_text())
    ledger = Ledger(out / "ledger.jsonl")
    if ledger.verify() is not None:
        raise RuntimeError(f"ledger does not verify: {ledger.verify()}")
    results, sessions = _ledger_results(ledger), _sessions(ledger)
    done = [results[j] for j in state["job_ids"] if results.get(j, {}).get("status") == "DONE"]
    summary: Dict[str, Any] = {"jobs_submitted": len(state["job_ids"]), "jobs_done": len(done),
                               "qpu_seconds": sum(r.get("qpu_seconds") or 0 for r in done)}
    if len(done) < design.min_jobs:
        summary["verdict"] = f"INCOMPLETE: {len(done)} jobs finished, fewer than {design.min_jobs}"
    else:
        source = {"kind": "HARDWARE", "node": state["backend"],
                  "generator": f"bridge.drift_check {code_hashes()['bridge/drift_check.py'][:12]}"}
        lines = write_telemetry([{"window": r["window"], "counts": r["counts"]} for r in done],
                                state["layout"], design.shots, f"drift-hw-{state['backend']}", source, "MEASURED")
        (out / "telemetry.jsonl").write_bytes(b"".join(line + b"\n" for line in lines))
        labels = [sessions[r["job_id"]][:12] for r in sorted(done, key=lambda r: r["window"][0])]
        summary["sessions"] = sorted(set(labels))
        summary["analysis"] = analyze_run(lines, state["layout"], design.shots, design.replicas,
                                          design.analysis_seed, design.alpha, labels)
    (out / "analysis.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def hardware(out: Path, key_file: Path, prereg: Path, design: HardwareDesign = HardwareDesign()) -> Dict[str, Any]:
    """Control, then design.jobs ledger-first submissions spaced design.spacing_s apart, then results."""
    from dnalang.backends import ibm
    out.mkdir(parents=True, exist_ok=True)
    ledger = Ledger(out / "ledger.jsonl")
    _use_key_file(key_file)
    svc = ibm.service()
    backend = svc.least_busy(operational=True, simulator=False)
    layout, readout = choose_layout(backend, design.qubits)
    circuits = probe_circuits(len(layout))
    ledger.append("preregistration", {"sha256": sha256_file(prereg), "code_sha256": code_hashes()})
    ledger.append("design", {**asdict(design), "backend": backend.name, "layout": layout,
                             "readout_error": {str(q): e for q, e in readout.items()},
                             "circuits": {c.name: c.sha256() for c in circuits}})

    control_dir = out / "control"
    control_dir.mkdir(exist_ok=True)
    source = {"kind": "SIMULATOR", "node": "aer", "seed": design.control_seed,
              "generator": f"qiskit-aer AerSimulator.from_backend({backend.name})"}
    lines = write_telemetry(control_run(backend, layout, design), layout, design.shots,
                            f"drift-control-{backend.name}", source, "SIMULATED")
    (control_dir / "telemetry.jsonl").write_bytes(b"".join(line + b"\n" for line in lines))
    control = analyze_run(lines, layout, design.shots, design.replicas, design.analysis_seed, design.alpha)
    (control_dir / "analysis.json").write_text(json.dumps(control, indent=1) + "\n")

    state = {"backend": backend.name, "layout": layout, "job_ids": []}
    for j in range(design.jobs):
        began = time.monotonic()
        job_id = ibm.submit(circuits, backend.name, layout, design.shots, ledger,
                            tags={"experiment": "drift_checks", "job": j})
        state["job_ids"].append(job_id)
        (out / "state.json").write_text(json.dumps(state, indent=1) + "\n")
        print(f"submitted {j + 1}/{design.jobs}: {job_id}", flush=True)
        if j + 1 < design.jobs:
            time.sleep(max(0.0, design.spacing_s - (time.monotonic() - began)))

    deadline = time.monotonic() + 8 * 3600
    while time.monotonic() < deadline:
        collect(svc, ledger, state["job_ids"])
        if len(_ledger_results(ledger)) >= len(state["job_ids"]):
            break
        time.sleep(60)
    return finish(out, design)


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(prog="python -m bridge.drift_check")
    sub = parser.add_subparsers(dest="command", required=True)
    sim = sub.add_parser("simulate")
    sim.add_argument("--out", type=Path, required=True)
    sim.add_argument("--prereg", type=Path, required=True)
    pre = sub.add_parser("preflight")
    pre.add_argument("--key-file", type=Path, required=True)
    hw = sub.add_parser("hardware")
    hw.add_argument("--key-file", type=Path, required=True)
    hw.add_argument("--out", type=Path, required=True)
    hw.add_argument("--prereg", type=Path, required=True)
    fin = sub.add_parser("finish")
    fin.add_argument("--out", type=Path, required=True)
    fin.add_argument("--key-file", type=Path)
    args = parser.parse_args(argv)
    if args.command == "simulate":
        print(json.dumps(simulate(args.out, args.prereg), indent=1))
    elif args.command == "preflight":
        print(json.dumps(preflight(args.key_file), indent=1, default=str))
    elif args.command == "hardware":
        print(json.dumps(hardware(args.out, args.key_file, args.prereg), indent=1, default=str))
    else:
        if args.key_file:
            from dnalang.backends import ibm
            _use_key_file(args.key_file)
            state = json.loads((args.out / "state.json").read_text())
            collect(ibm.service(), Ledger(args.out / "ledger.jsonl"), state["job_ids"])
        print(json.dumps(finish(args.out), indent=1, default=str))


if __name__ == "__main__":
    main()
