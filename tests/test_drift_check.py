"""bridge.drift_check: fast statistics against the exact ones, circuits, ledger-to-analysis path.

Synthetic counts and fake backends only; nothing here touches IBM hardware.
"""
import json
from fractions import Fraction

import numpy as np
import pytest
from osiris_governance import telemetry

from bridge import drift_check as dc
from dnalang.backends import aer
from dnalang.ledger import Ledger


def test_np_checks_match_the_exact_governance_checks():
    rng = np.random.default_rng(5)
    x = rng.binomial(4000, 0.02, (30, 12))
    t = np.sort(rng.choice(10_000, 12, replace=False)).astype(float)
    sessions = np.array([0] * 5 + [1] * 4 + [2] * 3)
    fast = dc.np_checks(x.astype(float), t, sessions)
    times = [int(v) for v in t]
    for i, row in enumerate(x.tolist()):
        exact = {"dispersion": telemetry.dispersion(row),
                 "lag1_autocorrelation": telemetry.lag1_autocorrelation(row),
                 "trend_r2": telemetry.trend_r2(times, row),
                 "cusum_peak": telemetry.cusum_peak(row),
                 "session_eta2": telemetry.session_eta2(row, sessions.tolist()),
                 "slope_sq": telemetry.trend_slope(times, row) ** 2}
        for name, value in exact.items():
            assert fast[name][i] == pytest.approx(float(value), rel=1e-9, abs=1e-15), name


def test_np_checks_give_zero_for_a_constant_series():
    stats = dc.np_checks(np.full((1, 6), 80.0), np.arange(6.0), np.array([0, 0, 0, 1, 1, 1]))
    assert all(float(v[0]) == 0 for k, v in stats.items() if k != "slope_sq")


def test_np_p_values_match_the_exact_ones_with_ties():
    stats = np.random.default_rng(8).integers(0, 4, (200, 9))             # small integers: many ties
    fast = dc.np_p_values(stats.astype(float))
    for row, p in zip(stats.tolist(), fast):
        assert p == float(telemetry.predictive_p_value(row[0], row[1:]))
    assert dc.np_p_values(np.array([5.0, 5, 1, 6])) == 0.75


def test_posterior_predictive_is_seeded_integer_counts():
    y = np.array([[80, 75, 90], [10, 12, 8]])
    a = dc.posterior_predictive(y, 4000, 7, np.random.default_rng(1))
    b = dc.posterior_predictive(y, 4000, 7, np.random.default_rng(1))
    assert a.shape == (2, 7, 3) and a.dtype.kind == "i"
    assert (a == b).all()
    assert abs(a[0].mean() - 81.7) < 10 and abs(a[1].mean() - 10) < 5      # replicas sit near each series' rate


def test_error_counts_read_bits_from_the_right():
    counts = {"0001": 10, "1000": 5, "0000": 85}
    assert dc.error_counts(counts, 4, "prep0") == [10, 0, 0, 5]
    assert dc.error_counts(counts, 4, "prep1") == [90, 100, 100, 95]


def test_noiseless_baseline_has_no_errors():
    circuits = dc.probe_circuits(2)
    assert [c.name for c in circuits] == list(dc.CIRCUITS)
    assert circuits[2].count("x") == 200 and circuits[1].count("x") == 2
    for circuit, counts in zip(dc.CIRCUITS, aer.run(circuits, shots=200, seed=3)):
        assert dc.error_counts(counts, 2, circuit) == [0, 0]


def test_holm_steps_down_and_stops():
    p = {"a": 0.001, "b": 0.012, "c": 0.02, "d": 0.5}
    assert dc.holm(p, 0.05) == {"a": True, "b": True, "c": True, "d": False}
    # c fails 0.05 / 2, so b is not rejected although 0.04 <= 0.05 / 1
    assert dc.holm({"a": 0.001, "b": 0.04, "c": 0.03}, 0.05) == {"a": True, "c": False, "b": False}


def test_fixed_cusum_alarms_on_a_large_step_only():
    flat = np.full((1, 20), 80)
    step = np.array([[60] * 10 + [100] * 10])          # about +-2.3 binomial SDs, 10 jobs each side
    assert not dc.fixed_cusum_alarm(flat, 4000, 5.0)[0]
    assert dc.fixed_cusum_alarm(step, 4000, 5.0)[0]
    assert not dc.fixed_cusum_alarm(step, 4000, 27.0)[0]


def test_true_rates_shapes_and_sizes():
    design = dc.SimDesign(jobs=10)
    sigma = dc.binomial_sd(0.02, design.shots)
    rng = np.random.default_rng(0)
    linear = dc.true_rates("linear", 4, 0.02, design, rng, 3)
    assert linear.shape == (3, 10)
    assert linear[0, -1] - linear[0, 0] == pytest.approx(4 * sigma)
    step = dc.true_rates("step", 2, 0.02, design, rng, 1)[0]
    assert step[5] - step[4] == pytest.approx(2 * sigma) and len(set(step[:5])) == 1
    assert dc.true_rates("walk", 1, 0.02, design, rng, 4).mean(axis=1) == pytest.approx(0.02)
    with pytest.raises(ValueError, match="unknown scenario"):
        dc.true_rates("sine", 1, 0.02, design, rng, 1)


def test_operating_characteristics_smoke():
    design = dc.SimDesign(jobs=10, trials=30, replicas=39, rates=(0.02,))
    rows = dc.operating_characteristics(design, [("null", 0), ("step", 8)], chunk=7)
    rates = {r["scenario"]: r["alarm_rate"] for r in rows}
    assert set(rates["null"]) == set(dc.PPC_CHECKS) | set(dc.FIXED_LIMITS)
    assert rates["step"]["session_eta2"] == 1 and rates["step"]["cusum_peak"] == 1
    assert rates["null"]["session_eta2"] < 0.3


def test_e1_verdict_applies_the_stated_rules():
    design = dc.SimDesign()
    rows = [{"rate": rate, "scenario": s, "size_sd": size,
             "alarm_rate": dict.fromkeys(dc.PPC_CHECKS + tuple(dc.FIXED_LIMITS), 0.04)}
            for rate in design.rates for s, size in dc.SCENARIOS]
    table = {(r["rate"], r["scenario"], r["size_sd"]): r["alarm_rate"] for r in rows}
    for rate in design.rates:
        for check in ("trend_r2", "cusum_peak"):
            table[(rate, "linear", 8)][check] = 0.95
        table[(rate, "step", 8)].update(cusum_peak=0.99, session_eta2=0.99)
        table[(rate, "walk", 2)]["lag1_autocorrelation"] = 0.85
        table[(rate, "jitter", 4)].update(dispersion=1.0, slope_sq=0.6)
    verdict = dc.e1_verdict(rows, design)
    assert all(v["pass"] for k, v in verdict.items() if k.startswith("P"))
    assert len(verdict["P2_sensitive_to_drift"]["rates"]) == 5 * len(design.rates)
    table[(0.005, "jitter", 4)]["trend_r2"] = 0.07          # above alpha + 2 SE = 0.0638, at one rate only
    assert not dc.e1_verdict(rows, design)["P3_time_order_checks_ignore_jitter"]["pass"]
    table[(0.02, "linear", 8)]["cusum_h27"] = 0.5
    assert not dc.e1_verdict(rows, design)["P5_h27_misses_linear_drift"]["pass"]


def fake_jobs(n, layout, rng, shots=4000):
    """Job records as collect() stores them: counts per circuit and an execution window."""
    jobs = []
    for j in range(n):
        counts = []
        for circuit in dc.CIRCUITS:
            errors = rng.binomial(shots, 0.02, len(layout))
            bits = {}
            for q, e in enumerate(errors):          # one error pattern per qubit keeps marginals exact
                key = ["1" if circuit == "prep1" else "0"] * len(layout)
                key[-1 - q] = dc.ERROR_BIT[circuit]
                bits["".join(key)] = int(e)
            ideal = "".join(["1" if circuit == "prep1" else "0"] * len(layout))
            bits[ideal] = shots - int(errors.sum())
            counts.append(bits)
        jobs.append({"window": (f"2026-10-26T15:{j:02d}:00.000000Z", f"2026-10-26T15:{j:02d}:03.500000Z"),
                     "counts": counts, "errors": None})
    return jobs


def test_telemetry_round_trip_and_analysis():
    layout = [3, 7]
    jobs = fake_jobs(12, layout, np.random.default_rng(2))
    source = {"kind": "HARDWARE", "node": "fake", "generator": "test"}
    lines = dc.write_telemetry(jobs, layout, 4000, "drift-test", source, "MEASURED")
    assert telemetry.verify_run(lines).ok
    times, values = telemetry.metric_series(lines, "q7_prep1_errors")
    assert times[1] == 60 and len(values) == 12
    plain = dc.analyze_run(lines, layout, 4000, 99, 1, 0.05)
    assert plain["tests"] == len(layout) * len(dc.CIRCUITS) * 4
    sessions = ["a"] * 6 + ["b"] * 6
    split = dc.analyze_run(lines, layout, 4000, 99, 1, 0.05, sessions)
    assert split["tests"] == len(layout) * len(dc.CIRCUITS) * 5
    series = split["series"]["q3_x100_errors"]
    assert Fraction(series["p_exact"]["dispersion"]) == Fraction(series["p"]["dispersion"]).limit_denominator(100)


def write_ledger(out, jobs, done):
    """The rows dnalang.backends.ibm.submit and collect() write, for len(jobs) jobs."""
    ledger = Ledger(out / "ledger.jsonl")
    ids = []
    for j, job in enumerate(jobs):
        intent = ledger.append("submit-intent", {"calibration_hash": "cal-a" if j < 6 else "cal-b"})
        ledger.append("submitted", {"intent_hash": intent["hash"], "job_id": f"job{j}"})
        ids.append(f"job{j}")
    for j, job in enumerate(jobs[:done]):
        ledger.append("result", {"job_id": f"job{j}", "status": "DONE", "counts": job["counts"],
                                 "window": list(job["window"]), "qpu_seconds": 2})
    (out / "state.json").write_text(json.dumps({"backend": "fake", "layout": [3, 7], "job_ids": ids}))


def test_finish_builds_telemetry_and_sessions_from_the_ledger(tmp_path):
    write_ledger(tmp_path, fake_jobs(12, [3, 7], np.random.default_rng(4)), done=12)
    summary = dc.finish(tmp_path, dc.HardwareDesign(replicas=49))
    assert summary["jobs_done"] == 12 and summary["qpu_seconds"] == 24
    assert summary["sessions"] == ["cal-a", "cal-b"]
    assert summary["analysis"]["tests"] == 2 * 3 * 5
    assert telemetry.verify_run((tmp_path / "telemetry.jsonl").read_bytes().splitlines()).ok


def test_finish_reports_too_few_jobs(tmp_path):
    write_ledger(tmp_path, fake_jobs(12, [3, 7], np.random.default_rng(4)), done=9)
    summary = dc.finish(tmp_path, dc.HardwareDesign(replicas=49))
    assert summary["verdict"].startswith("INCOMPLETE: 9 jobs")
    assert not (tmp_path / "telemetry.jsonl").exists()


def test_control_run_on_a_fake_backend():
    from qiskit_ibm_runtime.fake_provider import FakeFez
    backend = FakeFez()
    layout, readout = dc.choose_layout(backend, 2)
    assert len(layout) == 2 and all(0 <= e < 0.5 for e in readout.values())
    design = dc.HardwareDesign(jobs=3, shots=500)
    jobs = dc.control_run(backend, layout, design)
    assert len(jobs) == 3
    for job in jobs:
        for circuit, counts in zip(dc.CIRCUITS, job["counts"]):
            assert sum(counts.values()) == 500
            assert all(e < 250 for e in dc.error_counts(counts, 2, circuit))
