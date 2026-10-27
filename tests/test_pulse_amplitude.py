"""bridge.pulse_amplitude: model, change estimator, H1-H4 tests, ledger path.

Model data, a noiseless simulator and fake IBM objects only; nothing here touches IBM hardware.
"""
import json
from collections import Counter
from datetime import datetime, timezone
from types import SimpleNamespace

import numpy as np
import pytest
from osiris_governance import telemetry

from bridge import drift_check as dc
from bridge import pulse_amplitude as pa
from dnalang.backends import aer, ibm
from dnalang.ledger import Ledger

SMALL = pa.E4Design(replicas=199)        # 2 qubits: Holm's smallest H1 bar 0.05/4 > 1/200


def test_circuits():
    circuits = {c.name: c for c in pa.pulse_circuits(2)}
    assert tuple(circuits) == pa.CIRCUIT_NAMES
    assert circuits["prep1"].count("x") == 2 and circuits["prep0"].count("x") == 0
    assert circuits["amp48"].count("sx") == 2 and circuits["amp48"].count("x") == 96
    assert circuits["sx49"].count("sx") == 98 and circuits["sx49"].count("x") == 0
    assert circuits["ram0"].count("delay") == 2 and circuits["ram0"].count("rz") == 0
    assert circuits["ram90"].count("rz") == 2 and circuits["ram90"].count("sx") == 4
    assert circuits["x100"].count("x") == 200


def test_noiseless_baseline():
    counts = aer.run(pa.pulse_circuits(2), shots=2000, seed=11)
    by_name = dict(zip(pa.CIRCUIT_NAMES, (pa.ones(c, 2) for c in counts)))
    assert by_name["prep0"] == [0, 0] and by_name["prep1"] == [2000, 2000] and by_name["x100"] == [0, 0]
    assert by_name["ram0"] == [2000, 2000]                     # sx, delay, sx is x without detuning
    for name in [f"amp{n}" for n in pa.X_REPS] + [f"sx{m}" for m in pa.SX_REPS] + ["ram90"]:
        assert all(abs(v - 1000) < 120 for v in by_name[name])


def _rx(theta):
    return np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)], [-1j * np.sin(theta / 2), np.cos(theta / 2)]])


def test_model_matches_exact_rotations():
    delta, eps = 0.031, 0.0072
    sx, x = _rx(np.pi / 2 + delta), _rx(np.pi + eps)
    zero = np.array([1, 0], complex)
    sequences = [[], [x]] + [[sx] + [x] * n for n in pa.X_REPS] + [[sx] * m for m in pa.SX_REPS]
    exact = []
    for gates in sequences:
        state = zero
        for g in gates:
            state = g @ state
        exact.append(abs(state[1]) ** 2)
    assert pa.p_one(0.5, 1.0, 0.0, delta, eps) == pytest.approx(exact, abs=1e-9)


def test_jacobian_matches_finite_differences():
    point = np.array([0.497, 0.98, 0.001, -0.05, 0.008])
    numeric = []
    for i in (0, 1, 3, 4):
        h = np.zeros(5)
        h[i] = 1e-6
        numeric.append((pa.p_one(*(point + h)) - pa.p_one(*(point - h))) / 2e-6)
    assert pa.jacobian(*point) == pytest.approx(np.array(numeric).T, abs=1e-6)


def test_pooled_fit_recovers_the_parameters():
    rng = np.random.default_rng(3)
    p = pa.p_one(0.497, 0.98, 0.0, -0.05, 0.008)
    y = rng.binomial(200_000, np.broadcast_to(p, (10, p.size)))
    c, a, _, delta, eps = pa.fit_pooled(y, 200_000)
    assert (c, a) == (pytest.approx(0.497, abs=2e-3), pytest.approx(0.98, abs=3e-3))
    assert (delta, eps) == (pytest.approx(-0.05, abs=5e-4), pytest.approx(0.008, abs=5e-5))


@pytest.mark.parametrize("i", range(4))
def test_changes_and_their_ratio_survive_the_static_misfit(i):
    design = pa.E4Design(shots=10 ** 6, jobs=10)
    _, _, _, delta0, eps0 = pa.sim_theta(i)
    late = np.arange(10) >= 5
    y = pa.synthetic_counts(i, eps0 + 0.005 * late, delta0 + 0.0025 * late, np.zeros(10), design,
                            np.random.default_rng(4))
    fitted = y[:, :len(pa.FITTED)]
    phi = pa.per_job(fitted / design.shots, pa.fit_pooled(fitted, design.shots), design.shots, design.iterations)
    step_eps = phi[late, 3].mean() - phi[~late, 3].mean()
    step_delta = phi[late, 2].mean() - phi[~late, 2].mean()
    assert step_eps == pytest.approx(0.005, abs=3e-4)
    assert step_delta / step_eps == pytest.approx(0.5, abs=0.1)        # a gain change reads as gain


def test_pooled_fit_keeps_every_probability_physical():
    # perfect readout plus E3's unscaled residuals, which push sx33 (P near 0) below zero before clipping
    misfit = np.array(pa.E3_MISFIT[2] + pa.E3_MISFIT[2][3:6])
    p = np.clip(pa.p_one(0.5, 1.0, 0.0, -0.048, 0.0123) + misfit, 1e-6, 1 - 1e-6)
    y = np.random.default_rng(20).binomial(1000, np.broadcast_to(p, (10, p.size)))
    c, a, *_ = pa.fit_pooled(y, 1000)
    assert c - a / 2 >= 0 and c + a / 2 <= 1


def test_ramsey_frequency_shift_is_recovered():
    design = pa.E4Design(shots=200_000, jobs=10)
    shift = np.where(np.arange(10) >= 5, 300e3, 0.0)
    y = pa.synthetic_counts(0, np.full(10, 0.008), np.full(10, -0.04), shift, design, np.random.default_rng(5))
    c = pa.sim_theta(0)[0]
    phase = pa.ramsey_phase(y[:, 9] / design.shots, y[:, 10] / design.shots, c)
    found = pa.frequency_shift(phase, design.ramsey_delay_s)
    assert np.abs(found - (shift - 150e3)) == pytest.approx(0, abs=5e3)


def test_deming_slope_and_error_ratio():
    rng = np.random.default_rng(6)
    x = np.linspace(0, 0.01, 30)
    xs, ys = x + rng.normal(0, 1e-4, 30), 0.5 * x + rng.normal(0, 3e-4, 30)
    assert pa.deming(xs, ys, 9.0) == pytest.approx(0.5, abs=0.03)
    assert abs(pa.deming(xs, ys, 1.0) - pa.deming(xs, ys, 9.0)) > 1e-3     # the ratio matters
    assert not np.isfinite(pa.deming(np.array([0.0, 1.0, 0.0, 1.0]), np.array([0.0, 0.0, 1.0, 1.0]), 1.0))


@pytest.mark.parametrize("low, high, region", [
    (0.3, 0.8, "gain"), (0.25, 1.0, "gain"), (0.2, 0.8, "inconclusive"), (0.3, 1.1, "inconclusive"),
    (-0.1, 0.2, "x_only"), (-0.1, 0.25, "inconclusive"), (0.05, 0.2, "inconclusive"),
    (1.3, 2.5, "detuning_like"), (1.25, 2.5, "inconclusive"),
])
def test_ratio_regions(low, high, region):
    assert pa.ratio_region(low, high, pa.E4Design()) == region


def test_ratio_test_tells_gain_from_x_only():
    rng = np.random.default_rng(7)
    eps = {q: 0.008 + np.linspace(0, 0.005, 30) + rng.normal(0, 5e-4, 30) for q in range(3)}
    gain = {q: -0.05 + (eps[q] - 0.008) / 2 + rng.normal(0, 5e-4, 30) for q in range(3)}
    still = {q: -0.05 + rng.normal(0, 5e-4, 30) for q in range(3)}
    var = {q: 2.5e-7 for q in range(3)}
    assert pa.ratio_test(eps, gain, var, var, rng, SMALL)["region"] == "gain"
    assert pa.ratio_test(eps, still, var, var, rng, SMALL)["region"] == "x_only"


def test_frequency_test_catches_a_detuning():
    rng = np.random.default_rng(8)
    eps = {0: np.r_[np.zeros(15), np.full(15, 0.005)] + rng.normal(0, 5e-4, 30)}
    quiet = pa.frequency_test({0: rng.normal(0, 25e3, 30)}, eps, rng, SMALL)
    moved = pa.frequency_test({0: np.r_[np.zeros(15), np.full(15, 0.9e6)] + rng.normal(0, 25e3, 30)}, eps, rng, SMALL)
    flat = pa.frequency_test({0: np.r_[np.zeros(15), np.full(15, 0.9e6)] + rng.normal(0, 25e3, 30)},
                             {0: rng.normal(0, 5e-4, 30)}, rng, SMALL)
    assert quiet["passed"] == {0: True} and moved["passed"] == {0: False}
    assert flat["passed"] == {0: False} and flat["max_abs_hz"][0] > 0.8e6     # the limit alone fails it


def test_jump_test_separates_a_step_from_a_walk():
    rng = np.random.default_rng(9)
    t = 150.0 * np.arange(30)
    v = np.full(30, 2.9e-7)
    walk = np.r_[0, np.cumsum(rng.normal(0, np.sqrt(pa.E3_WALK_Q * 150), 29))] + rng.normal(0, np.sqrt(v))
    step = walk + 0.01 * (np.arange(30) >= 17)                  # qubit 82's drop in E3 was 0.011
    quiet = pa.jump_test(walk, t, v, rng, 999)
    found = pa.jump_test(step, t, v, rng, 999)
    assert quiet["p"] > 0.05
    assert found["p"] < 0.01 and found["interval"] == 16


def test_shared_test_needs_common_changes():
    rng = np.random.default_rng(10)
    common = np.cumsum(rng.normal(0, 1e-3, 30))
    shared = np.array([common + rng.normal(0, 3e-4, 30) for _ in range(4)])
    separate = np.array([np.cumsum(rng.normal(0, 1e-3, 30)) for _ in range(4)])
    assert pa.shared_test(shared, rng, 999)["p"] < 0.01
    assert pa.shared_test(separate, rng, 999)["p"] > 0.05


def _result(h1=True, region="gain", passed=None, jumps=False, shared=True):
    return {"H1_rotation_error_varies": h1, "H2": {"region": region},
            "H2F": {"passed": {0: True, 1: True} if passed is None else passed}, "H3_jumps": jumps, "H4_shared": shared}


@pytest.mark.parametrize("kwargs, start", [
    ({"h1": False}, "not evaluable"),
    ({"region": "detuning_like"}, "frequency, not amplitude"),
    ({"passed": {0: False, 1: False}}, "frequency, not amplitude"),
    ({}, "shared analog gain drift"),
    ({"jumps": True}, "step changes in a shared gain"),
    ({"shared": False}, "gain change, not shown to be shared"),
    ({"region": "x_only", "jumps": True}, "x-amplitude re-tunes"),
    ({"region": "x_only", "jumps": True, "shared": False}, "x-amplitude jumps, not shown"),
    ({"region": "x_only"}, "an x-only change"),
    ({"region": "inconclusive"}, "amplitude not excluded"),
])
def test_reading_follows_the_table(kwargs, start):
    assert pa.reading(_result(**kwargs)).startswith(start)


def test_reading_names_partly_falsified_qubits():
    assert pa.reading(_result(passed={0: True, 1: False})).endswith("H2F falsifies amplitude on qubits [1]")


def test_analyze_flags_a_rotation_change_on_its_own_qubit_and_parameter():
    rng = np.random.default_rng(11)
    k = SMALL.jobs
    late = np.arange(k) >= k // 2
    y = {}
    for i in range(2):
        _, _, _, delta, eps = pa.sim_theta(i)
        y[i] = pa.synthetic_counts(i, eps + (0.005 * late if i == 0 else 0), np.full(k, delta), np.zeros(k), SMALL, rng)
    result = pa.analyze([150.0 * j for j in range(k)], y, SMALL.shots, SMALL, 12, SMALL.ramsey_delay_s)
    assert result["flagged"] == ["q0/rotation_error/dispersion"]
    assert result["H1_rotation_error_varies"] and not result["HC_contrast_varies"]
    assert set(result["H2F"]["passed"]) == {0} and result["H3_jumps"]
    assert result["H3"]["per_qubit"][0]["interval"] == k // 2 - 1


def test_analyze_stops_when_nothing_varies():
    rng = np.random.default_rng(13)
    k = SMALL.jobs
    y = {i: pa.synthetic_counts(i, np.full(k, pa.sim_theta(i)[4]), np.full(k, pa.sim_theta(i)[3]), np.zeros(k),
                                SMALL, rng) for i in range(2)}
    result = pa.analyze([150.0 * j for j in range(k)], y, SMALL.shots, SMALL, 14, SMALL.ramsey_delay_s)
    assert not result["H1_rotation_error_varies"] and "H2" not in result
    assert result["reading"].startswith("not evaluable")


def test_h3_reports_underpowered_when_the_queue_spreads_jobs():
    rng = np.random.default_rng(15)
    k = SMALL.jobs
    late = np.arange(k) >= k // 2
    _, _, _, delta, eps = pa.sim_theta(0)
    y = {0: pa.synthetic_counts(0, eps + 0.005 * late, np.full(k, delta), np.zeros(k), SMALL, rng)}
    result = pa.analyze([600.0 * j for j in range(k)], y, SMALL.shots, SMALL, 16, SMALL.ramsey_delay_s)
    assert result["H3"] == {"underpowered": True, "short_intervals": 0} and not result["H3_jumps"]


def test_software_control_structure():
    control = pa.software_control(pa.E4Design(replicas=18), qubits=2)
    assert control["H1_tests"] == 4 and control["dispatch_allowed"] is True     # 18 replicas: p >= 1/19


# --------------------------------------------------------------------------------------------
# Hardware path, offline
# --------------------------------------------------------------------------------------------

DATE = datetime(2026, 10, 26, 9, 0, tzinfo=timezone.utc)


class FakeProps:
    last_update_date = DATE

    def to_dict(self):
        return {"backend_name": "fake_fez", "last_update_date": str(DATE)}

    def qubit_property(self, q):
        return {"T1": (1e-4, DATE), "frequency": (5e9 + q, DATE), "readout_error": (0.01, DATE)}

    def gate_property(self, gate, qubits):
        if gate not in ("x", "sx"):
            raise KeyError(gate)
        return {"gate_error": (2e-4, DATE), "gate_length": (3.2e-8, DATE)}


class FakeBackend:
    name = "fake_fez"

    def __init__(self, operational=True):
        self.operational = operational
        measure = {(q,): SimpleNamespace(error=0.01 * (q + 1)) for q in range(9)}
        self.target = type("Target", (), {"dt": 0.5e-9, "granularity": 16,
                                          "__getitem__": lambda self, key: measure})()
        self.coupling_map = SimpleNamespace(get_edges=lambda: [(q, q + 1) for q in range(8)])

    def status(self):
        return SimpleNamespace(operational=self.operational)

    def properties(self, refresh=False, datetime=None):
        return FakeProps()


def test_choose_extras_skips_base_and_neighbours():
    assert pa.choose_extras(FakeBackend(), [0, 5], 2) == [2, 7]


def test_aligned_delay_rounds_to_granularity():
    backend = FakeBackend()
    assert pa.aligned_delay(backend, 200e-9) == pytest.approx(200e-9)
    assert pa.aligned_delay(backend, 210e-9) == pytest.approx(208e-9)


def test_snapshot_and_dated_inside():
    snap = pa.snapshot(FakeProps(), [3])
    assert snap["qubits"]["3"]["frequency"] == [5e9 + 3, str(DATE)]
    assert snap["qubits"]["3"]["x_gate_error"] == [2e-4, str(DATE)] and "sx_gate_length" in snap["qubits"]["3"]
    assert len(snap["properties_sha256"]) == 64
    inside = pa.dated_inside([snap], "2026-10-26T08:59:00Z", "2026-10-26T09:01:00Z")
    assert {(f["qubit"], f["field"]) for f in inside} == {("3", n) for n in snap["qubits"]["3"]}
    assert pa.dated_inside([snap], "2026-10-26T09:00:01Z", "2026-10-26T09:01:00Z") == []
    assert pa.dated_inside([snap], "2026-10-26T08:59:00Z", "2026-10-26T08:59:59Z") == []
    assert pa.dated_inside([snap], "2026-10-26T08:59:00Z", "2026-10-26T09:01:00Z", [4]) == []


def _counts(ones_per_qubit, shots, rng):
    bits = np.zeros((shots, len(ones_per_qubit)), int)
    for q, n in enumerate(ones_per_qubit):
        bits[rng.permutation(shots)[:n], q] = 1
    return dict(Counter("".join(str(b) for b in row[::-1]) for row in bits))


def _fake_ibm(monkeypatch, remaining=404.0, operational=True):
    rng = np.random.default_rng(17)
    backend = FakeBackend(operational)
    svc = SimpleNamespace(backend=lambda name: backend, usage=lambda: {"usage_remaining_seconds": remaining})
    submitted = []

    def submit(circuits, backend_name, layout, shots, ledger, tags=None, optimization_level=0):
        intent = ledger.append("submit-intent", {"calibration_hash": "cal-a", "layout": list(layout)})
        job_id = f"job{len(submitted)}"
        ledger.append("submitted", {"intent_hash": intent["hash"], "job_id": job_id})
        submitted.append((job_id, list(layout), len(circuits)))
        return job_id

    def collect(svc, ledger, job_ids):
        k = len(job_ids)
        late = np.arange(k) >= k // 2
        per_qubit = [pa.synthetic_counts(i, pa.sim_theta(i)[4] + (0.005 * late if i == 0 else 0),
                                         np.full(k, pa.sim_theta(i)[3]), np.zeros(k), pa.E4Design(jobs=k), rng)
                     for i in range(4)]
        for j, job_id in enumerate(job_ids):
            counts = [_counts([per_qubit[i][j, c] for i in range(4)], 1000, rng) for c in range(len(pa.CIRCUIT_NAMES))]
            window = [f"2026-10-26T10:{(150 * j) // 60:02d}:{(150 * j) % 60:02d}.000000Z",
                      f"2026-10-26T10:{(150 * j + 5) // 60:02d}:{(150 * j + 5) % 60:02d}.000000Z"]
            ledger.append("result", {"job_id": job_id, "status": "DONE", "counts": counts, "window": window,
                                     "qpu_seconds": 5})

    monkeypatch.setattr(ibm, "service", lambda: svc)
    monkeypatch.setattr(ibm, "submit", submit)
    monkeypatch.setattr(dc, "_use_key_file", lambda path: None)
    monkeypatch.setattr(dc, "collect", collect)
    monkeypatch.setattr(pa.time, "sleep", lambda seconds: None)
    return submitted


def _paths(tmp_path):
    (tmp_path / "prereg.md").write_text("frozen\n")
    (tmp_path / "e3.json").write_text(json.dumps({"backend": "fake_fez", "layout": [0, 5]}))
    return tmp_path / "run", tmp_path / "key.json", tmp_path / "prereg.md", tmp_path / "e3.json"


HW = pa.E4Design(jobs=20, replicas=18, extra_qubits=2)      # 18 replicas: no test can reach alpha


def test_hardware_writes_a_snapshot_before_every_submission(monkeypatch, tmp_path):
    submitted = _fake_ibm(monkeypatch)
    out, key, prereg, after = _paths(tmp_path)
    summary = pa.hardware(out, key, prereg, after, HW)
    assert summary["jobs_done"] == 20 and summary["qpu_seconds"] == 100
    assert all(layout == [0, 5, 2, 7] and n == 12 for _, layout, n in submitted)
    rows = list(Ledger(out / "ledger.jsonl"))
    assert Ledger(out / "ledger.jsonl").verify() is None
    assert [r["kind"] for r in rows[:2]] == ["preregistration", "design"]
    intents = [i for i, r in enumerate(rows) if r["kind"] == "submit-intent"]
    assert len(intents) == 20
    for j, i in enumerate(intents):
        assert rows[i - 1]["kind"] == "calibration-snapshot" and rows[i - 1]["job"] == j
    state = json.loads((out / "state.json").read_text())
    assert state["ramsey_delay_s"] == pytest.approx(200e-9) and state["layout"] == [0, 5, 2, 7]
    lines = (out / "telemetry.jsonl").read_bytes().splitlines()
    assert telemetry.verify_run(lines).ok
    analysis = summary["analysis"]
    assert set(analysis["qubits"]) == {"0", "5", "2", "7"} and analysis["H1_tests"] == 8
    assert analysis["qubits"]["0"]["x100_errors"] == [int(v) for v in telemetry.metric_series(lines, "q0_x100_errors")[1]]


@pytest.mark.parametrize("remaining, operational, reason", [
    (199.0, True, "QPU budget below the minimum or unknown"),
    (None, True, "QPU budget below the minimum or unknown"),
    (404.0, False, "fake_fez not operational"),
])
def test_hardware_fails_closed(monkeypatch, tmp_path, remaining, operational, reason):
    submitted = _fake_ibm(monkeypatch, remaining, operational)
    out, key, prereg, after = _paths(tmp_path)
    summary = pa.hardware(out, key, prereg, after, HW)
    assert summary["verdict"].startswith("NOT RUN") and submitted == []
    assert list(Ledger(out / "ledger.jsonl"))[-1]["reason"] == reason
