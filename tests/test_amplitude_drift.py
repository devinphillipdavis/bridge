"""bridge.amplitude_drift: model fits, separation of rotation error from contrast, ledger path.

Model data and a noiseless simulator only; nothing here touches IBM hardware.
"""
import json

import numpy as np
import pytest
from osiris_governance import telemetry

from bridge import amplitude_drift as ad
from dnalang.backends import aer
from dnalang.ledger import Ledger

THETA = (0.5, 0.95, 0.0005, 0.02, 0.0068)           # (c, a, gamma, delta, eps)
ONE_QUBIT = ad.AmpDesign(replicas=199)               # 8 tests: Holm's smallest bar 0.05/8 > 1/200


def test_circuits():
    circuits = ad.amplitude_circuits(2)
    assert tuple(c.name for c in circuits) == ad.CIRCUIT_NAMES
    by_name = {c.name: c for c in circuits}
    assert by_name["prep0"].count("x") == 0 and by_name["prep1"].count("x") == 2
    assert by_name["amp0"].count("sx") == 2 and by_name["amp0"].count("x") == 0
    assert by_name["amp64"].count("sx") == 2 and by_name["amp64"].count("x") == 128
    assert by_name["x100"].count("sx") == 0 and by_name["x100"].count("x") == 200


def test_noiseless_baseline():
    counts = aer.run(ad.amplitude_circuits(2), shots=2000, seed=11)
    by_name = dict(zip(ad.CIRCUIT_NAMES, (ad.ones(c, 2) for c in counts)))
    assert by_name["prep0"] == [0, 0] and by_name["prep1"] == [2000, 2000] and by_name["x100"] == [0, 0]
    for n in ad.REPS:                                 # sx then an even number of x: an equator state
        assert all(abs(v - 1000) < 120 for v in by_name[f"amp{n}"])


def test_model_at_known_points():
    p = ad.p_one(0.5, 1.0, 0.0, 0.0, 0.0)
    assert p[0] == pytest.approx(0, abs=1e-8) and p[1] == pytest.approx(1) and p[2:] == pytest.approx(0.5)
    p = ad.p_one(0.5, 1.0, 0.0, 0.0, np.pi / 200)     # 100 x gates then add a quarter turn
    assert p[-1] == pytest.approx(1.0)


def test_pooled_fit_recovers_the_parameters():
    rng = np.random.default_rng(3)
    y = ad.synthetic(THETA, ad.AmpDesign(shots=200_000), rng)
    c, a, gamma, delta, eps = ad.fit_pooled(y, 200_000)
    assert (c, a, eps) == (pytest.approx(0.5, abs=2e-3), pytest.approx(0.95, abs=3e-3), pytest.approx(0.0068, abs=5e-5))
    assert delta == pytest.approx(0.02, abs=0.01)


def test_pooled_fit_recovers_contrast_decay():
    rng = np.random.default_rng(8)
    y = ad.synthetic((0.5, 0.95, 0.004, 0.02, 0.0068), ad.AmpDesign(shots=200_000), rng)
    gamma = ad.fit_pooled(y, 200_000)[2]
    assert gamma == pytest.approx(0.004, abs=3e-4)


def test_hypotheses_separate_dispersion_from_time_order():
    assert ad.hypotheses([]) == dict.fromkeys(
        ("HA_rotation_error_varies", "HAt_rotation_error_time_order", "HC_contrast_varies", "HCt_contrast_time_order"),
        False)
    h = ad.hypotheses(["q4/rotation_error/cusum_peak", "q54/contrast/dispersion"])
    assert h == {"HA_rotation_error_varies": False, "HAt_rotation_error_time_order": True,
                 "HC_contrast_varies": True, "HCt_contrast_time_order": False}


def test_per_job_coarse_to_fine_equals_the_full_grid():
    y = ad.synthetic(THETA, ad.AmpDesign(), np.random.default_rng(4))
    theta = ad.fit_pooled(y, 2000)
    c, _, gamma, delta, eps = theta
    full = ad._grid_best(y, 2000, c, gamma, delta,
                         eps + np.arange(-ad.EPS_WINDOW, ad.EPS_WINDOW + 1e-12, ad.EPS_STEP),
                         np.arange(0.0, 1.2 + 1e-12, ad.CONTRAST_STEP))
    eps_k, a_k = ad.per_job(y, 2000, theta)
    assert np.allclose(eps_k, full[0]) and np.allclose(a_k, full[1])


@pytest.mark.parametrize("steps, expect_rotation, expect_contrast", [
    ((0.0, 0.0), False, False),
    ((0.003, 0.0), True, False),
    ((0.0, -0.1), False, True),
])
def test_each_change_is_flagged_on_its_own_parameter(steps, expect_rotation, expect_contrast):
    rng = np.random.default_rng(5)
    y = {0: ad.synthetic(THETA, ONE_QUBIT, rng, *steps)}
    times = [300 * k for k in range(ONE_QUBIT.jobs)]
    result = ad.analyze(times, y, 2000, ONE_QUBIT, seed=6)
    assert result["tests"] == 8
    assert result["HA_rotation_error_varies"] is expect_rotation
    assert result["HC_contrast_varies"] is expect_contrast
    assert not result["qubits"]["0"]["at_window_edge"]


def test_software_control_structure():
    control = ad.software_control(ad.AmpDesign(jobs=10, replicas=19))
    assert set(control) == {"constant", "eps_step", "contrast_step", "dispatch_allowed"}
    assert control["dispatch_allowed"] is True               # 19 replicas cannot reach Holm's bar


def two_qubit_counts(n0, n1, shots):
    """A counts dict whose qubit marginals are exactly n0 and n1 ones."""
    both = max(0, n0 + n1 - shots)
    return {k: v for k, v in {"11": both, "10": n1 - both, "01": n0 - both,
                              "00": shots - n0 - n1 + both}.items() if v}


def test_finish_from_ledger_rows(tmp_path):
    rng = np.random.default_rng(7)
    design = ad.AmpDesign(replicas=19)
    y = {q: ad.synthetic(THETA, design, rng) for q in (3, 9)}
    ledger = Ledger(tmp_path / "ledger.jsonl")
    ids = []
    for k in range(design.jobs):
        intent = ledger.append("submit-intent", {"calibration_hash": "cal-a"})
        ledger.append("submitted", {"intent_hash": intent["hash"], "job_id": f"job{k}"})
        ids.append(f"job{k}")
    for k in range(design.jobs):
        counts = [two_qubit_counts(int(y[3][k, i]), int(y[9][k, i]), 2000) for i in range(len(ad.FITTED))]
        counts.append(two_qubit_counts(40 + k, 30, 2000))
        window = [f"2026-10-26T16:{k:02d}:00.000000Z", f"2026-10-26T16:{k:02d}:05.000000Z"]
        ledger.append("result", {"job_id": f"job{k}", "status": "DONE", "counts": counts, "window": window,
                                 "qpu_seconds": 8})
    (tmp_path / "state.json").write_text(json.dumps({"backend": "fake", "layout": [3, 9], "job_ids": ids}))
    summary = ad.finish(tmp_path, design)
    assert summary["jobs_done"] == 12 and summary["qpu_seconds"] == 96
    assert summary["analysis"]["tests"] == 2 * 2 * 4
    assert set(summary["analysis"]["qubits"]) == {"3", "9"}
    lines = (tmp_path / "telemetry.jsonl").read_bytes().splitlines()
    assert telemetry.verify_run(lines).ok
    assert [int(v) for v in telemetry.metric_series(lines, "q3_amp16_ones")[1]] == y[3][:, ad.FITTED.index("amp16")].tolist()
    assert [int(v) for v in telemetry.metric_series(lines, "q9_x100_errors")[1]] == [30] * 12
