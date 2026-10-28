"""bridge.gain_drift: E5's tests and reading, the one-in-flight scheduler, resume, ledger path.

Model data and fake IBM objects only; nothing here touches IBM hardware.
"""
import json
import time
from collections import Counter
from dataclasses import replace
from datetime import datetime, timedelta, timezone
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.stats import chi2, spearmanr

from bridge import drift_check as dc
from bridge import gain_drift as gd
from bridge import pulse_amplitude as pa
from dnalang.backends import ibm
from dnalang.ledger import Ledger

SMALL = gd.E5Design(replicas=199)


# --------------------------------------------------------------------------------------------
# Synthetic data
# --------------------------------------------------------------------------------------------

def test_synthetic_counts_match_e4_without_a_contrast_change():
    args = (np.full(6, 0.008), np.full(6, -0.04), np.linspace(0, 1e5, 6))
    ours = gd.synthetic_counts(1, *args, gd.E5Design(), np.random.default_rng(1))
    theirs = pa.synthetic_counts(1, *args, pa.E4Design(jobs=6), np.random.default_rng(1))
    assert (ours == theirs).all()


def test_a_rise_in_e1_lowers_prep1_only():
    design = gd.E5Design(shots=10 ** 6)
    rise = np.r_[np.zeros(3), np.full(3, 0.1)]
    y = gd.synthetic_counts(0, np.full(6, 0.008), np.full(6, -0.04), np.zeros(6), design,
                            np.random.default_rng(2), rise) / design.shots
    assert y[3:, 0].mean() == pytest.approx(y[:3, 0].mean(), abs=1e-3)             # prep0: c - a/2 unchanged
    assert y[:3, 1].mean() - y[3:, 1].mean() == pytest.approx(0.1, abs=2e-3)        # prep1 falls by the rise


# --------------------------------------------------------------------------------------------
# Joint gamma
# --------------------------------------------------------------------------------------------

def _drift(i, k=24, span=0.032, gamma=0.0):
    c, a, _, delta0, eps0 = pa.sim_theta(i)
    path = np.linspace(0, span, k)
    return path, eps0, pa.p_one(c, a, gamma, delta0 + path / 2, eps0 + path)


def test_fit_joint_does_not_take_a_drift_for_decay():
    shots = 10 ** 9
    path, eps0, p = _drift(0)
    theta, pooled = gd.fit_joint(p * shots, shots, gd.E5Design())
    assert pooled > 1e-4 and theta[2] < 1e-6                             # E4's pooled fit reads the drift as decay
    phi = pa.per_job(p, theta, shots, 5)
    assert np.abs(phi[:, 3] - eps0 - path).max() < 5e-6
    assert np.polyfit(phi[:, 3], phi[:, 2], 1)[0] == pytest.approx(0.5, abs=1e-3)
    bent = pa.per_job(p, np.r_[theta[:2], pooled, theta[3:]], shots, 5)
    assert np.abs(bent[:, 3] - eps0 - path).max() > 5e-4


@pytest.mark.parametrize("gamma", [0.0, 2e-4, 1e-3])
def test_fit_joint_keeps_a_real_decay(gamma):
    shots = 10 ** 9
    _, _, p = _drift(1, span=0.0, gamma=gamma)
    theta, _ = gd.fit_joint(p * shots, shots, gd.E5Design())
    assert theta[2] == pytest.approx(gamma, abs=1e-7)
    assert gamma > 0 or theta[2] == 0.0                                  # the bound itself, not 1e-9 above it


def test_estimate_uses_the_joint_gamma_for_jobs_and_replicas():
    design = gd.E5Design(shots=10 ** 6, replicas=9)
    path, eps0, p = _drift(2)
    rng = np.random.default_rng(21)
    c, a = pa.sim_theta(2)[:2]
    ram = np.broadcast_to([c + a / 2, c], (24, 2))
    y = rng.binomial(design.shots, np.concatenate([p, ram, np.full((24, 1), c - a / 2)], axis=1))
    found = gd.estimate(y, design.shots, design, design.ramsey_delay_s, np.random.default_rng(22))
    theta, pooled = gd.fit_joint(y[:, :len(pa.FITTED)], design.shots, design)
    assert found["theta"][2] == theta[2] < 1e-5 < pooled == found["pooled_gamma"]
    assert np.abs(found["phi"][:, 3] - eps0 - path).max() < 2.5e-4      # binomial noise at 10^6 shots
    freq = y[:, :len(pa.FITTED)] / design.shots
    draws = np.random.default_rng(22).binomial(design.shots, freq.mean(0), size=(9,) + freq.shape) / design.shots
    assert np.array_equal(found["replicas"], pa.per_job(draws, theta, design.shots, design.iterations))
    assert found["replicas"].shape == (9, 24, 4)


# --------------------------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------------------------

def test_frequency_test_correlates_increments_not_levels():
    rng = np.random.default_rng(3)
    k = 24
    trend = np.linspace(0, 1, k)
    eps = {0: 0.004 * trend + rng.normal(0, 5e-4, k)}
    shift = {0: 150e3 * trend - 75e3 + rng.normal(0, 10e3, k)}
    assert spearmanr(eps[0], shift[0])[0] > 0.6                       # levels share the trend
    assert gd.frequency_test(shift, eps, rng, SMALL)["failed"] == []
    walk = np.cumsum(rng.normal(0, 1, k))
    together = gd.frequency_test({0: 10e3 * (walk - walk.mean())}, {0: 5e-4 * walk + rng.normal(0, 1e-4, k)},
                                 rng, SMALL)
    assert together["failed"] == [0] and together["max_abs_hz"][0] < 200e3        # the correlation alone fails it
    assert gd.frequency_test({0: shift[0] + 250e3}, eps, rng, SMALL)["failed"] == [0]  # the limit alone fails it


def test_sharing_test_sees_a_common_trend_that_e4s_test_could_not():
    rng = np.random.default_rng(4)
    k = 24
    common = np.linspace(0, 0.008, k)
    x = np.array([common + rng.normal(0, 5e-4, k) for _ in range(10)])
    assert gd.sharing_test(np.diff(x, axis=-1), rng, 999)["p"] < 0.01
    assert pa.shared_test(x, rng, 999)["p"] > 0.05                    # centring removes a common trend
    separate = np.array([np.cumsum(rng.normal(0, 1e-3, k)) for _ in range(10)])
    assert gd.sharing_test(np.diff(separate, axis=-1), rng, 999)["p"] > 0.05
    assert gd.sharing_test(np.diff(x[:1], axis=-1), rng, 999) == {"mean_correlation": None, "p": 1.0}


def _noise(rng, *shape, k=24, sd=5e-4):
    return np.diff(rng.normal(0, sd, shape + (k,)), axis=-1)


def test_equality_test_ignores_a_common_component_and_sees_unequal_loadings():
    rng = np.random.default_rng(5)
    rep = _noise(rng, 499, 8)
    v = rep.var(0, ddof=1).mean(-1)
    common = rng.normal(0, 3e-3, 23)
    base = _noise(rng, 8)
    equal = base + common
    assert gd.equality_test(equal, rep, v)["statistic"] == pytest.approx(gd.equality_test(base, rep, v)["statistic"])
    assert gd.equality_test(equal, rep, v)["p"] > 0.05
    spread = base + np.linspace(0.5, 1.5, 8)[:, None] * common
    assert gd.equality_test(spread, rep, v)["p"] < 0.01
    assert gd.loading_test(spread, rep, v)["p"] < 0.01 and gd.loading_test(equal, rep, v)["p"] > 0.05


def test_equality_statistic_and_p_follow_the_formula():
    rng = np.random.default_rng(19)
    v = np.array([1.0, 4.0, 0.25])
    d = rng.normal(0, np.sqrt(v)[:, None], (3, 5))
    rep = rng.normal(0, np.sqrt(v)[None, :, None], (60, 3, 5))

    def t(x):
        mean = [sum(x[q, k] / v[q] for q in range(3)) / sum(1 / v) for k in range(5)]
        return sum((x[q, k] - mean[k]) ** 2 / v[q] for q in range(3) for k in range(5))

    found = gd.equality_test(d, rep, v)
    null = [t(r) for r in rep]
    assert found["statistic"] == pytest.approx(t(d))
    assert found["p"] == pytest.approx((1 + sum(n >= t(d) for n in null)) / 61)


def test_loading_test_holds_its_level_along_a_smooth_common_path():
    ps = []
    for seed in range(20):
        rng = np.random.default_rng(100 + seed)
        rep = _noise(rng, 499, 8)
        ps.append(gd.loading_test(_noise(rng, 8) + 0.016 / 23, rep, rep.var(0, ddof=1).mean(-1))["p"])
    assert 0.35 < np.mean(ps) < 0.65


def test_loadings_are_shares_of_the_others_mean():
    rng = np.random.default_rng(6)
    load = np.array([0.5, 1.0, 1.5, 1.0])
    common = rng.normal(0, 3e-3, 23)
    d = load[:, None] * common + rng.normal(0, 1e-6, (4, 23))
    found = gd.loadings(d, np.ones(4))
    for q in range(4):
        assert found[q] == pytest.approx(load[q] / np.delete(load, q).mean(), rel=1e-3)


def test_common_series_starts_at_zero_and_carries_the_replica_variance():
    rng = np.random.default_rng(7)
    levels = np.array([0.008, 0.005, 0.012])[:, None] + np.linspace(0, 0.004, 24)
    reps = rng.normal(0, 5e-4, (2000, 3, 24))
    series, var = gd.common_series(levels, reps, np.ones(3))
    assert series[0] == 0 and series[-1] == pytest.approx(0.004)
    assert var.mean() == pytest.approx(5e-4 ** 2 / 3, rel=0.05)


def test_cochran_tells_equal_betas_from_unequal_ones():
    same = {"per_qubit": {0: 0.5, 1: 0.52, 2: 0.48}, "weights": {0: 100.0, 1: 100.0, 2: 100.0}, "pooled": 0.5}
    differ = {"per_qubit": {0: 0.2, 1: 0.5, 2: 0.8}, "weights": {0: 1e3, 1: 1e3, 2: 1e3}, "pooled": 0.5}
    assert gd.cochran(same, 0.05)["homogeneous"] and gd.cochran(same, 0.05)["cochran_df"] == 2
    assert gd.cochran(same, 0.05)["cochran_p"] == pytest.approx(chi2.sf(0.08, 2))
    near = {"per_qubit": {0: 0.45, 1: 0.5, 2: 0.55}, "weights": {0: 1e3, 1: 1e3, 2: 1e3}, "pooled": 0.5}
    # Q = 5.0: p = 0.082 on 2 df
    assert gd.cochran(near, 0.10)["homogeneous"] is False and gd.cochran(near, 0.08)["homogeneous"] is True
    assert not gd.cochran(differ, 0.05)["homogeneous"]
    assert gd.cochran({"per_qubit": {0: 0.5}, "weights": {0: 1.0}, "pooled": 0.5}, 0.05)["cochran_q"] is None


def test_net_change_slope_separates_re_tunes_from_a_common_change():
    rng = np.random.default_rng(8)
    start = np.linspace(0.002, 0.012, 8)

    def series(change):
        return {q: np.r_[np.full(12, s), np.full(12, s + c)] + rng.normal(0, 1e-4, 24)
                for q, (s, c) in enumerate(zip(start, change))}

    retune = gd.net_change(series(-0.6 * start), rng, 999)
    common = gd.net_change(series(np.full(8, 0.0017)), rng, 999)
    assert retune["slope"] == pytest.approx(-0.6, abs=0.05) and retune["interval"][1] < -0.4
    assert common["slope"] == pytest.approx(0, abs=0.05) and -0.1 < common["interval"][0] < common["interval"][1] < 0.1


def test_large_moves_finds_the_jump():
    rng = np.random.default_rng(9)
    series = np.r_[0, np.cumsum(rng.normal(0, 1e-3, 23))]
    series[12:] += 0.02
    assert gd.large_moves(series, 3.0) == [11]
    assert gd.large_moves(np.zeros(24), 3.0) == []


def test_distances_and_the_distance_check():
    edges = [(0, 1), (1, 2), (2, 3), (3, 4)]
    assert gd._distances(edges, [0, 3])[(0, 3)] == 3
    rng = np.random.default_rng(10)
    parts = rng.normal(0, 1, (6, 200))
    residuals = {q: parts[q] + parts[q + 1] for q in range(5)}          # neighbours share a part
    found = gd.distance_check(residuals, edges)
    assert found["pairs"] == 10 and found["spearman"] < -0.5


def _result(gate=True, h5=(), hc=(), n=16, region="gain", homogeneous=True, h2=True, h3=False, under=False, h4=False):
    return {"gate": gate, "n_qubits": n, "excluded": {"H5": list(h5), "HC": list(hc)},
            "H1": {"region": region, "homogeneous": homogeneous}, "H2": {"flagged": h2},
            "H3": {"rejected": h3, "underpowered": under}, "H4": {"flagged": h4}}


@pytest.mark.parametrize("kwargs, key", [
    ({"gate": False, "h5": range(9)}, "not_evaluable"),
    ({"h5": range(9), "region": "x_only"}, "frequency"),
    ({"h5": range(8)}, "fully_common"),                       # exactly half is not more than half
    ({"region": "x_only"}, "not_gain"),
    ({"region": "inconclusive", "homogeneous": False}, "not_gain"),
    ({"homogeneous": False, "h2": False}, "gain_heterogeneous"),
    ({"h2": False}, "independent"),
    ({"h3": True, "under": True}, "partly_shared"),           # a rejection needs no power condition
    ({"under": True, "h4": True}, "shared_untested"),
    ({"h4": True}, "fully_common_steps"),
    ({}, "fully_common"),
])
def test_reading_follows_the_table(kwargs, key):
    result = _result(**kwargs)
    assert gd.category(result) == key
    assert gd.reading(result).startswith(gd.READINGS[key])


def test_reading_names_the_region_and_the_excluded_qubits():
    assert gd.reading(_result(region="detuning_like")) == "not a gain change: the ratio is detuning-like"
    text = gd.reading(_result(h5=[82], hc=[4, 22]))
    assert text.endswith("; excluded by H5: qubits [82]; excluded by HC: qubits [4, 22]")


# --------------------------------------------------------------------------------------------
# Analysis, control and simulation
# --------------------------------------------------------------------------------------------

def _walk_counts(qubits, design, seed, extra=None, k=24, readout=None):
    rng = np.random.default_rng(seed)
    walk = np.r_[0, np.cumsum(rng.normal(0, 1, k - 1))]
    walk *= 0.010 / np.ptp(walk)
    y = {}
    for i in range(qubits):
        _, _, _, delta, eps = pa.sim_theta(i)
        shift = np.zeros(k) if extra is None or extra[0] != i else extra[1]
        rise = 0.0 if readout is None or readout[0] != i else readout[1]
        y[i] = gd.synthetic_counts(i, eps + walk, delta + walk / 2, shift, design, rng, rise)
    return y


def test_analyze_excludes_frequency_and_contrast_changes_and_reads_a_shared_gain():
    design = gd.E5Design(replicas=299, r_min_rad=None)       # 7 qubits: Holm's first bar 0.05 / 14 > 1 / 300
    step = np.where(np.arange(24) >= 12, 200e3, -200e3)
    rise = np.where(np.arange(24) >= 6, np.random.default_rng(20).uniform(0.02, 0.1, 24), 0.0)
    y = _walk_counts(7, design, 11, extra=(5, step), readout=(6, rise))
    result = gd.analyze(600.0 * np.arange(24), y, design.shots, design, 12, design.ramsey_delay_s,
                        edges=[(i, i + 1) for i in range(6)])
    assert result["excluded"] == {"H5": [5], "HC": [6]} and result["primary"] == [0, 1, 2, 3, 4]
    assert result["gate"] and result["varying"] == [0, 1, 2, 3, 4]
    assert result["H1"]["region"] == "gain" and result["H1"]["homogeneous"]
    assert result["H2"]["flagged"] and not result["H3"]["rejected"] and result["H3"]["underpowered"]
    assert result["H3"]["common_range"] == pytest.approx(0.010, abs=0.002)
    assert result["reading"] == ("shared gain; whether it is fully common is untested (H3 underpowered); "
                                 "excluded by H5: qubits [5]; excluded by HC: qubits [6]")
    assert result["H3"]["statistic"] != result["H3"]["spread_test"]["statistic"]   # H3 is the loadings test
    assert result["qubits"]["0"]["pooled"]["gamma"] <= result["qubits"]["0"]["pooled_fit_gamma"]
    secondary = result["secondary"]
    assert set(secondary["per_qubit_jumps"]) == {0, 1, 2, 3, 4} and "distance" in secondary
    assert secondary["all_qubits"]["H2"]["flagged"] and "slope" in secondary["net_change"]
    powered = gd.analyze(600.0 * np.arange(24), y, design.shots, replace(design, r_min_rad=0.005), 12,
                         design.ramsey_delay_s, secondary=False)
    assert powered["reading"].startswith("fully common gain drift")


def test_analyze_stops_at_the_gate_when_nothing_varies():
    k = 24
    rng = np.random.default_rng(17)
    y = {i: gd.synthetic_counts(i, np.full(k, pa.sim_theta(i)[4]), np.full(k, pa.sim_theta(i)[3]), np.zeros(k),
                                SMALL, rng) for i in range(4)}
    result = gd.analyze(600.0 * np.arange(k), y, SMALL.shots, SMALL, 14, SMALL.ramsey_delay_s)
    assert not result["gate"] and "H1" not in result and result["reading"].startswith("not evaluable")
    raw = [record["p"][name]["dispersion"] for record in result["qubits"].values() for name in record["p"]]
    assert min(raw) <= 0.05 and result["flagged"] == []                  # Holm over all 8 dispersion tests


def test_analyze_does_not_depend_on_qubit_order():
    design = gd.E5Design(replicas=99)
    y = _walk_counts(3, design, 15)
    one = gd.analyze(600.0 * np.arange(24), y, design.shots, design, 16, design.ramsey_delay_s, gate=False,
                     secondary=False)
    two = gd.analyze(600.0 * np.arange(24), dict(reversed(list(y.items()))), design.shots, design, 16,
                     design.ramsey_delay_s, gate=False, secondary=False)
    assert one["H3"]["p"] == two["H3"]["p"] and one["qubits"] == two["qubits"]


def test_control_flags_and_software_control():
    flagged = {"flagged": ["q1/contrast/dispersion"], "excluded": {"H5": [2], "HC": [1]},
               "H1": {"region": "gain", "homogeneous": False}, "H2": {"flagged": True},
               "H3": {"rejected": True}, "H4": {"flagged": True}}
    assert gd.control_flags(flagged) == ["q1/contrast/dispersion", "q2/H5", "H1/gain", "H1/heterogeneous",
                                         "H2", "H3", "H4"]
    control = gd.software_control(gd.E5Design(replicas=18), qubits=4)
    assert control["dispersion_tests"] == 8 and control["dispatch_allowed"] == (control["control_flags"] == [])


def _seeds(n):
    return np.random.SeedSequence([1, n]).spawn(3)


def test_sim_truth_scenarios():
    design = gd.E5Design()
    t, s1, _ = gd.sim_truth("S1", design, _seeds(0), 4)
    assert t[1] - t[0] == 600
    inc = np.array([np.diff(s1[i][0]) for i in range(4)])
    assert np.allclose(inc, inc[0]) and np.allclose([np.diff(s1[i][1]) for i in range(4)], inc / 2)
    _, s2, _ = gd.sim_truth("S2", design, _seeds(0), 4)
    assert not np.allclose(np.diff(s2[0][0]), np.diff(s2[1][0]))
    for path in ("walk", "linear"):
        _, s3, _ = gd.sim_truth("S3", design, _seeds(0), 4, 0.010, path)
        share = [np.ptp(s3[i][0] - pa.sim_theta(i)[4]) / 0.010 for i in range(4)]
        assert all(0.5 <= s <= 1.5 for s in share) and len(set(np.round(share, 6))) == 4
    assert np.allclose(np.diff(s3[0][0]), np.diff(s3[0][0])[0])        # linear: equal steps
    _, s4, _ = gd.sim_truth("S4", design, _seeds(0), 4)
    for i in range(4):
        steps = np.flatnonzero(np.diff(s4[i][0]))
        assert 1 <= steps.size <= 2 and abs(s4[i][0][-1]) < abs(s4[i][0][0])
    _, s5, event = gd.sim_truth("S5", design, _seeds(0), 4)
    assert np.ptp(s5[event][2]) == 340e3 and not any(np.ptp(s5[i][2]) for i in range(4) if i != event)
    assert np.ptp((s5[event][0] - s1[event][0])) == pytest.approx(0.005)
    _, s6, event = gd.sim_truth("S6", design, _seeds(0), 4)
    rise = s6[event][3]
    assert rise[0] == 0 and rise[-1] >= 0.02 and np.allclose(s6[event][0], s1[event][0])
    _, s7, _ = gd.sim_truth("S7", design, _seeds(0), 4)
    for i in range(4):
        assert np.flatnonzero(np.diff(s7[i][0])).size == 1 and np.ptp(s7[i][0]) == pytest.approx(0.002)


def test_sim_trial_runs_one_scenario_end_to_end():
    row = gd.sim_trial("S1", 0, gd.E5Design(replicas=18, jobs=12), qubits=3)
    assert row["trial"] == 0 and row["category"] == gd.category(row) and row["event_qubit"] is None


@pytest.mark.parametrize("power, expected", [
    ({0.006: 0.5, 0.008: 0.85, 0.010: 0.79, 0.012: 0.9}, 0.012),
    ({0.006: 0.8, 0.008: 0.9}, 0.006),
    ({0.006: 0.1, 0.008: 0.5}, None),
])
def test_r_min_is_the_smallest_range_with_power_there_and_beyond(power, expected):
    assert gd.r_min_from(power) == expected


def _row(trial, category_inputs, event=None, excluded=None, span=0.009):
    row = {"trial": trial, "gate": True, "n_qubits": 16, "excluded": excluded or {"H5": [], "HC": []}, "flagged": [],
           "event_qubit": event,
           "H1": {"region": "gain", "homogeneous": True}, "H2": {"flagged": True},
           "H3": {"rejected": False, "common_range": span, "underpowered": True, "spread_p": 0.5},
           "H4": {"flagged": False}}
    row["H3"]["rejected"] = category_inputs.get("h3", False)
    return row


def test_summarize_sets_r_min_and_rereads_every_trial(tmp_path):
    def write(name, key, rows, scenario=None, common_range=None):
        (tmp_path / name).write_text(json.dumps({"key": key, "scenario": scenario or key, "common_range": common_range,
                                                 "seconds": 1.0, "rows": rows, "code_sha256": "same"}))

    (tmp_path / "S0.json").write_text(json.dumps({"key": "S0", "scenario": "S0", "seconds": 1.0,
                                                  "control": {"dispatch_allowed": True}, "code_sha256": "same"}))
    for r, rejected in ((0.006, 2), (0.008, 9), (0.010, 10)):
        rows = [_row(n, {"h3": n < rejected}) for n in range(10)]
        write(f"S3_{r}.json", gd.sim_key("S3", r), rows, "S3", r)
    write("S1_0-2.json", "S1", [_row(0, {}, span=0.012), _row(1, {}, span=0.007)])
    write("S1_2-3.json", "S1", [_row(2, {}, span=0.012)])
    write("S5.json", "S5", [_row(0, {}, 3, {"H5": [3], "HC": []}, 0.012), _row(1, {}, 4, span=0.012)])
    summary = gd.summarize(tmp_path)
    assert summary["r_min_rad"] == 0.008 and summary["S3_power"]["walk"]["0.008"] == 0.9
    assert summary["rates"]["S1"]["trials"] == 3
    assert summary["rates"]["S1"]["categories"]["fully_common"] == 2
    assert summary["rates"]["S5"]["event_excluded"] == 0.5
    assert summary["rates"]["S5"]["reading_unchanged"] == 0.5        # trial 1 moved from untested to fully common
    assert summary["verdict"]["S0_dispatch_allowed"] is True and summary["verdict"]["S7_H4_power_ge_0.80"] is None
    assert json.loads((tmp_path / "oc.json").read_text())["r_min_rad"] == 0.008
    (tmp_path / "S7.json").write_text(json.dumps({"key": "S7", "scenario": "S7", "seconds": 1.0, "rows": [],
                                                  "code_sha256": "other"}))
    with pytest.raises(RuntimeError, match="2 code versions"):
        gd.summarize(tmp_path)


def test_simulate_writes_one_file_per_invocation(tmp_path):
    record = gd.simulate(tmp_path, "S7", (0, 1), gd.E5Design(sim_replicas=18, jobs=8))
    assert record["key"] == "S7" and len(record["rows"]) == 1 and (tmp_path / "S7_0-1.json").exists()
    with pytest.raises(SystemExit):
        gd.main(["simulate", "--out", str(tmp_path), "--scenario", "S3"])


# --------------------------------------------------------------------------------------------
# Hardware path, offline
# --------------------------------------------------------------------------------------------

DATE = datetime(2026, 10, 26, 9, 0, tzinfo=timezone.utc)


class FakeProps:
    last_update_date = DATE

    def to_dict(self):
        return {"backend_name": "fake_fez", "last_update_date": str(DATE)}

    def qubit_property(self, q):
        return {"T1": (1e-4, DATE), "frequency": (5e9 + q, DATE)}

    def gate_property(self, gate, qubits):
        return {"gate_error": (2e-4, DATE)}


class FakeBackend:
    name = "fake_fez"

    def __init__(self, operational=True, dated=None):
        self.operational, self.dated = operational, dated
        measure = {(q,): SimpleNamespace(error=0.01 * (q + 1)) for q in range(9)}
        self.target = type("Target", (), {"dt": 0.5e-9, "granularity": 16,
                                          "__getitem__": lambda self, key: measure})()
        self.coupling_map = SimpleNamespace(get_edges=lambda: [(q, q + 1) for q in range(8)])

    def status(self):
        return SimpleNamespace(operational=self.operational)

    def properties(self, refresh=False, datetime=None):
        if datetime is not None and self.dated is not None:
            return self.dated(datetime)
        return FakeProps()


class Clock:
    def __init__(self):
        self.now = time.time()
        self.start = self.now

    def time(self):
        return self.now

    def sleep(self, seconds):
        self.now += seconds


def _counts(ones_per_qubit, shots, rng):
    bits = np.zeros((shots, len(ones_per_qubit)), int)
    for q, n in enumerate(ones_per_qubit):
        bits[rng.permutation(shots)[:n], q] = 1
    return dict(Counter("".join(str(b) for b in row[::-1]) for row in bits))


class FakeIBM:
    """Each job runs queue(j) seconds after submission and uses 5 QPU-seconds."""

    def __init__(self, clock, jobs, queue=lambda j: 60.0, budget=400.0, fail=None):
        self.clock, self.queue, self.budget, self.fail = clock, queue, budget, fail or {}
        self.submitted, self.usage_calls = {}, 0
        rng = np.random.default_rng(17)
        per_qubit = _walk_counts(4, gd.E5Design(), 18, k=jobs)
        self.counts = [[_counts([per_qubit[i][j, c] for i in range(4)], 1000, rng) for c in range(len(pa.CIRCUIT_NAMES))]
                       for j in range(jobs)]

    def remaining(self):
        self.usage_calls += 1
        value = self.budget(self) if callable(self.budget) else self.budget
        if value == "raise":
            raise ConnectionError("usage")
        return value

    def submit(self, circuits, backend_name, layout, shots, ledger, tags=None, optimization_level=0):
        j = tags["job"]
        if self.fail.get(j) == "before":
            del self.fail[j]
            raise ConnectionError("before the intent")
        intent = ledger.append("submit-intent", {"calibration_hash": "cal-a", "layout": list(layout), "tags": tags})
        if self.fail.get(j) in ("after", "crash"):
            raise (KeyboardInterrupt if self.fail.pop(j) == "crash" else ConnectionError)("after the intent")
        job_id = f"job{j}"
        ledger.append("submitted", {"intent_hash": intent["hash"], "job_id": job_id})
        self.submitted[job_id] = (j, self.clock.now, list(layout), len(circuits))
        return job_id

    def collect(self, svc, ledger, job_ids):
        done = dc._ledger_results(ledger)
        for job_id in job_ids:
            j, at, _, _ = self.submitted[job_id]
            end = at + self.queue(j)
            if job_id in done or self.clock.now < end:
                continue
            stamp = [dc._stamp(DATE + timedelta(seconds=s - self.clock.start)) for s in (end - 5, end)]   # in the past
            ledger.append("result", {"job_id": job_id, "status": "DONE", "counts": self.counts[j], "window": stamp,
                                     "qpu_seconds": 5})


def _fake_ibm(monkeypatch, jobs=20, operational=True, dated=None, **kwargs):
    clock = Clock()
    fake = FakeIBM(clock, jobs, **kwargs)
    backend = FakeBackend(operational, dated)
    svc = SimpleNamespace(backend=lambda name: backend, usage=lambda: {"usage_remaining_seconds": fake.remaining()})
    monkeypatch.setattr(ibm, "service", lambda: svc)
    monkeypatch.setattr(ibm, "submit", fake.submit)
    monkeypatch.setattr(dc, "_use_key_file", lambda path: None)
    monkeypatch.setattr(dc, "collect", fake.collect)
    monkeypatch.setattr(gd, "_clock", clock.time)
    monkeypatch.setattr(gd, "_sleep", clock.sleep)
    return fake


def _paths(tmp_path):
    (tmp_path / "prereg.md").write_text("frozen\n")
    (tmp_path / "e4.json").write_text(json.dumps({"backend": "fake_fez", "layout": [0, 5]}))
    return tmp_path / "run", tmp_path / "key.json", tmp_path / "prereg.md", tmp_path / "e4.json"


HW = gd.E5Design(jobs=20, replicas=18, extra_qubits=2, poll_s=60.0)     # 18 replicas: no replica test reaches alpha


def _submissions(fake):
    return sorted((j, at) for j, at, _, _ in fake.submitted.values())


def test_hardware_keeps_one_job_in_flight_with_snapshots_on_both_sides(monkeypatch, tmp_path):
    fake = _fake_ibm(monkeypatch, queue=lambda j: 1500.0 if j == 3 else 60.0)
    out, key, prereg, after = _paths(tmp_path)
    summary = gd.hardware(out, key, prereg, after, HW)
    assert summary["jobs_done"] == 20 and summary["qpu_seconds"] == 100 and summary["stopped"] == "all jobs submitted"
    assert summary["long_intervals"] == 18 and "analysis" in summary      # job 4 ran soon after the delayed job 3
    assert all(layout == [0, 5, 2, 7] and n == 12 for _, _, layout, n in fake.submitted.values())
    times = _submissions(fake)
    for (j, at), (_, nxt) in zip(times, times[1:]):
        assert nxt - at >= HW.spacing_s and nxt >= at + fake.queue(j)          # spaced, and after the result
    assert times[4][1] - times[3][1] >= 1500
    rows = list(Ledger(out / "ledger.jsonl"))
    assert Ledger(out / "ledger.jsonl").verify() is None
    assert [r["kind"] for r in rows[:2]] == ["preregistration", "design"]
    assert rows[0]["code_sha256"] == gd.code_hashes()
    kinds = [r["kind"] for r in rows]
    for i, r in enumerate(rows):
        if r["kind"] == "submit-intent":
            before = rows[i - 1]
            assert before["kind"] == "calibration-snapshot" and before["when"] == "before"
            assert before["job"] == r["tags"]["job"] and before["usage_remaining_seconds"] == 400.0
        if r["kind"] == "result":
            nxt = rows[i + 1]
            assert nxt["kind"] == "calibration-snapshot" and nxt["when"] == "after" and nxt["job_id"] == r["job_id"]
    assert kinds.count("submit-intent") == 20 and kinds.count("stopped") == 1
    state = json.loads((out / "state.json").read_text())
    assert state["coupling_edges"][0] == [0, 1] and state["next_job"] == 20 and state["lost"] == []
    assert (out / "control.json").exists() and (out / "telemetry.jsonl").exists()
    analysis = json.loads((out / "analysis.json").read_text())["analysis"]
    assert set(analysis["qubits"]) == {"0", "5", "2", "7"} and analysis["dispersion_tests"] == 8


def test_hardware_stops_at_the_wall_clock_cap(monkeypatch, tmp_path):
    fake = _fake_ibm(monkeypatch, queue=lambda j: 3000.0)
    out, key, prereg, after = _paths(tmp_path)
    summary = gd.hardware(out, key, prereg, after, HW)
    times = [at for _, at in _submissions(fake)]
    assert times[-1] - times[0] <= HW.cap_s and len(times) < 20
    assert summary["stopped"] == "wall-clock cap" and summary["jobs_pending"] == []
    assert summary["verdict"].startswith(f"INCOMPLETE: {len(times)} jobs finished (minimum 16)")


def test_hardware_stops_when_the_budget_runs_low(monkeypatch, tmp_path):
    fake = _fake_ibm(monkeypatch, budget=lambda f: 200.0 if f.usage_calls == 1 else 60.0 - 5 * len(f.submitted))
    out, key, prereg, after = _paths(tmp_path)
    summary = gd.hardware(out, key, prereg, after, HW)
    assert len(fake.submitted) == 9 and summary["stopped"] == "QPU budget below the stop level"
    snaps = [r for r in Ledger(out / "ledger.jsonl") if r.get("when") == "before"]
    assert min(s["usage_remaining_seconds"] for s in snaps) >= HW.stop_budget_s


def test_hardware_waits_while_the_budget_is_unknown(monkeypatch, tmp_path):
    fake = _fake_ibm(monkeypatch, budget=lambda f: "raise" if 2 <= f.usage_calls <= 5 else 400.0)
    out, key, prereg, after = _paths(tmp_path)
    gd.hardware(out, key, prereg, after, HW)
    times = _submissions(fake)
    assert len(times) == 20 and times[0][1] - fake.clock.start >= 4 * HW.poll_s


def test_hardware_retries_before_the_intent_and_loses_the_job_after_it(monkeypatch, tmp_path):
    fake = _fake_ibm(monkeypatch, fail={2: "before", 4: "after"})
    out, key, prereg, after = _paths(tmp_path)
    summary = gd.hardware(out, key, prereg, after, HW)
    rows = list(Ledger(out / "ledger.jsonl"))
    intents = Counter(r["tags"]["job"] for r in rows if r["kind"] == "submit-intent")
    assert intents == Counter(range(20)) and "job2" in fake.submitted and "job4" not in fake.submitted
    assert [r["job"] for r in rows if r["kind"] == "lost"] == [4]
    assert summary["jobs_lost"] == 1 and summary["jobs_done"] == 19


@pytest.mark.parametrize("budget, operational, reason", [
    (149.0, True, "QPU budget below the minimum or unknown"),
    ("raise", True, "QPU budget below the minimum or unknown"),
    (400.0, False, "fake_fez not operational"),
])
def test_hardware_fails_closed(monkeypatch, tmp_path, budget, operational, reason):
    fake = _fake_ibm(monkeypatch, budget=budget, operational=operational)
    out, key, prereg, after = _paths(tmp_path)
    summary = gd.hardware(out, key, prereg, after, HW)
    assert summary["verdict"].startswith("NOT RUN") and fake.submitted == {}
    assert list(Ledger(out / "ledger.jsonl"))[-1]["reason"] == reason


def test_hardware_does_not_dispatch_when_the_control_flags(monkeypatch, tmp_path):
    fake = _fake_ibm(monkeypatch)
    monkeypatch.setattr(gd, "software_control", lambda design, qubits: {"dispatch_allowed": False,
                                                                        "control_flags": ["H2"]})
    out, key, prereg, after = _paths(tmp_path)
    assert gd.hardware(out, key, prereg, after, HW)["control"] == ["H2"] and fake.submitted == {}


def test_hardware_refuses_a_ledger_that_has_rows(monkeypatch, tmp_path):
    _fake_ibm(monkeypatch, budget=149.0)
    out, key, prereg, after = _paths(tmp_path)
    gd.hardware(out, key, prereg, after, HW)
    with pytest.raises(RuntimeError, match="use resume"):
        gd.hardware(out, key, prereg, after, HW)


def test_resume_records_a_lost_job_and_finishes_the_schedule(monkeypatch, tmp_path):
    fake = _fake_ibm(monkeypatch, fail={5: "crash"})
    out, key, prereg, after = _paths(tmp_path)
    with pytest.raises(KeyboardInterrupt):
        gd.hardware(out, key, prereg, after, HW)
    assert len(fake.submitted) == 5
    summary = gd.resume(out, key, HW)
    rows = list(Ledger(out / "ledger.jsonl"))
    intents = Counter(r["tags"]["job"] for r in rows if r["kind"] == "submit-intent")
    assert intents == Counter(range(20))                                    # job 5 never resubmitted
    assert [r["job"] for r in rows if r["kind"] == "lost"] == [5]
    assert [r["kind"] for r in rows].count("resumed") == 1
    assert summary["jobs_done"] == 19 and summary["jobs_lost"] == 1 and summary["stopped"] == "all jobs submitted"
    gd.finish(out, HW)                                                      # idempotent from the ledger


def test_resume_refuses_changed_code(monkeypatch, tmp_path):
    _fake_ibm(monkeypatch, fail={1: "crash"})
    out, key, prereg, after = _paths(tmp_path)
    with pytest.raises(KeyboardInterrupt):
        gd.hardware(out, key, prereg, after, HW)
    monkeypatch.setattr(gd, "code_hashes", lambda: {"bridge/gain_drift.py": "0" * 64})
    with pytest.raises(RuntimeError, match="code differs"):
        gd.resume(out, key, HW)


def test_finish_is_incomplete_when_jobs_run_too_close(monkeypatch, tmp_path):
    _fake_ibm(monkeypatch)
    out, key, prereg, after = _paths(tmp_path)
    summary = gd.hardware(out, key, prereg, after, replace(HW, spacing_s=120.0))
    assert summary["jobs_done"] == 20 and summary["long_intervals"] < 12 and summary["verdict"].startswith("INCOMPLETE")


def test_retrospective_queries_each_job_midpoint(monkeypatch, tmp_path):
    asked = []

    def dated(moment):
        asked.append(moment)
        return FakeProps()

    _fake_ibm(monkeypatch, dated=dated)
    out, key, prereg, after = _paths(tmp_path)
    gd.hardware(out, key, prereg, after, HW)
    jumps = {"0": {"flagged": True, "interval": 5}, "2": {"flagged": False, "interval": 9}}
    (out / "analysis.json").write_text(json.dumps({"analysis": {
        "H4": {"flagged": True, "interval": 7},
        "secondary": {"large_common_moves": [3], "per_qubit_jumps": jumps}}}))
    record = gd.retrospective(tmp_path / "retro", key, out)
    assert record["available"] and len(record["snapshots"]) == 20 and len(asked) == 20
    windows = [s["at"] for s in record["snapshots"]]
    assert asked[0] == pa._midpoint(windows[0]) and windows == sorted(windows)
    assert set(record["dated_inside"]) == {"3", "5", "7"}
    assert (tmp_path / "retro" / "retrospective.json").exists()


def test_retrospective_stops_at_a_failed_query(monkeypatch, tmp_path):
    def dated(moment):
        raise ValueError("no properties")

    _fake_ibm(monkeypatch, dated=dated)
    out, key, prereg, after = _paths(tmp_path)
    gd.hardware(out, key, prereg, after, HW)
    record = gd.retrospective(tmp_path / "retro", key, out)
    assert record == {**record, "available": False, "error": "ValueError", "snapshots": []}


def test_after_snapshot_is_taken_late_when_a_poll_failed(monkeypatch, tmp_path):
    fake = _fake_ibm(monkeypatch)
    calls = {"n": 0}
    real = gd._after_snapshots

    def flaky(svc, ledger, state):
        calls["n"] += 1
        if calls["n"] == 3:
            raise ConnectionError("properties")
        real(svc, ledger, state)

    monkeypatch.setattr(gd, "_after_snapshots", flaky)
    out, key, prereg, after = _paths(tmp_path)
    gd.hardware(out, key, prereg, after, HW)
    rows = list(Ledger(out / "ledger.jsonl"))
    afters = Counter(r["job"] for r in rows if r.get("when") == "after")
    assert afters == Counter(range(20)) and len(fake.submitted) == 20


def test_timestamps_parse_as_utc():
    assert gd._ts_epoch({"ts": "2026-10-28T00:00:00Z"}) == datetime(2026, 10, 28, tzinfo=timezone.utc).timestamp()
    assert timedelta(seconds=gd.E5Design().cap_s) == timedelta(hours=8)
