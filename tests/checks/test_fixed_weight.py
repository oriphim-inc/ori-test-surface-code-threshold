"""Checks on the fixed-weight (stratified) code-capacity estimator, added by the Sceptic in the
second review (after commit 6fdd6d9). These were written AFTER the stratified results existed, so they
are consistency checks, not pre-registered tolerances. They test (1) the stored f_w against exact
enumeration, (2) agreement with the direct Monte Carlo, (3) the logic of the reported d at p = 0.1 %,
and (4) whether the global-fit Part D entries are compatible with the stratified data.
"""
import itertools
import math

import numpy as np
import pytest

TARGET = 1e-12


@pytest.fixture(scope="module")
def fw(summary):
    o = summary["code_capacity"]["overhead"]
    if "fixed_weight" not in o:
        pytest.skip("no fixed-weight results in summary.json")
    return o["fixed_weight"]


def _exact_f(d, w):
    from surfcode.sim import code_capacity_decoder
    m, H, L = code_capacity_decoder(d, 0.01)
    n = H.shape[1]
    sub = np.array(list(itertools.combinations(range(n), w)))
    err = np.zeros((len(sub), n), np.uint8)
    np.put_along_axis(err, sub, 1, axis=1)
    syn = (err.astype(np.int64) @ H.T.astype(np.int64)) % 2
    act = (err.astype(np.int64) @ L.T.astype(np.int64)) % 2
    pred = m.decode_batch(syn.astype(np.uint8))
    return float(np.any(pred != act, axis=1).mean())


@pytest.mark.parametrize("d,w", [(5, 3), (5, 4), (7, 4)])
def test_stored_fw_matches_exact_enumeration(fw, d, w):
    """Stored f_w (sampled) vs. exact enumeration of every weight-w error; 4 sigma."""
    k, N = fw["f_w"][str(d)][str(w)]
    f = _exact_f(d, w)
    sig = math.sqrt(max(f * (1 - f), 1.0 / N) / N)
    assert abs(k / N - f) <= 4 * sig, f"d={d} w={w}: sampled {k/N:.4e} vs exact {f:.4e}"


def test_fw_zero_below_t_plus_1_is_assumed_not_sampled(fw):
    """Lowest stored weight is t+1 for every d (weights <= t are corrected by MWPM; tested in tests/)."""
    for d, v in fw["f_w"].items():
        assert min(int(w) for w in v) == (int(d) - 1) // 2 + 1


def test_validation_vs_monte_carlo(fw):
    """Every direct-MC point lies inside the stratified interval [lo, hi(+truncation)] widened by 4 MC sigma."""
    assert len(fw["validation_vs_mc"]) >= 10
    bad = []
    for r in fw["validation_vs_mc"]:
        sig = (r["mc_hi"] - r["mc_lo"]) / (2 * 1.96)
        if not (r["fw_lo"] - 4 * sig <= r["mc"] <= r["fw_hi"] + 4 * sig):
            bad.append(r)
    assert not bad, bad


def test_reported_d_at_0p1_percent(fw):
    """At p = 0.1 %: reported d meets 1e-12 at the 95 % upper end with negligible truncation,
    and d-2 does not meet it at the central value."""
    req = fw["required"]["0.001"]
    d = req["d_upper95"]
    e = fw["p_L"]["0.001"]
    assert d is not None
    assert e[str(d)]["hi"] <= TARGET and e[str(d)]["trunc"] < 0.01 * TARGET
    assert e[str(d - 2)]["p_L"] > TARGET
    assert req["full"] == 2 * d * d - 1


def _d_from_stratified_extrapolation(fw, p):
    """d obtained by extrapolating the stratified p_L linearly in d (constant Lambda) from the
    d >= 7 points whose truncation is < 10 % of p_L."""
    e = fw["p_L"][p]
    ds = [int(x) for x, v in e.items() if int(x) >= 7 and v["trunc"] < 0.1 * v["p_L"]]
    assert len(ds) >= 3
    y = np.log([e[str(x)]["p_L"] for x in ds])
    b, a = np.polyfit(ds, y, 1)
    return next(D for D in range(3, 401, 2) if math.exp(a + b * D) <= TARGET)


@pytest.mark.parametrize("p", [
    "0.005",
    pytest.param("0.01", marks=pytest.mark.xfail(strict=True, reason=(
        "Finding, not a bug in the check: the intermediate global-fit entry (d=25) is too small. "
        "It is superseded by final_overhead (d=27), tested below. strict: if this ever passes, "
        "the data changed and the final rule must be re-examined."))),
])
def test_global_fit_entries_compatible_with_stratified(summary, fw, p):
    """The global-fit Part D d (main fit) must not be smaller than the d from the stratified
    extrapolation. Fails if the project's own low-p data say the global-fit entry is too small."""
    d_fw = _d_from_stratified_extrapolation(fw, p)
    d_fit = summary["code_capacity"]["overhead"]["main"]["targets"][p]["d"]
    assert d_fit >= d_fw, f"p={p}: global fit gives d={d_fit}, stratified extrapolation gives d={d_fw}"


@pytest.mark.parametrize("p", ["0.005", "0.01"])
def test_final_entries_compatible_with_stratified(summary, fw, p):
    """Added in the third review: the consolidated answer summary.json["final_overhead"] for code
    capacity must not be smaller than the stratified-extrapolation d, and must be self-consistent."""
    if "final_overhead" not in summary:
        pytest.skip("no final_overhead in summary.json")
    f = summary["final_overhead"]["code_capacity"][p]
    d_fw = _d_from_stratified_extrapolation(fw, p)
    assert f["d"] >= d_fw, f"p={p}: final d={f['d']}, stratified extrapolation gives d={d_fw}"
    assert f["full"] == 2 * f["d"] ** 2 - 1


def test_final_overhead_consistent(summary):
    """final_overhead: one odd d per (model, p), d >= every global-fit candidate, monotone in p,
    phenomenological >= code capacity, and the converged fixed-weight d is used at 0.1 %."""
    if "final_overhead" not in summary:
        pytest.skip("no final_overhead in summary.json")
    fo = summary["final_overhead"]
    for model, rows in fo.items():
        ds = []
        for p in sorted(rows, key=float):
            f = rows[p]
            assert f["d"] is not None and f["d"] % 2 == 1 and f["d"] >= 3
            assert f["full"] == 2 * f["d"] ** 2 - 1
            assert f["z_only"] == f["d"] ** 2 + (f["d"] ** 2 - 1) // 2
            fits = [v for k, v in f["candidates"].items() if k.startswith("global fit") and v]
            assert f["d"] >= max(fits)
            ds.append(f["d"])
        assert ds == sorted(ds), f"{model}: d not monotone in p: {ds}"
    for p in fo["code_capacity"]:
        assert fo["phenomenological"][p]["d"] >= fo["code_capacity"][p]["d"]
    fwreq = summary["code_capacity"]["overhead"]["fixed_weight"]["required"]["0.001"]["d_upper95"]
    assert fo["code_capacity"]["0.001"]["d"] == fwreq
