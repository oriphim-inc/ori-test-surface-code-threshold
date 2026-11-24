"""Overhead checks (PLAN part D): p_L = A (p/p_*)^{(d+1)/2}, smallest odd d with p_L <= 1e-12.

Fixed before the run:
- every target has an odd d >= 3, qubits full = 2d^2 - 1, z_only = d^2 + (d^2-1)/2;
- p_L_pred recomputed here from (A, p_*) equals the stored value and is <= 1e-12, while d-2
  would give > 1e-12 (d is the *smallest* such odd distance);
- d is non-decreasing in p (0.1 % <= 0.5 % <= 1 %) and phenomenological d >= code-capacity d;
- the fit window must lie below threshold (p_max < 0.6 p_th main, < 0.8 p_th alt) and p_* must be within a factor 3
  of the fitted threshold (p_* plays the role of an effective threshold; a value far off means
  the fixed-exponent form does not describe the data);
- robustness of the extrapolation: main vs alt fit windows must give d within 4 of each other
  (two odd steps) at every target p; larger spread means the 1e-12 answer is fit-window driven.
- global fit reduced chi2 <= 5 (generous: the fixed (d+1)/2 exponent ignores entropic
  corrections, so some misfit is expected; beyond 5 the form is not describing the data).
"""
import pytest

from conftest import MODELS

TARGET = 1e-12
P_TARGETS = ["0.001", "0.005", "0.01"]


def pred(g, p, d):
    return g["A"] * (p / g["p_star"]) ** ((d + 1) / 2)


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.parametrize("tag", ["main", "alt"])
def test_targets_consistent(summary, model, tag):
    g = summary[model]["overhead"][tag]
    assert sorted(g["targets"], key=float) == P_TARGETS
    prev = 0
    for ps in P_TARGETS:
        t = g["targets"][ps]
        p, d = float(ps), t["d"]
        assert d is not None, f"{model}/{tag} p={ps}: no distance found"
        assert d >= 3 and d % 2 == 1, f"{model}/{tag} p={ps}: d = {d}"
        assert t["full"] == 2 * d * d - 1
        assert t["z_only"] == d * d + (d * d - 1) // 2
        pl = pred(g, p, d)
        assert pl == pytest.approx(t["p_L_pred"], rel=1e-9)
        assert pl <= TARGET
        if d > 3:
            assert pred(g, p, d - 2) > TARGET, f"{model}/{tag} p={ps}: d={d} is not the smallest"
        assert d >= prev, f"{model}/{tag}: d decreases with p"
        prev = d


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.parametrize("tag", ["main", "alt"])
def test_fit_sane(summary, model, tag):
    g = summary[model]["overhead"][tag]
    pth = summary[model]["threshold"]["fss"]["p_th"]
    limit = 0.6 * pth if tag == "main" else 0.8 * pth
    assert g["p_max"] < limit, f"{model}/{tag}: fit window p_max={g['p_max']} too close to p_th={pth:.4f}"
    assert pth / 3 <= g["p_star"] <= 3 * pth, f"{model}/{tag}: p_* = {g['p_star']:.4g} vs p_th = {pth:.4g}"
    assert g["n_points"] >= 6 and g["dof"] >= 4
    assert g["chi2"] / g["dof"] <= 5, f"{model}/{tag}: global fit chi2/dof = {g['chi2']:.1f}/{g['dof']}"


@pytest.mark.parametrize("model", MODELS)
def test_main_alt_robust(summary, model):
    o = summary[model]["overhead"]
    for ps in P_TARGETS:
        dm, da = o["main"]["targets"][ps]["d"], o["alt"]["targets"][ps]["d"]
        assert abs(dm - da) <= 4, f"{model} p={ps}: d main={dm} vs alt={da} -- extrapolation not robust"


def test_phenom_needs_more(summary):
    for tag in ("main", "alt"):
        for ps in P_TARGETS:
            dc = summary["code_capacity"]["overhead"][tag]["targets"][ps]["d"]
            dp = summary["phenomenological"]["overhead"][tag]["targets"][ps]["d"]
            assert dp >= dc, f"{tag} p={ps}: phenomenological d={dp} < code-capacity d={dc}"
