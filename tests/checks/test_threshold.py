"""Threshold checks (PLAN parts A and C).

Two tiers, fixed before the run:

1. "Not a bug" band (hard). Outside it the pipeline is presumed wrong:
   - code capacity MWPM: [0.095, 0.111]. MWPM on the 2D bit-flip surface code has
     p_th ~ 10.3 % (Wang-Harrington-Preskill 2003; Dennis et al. 2002); published finite-size
     estimates on planar/rotated patches sit within ~10.0-10.5 %. The optimal (ML) threshold is
     ~10.9 %, so anything above ~11 % would be beating ML-level performance with MWPM -> bug.
     Lower edge symmetric (+-0.8 %).
   - phenomenological (q = p): [0.026, 0.034]. WHP 2003 MWPM 2.93 %; other MWPM studies give
     ~2.9-3.0 % (optimal 3D RPGM ~3.3 %). d-round blocks with a perfect final round shift
     small-d crossings somewhat; +-0.35 % is allowed before calling it a bug.
2. PLAN agreement (the stated success criterion "within a few tenths of a percent"):
   |p_th - ref| <= 0.003 (code capacity) or <= 0.0015 (phenomenological, ~5 % relative),
   each widened by 2 x bootstrap std. Failing tier 2 but passing tier 1 is a *reported
   disagreement*, not a reason to tune.

Also: the fit must be adequate (chi2/dof <= 3, enough bootstrap successes), precise enough to
say something (p_th_err <= 0.003 / 0.001), nu in a physically sensible range [0.7, 2.5]
(2D RBIM nu ~1.5, 3D RPGM nu ~1.0; this only catches a broken fit), and all adjacent-d
pairwise crossings must exist and agree with the FSS value within 0.01 / 0.004.
"""
import pytest

from conftest import LIT, MODELS

HARD = {"code_capacity": (0.095, 0.111), "phenomenological": (0.026, 0.034)}
AGREE_TOL = {"code_capacity": 0.003, "phenomenological": 0.0015}
MAX_ERR = {"code_capacity": 0.003, "phenomenological": 0.001}
CROSS_TOL = {"code_capacity": 0.010, "phenomenological": 0.004}


def fss(summary, model):
    return summary[model]["threshold"]["fss"]


@pytest.mark.parametrize("model", MODELS)
def test_threshold_not_a_bug(summary, model):
    lo, hi = HARD[model]
    pth = fss(summary, model)["p_th"]
    assert lo <= pth <= hi, f"{model}: p_th = {pth:.4f} outside sanity band [{lo}, {hi}]"


@pytest.mark.parametrize("model", MODELS)
def test_threshold_agrees_with_plan(summary, model):
    f = fss(summary, model)
    tol = AGREE_TOL[model] + 2 * f["p_th_err"]
    assert abs(f["p_th"] - LIT[model]) <= tol, (
        f"{model}: p_th = {f['p_th']:.4f} +- {f['p_th_err']:.4f} vs literature {LIT[model]} "
        f"(allowed {tol:.4f}) -- report as a disagreement, do not tune")


@pytest.mark.parametrize("model", MODELS)
def test_fss_fit_quality(summary, model):
    f = fss(summary, model)
    assert f["dof"] >= 10, f"{model}: too few points in FSS window (dof={f['dof']})"
    assert f["chi2"] / f["dof"] <= 3.0, f"{model}: FSS chi2/dof = {f['chi2']:.1f}/{f['dof']}"
    assert f["n_boot_ok"] >= 250, f"{model}: only {f['n_boot_ok']} bootstrap fits converged"
    assert 0 < f["p_th_err"] <= MAX_ERR[model], f"{model}: p_th_err = {f['p_th_err']}"
    assert 0.7 <= f["nu"] <= 2.5, f"{model}: nu = {f['nu']:.2f}"


@pytest.mark.parametrize("model", MODELS)
def test_fss_threshold_inside_window(summary, model):
    th = summary[model]["threshold"]
    lo, hi = th["window"]
    assert lo < th["fss"]["p_th"] < hi, f"{model}: p_th at/outside fit window {th['window']}"


@pytest.mark.parametrize("model", MODELS)
def test_pairwise_crossings(summary, model):
    th = summary[model]["threshold"]
    pth = th["fss"]["p_th"]
    bad = {}
    for pair, x in th["pairwise_crossings"].items():
        if x is None or abs(x - pth) > CROSS_TOL[model]:
            bad[pair] = x
    assert not bad, f"{model}: crossings missing or far from FSS p_th={pth:.4f}: {bad}"
