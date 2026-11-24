"""Below-threshold scaling checks (PLAN part B).

beta = fitted slope of ln p_L vs ln p at fixed d over p < 0.8 p_th (points with >= 20 fails).
Asymptotically (p -> 0) beta -> (d+1)/2 for odd d (smallest uncorrectable X error has weight
(d+1)/2, in space-time too for the phenomenological model). Over a finite p window the local
slope is *lower*: for d = 3 code capacity the exact MWPM p_L (enumeration of all 2^9 errors,
see test_independent.py) gives a fitted slope of 1.89 over p in [0.01, 0.06], i.e. ~5 % below 2.
For large d only high-p points have >= 20 failures (2e6 shot cap), so beta is expected to fall
well short of (d+1)/2 there; that is a limit of the data, not a failure. Tolerances:
- d = 3: beta in [1.70, 2.15] (+- 3 sigma); d = 5: beta in [2.40, 3.25] (+- 3 sigma);
- every d: beta <= (d+1)/2 + 0.3 + 3 sigma (a slope clearly steeper than the minimum-weight
  exponent would mean mis-counted failures or wrong p);
- beta increases with d (beta(d+2) > beta(d) - 3 sigma_combined), and beta(11) > beta(3) + 1.
Lambda (p_L ~ Lambda^{-(d+1)/2} at fixed p): Lambda > 1 at every fitted p below threshold
(and > 1 by 3 sigma at at least half of them), and non-increasing in p (Lambda(p_hi) <= Lambda(p_lo) + 3 sigma).
"""
import math

import pytest

from conftest import MODELS

BETA_BAND = {3: (1.70, 2.15), 5: (2.40, 3.25)}


def vs_p(summary, model):
    return {int(k): v for k, v in summary[model]["scaling"]["vs_p"].items()}


@pytest.mark.parametrize("model", MODELS)
def test_beta_present_for_all_d(summary, model):
    assert sorted(vs_p(summary, model)) == [3, 5, 7, 9, 11], f"{model}: beta fit missing for some d"


@pytest.mark.parametrize("model", MODELS)
def test_beta_small_d_matches_min_weight(summary, model):
    v = vs_p(summary, model)
    for d, (lo, hi) in BETA_BAND.items():
        b, e = v[d]["beta"], v[d]["beta_err"]
        assert lo - 3 * e <= b <= hi + 3 * e, f"{model} d={d}: beta = {b:.2f} +- {e:.2f}, band [{lo}, {hi}]"
        assert v[d]["expected"] == (d + 1) / 2


@pytest.mark.parametrize("model", MODELS)
def test_beta_not_steeper_than_min_weight(summary, model):
    for d, r in vs_p(summary, model).items():
        assert r["beta"] <= (d + 1) / 2 + 0.3 + 3 * r["beta_err"], (
            f"{model} d={d}: beta = {r['beta']:.2f} steeper than (d+1)/2 = {(d+1)/2}")


@pytest.mark.parametrize("model", MODELS)
def test_beta_increases_with_d(summary, model):
    v = vs_p(summary, model)
    ds = sorted(v)
    for d1, d2 in zip(ds, ds[1:]):
        s = math.hypot(v[d1]["beta_err"], v[d2]["beta_err"])
        assert v[d2]["beta"] > v[d1]["beta"] - 3 * s, f"{model}: beta({d2}) < beta({d1})"
    assert v[ds[-1]]["beta"] > v[ds[0]]["beta"] + 1, f"{model}: beta barely grows with d"


@pytest.mark.parametrize("model", MODELS)
def test_lambda_above_one_and_decreasing(summary, model):
    sc = summary[model]["scaling"]
    vd = sorted(((float(p), r) for p, r in sc["vs_d"].items()), key=lambda t: t[0])
    assert len(vd) >= 3, f"{model}: Lambda fitted at only {len(vd)} p values"
    sig = 0
    for p, r in vd:
        assert p < sc["p_cut"]
        assert r["n_points"] >= 2
        assert r["Lambda"] > 1, f"{model} p={p}: Lambda = {r['Lambda']:.2f} <= 1 below threshold"
        sig += r["Lambda"] - 3 * r["Lambda_err"] > 1
    assert sig >= len(vd) / 2, f"{model}: Lambda > 1 significant (3 sigma) at only {sig}/{len(vd)} p values"
    for (p1, r1), (p2, r2) in zip(vd, vd[1:]):
        s = math.hypot(r1["Lambda_err"], r2["Lambda_err"])
        assert r2["Lambda"] <= r1["Lambda"] + 3 * s, f"{model}: Lambda rises from p={p1} to p={p2}"
