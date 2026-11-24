"""Checks on results/raw.csv: completeness, sample size, internal consistency, monotonicity, ordering.

Tolerances (fixed before the run):
- every point >= 1e5 shots (PLAN floor) and <= 2e6 + 1e5 (top-up cap + one chunk);
- a point with < 200 failures must have been topped up to the 2e6 cap (run.py defaults);
- monotonicity in p and ordering in d are tested as "no reversal beyond 3 sigma"
  (combined binomial std); with ~180 comparisons a 3-sigma false alarm is ~0.25 expected,
  so one marginal reversal would be worth a look, not proof of a bug -- we fail at 3 sigma anyway.
- ordering windows use the literature thresholds (not the fitted one, to avoid circularity):
  code capacity: below p <= 0.08 (~0.78 p_th), above p >= 0.12 (~1.17 p_th);
  phenomenological: below p <= 0.024 (~0.82 p_th), above p >= 0.034 (~1.16 p_th).
"""
import math

import numpy as np
import pytest

from conftest import DS, MIN_SHOTS, MODELS, binom_sigma, grid

# Grids as agreed (PLAN.md 'Numerics'; values in run.py at commit c6a9494). A cut grid fails here.
P_CC = [0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06,
        0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.10, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13]
P_PH = [0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02,
        0.022, 0.024, 0.026, 0.027, 0.028, 0.029, 0.030, 0.031, 0.032, 0.034, 0.036]
PGRID = {"code_capacity": P_CC, "phenomenological": P_PH}
BELOW = {"code_capacity": 0.08, "phenomenological": 0.024}
ABOVE = {"code_capacity": 0.12, "phenomenological": 0.034}
MAX_SHOTS = 2_000_000
TARGET_FAILS = 200


@pytest.mark.parametrize("model", MODELS)
def test_grid_complete(raw, model):
    g = grid(raw, model)
    want = {(d, p) for d in DS for p in PGRID[model]}
    have = {(d, round(p, 6)) for d, p in g}
    missing = {(d, round(p, 6)) for d, p in want} - have
    assert not missing, f"{model}: missing (d,p) points {sorted(missing)}"
    assert len(g) == len(want), f"{model}: duplicate or extra points ({len(g)} vs {len(want)})"


def test_no_nan_and_consistent(raw):
    for r in raw:
        for k in ("p_L", "ci_lo", "ci_hi", "seconds"):
            assert math.isfinite(r[k]), f"non-finite {k} in {r}"
        assert 0 <= r["fails"] <= r["shots"]
        assert r["p_L"] == pytest.approx(r["fails"] / r["shots"], rel=1e-9, abs=1e-15)
        assert r["ci_lo"] - 1e-12 <= r["p_L"] <= r["ci_hi"] + 1e-12, r
        assert 0 <= r["ci_lo"] and r["ci_hi"] <= 1


def test_min_shots(raw):
    bad = [(r["model"], r["d"], r["p"], r["shots"]) for r in raw if r["shots"] < MIN_SHOTS]
    assert not bad, f"points below the agreed 1e5-shot floor: {bad}"


def test_topup_rule(raw):
    """Points with < 200 failures must sit at the 2e6 cap; no point exceeds cap + one chunk."""
    bad = []
    for r in raw:
        if r["shots"] > MAX_SHOTS + MIN_SHOTS:
            bad.append(("over cap", r["model"], r["d"], r["p"], r["shots"]))
        if r["fails"] < TARGET_FAILS and r["shots"] < MAX_SHOTS:
            bad.append(("not topped up", r["model"], r["d"], r["p"], r["shots"], r["fails"]))
    assert not bad, bad


def test_pL_below_half(raw):
    """A decoder never does worse than random guessing on average: p_L <= 0.5 + 3 sigma."""
    for r in raw:
        assert r["p_L"] <= 0.5 + 3 * binom_sigma(r["fails"], r["shots"]), r


@pytest.mark.parametrize("model", MODELS)
def test_monotone_in_p(raw, model):
    g = grid(raw, model)
    bad = []
    for d in DS:
        ps = sorted(p for dd, p in g if dd == d)
        for p1, p2 in zip(ps, ps[1:]):
            a, b = g[(d, p1)], g[(d, p2)]
            s = math.hypot(binom_sigma(a["fails"], a["shots"]), binom_sigma(b["fails"], b["shots"]))
            z = (b["p_L"] - a["p_L"]) / s
            if z < -3:
                bad.append((d, p1, p2, a["p_L"], b["p_L"], round(z, 2)))
    assert not bad, f"{model}: p_L decreases with p beyond 3 sigma: {bad}"


@pytest.mark.parametrize("model", MODELS)
def test_ordering_in_d(raw, model):
    """Below threshold larger d is better; above threshold larger d is worse.
    Fails on a reversal > 3 sigma; also requires the correct ordering to be *significant*
    (> 3 sigma) for every adjacent pair at the lowest p where both have >= 20 failures and at
    the highest p of the grid."""
    g = grid(raw, model)
    ps = sorted({p for _, p in g})
    bad = []
    for p in ps:
        if BELOW[model] < p < ABOVE[model]:
            continue
        sign = -1 if p <= BELOW[model] else +1        # expected sign of p_L(d+2) - p_L(d)
        for d1, d2 in zip(DS, DS[1:]):
            a, b = g[(d1, p)], g[(d2, p)]
            if a["fails"] == 0 and b["fails"] == 0:
                continue
            s = math.hypot(binom_sigma(a["fails"], a["shots"]), binom_sigma(b["fails"], b["shots"]))
            z = sign * (b["p_L"] - a["p_L"]) / s
            if z < -3:
                bad.append(("reversed", p, d1, d2, a["p_L"], b["p_L"], round(z, 2)))
    # significance at the extremes
    pmax = max(ps)
    for d1, d2 in zip(DS, DS[1:]):
        a, b = g[(d1, pmax)], g[(d2, pmax)]
        s = math.hypot(binom_sigma(a["fails"], a["shots"]), binom_sigma(b["fails"], b["shots"]))
        if (b["p_L"] - a["p_L"]) / s < 3:
            bad.append(("above-threshold ordering not significant", pmax, d1, d2, a["p_L"], b["p_L"]))
        low = [p for p in ps if p <= BELOW[model] and g[(d1, p)]["fails"] >= 20 and g[(d2, p)]["fails"] >= 20]
        if not low:
            bad.append(("no below-threshold p with >=20 failures for both", d1, d2))
            continue
        p = min(low)
        a, b = g[(d1, p)], g[(d2, p)]
        s = math.hypot(binom_sigma(a["fails"], a["shots"]), binom_sigma(b["fails"], b["shots"]))
        if (a["p_L"] - b["p_L"]) / s < 3:
            bad.append(("below-threshold ordering not significant", p, d1, d2, a["p_L"], b["p_L"]))
    assert not bad, f"{model}: {bad}"


@pytest.mark.parametrize("model", MODELS)
def test_low_p_points_have_data(raw, model):
    """Every d must have >= 3 below-threshold points with >= 20 failures, else the scaling fits
    (PLAN part B) for that d rest on too little."""
    g = grid(raw, model)
    for d in DS:
        n = sum(1 for (dd, p), r in g.items() if dd == d and p <= BELOW[model] and r["fails"] >= 20)
        assert n >= 3, f"{model} d={d}: only {n} fit-usable below-threshold points"
