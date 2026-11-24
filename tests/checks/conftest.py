"""Sceptic's acceptance checks (written before the main run, 2026-11-13).

Fixtures load results/raw.csv and results/summary.json; every test that needs them skips
cleanly if they are absent. Tolerances are fixed here, before any result exists, and must not
be changed to make a run pass (only the Sceptic edits tests/checks/).
"""
import csv
import json
import math
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS = os.path.join(ROOT, "results")
RAW = os.path.join(RESULTS, "raw.csv")
SUMMARY = os.path.join(RESULTS, "summary.json")

MODELS = ("code_capacity", "phenomenological")
DS = (3, 5, 7, 9, 11)
MIN_SHOTS = 100_000          # PLAN.md: >= 1e5 trials per point, never fewer
LIT = {"code_capacity": 0.103, "phenomenological": 0.0293}   # PLAN.md references (WHP 2003)


@pytest.fixture(scope="session")
def raw():
    if not os.path.exists(RAW):
        pytest.skip("results/raw.csv not present (main run not done)")
    rows = []
    with open(RAW) as fh:
        for r in csv.DictReader(fh):
            rows.append(dict(model=r["model"], d=int(r["d"]), p=float(r["p"]),
                             shots=int(r["shots"]), fails=int(r["fails"]),
                             p_L=float(r["p_L"]), ci_lo=float(r["ci_lo"]),
                             ci_hi=float(r["ci_hi"]), seconds=float(r["seconds"])))
    return rows


@pytest.fixture(scope="session")
def summary():
    if not os.path.exists(SUMMARY):
        pytest.skip("results/summary.json not present (main run not done)")
    with open(SUMMARY) as fh:
        return json.load(fh)


def grid(rows, model):
    """{(d, p): row} for one model."""
    return {(r["d"], r["p"]): r for r in rows if r["model"] == model}


def binom_sigma(k, n):
    """Std error of k/n, floored so that 0-failure points are not infinitely precise."""
    ph = max(k, 1) / n
    return math.sqrt(ph * (1 - ph) / n)
