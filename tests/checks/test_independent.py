"""Independent correctness checks of the simulation pipeline.

1. Independent MWPM (networkx blossom on the detector graph, BFS distances, one boundary copy
   per fired detector): on random errors, PyMatching's correction must (a) reproduce the
   syndrome and (b) have exactly the minimum total weight. Code capacity d = 3, 5, 7 and
   phenomenological d = 3, 5 (q = p, so all edges have equal weight). Does not need results/.
2. Exact d = 3 code capacity: enumerate all 2^9 X errors through the project's decoder ->
   exact p_L(p) polynomial. Every d = 3 code-capacity point in results/raw.csv must agree
   within 4 sigma (binomial). Also: exact ML (optimal coset) p_L is a lower bound that MWPM
   must respect, and the MWPM threshold cannot exceed ML performance.
3. d = 3 phenomenological: enumerate all space-time faults of weight <= 4 (39 fault locations)
   -> lower bound L on p_L; upper bound U = L + P(weight >= 5). Each raw d = 3 phenomenological
   point at p <= 0.01 must satisfy n L - 4 sigma <= fails <= n U + 4 sigma.
"""
import itertools
import math

import networkx as nx
import numpy as np
import pymatching
import pytest
from scipy.stats import binom

from conftest import grid
from surfcode.code import hz, logical_z, spacetime_matrices
from surfcode.sim import code_capacity_decoder, phenomenological_decoder


# ---------------------------------------------------------------- 1. independent MWPM
def _min_weight_independent(H, syn):
    """Minimum number of columns of H whose sum is syn, assuming every column has weight 1 or 2
    (graph-like), via blossom on the complete graph of fired detectors + boundary copies."""
    m = H.shape[0]
    G = nx.Graph()
    G.add_nodes_from(range(m + 1))
    B = m
    for c in range(H.shape[1]):
        rows = np.flatnonzero(H[:, c])
        assert len(rows) in (1, 2)
        G.add_edge(rows[0], rows[1] if len(rows) == 2 else B)
    fired = list(np.flatnonzero(syn))
    if not fired:
        return 0
    dist = {u: nx.single_source_shortest_path_length(G, u) for u in fired}
    K = nx.Graph()
    big = 10 * (m + 1)
    for i, u in enumerate(fired):
        K.add_edge(("d", u), ("b", u), weight=big - dist[u][B])
        for v in fired[i + 1:]:
            K.add_edge(("d", u), ("d", v), weight=big - dist[u][v])
            K.add_edge(("b", u), ("b", v), weight=big)      # boundary copies pair for free
    M = nx.max_weight_matching(K, maxcardinality=True)
    assert 2 * len(M) == K.number_of_nodes()
    return sum(big - K[a][b]["weight"] for a, b in M)


def _compare(H, p, n_trials, seed):
    H = np.asarray(H, dtype=np.uint8)
    mt = pymatching.Matching.from_check_matrix(H)          # uniform weights, returns correction
    rng = np.random.default_rng(seed)
    checked = 0
    for _ in range(n_trials):
        e = (rng.random(H.shape[1]) < p).astype(np.uint8)
        syn = (H.astype(int) @ e) % 2
        corr = mt.decode(syn)
        assert np.array_equal((H.astype(int) @ corr) % 2, syn), "correction does not reproduce syndrome"
        w_ind = _min_weight_independent(H, syn)
        assert int(corr.sum()) == w_ind, f"PyMatching weight {corr.sum()} != independent min {w_ind}"
        assert w_ind <= e.sum()
        checked += syn.any()
    assert checked > n_trials // 3


@pytest.mark.parametrize("d", [3, 5, 7])
def test_independent_mwpm_code_capacity(d):
    _compare(hz(d), 0.08, 150, seed=100 + d)


@pytest.mark.parametrize("d", [3, 5])
def test_independent_mwpm_phenomenological(d):
    H, *_ = spacetime_matrices(d, d)
    _compare(H.toarray(), 0.02, 100, seed=200 + d)


def test_phenom_decoder_weights_uniform_when_q_eq_p():
    """The decoder used in the run must weight data and measurement edges equally when q = p."""
    m, H, L, nD, nM = phenomenological_decoder(3, 0.02, 0.02)
    ws = {round(e[2]["weight"], 9) for e in m.edges()}
    assert len(ws) == 1 and abs(ws.pop() - math.log(0.98 / 0.02)) < 1e-6


# ---------------------------------------------------------------- 2. exact d = 3 code capacity
def _d3_tables():
    m, H, L = code_capacity_decoder(3, 0.05)
    E = np.array(list(itertools.product([0, 1], repeat=9)), dtype=np.uint8)
    syn = (E.astype(int) @ H.T) % 2
    act = (E.astype(int) @ L.T)[:, 0] % 2
    fail = m.decode_batch(syn.astype(np.uint8))[:, 0] != act
    w = E.sum(1)
    skey = syn @ (1 << np.arange(syn.shape[1]))
    return fail, w, skey, act


def exact_mwpm_d3(p):
    fail, w, _, _ = _d3_tables()
    return float(np.sum(fail * p**w * (1 - p) ** (9 - w)))


def exact_ml_d3(p):
    _, w, skey, act = _d3_tables()
    pr = p**w * (1 - p) ** (9 - w)
    tot = 0.0
    for s in np.unique(skey):
        sel = skey == s
        tot += min(pr[sel & (act == 0)].sum(), pr[sel & (act == 1)].sum())
    return tot


def test_d3_exact_sanity():
    fail, w, _, _ = _d3_tables()
    assert not fail[w <= 1].any()               # corrects all single errors
    assert fail[w == 2].sum() > 0               # (d+1)/2 = 2 errors can fail
    for p in (0.01, 0.05, 0.1):
        assert exact_ml_d3(p) <= exact_mwpm_d3(p) + 1e-15


def test_d3_code_capacity_matches_exact(raw):
    g = grid(raw, "code_capacity")
    pts = [(p, r) for (d, p), r in g.items() if d == 3]
    assert pts
    bad = []
    for p, r in pts:
        P = exact_mwpm_d3(p)
        n, k = r["shots"], r["fails"]
        z = (k - n * P) / math.sqrt(n * P * (1 - P))
        if abs(z) > 4:
            bad.append((p, k / n, P, round(z, 2)))
    assert not bad, f"d=3 code capacity disagrees with exact enumeration: {bad}"


# ---------------------------------------------------------------- 3. d = 3 phenomenological
def _phenom_d3_counts(max_w=4):
    m, H, L, nD, nM = phenomenological_decoder(3, 0.02, 0.02)
    H = H.toarray().astype(int); L = L.toarray().astype(int)
    N = nD + nM
    counts = np.zeros(max_w + 1, dtype=np.int64)
    for w in range(1, max_w + 1):
        F = np.zeros((math.comb(N, w), N), dtype=np.uint8)
        for i, supp in enumerate(itertools.combinations(range(N), w)):
            F[i, list(supp)] = 1
        syn = ((F.astype(int) @ H.T) % 2).astype(np.uint8)
        act = (F.astype(int) @ L.T)[:, 0] % 2
        counts[w] = int((m.decode_batch(syn)[:, 0] != act).sum())
    return counts, N


def test_d3_phenomenological_matches_enumeration(raw):
    g = grid(raw, "phenomenological")
    pts = [(p, r) for (d, p), r in g.items() if d == 3 and p <= 0.01]
    assert pts
    counts, N = _phenom_d3_counts(4)
    assert counts[1] == 0 and counts[2] > 0
    bad = []
    for p, r in pts:
        lo = sum(counts[w] * p**w * (1 - p) ** (N - w) for w in range(len(counts)))
        hi = lo + binom.sf(len(counts) - 1, N, p)
        n, k = r["shots"], r["fails"]
        s = math.sqrt(n * hi * (1 - hi))
        if not (n * lo - 4 * s <= k <= n * hi + 4 * s):
            bad.append((p, k / n, lo, hi))
    assert not bad, f"d=3 phenomenological disagrees with low-weight enumeration: {bad}"
