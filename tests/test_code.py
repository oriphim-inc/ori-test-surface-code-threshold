import itertools
import numpy as np
import pytest
from surfcode.code import hz, faces_to_matrix, x_faces, logical_x, logical_z, spacetime_matrices
from surfcode.sim import code_capacity_decoder, phenomenological_decoder


def gf2_rank(M):
    M = M.copy() % 2
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = [i for i in range(r, rows) if M[i, c]]
        if not piv:
            continue
        M[[r, piv[0]]] = M[[piv[0], r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
    return r


@pytest.mark.parametrize("d", [3, 5, 7, 9, 11])
def test_code_structure(d):
    Hz = hz(d); Hx = faces_to_matrix(x_faces(d), d * d)
    assert Hz.shape == ((d * d - 1) // 2, d * d)
    assert Hx.shape == ((d * d - 1) // 2, d * d)
    assert not ((Hz.astype(int) @ Hx.T.astype(int)) % 2).any()
    assert gf2_rank(Hz) + gf2_rank(Hx) == d * d - 1   # one logical qubit
    assert not ((Hz.astype(int) @ logical_x(d).T) % 2).any()
    assert not ((Hx.astype(int) @ logical_z(d).T) % 2).any()
    assert int((logical_x(d) @ logical_z(d).T)[0, 0]) % 2 == 1
    # every data qubit is in 1 or 2 Z checks (graph-like for matching)
    assert set(Hz.sum(0)) <= {1, 2}


def test_d3_distance_bruteforce():
    d = 3; Hz = hz(d); Lz = logical_z(d)[0]
    best = 99
    for bits in itertools.product([0, 1], repeat=d * d):
        e = np.array(bits)
        if not ((Hz @ e) % 2).any() and (Lz @ e) % 2:
            best = min(best, e.sum())
    assert best == 3


@pytest.mark.parametrize("d", [3, 5, 7])
def test_mwpm_corrects_up_to_t(d):
    m, Hz, L = code_capacity_decoder(d, 0.05)
    n = d * d; t = (d - 1) // 2
    for w in range(t + 1):
        for supp in itertools.combinations(range(n), w):
            e = np.zeros(n, np.uint8); e[list(supp)] = 1
            pred = m.decode((Hz @ e) % 2)
            assert pred[0] == (L[0] @ e) % 2


@pytest.mark.parametrize("d", [3, 5])
def test_spacetime_matches_roundwise_simulation(d):
    """Explicit round-by-round simulation gives the same detectors as H @ faults."""
    T = d
    H, L, nD, nM = spacetime_matrices(d, T)
    H = H.toarray(); L = L.toarray()
    Hz = hz(d); m, n = Hz.shape
    rng = np.random.default_rng(0)
    for _ in range(50):
        f = (rng.random(nD + nM) < 0.1).astype(np.uint8)
        D = f[:nD].reshape(T, n); M = f[nD:].reshape(T, m)
        cum = np.zeros(n, np.uint8); prev = np.zeros(m, np.uint8); dets = []
        for t in range(T):
            cum ^= D[t]
            s = ((Hz @ cum) % 2).astype(np.uint8) ^ M[t]
            dets.append(s ^ prev); prev = s
        s = ((Hz @ cum) % 2).astype(np.uint8)
        dets.append(s ^ prev)
        assert np.array_equal(np.concatenate(dets), (H @ f) % 2)
        assert (L @ f)[0] % 2 == (logical_z(d)[0] @ cum) % 2


@pytest.mark.parametrize("d", [3, 5])
def test_phenom_corrects_up_to_t(d):
    m, H, L, nD, nM = phenomenological_decoder(d, 0.02, 0.02)
    H = H.toarray(); L = L.toarray(); N = nD + nM; t = (d - 1) // 2
    for w in range(t + 1):
        for supp in itertools.combinations(range(N), w):
            f = np.zeros(N, np.uint8); f[list(supp)] = 1
            assert m.decode((H @ f) % 2)[0] == (L @ f)[0] % 2


def test_fixedweight_matches_exact_d3():
    """Stratified estimator equals exact enumeration of all 2^9 errors at d = 3."""
    from scipy.stats import binom
    from surfcode.fixedweight import failure_fraction, stratified_pl
    d = 3; n = 9
    m, Hz, L = code_capacity_decoder(d, 0.05)
    nfail = np.zeros(n + 1); ntot = np.zeros(n + 1)
    for bits in itertools.product([0, 1], repeat=n):
        e = np.array(bits, np.uint8); w = e.sum()
        ntot[w] += 1
        nfail[w] += m.decode((Hz @ e) % 2)[0] != (L[0] @ e) % 2
    assert nfail[0] == nfail[1] == 0
    rng = np.random.default_rng(3)
    fw = {w: failure_fraction(d, w, 200_000, rng) for w in range(2, 6)}
    for w in range(2, 6):
        f_exact = nfail[w] / ntot[w]
        k, N = fw[w]
        assert abs(k / N - f_exact) < 5 * np.sqrt(f_exact * (1 - f_exact) / N) + 1e-9
    p = 0.01
    exact = sum(binom.pmf(w, n, p) * nfail[w] / ntot[w] for w in range(n + 1))
    r = stratified_pl(d, p, fw)
    assert r["lo"] <= exact <= r["hi"]
