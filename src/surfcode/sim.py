"""Monte-Carlo sampling and MWPM decoding."""
import numpy as np
import pymatching

from .code import hz, logical_z, spacetime_matrices


def _weight(p):
    return float(np.log((1 - p) / p))


def code_capacity_decoder(d, p):
    H = hz(d)
    L = logical_z(d)
    m = pymatching.Matching.from_check_matrix(H, weights=np.full(H.shape[1], _weight(p)),
                                              faults_matrix=L)
    return m, H, L


def sample_code_capacity(d, p, shots, rng, batch=200_000):
    """Return number of logical failures in `shots` trials (one perfect syndrome round)."""
    m, H, L = code_capacity_decoder(d, p)
    n = H.shape[1]
    fails = 0
    done = 0
    while done < shots:
        b = min(batch, shots - done)
        err = (rng.random((b, n)) < p).astype(np.uint8)
        syn = (err @ H.T) % 2
        actual = (err @ L.T) % 2
        pred = m.decode_batch(syn.astype(np.uint8))
        fails += int(np.any(pred != actual, axis=1).sum())
        done += b
    return fails


def phenomenological_decoder(d, p, q, rounds=None):
    T = d if rounds is None else rounds
    H, L, nD, nM = spacetime_matrices(d, T)
    w = np.concatenate([np.full(nD, _weight(p)), np.full(nM, _weight(q))])
    m = pymatching.Matching.from_check_matrix(H, weights=w, faults_matrix=L)
    return m, H, L, nD, nM


def sample_phenomenological(d, p, shots, rng, q=None, rounds=None, batch=None):
    """Data errors p and measurement errors q (=p) in each of `rounds` (=d) noisy rounds,
    followed by a perfect round. Faults sampled independently per space-time column."""
    q = p if q is None else q
    m, H, L, nD, nM = phenomenological_decoder(d, p, q, rounds)
    Hcsr = H.tocsr().T.tocsr()  # (cols x dets) for err @ H.T
    Lt = L.T.tocsr()
    probs = np.concatenate([np.full(nD, p), np.full(nM, q)])
    if batch is None:  # keep each random block ~5e6 entries (memory)
        batch = max(1000, int(5e6 // (nD + nM)))
    fails = 0
    done = 0
    while done < shots:
        b = min(batch, shots - done)
        err = (rng.random((b, nD + nM)) < probs).astype(np.uint8)
        syn = np.asarray((sp_mul(err, Hcsr)) % 2, dtype=np.uint8)
        actual = np.asarray(sp_mul(err, Lt) % 2, dtype=np.uint8)
        pred = m.decode_batch(syn)
        fails += int(np.any(pred != actual, axis=1).sum())
        done += b
    return fails


def sp_mul(dense, sparse):
    """dense (b x k) uint8 @ sparse (k x j) -> dense int array."""
    return np.asarray((sparse.T @ dense.T.astype(np.int32)).T)
