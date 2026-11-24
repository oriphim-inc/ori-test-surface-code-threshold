"""Fixed-weight (stratified) estimate of the code-capacity logical error rate at low p.

p_L(p) = sum_w P(W = w) f_w, W ~ Binomial(n, p), f_w = MWPM failure fraction over uniformly random
weight-w X errors. f_w = 0 for w <= (d-1)/2 (MWPM corrects up to t errors). Truncating at w_max leaves a
remainder in [0, P(W > w_max)], so the estimate comes with a rigorous truncation bound.
"""
import numpy as np
from scipy.stats import binom

from .sim import code_capacity_decoder


def failure_fraction(d, w, samples, rng, batch=200_000):
    m, H, L = code_capacity_decoder(d, 0.01)  # uniform weights: decoder independent of p
    n = H.shape[1]
    fails = 0; done = 0
    while done < samples:
        b = min(batch, samples - done)
        idx = np.argsort(rng.random((b, n)), axis=1)[:, :w]
        err = np.zeros((b, n), np.uint8)
        np.put_along_axis(err, idx, 1, axis=1)
        syn = (err @ H.T) % 2
        actual = (err @ L.T) % 2
        pred = m.decode_batch(syn.astype(np.uint8))
        fails += int(np.any(pred != actual, axis=1).sum())
        done += b
    return fails, samples


def stratified_pl(d, p, fw):
    """fw: dict w -> (fails, samples). Returns estimate, 95% stat interval, truncation bound."""
    n = d * d
    est = 0.0; var = 0.0
    for w, (k, N) in fw.items():
        pw = binom.pmf(w, n, p)
        f = k / N
        est += pw * f
        var += pw**2 * max(f * (1 - f), 1.0 / N**2) / N
    wmax = max(fw)
    trunc = float(binom.sf(wmax, n, p))
    sd = np.sqrt(var)
    return dict(p_L=float(est), lo=float(est - 1.96 * sd), hi=float(est + 1.96 * sd + trunc), trunc=trunc)
