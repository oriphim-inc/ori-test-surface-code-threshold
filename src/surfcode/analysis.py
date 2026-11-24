"""Statistics and fits: Wilson intervals, finite-size threshold fit, scaling fits, overhead."""
import numpy as np
from scipy.optimize import curve_fit, brentq


def wilson(k, n, z=1.96):
    k = np.asarray(k, float); n = np.asarray(n, float)
    ph = k / n
    den = 1 + z**2 / n
    c = (ph + z**2 / (2 * n)) / den
    h = z * np.sqrt(ph * (1 - ph) / n + z**2 / (4 * n**2)) / den
    return c - h, c + h


def _fss(X, pth, nu, a, b, c):
    p, d = X
    x = (p - pth) * d ** (1 / nu)
    return a + b * x + c * x**2


def fss_fit(p, d, pl, n, p0):
    """Fit p_L = a + b x + c x^2, x = (p - p_th) d^(1/nu). Returns params, cov."""
    sig = np.sqrt(np.maximum(pl * (1 - pl), 1e-12) / n)
    popt, pcov = curve_fit(_fss, (p, d), pl, p0=[p0, 1.0, 0.1, 1.0, 0.0],
                           sigma=sig, absolute_sigma=True, maxfev=20000)
    return popt, pcov


def fss_threshold(p, d, k, n, p0, n_boot=300, seed=1):
    """Threshold from FSS fit with parametric bootstrap (binomial resampling)."""
    p = np.asarray(p, float); d = np.asarray(d, float)
    k = np.asarray(k, float); n = np.asarray(n, float)
    pl = k / n
    popt, pcov = fss_fit(p, d, pl, n, p0)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        kb = rng.binomial(n.astype(np.int64), pl)
        try:
            pb, _ = fss_fit(p, d, kb / n, n, popt[0])
            boots.append(pb[0])
        except RuntimeError:
            pass
    boots = np.array(boots)
    resid = (pl - _fss((p, d), *popt)) / np.sqrt(np.maximum(pl * (1 - pl), 1e-12) / n)
    chi2 = float(np.sum(resid**2)); dof = len(p) - 5
    return dict(p_th=float(popt[0]), p_th_err=float(np.std(boots)), nu=float(popt[1]),
                params=popt.tolist(), chi2=chi2, dof=dof, n_boot_ok=int(len(boots)))


def pairwise_crossing(pgrid, pl1, pl2):
    """Crossing of two curves by linear interpolation of log p_L difference; None if no sign change."""
    diff = np.log(pl1) - np.log(pl2)
    for i in range(len(pgrid) - 1):
        if np.isfinite(diff[i]) and np.isfinite(diff[i + 1]) and diff[i] * diff[i + 1] < 0:
            t = diff[i] / (diff[i] - diff[i + 1])
            return float(pgrid[i] + t * (pgrid[i + 1] - pgrid[i]))
    return None


def fit_vs_d(ds, k, n, min_fail=20):
    """Weighted fit ln p_L = a - alpha d. Returns alpha, err, Lambda = exp(2 alpha)."""
    ds = np.asarray(ds, float); k = np.asarray(k, float); n = np.asarray(n, float)
    ok = k >= min_fail
    if ok.sum() < 2:
        return None
    y = np.log(k[ok] / n[ok]); s = 1 / np.sqrt(k[ok])  # rel. error of p_L ~ 1/sqrt(k)
    A = np.vstack([np.ones(ok.sum()), -ds[ok]]).T
    W = np.diag(1 / s**2)
    cov = np.linalg.inv(A.T @ W @ A)
    a, alpha = cov @ A.T @ W @ y
    return dict(alpha=float(alpha), alpha_err=float(np.sqrt(cov[1, 1])),
                Lambda=float(np.exp(2 * alpha)), Lambda_err=float(2 * np.exp(2 * alpha) * np.sqrt(cov[1, 1])),
                a=float(a), n_points=int(ok.sum()), ds=ds[ok].tolist())


def fit_vs_p(ps, k, n, min_fail=20):
    """Weighted fit ln p_L = c + beta ln p at fixed d."""
    ps = np.asarray(ps, float); k = np.asarray(k, float); n = np.asarray(n, float)
    ok = k >= min_fail
    if ok.sum() < 2:
        return None
    y = np.log(k[ok] / n[ok]); s = 1 / np.sqrt(k[ok])
    A = np.vstack([np.ones(ok.sum()), np.log(ps[ok])]).T
    W = np.diag(1 / s**2)
    cov = np.linalg.inv(A.T @ W @ A)
    c, beta = cov @ A.T @ W @ y
    return dict(beta=float(beta), beta_err=float(np.sqrt(cov[1, 1])), n_points=int(ok.sum()),
                p_range=[float(ps[ok].min()), float(ps[ok].max())])


def global_fit(p, d, k, n, min_fail=20, p_max=None):
    """Fit ln p_L = ln A + (d+1)/2 * (ln p - ln p_star). Returns A, p_star and covariance."""
    p = np.asarray(p, float); d = np.asarray(d, float); k = np.asarray(k, float); n = np.asarray(n, float)
    ok = k >= min_fail
    if p_max is not None:
        ok &= p <= p_max
    e = (d[ok] + 1) / 2
    y = np.log(k[ok] / n[ok]) - e * np.log(p[ok])
    # y = ln A - e ln p_star  -> linear in (lnA, ln p_star)
    A = np.vstack([np.ones(ok.sum()), -e]).T
    W = np.diag(k[ok])
    cov = np.linalg.inv(A.T @ W @ A)
    lnA, lnps = cov @ A.T @ W @ y
    resid = (y - A @ np.array([lnA, lnps])) * np.sqrt(k[ok])
    return dict(A=float(np.exp(lnA)), p_star=float(np.exp(lnps)), cov=cov.tolist(),
                chi2=float(np.sum(resid**2)), dof=int(ok.sum() - 2), n_points=int(ok.sum()))


def predict(fit, p, d):
    return fit["A"] * (p / fit["p_star"]) ** ((d + 1) / 2)


def required_distance(fit, p, target=1e-12, d_max=401):
    for d in range(3, d_max + 1, 2):
        if predict(fit, p, d) <= target:
            return d
    return None


def qubits(d):
    return dict(full=2 * d * d - 1, z_only=d * d + (d * d - 1) // 2)
