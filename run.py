"""Reproduce every result: python run.py [--shots N] [--max-shots M] [--target-fails K] [--out DIR]

Rotated surface code, X noise, MWPM (PyMatching). Parts A-D of PLAN.md.
"""
import argparse
import csv
import json
import os
import sys
import time
from multiprocessing import Pool

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), ".mplcache"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pymatching
from surfcode.sim import sample_code_capacity, sample_phenomenological
from surfcode import analysis as an
from surfcode.fixedweight import failure_fraction, stratified_pl

DS = [3, 5, 7, 9, 11]
P_CC = [0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06,
        0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.10, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13]
P_PH = [0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02,
        0.022, 0.024, 0.026, 0.027, 0.028, 0.029, 0.030, 0.031, 0.032, 0.034, 0.036]
FSS_WINDOW = {"code_capacity": (0.085, 0.12), "phenomenological": (0.024, 0.034)}
GLOBAL_PMAX = {"code_capacity": 0.05, "phenomenological": 0.015}
GLOBAL_PMAX_ALT = {"code_capacity": 0.07, "phenomenological": 0.02}
TARGET = 1e-12
FW_DS = [3, 5, 7, 9, 11, 13, 15]
FW_EXTRA = 4            # weights t+1 .. t+1+FW_EXTRA
FW_MIN, FW_MAX, FW_FAILS = 200_000, 20_000_000, 400
P_TARGETS = [0.001, 0.005, 0.01]


def run_point(args):
    model, d, p, shots_min, shots_max, target_fails, seed = args
    rng = np.random.default_rng(seed)
    f = sample_code_capacity if model == "code_capacity" else sample_phenomenological
    t0 = time.time()
    fails, shots = 0, 0
    chunk = shots_min
    while True:
        fails += f(d, p, chunk, rng)
        shots += chunk
        if fails >= target_fails or shots >= shots_max:
            break
        chunk = min(shots_min, shots_max - shots)
    return dict(model=model, d=d, p=p, shots=shots, fails=fails, seconds=time.time() - t0)


def run_fw(args):
    d, w, seed = args
    rng = np.random.default_rng(seed)
    k = N = 0
    while k < FW_FAILS and N < FW_MAX:
        dk, dn = failure_fraction(d, w, min(FW_MIN * max(1, N // FW_MIN), FW_MAX - N), rng)
        k += dk; N += dn
    return d, w, k, N


def stratified_part(a, pool, rows):
    jobs = []
    for d in FW_DS:
        t = (d - 1) // 2
        for w in range(t + 1, t + 2 + FW_EXTRA):
            jobs.append((d, w, np.random.SeedSequence([a.seed, 7, d, w]).generate_state(1)[0]))
    jobs.sort(key=lambda j: -j[0])
    fw = {}
    for d, w, k, N in pool.imap_unordered(run_fw, jobs):
        fw.setdefault(d, {})[w] = (k, N)
    out = dict(f_w={str(d): {str(w): list(v) for w, v in sorted(fw[d].items())} for d in sorted(fw)})
    est = {}
    for p in P_TARGETS + [0.01, 0.015, 0.02]:
        est[f"{p:g}"] = {str(d): stratified_pl(d, p, fw[d]) for d in sorted(fw)}
    out["p_L"] = est
    # validation against direct Monte Carlo where both exist and truncation is small
    val = []
    for r_ in rows:
        if r_["model"] == "code_capacity" and r_["p"] in (0.01, 0.015, 0.02) and r_["d"] in fw:
            s = est[f"{r_['p']:g}"][str(r_["d"])]
            lo, hi = an.wilson(r_["fails"], r_["shots"])
            val.append(dict(d=r_["d"], p=r_["p"], mc=r_["fails"] / r_["shots"], mc_lo=float(lo), mc_hi=float(hi),
                            fw=s["p_L"], fw_lo=s["lo"], fw_hi=s["hi"], trunc=s["trunc"]))
    out["validation_vs_mc"] = val
    need = {}
    for p in P_TARGETS:
        e = est[f"{p:g}"]
        conv = {d: v for d, v in e.items() if v["trunc"] < 0.1 * TARGET}
        d_c = next((int(d) for d, v in conv.items() if v["p_L"] <= TARGET), None)
        d_u = next((int(d) for d, v in conv.items() if v["hi"] <= TARGET), None)
        need[f"{p:g}"] = dict(d_central=d_c, d_upper95=d_u, converged_ds=[int(x) for x in conv],
                              **({"full": an.qubits(d_u)["full"], "z_only": an.qubits(d_u)["z_only"]} if d_u else {}))
    out["required"] = need
    return out


def simulate(model, ps, ds, a, pool):
    jobs = []
    for i, d in enumerate(ds):
        for j, p in enumerate(ps):
            seed = [a.seed, 0 if model == "code_capacity" else 1, d, j]
            jobs.append((model, d, p, a.shots, max(a.max_shots, a.shots), a.target_fails,
                         np.random.SeedSequence(seed).generate_state(1)[0]))
    # biggest jobs first
    jobs.sort(key=lambda j: (-j[1], j[2]))
    return list(pool.imap_unordered(run_point, jobs))


def table(rows, model):
    rs = [r for r in rows if r["model"] == model]
    ds = sorted({r["d"] for r in rs}); ps = sorted({r["p"] for r in rs})
    K = np.zeros((len(ds), len(ps))); N = np.zeros_like(K)
    for r in rs:
        K[ds.index(r["d"]), ps.index(r["p"])] = r["fails"]
        N[ds.index(r["d"]), ps.index(r["p"])] = r["shots"]
    return np.array(ds), np.array(ps), K, N


def threshold_part(model, rows, out, label, ref):
    ds, ps, K, N = table(rows, model)
    PL = K / N
    lo, hi = FSS_WINDOW[model]
    sel = (ps >= lo) & (ps <= hi)
    P, D = np.meshgrid(ps[sel], ds)
    fss = an.fss_threshold(P.ravel(), D.ravel(), K[:, sel].ravel(), N[:, sel].ravel(),
                           p0=(lo + hi) / 2)
    big = D.ravel() >= 5  # systematics: drop d = 3
    fss_no3 = an.fss_threshold(P.ravel()[big], D.ravel()[big], K[:, sel].ravel()[big],
                               N[:, sel].ravel()[big], p0=(lo + hi) / 2)
    cross = {}
    for i in range(len(ds) - 1):
        cross[f"{ds[i]}-{ds[i+1]}"] = an.pairwise_crossing(ps[sel], PL[i, sel], PL[i + 1, sel])
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
    for ax, zoom in zip(axs, [False, True]):
        for i, d in enumerate(ds):
            m = (K[i] > 0) & (sel if zoom else np.ones_like(sel))
            lo95, hi95 = an.wilson(K[i, m], N[i, m])
            ax.errorbar(ps[m], PL[i, m], yerr=[PL[i, m] - lo95, hi95 - PL[i, m]], marker="o", ms=3,
                        capsize=2, label=f"d={d}")
        for f, name, ls in [(fss, "all d", "--"), (fss_no3, "d ≥ 5", "-")]:
            q = f["chi2"] / f["dof"]
            tag = " REJECTED" if q > 3 else ""
            ax.axvline(f["p_th"], color="k", ls=ls, lw=1,
                       label=f"FSS {name}: {100*f['p_th']:.2f}% (χ²/dof {q:.1f}){tag}")
        cv = [c for c in cross.values() if c is not None]
        if cv:
            ax.axvspan(min(cv), max(cv), color="orange", alpha=0.15, label="adjacent-d crossings (" + ", ".join(k for k, c in cross.items() if c is not None) + ")")
        ax.axvline(ref, color="grey", ls=":", lw=1, label=f"literature {100*ref:.2f}%")
        ax.set_xlabel("physical error rate p"); ax.set_ylabel("logical error rate p_L")
        if not zoom:
            ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(label + (" (zoom near threshold)" if zoom else ""))
        ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(out, f"threshold_{model}.png"), dpi=150); plt.close(fig)
    for f in (fss, fss_no3):
        f["accepted"] = bool(f["chi2"] / f["dof"] <= 3)
    cv = [c for c in cross.values() if c is not None]
    headline = dict(fss_d_ge_5=fss_no3["p_th"], fss_d_ge_5_accepted=fss_no3["accepted"],
                    crossing_range=[min(cv), max(cv)] if cv else None,
                    largest_d_crossing=cv[-1] if cv else None,
                    note="all-d FSS kept for reference; rejected if chi2/dof > 3. Bootstrap error is "
                         "statistical only; the spread of crossings is the systematic scale.")
    return dict(fss=fss, fss_without_d3=fss_no3, pairwise_crossings=cross, window=[lo, hi], headline=headline)


def scaling_part(model, rows, out, label, p_th):
    ds, ps, K, N = table(rows, model)
    below = ps < 0.8 * p_th
    vs_d = {}
    for j in np.where(below)[0]:
        r = an.fit_vs_d(ds, K[:, j], N[:, j])
        if r:
            vs_d[f"{ps[j]:g}"] = r
    vs_p = {}
    for i, d in enumerate(ds):
        r = an.fit_vs_p(ps[below], K[i, below], N[i, below])
        if r:
            vs_p[str(int(d))] = dict(r, expected=(d + 1) / 2)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
    ax = axs[0]
    for key, r in vs_d.items():
        p = float(key); j = list(ps).index(p)
        ok = K[:, j] >= 20
        ax.semilogy(ds[ok], K[ok, j] / N[ok, j], "o", ms=4)
        dd = np.linspace(3, 11, 50)
        ax.semilogy(dd, np.exp(r["a"] - r["alpha"] * dd), "-", lw=1, color=ax.lines[-1].get_color(),
                    label=f"p={p:g}: Λ={r['Lambda']:.2f}")
    ax.set_xlabel("code distance d"); ax.set_ylabel("p_L"); ax.set_title(label + ": p_L ∝ Λ^{-(d+1)/2}", fontsize=10)
    ax.legend(fontsize=7)
    ax = axs[1]
    ks = list(vs_p.keys())
    ax.errorbar([int(k) for k in ks], [vs_p[k]["beta"] for k in ks], yerr=[vs_p[k]["beta_err"] for k in ks],
                fmt="o", label="fitted exponent β (p_L ∝ p^β)")
    ax.plot(ds, (ds + 1) / 2, "k--", label="(d+1)/2")
    ax.set_xlabel("d"); ax.set_ylabel("β"); ax.set_title("local exponent β (p_L ∝ p^β), p < 0.8 p_th", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out, f"scaling_{model}.png"), dpi=150); plt.close(fig)
    return dict(vs_d=vs_d, vs_p=vs_p, p_cut=float(0.8 * p_th))


def overhead_part(model, rows):
    ds, ps, K, N = table(rows, model)
    P, D = np.meshgrid(ps, ds)
    res = {}
    for tag, pmax in [("main", GLOBAL_PMAX[model]), ("alt", GLOBAL_PMAX_ALT[model])]:
        g = an.global_fit(P.ravel(), D.ravel(), K.ravel(), N.ravel(), p_max=pmax)
        g["p_max"] = pmax
        g["targets"] = {}
        for p in P_TARGETS:
            d = an.required_distance(g, p, TARGET)
            g["targets"][f"{p:g}"] = dict(d=d, p_L_pred=an.predict(g, p, d) if d else None,
                                         **(an.qubits(d) if d else {}))
        res[tag] = g
    # cross-check: direct fixed-p extrapolation in d, ln p_L = a - alpha d, where p is on the grid
    direct = {}
    for p in P_TARGETS:
        if p in list(ps):
            j = list(ps).index(p)
            f = an.fit_vs_d(ds, K[:, j], N[:, j])
            if f:
                d = next((d for d in range(3, 402, 2) if np.exp(f["a"] - f["alpha"] * d) <= TARGET), None)
                direct[f"{p:g}"] = dict(d=d, Lambda=f["Lambda"], fitted_ds=f["ds"],
                                        p_L_pred=float(np.exp(f["a"] - f["alpha"] * d)) if d else None,
                                        **(an.qubits(d) if d else {}))
    res["direct_fixed_p"] = direct
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shots", type=int, default=100_000, help="min trials per point (>=1e5 for the agreed run)")
    ap.add_argument("--max-shots", type=int, default=2_000_000, help="cap when topping up low-p points")
    ap.add_argument("--target-fails", type=int, default=200, help="top up until this many failures (or cap)")
    ap.add_argument("--models", default="code_capacity,phenomenological")
    ap.add_argument("--ds", default=",".join(map(str, DS)))
    ap.add_argument("--seed", type=int, default=20261113)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--out", default="results")
    ap.add_argument("--skip-fixed-weight", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    ds = [int(x) for x in a.ds.split(",")]
    models = a.models.split(",")
    t0 = time.time()
    rows = []
    with Pool(a.workers) as pool:
        for model in models:
            t1 = time.time()
            rows += simulate(model, P_CC if model == "code_capacity" else P_PH, ds, a, pool)
            print(f"{model}: {time.time()-t1:.1f} s", flush=True)
        strat = None
        if "code_capacity" in models and not a.skip_fixed_weight:
            t1 = time.time()
            strat = stratified_part(a, pool, rows)
            print(f"fixed-weight: {time.time()-t1:.1f} s", flush=True)
    rows.sort(key=lambda r: (r["model"], r["d"], r["p"]))
    with open(os.path.join(a.out, "raw.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["model", "d", "p", "shots", "fails", "p_L", "ci_lo", "ci_hi", "seconds"])
        w.writeheader()
        for r in rows:
            lo, hi = an.wilson(r["fails"], r["shots"])
            w.writerow(dict(r, p_L=r["fails"] / r["shots"], ci_lo=float(lo), ci_hi=float(hi)))

    summary = dict(settings=vars(a), pymatching=pymatching.__version__, numpy=np.__version__,
                   min_shots_per_point=min(r["shots"] for r in rows))
    refs = {"code_capacity": 0.103, "phenomenological": 0.0293}
    labels = {"code_capacity": "Code capacity (perfect syndrome)",
              "phenomenological": "Phenomenological (q = p, d rounds)"}
    for model in models:
        th = threshold_part(model, rows, a.out, labels[model], refs[model])
        sc = scaling_part(model, rows, a.out, labels[model], th["fss"]["p_th"])
        ov = overhead_part(model, rows)
        summary[model] = dict(threshold=th, scaling=sc, overhead=ov, literature_threshold=refs[model])
    if strat is not None:
        summary["code_capacity"]["overhead"]["fixed_weight"] = strat

    # consolidated final answer per (model, p)
    final = {}
    for model in models:
        ov = summary[model]["overhead"]
        for p in P_TARGETS:
            key = f"{p:g}"
            cand = {f"global fit ({t}, p ≤ {ov[t]['p_max']:g})": ov[t]["targets"][key]["d"] for t in ("main", "alt")}
            if key in ov.get("direct_fixed_p", {}):
                cand["direct fixed-p extrapolation in d"] = ov["direct_fixed_p"][key]["d"]
            fw = ov.get("fixed_weight", {}).get("required", {}).get(key)
            if fw and fw["d_upper95"]:
                d, method = fw["d_upper95"], "fixed-weight sum (converged; supersedes fits)"
                cand["fixed-weight sum"] = d
            else:
                ok = [v for v in cand.values() if v]
                d = max(ok) if ok else None
                method = "largest d among the candidate extrapolations (conservative)"
            final.setdefault(model, {})[key] = dict(d=d, method=method, candidates=cand,
                                                     **(an.qubits(d) if d else {}))
    summary["final_overhead"] = final

    # overhead table
    lines = ["## Final answer: smallest odd d with p_L ≤ 1e-12 (per noise model and p)", "",
             "| noise model | p | d | physical qubits 2d²−1 | Z-only patch d²+(d²−1)/2 | method | candidates |",
             "|---|---|---|---|---|---|---|"]
    for model in models:
        for key, f in final[model].items():
            lines.append(f"| {model} | {float(key)*100:g}% | {f['d']} | {f.get('full')} | {f.get('z_only')} | "
                         f"{f['method']} | {'; '.join(f'{k}: {v}' for k, v in f['candidates'].items())} |")
    lines += ["", "## Intermediate: global fit p_L = A (p/p*)^((d+1)/2); superseded where the final table differs", "",
              "| noise model | p | fit window | d needed | p_L predicted | qubits (2d²−1) | qubits Z-only (d²+(d²−1)/2) |",
             "|---|---|---|---|---|---|---|"]
    for model in models:
        for tag in ["main", "alt"]:
            g = summary[model]["overhead"][tag]
            for p, t in g["targets"].items():
                if t["d"] is None:
                    lines.append(f"| {model} | {float(p)*100:g}% | p ≤ {g['p_max']:g} ({tag}) | none ≤ 401 | – | – | – |")
                    continue
                lines.append(f"| {model} | {float(p)*100:g}% | p ≤ {g['p_max']:g} ({tag}) | {t['d']} | "
                             f"{t['p_L_pred']:.2e} | {t['full']} | {t['z_only']} |")
    lines += ["", "Cross-check: direct extrapolation in d of ln p_L = a − αd measured at the same p (only where p is on the grid)", "",
              "| noise model | p | d fitted | Λ | d needed | qubits (2d²−1) |", "|---|---|---|---|---|---|"]
    for model in models:
        for p, t in summary[model]["overhead"]["direct_fixed_p"].items():
            lines.append(f"| {model} | {float(p)*100:g}% | {','.join(str(int(x)) for x in t['fitted_ds'])} | "
                         f"{t['Lambda']:.2f} | {t['d']} | {t.get('full')} |")
    if strat is not None:
        lines += ["", "Code capacity, fixed-weight (stratified) estimate p_L = Σ_w P(W=w) f_w, w = t+1..t+5; "
                  "only d with truncation bound < 1e-13 are used", "",
                  "| p | d (central ≤ 1e-12) | d (95% upper ≤ 1e-12) | qubits (2d²−1) | p_L at d−2 | p_L at d |",
                  "|---|---|---|---|---|---|"]
        for p, t in strat["required"].items():
            e = strat["p_L"][p]
            du = t["d_upper95"]
            f = lambda d: (f"{e[str(d)]['p_L']:.2e} [{e[str(d)]['lo']:.2e}, {e[str(d)]['hi']:.2e}]"
                           if d and str(d) in e else "–")
            conv = "" if du else f" (not resolved; converged d: {t['converged_ds']})"
            lines.append(f"| {float(p)*100:g}% | {t['d_central']} | {du}{conv} | {t.get('full', '–')} | "
                         f"{f(du - 2) if du else '–'} | {f(du)} |")
    with open(os.path.join(a.out, "overhead_table.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    summary["total_seconds"] = time.time() - t0
    with open(os.path.join(a.out, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=float)
    for model in models:
        f = summary[model]["threshold"]["fss"]
        print(f"{model}: p_th = {100*f['p_th']:.3f} ± {100*f['p_th_err']:.3f} %  (chi2/dof {f['chi2']:.1f}/{f['dof']}, nu {f['nu']:.2f})")
        f3 = summary[model]["threshold"]["fss_without_d3"]
        print(f"  without d=3: p_th = {100*f3['p_th']:.3f} ± {100*f3['p_th_err']:.3f} %  (chi2/dof {f3['chi2']:.1f}/{f3['dof']})")
        print("  headline:", summary[model]["threshold"]["headline"])
        print("  crossings:", summary[model]["threshold"]["pairwise_crossings"])
        print("  Lambda:", {k: round(v["Lambda"], 2) for k, v in summary[model]["scaling"]["vs_d"].items()})
        print("  beta:", {k: round(v["beta"], 2) for k, v in summary[model]["scaling"]["vs_p"].items()})
    print("\n".join(lines))
    print(f"total {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()
