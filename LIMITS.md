# LIMITS: what this project does NOT show

Written by the Sceptic after the main run (commit 4bfc496); updated after the fix (commit 6fdd6d9), which
added the fixed-weight (stratified) code-capacity estimator and the threshold headline, and again after
commit cafdd19, which added the consolidated final answer (summary.json["final_overhead"]). The raw Monte
Carlo data (results/raw.csv) did not change between the two commits; only timings did. Inputs:
results/raw.csv, results/summary.json, results/overhead_table.md, and short targeted checks described
below. Those checks' scripts are not part of run.py; the checks on the stratified estimator are in
tests/checks/test_fixed_weight.py.

## Noise model and decoder
- **Only independent X (bit-flip) noise.** No Z or Y errors, no correlations between X and Z (as in
  depolarizing noise), and no circuit-level noise: no CNOT faults, hook errors, leakage or idle errors.
  Thresholds under circuit-level noise are roughly 0.5–1 %, far below the numbers here. The qubit counts
  in Part D are therefore **not** a hardware estimate.
- **Phenomenological model:** q = p, d noisy rounds, then one perfect readout round. Failure is counted
  **per block of d rounds**, not per round. At these sizes the crossings depend strongly on that choice.
  The block crossings are 3-5: 2.12 %, 5-7: 2.74 %, 7-9: 2.83 %, 9-11: 2.89 %. If p_L is converted to a
  per-round rate, no adjacent pair crosses anywhere on the grid, which runs up to 3.6 %. Neither
  normalisation is the asymptotic threshold.
- **The decoder is MWPM with uniform weights ln((1−p)/p).** Two boundary qubits share each boundary
  detector, so those edges really carry about 2p, but they get the same weight as bulk edges. This is
  "minimum weight" as the plan defines it, not the best weighting. It lowers small-d crossings slightly.
  This project does not show what a re-weighted, correlated or maximum-likelihood decoder would give.
- Only memory of one logical (Z̄ parity against X errors) is tested. Logical Z errors are not simulated.

## Thresholds (Parts A and C)
- **d ≤ 11 is not the asymptotic regime.** The pre-agreed all-d FSS fits (code capacity 9.535 ± 0.025 %,
  phenomenological 2.648 ± 0.006 %) both fail the plan's agreement window. The bootstrap errors cover
  only binomial noise. The true uncertainty is dominated by finite-size systematics, which are 10–40×
  larger:
  - Code capacity, pairwise crossings for d = 3–11: 8.99, 9.38, 9.80, 9.55 %. The d ≥ 5 FSS fit gives
    9.63 % (χ²/dof 1.2). A targeted check at larger d (3e5 shots per point, p = 9.5, 10.0 and 10.3 %)
    put the 11-17 crossing at about 10.05 % and the 17-25 crossing at about 10.1 %. The crossings drift
    up towards the literature 10.3 % as d grows. That supports a finite-size explanation over a
    lattice bug, but this project does not measure 10.3 %.
  - Phenomenological: the all-d FSS fit is rejected (χ²/dof = 1126/40 ≈ 28, because d = 3 crosses d = 5
    near 2.1 %, below the fit window, so the 3-5 crossing is reported as "none"). The figure and
    summary.json now label it REJECTED. The fit without d = 3 gives 2.806 ± 0.006 % (χ²/dof 64/31 ≈ 2,
    accepted by the χ²/dof ≤ 3 rule, but its error bar is still too small). The d ≥ 5 crossings rise
    2.74 → 2.83 → 2.89 %, approaching 2.93 % but not reaching it.
  - The defensible statement is a range plus a trend: code capacity ≈ 9.0–9.8 % at d ≤ 11 (≈ 10.1 % by
    d = 25 in the targeted check), rising with d; phenomenological ≈ 2.7–2.9 % for d ≥ 5, rising with d.
    A single number with a ±0.01–0.03 % error is not defensible. The `headline` field in summary.json
    (d ≥ 5 fit plus crossing range) follows this framing; the orange band in the phenomenological figure
    spans the d ≥ 5 crossings only, because 3-5 has no crossing in the window.
- Excluding d = 3 was proposed as a systematics check before any data existed (Sceptic checks,
  notes/review.md). It is a systematic estimate, not a replacement for the pre-agreed all-d fit. The
  all-d fit's disagreement with the literature stands and must be reported as found.

## Scaling (Part B)
- The fitted exponents β are **local slopes** over p ≈ 1–7.5 % (code capacity) and 0.5–2 %
  (phenomenological), using only points with ≥ 20 failures. For code capacity they fall short of
  (d+1)/2 (d = 11: 4.83 against 6). That is expected pre-asymptotic behaviour, not a measurement of the
  asymptotic exponent.
- Λ is fitted on d = 3–11, and at the lowest p only on d = 3–7 (or 3–9). Λ still changes with d at these
  sizes. The values are not asymptotic Λ.

## Overhead (Part D)
- **The reported answer is summary.json["final_overhead"] and the "Final answer" table at the top of
  results/overhead_table.md**: code capacity d = 13 / 19 / 27 (337 / 721 / 1457 qubits) and
  phenomenological d = 15 / 31 / 53 (449 / 1921 / 5617 qubits) at p = 0.1 / 0.5 / 1 %. Only the
  code-capacity 0.1 % entry is a direct, converged estimate (no extrapolation). Every other entry is an
  extrapolation of 6–10 decades and should be read as ±2 in d (see below). The global-fit table below it
  is intermediate; where it differs (code capacity 0.1 % and 1 %), it is superseded.
- How the rule picks the answer: if the fixed-weight sum converged, its 95 %-upper d; otherwise the
  largest d among the global fits (main, alt) and the direct fixed-p extrapolation. At code capacity 1 %
  the 27 comes from the direct fixed-p fit, which is not independent of the global fit. The value is
  right for a separate reason: a constant-Λ extrapolation of the converged stratified points also gives
  27 (tests/checks/test_fixed_weight.py::test_final_entries_compatible_with_stratified). At code capacity
  0.5 % no p-grid point exists, so the rule falls back on the global fits alone (19). The stratified
  extrapolation agrees (19, margin about 2×), but the rule did not use it.
- The global ansatz A·(p/p*)^{(d+1)/2} fits the data badly (χ²/dof 12–47) and is extrapolated over 6–10
  decades in p_L, from measured p_L ≥ ~1e-6 down to 1e-12. The "direct" fixed-p cross-check is **not
  independent**: at fixed p the global ansatz is also a pure exponential in d, fitted to the same data.
- **Code capacity, p = 0.1 %: d = 13 (337 qubits), not 11.** This now comes from the fixed-weight
  (stratified) estimator p_L = Σ_w P(W = w)·f_w over w = t+1..t+5, which needs no extrapolation in p
  (truncation bound ≤ 1e-18 here). It gives p_L(d = 11) = 1.06e-12 [0.97, 1.14]e-12 and
  p_L(d = 13) = 1.5e-14. The estimator was checked against exact enumeration (all 512 errors at d = 3;
  every weight-3/4 error at d = 5 and weight-4 error at d = 7), and against all 15 direct Monte Carlo
  points at p = 1–2 %, d = 3–11. The stored f_6 at d = 11 (2.781e-4) agrees with the exact count of
  weight-6 subsets of minimum-weight logical paths (1 068 280 / C(121,6) = 2.780e-4). For w = t+1 that count
  is exact, not just a bound, because a weight-(t+1) error fails only if it lies on a minimum logical.
  With f_6 exact, p_L(11) ≈ 1.05e-12 and the remaining error comes from the weight-7 term (~4 %).
  So d = 11 misses the target by about 5 %, which is narrow but firm. Since p_L ∝ p⁶ there, d = 11 would suffice
  at p ≈ 0.099 %. The answer is sharp in this sense: it sits right at a step in d.
- **The global fit underestimates p_L at every p the stratified data can check, and the bias grows
  with d.** Ratio of stratified p_L to the main-fit prediction at d = 9, 11, 13, 15: 1.7–2.8 at 0.1 %,
  1.7–2.3 at 0.5 %, 1.5–1.9 at 1 % (alt fit: 2.2–5.5). At 0.5 % and 1 %, the stratified sum is not
  converged for d ≥ 13 (truncation too large), so those Part D entries still come from the global fit.
  Constant-Λ extrapolation of the converged stratified points (d = 7–11 at 1 %, Λ ≈ 6.7 per step; d = 7–13
  at 0.5 %, Λ ≈ 13.4) gives:
  - 0.5 %: d = 19 (p_L ≈ 5e-13), the same as the global fit, but the margin is about 2×, not 5×.
  - 1 %: d = 27 (p_L(25) ≈ 1.5e-12), against 25 from the global fit. The direct fixed-p fit also gives 27.
    A Λ that keeps rising slowly with d (6.4 → 6.9 over d = 3–13) could bring d = 25 back to the target,
    but it would take a rise that the data do not show. Read code capacity at 1 % as **27 (1457 qubits),
    with 25 possible**. The final answer now reports 27. The intermediate global-fit entry (25) is kept as
    a strict expected failure in tests/checks/test_fixed_weight.py, because its failure is a finding.
- **Phenomenological Part D (15 / 31 / 53) has no stratified check.** It uses the same global ansatz,
  which underestimates p_L for code capacity by 1.5–3×. Phenomenological at 0.1 % (d = 15) has only about
  2.5× margin over d = 13 in the fit. Each phenomenological entry could be one step (+2) too small.
  Read them as lower-leaning estimates, ±2 in d. The "conservative" rule (largest candidate) does not
  help here: all three candidates agree (15/15, 31/31/31, 53/53/53) because they come from the same data
  and the same exponential-in-d form. Agreement among them is not evidence that the entries are right.
- Qubit counts 2d² − 1 assume a single rotated patch with both check types and no routing, factory or
  lattice-surgery overhead.

## Statistics
- Points stop at 200 failures or 2e6 shots. That biases p_L upward by well under 1 %, which is
  negligible. Low-p points at large d with few failures are excluded from the fits (< 20 failures). The
  lowest-p behaviour of d = 9 and 11 is therefore not measured by direct Monte Carlo. For code capacity
  the fixed-weight estimator covers it; for the phenomenological model it remains unmeasured.
- The stratified estimator's lower 95 % limit ignores the truncation (it is a partial sum, so it is a
  lower bound up to statistics); the upper limit adds the full truncation bound P(W > w_max). Strata with
  0 failures would get an almost-zero variance (floor 1/N² instead of a rule-of-three bound). No stratum
  used at p = 0.1 % has fewer than 46 failures, so this does not affect the result.
