# PLAN — rotated surface code under bit-flip noise, MWPM decoding

Requested by held-out runner.

## Question
For the rotated surface code (d = 3, 5, 7, 9, 11) with independent X (bit-flip) noise, decoded with
minimum-weight perfect matching (PyMatching 2):
- **A** logical error rate p_L vs physical p for each d; threshold p_th from the crossing.
- **B** below threshold, scaling of p_L with d; extract the exponent.
- **C** phenomenological noise (data errors p and syndrome-measurement errors q = p each round, d noisy
  rounds + one final perfect readout); new threshold.
- **D** physical qubits needed for p_L = 1e-12 at p = 0.1 %, 0.5 %, 1 %.

## System (fixed)
- Rotated surface code, d×d data qubits at (r, c). Faces (i, j), 0 ≤ i, j ≤ d, act on the existing data
  qubits among (i−1, j−1), (i−1, j), (i, j−1), (i, j). Z-type faces: (i+j) even, bulk or on the top/bottom
  rows (i = 0 or d). This gives (d²−1)/2 Z checks; X checks are not needed for X noise.
- Logical X̄ = X on a row; logical Z̄ = Z on column 0. A residual X error (error ⊕ correction, zero
  syndrome) is a logical failure iff it has odd overlap with column 0.
- Code capacity (A, B): one perfect syndrome round. Edge weights ln((1−p)/p).
- Phenomenological (C): rounds t = 0..d−1: X error with prob p on each data qubit, then Z-check
  measurement flipped with prob q = p; then a perfect round (data readout). Detectors = differences of
  consecutive syndrome rounds, (d+1) layers. Space-time check matrix built explicitly; weights ln((1−p)/p)
  for both edge types. Failure rate is per block of d rounds.
- Physical qubit count of a distance-d rotated patch: n = d² data + (d²−1) ancillas = 2d² − 1 (both
  check types, since a real memory needs both). (d²−1)/2 + d² for a Z-only patch is also listed.

## Numerics
- ≥ 1e5 trials per (d, p) point (the requested minimum; never fewer). Low-p points may use more shots.
- Grids: code capacity p ∈ {0.01..0.05} ∪ dense 0.07–0.13; phenomenological p ∈ {0.005..0.02} ∪
  dense 0.022–0.036. Exact values in run.py.
- Error bars: Wilson 95 % intervals.

## Analysis
- **Threshold (A, C)**: finite-size scaling fit p_L = a + b x + c x², x = (p − p_th) d^{1/ν}, on points
  near the crossing; uncertainty by bootstrap over binomial resampling. Pairwise crossings of adjacent-d
  curves are reported too.
- **Scaling (B)**: at fixed p below threshold, fit ln p_L = a − α d (equivalently p_L ∝ Λ^{−(d+1)/2},
  Λ = e^{2α}); also the power-law exponent of p_L in p at each d, compared to the expected (d+1)/2.
  Only points with ≥ 20 failures are fitted.
- **Overhead (D)**: global fit p_L = A (p / p_*)^{(d+1)/2} on below-threshold points with ≥ 20 failures,
  for each noise model; smallest odd d with p_L ≤ 1e-12 at each p; n = 2d² − 1. This is an
  extrapolation over many decades and is labelled so.

## What counts as the answer
- Code-capacity threshold within ~10.3 % (Wang, Harrington, Preskill, Ann. Phys. 303, 31 (2003);
  Dennis et al., J. Math. Phys. 43, 4452 (2002)) to within a few tenths of a percent.
- Phenomenological threshold near ~2.9–3.0 % (Wang–Harrington–Preskill 2003, 2.93 %).
- A disagreement is reported as found, not tuned away.

## Scope
- Bit-flip only (no Y/Z, no circuit-level noise). If runtime is too large, shots may be reduced only to
  the 1e5 floor; d range and p grids are not cut without asking.
