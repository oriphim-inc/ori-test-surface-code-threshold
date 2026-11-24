## Final answer: smallest odd d with p_L ≤ 1e-12 (per noise model and p)

| noise model | p | d | physical qubits 2d²−1 | Z-only patch d²+(d²−1)/2 | method | candidates |
|---|---|---|---|---|---|---|
| code_capacity | 0.1% | 13 | 337 | 253 | fixed-weight sum (converged; supersedes fits) | global fit (main, p ≤ 0.05): 11; global fit (alt, p ≤ 0.07): 11; fixed-weight sum: 13 |
| code_capacity | 0.5% | 19 | 721 | 541 | largest d among the candidate extrapolations (conservative) | global fit (main, p ≤ 0.05): 19; global fit (alt, p ≤ 0.07): 19 |
| code_capacity | 1% | 27 | 1457 | 1093 | largest d among the candidate extrapolations (conservative) | global fit (main, p ≤ 0.05): 25; global fit (alt, p ≤ 0.07): 25; direct fixed-p extrapolation in d: 27 |
| phenomenological | 0.1% | 15 | 449 | 337 | largest d among the candidate extrapolations (conservative) | global fit (main, p ≤ 0.015): 15; global fit (alt, p ≤ 0.02): 15 |
| phenomenological | 0.5% | 31 | 1921 | 1441 | largest d among the candidate extrapolations (conservative) | global fit (main, p ≤ 0.015): 31; global fit (alt, p ≤ 0.02): 31; direct fixed-p extrapolation in d: 31 |
| phenomenological | 1% | 53 | 5617 | 4213 | largest d among the candidate extrapolations (conservative) | global fit (main, p ≤ 0.015): 53; global fit (alt, p ≤ 0.02): 53; direct fixed-p extrapolation in d: 53 |

## Intermediate: global fit p_L = A (p/p*)^((d+1)/2); superseded where the final table differs

| noise model | p | fit window | d needed | p_L predicted | qubits (2d²−1) | qubits Z-only (d²+(d²−1)/2) |
|---|---|---|---|---|---|---|
| code_capacity | 0.1% | p ≤ 0.05 (main) | 11 | 5.61e-13 | 241 | 181 |
| code_capacity | 0.5% | p ≤ 0.05 (main) | 19 | 2.02e-13 | 721 | 541 |
| code_capacity | 1% | p ≤ 0.05 (main) | 25 | 5.51e-13 | 1249 | 937 |
| code_capacity | 0.1% | p ≤ 0.07 (alt) | 11 | 3.56e-13 | 241 | 181 |
| code_capacity | 0.5% | p ≤ 0.07 (alt) | 19 | 8.55e-14 | 721 | 541 |
| code_capacity | 1% | p ≤ 0.07 (alt) | 25 | 1.72e-13 | 1249 | 937 |
| phenomenological | 0.1% | p ≤ 0.015 (main) | 15 | 3.95e-13 | 449 | 337 |
| phenomenological | 0.5% | p ≤ 0.015 (main) | 31 | 3.13e-13 | 1921 | 1441 |
| phenomenological | 1% | p ≤ 0.015 (main) | 53 | 6.24e-13 | 5617 | 4213 |
| phenomenological | 0.1% | p ≤ 0.02 (alt) | 15 | 3.70e-13 | 449 | 337 |
| phenomenological | 0.5% | p ≤ 0.02 (alt) | 31 | 2.87e-13 | 1921 | 1441 |
| phenomenological | 1% | p ≤ 0.02 (alt) | 53 | 5.54e-13 | 5617 | 4213 |

Cross-check: direct extrapolation in d of ln p_L = a − αd measured at the same p (only where p is on the grid)

| noise model | p | d fitted | Λ | d needed | qubits (2d²−1) |
|---|---|---|---|---|---|
| code_capacity | 1% | 3,5,7 | 6.62 | 27 | 1457 |
| phenomenological | 0.5% | 3,5,7,9 | 5.37 | 31 | 1921 |
| phenomenological | 1% | 3,5,7,9,11 | 2.59 | 53 | 5617 |

Code capacity, fixed-weight (stratified) estimate p_L = Σ_w P(W=w) f_w, w = t+1..t+5; only d with truncation bound < 1e-13 are used

| p | d (central ≤ 1e-12) | d (95% upper ≤ 1e-12) | qubits (2d²−1) | p_L at d−2 | p_L at d |
|---|---|---|---|---|---|
| 0.1% | 13 | 13 | 337 | 1.06e-12 [9.67e-13, 1.14e-12] | 1.54e-14 [1.42e-14, 1.66e-14] |
| 0.5% | None | None (not resolved; converged d: [3]) | – | – | – |
| 1% | None | None (not resolved; converged d: []) | – | – | – |
