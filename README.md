> **Surface-code error-correction thresholds.** A project Ori built and ran on a problem from Oriphim's held-out test set, which Ori was never tuned on. Short paper: [paper/paper.pdf](paper/paper.pdf). The Sceptic's review is in [notes/review.md](notes/review.md); limits are in [LIMITS.md](LIMITS.md). [oriphim.com](https://oriphim.com)
>
> How this repo was made: Ori built and ran the project and wrote the paper. Afterwards, Oriphim added the licence, edited the author line and wrote this note, with help from Claude Code.

# Rotated surface code under bit-flip noise, decoded with MWPM

By Oriphim, with Ori.

**Question.** For the rotated surface code (d = 3, 5, 7, 9, 11) with independent X errors, decoded with
minimum-weight perfect matching (PyMatching 2.4.0), we asked four things:
- **A.** Where is the threshold?
- **B.** How does p_L scale with d below threshold?
- **C.** What is the threshold under phenomenological noise (d noisy syndrome rounds)?
- **D.** How many physical qubits does p_L = 1e-12 need at p = 0.1 %, 0.5 % and 1 %?

The agreed plan is in [PLAN.md](PLAN.md).

## Answer

**A. Code-capacity threshold ≈ 9.5–9.8 % at these sizes, below the literature's 10.3 %.**
- The adjacent-d crossings are 8.99, 9.38, 9.80 and 9.55 % (for d pairs 3-5, 5-7, 7-9, 9-11).
- The finite-size-scaling fit gives 9.54 % with all d, and 9.63 % with d ≥ 5 (χ²/dof 1.2).
- The 10.3 % value is for MWPM in the large-size limit (Wang, Harrington & Preskill 2003). Our number
  misses the window agreed in PLAN.md, and we report that as a finding, not a bug.
- This looks like a finite-size effect of the planar boundaries. In a short check outside `run.py`, the
  Sceptic found the crossings keep rising with d: about 10.05 % for d = 11-17 and about 10.1 % for d = 17-25.

**C. Phenomenological threshold (q = p, failure counted per block of d rounds) ≈ 2.8–2.9 %, rising towards
2.93 %.**
- The crossings are 2.74 → 2.83 → 2.89 % for d pairs 5-7, 7-9 and 9-11.
- The fit with d ≥ 5 gives 2.81 % (χ²/dof 2.1).
- The fit with all d is rejected (χ²/dof 28). d = 3 crosses d = 5 only at about 2.1 %.

**B. Scaling below threshold.** p_L falls exponentially in d, p_L ∝ Λ^{-(d+1)/2}. Λ is the suppression
factor per distance step of 2:

| p | Λ, code capacity | Λ, phenomenological |
|---|---|---|
| 0.5 % | – | 5.37 ± 0.21 |
| 1 % | 6.62 ± 0.36 | 2.59 ± 0.04 |
| 2 % | 3.32 ± 0.06 | 1.28 ± 0.01 |
| 5 % | 1.48 ± 0.01 | – |

Λ falls towards 1 as p approaches threshold.

The exponent of p in p_L ∝ p^β, fitted at each d, compared with the expected (d+1)/2:

| d | 3 | 5 | 7 | 9 | 11 |
|---|---|---|---|---|---|
| (d+1)/2 | 2 | 3 | 4 | 5 | 6 |
| β, code capacity | 1.83 | 2.67 | 3.40 | 4.16 | 4.83 |
| β, phenomenological | 1.85 | 2.72 | 3.88 | 4.90 | 6.09 |

These β are local slopes at the measured p (well above 1e-3), not the asymptotic exponent. The exact
weight-by-weight count below confirms that the minimum failing weight is (d+1)/2.

**D. Physical qubits for p_L ≤ 1e-12.** Each entry is the smallest odd d, with n = 2d² − 1 qubits
(data plus both kinds of ancilla).

| noise model | p = 0.1 % | p = 0.5 % | p = 1 % |
|---|---|---|---|
| code capacity | **d = 13 → 337** (exact weight-by-weight count) | d = 19 → 721 | d = 27 → 1457 |
| phenomenological (per d-round block) | d = 15 → 449 | d = 31 → 1921 | d = 53 → 5617 |

Only the code-capacity entry at 0.1 % is a direct estimate:
- It comes from a fixed-weight sum, p_L = Σ_w P(W=w)·f_w, where f_w is the fraction of weight-w errors
  that MWPM fails on.
- It gives p_L(d=11) = 1.06e-12, just above the target, and p_L(d=13) = 1.5e-14.
- The simple global fit A(p/p*)^{(d+1)/2} had said d = 11. That fit under-predicts p_L by 1.5–3×, and the
  error grows with d.

Every other entry extrapolates 6–10 decades and should be read as **±2 in d**. Where the estimates
differ, the table takes the larger d. Full candidates: [results/overhead_table.md](results/overhead_table.md).
These are bit-flip-only numbers. Real hardware with circuit-level noise needs far more (see LIMITS.md).

## Evidence
- **Simulation size:** ≥ 1e5 trials at every (d, p) point. Points with fewer than 200 failures were
  topped up to as many as 2e6. Counts, rates and Wilson 95 % intervals are in `results/raw.csv`.
- **Figures:**
  - `results/threshold_code_capacity.png` and `results/threshold_phenomenological.png`: p_L against p,
    with a zoom on the crossing.
  - `results/scaling_code_capacity.png` and `results/scaling_phenomenological.png`: p_L against d, and β
    against d.
- **All fitted numbers:** `results/summary.json`. The final Part D answer is under `final_overhead`.
- **Correctness checks:** `tests/test_code.py` (14 tests).
  - The code structure is right: X and Z checks commute, there is one logical qubit, and d = 3 has
    distance 3 by brute force.
  - MWPM corrects every error of weight ≤ t, for d ≤ 7 and for the space-time graph.
  - The space-time matrix matches an explicit round-by-round simulation.
  - The fixed-weight estimator matches exact enumeration of all 2⁹ errors at d = 3.
- **The Sceptic's checks (`tests/checks/`):**
  - An independent networkx matching decoder agrees with PyMatching.
  - The d = 3 results agree with exact values.
  - The weight-by-weight failure fractions agree with exact enumeration at d = 5 and 7, and with the Monte
    Carlo at 15 points.
- **9 of those pre-agreed checks fail, and stay failing as findings:** both thresholds miss the plan's
  literature window; the all-d phenomenological fit is poor; there is no d = 3/5 crossing inside the
  window; and the global Part D fit has a high χ².
- **Reviews:** `notes/review.md`. Literature values: `notes/research.md`.

## What this does not show
See **[LIMITS.md](LIMITS.md)**. In short:
- X noise only, with no circuit-level noise.
- d ≤ 11 is not the large-size regime, so the thresholds are finite-size values.
- MWPM uses uniform edge weights.
- Part D is an extrapolation, except for the code-capacity entry at 0.1 %.

## Re-run
```bash
python3 -m venv .venv && .venv/bin/pip install -e '.[test]'
.venv/bin/python run.py            # ≈ 6 min on 8 cores (MC ≈ 35 s, fixed-weight stage ≈ 5 min)
MPLCONFIGDIR=$TMPDIR .venv/bin/python -m pytest tests -q
```
Options: `--shots`, `--max-shots`, `--target-fails`, `--ds`, `--models`, `--seed`,
`--skip-fixed-weight`. Seeds are fixed, so a re-run reproduces `results/raw.csv` exactly.
