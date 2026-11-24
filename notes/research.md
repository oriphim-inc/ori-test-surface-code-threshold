# Reference values (standard literature; written by Ori, no web lookup)

- **Code-capacity (perfect syndrome) bit-flip threshold, MWPM: ≈ 10.3 %.**
  C. Wang, J. Harrington, J. Preskill, "Confinement-Higgs transition in a disordered gauge theory and the
  accuracy threshold for quantum memory", Ann. Phys. 303, 31 (2003), arXiv:quant-ph/0207088 (toric code,
  large L). Also E. Dennis, A. Kitaev, A. Landahl, J. Preskill, J. Math. Phys. 43, 4452 (2002),
  arXiv:quant-ph/0110143.
- **Optimal (maximum-likelihood) decoder, code capacity: ≈ 10.9 %** (Nishimori point of the random-bond
  Ising model), Dennis et al. 2002; Honecker, Picco, Pujol, PRL 87, 047201 (2001). MWPM cannot exceed it.
- **Phenomenological (data p, measurement q = p) MWPM threshold: ≈ 2.9–3.0 %** (Wang–Harrington–Preskill
  2003 quote 2.93 %; Dennis et al. 2002).
- **Finite-size scaling ansatz** p_L = a + b x + c x², x = (p − p_th) d^{1/ν}: Wang–Harrington–Preskill 2003;
  ν ≈ 1.5 (2D RBIM, code capacity), ≈ 1.0 (3D RPGM, phenomenological) are typical fitted values.
- **Rotated surface code layout and qubit count 2d² − 1** (d² data + d² − 1 measure qubits): Horsman, Fowler,
  Devitt, Van Meter, NJP 14, 123011 (2012); Tomita & Svore, PRA 90, 062320 (2014).
- **Low-p scaling p_L ≈ A (p/p_th)^{(d+1)/2}**: Fowler, Mariantoni, Martinis, Cleland, PRA 86, 032324 (2012).
- **PyMatching**: O. Higgott & C. Gidney, "Sparse Blossom", Quantum 9, 1600 (2025), arXiv:2303.15933.
