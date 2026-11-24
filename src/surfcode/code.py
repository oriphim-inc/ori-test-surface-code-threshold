"""Rotated surface code: Z-type check matrix, logical operators, space-time matrices."""
import numpy as np
import scipy.sparse as sp


def qubit_index(r, c, d):
    return r * d + c


def z_faces(d):
    """List of Z-type faces; each face is a list of data-qubit indices.

    Faces (i, j), 0 <= i, j <= d, touch data qubits (i-1, j-1), (i-1, j), (i, j-1), (i, j)
    that exist. Z-type: (i + j) even, and either in the bulk (1 <= i, j <= d-1) or a
    weight-2 face on the top/bottom boundary rows (i = 0 or i = d, 1 <= j <= d-1).
    """
    if d < 3 or d % 2 == 0:
        raise ValueError("d must be odd and >= 3")
    faces = []
    for i in range(d + 1):
        for j in range(d + 1):
            if (i + j) % 2:
                continue
            bulk = 1 <= i <= d - 1 and 1 <= j <= d - 1
            tb = (i == 0 or i == d) and 1 <= j <= d - 1
            if not (bulk or tb):
                continue
            qs = [qubit_index(r, c, d) for r, c in
                  ((i - 1, j - 1), (i - 1, j), (i, j - 1), (i, j))
                  if 0 <= r < d and 0 <= c < d]
            faces.append(qs)
    return faces


def x_faces(d):
    """X-type faces (only used for tests of commutation): (i+j) odd, bulk or left/right boundary."""
    faces = []
    for i in range(d + 1):
        for j in range(d + 1):
            if (i + j) % 2 == 0:
                continue
            bulk = 1 <= i <= d - 1 and 1 <= j <= d - 1
            lr = (j == 0 or j == d) and 1 <= i <= d - 1
            if not (bulk or lr):
                continue
            qs = [qubit_index(r, c, d) for r, c in
                  ((i - 1, j - 1), (i - 1, j), (i, j - 1), (i, j))
                  if 0 <= r < d and 0 <= c < d]
            faces.append(qs)
    return faces


def faces_to_matrix(faces, n):
    H = np.zeros((len(faces), n), dtype=np.uint8)
    for k, qs in enumerate(faces):
        H[k, qs] = 1
    return H


def hz(d):
    """Z-check matrix ((d^2-1)/2 x d^2) that detects X errors."""
    return faces_to_matrix(z_faces(d), d * d)


def logical_z(d):
    """Logical Z = Z on column 0. An X residual is a logical error iff odd overlap."""
    L = np.zeros((1, d * d), dtype=np.uint8)
    L[0, [qubit_index(r, 0, d) for r in range(d)]] = 1
    return L


def logical_x(d):
    """Logical X = X on row 0."""
    L = np.zeros((1, d * d), dtype=np.uint8)
    L[0, [qubit_index(0, c, d) for c in range(d)]] = 1
    return L


def spacetime_matrices(d, rounds):
    """Phenomenological space-time check matrix.

    Rounds t = 0..rounds-1: data X errors (column block D_t) then noisy measurement
    (column block M_t). Then one perfect round. Detector layers k = 0..rounds:
    det_k = s_k xor s_{k-1} (s_{-1} = 0), s_rounds perfect.
    D_t flips Hz in layer t (subsequent layers see it in both, so the difference cancels).
    M_t flips layer t and layer t+1.
    Returns (H, L, n_data_cols, n_meas_cols): H is (m*(rounds+1)) x (n*rounds + m*rounds).
    """
    H0 = sp.csr_matrix(hz(d))
    m, n = H0.shape
    T = rounds
    blocks_D = sp.kron(sp.eye(T + 1, T, format="csr"), H0)          # layer t <- D_t
    # M_t -> layers t and t+1
    S = sp.eye(T + 1, T, format="csr") + sp.eye(T + 1, T, k=-1, format="csr")
    blocks_M = sp.kron(S, sp.eye(m))
    H = sp.hstack([blocks_D, blocks_M]).tocsc()
    H.data = H.data % 2
    H.eliminate_zeros()
    L0 = sp.csr_matrix(logical_z(d))
    L = sp.hstack([sp.kron(np.ones((1, T)), L0), sp.csr_matrix((1, m * T))]).tocsc()
    return H.astype(np.uint8), L.astype(np.uint8), n * T, m * T
