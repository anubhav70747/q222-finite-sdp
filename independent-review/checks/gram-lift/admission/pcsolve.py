"""pcsolve: congruence ('shear') map used by the affine-invariant (G-metric) Newton projection in exact_primal.py.
Only tmap is kept: for the Cholesky factor L of a nearby primal solution, the matrix T with  triu-vec(L H L^T) = T triu-vec(H)."""
import numpy as np


def tmap(L):
    """matrix T (n x n), n = m(m+1)/2, with  triu-vec(L H L^T) = T triu-vec(H)  (unscaled triu coords, column-major, idx = q(q+1)/2 + p, p<=q)."""
    m = L.shape[0]; q_, p_ = np.array([(q, p) for q in range(m) for p in range(q + 1)]).T      # column-major triu: q>=p
    P = p_; Q = q_            # target (p,q), p<=q ; source (r,s) same list
    LP = L[P]; LQ = L[Q]      # (n, m)
    T = LP[:, P]*LQ[:, Q] + np.where(P[None, :] != Q[None, :], LP[:, Q]*LQ[:, P], 0.0)
    return T

