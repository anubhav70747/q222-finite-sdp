"""centre.py: double-precision helpers on a rs_core.Structure (svec scaling, cone shift, dependent-row detection, Gram unpacking).
Only the helpers used by the admission check are kept.  Everything computed here is a *candidate*; the certificate re-verifies
every claim that matters in exact rational / ball arithmetic (exact_primal.py)."""
import math
import numpy as np, scipy.linalg as sla

SQ2 = math.sqrt(2.0)


def cx_vec(S):
    """x = Cx y : x_pp = y_pp, x_pq = y_pq/sqrt2."""
    cx = np.ones(S.ng)
    for bi, m in enumerate(S.ms):
        p_, q_ = np.triu_indices(m); pos = q_*(q_ + 1)//2 + p_
        v = np.where(p_ == q_, 1.0, 1/SQ2); cx[S.offs[bi] + pos] = v
    return cx


def shift_y(S):
    """svec (y-coordinates) of Sh_b = (V^T V)^-1 = I - a a^T/(1+a.a) (so G' - t Sh >= 0 <=> complement-orthonormal Gram >= t)."""
    sh = np.zeros(S.ng)
    for bi, m in enumerate(S.ms):
        a = np.array([float(x) for x in S.info[bi]['a']]); M = np.eye(m) - np.outer(a, a)/(1 + a @ a)
        p_, q_ = np.triu_indices(m); pos = q_*(q_ + 1)//2 + p_
        sh[S.offs[bi] + pos] = M[p_, q_]*np.where(p_ == q_, 1.0, SQ2)
    return sh


def drop_rows(S, An, drop=6):
    """numerically dependent rows: returns (keep index array, dropped index array)."""
    U, sv, Vt = np.linalg.svd(An.toarray(), full_matrices=False); rk = int(np.sum(sv > 1e-10*sv[0]))
    assert An.shape[0] - rk == drop, ('rank', rk, An.shape[0])
    _, _, pv = sla.qr(U[:, rk:].T, mode='economic', pivoting=True); gone = np.sort(pv[:drop])
    return np.setdiff1d(np.arange(An.shape[0]), gone), gone, sv, rk


def y_to_Gs(S, y):
    return S.unvec(y*cx_vec(S))


def Gs_scaled(S, y):
    """Gs = V G' V^T (current scaling) diagonals"""
    out = []
    for bi, G in enumerate(y_to_Gs(S, y)):
        inf = S.info[bi]; n = inf['n']; V = np.zeros((n, n - 1)); a = np.array([float(x) for x in inf['a']])
        for k in range(n):
            if k != inf['i0']: V[k, inf['pos'][k]] = 1.0
        V[inf['i0'], :] = -a; out.append(V @ G @ V.T)
    return out


