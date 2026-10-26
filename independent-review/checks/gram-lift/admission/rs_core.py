"""rs_core: structure of the contact-rescaled (C1/C2) Gram problem.  Only Structure (and what it needs) is kept from the covering pilot.

KEY FACT (exact, used throughout).  Write the Gram identity F_c(s^6,t^6,psi) = Z^T G Z in the basis rescaled by the contact
values:  D = diag(w_k(theta) s_k),  w_k = |u_c^{i/6} v_c^{j/6} psi_{l,c}|,  s_k = 2^{e_k} a FIXED dyadic constant,  Gs = D G D,
Gs = V G' V^T with V^T zh = 0, zh_k = sigma_l / s_k  (sigma_l = sign psi_{l,c}, constant on the whole C2 chart).  Then every
coefficient equation (row `key` = (p,q,l,m)) reads
        sum_{(r,c) in key} Gs[r,c]/(s_r s_c)  =  beta_key(theta) := mu_key(theta) f_key(theta),
        mu_key = u_c^{p/6} v_c^{q/6} |psi_l psi_m|,     f_key = coefficient of u^{p/6} v^{q/6} psi_l psi_m in F_c (0 unless 6|p,q).
So the coefficient matrix A_0 is a CONSTANT dyadic matrix independent of theta; only the right-hand side beta(theta) (a RATIONAL
function of the chart parameters, supported on <= 90 keys) moves.
"""
import os, sys
from fractions import Fraction as Fr
import numpy as np, scipy.sparse as sps

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_core import NAMES, _TFR                    # T tensor

ACTIVE = ['d', 'n', 'm', 'C', 'A', 'D', 'alpha']   # C2 chart: lam = beta = 0 identically
SIGMA = (1, -1, 1, 1)                               # sign of psi_c = (b/a, -c r, c h/a, 1) on the whole C2 chart


# ---------------------------------------------------------------------------------------------- exact data on weights
def weight_tensor():
    """TW[(a,b,l,m)] = {w: Fraction}: coefficient of u^a v^b psi_l psi_m (l<=m) in F_c is sum_w TW[..][w]*c_w."""
    TW = {}
    for (p, a, b, l, m), tv in _TFR.items():
        if l > m or NAMES[p] not in ACTIVE: continue
        d = TW.setdefault((a, b, l, m), {}); d[NAMES[p]] = d.get(NAMES[p], Fr(0)) + tv*(1 if l == m else 2)
    return {k: {w: x for w, x in d.items() if x != 0} for k, d in TW.items() if any(x != 0 for x in d.values())}


TW = weight_tensor()


def beta_keys():
    """keys (p,q,l,m) = (6a,6b,l,m) with theta-dependent right-hand side."""
    return sorted((6*a, 6*b, l, m) for (a, b, l, m) in TW)


# ---------------------------------------------------------------------------------------------- rows of the Gram map
def build_rows(blocks):
    rows = {}
    for bi, blk in enumerate(blocks):
        for r, (ir, jr, lr) in enumerate(blk):
            for c, (ic, jc, lc) in enumerate(blk):
                rows.setdefault((ir + ic, jr + jc, min(lr, lc), max(lr, lc)), []).append((bi, r, c))
    return rows


class Structure:
    """Constant dyadic coefficient matrix A_0 (rows = keys) for fixed blocks and fixed dyadic scalings 2^e."""
    def __init__(self, blocks, e):
        self.blocks = blocks; self.e = e; nb = len(blocks)
        self.rows = build_rows(blocks)
        self.keys = sorted(self.rows)                       # rows
        self.kidx = {k: i for i, k in enumerate(self.keys)}
        missing = [k for k in beta_keys() if k not in self.rows]
        assert not missing, 'F has monomials not reachable by the basis: %s' % missing
        self.info = []; ms = []
        for bi, blk in enumerate(blocks):
            n = len(blk); zh = [Fr(SIGMA[l])/Fr(2)**e[bi][k] for k, (_, _, l) in enumerate(blk)]
            i0 = max(range(n), key=lambda k: abs(zh[k]))     # pivot = largest |zh| (smallest s_k)
            pos = [-1]*n; cnt = 0
            for k in range(n):
                if k != i0: pos[k] = cnt; cnt += 1
            a = [zh[k]/zh[i0] for k in range(n) if k != i0]
            self.info.append(dict(n=n, i0=i0, pos=pos, a=a, zh=zh)); ms.append(n - 1)
        self.ms = ms; self.dims = [m*(m + 1)//2 for m in ms]; self.offs = np.concatenate([[0], np.cumsum(self.dims)]).astype(int)
        self.ng = int(self.offs[-1]); self._build()

    def _build(self):
        ent = [dict() for _ in self.keys]                    # per row: column -> Fraction
        for ki, key in enumerate(self.keys):
            d = ent[ki]
            for (bi, r, c) in self.rows[key]:
                inf = self.info[bi]; i0 = inf['i0']; pos = inf['pos']; a = inf['a']; e = self.e[bi]
                coef = Fr(1, 1)/(Fr(2)**(e[r] + e[c]))
                vr = {pos[r]: Fr(1)} if r != i0 else {p: -a[p] for p in range(len(a))}
                vc = {pos[c]: Fr(1)} if c != i0 else {p: -a[p] for p in range(len(a))}
                for p, x in vr.items():
                    for q, y in vc.items():
                        lo, hi = (p, q) if p <= q else (q, p); col = int(self.offs[bi]) + hi*(hi + 1)//2 + lo
                        d[col] = d.get(col, Fr(0)) + coef*x*y
        self.ent = [{c: x for c, x in d.items() if x != 0} for d in ent]
        I, J, X = [], [], []
        for ri, d in enumerate(self.ent):
            for c, x in d.items(): I.append(ri); J.append(c); X.append(float(x))
        self.A = sps.csr_matrix((X, (I, J)), shape=(len(self.keys), self.ng))   # x = triu entries of G' (unscaled), column-major triu

    def unvec(self, x):
        """list of symmetric matrices G'_b from x (triu entries, unscaled)."""
        Gs = []
        for bi, m in enumerate(self.ms):
            seg = x[self.offs[bi]:self.offs[bi + 1]]; G = np.zeros((m, m)); p_, q_ = np.triu_indices(m)
            G[p_, q_] = seg[q_*(q_ + 1)//2 + p_]; Gs.append(G + G.T - np.diag(np.diag(G)))
        return Gs


