"""hcore: vectorised native C1/C2 chart, constant-A_0 beta(theta), margin SDP with dual extraction (the float stage of the admission check).

mu convention (same as rs-cover/nd-survey): in contact-rescaled coordinates (w_k = |contact value of feature k|, then dyadic scaling 2^e),
  Gs = V G' V^T,  A_e x = beta(theta),  mu_e(theta) = max{t : G' - t (V^T V)^-1 >= 0}  (= lambda_min of the Gram on the orthonormal contact complement),
  dimensionless version  t_rel = lambda_min / max-diagonal of that Gram.
The dual y (key space, 412 rows) satisfies  <y, beta> = mu_e, and y in D (cone is e-independent: A_e^* y = S^-1 A_0^* y S^-1).
"""
import os, sys, math
import numpy as np, scipy.sparse as sps, clarabel
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rs_core, nd_core, gramtools
from rs_core import Structure
from centre import cx_vec, shift_y, y_to_Gs, Gs_scaled, drop_rows
from nd_core import _TFR

# ------------------------------------------------------------------------------------------ coefficient tensor (all 9 weights)
def weight_tensor_all():
    TW = {}
    for (p, a, b, l, m), tv in _TFR.items():
        if l > m: continue
        d = TW.setdefault((a, b, l, m), {}); d[p] = d.get(p, Fr(0)) + tv*(1 if l == m else 2)
    return {k: {w: x for w, x in d.items() if x != 0} for k, d in TW.items() if any(x != 0 for x in d.values())}
TWA = weight_tensor_all()
TKEYS = sorted(TWA)                                              # (a,b,l,m)
TMAT = np.zeros((len(TKEYS), 9))
for i, k in enumerate(TKEYS):
    for w, x in TWA[k].items(): TMAT[i, w] = float(x)
TA = np.array([k[0] for k in TKEYS]); TB = np.array([k[1] for k in TKEYS]); TL = np.array([k[2] for k in TKEYS]); TM = np.array([k[3] for k in TKEYS])
ROWKEYS = [(6*a, 6*b, l, m) for (a, b, l, m) in TKEYS]

# ------------------------------------------------------------------------------------------ vectorised native sources
def c2_vec(P):
    """P: (5,N) c,a,b,t,r  ->  W (9,N) weights (d,lam,n,m,C,A,D,alpha,beta), CT (6,N) contact (u,v,W,X,Y,Z), chart mask (N,), gaps dict"""
    c, a, b, t, r = P; z = c*c; Ac = 1 + z*r; h = b - a*a*t; g = (1 + z)*b - (1 - r + a*a*z*t); j = 1 + z*r - t
    Dl = ((1 + z)*z*b*b + 2*z*(1 + z)*b*r + (1 + z)*(1 + a*a*z)*b + z*z*(1 - a*a)*r*r + 2*z*(1 + a*a*z)*r + 1 + a*a*z)
    N = (1 + z)*(b + r)*(b + a*a*(1 - t)) - a*a*(1 + a*a*z)*(1 - r)*t
    one = np.ones_like(c); zero = np.zeros_like(c)
    with np.errstate(all='ignore'):
        W = np.array([one, zero, ((b + r)*Ac*(Ac + a*a*z*t) - r*Dl)/(2*r*Dl), (1 + a*a*z)*(1 - r)*Ac*(Ac + a*a*z*t)/(2*g*Dl), Ac*N/(2*a*a*t*Dl),
                      a*a*Ac*(b + r)*(Ac + a*a*z*t)/(2*h*Dl), (1 + z)*(b + r)*(b - a*a*z*r)*(Ac + a*a*z*t)/(2*a*a*j*Dl),
                      z*z*(b + r)*N/(2*(1 - z*h)*Dl), zero])
        CT = np.array([c, a*c, b/a, -c*r, c*(b/a - a*t), one])
        gp = {'h': h/b, 'g': g/((1 + z)*b), 'j': j/(1 + z*r), 'zh': 1 - z*h, 'bz': (b - a*a*z*r)/b, 'r0': r, 'r1': 1 - r}
        ok = (a > 0) & (c > 0) & (b > 0) & (t > 0) & (r > 0) & (r < 1) & (h > 0) & (g > 0) & (j > 0) & (1 - z*h > 0) & (b - a*a*z*r > 0)
        ok &= np.all(W[1:] >= 0, axis=0) & np.all(np.isfinite(W), axis=0) & np.all(np.isfinite(CT), axis=0)
    return W, CT, ok, gp


def c1_vec(P):
    k, f, g, a, t = P; kap, eta, zeta = k*k, f*f, g*g; V = kap*eta*zeta; Lc = 1 - a*eta; Hc = 1 + (1 - a)*eta; Qc = 1 - kap*zeta + a*V; K = 1 + kap
    one = np.ones_like(k); zero = np.zeros_like(k)
    with np.errstate(all='ignore'):
        e = k*f/g; Dl = 1/e; R = 1 + kap*(1 - a*eta); T = 1 + kap - kap*zeta*t; M = R - t; J = T - a
        P4 = 1 + kap + eta*kap*(a - kap - 1) - zeta*kap*t + V*(-zeta*kap*t**2 + 2*kap*t + a*(1 - t) + t)
        Pal = 1 + kap - a*eta*kap + zeta*((1 + kap) - (2 + kap)*t) + V*(-zeta*t**2 + a*(t - 1) + t)
        Pd = (1 + kap)**2 - a*eta*kap**2 + zeta*kap*((kap + 1) - (kap + 2)*t) + V*(-zeta*kap*t**2 + kap*a*(t - 1) + kap*t + kap + 1)
        c2 = 2*K*(1 - t*zeta)*Lc*(1 + V)*Hc*Qc; cb = 2*K*Lc*(1 + (1 - t)*zeta)*Hc*Qc; c3 = 2*K*Lc*Hc*Qc*(1 - kap*eta + t*V)
        c5 = 2*K*Lc*Hc*Qc*(1 + V); c4 = 2*Qc*P4; cal = 2*Lc*Hc*Pal; cd = 4*Lc*Qc*Pd
        dB = 2*a*eta*cd/c3
        W = np.array([one, zero, zero, 2*a*e*e*c2/(M*c3)/dB, 2/dB, 2*a*eta**2*c4/c3/dB, 2*a*c5/(J*c3)/dB - 0.5,
                      2*a*Dl**2*k**6*eta**2*cal/c3/dB, 2*a*e*e*cb/(t*c3)/dB - 0.5])
        CT = np.array([k, Dl*k*k*eta, Dl*k, -k*(1 - t*zeta), Dl*k*k*Lc, one])
        gp = {'Lc': Lc/1, 'Hc': Hc/(1 + eta + a*eta), 'Qc': Qc/(1 + kap*zeta + a*V), 'M': M/R, 'J': J/T, 'e3': (1 - kap*eta + t*V)/(1 + kap*eta + t*V),
              'k1': 1 - k, 'k0': k, 't': t}
        ok = (k > 0) & (f > 0) & (g > 0) & (a > 0) & (t > 0) & (k <= 1) & (Lc > 0) & (Hc > 0) & (Qc > 0) & (M > 0) & (J > 0) & (c3 > 0) & (dB > 0)
        ok &= np.all(W[1:] >= 0, axis=0) & np.all(np.isfinite(W), axis=0) & np.all(np.isfinite(CT), axis=0)
    return W, CT, ok, gp

VEC = {'C2': c2_vec, 'C1': c1_vec}

def beta_all(fam, P):
    """P (5,N) -> beta (N, nTkeys) float64 (contact-rescaled rhs), chart mask, gaps"""
    P = np.atleast_2d(np.asarray(P, float).T).T if np.ndim(P) == 1 else np.asarray(P, float)
    W, CT, ok, gp = VEC[fam](P)
    f = (TMAT @ W).T                                        # (N, nkeys) coefficient of u^a v^b psi_l psi_m in F_c
    psa = np.abs(CT[2:]); uc, vc = CT[0], CT[1]
    with np.errstate(all='ignore'):
        mu = psa[TL].T*psa[TM].T*uc[:, None]**TA[None, :]*vc[:, None]**TB[None, :]
    return mu*f, ok, gp, W, CT


# ------------------------------------------------------------------------------------------ fixed structures per family
_STRUCT = {}
REF = {'C2': (Fr(1), Fr(1, 2), Fr(13, 20), Fr(2, 5), Fr(3, 10)), 'C1': None}

def ref_point(fam):
    if fam == 'C2': return REF['C2']
    # find a rational interior C1 point
    rng = np.random.default_rng(1)
    while True:
        p = (10**rng.uniform(-1, 0), 10**rng.uniform(-1, 0.5), 10**rng.uniform(-1, 0.5), 10**rng.uniform(-1, 0.5), 10**rng.uniform(-1, 0.5))
        _, _, ok, _ = c1_vec(np.array(p)[:, None])
        if ok[0]: return tuple(Fr(x).limit_denominator(1000) for x in p)


def blocks_at(fam, th):
    src, ct, chart = nd_core.source(fam, th); assert all(chart)
    return gramtools.make_blocks(nd_core.fcoef_exact(src)), src, ct


def get_blocks(fam):
    if fam not in _STRUCT:
        th = ref_point(fam); blocks, src, ct = blocks_at(fam, th); _STRUCT[fam] = dict(blocks=blocks, th=th)
    return _STRUCT[fam]['blocks']


def beta_rows(S, bt):
    """bt: (nkeys,) beta over TKEYS -> vector over S.keys (zeros elsewhere); keys must be rows"""
    out = np.zeros(len(S.keys))
    for i, k in enumerate(ROWKEYS):
        if bt[i] != 0.0: out[S.kidx[k]] += bt[i]
    return out


# ------------------------------------------------------------------------------------------ margin SDP with duals
def solve_dual(S, beta_f, tol=1e-11, maxit=400, keep=None, want_dual=True, delta=0.0):
    rn = np.sqrt(np.asarray(S.A.multiply(S.A).sum(axis=1)).ravel()); An = sps.diags(1/rn) @ S.A; bn = beta_f/rn
    if keep is None: keep, gone, sv, rk = drop_rows(S, An)
    cx = cx_vec(S); Ay = (An[keep] @ sps.diags(cx)).tocsc(); by = bn[keep]
    Ad = Ay.toarray(); Uu, ss, Vt = np.linalg.svd(Ad, full_matrices=False); by_c = Uu @ (Uu.T @ by)
    incons = float(np.linalg.norm(by - by_c)/max(np.linalg.norm(by), 1e-300)); by = by_c
    sh = shift_y(S); ng = S.ng
    Afull = sps.vstack([sps.hstack([sps.csc_matrix((len(by), 1)), Ay]), sps.hstack([sps.csc_matrix(sh.reshape(-1, 1)), -sps.identity(ng)])]).tocsc()
    qv = np.zeros(1 + ng); qv[0] = -1.0
    if delta:                                                  # maximise t + delta*tr(G'):  dual y then has A^*y >= delta*I (strictly feasible, rounds to an exact PD certificate)
        for bi, m in enumerate(S.ms): qv[1 + int(S.offs[bi]) + np.arange(m)*(np.arange(m) + 1)//2 + np.arange(m)] = -delta
    st = clarabel.DefaultSettings(); st.verbose = False; st.max_iter = maxit; st.tol_gap_abs = st.tol_gap_rel = st.tol_feas = tol
    cones = [clarabel.ZeroConeT(len(by))] + [clarabel.PSDTriangleConeT(m) for m in S.ms]
    sol = clarabel.DefaultSolver(sps.csc_matrix((1 + ng, 1 + ng)), qv, Afull, np.concatenate([by, np.zeros(ng)]), cones, st).solve()
    x = np.array(sol.x); out = dict(t=float(x[0]), y=x[1:], status=str(sol.status), keep=keep, incons=incons, iters=sol.iterations)
    if want_dual:
        z = np.array(sol.z); yk = np.zeros(len(S.keys)); yk[keep] = z[:len(by)]/rn[keep]; out['ydual'] = yk     # key-space dual:  <ydual, beta> ~ t
        out['zpsd'] = z[len(by):]
    return out


def gram_stats(S, out):
    """lambda_min of the orthonormal-complement Gram and max diagonal of Gs (current scaling)."""
    Gps = y_to_Gs(S, out['y']); mn = min(np.linalg.eigvalsh(G)[0] for G in Gps)
    Gs = Gs_scaled(S, out['y']); md = max(np.max(np.diag(G)) for G in Gs)
    return float(mn), float(md)


# ------------------------------------------------------------------------------------------ family-specific structure (assertion fix)
ACTIVE_IDX = {'C2': [0, 2, 3, 4, 5, 6, 7], 'C1': [0, 3, 4, 5, 6, 7, 8]}      # weights that are not identically zero in the chart


def active_keys(fam):
    """row keys (6a,6b,l,m) whose coefficient involves a non-identically-zero weight"""
    return [ROWKEYS[i] for i, k in enumerate(TKEYS) if any(TMAT[i, w] != 0 for w in ACTIVE_IDX[fam])]


def make_struct(fam, e=None):
    blocks = get_blocks(fam); e = e or [[0]*len(b) for b in blocks]
    old = rs_core.beta_keys; rs_core.beta_keys = lambda: active_keys(fam)
    try: S = Structure(blocks, e)
    finally: rs_core.beta_keys = old
    return S


def equilibrated(fam, P, passes=3, e0=None, tol=1e-10, verbose=False, blocks=None, check_ok=True):
    """P: single parameter vector.  iterate: solve, rescale basis by power-of-two approx of 1/sqrt(diag Gs).  Returns dict (best Solved pass)."""
    bt, ok, gp, W, CT = beta_all(fam, np.asarray(P, float)[:, None]); assert ok[0] or not check_ok
    blocks = blocks or get_blocks(fam); e = [list(x) for x in e0] if e0 else [[0]*len(b) for b in blocks]; hist = []; best = None; o = None
    for it in range(passes + 1):
        S = make_struct(fam, e); bf = beta_rows(S, bt[0])
        o = solve_dual(S, bf, tol=tol); mn, md = gram_stats(S, o); o.update(mineig=mn, maxdiag=md, trel=mn/md if md > 0 else float('nan'), e=[list(x) for x in e], it=it)
        hist.append((o['t'], o['trel'], o['status']))
        if verbose: print('   pass', it, 't=%.3e trel=%.3e %s' % (o['t'], o['trel'], o['status']), flush=True)
        if o['status'] in ('Solved', 'AlmostSolved') and np.isfinite(o['trel']) and (best is None or o['trel'] > best['trel']): best = o
        if it == passes: break
        Gs = Gs_scaled(S, o['y']); dg = [np.maximum(np.diag(G), 1e-10*md) for G in Gs]
        e = [[e[bi][k] - int(round(0.5*math.log2(dg[bi][k]))) for k in range(len(e[bi]))] for bi in range(len(blocks))]
    res = best or o; res['hist'] = hist; res['beta'] = bt[0]
    return res
