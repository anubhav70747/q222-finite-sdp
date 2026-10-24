#!/usr/bin/env python
"""Sextic Puiseux-Gram lift test for the parent form F_c (psi = (W,X,Y,Z), D_u = 1+u^2, D_v = 1+v^2).

For u,v > 0 and basis Z = {u^(i/k) v^(j/k) psi_l : 0 <= i,j <= k}, d_k(e) = max d such that
F_c = Z^T G Z with G PSD (coefficients of u^(p/k) v^(q/k) psi_l psi_m are matched; those with p or q not a
multiple of k must vanish). The sign symmetries (u,Y,Z)->(-u,-Y,-Z), (v,X,Z)->(-v,-X,-Z) (see symcheck.py)
cover the other quadrants and let G be split into 4 blocks (exact). Zero events allow Newton-polytope pruning.

Prints d* (inlined validity oracle: det P / (16 Du^2 Dv^2) = h + sum_e q_e U^e1 V^e2, U=u^2, V=v^2),
d_k, d_k/d*, and the delta-margin t at delta = 1e-3: max t s.t. the Gram identity holds at d = (1-delta) d*
and G - t I >= 0, in a basis equilibrated by the diagonal of the max-d Gram. t > 0 (scaling like delta) means
strictly feasible at (1-delta) d*; t < 0 means a real gap.
NOTE: Clarabel can be noise-limited on face (zero events) or multiscale cases (ratios far from 1, even > 1,
and tiny |t| ~ 1e-9). The delta-margin (--robust: iterated equilibration, delta = 1e-2,1e-3,1e-4, t/delta
constant) is the robust test; confirm borderline cases with a high-precision solver.
Usage: python gram_lift.py --k 2 --e "1,123/100000,..." | --case random8 [--robust] [--tol 1e-9]"""
import argparse, json, math, os, time
from fractions import Fraction
import numpy as np, scipy.sparse as sps, sympy as sp, mpmath as mp
from scipy.optimize import linprog
import clarabel

u, v, d = sp.symbols('u v d'); W, X, Y, Z = sp.symbols('W X Y Z')
PSI = [W, X, Y, Z]; PAR = list(sp.symbols('d lam n m C A D alpha beta'))   # c = (d, lam, n, m, C, A, D, alpha, beta)

def form(c):
    dd, lam, n, m, C, A, D, al, be = c; Du, Dv = 1 + u**2, 1 + v**2
    return (4*lam*Du*(W + v*X)**2 + 4*n*Du*Dv*X**2 + 4*m*(v*W - X + u*v*Y - u*Z)**2 + 4*C*Dv*(u*W - Y)**2
            + 4*A*Du*Dv*Y**2 + (4*D + 2*dd)*(u*W + u*v*X - Y - v*Z)**2 + (4*al + 2*dd)*Du*(v*Y - Z)**2
            + (4*be + 2*dd)*Dv*(X + u*Z)**2 - 2*dd*Du*Dv*Z**2)
def pmat(c): F = form(c); return sp.Matrix(4, 4, lambda i, j: sp.Rational(1, 2)*sp.diff(F, PSI[i], PSI[j]))

def tensors():     # T[param, a, b, l, m] = coefficient of u^a v^b in P_lm (linear in the 9 parameters)
    P = pmat(PAR); T = np.zeros((9, 3, 3, 4, 4))
    for i in range(4):
        for j in range(4):
            for mon, co in sp.Poly(sp.expand(P[i, j]), u, v, *PAR).terms():
                T[list(mon[2:]).index(1), mon[0], mon[1], i, j] += float(co)
    return T

# ---------------- validity oracle (exact rational e/s, d symbolic) ----------------
def oracle(e_s):
    P = pmat([d] + [sp.Rational(x.numerator, x.denominator) for x in e_s])
    Q, r = sp.div(sp.Poly(P.det(method='berkowitz'), u, v), sp.Poly(sp.expand(16*(1 + u**2)**2*(1 + v**2)**2), u, v))
    assert r.is_zero
    co = {(a//2, b//2): sp.Poly(c, d) for (a, b), c in Q.terms()}      # coefficient of U^alpha V^beta
    E = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0)]
    hp = co.get((1, 1), sp.Poly(0, d)); qp = [co.get((a + 1, b + 1), sp.Poly(0, d)) for a, b in E]
    mpv = lambda t: mp.mpf(t) if isinstance(t, float) else mp.mpf(int(t.p))/int(t.q)
    pev = lambda p, x: sum((mpv(c)*mpv(x)**k for (k,), c in p.terms()), mp.mpf(0))
    mp.mp.dps = 60
    def margin(dd):
        q = np.array([float(pev(p, dd)) for p in qp]); lf, z = phi_log(q, np.array(E, float))
        return pev(hp, dd) + (0 if lf == -np.inf else mp.e**mp.mpf(lf)), z
    return margin

def phi_log(q, E):   # min_z log sum_e q_e exp(e.z)  (damped Newton); -inf if q == 0
    idx = np.where(q > 0)[0]
    if len(idx) == 0: return -np.inf, np.zeros(2)
    lq = np.log(q[idx]); Ee = E[idx]; z = np.zeros(2)
    def fg(z):
        a = lq + Ee @ z; mx = a.max(); w = np.exp(a - mx); s = w.sum(); p = w/s; g = Ee.T @ p
        return mx + math.log(s), g, (Ee.T*p) @ Ee - np.outer(g, g)
    f, g, H = fg(z)
    for _ in range(200):
        mu = 1e-12*max(1.0, np.trace(H)) + 1e-14
        dz = -np.linalg.solve(H + mu*np.eye(2), g)
        if np.linalg.norm(dz) > 20: dz *= 20/np.linalg.norm(dz)
        t = 1.0
        while t > 1e-12:
            f2, g2, H2 = fg(z + t*dz)
            if f2 <= f + 1e-4*t*(g @ dz): break
            t *= 0.5
        else: break
        z = z + t*dz; df = f - f2; f, g, H = f2, g2, H2
        if np.linalg.norm(g) < 1e-13 or (df < 1e-16*max(1, abs(f)) and np.linalg.norm(g) < 1e-9) or np.max(np.abs(z)) > 400: break
    return f, z

def dstar(margin):    # largest d with margin >= 0 (bisection, rtol 1e-13); e already scaled to max 1
    ok = lambda dd: margin(dd)[0] >= 0
    lo, hi = 0.0, 1.0
    if ok(hi):
        while ok(hi): lo, hi = hi, hi*2
    else:
        while not ok(hi): hi /= 2
        lo, hi = hi, hi*2
    while hi - lo > 1e-13*hi:
        mid = 0.5*(lo + hi); lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    dd = 0.5*(lo + hi); return dd, margin(dd)[1]

# ---------------- Gram SDP ----------------
def in_hull(pts, t):
    pts = np.array(pts, float); A = np.vstack([pts.T, np.ones(len(pts))])
    return linprog(np.zeros(len(pts)), A_eq=A, b_eq=[t[0], t[1], 1.0], bounds=(0, None), method='highs').status == 0

def gram(k, T, e, lam, psis, fixed_d=None, bscale=None, tol=1e-9):
    """max d (fixed_d None) or max t with G - tI >= 0 at d=fixed_d. e scaled to max 1; basis rescaled u=lam_u u', v=lam_v v'."""
    act = np.concatenate([[True], e != 0]); sets = []
    for l in range(4):     # Newton-polytope pruning per psi_l from the support of P_ll (d free)
        S = [(a, b) for a in range(3) for b in range(3) if np.any(T[act, a, b, l, l] != 0)]
        sets.append([(i, j) for i in range(k + 1) for j in range(k + 1) if in_hull(S, (2*i/k, 2*j/k))])
    sc = np.array([[lam[0]**a*lam[1]**b for b in range(3)] for a in range(3)])
    f = sc[:, :, None, None]*psis[None, None, :, None]*psis[None, None, None, :]
    P0 = np.einsum('p,pabij->abij', np.concatenate([[0.0], e]), T)*f; P1 = T[0]*f
    cls = {}   # sign-symmetry blocks: parity of i (+ u-odd psi if k odd), same for j
    for l in range(4):
        for (i, j) in sets[l]:
            key = ((i + (l >= 2)*(k % 2)) % 2, (j + (l in (1, 3))*(k % 2)) % 2); cls.setdefault(key, []).append((i, j, l))
    blocks = [cls[c] for c in sorted(cls)]; sizes = [len(b) for b in blocks]
    offs = np.cumsum([0] + [n*(n + 1)//2 for n in sizes]); ng = offs[-1]; rows = {}; ent = []
    row = lambda key: rows.setdefault(key, len(rows))
    for bi, blk in enumerate(blocks):
        for r in range(len(blk)):
            for c in range(len(blk)):
                (ir, jr, lr), (ic, jc, lc) = blk[r], blk[c]
                w = (1.0 if r == c else 2**-0.5)*(2 if lr == lc else 1)*(1.0 if bscale is None else bscale[bi][r]*bscale[bi][c])
                ent.append((row((ir + ic, jr + jc, min(lr, lc), max(lr, lc))), 1 + offs[bi] + max(r, c)*(max(r, c) + 1)//2 + min(r, c), w))
    nz = np.any(T != 0, axis=0)
    for a, b, l, m in zip(*np.where(nz)):
        if l <= m: row((k*a, k*b, l, m))
    nr = len(rows); r0 = np.zeros(nr); r1 = np.zeros(nr)
    for (p, q, l, m), r in rows.items():
        if p % k == 0 and q % k == 0 and p//k <= 2 and q//k <= 2: r0[r] = 2*P0[p//k, q//k, l, m]; r1[r] = 2*P1[p//k, q//k, l, m]
    if fixed_d is not None: r0 = r0 + fixed_d*r1; r1 = 0*r1
    Aeq = sps.coo_matrix(([x[2] for x in ent], ([x[0] for x in ent], [x[1] for x in ent])), shape=(nr, 1 + ng)).tocsc() \
        + sps.csc_matrix((-r1, (np.arange(nr), np.zeros(nr, int))), shape=(nr, 1 + ng))
    iv = np.zeros(ng)
    if fixed_d is not None:
        for bi, n in enumerate(sizes): iv[offs[bi] + np.arange(n)*(np.arange(n) + 1)//2 + np.arange(n)] = 1.0
    A = sps.vstack([Aeq, sps.hstack([sps.csc_matrix(iv.reshape(-1, 1)), -sps.identity(ng)])]).tocsc()
    q = np.zeros(1 + ng); q[0] = -1.0; st = clarabel.DefaultSettings(); st.verbose = False; st.max_iter = 300
    st.tol_gap_abs = st.tol_gap_rel = st.tol_feas = tol
    sol = clarabel.DefaultSolver(sps.csc_matrix((1 + ng, 1 + ng)), q, A, np.concatenate([r0, np.zeros(ng)]),
                                 [clarabel.ZeroConeT(nr)] + [clarabel.PSDTriangleConeT(n) for n in sizes], st).solve()
    x = np.array(sol.x)
    diag = [x[1 + offs[bi] + np.arange(n)*(np.arange(n) + 1)//2 + np.arange(n)] for bi, n in enumerate(sizes)]
    return dict(val=x[0], status=str(sol.status), iters=sol.iterations, sizes=sizes, nbasis=[len(s) for s in sets], diag=diag)

def equilibrate(diag, old=None):
    g = [np.maximum(np.abs(x), 1e-10*np.abs(np.concatenate(diag)).max()) for x in diag]; gm = np.sqrt(max(x.max() for x in g))
    new = [np.sqrt(x)/gm for x in g]; return new if old is None else [a*b for a, b in zip(old, new)]

def parse_e(items): return [Fraction(str(x).strip()) if isinstance(x, str) else Fraction(float(x)) for x in items]

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--k', type=int, default=2); ap.add_argument('--e'); ap.add_argument('--case')
    ap.add_argument('--cases', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cases.json'))
    ap.add_argument('--tol', type=float, default=1e-9); ap.add_argument('--robust', action='store_true'); a = ap.parse_args()
    ev = parse_e(a.e.split(',') if a.e else json.load(open(a.cases))[a.case]['e']); assert len(ev) == 8 and min(ev) >= 0
    s = max(ev); e_s = [x/s for x in ev]; T = tensors(); t0 = time.time()
    dsn, z = dstar(oracle(e_s)); ds = float(s)*dsn; lam = np.exp(z/2)
    e = np.array([float(x) for x in e_s]); c = np.concatenate([[dsn], e])
    Pc = np.einsum('p,pabij->abij', c, T)
    psis = 1/np.sqrt(np.maximum([np.sum(np.abs(Pc[:, :, l, l])*np.array([[lam[0]**i*lam[1]**j for j in range(3)] for i in range(3)])) for l in range(4)], 1e-300))
    print('e/max(e) =', np.array2string(e, precision=4), '\nd* = %.9e   contact (lam_u, lam_v) = (%.4g, %.4g)   [%.1fs]' % (ds, lam[0], lam[1], time.time() - t0))
    r = gram(a.k, T, e, lam, psis, tol=a.tol); dk = r['val']*float(s)
    print('k=%d  Gram blocks %s (order %d)  basis per psi %s' % (a.k, r['sizes'], sum(r['sizes']), r['nbasis']))
    print('d_k = %.9e   d_k/d* = %.9f   [%s, %d iters]' % (dk, dk/ds, r['status'], r['iters']))
    def margin(dl, bs): m = gram(a.k, T, e, lam, psis, fixed_d=(1 - dl)*dsn, bscale=bs, tol=1e-10); return m['val'], m['status']
    if not a.robust:
        t, stt = margin(1e-3, equilibrate(r['diag']))
        print('delta-margin t(1e-3) = %.3e [%s]   (t>0: feasible at 0.999 d*; t<0: gap)' % (t, stt))
    else:
        bs = None
        for _ in range(2): bs = equilibrate(gram(a.k, T, e, lam, psis, fixed_d=0.99*dsn, bscale=bs)['diag'], bs)
        for dl in (1e-2, 1e-3, 1e-4):
            t, stt = margin(dl, bs); print('delta=%g  t = %.3e  t/delta = %.3e [%s]' % (dl, t, t/dl, stt))
