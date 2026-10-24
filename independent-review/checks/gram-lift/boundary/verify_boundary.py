"""Exact verification of boundary (extreme) source certificates in the sextic Puiseux-Gram lift L_6.

Usage:  python verify_boundary.py [DATA.json ...]      (default: every *.json in ./data next to this script)

Source order (d, lambda, n, m, C, A, D, alpha, beta); parent
  F_c = 4 lam Du (W+vX)^2 + 4 n Du Dv X^2 + 4 m (vW - X + uvY - uZ)^2 + 4 C Dv (uW - Y)^2 + 4 A Du Dv Y^2
        + (4D+2d)(uW + uvX - Y - vZ)^2 + (4 alpha + 2d) Du (vY - Z)^2 + (4 beta + 2d) Dv (X + uZ)^2 - 2d Du Dv Z^2.
Everything is exact (fractions + sympy); no floating point enters any decision.
Checks per source
 (1) rational data: native-chart parameters reproduce the stored source (C1 or C2 native formulas), d = 1,
     all event weights >= 0, chart positivity, u_c and v_c are sixth powers of rationals.
 (2) exact zero and zero gradient of F_c at the contact.
 (3) F_c >= 0 (determinant / weighted AM-GM route): det P / (16 Du^2 Dv^2) = U V L(U,V), L = h + sum_e q_e U^e1 V^e2,
     q_e >= 0, L(U_c,V_c) = 0 and grad_log L = 0 at the contact => L >= 0 on U,V > 0 by weighted AM-GM (rational weights);
     the weight support spans R^2 (so the zero is isolated) and P is positive definite at rational points of all four
     open quadrants (exact LDL^T) => P >= 0 on R^2.
 (4) F_c(s^6, t^6, psi) = Z^T (D^-1 V G' V^T D^-1) Z as polynomials in s, t, psi (s = u^(1/6), t = v^(1/6)),
     block-diagonal in the parity classes of the Puiseux exponents, V = pivot form of the complement of the contact vector.
 (5) every G'_b is positive definite (exact LDL^T): G = D^-1 V G' V^T D^-1 is PSD with rank sum(n_b - 1)."""
import sys, os, json, glob
from fractions import Fraction as Fr
import sympy as sp

u, v, W, X, Y, Z = sp.symbols('u v W X Y Z'); PSI = [W, X, Y, Z]
results = []


def report(name, ok, detail=""):
    results.append(bool(ok)); print("  %-4s %s%s" % ("PASS" if ok else "FAIL", name, ("  [" + detail + "]") if detail else ""))
    return bool(ok)


def sr(x): x = Fr(x); return sp.Rational(x.numerator, x.denominator)


def iroot6(n):
    lo, hi = 0, 1
    while hi ** 6 <= n: hi *= 2
    while lo < hi - 1:
        mid = (lo + hi) // 2
        if mid ** 6 <= n: lo = mid
        else: hi = mid
    return lo if lo ** 6 == n else None


def sixth_root(x):
    p, q = iroot6(x.numerator), iroot6(x.denominator)
    return None if p is None or q is None else Fr(p, q)


def parent(c):
    dd, lam, n, m, C, A, D, al, be = c; Du, Dv = 1 + u**2, 1 + v**2
    return (4*lam*Du*(W + v*X)**2 + 4*n*Du*Dv*X**2 + 4*m*(v*W - X + u*v*Y - u*Z)**2 + 4*C*Dv*(u*W - Y)**2
            + 4*A*Du*Dv*Y**2 + (4*D + 2*dd)*(u*W + u*v*X - Y - v*Z)**2 + (4*al + 2*dd)*Du*(v*Y - Z)**2
            + (4*be + 2*dd)*Dv*(X + u*Z)**2 - 2*dd*Du*Dv*Z**2)


# ---------------- native formulas (appendix C2 source; C1 native source SD1 normalised to d = 1) ----------------
def c2_source(c, a, b, t, r):
    z = c*c; Ac = 1 + z*r; h = b - a*a*t; g = (1 + z)*b - (1 - r + a*a*z*t); j = 1 + z*r - t
    Dl = ((1 + z)*z*b*b + 2*z*(1 + z)*b*r + (1 + z)*(1 + a*a*z)*b + z*z*(1 - a*a)*r*r + 2*z*(1 + a*a*z)*r + 1 + a*a*z)
    N = (1 + z)*(b + r)*(b + a*a*(1 - t)) - a*a*(1 + a*a*z)*(1 - r)*t
    src = dict(d=Fr(1), lam=Fr(0), n=((b + r)*Ac*(Ac + a*a*z*t) - r*Dl)/(2*r*Dl),
               m=(1 + a*a*z)*(1 - r)*Ac*(Ac + a*a*z*t)/(2*g*Dl), C=Ac*N/(2*a*a*t*Dl),
               A=a*a*Ac*(b + r)*(Ac + a*a*z*t)/(2*h*Dl), D=(1 + z)*(b + r)*(b - a*a*z*r)*(Ac + a*a*z*t)/(2*a*a*j*Dl),
               alpha=z*z*(b + r)*N/(2*(1 - z*h)*Dl), beta=Fr(0))
    chart = [a > 0, c > 0, b > 0, t > 0, 0 < r < 1, h > 0, g > 0, j > 0, 1 - z*h > 0, b - a*a*z*r > 0]
    return src, (c, a*c, b/a, -c*r, c*(b/a - a*t), Fr(1)), chart


def c1_source(k, f, g, a, t):
    kap, eta, zeta = k*k, f*f, g*g; V = kap*eta*zeta; Lc = 1 - a*eta; Hc = 1 + (1 - a)*eta; Qc = 1 - kap*zeta + a*V; K = 1 + kap
    e = k*f/g; Dl = 1/e; R = 1 + kap*(1 - a*eta); T = 1 + kap - kap*zeta*t; M = R - t; J = T - a
    P4 = 1 + kap + eta*kap*(a - kap - 1) - zeta*kap*t + V*(-zeta*kap*t**2 + 2*kap*t + a*(1 - t) + t)
    Pal = 1 + kap - a*eta*kap + zeta*((1 + kap) - (2 + kap)*t) + V*(-zeta*t**2 + a*(t - 1) + t)
    Pd = (1 + kap)**2 - a*eta*kap**2 + zeta*kap*((kap + 1) - (kap + 2)*t) + V*(-zeta*kap*t**2 + kap*a*(t - 1) + kap*t + kap + 1)
    c2 = 2*K*(1 - t*zeta)*Lc*(1 + V)*Hc*Qc; cb = 2*K*Lc*(1 + (1 - t)*zeta)*Hc*Qc; c3 = 2*K*Lc*Hc*Qc*(1 - kap*eta + t*V)
    c5 = 2*K*Lc*Hc*Qc*(1 + V); c4 = 2*Qc*P4; cal = 2*Lc*Hc*Pal; cd = 4*Lc*Qc*Pd
    dB = 2*a*eta*cd/c3
    src = dict(d=Fr(1), lam=Fr(0), n=Fr(0), m=2*a*e*e*c2/(M*c3)/dB, C=2/dB, A=2*a*eta**2*c4/c3/dB,
               D=2*a*c5/(J*c3)/dB - Fr(1, 2), alpha=2*a*Dl**2*k**6*eta**2*cal/c3/dB, beta=2*a*e*e*cb/(t*c3)/dB - Fr(1, 2))
    chart = [k > 0, f > 0, g > 0, a > 0, t > 0, k <= 1, Lc > 0, Hc > 0, Qc > 0, M > 0, J > 0, c3 > 0, dB > 0]
    return src, (k, Dl*k*k*eta, Dl*k, -k*(1 - t*zeta), Dl*k*k*Lc, Fr(1)), chart


def ldl_pd(M):
    n = len(M); A = [row[:] for row in M]
    for k in range(n):
        if A[k][k] <= 0: return False
        for i in range(k + 1, n):
            f = A[i][k]/A[k][k]
            if f:
                for j in range(k + 1, n): A[i][j] -= f*A[k][j]
    return True


def verify(path):
    D_ = json.load(open(path)); print("== %s  (%s source)" % (D_['name'], D_['family'])); start = len(results)
    S = D_['source']; names = ['d', 'lam', 'n', 'm', 'C', 'A', 'D', 'alpha', 'beta']
    src = {k: Fr(S[k]) for k in names}; cvec = [src[k] for k in names]
    P_ = {k: Fr(x) for k, x in D_['params'].items()}
    if D_['family'] == 'C2': s2, ct2, chart = c2_source(P_['c'], P_['a'], P_['b'], P_['t'], P_['r'])
    else: s2, ct2, chart = c1_source(P_['k'], P_['f'], P_['g'], P_['a'], P_['t'])
    contact = tuple(Fr(x) for x in D_['contact'])
    # (1) rational data
    report("(1a) stored source = native-formula source for the stored parameters; stored contact = native contact",
           all(s2[k] == src[k] for k in names) and tuple(ct2) == contact)
    report("(1b) native chart conditions hold", all(chart))
    report("(1c) d = 1, lambda/n/beta zero pattern, all event weights rational and >= 0",
           src['d'] == 1 and all(src[k] >= 0 for k in names[1:]) and ((src['lam'] == 0 and src['beta'] == 0) if D_['family'] == 'C2' else (src['lam'] == 0 and src['n'] == 0)))
    sc, tc = sixth_root(contact[0]), sixth_root(contact[1])
    report("(1d) u_c and v_c are sixth powers of rationals", sc is not None and tc is not None and sc**6 == contact[0] and tc**6 == contact[1],
           "u_c^(1/6) = %s, v_c^(1/6) = %s" % (sc, tc))
    # (2) contact
    F = parent([sr(x) for x in cvec]); sub = dict(zip((u, v, W, X, Y, Z), [sr(x) for x in contact]))
    val = sp.expand(F).subs(sub); gr = [sp.diff(F, w).subs(sub) for w in (u, v, W, X, Y, Z)]
    report("(2) F_c(contact) = 0 and grad F_c(contact) = 0 exactly", val == 0 and all(g == 0 for g in gr))
    # (3) determinant / AM-GM
    Pm = sp.Matrix(4, 4, lambda i, j: sp.Rational(1, 2)*sp.diff(sp.expand(F), PSI[i], PSI[j]))
    det = sp.expand(Pm.det(method='berkowitz')); Du, Dv = 1 + u**2, 1 + v**2
    qq, rr = sp.div(sp.Poly(det, u, v), sp.Poly(sp.expand(16*Du**2*Dv**2), u, v))
    shape = rr.is_zero and all(a % 2 == 0 and b % 2 == 0 and 0 <= a <= 4 and 0 <= b <= 4 for (a, b), _ in qq.terms())
    report("(3a) det P = 16 Du^2 Dv^2 U V L(U,V) with L a Laurent polynomial in U = u^2, V = v^2, exponents in {-1,0,1}^2", shape)
    co = {(a//2 - 1, b//2 - 1): Fr(int(c.p), int(c.q)) for (a, b), c in qq.terms()}
    h = co.get((0, 0), Fr(0)); q = {(e1, e2): co.get((e1, e2), Fr(0)) for e1 in (-1, 0, 1) for e2 in (-1, 0, 1) if (e1, e2) != (0, 0)}
    Uc, Vc = contact[0]**2, contact[1]**2
    report("(3b) all 8 off-centre q_e >= 0 and h < 0", all(x >= 0 for x in q.values()) and h < 0, "h = %.4e" % float(h))
    mono = lambda e: Uc**e[0]*Vc**e[1]
    L0 = h + sum(x*mono(e) for e, x in q.items())
    g1 = sum(x*mono(e)*e[0] for e, x in q.items()); g2 = sum(x*mono(e)*e[1] for e, x in q.items())
    report("(3c) L(U_c,V_c) = 0 and U dL/dU = V dL/dV = 0 at the contact (exact)", L0 == 0 and g1 == 0 and g2 == 0)
    beta = {e: x*mono(e)/(-h) for e, x in q.items() if x > 0}
    supp = list(beta); span = len(supp) >= 3 and any(e1[0]*e2[1] - e1[1]*e2[0] != 0 for e1 in supp for e2 in supp)
    report("(3d) AM-GM weights beta_e = q_e x_c^e/|h| >= 0 sum to 1 with sum beta_e e = 0, support spans R^2 (L >= 0, zero only at contact)",
           sum(beta.values()) == 1 and sum(b*e[0] for e, b in beta.items()) == 0 and sum(b*e[1] for e, b in beta.items()) == 0 and span,
           "support size %d" % len(supp))
    pts = [Fr(1, 2), Fr(2), Fr(1), Fr(1, 4), Fr(4), Fr(3), Fr(1, 3)]; quad_ok = True
    for su in (1, -1):
        for sv in (1, -1):
            hit = False
            for pu in pts:
                for pv in pts:
                    M = Pm.subs({u: sr(su*pu), v: sr(sv*pv)}); Mf = [[Fr(int(M[i, j].p), int(M[i, j].q)) for j in range(4)] for i in range(4)]
                    if ldl_pd(Mf): hit = True; break
                if hit: break
            quad_ok &= hit
    report("(3e) P(u,v) positive definite at a rational point of each open quadrant (exact LDL^T)  => F_c >= 0 on R^2", quad_ok)
    # (4) Gram identity
    blocks = [[tuple(e) for e in b] for b in D_['basis_blocks']]; i0s = D_['pivot_index']
    Gp = [[[Fr(x) for x in row] for row in M] for M in D_['gram_reduced']]
    zc = [[sc**i * tc**j * contact[2 + l] for (i, j, l) in blk] for blk in blocks]
    report("(4a) contact basis vector z_c is nonzero entrywise (so the scaling z/z_c is exact) and shapes are consistent",
           all(all(x != 0 for x in z) for z in zc) and all(len(Gp[b]) == len(blocks[b]) - 1 for b in range(len(blocks))))
    rhs = {}
    for bi, blk in enumerate(blocks):
        n = len(blk); i0 = i0s[bi]; others = [r for r in range(n) if r != i0]; a_ = [1/x for x in zc[bi]]
        def key(r1, r2):
            (i1, j1, l1), (i2, j2, l2) = blk[r1], blk[r2]; return (i1 + i2, j1 + j2, min(l1, l2), max(l1, l2))
        for p, rp in enumerate(others):
            for q_, rq in enumerate(others):
                g = Gp[bi][p][q_]
                if g == 0: continue
                # (z_rp/zc_rp - z_i0/zc_i0)(z_rq/zc_rq - z_i0/zc_i0)
                for (r1, c1) in ((rp, a_[rp]), (i0, -a_[i0])):
                    for (r2, c2) in ((rq, a_[rq]), (i0, -a_[i0])):
                        k = key(r1, r2); rhs[k] = rhs.get(k, Fr(0)) + g*c1*c2
    lhs = {}
    Fp = sp.Poly(sp.expand(F), u, v, *PSI)
    for mon, c in Fp.terms():
        a, b = mon[0], mon[1]; ex = mon[2:]; idxs = [i for i in range(4) for _ in range(ex[i])]
        lhs[(6*a, 6*b, idxs[0], idxs[1])] = lhs.get((6*a, 6*b, idxs[0], idxs[1]), Fr(0)) + Fr(int(c.p), int(c.q))
    rhs = {k: x for k, x in rhs.items() if x != 0}; lhs = {k: x for k, x in lhs.items() if x != 0}
    report("(4b) F_c(s^6,t^6,psi) = sum_b p_b^T G'_b p_b as Puiseux polynomials (exact; %d monomials), p_i = z_i/zc_i - z_i0/zc_i0" % len(lhs), lhs == rhs)
    parity_ok = all(len({(i % 2, j % 2) for (i, j, l) in blk}) == 1 for blk in blocks) and len({(blk[0][0] % 2, blk[0][1] % 2) for blk in blocks}) == len(blocks)
    report("(4c) blocks are the 4 distinct parity classes (i mod 2, j mod 2) of the Puiseux exponents; Gram is zero across blocks and on pruned monomials", parity_ok,
           "block orders %s of 196 monomials" % [len(b) for b in blocks])
    # (5) positive definiteness
    pd = [ldl_pd(M) for M in Gp]
    report("(5) every reduced Gram block G'_b is positive definite (exact LDL^T)", all(pd),
           "rank G = %d, kernel dim within blocks = %d (the contact vectors), %d monomials pruned" % (sum(len(b) - 1 for b in blocks), len(blocks), 196 - sum(len(b) for b in blocks)))
    ok = all(results[start:]); print("  ==> %s %s" % (D_['name'], "PASS" if ok else "FAIL")); return ok


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    files = sys.argv[1:] or sorted(glob.glob(os.path.join(here, 'data', '*.json')))
    allok = [verify(f) for f in files]
    print("\nOVERALL:", "PASS" if all(allok) else "FAIL", "(%d/%d sources, %d/%d checks)" % (sum(allok), len(allok), sum(results), len(results)))
