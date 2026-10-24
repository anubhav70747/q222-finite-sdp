"""Exact verification of the separation certificate in cert_data.json.

Usage:  python verify_cert.py PENCIL_JSON [CERT_JSON]
  PENCIL_JSON : ramified_order49_pencil.json (the order-49 affine pencil)
  CERT_JSON   : defaults to cert_data.json next to this script

Everything is checked in exact rational arithmetic (fractions + sympy for the
symbolic determinant).  The source order is (d, lambda, n, m, C, A, D, alpha, beta).

(A) c_sep is valid, i.e. P_c(u,v) is PSD for all real u, v, hence F_c >= 0.
    * det P_c = 16 Du^2 Dv^2 u^2 v^2 L(u^2, v^2) with L = h + sum_{e != 0} q_e U^e1 V^e2
      (exact polynomial division is performed here).
    * L >= h + sum_k m_k > 0 on U, V > 0, via 12-circuit weighted AM-GM bounds
      m_k <= prod (c_e/beta_e)^beta_e, certified by integer-power comparison.
    * P_c is positive definite at one rational point per sign quadrant (exact LDL^T);
      since det P_c never vanishes on an open quadrant the inertia is constant there,
      so P_c > 0 on all four open quadrants and PSD everywhere by closure.
(B) c_sep is not in C_ram + pi(C_ram): rational PSD matrices Y, Y' (49 x 49, pencil
    block structure) annihilate every auxiliary matrix, define one functional ell with
    ell_k = Tr(A_k Y) and ell_{pi(k)} = Tr(A_k Y'), and ell(c_sep) < 0.
"""
import sys, os, json
from fractions import Fraction as Fr
import sympy as sp

PI = [0, 4, 5, 3, 1, 2, 6, 8, 7]   # mirror: swaps lambda<->C, n<->A, alpha<->beta; fixes d, m, D
results = []


def report(name, ok, detail=""):
    results.append(bool(ok))
    print("%-4s %s%s" % ("PASS" if ok else "FAIL", name, ("  [" + detail + "]") if detail else ""))


def ldl_positive_definite(M):
    n = len(M)
    A = [row[:] for row in M]
    for k in range(n):
        if A[k][k] <= 0:
            return False
        for i in range(k + 1, n):
            f = A[i][k] / A[k][k]
            if f:
                for j in range(k + 1, n):
                    A[i][j] -= f * A[k][j]
    return True


# ------------------------------------------------------------------ inputs
here = os.path.dirname(os.path.abspath(__file__))
pencil_path = sys.argv[1]
cert_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "cert_data.json")
pencil = json.load(open(pencil_path))
cert = json.load(open(cert_path))
names = pencil["source_names"]
assert names == cert["source_names"] == ['d', 'lambda', 'n', 'm', 'C', 'A', 'D', 'alpha', 'beta']
csep = {k: Fr(v) for k, v in cert["c_sep"].items()}
cvec = [csep[k] for k in names]
print("c_sep =", {k: str(v) for k, v in csep.items()})
report("c_sep has d > 0 and all event weights >= 0", csep['d'] > 0 and all(csep[k] >= 0 for k in names[1:]))

# ------------------------------------------------------------------ part A
u, v, W, X, Y, Z = sp.symbols('u v W X Y Z')
w = sp.symbols('lam n m C A D alpha beta d')
lam, n, m, C, A, D, al, be, d = w
Du, Dv = 1 + u**2, 1 + v**2
F = (4*lam*Du*(W + v*X)**2 + 4*n*Du*Dv*X**2 + 4*m*(v*W - X + u*v*Y - u*Z)**2
     + 4*C*Dv*(u*W - Y)**2 + 4*A*Du*Dv*Y**2 + (4*D + 2*d)*(u*W + u*v*X - Y - v*Z)**2
     + (4*al + 2*d)*Du*(v*Y - Z)**2 + (4*be + 2*d)*Dv*(X + u*Z)**2 - 2*d*Du*Dv*Z**2)
psi = [W, X, Y, Z]
P = sp.Matrix(4, 4, lambda i, j: sp.Rational(1, 2) * sp.diff(F, psi[i], psi[j]))
det = sp.expand(P.det(method='berkowitz'))
qq, rr = sp.div(sp.Poly(det, u, v, domain=sp.QQ[w]), sp.Poly(sp.expand(16 * Du**2 * Dv**2), u, v, domain=sp.QQ[w]))
report("det P_c divisible by 16 Du^2 Dv^2 (exact)", rr.is_zero)
quot = sp.Poly(qq.as_expr(), u, v)
laurent = {}
shape_ok = True
for (a, b), c in quot.terms():
    if a % 2 or b % 2 or not (0 <= a <= 4 and 0 <= b <= 4):
        shape_ok = False
    laurent[(a // 2 - 1, b // 2 - 1)] = sp.Poly(c, *w)
report("quotient = u^2 v^2 L(u^2,v^2), exponents in {-1,0,1}^2", shape_ok)

# order of generators (lam, n, m, C, A, D, alpha, beta, d) <- source names
order = ['lambda', 'n', 'm', 'C', 'A', 'D', 'alpha', 'beta', 'd']
vals = [csep[k] for k in order]


def eval_poly(poly):
    tot = Fr(0)
    for mon, co in poly.terms():
        t = Fr(int(co.p), int(co.q))
        for x, k in zip(vals, mon):
            if k:
                t *= x**k
        tot += t
    return tot


h = eval_poly(laurent[(0, 0)])
q = {e: eval_poly(p) for e, p in laurent.items() if e != (0, 0)}
report("all 8 off-centre q_e(c_sep) >= 0", len(q) == 8 and all(x >= 0 for x in q.values()),
       "h = %.6e" % float(h))

circuits = cert["partA"]["circuits"]
tot = {e: Fr(0) for e in q}
msum = Fr(0)
circ_ok = True
for k, ci in enumerate(circuits):
    S = [tuple(e) for e in ci["support"]]
    p = ci["p"]
    N = ci["N"]
    cs = [Fr(x) for x in ci["c"]]
    mk = Fr(ci["m"])
    ok = (sum(p) == N and all(pe > 0 for pe in p) and all(e in q for e in S)
          and sum(pe * e[0] for pe, e in zip(p, S)) == 0 and sum(pe * e[1] for pe, e in zip(p, S)) == 0
          and all(c >= 0 for c in cs) and mk >= 0 and len(S) == len(p) == len(cs))
    prod = Fr(1)
    for c, pe in zip(cs, p):
        prod *= (c * N / pe) ** pe
    ok = ok and (mk ** N <= prod)   # m_k <= prod (c_e/beta_e)^beta_e, beta_e = p_e/N
    circ_ok = circ_ok and ok
    for e, c in zip(S, cs):
        tot[e] += c
    msum += mk
report("each circuit balanced, c >= 0, and m_k^N <= prod (N c_e/p_e)^p_e (integer powers)", circ_ok,
       "%d circuits used" % len(circuits))
report("sum_k c^(k)_e <= q_e for all 8 exponents", all(tot[e] <= q[e] for e in q))
report("h + sum_k m_k > 0  (hence L(U,V) > 0 for all U,V > 0)", h + msum > 0,
       "h + sum m = %.6e" % float(h + msum))

pts_ok = True
for su in (1, -1):
    for sv in (1, -1):
        sub = {w[i]: sp.Rational(vals[i].numerator, vals[i].denominator) for i in range(9)}
        sub[u] = sp.Rational(su, 4)
        sub[v] = sp.Rational(4 * sv)
        Pm = P.subs(sub)
        M = [[Fr(int(Pm[i, j].p), int(Pm[i, j].q)) for j in range(4)] for i in range(4)]
        pts_ok = pts_ok and ldl_positive_definite(M)
report("P_c(+-1/4, +-4) positive definite in all four sign quadrants (exact LDL^T)", pts_ok)
report("(A) c_sep is valid: P_c PSD for all real (u,v), so F_c >= 0 and c_sep in C_can",
       all(results))

# ------------------------------------------------------------------ part B
nA = len(pencil["auxiliary_matrices"])
blocks = [(b["offset"], b["order"]) for b in pencil["blocks"] if b["order"] > 1]
assert pencil["order"] == 49 and pencil["source_dimension"] == 9 and nA == 148


def block_of(i):
    for bi, (o, nn) in enumerate(blocks):
        if o <= i < o + nn:
            return bi
    return None


def build_Y(scalars, blks):
    Ymat = {}
    for i, x in enumerate(scalars):
        Ymat[(i, i)] = Fr(x)
    for bi, ent in enumerate(blks):
        o = blocks[bi][0]
        for p, qd, x in ent:
            Ymat[(o + p, o + qd)] = Fr(x)
    return Ymat


def trace(sparse, Ymat):
    """Tr(M Y) for symmetric M given by upper-triangular entries and symmetric Y (upper stored)."""
    t = Fr(0)
    for i, j, x in sparse:
        yv = Ymat.get((i, j), Fr(0))
        t += Fr(x) * yv * (1 if i == j else 2)
    return t


pb = cert["partB"]
Ymat = build_Y(pb["Y_scalar"], pb["Y_blocks"])
Ypmat = build_Y(pb["Yp_scalar"], pb["Yp_blocks"])
ell = [Fr(x) for x in pb["ell"]]


def psd_check(Ym):
    inside = all(block_of(i) is not None and block_of(i) == block_of(j) for (i, j) in Ym if i >= 9 or j >= 9)
    scal = all(Ym.get((i, i), Fr(0)) >= 0 for i in range(9))
    only_diag_scalars = all(i == j for (i, j) in Ym if i < 9 or j < 9)
    pd = True
    for bi, (o, nn) in enumerate(blocks):
        M = [[Ym.get((o + min(p, qd), o + max(p, qd)), Fr(0)) for qd in range(nn)] for p in range(nn)]
        pd = pd and ldl_positive_definite(M)
    return inside and scal and only_diag_scalars and pd


report("Y PSD: block structure, 9 scalar entries >= 0, four blocks positive definite (exact LDL^T)", psd_check(Ymat))
report("Y' PSD: block structure, 9 scalar entries >= 0, four blocks positive definite (exact LDL^T)", psd_check(Ypmat))
aux = pencil["auxiliary_matrices"]
report("Tr(B_j Y) = 0 exactly for all 148 auxiliaries", all(trace(a, Ymat) == 0 for a in aux))
report("Tr(B_j Y') = 0 exactly for all 148 auxiliaries", all(trace(a, Ypmat) == 0 for a in aux))
src = pencil["source_matrices"]
report("ell_k = Tr(A_k Y) for k = 0..8", all(trace(src[k], Ymat) == ell[k] for k in range(9)))
report("ell_{pi(k)} = Tr(A_k Y') for k = 0..8", all(trace(src[k], Ypmat) == ell[PI[k]] for k in range(9)))
val = sum(l * c for l, c in zip(ell, cvec))
report("ell(c_sep) < 0", val < 0, "ell(c_sep) = %.6e" % float(val))
report("(B) ell >= 0 on C_ram and on pi(C_ram), ell(c_sep) < 0: c_sep not in C_ram + pi(C_ram)",
       all(results[-7:]))

# ------------------------------------------------------------------ sanity (informational)
for fn, label in [("ramified_order49_cbar_auxiliaries.json", "cbar"),
                  ("ramified_order49_source_auxiliaries.json", "source")]:
    p = os.path.join(os.path.dirname(os.path.abspath(pencil_path)), fn)
    if os.path.exists(p):
        pt = [Fr(x) for x in json.load(open(p))["source"]]
        a = sum(l * c for l, c in zip(ell, pt))
        b = sum(l * c for l, c in zip(ell, [pt[PI[k]] for k in range(9)]))
        print("INFO ell(%s) = %.4e >= 0: %s ; ell(pi %s) = %.4e >= 0: %s"
              % (label, float(a), a >= 0, label, float(b), b >= 0))

print()
print("OVERALL:", "PASS" if all(results) else "FAIL", "(%d/%d checks)" % (sum(results), len(results)))
