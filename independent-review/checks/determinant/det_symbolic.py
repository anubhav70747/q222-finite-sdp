"""Symbolic check of Proof.tex eq. (insight-determinant) / 04_geometry.tex SG5-SG6.

det P_c = 16 D_u^2 D_v^2 u^2 v^2 L_c(u^2, v^2), where L_c is a Laurent polynomial
with exponents in {-1,0,1}^2 and every off-centre coefficient is a polynomial in the
event weights with nonnegative coefficients. Only the centre h_c can be negative.
"""
import sympy as sp

u, v, W, X, Y, Z = sp.symbols('u v W X Y Z')
w = sp.symbols('lam n m C A D alpha beta d')
lam, n, m, C, A, D, al, be, d = w
Du, Dv = 1 + u**2, 1 + v**2
F = (4*lam*Du*(W + v*X)**2 + 4*n*Du*Dv*X**2 + 4*m*(v*W - X + u*v*Y - u*Z)**2
     + 4*C*Dv*(u*W - Y)**2 + 4*A*Du*Dv*Y**2 + (4*D + 2*d)*(u*W + u*v*X - Y - v*Z)**2
     + (4*al + 2*d)*Du*(v*Y - Z)**2 + (4*be + 2*d)*Dv*(X + u*Z)**2 - 2*d*Du*Dv*Z**2)
psi = [W, X, Y, Z]
P = sp.Matrix(4, 4, lambda i, j: sp.Rational(1, 2)*sp.diff(F, psi[i], psi[j]))
det = sp.expand(P.det(method='berkowitz'))

# Divide by 16 Du^2 Dv^2 (univariate-in-each-variable, so the quotient is exact).
qq, rr = sp.div(sp.Poly(det, u, v, domain=sp.QQ[w]),
                sp.Poly(sp.expand(16*Du**2*Dv**2), u, v, domain=sp.QQ[w]))
print("det divisible by 16 Du^2 Dv^2:", rr.is_zero)

# The quotient is u^2 v^2 * L(u^2, v^2): only even powers, each between 0 and 4.
quot = sp.Poly(qq.as_expr(), u, v)
ok = True
rows = []
for (a, b), c in sorted(quot.terms()):
    if a % 2 or b % 2 or not (0 <= a <= 4) or not (0 <= b <= 4):
        ok = False
    e = (a//2 - 1, b//2 - 1)
    nonneg = all(co >= 0 for co in sp.Poly(c, *w).coeffs())
    rows.append((e, nonneg, sp.factor(c)))
print("quotient = u^2 v^2 * Laurent polynomial in U=u^2, V=v^2 on {-1,0,1}^2:", ok)
for e, nonneg, c in rows:
    tag = "h_c (centre)" if e == (0, 0) else f"q_{e}  nonnegative coefficients: {nonneg}"
    print(f"  {tag}: {c}")
print("all 8 off-centre q_e have only nonnegative coefficients in the weights:",
      all(nn for e, nn, _ in rows if e != (0, 0)))
