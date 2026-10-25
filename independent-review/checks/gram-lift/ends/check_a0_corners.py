"""Two corner blocks of every Gram certificate, and why the a -> 0 end of C2 has margin ~ a^2.

In any Gram certificate F_c = Z^T G Z (any k), the entries on the pairs (uZ, uvY) and (uvX, vZ)
are fixed by single monomials of F_c, so the 2x2 blocks B1, B2 below are the same for every G.
On C2 (lam = beta = 0) their relative determinants are
    rho1 = 4A/(4m + 4alpha + 2d + 4A),   rho2 = (d + 2n)/(2D + d + 2n).
Along the native C2 chart with v_c = a c -> 0, A ~ a^2 and D ~ a^-2, so both are ~ a^2:
no Gram certificate has a normalised margin bounded below at this end, and no diagonal
rescaling changes that.  In the sheared basis (uZ - uvY, uvY) B1 = diag(K, 4A) is positive
for every a > 0, so this is a loss of scale, not of positivity.
Usage: python check_a0_corners.py   (a few seconds)
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import sympy as sp
import gram_lift as gl

u, v = gl.u, gl.v
W, X, Y, Z = gl.PSI
d, lam, n, m, C, A, D, al, be = gl.PAR
P = sp.Poly(sp.expand(gl.form(gl.PAR)), u, v, W, X, Y, Z)
co = P.coeff_monomial
B1 = sp.Matrix([[co(u**2*Z**2), co(u**2*v*Y*Z)/2], [co(u**2*v*Y*Z)/2, co(u**2*v**2*Y**2)]])
B2 = sp.Matrix([[co(u**2*v**2*X**2), co(u*v**2*X*Z)/2], [co(u*v**2*X*Z)/2, co(v**2*Z**2)]])
K = 4*m + 4*al + 2*d
assert sp.expand(B1 - sp.Matrix([[K + 4*be, -K], [-K, K + 4*A]])) == sp.zeros(2)
assert sp.expand(B2 - sp.Matrix([[4*D + 2*d + 4*lam + 4*n, 4*be - 4*D], [4*be - 4*D, 4*D]])) == sp.zeros(2)
assert sp.expand(B1.det() - 8*(A*(2*al + 2*be + d + 2*m) + be*(2*al + d + 2*m))) == 0
assert sp.expand(B2.det() - 8*(D*(d + 2*lam + 2*n + 4*be) - 2*be**2)) == 0
c2 = {lam: 0, be: 0}
rho = lambda B: sp.simplify((B.det()/(B[0, 0]*B[1, 1])).subs(c2))
assert sp.simplify(rho(B1) - 4*A/(K + 4*A)) == 0
assert sp.simplify(rho(B2) - (d + 2*n)/(2*D + d + 2*n)) == 0
print('PASS: B1, B2 and their determinants; on C2 rho1 = 4A/(K+4A), rho2 = (d+2n)/(2D+d+2n)')

# native C2 chart (c2_short.tex:275-292, as in ../../physics/check_c2_to_c1_counterexample.py), d = 1
a = sp.symbols('a', positive=True)
def c2_native(c, b, t, r):
    z = c**2
    h = b - a**2*t; g = (1 + z)*b - (1 - r + a**2*z*t); j = 1 + z*r - t; Ac = 1 + z*r
    Dl = ((1 + z)*z*b**2 + 2*z*(1 + z)*b*r + (1 + z)*(1 + a**2*z)*b + z**2*(1 - a**2)*r**2
          + 2*z*(1 + a**2*z)*r + 1 + a**2*z)
    w = {n: ((b + r)*Ac*(Ac + a**2*z*t) - r*Dl)/(2*r*Dl),
         m: (1 + a**2*z)*(1 - r)*Ac*(Ac + a**2*z*t)/(2*g*Dl),
         A: a**2*Ac*(b + r)*(Ac + a**2*z*t)/(2*h*Dl),
         D: (1 + z)*(b + r)*(b - a**2*z*r)*(Ac + a**2*z*t)/(2*a**2*j*Dl),
         al: z**2*(b + r)*((1 + z)*(b + r)*h + a**2*(g + j + r*t*(1 + a**2*z)))/(2*(1 - z*h)*Dl),
         d: 1}
    contact_vY_over_Z = (a*c)*c*(b/a - a*t)       # v_c Y_c / Z_c at the contact (Z_c = 1)
    return w, contact_vY_over_Z

R = sp.Rational
for name, (c, b, t, r), want in (('P1', (R(2, 5), R(27, 20), R(9, 20), R(1, 2)), (0.19766, 0.43448)),
                                 ('P2', (R(1), R(13, 20), R(2, 5), R(3, 10)), (0.12426, 1.50000)),
                                 ('P3', (R(3, 2), R(77, 200), R(1, 2), R(1, 5)), (0.06264, 2.75225))):
    w, ratio = c2_native(c, b, t, r)
    r1 = sp.limit(rho(B1).subs(w)/a**2, a, 0); r2 = sp.limit(rho(B2).subs(w)/a**2, a, 0)
    A0 = sp.limit(w[A]/a**2, a, 0); D0 = sp.limit(w[D]*a**2, a, 0); nr = sp.limit(ratio, a, 0)
    print('%s  A/a^2 -> %s, a^2 D -> %s,  rho1/(2a^2) -> %.5f, rho2/(2a^2) -> %.5f,  v_cY_c/Z_c -> %s (c^2 b = %s)'
          % (name, A0, D0, float(r1/2), float(r2/2), nr, c**2*b))
    assert A0 > 0 and D0 > 0 and nr == c**2*b
    assert abs(float(r1/2) - want[0]) < 1e-4 and abs(float(r2/2) - want[1]) < 1e-4
print('PASS: at P1-P3 both relative determinants are exactly of order a^2 (limits match the SDP survey)')
