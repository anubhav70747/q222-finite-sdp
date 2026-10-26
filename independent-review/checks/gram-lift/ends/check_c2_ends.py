"""Exact facts behind the r and c ends of the native C2 chart (c, a, b, t, r), d = 1, lam = beta = 0.

1. r -> 1 is not an end of the weight-nonnegative chart: n(r=1) = -c^2 (h + a^2) / (2(c^2(a^2+b+1)+1)) < 0
   whenever h = b - a^2 t > 0, so n reaches 0 at some r_n < 1 first (an n = 0 wall, not r = 1).
2. r -> 0: n = n0/r + O(1) with n0 > 0, and m, C, A, D, alpha have finite limits.
3. c -> 0 and c -> infinity: along the torus orbit O(L) = (L u_c, L v_c, W_c/L, X_c/L^2, Y_c/L^2, Z_c/L)
   through the contact, every C2 parent satisfies
       F_c(O(L)) = 4C (u_c v_c W_c)^2 (L - 1/L)^2 + [(4n+4m+2) X_c^2 + (4A+4C+4D+2) Y_c^2] (1 - 1/L^2)^2,
   a positive first-order lift of the orbit that the c-end limit forms vanish on.
All three are identities of rational functions in (c, a, b, t, r) (and L).
Usage: python check_c2_ends.py   (about a minute)
"""
import os, sys, random
from fractions import Fraction as Fr
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(HERE, '..', 'boundary'))
import gram_lift as gl
import verify_boundary as vb

c, a, b, t, r, L = sp.symbols('c a b t r L', positive=True)
u, v = gl.u, gl.v
W, X, Y, Z = gl.PSI
d, lam, n, m, C, A, D, al, be = gl.PAR


def c2_weights(c, a, b, t, r):
    z = c*c; Ac = 1 + z*r; h = b - a*a*t; g = (1 + z)*b - (1 - r + a*a*z*t); j = 1 + z*r - t
    Dl = ((1 + z)*z*b*b + 2*z*(1 + z)*b*r + (1 + z)*(1 + a*a*z)*b + z*z*(1 - a*a)*r*r + 2*z*(1 + a*a*z)*r + 1 + a*a*z)
    N = (1 + z)*(b + r)*(b + a*a*(1 - t)) - a*a*(1 + a*a*z)*(1 - r)*t
    return {d: 1, lam: 0, be: 0, n: ((b + r)*Ac*(Ac + a*a*z*t) - r*Dl)/(2*r*Dl),
            m: (1 + a*a*z)*(1 - r)*Ac*(Ac + a*a*z*t)/(2*g*Dl), C: Ac*N/(2*a*a*t*Dl),
            A: a*a*Ac*(b + r)*(Ac + a*a*z*t)/(2*h*Dl), D: (1 + z)*(b + r)*(b - a*a*z*r)*(Ac + a*a*z*t)/(2*a*a*j*Dl),
            al: z*z*(b + r)*N/(2*(1 - z*h)*Dl)}


Wt = c2_weights(c, a, b, t, r)
contact = (c, a*c, b/a, -c*r, c*(b/a - a*t), sp.Integer(1))
# the transcription agrees with the repository's native C2 source at random rational points
random.seed(5)
names = {n: 'n', m: 'm', C: 'C', A: 'A', D: 'D', al: 'alpha'}
for _ in range(5):
    P = [Fr(random.randint(1, 30), 10), Fr(random.randint(1, 30), 10), Fr(random.randint(1, 60), 10),
         Fr(random.randint(1, 30), 10), Fr(random.randint(1, 99), 100)]
    s_repo, ct_repo, _ = vb.c2_source(*P)
    sub = dict(zip((c, a, b, t, r), [sp.Rational(x.numerator, x.denominator) for x in P]))
    for s_, nm in names.items():
        assert sp.Rational(Wt[s_].subs(sub)) == sp.Rational(s_repo[nm].numerator, s_repo[nm].denominator)
    for x_, y_ in zip(contact, ct_repo):
        assert sp.Rational(x_.subs(sub)) == sp.Rational(y_.numerator, y_.denominator)
print('PASS: weights and contact agree with verify_boundary.c2_source')

# 1. n at r = 1
n1 = sp.cancel(Wt[n].subs(r, 1))
assert sp.simplify(n1 + c**2*((b - a**2*t) + a**2)/(2*(c**2*(a**2 + b + 1) + 1))) == 0
print('PASS: n(r=1) = -c^2 (h + a^2)/(2(c^2(a^2+b+1)+1)) < 0 for h > 0')

# 2. r -> 0
n0 = sp.limit(Wt[n]*r, r, 0)
assert sp.simplify(n0 - b*(1 + a**2*c**2*t)/(2*(a**2*b*c**4 + a**2*b*c**2 + a**2*c**2 + b**2*c**4 + b**2*c**2
                                                    + b*c**2 + b + 1))) == 0
def no_pole_at_r0(expr):
    # the denominator does not vanish identically at r = 0; what remains are the chart gaps g, h, 1 - t, 1 - z h
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    return sp.expand(den.subs(r, 0)) != 0
assert no_pole_at_r0(Wt[n] - n0/r)
for s_ in (m, C, A, D, al):
    assert no_pole_at_r0(Wt[s_])
print('PASS: n = n0/r + O(1) with n0 > 0; m, C, A, D, alpha have no pole at r = 0')

# 3. orbit identity, first with symbolic weights, then on the C2 chart
uc, vc, Wc, Xc, Yc, Zc = contact
F = gl.form(gl.PAR).subs({d: 1, lam: 0, be: 0})
orbit = {u: L*uc, v: L*vc, W: Wc/L, X: Xc/L**2, Y: Yc/L**2, Z: Zc/L}
rhs = 4*C*(uc*vc*Wc)**2*(L - 1/L)**2 + ((4*n + 4*m + 2)*Xc**2 + (4*A + 4*C + 4*D + 2)*Yc**2)*(1 - L**-2)**2
diff = sp.Poly(sp.expand(sp.expand(F.subs(orbit) - rhs)*L**4), L)
for _, co in diff.terms():
    assert sp.cancel(sp.together(co.subs(Wt))) == 0
print('PASS: orbit identity F_c(O(L)) holds for every C2 source')
