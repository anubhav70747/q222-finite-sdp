"""One-source exact obstruction to the recovered C1 three-channel allocation.

Python + SymPy. No solver, sampling, networking or file writes.
Reconstruct SG2, then precisely the displayed B1,B2,w allocation.
"""
import json
import sympy as S

u, v, s, t, W, X, Y, Z = S.symbols("u v s t W X Y Z")
d, m, C, A, D, alpha, beta = map(
    S.Rational,
    ["7301/540", "40/9", "2", "2891/135",
     "4699/1080", "81/40", "49/1080"],
)
assert all(a > 0 for a in (d, m, C, A, D, alpha, beta))
# Native SD1 parameters and the explicit strict gaps in the public example.
k = f_native = g_native = S.Integer(1)
a_native, t_native = S.Rational(7, 8), S.Rational(3, 5)
native_gaps = {"M": S.Rational(21, 40), "J": S.Rational(21, 40),
               "E": S.Rational(3, 5), "Q": S.Rational(7, 8)}
assert all(a > 0 for a in native_gaps.values())

Du, Dv = 1 + u*u, 1 + v*v
# Full cleared SG2 parent, with lambda=n=0.
F = (
    4*m*(v*W-X+u*v*Y-u*Z)**2
    + 4*C*Dv*(u*W-Y)**2 + 4*A*Du*Dv*Y**2
    + (4*D+2*d)*(u*W+u*v*X-Y-v*Z)**2
    + (4*alpha+2*d)*Du*(v*Y-Z)**2
    + (4*beta+2*d)*Dv*(X+u*Z)**2
    - 2*d*Du*Dv*Z**2
)
contact = {u: S.Integer(1), v: S.Integer(1), W: S.Integer(1),
           X: S.Rational(-2, 5), Y: S.Rational(1, 8), Z: S.Integer(1)}
assert contact[W] > 0 and contact[X] < 0 and contact[Y] > 0
assert F.subs(contact) == 0
assert all(S.diff(F, a).subs(contact) == 0 for a in (u, v, W, X, Y, Z))
# EN2 contact amplitude signs, with Zc=1.
contact_rows = {
    "uW-Y": u*W-Y, "1-vY": 1-v*Y,
    "Rm": v*W-X+u*v*Y-u*Z,
    "minus_RD": -(u*W+u*v*X-Y-v*Z),
    "X+u": X+u, "W+vX": W+v*X,
}
assert all(row.subs(contact) > 0 for row in contact_rows.values())

pc = {W: contact[W], X: contact[X], Y: contact[Y], Z: contact[Z]}
shift = {u: 1+s, v: 1+t}
gamma = S.Poly(S.expand(F.subs(pc).subs(shift)), s, t)
# psi=(xi+xi_c Z,Z), so zeta_W is the W derivative at the fixed contact state.
zeta_W = S.Poly(S.expand((S.diff(F, W)/2).subs(pc).subs(shift)), s, t)
gc = lambda i, j: gamma.coeff_monomial(s**i*t**j)
zc = lambda i, j: zeta_W.coeff_monomial(s**i*t**j)
K = S.Matrix([[gc(2, 0), gc(1, 1)/2],
              [gc(1, 1)/2, gc(0, 2)]])
g = S.Matrix([gc(2, 1)/2, gc(1, 2)/2])
f = gc(2, 2)
assert K[0, 0] > 0 and K.det() > 0
Kinv = K.inv()
rho = Kinv*g
sigma = S.factor(f - (g.T*Kinv*g)[0])
assert sigma > 0

Js, Jt, aa, bb, cc, dd, ee, jj = [
    zc(i, j) for i, j in
    [(1, 0), (0, 1), (2, 0), (0, 2),
     (1, 1), (2, 1), (1, 2), (2, 2)]
]
cprime = cc-rho[0]*Js-rho[1]*Jt
B1 = Js+aa*s+cprime*t
B2 = Jt+bb*t
w = (dd-rho[0]*aa)*s+(ee-rho[0]*cprime-rho[1]*bb)*t+jj*s*t
h1, h2 = s*(1+rho[0]*t), t*(1+rho[1]*s)
assert S.expand(B1*h1+B2*h2+w*s*t-zeta_W.as_expr()) == 0
assert S.expand(
    (S.Matrix([h1, h2]).T*K*S.Matrix([h1, h2]))[0]
    + sigma*s*s*t*t-gamma.as_expr()
) == 0

corner = {s: -1, t: -1}  # u=v=0
Bval = S.Matrix([S.factor(B1.subs(corner)), S.factor(B2.subs(corner))])
wval = S.factor(w.subs(corner))
H_WW = S.diff(F, W, W)/2
H_corner = H_WW.subs({u: 0, v: 0})
F_corner = F.subs({u: 0, v: 0, W: 1, X: 0, Y: 0, Z: 0})
R_corner = S.factor(H_corner-(Bval.T*Kinv*Bval)[0]-wval*wval/sigma)
assert F_corner == H_corner == 0
assert Bval[1] == -S.Rational(413, 9)
assert R_corner == -S.Rational(
    136060499747820067655375429183840,
    4375356724835654725138132413567,
)
assert R_corner < 0

print(json.dumps({
    "source_order": ["d", "m", "C", "A", "D", "alpha", "beta"],
    "source": list(map(str, (d, m, C, A, D, alpha, beta))),
    "contact": ["1", "1", "1", "-2/5", "1/8", "1"],
    "all_source_weights_and_native_gaps_positive": True,
    "contact_state_and_amplitude_signs": True,
    "full_parent_value_and_six_first_derivatives_zero": True,
    "K": [[str(a) for a in row] for row in K.tolist()],
    "g": list(map(str, g)), "rho": list(map(str, rho)),
    "sigma": str(sigma), "B_W_at_00": list(map(str, Bval)),
    "w_W_at_00": str(wval), "H_WW_at_00": str(H_corner),
    "F_at_00_eW": str(F_corner), "R_WW_at_00": str(R_corner),
    "scope": "Excludes this three-channel allocation, not the sixth-root cone.",
}, indent=2))

