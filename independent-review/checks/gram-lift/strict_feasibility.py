"""Exact strictly feasible point of the full sextic Gram system; no solver.

Run from any folder with Python and SymPy. Writes the rational witness beside
this script. This proves strict feasibility, not exactness of L_6.
"""
from pathlib import Path
import json
import sympy as s

HERE = Path(__file__).resolve().parent
u, v = s.symbols('u v')
W, X, Y, Z = PSI = s.symbols('W X Y Z')
corners = [(0, 0), (1, 0), (0, 1), (1, 1)]
z1 = s.Matrix([u**i*v**j*p for i, j in corners for p in PSI])

def vector(f):
    f = s.expand(f)
    return s.Matrix([f.coeff(p).coeff(u, i).coeff(v, j)
                     for i, j in corners for p in PSI])

def gram(roots):
    return sum((vector(f)*vector(f).T for f in roots), s.zeros(16))

event_roots = [
    W+v*X, u*(W+v*X),
    X, u*X, v*X, u*v*X,
    v*W-X+u*v*Y-u*Z,
    u*W-Y, v*(u*W-Y),
    Y, u*Y, v*Y, u*v*Y,
    u*W+u*v*X-Y-v*Z,
    v*Y-Z, u*(v*Y-Z),
    X+u*Z, v*(X+u*Z),
]
Q0 = 4*gram(event_roots)
Qd_plus = 2*gram([u*W+u*v*X-Y-v*Z,
                  v*Y-Z, u*(v*Y-Z), X+u*Z, v*(X+u*Z)])
Qd_minus = 2*gram([Z, u*Z, v*Z, u*v*Z])
Qd = Qd_plus-Qd_minus
assert Q0.rank() == 16
assert s.trace(Q0.inv()) == s.Rational(21, 4)
assert Qd_minus == s.diag(*([0, 0, 0, 2]*4))
eps = s.Rational(1, 21)
Qc = Q0+eps*Qd-eps*s.eye(16)
assert all(Qc[:i, :i].det() > 0 for i in range(1, 17))

# Each interior monomial's diagonal is cancelled by one off-diagonal.
H = s.diag(1, *([s.Rational(1, 2)]*5), 1)
for i in range(1, 6):
    H[i-1, i+1] = H[i+1, i-1] = -s.Rational(1, 4)
r = s.symbols('r')
zr = s.Matrix([r**i for i in range(7)])
assert s.expand((zr.T*H*zr)[0]-1-r**12) == 0
assert all(H[:i, :i].det() > 0 for i in range(1, 8))

basis = [(i, j, l) for i in range(7) for j in range(7) for l in range(4)]
indices = {b: i for i, b in enumerate(basis)}
J = s.zeros(196, 16)
for col, (i, j, l) in enumerate((i, j, l) for i, j in corners for l in range(4)):
    J[indices[(6*i, 6*j, l)], col] = 1
G = J*Qc*J.T+eps*s.kronecker_product(H, H, s.eye(4))

# Exact coefficient identity, with u=r^6 and v=t^6.
coefficients = {}
records = []
for i in range(196):
    for j in range(i, 196):
        value = G[i, j]
        if value == 0:
            continue
        records.append([i, j, str(value)])
        ai, bi, li = basis[i]
        aj, bj, lj = basis[j]
        key = (ai+aj, bi+bj, min(li, lj), max(li, lj))
        coefficients[key] = coefficients.get(key, 0)+(1 if i == j else 2)*value
coefficients = {k: x for k, x in coefficients.items() if x != 0}
target = s.Poly(s.expand((z1.T*(Q0+eps*Qd)*z1)[0]), u, v, *PSI)
expected = {}
for exponents, value in target.terms():
    a, b, *state = exponents
    ls = [l for l, power in enumerate(state) for _ in range(power)]
    expected[(6*a, 6*b, *ls)] = value
assert coefficients == expected

out = {
    'claim': 'Strict feasibility of L6 and all nine source signs; not exactness of L6.',
    'source_order': ['d', 'lambda', 'n', 'm', 'C', 'A', 'D', 'alpha', 'beta'],
    'source': ['1/21']+['1']*8,
    'basis': basis,
    'gram_upper_triangle': records,
    'symmetry': 'Unlisted entries are zero; reflect off-diagonal records.',
    'epsilon': str(eps),
    'corner_Q0': [[str(x) for x in row] for row in Q0.tolist()],
    'corner_Qd': [[str(x) for x in row] for row in Qd.tolist()],
    'H': [[str(x) for x in row] for row in H.tolist()],
    'trace_Q0_inverse': '21/4',
    'strict_bound': 'G >= (H tensor H tensor I4)/21 > 0',
    'coefficient_identity_exact': True,
}
(HERE/'strict_feasibility.json').write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8', newline='\n')
print(json.dumps({k: out[k] for k in ['claim', 'source', 'trace_Q0_inverse', 'strict_bound', 'coefficient_identity_exact']}, indent=2))
print('PSD order:', len(basis), '; nonzero upper-triangle records:', len(records))
