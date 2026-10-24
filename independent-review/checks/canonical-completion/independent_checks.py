"""Independent checks of claims in the canonical-completion notes.

Reads repo data read-only (path given as argv[1] = canonical-completion dir).
Writes nothing into the repo.
"""
import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

CC = Path(sys.argv[1])

u, v, U = sp.symbols('u v U')
W, X, Y, Z = sp.symbols('W X Y Z')
d, lam, n, m, C, A, D, al, be = sp.symbols('d lambda n m C A D alpha beta')
src = [d, lam, n, m, C, A, D, al, be]
Du, Dv = 1 + u**2, 1 + v**2


def parent(c, u, v, W, X, Y, Z):
    d, lam, n, m, C, A, D, al, be = c
    Du, Dv = 1 + u**2, 1 + v**2
    return (4*lam*Du*(W + v*X)**2 + 4*n*Du*Dv*X**2
            + 4*m*(v*W - X + u*v*Y - u*Z)**2
            + 4*C*Dv*(u*W - Y)**2 + 4*A*Du*Dv*Y**2
            + (4*D + 2*d)*(u*W + u*v*X - Y - v*Z)**2
            + (4*al + 2*d)*Du*(v*Y - Z)**2
            + (4*be + 2*d)*Dv*(X + u*Z)**2 - 2*d*Du*Dv*Z**2)


F = parent(src, u, v, W, X, Y, Z)
P = sp.Matrix([[4*m*v**2, -4*m*v], [-4*m*v, 4*m + (4*be + 2*d)*Dv]])
N = sp.Matrix([[4*C*Dv + 4*D + 2*d, (4*D + 2*d)*v], [(4*D + 2*d)*v, (4*D + 2*d)*v**2]])
E = sp.Matrix([[4*A*Dv + (4*al + 2*d)*v**2, -(4*al + 2*d)*v],
               [-(4*al + 2*d)*v, 4*al - 2*d*v**2]])
K = 4*lam*sp.Matrix([[1, v], [v, v**2]]) + 4*n*Dv*sp.Matrix([[0, 0], [0, 1]])

S, T, P0, Q0 = sp.symbols('S T P0 Q0')
a = sp.Matrix([S, T]); b = sp.Matrix([P0, Q0])
Phi = ((1 + U)*(a.T*P*a)[0] + U*(1 + U)*(b.T*N*b)[0]
       + U*((a - b).T*E*(a - b))[0] + ((a + U*b).T*K*(a + U*b))[0])

# 1. F_c = D_u * Phi_c under the change of variables.
sub = {W: S + u**2*P0, X: T + u**2*Q0, Y: u*(S - P0), Z: u*(T - Q0)}
r1 = sp.expand(F.subs(sub, simultaneous=True) - Du*Phi.subs(U, u**2))
print('F_c = D_u Phi_c identity:', r1 == 0)

# 2. Split identity with arbitrary J.
j = sp.symbols('j0:4'); J = sp.Matrix([[j[0], j[1]], [j[2], j[3]]])
G = sp.Matrix(sp.BlockMatrix([[P + K, J], [J.T, N + K]]))
H = sp.Matrix(sp.BlockMatrix([[P + E, K - E - J], [K - E - J.T, N + E]]))
x1 = sp.Matrix([S, T, U*P0, U*Q0]); x2 = sp.Matrix([S, T, P0, Q0])
r2 = sp.expand((x1.T*G*x1)[0] + U*(x2.T*H*x2)[0] - Phi)
print('split identity for every J:', r2 == 0)

# 3. Alice-Bob transport: u^2 v^2 F_{pi c}(-1/v,-1/u; S psi) = F_c(u,v;psi).
pic = [d, C, A, m, lam, n, D, be, al]
lhs = sp.expand(u**2*v**2*parent(pic, -1/v, -1/u, W, Y, X, Z))
print('pi transport identity:', sp.simplify(lhs - sp.expand(F)) == 0)
cbar = [1, Fraction(1, 2**34), Fraction(1, 2**14), Fraction(1, 2**24), Fraction(1, 576), 161,
        Fraction(1, 2), 76, Fraction(1, 250)]
perm = [0, 4, 5, 3, 1, 2, 6, 8, 7]
cdag = [cbar[i] for i in perm]
print('pi(cbar) =', [str(x) for x in cdag])

# 4. Motzkin square-root identity.
x, s = sp.symbols('x s')
M = lambda x, y: x**4*y**2 + x**2*y**4 - 3*x**2*y**2 + 1
rhs = (x**2*s**2 - 1)**2 + (x*s**2*(s**2 - 1))**2 + 2*(x*s*(s**2 - 1))**2
print('Motzkin half-line identity:', sp.expand(M(x, s**2) - rhs) == 0)

# 5. Forced polynomial selector minor at v=1 for the simple C1 source.
vals = {d: 359, lam: 0, n: 0, m: sp.Rational(75, 16), C: 60, A: 416, D: sp.Rational(41, 2),
        al: sp.Rational(831, 2), be: sp.Rational(151, 2)}
Jp = sp.Matrix([[-100*v**2, 0], [2000*v, 1200*v**2]])
Gs = sp.Matrix(sp.BlockMatrix([[P + K, Jp], [Jp.T, N + K]])).subs(vals).subs(v, 1)
sub2 = Gs.extract([1, 2], [1, 2])
print('G_J(1) on (a2,b1):', sub2.tolist(), ' value at (1,-1):',
      (sp.Matrix([1, -1]).T*sub2*sp.Matrix([1, -1]))[0])

# 6. Kernel dimension of the auxiliary adjoint on the core of the order-49 pencil.
pen = json.loads((CC/'ramified-completion/data/ramified_order49_pencil.json').read_text())
blocks = pen['blocks']
core = [bl for bl in blocks if bl['order'] > 1]
coords = {}
for bl in core:
    o, k = bl['offset'], bl['order']
    for i in range(k):
        for jj in range(i, k):
            coords[(o + i, o + jj)] = len(coords)
print('core block orders:', [bl['order'] for bl in core], ' Gram coords:', len(coords))


def adjoint_row(entries):
    row = [0]*len(coords)
    for i, jj, val in entries:
        i, jj = min(i, jj), max(i, jj)
        if (i, jj) in coords:
            # trace pairing: off-diagonal upper entries count twice
            row[coords[(i, jj)]] += sp.Rational(val)*(1 if i == jj else 2)
    return row


Baux = sp.Matrix([adjoint_row(e) for e in pen['auxiliary_matrices']])
rk = Baux.rank()
print('order-49 aux adjoint rank:', rk, ' kernel dim per chart:', len(coords) - rk,
      ' assembled order:', 8*pen['order'], ' assembled aux:', 8*(len(coords) - rk))
