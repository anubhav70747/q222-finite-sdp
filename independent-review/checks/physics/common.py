"""Shared setup. g_k = sqrt(c_k) * h_k, so g g^T = c_k h h^T is polynomial (sqrt only appears squared)."""
import itertools, random
import sympy as sp
from sympy import Matrix, Rational, symbols, expand, zeros, eye

u, v = symbols('u v', real=True)
lam, n, m, C, A, D, d, al, be = symbols('lam n m C A D d alpha beta', real=True)
W, X, Y, Z = symbols('W X Y Z', real=True)
psi = Matrix([W, X, Y, Z])
Du, Dv = 1 + u**2, 1 + v**2

NAMES = ['lam', 'n', 'm', 'C', 'A', 'D', 'alpha', 'beta']
# weights a_k in the order of NAMES
a_orig = [4*lam, 4*n, 4*m, 4*C, 4*A, 4*D + 2*d, 4*al + 2*d, 4*be + 2*d]
# scalar c_k and polynomial vector h_k  (g_k = sqrt(c_k) h_k)
c = [Du, Du*Dv, 1, Dv, Du*Dv, 1, Du, Dv]
h = [Matrix(x) for x in ([1, v, 0, 0], [0, 1, 0, 0], [v, -1, u*v, -u], [u, 0, -1, 0],
                         [0, 0, 1, 0], [u, u*v, -1, -v], [0, 0, v, -1], [0, 1, 0, u])]
r = 2*d*Du*Dv                      # ell = sqrt(r) e4
e4 = Matrix([0, 0, 0, 1])
LL = r * e4 * e4.T                 # ell ell^T


def T_of(a):
    """T = sum a_k c_k h_k h_k^T for any weight list a."""
    return sum((a[k]*c[k]*h[k]*h[k].T for k in range(8)), zeros(4, 4))


F = (4*lam*Du*(W + v*X)**2 + 4*n*Du*Dv*X**2 + 4*m*(v*W - X + u*v*Y - u*Z)**2
     + 4*C*Dv*(u*W - Y)**2 + 4*A*Du*Dv*Y**2 + (4*D + 2*d)*(u*W + u*v*X - Y - v*Z)**2
     + (4*al + 2*d)*Du*(v*Y - Z)**2 + (4*be + 2*d)*Dv*(X + u*Z)**2 - 2*d*Du*Dv*Z**2)

T = T_of(a_orig)
P = T - LL


def cross3(a, b, cc):
    """Generalized cross product of three 4-vectors: w_i = det[a b cc e_i]."""
    return Matrix([Matrix.hstack(a, b, cc, Matrix([1 if j == i else 0 for j in range(4)])).det()
                   for i in range(4)])


def rand_q(rng, lo=1, hi=9):
    return Rational(rng.randint(lo, hi), rng.randint(1, 5))
