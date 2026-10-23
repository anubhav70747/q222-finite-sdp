#!/usr/bin/env python3
"""Exact verification of the ramified quartic obstruction (standard library).

Run: python3 verify.py
No optimizer, floating point, symbolic package, or regeneration of the
accepted order-49 pencil is used. Positivity is checked independently by
fraction-free leading principal determinants, rather than the producer's LDL.
The accepted validity/admission of cbar is an input; its all-source physical
transport is verified here as a homogeneous polynomial identity.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
data = json.loads((HERE / 'data/certificate.json').read_text())
N = [[[[int(x) for x in row] for row in mat] for mat in side]
     for side in data['moment_numerators']]
den = int(data['denominator'])
assert den == 15 * 2**60
assert len(N) == 2 and all(len(side) == 5 for side in N)
for side in N:
    for mat in side:
        assert len(mat) == 4 and all(len(row) == 4 for row in mat)
        assert all(mat[i][j] == mat[j][i] for i in range(4) for j in range(4))
assert all(N[0][k][i][j] == N[1][k][i][j]
           for k in range(5) for i in range(2) for j in range(2, 4))

def hankel(side, n, shift):
    return [[N[side][i//4 + j//4 + shift][i % 4][j % 4]
             for j in range(4*n)] for i in range(4*n)]

blocks = [hankel(0, 3, 0), hankel(0, 2, 1),
          hankel(1, 3, 0), hankel(1, 2, 1)]

def leading_determinants(mat):
    """Bareiss elimination, without pivoting: returns every leading minor."""
    a = [row[:] for row in mat]
    minors = []
    previous = 1
    for k in range(len(a)-1):
        pivot = a[k][k]
        assert pivot > 0
        minors.append(pivot)
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                numerator = pivot*a[i][j] - a[i][k]*a[k][j]
                value, remainder = divmod(numerator, previous)
                assert remainder == 0
                a[i][j] = value
        for i in range(k+1, len(a)):
            a[i][k] = a[k][i] = 0
        previous = pivot
    minors.append(a[-1][-1])
    assert min(minors) > 0
    return minors

minors = [leading_determinants(mat) for mat in blocks]
pivots = [[Q(a, b) for a, b in zip(ds, [1] + ds[:-1])] for ds in minors]
floor_pivots = [[p.numerator // p.denominator for p in ps] for ps in pivots]
stored_floors = [[str((p/15).numerator // (p/15).denominator) for p in ps]
                 for ps in pivots]
assert stored_floors == data['ldl_pivot_floors']

# Compute the functional from P,N,E,K directly. This calculation does not
# import the source map or the compiler used to build the accepted pencil.
def add(a, b):
    return [[x+y for x, y in zip(ar, br)] for ar, br in zip(a, b)]

def sub(a, b):
    return [[x-y for x, y in zip(ar, br)] for ar, br in zip(a, b)]

def block2(a, b, c, d):
    return [a[i] + b[i] for i in range(2)] + [c[i] + d[i] for i in range(2)]

def source_coefficients(c):
    d, la, n, m, C, A, D, al, be = c
    P = [[[0, 0], [0, 4*m+4*be+2*d]],
         [[0, -4*m], [-4*m, 0]],
         [[4*m, 0], [0, 4*be+2*d]]]
    NN = [[[4*C+4*D+2*d, 0], [0, 0]],
          [[0, 4*D+2*d], [4*D+2*d, 0]],
          [[4*C, 0], [0, 4*D+2*d]]]
    E = [[[4*A, 0], [0, 4*al]],
         [[0, -4*al-2*d], [-4*al-2*d, 0]],
         [[4*A+4*al+2*d, 0], [0, -2*d]]]
    K = [[[4*la, 0], [0, 4*n]],
         [[0, 4*la], [4*la, 0]],
         [[0, 0], [0, 4*la+4*n]]]
    zero = [[0, 0], [0, 0]]
    G = [block2(add(P[k], K[k]), zero, zero, add(NN[k], K[k]))
         for k in range(3)]
    H = [block2(add(P[k], E[k]), sub(K[k], E[k]),
                sub(K[k], E[k]), add(NN[k], E[k])) for k in range(3)]
    return G, H

def tracepair(a, b):
    return sum(a[i][j]*b[j][i] for i in range(len(a)) for j in range(len(a)))

ell = list(map(int, data['integer_source_functional']))
for j in range(9):
    c = [int(i == j) for i in range(9)]
    G, H = source_coefficients(c)
    val = sum(16**k*(tracepair(N[0][2*k], G[k]) +
                    tracepair(N[1][2*k], H[k])) for k in range(3))
    assert val == 12*ell[j]

# Independent comparison with the unchanged accepted COO arrays.
pencil_file = HERE / 'reused/ramified_order49_pencil.json'
pencil = json.loads(pencil_file.read_text())
assert pencil['source_names'] == data['source_order']
assert pencil['order'] == 49
offsets = [9, 21, 29, 41]

def entry(i, j):
    for off, mat in zip(offsets, blocks):
        if off <= i < off+len(mat) and off <= j < off+len(mat):
            return mat[i-off][j-off]
    return 0

def pair_coo(es):
    return sum((1 if i == j else 2)*Q(z)*entry(i, j) for i, j, z in es)

assert len(pencil['auxiliary_matrices']) == 148
assert all(pair_coo(A) == 0 for A in pencil['auxiliary_matrices'])
assert [pair_coo(A) for A in pencil['source_matrices']] == [12*x for x in ell]
cbar = list(map(Q, data['established_cbar']))
c = list(map(Q, data['source']))
pi = [0, 4, 5, 3, 1, 2, 6, 8, 7]
assert c == [cbar[i] for i in pi]
assert all(x > 0 for x in c)
negative = sum(a*b for a, b in zip(ell, c))
assert negative == Q(data['integer_functional_value']) < 0
assert Q(data['functional_scale']) == Q(12, den)
assert Q(data['negative_value']) == Q(12, den)*negative

# Homogeneous polynomial transport, including axes and projective endpoints.
# Variables (a0,a1,b0,b1,W,X,Y,Z); no divisions occur in these checks.
ZERO = (0,)*8
def padd(*ps):
    out = {}
    for p in ps:
        for e, x in p.items():
            out[e] = out.get(e, 0) + x
    return {e: x for e, x in out.items() if x}

def pscale(k, p):
    return {e: k*x for e, x in p.items() if k*x}

def pmul(p, q):
    out = {}
    for e, x in p.items():
        for f, y in q.items():
            g = tuple(a+b for a, b in zip(e, f))
            out[g] = out.get(g, 0) + x*y
    return {e: x for e, x in out.items() if x}

def psq(p):
    return pmul(p, p)

variables = [{tuple(int(i == j) for i in range(8)): 1} for j in range(8)]
def homogeneous_packets(args):
    a0, a1, b0, b1, W, X, Y, Z = args
    Da, Db = padd(psq(a0), psq(a1)), padd(psq(b0), psq(b1))
    dab = pmul(Da, Db)
    lroot = padd(pmul(b0, W), pmul(b1, X))
    mroot = padd(pmul(pmul(a0, b1), W), pscale(-1, pmul(pmul(a0, b0), X)),
                 pmul(pmul(a1, b1), Y), pscale(-1, pmul(pmul(a1, b0), Z)))
    croot = padd(pmul(a1, W), pscale(-1, pmul(a0, Y)))
    droot = padd(pmul(pmul(a1, b0), W), pmul(pmul(a1, b1), X),
                 pscale(-1, pmul(pmul(a0, b0), Y)), pscale(-1, pmul(pmul(a0, b1), Z)))
    aroot = padd(pmul(b1, Y), pscale(-1, pmul(b0, Z)))
    broot = padd(pmul(a0, X), pmul(a1, Z))
    events = [pmul(Da, psq(lroot)), pmul(dab, psq(X)), psq(mroot),
              pmul(Db, psq(croot)), pmul(dab, psq(Y)), psq(droot),
              pmul(Da, psq(aroot)), pmul(Db, psq(broot))]
    core = pscale(2, padd(events[5], events[6], events[7],
                          pscale(-1, pmul(dab, psq(Z)))))
    return [core] + [pscale(4, p) for p in events]

packets = homogeneous_packets(variables)
a0, a1, b0, b1, W, X, Y, Z = variables
transformed = homogeneous_packets([b1, pscale(-1, b0), a1, pscale(-1, a0), W, Y, X, Z])
assert all(packets[j] == transformed[pi[j]] for j in range(9))

# The transported exact contact retains the value and all six derivatives.
star = json.loads((HERE / 'reused/ESTABLISHED_VALID_SOURCE.json').read_text())
assert star['source_order'] == ['d', 'm', 'C', 'A', 'D', 'alpha', 'beta']
s0 = list(map(Q, star['source']))
source_star = [s0[0], Q(0), Q(0)] + s0[1:]
source_transposed = [source_star[i] for i in pi]
assert all(a >= b for a, b in zip(c, source_transposed))
Fstar = padd(*(pscale(a, p) for a, p in zip(source_transposed, packets)))
point = [Q(1), Q(1, 16), Q(1), Q(13), Q(140), Q(-1, 50), Q(1, 32), Q(1)]
def peval(p):
    total = Q(0)
    for e, a in p.items():
        value = a
        for x, n in zip(point, e):
            value *= x**n
        total += value
    return total

assert peval(Fstar) == 0
for i in [1, 3, 4, 5, 6, 7]:
    deriv = {tuple(n-int(j == i) for j, n in enumerate(e)): a*e[i]
             for e, a in Fstar.items() if e[i]}
    assert peval(deriv) == 0
star_negative = sum(a*b for a, b in zip(ell, source_transposed))
assert star_negative == Q(-4844338143654932884745484309348314410,
                         802811514300974151343419) < 0

print('PASS: all-source homogeneous transport, including every endpoint and zero source face.')
print('PASS: four positive-definite INTEGER Hankel matrices; all leading minors are positive.')
print('Minimum integer LDL-pivot floors:', [min(ps) for ps in floor_pivots])
print('PASS: dual trace functional derived directly from the accepted P,N,E,K formulas.')
print('PASS: unchanged accepted COO arrays: 148 auxiliary cancellations and all nine source terms.')
print('PASS: strictly positive source, and its exact contact predecessor, are both separated.')
print('ell(cdagger) =', negative)
print('Accepted pencil SHA256 =', hashlib.sha256(pencil_file.read_bytes()).hexdigest())
