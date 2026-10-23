"""Exact self-contained validity check for the completion counterexample.

Standard library only. Expands the original four-state determinant directly,
without importing a source map, Gram archive, builder, or symbolic package.
Only the two named result files in --output-dir are written.
"""
from fractions import Fraction as Q
from pathlib import Path
from itertools import permutations
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source', type=Path, default=HERE/'ESTABLISHED_VALID_SOURCE.json')
parser.add_argument('--public-source', type=Path, default=HERE.parent/'outer-candidates/strict_outer_physical_source.json')
parser.add_argument('--output-dir', type=Path, default=HERE)
args = parser.parse_args()
expected_hash = '0479aeabc84cb1e3a71073d5162df83af8603c5f94609742498b0659ad793b75'
raw = args.source.read_bytes()
public_raw = args.public_source.read_bytes()
assert raw == public_raw
assert hashlib.sha256(raw).hexdigest() == expected_hash
source = json.loads(raw)
assert source['source_order'] == ['d', 'm', 'C', 'A', 'D', 'alpha', 'beta']
d, m, C, A, D, al, be = map(Q, source['source'])
assert d == 1 and min(m, C, A, D, al, be) > 0

def add(*polys):
    out = {}
    for poly in polys:
        for exponent, coefficient in poly.items():
            out[exponent] = out.get(exponent, Q(0)) + coefficient
    return {e: c for e, c in out.items() if c}

def scale(a, poly):
    return {e: a*c for e, c in poly.items() if a*c}

def mul(left, right):
    out = {}
    for e, a in left.items():
        for f, b in right.items():
            g = (e[0]+f[0], e[1]+f[1])
            out[g] = out.get(g, Q(0)) + a*b
    return {e: c for e, c in out.items() if c}

one = {(0, 0): Q(1)}
u, v, zero = {(1, 0): Q(1)}, {(0, 1): Q(1)}, {}
uv = mul(u, v)
Du, Dv = add(one, mul(u, u)), add(one, mul(v, v))
DuDv = mul(Du, Dv)
# Original SG2 packets at lambda=n=0, with state order (W,X,Y,Z).
rows = [
    (scale(4*m, one), [v, scale(-1, one), uv, scale(-1, u)]),
    (scale(4*C, Dv), [u, zero, scale(-1, one), zero]),
    (scale(4*A, DuDv), [zero, zero, one, zero]),
    (scale(4*D+2*d, one), [u, uv, scale(-1, one), scale(-1, v)]),
    (scale(4*al+2*d, Du), [zero, zero, v, scale(-1, one)]),
    (scale(4*be+2*d, Dv), [zero, one, zero, u]),
    (scale(-2*d, DuDv), [zero, zero, zero, one]),
]
P = [[add(*(mul(w, mul(row[i], row[j])) for w, row in rows))
      for j in range(4)] for i in range(4)]
determinant = {}
for perm in permutations(range(4)):
    sign = (-1)**sum(perm[i] > perm[j] for i in range(4) for j in range(i+1, 4))
    term = one
    for i, j in enumerate(perm):
        term = mul(term, P[i][j])
    determinant = add(determinant, scale(sign, term))

k = Q(365997345792000000, 806534883104527170621493981142180739537045858457)
weights = {
    (-1, 0): 1347083645295592961374032847731,
    (-1, 1): 2737791187307806620394799996019,
    (0, -1): 765091750295568657441489249147602,
    (0, 1): 762274793296209851400745596627083,
    (1, 0): 4005709020552400161419980319250,
    (1, 1): 79165812050999420348852524500,
}
assert k > 0 and min(weights.values()) > 0
balances = [sum(n*e[i] for e, n in weights.items()) for i in [0, 1]]
assert balances == [0, 0]
# u^2*v^2*k*sum n_e[(169*u^2)^e1*(v^2/256)^e2-1].
core = {(2, 2): -k*sum(weights.values())}
for e, n in weights.items():
    core[(2*e[0]+2, 2*e[1]+2)] = k*n*Q(169)**e[0]*Q(256)**(-e[1])
rhs = scale(16, mul(mul(Du, Du), mul(mul(Dv, Dv), core)))
assert determinant == rhs

u0, v0 = Q(source['contact']['u']), Q(source['contact']['v'])
assert (u0, v0) == (Q(1, 13), Q(16))
psi = list(map(Q, source['contact']['state']))
def evaluate(poly):
    return sum(c*u0**a*v0**b for (a, b), c in poly.items())
kernel = [sum(evaluate(P[i][j])*psi[j] for j in range(4)) for i in range(4)]
F = add(*(scale(psi[i]*psi[j], P[i][j]) for i in range(4) for j in range(4)))
derivatives = []
for coordinate in [0, 1]:
    derivative = {tuple(a-int(i == coordinate) for i, a in enumerate(e)): c*e[coordinate]
                  for e, c in F.items() if e[coordinate]}
    derivatives.append(evaluate(derivative))
assert kernel == [0]*4 and derivatives == [0, 0]

# The three positive SG2 packets give T^T diag(4m,4C*Dv,4A*Du*Dv) T,
# T=[[v,-1,uv],[u,0,-1],[0,0,1]], whose determinant is u.
T = [[v, scale(-1, one), uv], [u, zero, scale(-1, one)], [zero, zero, one]]
det_T = {}
for perm in permutations(range(3)):
    sign = (-1)**sum(perm[i] > perm[j] for i in range(3) for j in range(i+1, 3))
    term = one
    for i, j in enumerate(perm):
        term = mul(term, T[i][j])
    det_T = add(det_T, scale(sign, term))
assert det_T == u
assert args.source.read_bytes() == raw and args.public_source.read_bytes() == public_raw

data = {
    'source_sha256': expected_hash, 'source_order': source['source_order'],
    'source': source['source'], 'zero_source_coordinates': {'lambda': '0', 'n': '0'},
    'contact': source['contact'],
    'k': str(k), 'normalized_testing': {'X': '169*u^2', 'Y': 'v^2/256'},
    'integer_weights': [{'exponent': list(e), 'weight': str(n)} for e, n in weights.items()],
    'sum_weights': str(sum(weights.values())),
    'identity': 'det(P)=16*(1+u^2)^2*(1+v^2)^2*u^2*v^2*k*sum_e n_e*(X^e1*Y^e2-1)',
    'strict_transverse_frame': [['v', '-1', 'u*v'], ['u', '0', '-1'], ['0', '0', '1']],
    'frame_determinant': 'u',
    'positive_frame_gram_determinant': '64*m*C*A*u^2*(1+u^2)*(1+v^2)^2',
    'scope': 'Self-contained all-state canonical validity of this one fixed rational source; no large certificate or global adequacy claim.',
}
checks = {
    'source_matches_public_bytes': True, 'source_sha256': expected_hash,
    'inputs_unchanged': True, 'direct_determinant_polynomial_identity': True,
    'determinant_nonzero_terms': len(determinant),
    'determinant_bidegree': [max(e[i] for e in determinant) for i in [0, 1]],
    'six_weights_strictly_positive': True, 'integer_balances': balances,
    'full_state_contact_kernel': [str(x) for x in kernel],
    'fixed_state_angular_derivatives': [str(x) for x in derivatives],
    'strict_transverse_frame_determinant': 'u',
    'method': 'Exact Fraction polynomial arithmetic; 24 determinant permutations; no imported certificate/map or solver.',
}
assert args.output_dir.is_dir()
(args.output_dir/'pro_completion_selfcontained_data.json').write_text(json.dumps(data, indent=2)+'\n')
(args.output_dir/'pro_completion_selfcontained_checks.json').write_text(json.dumps(checks, indent=2)+'\n')
print(json.dumps(checks, indent=2))
