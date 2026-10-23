"""Compile and check the quadratic half-line selector cone. Standard library only.

Writes only three local JSON files. No solver, source sampling, or adequacy
test is performed. All exported entries and polynomial identities are exact.
"""
from pathlib import Path
from fractions import Fraction as Q
import json

HERE = Path(__file__).resolve().parent
SOURCES = ['d', 'lambda', 'n', 'm', 'C', 'A', 'D', 'alpha', 'beta']
AUX = ([f'J{j}_{r}{s}' for j in range(3) for r in range(2) for s in range(2)]
       + [f'B{kind}_{r}{s}' for kind in ['G', 'H'] for r in range(4) for s in range(4)])
PARAMETERS = SOURCES + AUX
INDEX = {name: i for i, name in enumerate(PARAMETERS)}

def lin(**terms):
    return {INDEX[name]: value for name, value in terms.items() if value}

def add(*terms):
    out = {}
    for term in terms:
        for key, value in term.items():
            out[key] = out.get(key, 0) + value
    return {key: value for key, value in out.items() if value}

def neg(term):
    return {key: -value for key, value in term.items()}

def zero(n):
    return [[{} for _ in range(n)] for _ in range(n)]

def diagonal2(a, b):
    return [[a, {}], [{}, b]]

def offdiagonal2(a):
    return [[{}, a], [a, {}]]

def madd(A, B):
    return [[add(a, b) for a, b in zip(ar, br)] for ar, br in zip(A, B)]

def mneg(A):
    return [[neg(a) for a in row] for row in A]

def transpose(A):
    return [list(row) for row in zip(*A)]

def block(A, B, C, D):
    return [a+b for a, b in zip(A, B)] + [c+d for c, d in zip(C, D)]

# Coefficients of 1,v,v^2 in the authoritative P,N,E,K.
P = [diagonal2({}, lin(m=4, beta=4, d=2)),
     offdiagonal2(lin(m=-4)),
     diagonal2(lin(m=4), lin(beta=4, d=2))]
N = [diagonal2(lin(C=4, D=4, d=2), {}),
     offdiagonal2(lin(D=4, d=2)),
     diagonal2(lin(C=4), lin(D=4, d=2))]
E = [diagonal2(lin(A=4), lin(alpha=4)),
     offdiagonal2(lin(alpha=-4, d=-2)),
     diagonal2(lin(A=4, alpha=4, d=2), lin(d=-2))]
K = [diagonal2(lin(**{'lambda': 4}), lin(n=4)),
     offdiagonal2(lin(**{'lambda': 4})),
     diagonal2({}, lin(**{'lambda': 4, 'n': 4}))]
G, H = [], []
for j in range(3):
    J = [[lin(**{f'J{j}_{r}{s}': 1}) for s in range(2)] for r in range(2)]
    G.append(block(madd(P[j], K[j]), J, transpose(J), madd(N[j], K[j])))
    off = madd(madd(K[j], mneg(E[j])), mneg(J))
    H.append(block(madd(P[j], E[j]), off, transpose(off), madd(N[j], E[j])))

BG = [[lin(**{f'BG_{r}{s}': 1}) for s in range(4)] for r in range(4)]
BH = [[lin(**{f'BH_{r}{s}': 1}) for s in range(4)] for r in range(4)]
SG = block(G[0], BG, transpose(BG), G[2])
SH = block(H[0], BH, transpose(BH), H[2])
TG = madd(madd(G[1], mneg(BG)), mneg(transpose(BG)))
TH = madd(madd(H[1], mneg(BH)), mneg(transpose(BH)))

# Exact algebraic premises of the bounded-auxiliary/closed-projection lemma.
TRACE_G = [8, 16, 16, 8, 8, 0, 8, 0, 8]
TRACE_H = [8, 0, 0, 8, 8, 16, 8, 16, 8]
assert add(*(SG[i][i] for i in range(8))) == {i: x for i, x in enumerate(TRACE_G) if x}
assert add(*(SH[i][i] for i in range(8))) == {i: x for i, x in enumerate(TRACE_H) if x}
assert all(key < 9 for matrix in [SG, SH] for i in range(8) for key in matrix[i][i])
assert all(G[1][i][i] == H[1][i][i] == {} for i in range(4))
assert all(TG[i][i] == add(neg(BG[i][i]), neg(BG[i][i])) for i in range(4))
assert all(TH[i][i] == add(neg(BH[i][i]), neg(BH[i][i])) for i in range(4))
assert max(TRACE_G+TRACE_H) == 16

def export(blocks, orders, name, external_signs):
    entries, offset = [], 0
    for matrix, order in zip(blocks, orders):
        assert len(matrix) == order and all(len(row) == order for row in matrix)
        assert matrix == transpose(matrix)
        for r in range(order):
            for s in range(order):
                for parameter, coefficient in sorted(matrix[r][s].items()):
                    entries.append([parameter, offset+r, offset+s, coefficient])
        offset += order
    data = {
        'format': 'q222-halfline-quadratic-selector-affine-pencil-v1',
        'scope': 'Exact quadratic-selector cone on v>=0; equality with all canonical validity is NOT proved.',
        'order': sum(orders), 'blocks': orders,
        'source_order': SOURCES, 'auxiliary_order': AUX, 'auxiliary_count': len(AUX),
        'parameter_order': PARAMETERS, 'coefficient_equations_remaining': 0,
        'source_nonnegativity_external': external_signs,
        'indexing': 'zero-based; full symmetric COO; do not mirror or double entries',
        'entry_order': ['parameter_index', 'row', 'column', 'integer_coefficient'],
        'gram_block_order': ['SG', 'SH', 'TG', 'TH'] + ([] if external_signs else SOURCES),
        'trace_bounds': {
            'source_order': SOURCES, 'trace_SG_coefficients': TRACE_G,
            'trace_SH_coefficients': TRACE_H,
            'trace_TG_bound': '0 <= trace(TG) <= trace(SG)',
            'trace_TH_bound': '0 <= trace(TH) <= trace(SH)',
            'feasible_auxiliary_supnorm_bound': '16*sum(c_i) for c>=0; 16*sum(abs(c_i)) without external source signs',
            'closed_projection': True,
            'justification': 'Bounded-auxiliary lemma in halfline_gram.tex; trace and zero-diagonal identities checked exactly by compiler.',
        },
        'entries': entries,
    }
    (HERE/name).write_text(json.dumps(data, indent=2)+'\n')
    return data

core = export([SG, SH, TG, TH], [8, 8, 4, 4],
              'HALFLINE_QUADRATIC_SELECTOR_ORDER24.coo.json', True)
canonical = export([SG, SH, TG, TH] + [[[lin(**{name: 1})]] for name in SOURCES],
                   [8, 8, 4, 4]+[1]*9,
                   'HALFLINE_QUADRATIC_SELECTOR_ORDER33.coo.json', False)

# Independent direct-value reconstruction of the original matrices, not the
# coefficient arrays above. Three testing values determine each quadratic.
def numerical_matrices(parameters, v):
    d, la, n, m, C, A, D, al, be = parameters[:9]
    dv = 1+v*v
    p = [[4*m*v*v, -4*m*v], [-4*m*v, 4*m+(4*be+2*d)*dv]]
    nmat = [[4*C*dv+4*D+2*d, (4*D+2*d)*v], [(4*D+2*d)*v, (4*D+2*d)*v*v]]
    e = [[4*A*dv+(4*al+2*d)*v*v, -(4*al+2*d)*v], [-(4*al+2*d)*v, 4*al-2*d*v*v]]
    k = [[4*la, 4*la*v], [4*la*v, 4*la*v*v+4*n*dv]]
    j = [[sum(parameters[9+4*a+2*r+s]*v**a for a in range(3)) for s in range(2)] for r in range(2)]
    plus = lambda x, y: [[a+b for a, b in zip(ar, br)] for ar, br in zip(x, y)]
    minus = lambda x, y: [[a-b for a, b in zip(ar, br)] for ar, br in zip(x, y)]
    off = minus(minus(k, e), j)
    return (block(plus(p, k), j, transpose(j), plus(nmat, k)),
            block(plus(p, e), off, transpose(off), plus(nmat, e)))

def evaluate(data, parameters):
    out = [[Q(0) for _ in range(data['order'])] for _ in range(data['order'])]
    for parameter, row, column, coefficient in data['entries']:
        out[row][column] += coefficient*parameters[parameter]
    return out

def gram_polynomial(full, offsetS, offsetT, v):
    return [[full[offsetS+r][offsetS+s]
             + v*(full[offsetS+r][offsetS+4+s]+full[offsetS+4+r][offsetS+s])
             + v*v*full[offsetS+4+r][offsetS+4+s]
             + v*full[offsetT+r][offsetT+s]
             for s in range(4)] for r in range(4)]

checks = 0
for parameter in range(len(PARAMETERS)):
    values = [Q(int(i == parameter)) for i in range(len(PARAMETERS))]
    M = evaluate(canonical, values)
    assert M == transpose(M)
    assert [M[24+i][24+i] for i in range(9)] == values[:9]
    assert [row[:24] for row in M[:24]] == evaluate(core, values)
    for v in map(Q, [-1, 0, 2]):
        g, h = numerical_matrices(values, v)
        assert gram_polynomial(M, 0, 16, v) == g
        assert gram_polynomial(M, 8, 20, v) == h
        checks += 2
assert len(AUX) == 44
assert evaluate(canonical, [Q(0)]*len(PARAMETERS)) == [[Q(0)]*33 for _ in range(33)]
result = {
    'exact_coefficient_reconstruction_checks': checks,
    'method': 'Every one of 53 coefficient-basis vectors at three rational testing values; determines all quadratic coefficients.',
    'core_block_orders': [8, 8, 4, 4], 'core_order': 24,
    'canonical_block_orders': [8, 8, 4, 4]+[1]*9, 'canonical_order': 33,
    'auxiliary_count': 44, 'source_count': 9, 'coefficient_equations_remaining': 0,
    'core_full_symmetric_entries': len(core['entries']),
    'canonical_full_symmetric_entries': len(canonical['entries']),
    'zero_data_gives_zero_psd_pencil': True,
    'trace_SG_source_coefficients': TRACE_G,
    'trace_SH_source_coefficients': TRACE_H,
    'S_gram_diagonals_source_only': True,
    'middle_quadratic_coefficients_have_zero_diagonal': True,
    'T_diagonals_equal_minus_twice_B_diagonals': True,
    'feasible_auxiliary_supnorm_bound_constant': 16,
    'bound_uses': 'sum(c_i) for canonical orthant; sum(abs(c_i)) for unrestricted coefficient input',
    'scope': 'Exact half-line quadratic-selector cone only; no universal adequacy claim or solver run.',
}
(HERE/'halfline_gram_checks.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
