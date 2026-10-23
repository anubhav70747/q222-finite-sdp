#!/usr/bin/env python3
"""Export the fixed closed testing cover and all 80 linear Gram equations.

Standard-library exact arithmetic. The equations use the 228 independent
entries of A_G (12), B_G (8), A_H (12), B_H (8). They are equivalent to
the evaluated affine order-49 pencil; universal adequacy is not claimed.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import json
import sys
from verify_ramified import source_coefficients

HERE = Path(__file__).resolve().parent
CHECK_ONLY = '--check' in sys.argv
def emit(name, value):
    path = HERE/name
    if CHECK_ONLY:
        assert json.loads(path.read_text()) == json.loads(json.dumps(value)), name
    else:
        path.write_text(json.dumps(value, indent=2)+'\n')
source_names = ['d', 'lambda', 'n', 'm', 'C', 'A', 'D', 'alpha', 'beta']
block_names = ['A_G', 'B_G', 'A_H', 'B_H']
block_sizes = [12, 8, 12, 8]
gram_variables = [(name, i, j) for name, n in zip(block_names, block_sizes)
                  for i in range(n) for j in range(i, n)]
gram_index = {key: i for i, key in enumerate(gram_variables)}
assert len(gram_variables) == 228

def gram_coefficient(side, degree, i, j):
    terms = defaultdict(int)
    for letter, count, shift in [('A', 3, 0), ('B', 2, 1)]:
        for a in range(count):
            for b in range(count):
                if a+b+shift == degree:
                    row, col = sorted([4*a+i, 4*b+j])
                    terms[gram_index[(letter+'_'+side, row, col)]] += 1
    return terms

source_blocks = []
for i in range(9):
    c = [F(int(i == j)) for j in range(9)]
    source_blocks.append(source_coefficients(c))

rows = []
for k in range(5):
    for side in ['G', 'H']:
        for start in [0, 2]:
            for i in range(start, start+2):
                for j in range(i, start+2):
                    terms = gram_coefficient(side, k, i, j)
                    rhs = []
                    for src, (P, N, E, K) in enumerate(source_blocks):
                        base = P if start == 0 else N
                        extra = K if side == 'G' else E
                        x = base[k][i-start][j-start] + extra[k][i-start][j-start]
                        if x:
                            rhs.append([src, str(x)])
                    rows.append({'degree': k, 'entry': [i, j], 'condition': side+'_diagonal',
                                 'gram_terms': sorted(terms.items()), 'source_terms': rhs})
    for i in range(2):
        for j in range(2, 4):
            terms = gram_coefficient('G', k, i, j)
            for a, x in gram_coefficient('H', k, i, j).items():
                terms[a] += x
            rhs = []
            for src, (P, N, E, K) in enumerate(source_blocks):
                x = K[k][i][j-2]-E[k][i][j-2]
                if x:
                    rhs.append([src, str(x)])
            rows.append({'degree': k, 'entry': [i, j], 'condition': 'G_plus_H_cross',
                         'gram_terms': sorted(terms.items()), 'source_terms': rhs})
assert len(rows) == 80

# Independently compare this explicit equation map against every supplied
# affine basis matrix, and both supplied feasible auxiliary vectors.
pencil = json.loads((HERE/'ramified_order49_pencil.json').read_text())
offsets = dict(zip(block_names, [9, 21, 29, 41]))
for basis_index, basis in enumerate(pencil['source_matrices']+pencil['auxiliary_matrices']):
    sparse = {(i, j): F(x) for i, j, x in basis}
    values = [sparse.get((offsets[b]+i, offsets[b]+j), F(0)) for b, i, j in gram_variables]
    for row in rows:
        lhs = sum(F(x)*values[i] for i, x in row['gram_terms'])
        rhs = sum(F(x) for i, x in row['source_terms'] if i == basis_index)
        assert lhs == rhs

equations = {
    'meaning': 'sum(gram_terms)=sum(source_terms), for every row',
    'coordinate': 'v=16*s^2, s>=0', 'source_names': source_names,
    'gram_variables': [{'block': b, 'row': i, 'column': j} for b, i, j in gram_variables],
    'entry_convention': 'Each symmetric Gram entry is one variable, with zero-based indices.',
    'equation_count': 80, 'Gram_variable_count': 228,
    'source_independent_kernel_dimension': 148,
    'equations': rows,
    'eliminated_selector': 'J_k is the upper-right 2 by 2 block of the G polynomial coefficient.',
    'quartic_product_rule': 'A coefficient k sums A_(a,b) with a+b=k; B sums B_(a,b) with a+b+1=k.'
}
emit('linear_coefficient_equations.json', equations)

pg = ['a0^2*a1*W+a0*a1^2*Y', 'epsilon*(a0^2*a1*X+a0*a1^2*Z)',
      'a1^3*W-a0*a1^2*Y', 'epsilon*(a1^3*X-a0*a1^2*Z)']
ph = ['a0*a1^2*W+a1^3*Y', 'epsilon*(a0*a1^2*X+a1^3*Z)',
      'a0*a1^2*W-a0^2*a1*Y', 'epsilon*(a0*a1^2*X-a0^2*a1*Z)']
gn = ['tau^4', 'rho^2*tau^2', 'rho^4']
rn = ['rho*tau^3', 'rho^3*tau']
groups = []
for name, bob, alice in [('A_G', gn, pg), ('B_G', rn, pg), ('A_H', gn, ph), ('B_H', rn, ph)]:
    functions = [f'({x})*({y})/kappa' for x in bob for y in alice]
    groups.append({'block': name, 'size': len(functions), 'functions': functions})
assert sum(x['size'] for x in groups) == 40
cover = {
    'description': 'Forty evaluated algebraic functions on one fixed compact cover; no source coordinates occur.',
    'variables': ['a0', 'a1', 'rho', 'tau', 'epsilon', 'kappa', 'W', 'X', 'Y', 'Z'],
    'equations': ['a0^2+a1^2=1', 'rho^2+tau^2=1', 'epsilon^2=1',
                  'kappa^2=tau^8+256*rho^8'],
    'inequalities': ['a0>=0', 'rho>=0', 'tau>=0', 'kappa>0'],
    'state_domain': 'All real W,X,Y,Z; one may restrict W^2+X^2+Y^2+Z^2=1 by homogeneity. Every sign is retained.',
    'uniform_bounds': {'kappa_squared': '>=1/8', 'absolute_value_of_each_feature_on_unit_states': '<3'},
    'finite_chart': {'u': 'a1/a0', 'v': 'epsilon*16*(rho/tau)^4', 's': '(rho/tau)^2'},
    'Bob_unit_angle': ['tau^4/kappa', 'epsilon*16*rho^4/kappa'],
    'physical_slack': 'F_c^h=a0^2*b0^2*F_c(a1/a0,b1/b0;psi), interpreted by polynomial homogenization at endpoints.',
    'represented_target': 'a1^2*F_c^h', 'source_independent_multiplier': 'a1^2',
    'feature_blocks': groups,
    'Gram_product_equations': 'linear_coefficient_equations.json',
    'endpoint_values': {
        'Bob_zero': {'g': ['1', '0', '0'], 'r': ['0', '0']},
        'Bob_infinity': {'g': ['0', '0', '1/16'], 'r': ['0', '0']},
        'Alice_zero': 'Every feature equals zero; multiplier a1^2 is zero.',
        'Alice_infinity': {'P_G': ['0', '0', 'a1*W', 'a1*epsilon*X'],
                           'P_H': ['a1*Y', 'a1*epsilon*Z', '0', '0'], 'a1': 'either +1 or -1'}
    },
    'scope': 'The evaluated affine cone is a sufficient canonical certificate and strictly contains C33. Its equality with the full canonical valid cone remains unproved.'
}
emit('compact_fixed_functions.json', cover)
print('PASS: forty compact algebraic functions and all 80 exact linear Gram equations exported.')
print('PASS: equation map checked against all 157 evaluated affine basis matrices.')
