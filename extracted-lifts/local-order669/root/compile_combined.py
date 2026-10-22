#!/usr/bin/env python3
"""Combine only proved local fixed C1 certificate cones.

This is an evaluated source-cone pencil, NOT an exact lift of full Q222.
Python standard library only. Inputs and all arithmetic are exact rationals.
"""
import json
from pathlib import Path
from fractions import Fraction
from collections import defaultdict

HERE = Path(__file__).resolve().parent
POS = HERE.parent / 'uniform_tf1' / 'C1_BOTH_SIGNS_LOCAL_SOURCE_PENCIL.coo.json'
SIX = HERE.parent / 'candidate_boundary' / 'six_root_structured_source_pencil.coo.json'


def read_pencil(path):
    data = json.loads(path.read_text())
    if 'canonical_source_order' in data:
        # This source uses global matrix indices and k=0 for the constant.
        from bisect import bisect_right
        blocks = data['blocks']
        offsets = [b['offset'] for b in blocks]
        converted = []
        for k, i, j, value in data['entries']:
            bi, bj = bisect_right(offsets, i)-1, bisect_right(offsets, j)-1
            assert bi == bj and k > 0
            converted.append([bi, i-offsets[bi], j-offsets[bi], k-1, value])
        data = {'source_order': data['canonical_source_order'],
                'auxiliary_order': ['z_'+str(i) for i in range(data['auxiliary_count'])],
                'block_orders': data['block_orders'], 'coo': converted}
    elif 'auxiliary_order' not in data:
        data['auxiliary_order'] = ['y_'+str(i) for i in range(data['auxiliary_dimension'])]
    return data


def main():
    a, b = read_pencil(POS), read_pencil(SIX)
    # Sources are d,m,C,A,D,alpha,beta; the split source z is free.
    assert a['source_order'] == b['source_order']
    na, nb = len(a['auxiliary_order']), len(b['auxiliary_order'])
    entries = defaultdict(Fraction)
    for block, i, j, var, value in a['coo']:
        newvar = 7 + var if var < 7 else 14 + var - 7
        entries[block, i, j, newvar] += Fraction(value)
    offset = len(a['block_orders'])
    for block, i, j, var, value in b['coo']:
        v = Fraction(value)
        if var < 7:
            entries[block + offset, i, j, var] += v
            entries[block + offset, i, j, 7 + var] -= v
        else:
            entries[block + offset, i, j, 14 + na + var - 7] += v
    coo = [list(k) + [str(v)] for k, v in sorted(entries.items()) if v]
    orders = a['block_orders'] + b['block_orders']
    aux = ['split_' + s for s in a['source_order']]
    aux += ['plus_' + s for s in a['auxiliary_order']]
    aux += ['six_' + s for s in b['auxiliary_order']]
    output = {
        'schema': 'Q222.C1.fixed_local_sum_source_pencil.v1',
        'scope': 'Exact evaluated representation of (K_POS intersect K_NEG) + K_SIX, the printed LOCAL certificate cone; NOT a full-Q222 lift.',
        'equivalence': 'c belongs to the local sum certificate cone iff there exists z,y_plus,y_six with blockdiag(L_both(z,y_plus),L_six(c-z,y_six)) PSD.',
        'source_order': a['source_order'],
        'auxiliary_order': aux,
        'block_orders': orders,
        'total_matrix_order': sum(orders),
        'max_block_order': max(orders),
        'constant_entries': [],
        'coo': coo,
        'convention': '[block,i,j,variable,rational coefficient]; variable0..6 are canonical source coefficients, remaining variables are auxiliaries in auxiliary_order. Upper triangle stores actual matrix entries; mirror off-diagonal once.',
        'uniform_source_domain': {
            'formula_source': 'source_definitions.json / P1 in the proof note',
            'bounds': ['0<x<=2^-30', '0<sigma<=2^-30', '0<z<=1', '0<b<1', '0<h<1', 'r=b*z<1', 'E>0', 'M>0', 'J>0', 'D>=0', 'beta>=0'],
            'normalization': 'd=1; positive source scaling is also allowed',
            'lambda': 'No restriction: lambda>0 uses the first summand; lambda<=0 uses the second.',
            'H': 'No upper bound on H=h/sigma.'
        },
        'testing_domain': ['abs(u)<=2^-24', 'abs(w)<=2^-12', 'all real state vectors, and all complex state vectors by real/imaginary splitting'],
        'directions': {'canonical_source_to_pencil': 'Established on uniform_source_domain by the explicit fixed-family SOS construction.', 'pencil_to_local_nonnegativity': 'Established for every feasible source on testing_domain by the fixed Gram identities.', 'pencil_to_global_quantum_validity': 'Not established; not claimed.'}
    }
    # Verify the substitution exactly by pulling back to the input coefficient indices.
    inverse_a = defaultdict(Fraction)
    inverse_b = defaultdict(Fraction)
    for block, i, j, var, value in coo:
        v = Fraction(value)
        if block < offset:
            oldvar = var - 7 if var < 14 else var - 14 + 7
            inverse_a[block, i, j, oldvar] += v
        else:
            if var < 7:
                inverse_b[block-offset, i, j, var] += v
            elif var < 14:
                # Independent consistency check: coefficient of split is negative source.
                assert v == -entries[block, i, j, var-7]
            else:
                assert var >= 14 + na
                inverse_b[block-offset, i, j, var-(14+na)+7] += v
    expected_a = defaultdict(Fraction)
    expected_b = defaultdict(Fraction)
    for k, dest in ((a, expected_a), (b, expected_b)):
        for block, i, j, var, value in k['coo']:
            dest[block, i, j, var] += Fraction(value)
    assert dict(inverse_a) == dict(expected_a)
    assert dict(inverse_b) == dict(expected_b)
    assert len(aux) == 7 + na + nb
    assert all(0 <= k[0] < len(orders) and 0 <= k[1] <= k[2] < orders[k[0]] and 0 <= k[3] < 7+len(aux) for k in entries)
    (HERE / ('C1_UNIFORM_LOCAL_ORDER'+str(sum(orders))+'.coo.json')).write_text(json.dumps(output, indent=2)+'\n')
    receipt = {'status': 'passed', 'both_input_pencils_exactly_recovered': True, 'source_dimensions': 7, 'auxiliary_dimensions': len(aux), 'blocks': len(orders), 'total_matrix_order': sum(orders), 'max_block_order': max(orders), 'nonzero_upper_coefficient_records': len(coo), 'scope': output['scope']}
    (HERE / 'COMBINED_PENCIL_VERIFIED.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
