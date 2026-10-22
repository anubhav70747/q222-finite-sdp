"""Independently reconstruct all polynomial equations of an extracted pencil.

Usage: python verify_weighted_pencil.py FEATURES.json MATRIX_DIRECTORY
Requires SymPy. No numerical arithmetic or SDP solver is used.
"""
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json
import sympy as sp


def verify(features, directory):
    data = json.loads(features.read_text(encoding='utf-8'))
    amap = json.loads((directory/'COEFFICIENT_MAP.coo.json').read_text())
    pencil = json.loads((directory/'SOURCE_PENCIL.coo.json').read_text())
    p,w,W,X,Y,Z = variables = sp.symbols('p w W X Y Z')
    ns = dict(zip(('p','w','W','X','Y','Z'),variables))
    assert data['testing_variables'] == list(ns)
    polynomial = lambda expr:sp.Poly(expr,*variables,domain=sp.QQ)
    parse = lambda text:sp.sympify(text,locals=ns)
    # Define the complete canonical parent independently of the input target.
    dp,dw = 1+p*p,1+w*w
    ev_m = p*W-p*w*X+Y-w*Z
    ev_D = w*W+X-p*w*Y-p*Z
    ev_a = Y-w*Z
    ev_b = p*X+Z
    targets = [2*ev_D**2+2*dp*ev_a**2+2*dw*ev_b**2-2*dp*dw*Z**2,
               4*ev_m**2,4*dw*(W-p*Y)**2,4*dp*dw*Y**2,
               4*ev_D**2,4*dp*ev_a**2,4*dw*ev_b**2]
    assert data['canonical_coefficient_order'] == ['d','m','C','A','D','alpha','beta']
    for label, target in zip(data['canonical_coefficient_order'],targets):
        assert polynomial(parse(data['target_coefficient_polynomials'][label])-target).is_zero
    basis = [tuple(e) for e in amap['row_monomials']]
    assert len(basis) == len(set(basis))
    column_polynomials = []
    expected_columns = []
    for k, block in enumerate(data['blocks']):
        vector = list(map(parse,block['testing_vector']))
        gate = parse(block['scalar_multiplier'])
        assert gate.free_symbols <= set(variables)
        for i in range(len(vector)):
            for j in range(i,len(vector)):
                expected_columns.append([k,i,j])
                column_polynomials.append(polynomial(gate*vector[i]*vector[j]*(1 if i==j else 2)))
    assert amap['column_gram_entries'] == expected_columns
    assert pencil['block_orders'] == [len(b['testing_vector']) for b in data['blocks']]
    nr,nc = amap['shape']
    assert (nr,nc) == (len(basis),len(expected_columns))
    rows = [{} for _ in basis]
    for i,j,c in amap['coo']:
        assert j not in rows[i]
        rows[i][j] = sp.Rational(c)
    for j, poly in enumerate(column_polynomials):
        actual = {e:rows[i][j] for i,e in enumerate(basis) if j in rows[i]}
        assert actual == poly.as_dict(), ('feature coefficient',j)
    tr = [{} for _ in basis]
    for i,j,c in amap['target_coo']:
        assert j not in tr[i]
        tr[i][j] = sp.Rational(c)
    for j, target in enumerate(targets):
        actual = {e:tr[i][j] for i,e in enumerate(basis) if j in tr[i]}
        assert actual == polynomial(target).as_dict(), ('canonical coefficient',j)
    lookup = {tuple(col):i for i,col in enumerate(expected_columns)}
    entries = [{} for _ in expected_columns]
    for block,i,j,v,c in pencil['coo']:
        assert i<=j
        col = lookup[block,i,j]
        assert v not in entries[col]
        entries[col][v] = sp.Rational(c)
    for i,row in enumerate(rows):
        actual = {}
        for j,c in row.items():
            for v,a in entries[j].items():
                actual[v] = actual.get(v,sp.S.Zero)+c*a
        actual = {v:c for v,c in actual.items() if c}
        assert actual == tr[i], ('affine pencil equation',i)
    free = pencil['free_gram_columns']
    assert len(free) == len(set(free))
    assert all(entries[j] == {7+k:sp.S.One} for k,j in enumerate(free))
    # Independent elimination verifies that no additional free directions were lost.
    echelon = {}
    for row in rows:
        row = row.copy()
        while row:
            pivot = min(row)
            factor = row[pivot]
            if pivot not in echelon:
                echelon[pivot] = {j:c/factor for j,c in row.items()}
                break
            for j,c in echelon[pivot].items():
                new = row.get(j,sp.S.Zero)-factor*c
                if new: row[j] = new
                else: row.pop(j,None)
    assert len(echelon)+len(free) == nc
    assert set(echelon) == set(pencil['pivot_gram_columns'])
    assert not pencil['constant_entries']
    result = {'status':'passed','feature_sha256':hashlib.sha256(features.read_bytes()).hexdigest(),
              'polynomial_equations_checked':nr,'independent_rank':len(echelon),
              'free_directions_checked':len(free),'canonical_source_coefficients':7,
              'total_matrix_order':sum(pencil['block_orders']),
              'qualification':'Exact equivalence of the serialized pencil and the given weighted Gram cone. This does not assert full quantum-set coverage.'}
    (directory/'independent_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('features',type=Path)
    ap.add_argument('matrix_directory',type=Path)
    args=ap.parse_args()
    verify(args.features,args.matrix_directory)
