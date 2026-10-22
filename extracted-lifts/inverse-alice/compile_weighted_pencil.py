"""Compile a rational affine pencil from a fixed polynomial weighted Gram family.

Usage: python compile_weighted_pencil.py FEATURES.json OUTPUT_DIRECTORY
The compiler uses exact rational arithmetic and does not invoke an SDP solver.
Its conclusion concerns precisely the supplied certificate cone.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sympy as s


def add(dst, src, factor=F(1)):
    for key, val in src.items():
        v = dst.get(key, F(0)) + factor * val
        if v:
            dst[key] = v
        else:
            dst.pop(key, None)


def fs(x):
    return str(x.numerator) if x.denominator == 1 else f'{x.numerator}/{x.denominator}'


def compile_family(path, out):
    raw = path.read_bytes()
    data = json.loads(raw)
    vs = s.symbols(' '.join(data['testing_variables']))
    names = dict(zip(data['testing_variables'], vs))
    def poly(text):
        expr = s.sympify(text, locals=names)
        assert expr.free_symbols <= set(vs), 'The testing family depends on a source.'
        return {tuple(e): F(int(c.p), int(c.q))
                for e, c in s.Poly(expr, *vs, domain=s.QQ).terms() if c}
    columns, products = [], []
    ambient = set()
    for k, block in enumerate(data['blocks']):
        v = [s.sympify(a, locals=names) for a in block['testing_vector']]
        w = s.sympify(block['scalar_multiplier'], locals=names)
        for i in range(len(v)):
            for j in range(i, len(v)):
                columns.append([k, i, j])
                pp = poly(w*v[i]*v[j]*(1 if i == j else 2))
                products.append(pp)
                ambient.update(pp)
    source_order = data['canonical_coefficient_order']
    target = [poly(data['target_coefficient_polynomials'][v]) for v in source_order]
    ambient.update(e for pp in target for e in pp)
    ambient = sorted(ambient)
    rowid = {e:i for i,e in enumerate(ambient)}
    rows = [{} for _ in ambient]
    for j, pp in enumerate(products):
        for e, c in pp.items():
            rows[rowid[e]][j] = c
    rhs = [{} for _ in ambient]
    for j, pp in enumerate(target):
        for e, c in pp.items():
            rhs[rowid[e]][j] = c
    pivots = {}
    for row, tr in zip(rows, rhs):
        a, b = row.copy(), tr.copy()
        while a:
            p = min(a)
            if p not in pivots:
                div = a[p]
                pivots[p] = ({j:c/div for j,c in a.items()},
                             {j:c/div for j,c in b.items()})
                break
            pivot, pb = pivots[p]
            factor = a[p]
            add(a, pivot, -factor)
            add(b, pb, -factor)
        else:
            assert not b, f'Additional visible source constraint: {b}'
    nsource = len(source_order)
    free = [j for j in range(len(columns)) if j not in pivots]
    exprs = {j:{nsource+i:F(1)} for i,j in enumerate(free)}
    for p in sorted(pivots, reverse=True):
        a, b = pivots[p]
        expr = b.copy()
        for j, c in a.items():
            if j != p:
                add(expr, exprs[j], -c)
        exprs[p] = expr
    # Check every original polynomial coefficient, including dependent rows.
    for row, tr in zip(rows, rhs):
        check = {}
        for j, c in row.items():
            add(check, exprs[j], c)
        assert check == tr
    assert all(exprs[j] == {nsource+i:F(1)} for i,j in enumerate(free))
    orders = [len(b['testing_vector']) for b in data['blocks']]
    coo = [[*columns[j], v, fs(c)] for j, expr in sorted(exprs.items())
           for v, c in sorted(expr.items())]
    counts = {'blocks':len(orders), 'total_order':sum(orders),
              'maximum_block_order':max(orders), 'gram_entries':len(columns),
              'coefficient_equations':len(rows), 'coefficient_rank':len(pivots),
              'auxiliary_dimension':len(free), 'affine_entries':len(coo)}
    matrix = {'schema':'q222.weighted_gram_coefficient_map.v1',
              'row_monomials':ambient, 'column_gram_entries':columns,
              'shape':[len(rows),len(columns)],
              'coo':[[i,j,fs(c)] for i,row in enumerate(rows) for j,c in sorted(row.items())],
              'target_source_order':source_order,
              'target_coo':[[i,j,fs(c)] for i,row in enumerate(rhs) for j,c in sorted(row.items())]}
    pencil = {'schema':'q222.weighted_source_pencil.v1',
              'scope':data['scope'],
              'claim':'Exact representation of the supplied fixed weighted Gram cone; no global coverage assertion.',
              'source_order':source_order, 'block_orders':orders,
              'constant_entries':[], 'coo':coo, 'counts':counts,
              'free_gram_columns':free, 'pivot_gram_columns':sorted(pivots),
              'convention':'[block,i,j,variable,rational coefficient]; variables 0..6 are source coefficients; later variables are independent auxiliaries. Mirror each off-diagonal entry.'}
    receipt = {'status':'passed', 'input_sha256':hashlib.sha256(raw).hexdigest(),
               'counts':counts, 'verified_original_coefficient_equations':len(rows),
               'verified_free_coordinates':len(free),
               'qualification':'The compiler establishes both directions between the fixed weighted Gram cone and this pencil. Source admission and global coverage require separate proofs.'}
    out.mkdir(parents=True, exist_ok=True)
    for name, item in [('COEFFICIENT_MAP.coo.json',matrix),
                       ('SOURCE_PENCIL.coo.json',pencil),('compilation.json',receipt)]:
        (out/name).write_text(json.dumps(item,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt,indent=2))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('features', type=Path)
    ap.add_argument('output_directory', type=Path)
    args = ap.parse_args()
    compile_family(args.features,args.output_directory)
