"""Exact rank of the C1/C2 support-face Gram maps, without a numerical SDP.

An integer matrix is annihilated by six independent contact derivatives.
Gaussian elimination modulo the prime65521 supplies the matching lower
bound. Monomial scaling transfers the rank to every nonzero contact.
"""
from pathlib import Path
import json
import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
PRIME = 65521
assert sp.isprime(PRIME)

def keep(family, i, j, l):
    if l == 0:
        return i+j >= 6
    if family == 'C1' and l == 1:
        return j >= i
    if family == 'C2' and l == 3:
        return i+j <= 6
    return True

def product(a, b):
    i, j, l = a
    ii, jj, ll = b
    return i+ii, j+jj, min(l, ll), max(l, ll)

def modular_rank(matrix):
    a = matrix.copy() % PRIME
    r = 0
    pivots = []
    for c in range(a.shape[1]):
        nz = np.flatnonzero(a[r:, c])
        if not len(nz):
            continue
        pivot = r+int(nz[0])
        if pivot != r:
            a[[r, pivot]] = a[[pivot, r]]
        value = int(a[r, c])
        pivots.append([c, value])
        a[r, c:] = (a[r, c:]*pow(value, -1, PRIME)) % PRIME
        a[r+1:, c:] = (a[r+1:, c:]-a[r+1:, c:c+1]*a[r:r+1, c:]) % PRIME
        r += 1
        if r == a.shape[0]:
            break
    return r, pivots

out = {'prime': PRIME, 'families': {}}
for family in ['C1', 'C2']:
    blocks = [[(i, j, l) for l in range(4) for i in range(7) for j in range(7)
               if i % 2 == pi and j % 2 == pj and keep(family, i, j, l)]
              for pi, pj in [(0, 0), (0, 1), (1, 0), (1, 1)]]
    for suffix in ['A', 'B']:
        seed = json.loads((HERE/'boundary/data'/f'{family}{suffix}.json').read_text())
        assert [set(b) for b in blocks] == [set(map(tuple, b)) for b in seed['basis_blocks']]
    rows = sorted({product(a, b) for block in blocks for a in block for b in block})
    row_index = {key: i for i, key in enumerate(rows)}
    columns = []
    for block in blocks:
        ref = block[0]
        for p in range(1, len(block)):
            for q in range(p, len(block)):
                col = {}
                for a, b, sign in [(block[p], block[q], 1), (block[p], ref, -1),
                                   (ref, block[q], -1), (ref, ref, 1)]:
                    key = product(a, b)
                    col[key] = col.get(key, 0)+sign
                columns.append(col)
    matrix = np.zeros((len(rows), len(columns)), dtype=np.int64)
    for j, col in enumerate(columns):
        for key, value in col.items():
            matrix[row_index[key], j] = value
    derivatives = np.array([[int(l == i)+int(ll == i) for a, b, l, ll in rows]
                            for i in range(4)]+
                           [[a for a, b, l, ll in rows], [b for a, b, l, ll in rows]], dtype=np.int64)
    assert np.all(derivatives @ matrix == 0)
    assert sp.Matrix(derivatives.tolist()).rank() == 6
    rank, pivots = modular_rank(matrix)
    assert rank == len(rows)-6 == 406
    result = {'block_orders': list(map(len, blocks)), 'features': sum(map(len, blocks)),
              'coefficient_rows': len(rows), 'reduced_Gram_variables': len(columns),
              'independent_annihilators': 6, 'modular_rank_lower_bound': rank,
              'exact_rank': rank, 'pivot_columns_and_residues': pivots}
    out['families'][family] = result
    print(f'{family}: {len(rows)} rows; rank >= {rank} modulo{PRIME}, rank <= {len(rows)-6} over Q. Exact rank{rank}.', flush=True)
(HERE/'contact_rank.json').write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8', newline='\n')
