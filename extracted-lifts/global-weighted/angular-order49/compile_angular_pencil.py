"""Write the rational affine pencil for the nine-dimensional Laurent cone."""
from pathlib import Path
from fractions import Fraction as F
import json

HERE = Path(__file__).resolve().parent
E = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1) if (i, j) != (0, 0)]
pairs = [[e, (-e[0], -e[1])] for e in [(1, 0), (0, 1), (1, 1), (1, -1)]]
tri = [[(s, t), (-s, 0), (0, -t)] for s in (-1, 1) for t in (-1, 1)]
weighted = [[(-1, 1), (-1, -1), (1, 0)],
            [(1, 1), (1, -1), (-1, 0)],
            [(1, -1), (-1, -1), (0, 1)],
            [(1, 1), (-1, 1), (0, -1)]]
names = ['h'] + [f'q_{i}_{j}' for i, j in E]
names += [f'r_{k}_{j}' for k in range(8) for j in range(3)]
names += [f's_{k}' for k in range(12)]
names += [f'x_{k}' for k in range(8)] + [f'y_{k}' for k in range(4)]
assert len(names) == 57 and len(set(names)) == 57
ix = {name: j for j, name in enumerate(names)}
records, blocks = [], []
offset = 0

def linear(**terms):
    return {k: F(v) for k, v in terms.items() if v}

def entry(i, j, terms):
    for name, value in sorted(terms.items(), key=lambda kv: ix[kv[0]]):
        value = F(value)
        if value:
            records.append([offset+i, offset+j, ix[name], value.numerator, value.denominator])

def two(name, a, b, c):
    global offset
    blocks.append({'name': name, 'offset': offset, 'order': 2})
    entry(0, 0, a); entry(0, 1, b); entry(1, 1, c)
    offset += 2

def residual(e):
    ans = {f'q_{e[0]}_{e[1]}': F(2)}
    for k, circuit in enumerate(tri + weighted):
        for j, f in enumerate(circuit):
            if e == f:
                ans[f'r_{k}_{j}'] = F(-2)
    return ans

for k, (e, f) in enumerate(pairs):
    two(f'pair_{k}', residual(e), linear(**{f's_{k}': 1}), residual(f))
for k in range(4):
    t = linear(**{f's_{k+4}': F(1, 3)})
    x, y = linear(**{f'x_{k}': 1}), linear(**{f'y_{k}': 1})
    two(f'triple_{k}_a', linear(**{f'r_{k}_0': 1}), x, linear(**{f'r_{k}_1': 1}))
    two(f'triple_{k}_b', linear(**{f'r_{k}_2': 1}), y, t)
    two(f'triple_{k}_c', x, t, y)
for k in range(4, 8):
    x = linear(**{f'x_{k}': 1})
    two(f'weighted_{k-4}_a', linear(**{f'r_{k}_0': 1}), x, linear(**{f'r_{k}_1': 1}))
    two(f'weighted_{k-4}_b', x, linear(**{f's_{k+4}': F(1, 4)}), linear(**{f'r_{k}_2': F(1, 2)}))
blocks.append({'name': 'constant_remainder', 'offset': offset, 'order': 1})
entry(0, 0, {'h': 1, **{f's_{k}': 1 for k in range(12)}})
offset += 1
assert offset == 49 and len(blocks) == 25 and len(records) == 109
data = {
    'scope': 'Exact affine lift of the nine-dimensional nonnegative Laurent cone with nonnegative outer coefficients. Not an affine lift of Q222.',
    'input_variables': names[:9], 'auxiliary_variables': names[9:],
    'variables': names, 'constant_matrix': 'zero',
    'laurent_exponents': E, 'pairs': pairs, 'equal_triangles': tri, 'weighted_triangles': weighted,
    'order': offset, 'blocks': blocks,
    'coo_convention': '[row,column,variable_index,numerator,denominator], zero-based; mirror every off-diagonal entry, with no factor of one half',
    'entries': sorted(records),
}
(HERE/'ANGULAR_CONE_ORDER49.coo.json').write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8', newline='\n')
print(json.dumps({'order': offset, 'order2_blocks': 24, 'scalar_blocks': 1,
                  'inputs': 9, 'auxiliaries': 48, 'rational_entries': len(records)}))
