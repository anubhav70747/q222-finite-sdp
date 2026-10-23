"""Literal rational strict dual trace certificate, checked against order-24 COO.

Standard library only. Does not import the compiler or run an SDP solver.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
INPUT = HERE / 'HALFLINE_QUADRATIC_SELECTOR_ORDER24.coo.json'
raw = INPUT.read_bytes()
data = json.loads(raw)
h = Q(1, 2)
# Literal constant blocks, independent of every source and auxiliary.
WS = [[1,0,0,0,h,0,0,0], [0,1,0,0,0,h,0,0],
      [0,0,1,0,0,0,h,0], [0,0,0,1,0,0,0,h],
      [h,0,0,0,1,0,0,0], [0,h,0,0,0,1,0,0],
      [0,0,h,0,0,0,1,0], [0,0,0,h,0,0,0,1]]
WT = [[h,0,0,0], [0,h,0,0], [0,0,h,0], [0,0,0,h]]
W = [[Q(0) for j in range(24)] for i in range(24)]
offset = 0
for block in [WS, WS, WT, WT]:
    for i, row in enumerate(block):
        for j, value in enumerate(row):
            W[offset+i][offset+j] = Q(value)
    offset += len(block)
assert data['order'] == offset == 24
assert data['blocks'] == [8, 8, 4, 4]
assert data['indexing'] == 'zero-based; full symmetric COO; do not mirror or double entries'

# W = I/2 + R^T R/2, with R joining the paired four-dimensional blocks.
R = [[Q(0) for j in range(24)] for i in range(8)]
for k in range(4):
    R[k][k] = R[k][k+4] = 1
    R[k+4][k+8] = R[k+4][k+12] = 1
assert all(W[i][j] == Q(i == j, 2) + h*sum(R[k][i]*R[k][j] for k in range(8))
           for i in range(24) for j in range(24))
assert W[16][16] == h  # the lower eigenvalue is attained

coeff = [Q(0) for name in data['parameter_order']]
matrices = [[[Q(0) for j in range(24)] for i in range(24)] for name in coeff]
for k, i, j, value in data['entries']:
    matrices[k][i][j] += Q(value)
    # Full symmetric COO: sum every printed entry once; no mirroring.
    coeff[k] += W[j][i]*Q(value)
assert all(M == [list(row) for row in zip(*M)] for M in matrices)
assert coeff == [Q(16)]*9 + [Q(0)]*data['auxiliary_count']
assert data['auxiliary_count'] == 44

result = {
    'input': INPUT.name, 'input_sha256': hashlib.sha256(raw).hexdigest(),
    'literal_W_block_orders': [8,8,4,4],
    'literal_W_blocks': [[[str(Q(x)) for x in row] for row in b] for b in [WS,WS,WT,WT]],
    'W_lower_eigenvalue_exact': '1/2', 'trace_W': str(sum(W[i][i] for i in range(24))),
    'strict_PSD_proof': 'W=I24/2+R^T R/2; the TG coordinates attain eigenvalue1/2',
    'source_order': data['source_order'],
    'source_trace_coefficients': [str(x) for x in coeff[:9]],
    'auxiliary_trace_coefficients': [str(x) for x in coeff[9:]],
    'coefficient_count_checked': len(coeff),
    'identity': '<W,L(c,y)>=16*sum(c), for every real c,y',
    'PSD_trace_bound': 'L>=0 implies trace(L)<=32*sum(c)',
    'scope': 'Fixed-source matrix fibers are compact. The canonical nonnegative-source, sum(c)=1 matrix section is compact. Auxiliary representatives are bounded modulo the auxiliary kernel. The source projection is closed and infeasible fixed sources admit strict PSD separating certificates. W/20 is a strict feasible dual point for the specified phase-I dual; no arbitrary SDP dual-attainment assertion is made.',
    'phase_I': 'min t: L(c,y)+t*I24>=0. Dual max -<Y,L(c,0)>: Y>=0, <Y,Bj>=0, trace(Y)=1. Both sides are strictly feasible and finite; both optima are attained and equal. The optimal t is positive exactly for infeasible fixed-source fibers.'
}
(HERE/'strict_trace_certificate_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('literal_W_blocks','auxiliary_trace_coefficients')},indent=2))
