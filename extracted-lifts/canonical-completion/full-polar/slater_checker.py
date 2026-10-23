"""Exact rational Slater witness for the canonical order-33 half-line pencil.

Standard library only; literal source, selector, and cross-block matrices.
No source sweep, numerical SDP, or compiler import.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--pencil', type=Path, default=HERE.parent/'halfline'/'HALFLINE_QUADRATIC_SELECTOR_ORDER33.coo.json')
args = parser.parse_args()
raw = args.pencil.read_bytes()
data = json.loads(raw)
assert data['blocks'] == [8,8,4,4]+[1]*9
assert data['order'] == 33
assert data['indexing'] == 'zero-based; full symmetric COO; do not mirror or double entries'
source = [Q(1,100)] + [Q(1)]*8
J = [[[4,0],[0,4]], [[0,4],[4,0]], [[0,0],[0,8]]]
BG100 = [[-1,400,0,400],[-400,-1,0,0],[0,400,-1,801],[0,0,1,-1]]
BH100 = [[-1,-1,0,1],[-801,-1,401,0],[0,1,-1,400],[401,0,-400,-1]]
BG = [[Q(x,100) for x in row] for row in BG100]
BH = [[Q(x,100) for x in row] for row in BH100]
aux = [Q(x) for block in J for row in block for x in row]
aux += [x for block in [BG,BH] for row in block for x in row]
values = source+aux
assert len(aux)==44 and len(values)==len(data['parameter_order'])==53
M = [[Q(0) for j in range(33)] for i in range(33)]
for k,i,j,value in data['entries']:
    M[i][j] += Q(value)*values[k]
assert M == [list(row) for row in zip(*M)]
assert [M[24+i][24+i] for i in range(9)] == source
assert all(x>0 for x in source)

def pivots_and_minors(matrix):
    B = [row[:] for row in matrix]
    pivots, minors, running = [], [], Q(1)
    for k in range(len(B)):
        pivot = B[k][k]
        assert pivot > 0, (k,pivot)
        pivots.append(pivot)
        running *= pivot
        minors.append(running)
        for i in range(k+1,len(B)):
            for j in range(k+1,len(B)):
                B[i][j] -= B[i][k]*B[k][j]/pivot
    return pivots,minors

blocks, offset = [],0
for name,order in zip(['SG','SH','TG','TH']+data['source_order'],data['blocks']):
    block = [row[offset:offset+order] for row in M[offset:offset+order]]
    pivots,minors = pivots_and_minors(block)
    blocks.append({'name':name,'order':order,
                   'matrix':[[str(x) for x in row] for row in block],
                   'LDL_pivots':[str(x) for x in pivots],
                   'leading_principal_minors':[str(x) for x in minors]})
    offset += order
for index in [2,3]:
    assert blocks[index]['matrix']==[[str(Q(i==j,50)) for j in range(4)] for i in range(4)]

# A small hand bound is available as well as exact Sylvester minors.
# Coordinates are (a1,a2,b1,b2,va1,va2,vb1,vb2).
common = [[0,-1,0,0,1,0,0,0],[0,1,0,0,0,0,0,0],
          [0,0,0,0,0,1,0,0],[0,0,1,0,0,0,0,0],
          [0,0,0,0,0,0,1,0],[0,0,1,0,0,0,0,1]]
rowsG = common+[[1,0,1,0,0,1,0,1],[0,1,0,1,0,0,0,0],[0,0,0,0,0,1,0,1]]
rowsH = common+[[1,0,-1,0,0,0,0,0],[0,0,0,0,1,0,-1,0],[0,-1,0,1,1,0,-1,0]]
delta_bounds=[]
for index,rows in enumerate([rowsG,rowsH]):
    S0=[[4*Q(sum(row[i]*row[j] for row in rows)) for j in range(8)] for i in range(8)]
    S=[[Q(x) for x in row] for row in blocks[index]['matrix']]
    delta_bounds.append(max(sum(abs(S[i][j]-S0[i][j]) for j in range(8)) for i in range(8)))
assert delta_bounds == [Q(1,25),Q(7,100)]

result = {'input':args.pencil.name,'input_sha256':hashlib.sha256(raw).hexdigest(),
          'source_order':data['source_order'],'source':[str(x) for x in source],
          'selector':'J(v)=K(v)=[[4,4v],[4v,4+8v^2]]',
          'J_coefficient_matrices':J,
          'BG':[[str(x) for x in row] for row in BG],
          'BH':[[str(x) for x in row] for row in BH],
          'auxiliary_order':data['auxiliary_order'],'auxiliaries':[str(x) for x in aux],
          'all_source_coefficients_strictly_positive':True,
          'all_thirteen_PSD_blocks_positive_definite':True,
          'criterion':'All rational leading principal minors are strictly positive (Sylvester).',
          'hand_proof_baseline_lower_bounds':['4/5','2/3'],
          'exact_perturbation_row_bounds':[str(x) for x in delta_bounds],
          'hand_proof_final_lower_bounds':['19/25','179/300'],
          'blocks':blocks,
          'scope':'A rational Slater point of the evaluated candidate certificate cone. No universal adequacy or exact full-Q222 lift is asserted.'}
(HERE/'slater_witness.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({key:value for key,value in result.items() if key not in ('blocks','auxiliary_order','auxiliaries')},indent=2))
