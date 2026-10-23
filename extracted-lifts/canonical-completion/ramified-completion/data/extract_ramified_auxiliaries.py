"""Encode the independently certified ramified source in the affine pencil.

Run with a path to ramified_exact.json. Python standard library only.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent
path=Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'ramified_exact.json'
c=json.loads(path.read_text())
p=json.loads((HERE/'ramified_order49_pencil.json').read_text())
def mat(A):return [[F(x) for x in row] for row in A]
def tr(A):return list(map(list,zip(*A)))
def mul(A,B):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
def coeff(A,B,k):
 return [[sum(A[4*a+i][4*b+j] for a in range(3) for b in range(3) if a+b==k)
          +sum(B[4*a+i][4*b+j] for a in range(2) for b in range(2) if a+b+1==k)
          for j in range(4)] for i in range(4)]
blocks=[]
if 'Gram_matrices' in c:
 blocks=[mat(A) for A in c['Gram_matrices']]
else:
 for idx in range(4):
  frame=mat(c['face_matrices'][idx]);red=mat(c['reduced_Gram_blocks'][idx])
  scaled=mul(mul(frame,red),tr(frame))
  diagonal=[F(c['D_G' if idx<2 else 'D_H'][i][i]) for i in range(4)]
  native=[[x/diagonal[i%4]/diagonal[j%4] for j,x in enumerate(row)] for i,row in enumerate(scaled)]
  blocks.append(native)
vals={}
for k in range(5):
 Q=coeff(blocks[0],blocks[1],k)
 for i in range(2):
  for j in range(2):vals[f'J{k}_{i}{j}']=Q[i][j+2]
for side,A,B in [('G',*blocks[:2]),('H',*blocks[2:])]:
 for i in range(4):
  for j in range(i,4):
   vals[f'{side}_B00_{i}{j}']=B[i][j]
   vals[f'{side}_B11_{i}{j}']=B[i+4][j+4]
  for j in range(4):
   vals[f'{side}_B01_{i}{j}']=B[i][j+4]
   vals[f'{side}_A02_{i}{j}']=A[i][j+8]
  for j in range(i+1,4):
   vals[f'{side}_K01_{i}{j}']=(A[i][j+4]-A[j][i+4])/2
   vals[f'{side}_K12_{i}{j}']=(A[i+4][j+8]-A[j+4][i+8])/2
source=list(map(F,c['source']))
aux=[vals[name] for name in p['auxiliary_names']]
assert len(aux)==148
L=[[F(0) for j in range(49)] for i in range(49)]
for value,sparse in zip(source+aux,p['source_matrices']+p['auxiliary_matrices']):
 for i,j,x in sparse:
  L[i][j]+=value*F(x)
  if i!=j:L[j][i]+=value*F(x)
for i in range(9):assert L[i][i]==source[i]
for offset,A in zip([9,21,29,41],blocks):
 assert A==[row[offset:offset+len(A)] for row in L[offset:offset+len(A)]]
out={'description':'Exact auxiliary values for the excluded valid source in the ramified affine order49 cone.',
     'source_order':p['source_names'],'source':[str(x) for x in source],
     'auxiliary_order':p['auxiliary_names'],'auxiliary_values':[str(x) for x in aux],
     'block_orders':[1]*9+[12,8,12,8],
     'PSD_certificate':path.name}
if 'face_matrices' in c: out['non_scalar_block_ranks']=[7,3,11,7]
(HERE/('ramified_order49_cbar_auxiliaries.json' if 'Gram_matrices' in c else 'ramified_order49_source_auxiliaries.json')).write_text(json.dumps(out,indent=2)+'\n')
print('PASS: all148 auxiliary values reconstruct the four exact certified Gram blocks.')
