"""Independent exact verification of the emitted sparse order49 pencil.

Reads the evaluated coefficient matrices, not the builder's expressions.
Each degree <=4 polynomial identity is proved at five exact integer points.
Also recovers every one of the148 auxiliary coordinates from its basis matrix.
Python standard library only. Run python3 verify_ramified_pencil.py.
"""
from fractions import Fraction as F
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
p=json.loads((HERE/'ramified_order49_pencil.json').read_text())
assert p['order']==49 and p['source_dimension']==9 and p['auxiliary_dimension']==148
assert len(p['source_matrices'])==9 and len(p['auxiliary_matrices'])==148

def matrix(data):
 A=[[F(0) for j in range(49)] for i in range(49)]
 for i,j,x in data:A[i][j]=A[j][i]=F(x)
 return A

def part(A,o,n):return [row[o:o+n] for row in A[o:o+n]]

def polycoeff(A,B,k):
 return [[sum(A[4*a+i][4*b+j] for a in range(3) for b in range(3) if a+b==k)
          +sum(B[4*a+i][4*b+j] for a in range(2) for b in range(2) if a+b+1==k)
          for j in range(4)] for i in range(4)]

def evalgram(A,B,s):
 return [[sum(A[4*a+i][4*b+j]*s**(a+b) for a in range(3) for b in range(3))
          +sum(B[4*a+i][4*b+j]*s**(a+b+1) for a in range(2) for b in range(2))
          for j in range(4)] for i in range(4)]

def physical(c,s):
 d,la,n,m,C,AA,D,al,be=c;v=16*s*s;dv=1+v*v
 P=[[4*m*v*v,-4*m*v],[-4*m*v,4*m+(4*be+2*d)*dv]]
 N=[[4*C*dv+4*D+2*d,(4*D+2*d)*v],[(4*D+2*d)*v,(4*D+2*d)*v*v]]
 E=[[4*AA*dv+(4*al+2*d)*v*v,-(4*al+2*d)*v],[-(4*al+2*d)*v,4*al-2*d*v*v]]
 K=[[4*la,4*la*v],[4*la*v,4*la*v*v+4*n*dv]]
 return P,N,E,K

for h,sparse in enumerate(p['source_matrices']+p['auxiliary_matrices']):
 L=matrix(sparse)
 c=[F(h==i) for i in range(9)]
 assert all(L[i][i]==c[i] for i in range(9))
 blocks=[part(L,o,n) for o,n in zip([9,21,29,41],[12,8,12,8])]
 AG,BG,AH,BH=blocks
 # At five points, each of the16 degree4 constraints holds. This proves
 # all80 coefficient equations, with no floating-point tests.
 for s in map(F,range(5)):
  G=evalgram(AG,BG,s);H=evalgram(AH,BH,s)
  P,N,E,K=physical(c,s)
  for i in range(2):
   for j in range(2):
    assert G[i][j]==P[i][j]+K[i][j]
    assert G[i+2][j+2]==N[i][j]+K[i][j]
    assert H[i][j]==P[i][j]+E[i][j]
    assert H[i+2][j+2]==N[i][j]+E[i][j]
    assert G[i][j+2]+H[i][j+2]==K[i][j]-E[i][j]
 # Recover the twenty selector entries and128 independent Gram freedoms.
 recovered={}
 for k in range(5):
  Q=polycoeff(AG,BG,k)
  for i in range(2):
   for j in range(2):recovered[f'J{k}_{i}{j}']=Q[i][j+2]
 for side,A,B in [('G',AG,BG),('H',AH,BH)]:
  for i in range(4):
   for j in range(i,4):
    recovered[f'{side}_B00_{i}{j}']=B[i][j]
    recovered[f'{side}_B11_{i}{j}']=B[i+4][j+4]
   for j in range(4):
    recovered[f'{side}_B01_{i}{j}']=B[i][j+4]
    recovered[f'{side}_A02_{i}{j}']=A[i][j+8]
   for j in range(i+1,4):
    recovered[f'{side}_K01_{i}{j}']=(A[i][j+4]-A[j][i+4])/2
    recovered[f'{side}_K12_{i}{j}']=(A[i+4][j+8]-A[j+4][i+8])/2
 assert [recovered[name] for name in p['auxiliary_names']]==[F(h==9+i) for i in range(148)]
print('PASS: all80 coefficient equations for all157 emitted basis matrices.')
print('PASS: all148 auxiliary coordinates recover exactly; no coefficient equations remain.')
