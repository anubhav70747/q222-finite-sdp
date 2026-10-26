import sys,sympy as S,itertools
sys.path.insert(0,'/mnt/data/Q222_C1_fixed_admission_2026_10_26')
from c1_small_map import FEATS, ROWS, PAIR, column, key, sym_parent, u,v,W,X,Y,Z,psi
from eight_gauge_constructor import *
ODD4=[('halfuv'+n,S.Rational(1,2),S.Rational(1,2),i) for i,n in enumerate(('W','X','Y','Z'))]
OROWS=[[(i,S.S.One),(3,-S.S.One)] for i in range(3)]
OPAIRS=list(itertools.combinations_with_replacement(range(3),2))
def oct(i,j):
 d={}
 for p,a in OROWS[i]:
  for q,b in OROWS[j]:
   zz=key(ODD4[p][1]+ODD4[q][1],ODD4[p][2]+ODD4[q][2],ODD4[p][3],ODD4[q][3]);d[zz]=d.get(zz,0)+(1 if i==j else 2)*a*b
 return {zz:val for zz,val in d.items() if val}
Bcols=[column(i,j) for i,j in PAIR if (i<7)==(j<7)] + [oct(i,j) for i,j in OPAIRS]
keys=sorted({zz for dic in Bcols for zz in dic})
mat=S.MutableSparseMatrix(len(keys),len(Bcols),{(keys.index(k),j):val for j,col in enumerate(Bcols) for k,val in col.items()})
print('18 feature raw sign-averaged matrix',mat.shape, 'rank',mat.rank(),'null',len(Bcols)-mat.rank(),'coker',len(keys)-mat.rank())
print('new odd half columns',OPAIRS)
print('new keys not prior',sorted(set(keys)-set(KEYS)) if 'KEYS' in globals() else '')
# Two contact blocks of 7,5 and one half block of 3, iff off parity outputs zero
print('odd columns:')
for ij in OPAIRS: print(ij,oct(*ij))
