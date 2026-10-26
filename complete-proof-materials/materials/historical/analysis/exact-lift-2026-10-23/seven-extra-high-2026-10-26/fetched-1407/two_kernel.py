"""Exact two-kernel 14 even + 2 odd parameter-free linear coefficient map.
Root differences remain in the original Z_6 span.
"""
from map_core import *
GA=[0,1,2,3,4,12,13];PA=5
GB=[6,7,8,9,10];PB=11
ROWS=[[(i,S.Integer(1)),(PA,S.Integer(-1))] for i in GA]+[[(i,S.Integer(1)),(PB,S.Integer(-1))] for i in GB]
assert len(ROWS)==12
DIFFCOLS=[(i,j) for i in range(12) for j in range(i,12)]

def diff_col(i,j):
 o={}
 for p,c in ROWS[i]:
  for q,d in ROWS[j]:
   _,u1,v1,s1=FEATS[p];_,u2,v2,s2=FEATS[q]
   k=key(u1+u2,v1+v2,s1,s2)
   o[k]=o.get(k,0)+c*d*(1 if i==j else 2)
 return {k:v for k,v in o.items() if v}
M=S.MutableSparseMatrix(len(KEYS),79,{
 (KEY_IND[k],l):v for l,(i,j) in enumerate(DIFFCOLS) for k,v in diff_col(i,j).items()})
for k,v in gram_column(('odd',0,1),'laplacian').items():M[KEY_IND[k],78]=v
assert all(i in KEY_IND for col in range(M.cols) for i in [KEYS[row] for row in range(M.rows) if M[row,col]])
if __name__=='__main__':
 print('Two-kernel coefficient map shape',M.shape,'rank',M.rank(),'cokernel',M.rows-M.rank(),'Gram null',M.cols-M.rank())
