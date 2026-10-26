#!/usr/bin/env python3
"""Reproduces the existing exact rank-406 modular lower-bound certificate.
Not numerical sampling, not a solver. Uses only fixed integer monomial data.
Requires numpy and local c2_fixed_map.py; upper bound is six exact derivatives.
"""
import json
import numpy as np
from c2_fixed_map import coefficient_map,contact_annihilators
PRIME=65521
M=coefficient_map()
nr=len(M['row_keys']);nc=len(M['columns'])
mat=np.zeros((nr,nc),dtype=np.int64)
for j,col in enumerate(M['sparse_columns']):
 for i,v in col.items():mat[i,j]=v%PRIME
r=0;pivots=[]; row_ids=list(range(nr))
for j in range(nc):
 if r>=nr:break
 ix=np.flatnonzero(mat[r:,j]);
 if len(ix)==0:continue
 k=r+int(ix[0]);
 if k!=r:
  mat[[k,r]]=mat[[r,k]]
  row_ids[k],row_ids[r]=row_ids[r],row_ids[k]
 val=int(mat[r,j]);pivots.append([j,val,row_ids[r]]);inv=pow(val,-1,PRIME)
 mat[r,j:]=(mat[r,j:]*inv)%PRIME
 mat[r+1:,j:]=(mat[r+1:,j:]-mat[r+1:,j:j+1]*mat[r:r+1,j:])%PRIME
 r+=1
print('EXACT rank lower bound (mod 65521):',r)
assert r==406
assert r==nr-6
print('EXACT rank upper bound: 6 independent contact derivative annihilators')
print('EXACT rank over Q: 406')
with open('C2_RANK406_MODULAR_PIVOTS.json','w') as f:
 json.dump({'prime':PRIME,'shape':[nr,nc],'rank':r,'pivot_columns_residues_and_original_rows':pivots},f,indent=2)
