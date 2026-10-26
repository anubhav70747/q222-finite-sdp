#!/usr/bin/env python3
"""Strictly positive forced C1 diagonal masses and 3 forced 2x2 Schur margins.
These are partial principal minors, NOT full 7x7 or 5x5 PSD.
"""
from c1_small_map import *
uc,vc,wc,xc,yc,zc=s.symbols('uc vc wc xc yc zc',nonzero=True)
d,m,C,A,D,al,be=s.symbols('d m C A D alpha beta')
raw=coefficients(sym_parent((d,0,0,m,C,A,D,al,be)))
scales=(uc,vc,wc,xc,yc,zc)
sol_map=M[:,BLKCOL]
for i in range(12):
 p=i if i<7 else i+1
 _,uu,vv,st=FEATS[p];k0=key(2*uu,2*vv,st,st)
 idx=KEY_INDEX[k0];inz=[(j,sol_map[idx,j]) for j in range(44) if sol_map[idx,j]]
 assert inz==[(BLKCOL.index(PAIR.index((i,i))),s.S.One)],(i,k0,inz)
 expected=raw[idx]*uc**(2*uu)*vc**(2*vv)*(wc,xc,yc,zc)[st]**2
 assert expected!=0,(i,k0)
 print(f'DIAGONAL {i:2d}: {s.factor(expected)}')
# Three forced principal blocks, independent of 8 gauge parameters.
blocks=[
 ((1,2),4*m*vc**2*wc**2,2*(2*m+2*be+d)*xc**2,-4*m*vc*wc*xc,
     8*m*(2*be+d)*vc**2*wc**2*xc**2),
 ((0,4),2*(2*C+2*D+d)*uc**2*wc**2,2*(2*A+2*C+2*D+d)*yc**2,-2*(2*C+2*D+d)*uc*wc*yc,
     8*A*(2*C+2*D+d)*uc**2*wc**2*yc**2),
 ((0,3),2*(2*C+2*D+d)*uc**2*wc**2,2*(2*D+d)*uc**2*vc**2*xc**2,2*(2*D+d)*uc**2*vc*wc*xc,
     8*C*(2*D+d)*uc**4*vc**2*wc**2*xc**2),
]
for (i,j),a,b,h,det in blocks:
 assert s.expand(a*b-h*h-det)==0,(i,j)
 # Cross monomial has unique coefficient 2 in reduced map.
 fi=i if i<7 else i+1; fj=j if j<7 else j+1
 ki=key(FEATS[fi][1]+FEATS[fj][1],FEATS[fi][2]+FEATS[fj][2],FEATS[fi][3],FEATS[fj][3]);index=KEY_INDEX[ki]
 support=[(p,sol_map[index,p]) for p in range(44) if sol_map[index,p]]
 assert support==[(BLKCOL.index(PAIR.index((i,j))),s.Integer(2))]
 exp=raw[index]*uc**ki[0]*vc**ki[1]*(wc,xc,yc,zc)[ki[2]]*(wc,xc,yc,zc)[ki[3]]
 assert s.cancel(exp-2*h)==0,(i,j)
 print(f'FORCED {i},{j} Schur determinant: {s.factor(det)} strictly positive if d,m,C,A,beta>0 and contact entries nonzero')
print('PASS: 12 forced diagonal masses and 3 strictly positive source-uniform 2x2 principal blocks. NO FULL PSD CLAIM.')
