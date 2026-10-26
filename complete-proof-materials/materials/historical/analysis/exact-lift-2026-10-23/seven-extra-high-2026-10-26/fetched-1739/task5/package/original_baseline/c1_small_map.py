#!/usr/bin/env python3
"""C1 candidate 14+2 fixed integer and half-power feature map, exactly rational."""
import sympy as s
from itertools import combinations_with_replacement
u,v=s.symbols('u v');W,X,Y,Z=s.symbols('W X Y Z');psi=(W,X,Y,Z)
# Contact-parity groups (simultaneous u,v sign gauge): 8 odd, 6 even.
FEATS=[('uW',1,0,0),('vW',0,1,0),('X',0,0,1),('uvX',1,1,1),
 ('Y',0,0,2),('uvY',1,1,2),('uZ',1,0,3),('vZ',0,1,3),
 ('uvW',1,1,0),('vX',0,1,1),('uY',1,0,2),('vY',0,1,2),('Z',0,0,3),('uvZ',1,1,3)]
# Two homogeneous contact-kernel support groups, each normalized to evaluation 1
SETS=(list(range(8)),list(range(8,14)));ANCHORS=(7,13)
ROWS=[[(i,s.S.One),(anc,-s.S.One)] for group,anc in zip(SETS,ANCHORS) for i in group if i!=anc]
assert len(ROWS)==12
# odd sqrt(uv)W,Z squared difference (Wc,Zc normalized to 1)
ODD=[('sqrtuvW',s.Rational(1,2),s.Rational(1,2),0),('sqrtuvZ',s.Rational(1,2),s.Rational(1,2),3)]
PAIR=list(combinations_with_replacement(range(12),2))

def key(ai,bi,ci,di):return (ai,bi,min(ci,di),max(ci,di))

def column(i,j):
 d={}
 for p,a in ROWS[i]:
  for q,b in ROWS[j]:
   z=key(FEATS[p][1]+FEATS[q][1],FEATS[p][2]+FEATS[q][2],FEATS[p][3],FEATS[q][3]);d[z]=d.get(z,0)+(1 if i==j else 2)*a*b
 return {z:t for z,t in d.items() if t}

def odd_column():
 d={}
 for p,a in ((0,1),(1,-1)):
  for q,b in ((0,1),(1,-1)):
   z=key(ODD[p][1]+ODD[q][1],ODD[p][2]+ODD[q][2],ODD[p][3],ODD[q][3]);d[z]=d.get(z,0)+a*b
 return {z:t for z,t in d.items() if t}
COLS=[column(i,j) for i,j in PAIR]+[odd_column()]
# polynomial target (ordinary angle degrees 0...2). Some impossible terms excluded.
KEYS=sorted({key(i,j,l,m) for i in range(3) for j in range(3) for l in range(4) for m in range(l,4) if any((i,j,l,m) in c for c in COLS)},key=lambda k:k)
KEY_INDEX={k:i for i,k in enumerate(KEYS)}
M=s.MutableSparseMatrix(len(KEYS),len(COLS),{(KEY_INDEX[k],j):val for j,c in enumerate(COLS) for k,val in c.items()})

def sym_parent(coef):
 dd,la,n,m,C,A,D,al,be=coef;du=1+u*u;dv=1+v*v
 return s.expand(4*la*du*(W+v*X)**2+4*n*du*dv*X**2+4*m*(v*W-X+u*v*Y-u*Z)**2+4*C*dv*(u*W-Y)**2 + 4*A*du*dv*Y**2+(4*D+2*dd)*(u*W+u*v*X-Y-v*Z)**2+(4*al+2*dd)*du*(v*Y-Z)**2 +(4*be+2*dd)*dv*(X+u*Z)**2-2*dd*du*dv*Z**2)

def coefficients(F):
 co={};poly=s.Poly(F,u,v,W,X,Y,Z)
 for exp,val in poly.terms():
  ii,jj=exp[:2];z=exp[2:];inds=[i for i,c in enumerate(z) for _ in range(c)]
  if len(inds)!=2:raise ValueError(exp)
  name=key(ii,jj,*inds);co[name]=co.get(name,0)+val
 if not set(co).issubset(KEY_INDEX):raise ValueError('Unrepresented targets: '+str(set(co)-set(KEY_INDEX)))
 return s.Matrix([co.get(k,0) for k in KEYS])

if __name__=='__main__':
 print('C1 features', len(FEATS),len(ODD),'rows',len(ROWS),'matrix',M.shape,'rank',M.rank(),'null',M.cols-M.rank(),'cokernel',M.rows-M.rank(),flush=True)
 dd,la,n,m,C,A,D,al,be=s.symbols('d lambda n m C A D alpha beta')
 ff=sym_parent((dd,0,0,m,C,A,D,al,be)); fv=coefficients(ff)
 print('full native support is in coefficient map:',len(KEYS),'terms in parent:',sum(x!=0 for x in fv), flush=True)
 print('left null rank:',len(M.T.nullspace()),flush=True)

# Optional sign-gauge averaged block map: 7 odd-contact differences, five even-contact differences and one odd sqrtUV row.
BLKCOL=[j for j,(i,k) in enumerate(PAIR) if (i<7)==(k<7)]+[78]
BLOCK_MAP=M[:,BLKCOL]
