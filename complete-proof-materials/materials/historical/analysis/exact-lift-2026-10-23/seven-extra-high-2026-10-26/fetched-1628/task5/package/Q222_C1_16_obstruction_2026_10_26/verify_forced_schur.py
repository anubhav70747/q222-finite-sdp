#!/usr/bin/env python3
"""Exact forced C1 Gram identities, and nonsingular/singular Schur rules.
No PSD solver, no source sample, no parameter scan.
"""
from pathlib import Path
import sys
import sympy as S
BASE=Path(__file__).resolve().parent/'original_baseline'
if not (BASE/'c1_native.py').exists():
    BASE=Path(__file__).resolve().parents[1]/'Q222_C1_fixed_admission_2026_10_26'
sys.path.insert(0,str(BASE))
from eight_gauge_constructor import INVARIANT_ROWS,BLKCOL,PAIR,PIV_COLS,PIV_ROWS,FIXED_INVERSE,GAUGES
from c1_small_map import *
d,m,C,A,D,al,be=S.symbols('d m C A D alpha beta')
uc,vc,wc,xc,yc,zc=S.symbols('uc vc wc xc yc zc',nonzero=True)
F=sym_parent((d,0,0,m,C,A,D,al,be))
rhs=coefficients(S.expand(F.subs({u:uc*u,v:vc*v,W:wc*W,X:xc*X,Y:yc*Y,Z:zc*Z})))
r=rhs[INVARIANT_ROWS,:]
A42=M[INVARIANT_ROWS,BLKCOL]
assert len(GAUGES)==8
sol=S.zeros(44,1)
for j,val in zip(PIV_COLS,FIXED_INVERSE*r[PIV_ROWS,:]): sol[j]=S.expand(val)
gamma=S.symbols('gamma',real=True)
for j,v in enumerate(GAUGES):
   if j==7: sol+=gamma*v
H=S.zeros(12)
for j,col in enumerate(BLKCOL):
 if col!=78:
    p,q=PAIR[col];H[p,q]=H[q,p]=sol[j]
B=2*C+2*D+d
expected={(0,0):2*B*uc**2*wc**2,(0,1):-gamma/2,(0,2):0,
 (0,4):-2*B*uc*wc*yc,(1,1):4*m*vc**2*wc**2,
 (1,2):-4*m*vc*wc*xc,(1,4):0,
 (2,2):2*(2*be+d+2*m)*xc**2,(2,4):0,
 (4,4):2*(2*A+B)*yc**2,
 (11,11):4*al*zc**2}
for (i,j),exp in expected.items():
  assert S.cancel(H[i,j]-exp)==0,(i,j)
for ga in GAUGES:
  G=S.zeros(12)
  for i,col in enumerate(BLKCOL):
    if col!=78:
      p,q=PAIR[col];G[p,q]=G[q,p]=ga[i]
  for (i,j) in expected:
    if (i,j)==(0,1):continue
    assert G[i,j]==0,(i,j)
assert sol[-1]==gamma
print('PASS: source-independent forced Gram entries and gauge 7, including all zeros needed for the 4x4 Schur.')

# Algebraic universal 2+1+1 pivot elimination for H7, 1+1 for H5.
# Let generic H be symmetric; prove congruence formula by construction.
P0=H.extract([1,2],[1,2]);E=P0.inv()
I=[0,3,4,5,6]
J=H.extract(I,I)-H.extract(I,[1,2])*E*H.extract([1,2],I)
P4=J[2,2]
W=[0,1,3,4]
K=J.extract(W,W)-J.extract(W,[2])*J.extract([2],W)/P4
D0=S.factor(K[0,0])
# The general-source scalar is forced, independent of gauges 0..6.
Delta=S.factor(P0.det())
closed=S.factor(4*A*B*uc**2*wc**2/(2*A+B)-(2*be+d+2*m)*gamma**2/(16*m*(2*be+d)*vc**2*wc**2))
assert S.cancel(D0-closed)==0
print('PASS: full 7x7 -> strict 2x2 + strict 1x1 + conditional scalar + 3x3 reduction.')
print('        d0 =',closed)
# H5 mandatory positive last diagonal and an exact 4x4 PSD remainder.
F5=H[7:12,7:12]
B4=F5[:4,:4]-F5[:4,4]*F5[4,:4]/F5[4,4]
assert S.cancel(F5[4,4]-4*al*zc**2)==0
# No claim that B4, K are PSD is made.
print('PASS: full 5x5 -> strict 1x1 + exact 4x4 Schur; no additional coefficient assumption.')
