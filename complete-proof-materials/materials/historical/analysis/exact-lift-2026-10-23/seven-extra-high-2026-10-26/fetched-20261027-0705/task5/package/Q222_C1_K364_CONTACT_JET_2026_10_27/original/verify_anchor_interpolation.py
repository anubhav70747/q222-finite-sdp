#!/usr/bin/env python3
"""Exact algebraic source-only contact interpolation, with Bob-zero Schur safety.
Calderon 2xn SOS theorem supplies an algebraic PSD contact skew;
this checker verifies interpolation identities, not the external theorem or global PSD.
"""
import itertools,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from verify_full_parent_bob import (S,P0,P1,P2,S0,S4,d,m,C,A,D,al,be,v,Dv)
vc=S.symbols('vc', positive=True)
source=2*D+d-2*m
S1=S.zeros(4)
for i,j in ((0,3),(1,2)):
 S1[i,j]=source;S1[j,i]=-source
raw=S.symbols('z0:6',real=True)
Sc=S.zeros(4)
for (i,j),coef in zip(itertools.combinations(range(4),2),raw):
 Sc[i,j]=coef;Sc[j,i]=-coef
T=Sc-S0-vc*S1-vc**4*S4
WX=S.zeros(4);WX[0,1]=T[0,1];WX[1,0]=-T[0,1]
S2=(T-WX)/vc**2; S3=WX/vc**3
assert S2[0,1]==0 and S2[1,0]==0
Sinterp=S0+v*S1+v**2*S2+v**3*S3+v**4*S4
assert (Sinterp.subs(v,vc)-Sc).applyfunc(S.cancel)==S.zeros(4)
assert Sinterp.subs(v,0)==S0
assert Sinterp.applyfunc(lambda ee:S.expand(ee).coeff(v,4))==S4
P=(Dv*P0).row_join(Dv*P1/2+Sinterp).col_join((Dv*P1/2-Sinterp).row_join(Dv*P2))
Dc=1+vc**2
J=(Dc*P0.subs(v,vc)).row_join(Dc*P1.subs(v,vc)/2+Sc).col_join(
 (Dc*P1.subs(v,vc)/2-Sc).row_join(Dc*P2.subs(v,vc)))
assert (P.subs(v,vc)-J).applyfunc(S.cancel)==S.zeros(8)
assert (P-P.T).applyfunc(S.cancel)==S.zeros(8)
print('PASS one source-only degree-4 skew polynomial with exact endpoints/contact')
print('PASS forced S2[W,X]=0, preserving both exact strictly positive Bob-zero pivots')
print('PASS arbitrary source-algebraic contact skew is transported without testing denominators')
print('UNPROVED: resulting 8x8 polynomial is PSD at every intermediate Bob v>0.')
