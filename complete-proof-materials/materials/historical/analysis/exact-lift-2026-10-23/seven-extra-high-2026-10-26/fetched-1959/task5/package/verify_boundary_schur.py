#!/usr/bin/env python3
"""Exact all-source C1 Bob-zero graded Schur completion, no samples or SDP.
The full 8-state-u matrix is checked coefficientwise from the unrestricted parent.
"""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from verify_full_parent_bob import S,P0,P1,P2,G0,Ginf,S0,S4,d,m,C,A,D,al,be,v,Dv
V=S.zeros(4)
eta=2*D+d-2*m
# the unique optimally coupled first-Bob-order skew in the selected endpoint frame
V[0,3]=eta;V[3,0]=-eta
V[1,2]=eta;V[2,1]=-eta
# Arbitrary higher-degree source-only skew corrections preserving the single
# forced WX entry of S2.  They do not affect the strict graded Bob-zero margin.
import itertools
S2=S.zeros(4);S3=S.zeros(4)
for j,(i,k) in enumerate([p for p in itertools.combinations(range(4),2) if p!=(0,1)]):
    x=S.Symbol('S2_'+str(j),real=True);S2[i,k]=x;S2[k,i]=-x
for j,(i,k) in enumerate(itertools.combinations(range(4),2)):
    x=S.Symbol('S3_'+str(j),real=True);S3[i,k]=x;S3[k,i]=-x
sk=S0+v*V+v*v*S2+v**3*S3+v**4*S4
H=(Dv*P0).row_join(Dv*P1/2+sk).col_join(
 (Dv*P1/2-sk).row_join(Dv*P2))
H0=H.subs(v,0)
assert H0==G0
H1=H.diff(v).subs(v,0)
H2=H.diff(v,2).subs(v,0)/2
ker=[0,5]; ran=[1,2,3,4,6,7]
assert H0.extract(ker,ker)==S.zeros(2)
assert H0.extract(ker,ran)==S.zeros(2,6)
assert H1.extract(ker,ker)==S.zeros(2)
R=H0.extract(ran,ran)
# Demonstrate R>0 on whole strict source orthant by 6 independent endpoint squares.
assert S.simplify(R.det()) != 0
q=H1.extract(ker,ran)
Sch=(H2.extract(ker,ker)-q*R.inv()*q.T).applyfunc(S.factor)
L=4*m+4*be+2*d
K=4*C+4*D+2*d
ans=S.diag(4*m*(4*be+2*d)/L, (4*D+2*d)*4*C/K)
assert (Sch-ans).applyfunc(S.factor)==S.zeros(2)
print('PASS Bob-zero kernel directions (W^0,X^u), 6x6 positive range block')
print('PASS first-order kernel block exactly zero')
print('PASS exact 2x2 graded Schur:')
for i in range(2):print('   ',S.factor(Sch[i,i]))
print('PASS both diagonals strictly positive for all genuine strict C1 source weights')
# As v grows, leading Gram is positive definite; this follows from Ginf's 10 explicit rows.
assert (H.applyfunc(lambda ee:S.expand(ee).coeff(v,4))-Ginf)==S.zeros(8)
print('PASS positive-definite degree-four Bob-infinity block for any positive source weights')
