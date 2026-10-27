#!/usr/bin/env python3
"""One targeted exact rational regression for the universal analytic contact-jet formula.
This is NOT evidence of universal source positivity; that is proved in PROOF.md.
No numerical sweep, no SDP or approximate witness construction.
"""
import sys
from pathlib import Path
import sympy as s
sys.path.insert(0,str(Path(__file__).resolve().parent))
from source_only_contact_jet import source_contact_jet
sys.path.insert(0,str(Path(__file__).resolve().parent/'original'))
import c1_native as src
# A native strict already-published elementary analytic contact, NOT a search over points.
params={src.k:s.Integer(1),src.f:s.Integer(1),src.g:s.Integer(1),src.a:s.Rational(3,5),src.t:s.Rational(3,5)}
weights=[s.factor(w.subs(params)) for w in src.NATIVE]
assert weights == [s.Integer(8),0,0,2,2,3,1,s.Rational(8,3),s.Rational(2,3)]
uc,vc,*psi=[s.factor(x.subs(params)) for x in src.CONTACT]
assert [uc,vc,*psi]==[1,1,1,s.Rational(-2,5),s.Rational(2,5),1]
v=s.symbols('v',real=True)
u=s.symbols('u',real=True)
W,X,Y,Z=s.symbols('W X Y Z',real=True)
d,m,C,A,D,al,be=weights[0],weights[3],weights[4],weights[5],weights[6],weights[7],weights[8]
Du=1+u*u;Dv=1+v*v
F=(4*m*(v*W-X+u*v*Y-u*Z)**2+4*C*Dv*(u*W-Y)**2+4*A*Du*Dv*Y**2
  +(4*D+2*d)*(u*W+u*v*X-Y-v*Z)**2+(4*al+2*d)*Du*(v*Y-Z)**2
  +(4*be+2*d)*Dv*(X+u*Z)**2-2*d*Du*Dv*Z**2)
Q=s.hessian(F,(W,X,Y,Z))/2
P0=Q.subs(u,0);P1=Q.diff(u).subs(u,0);P2=Q.diff(u,2).subs(u,0)/2
r=source_contact_jet(v,P0,P1,P2,uc,vc,s.Matrix(psi),s.Integer(0))
print('contact complement rational Schur pivots source gauge tau=0')
IA=[0,1,6,7];IB=[2,3,4,5]
for name,ix,kx in [('A',IA,r['kerA']),('B',IB,r['kerB'])]:
    zv=s.Matrix([kx[i] for i in ix]); mat=r['Hcontact'].extract(ix,ix)
    Qj=s.Matrix.hstack(*[s.eye(4)[:,j]-zv[j]/zv[0]*s.eye(4)[:,0] for j in [1,2,3]])
    cm=(Qj.T*mat*Qj).applyfunc(s.cancel)
    piv=[s.factor(cm[:t,:t].det()) for t in [1,2,3]]
    assert all(e>0 for e in piv), (name,piv)
    print(name, list(map(str,piv)))
print('native strict physical contact second-order Schur curvature',s.factor(r['physical_schur']))
assert r['physical_schur']>0
print('PASS exact source-only first-derivative minimization and equal positive second-order soft pivots')
xx=s.Matrix([W,X,Y,Z]).col_join(u*s.Matrix([W,X,Y,Z]))
assert s.expand((xx.T*r['Hpol']*xx)[0]-Dv*F)==0
print('PASS whole u-v-four-state coefficient matching for the interpolated contact polynomial')
print('NOTE: one rational regression only; no numerical sweep or global-v PSD inference.')
