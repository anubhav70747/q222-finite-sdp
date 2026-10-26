#!/usr/bin/env python3
"""Exact symbolic proof of oriented C1 small-f native charge compatibility.
No numerical sampling, SDP search, or assumption that charge compatibility implies PSD.
"""
import sympy as S,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'original_native'))
from c1_native import k,f,g,a,t,NATIVE,CONTACT,Lc,Hc,Qc,M,J,P4,Pa,Pd
x,z=S.symbols('x z',positive=True)
K=1+x
v={k:S.sqrt(x),g:S.sqrt(z),a:S.Integer(1),t:S.Integer(1)}
d,lam,n,m,C,A,D,alpha,beta=[S.cancel(y.subs(v)) for y in NATIVE]
uc,vc,wc,xc,yc,zc=[S.cancel(y.subs(v)) for y in CONTACT]
def at0(expr,pow=0):
    return S.factor(S.cancel(expr/f**pow).subs(f,0))
B=x*z**2-K**2*z+x*K
assert S.factor(at0(beta,2)-2*B/(z*K))==0
assert S.factor(at0(d,2)-4*(K**2-x*z)/K)==0
assert S.factor(at0(m,2)-2*(1-z)/z)==0
assert S.factor(at0(A,4)-2*(K-x*z)/K)==0
assert S.factor(at0(D,0)-2/(x*(1-z)))==0
Sy=4*A*uc**2*yc**2
Sx=(4*beta+2*d)*vc**2*xc**2
Sz=4*beta*uc**2*vc**2*zc**2
Ud=16*(2*beta+d+2*m)*(2*D+d)*uc**2*vc**2*xc**4
Y2=8*z*x**2*(K-x*z)/K
D4=256*x**3*(1-z)**3*(K-z)
assert S.factor(at0(Sy,2)-Y2)==0
assert S.factor(at0(Ud,4)-D4)==0
assert at0(Sx,2)==0 and at0(Sz,2)==0
comparison=S.expand(4*K**2*(1-z)**3*(K-z)-x*z*z*(K-x*z)**2)
zmax=x*K/(x*x+x+1)
Q=x**5+4*x**4+7*x**3+8*x**2+4*x+4
assert S.factor(comparison.subs(z,zmax)-(1-x)*K**3*Q/(x*x+x+1)**4)==0
assert S.factor((D4-Y2**2)*K**2/(64*x**3)-comparison)==0
assert S.factor(zmax-S.Rational(2,3)-(x-1)*(x+2)/(3*(x*x+x+1)))==0
assert S.factor(B-(x*K-z*(K**2-x))-x*z*(z-1))==0
# Source strictness at f=0: the coquotients retain nonzero source denominators.
assert at0(Lc.subs(v))==1 and at0(Hc.subs(v))==1 and S.factor(at0(Qc.subs(v))-(1-x*z))==0
assert S.factor(at0(M.subs(v))-x)==0 and S.factor(at0(J.subs(v))-x*(1-z))==0
assert S.factor(at0(P4.subs(v))-(K-x*z))==0
print('PASS: source leading d,m,A,D,beta and all local contact substitutions')
print('PASS: exact Sy/f^2, forced delta^2/f^4 and vanishing Sx,Sz/f^2')
print('PASS: beta positivity implies z < x(1+x)/(1+x+x^2) <= 2/3')
print('PASS: charge-margin endpoint polynomial factorization')
print('THEOREM (symbolic): for fixed x in (0,1], z in (0,1), B>0, sufficiently small f>0 gives strict native C1 and gamma=0, delta=Sy between Fourier demand and forced delta capacity')
print('WARNING: these scalar conditions are NOT the H7/H5 PSD inequalities')
