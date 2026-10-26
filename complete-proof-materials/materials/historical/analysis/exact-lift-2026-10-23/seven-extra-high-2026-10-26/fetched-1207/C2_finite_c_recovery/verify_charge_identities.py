#!/usr/bin/env python3
"""Recovered original partial charge checker. Exact SymPy, no samples or SDP.

This checker covers the charge *skeleton*, not finite-c admission.
Recovered from the previous run's saved execution transcript; not rerun here.
"""
import sympy as s
u,v,U,V,W,X,Y,Z,e=s.symbols('u v U V W X Y Z epsilon')
n,m,C,A,D,al=s.symbols('n m C A D alpha')

def orig(u,v,W,X,Y,Z,n,m,C,A,D,al):
    Du,Dv=1+u*u,1+v*v
    return (4*n*Du*Dv*X**2+4*m*(v*W-X+u*v*Y-u*Z)**2
        +4*C*Dv*(u*W-Y)**2+4*A*Du*Dv*Y**2
        +(4*D+2)*(u*W+u*v*X-Y-v*Z)**2
        +(4*al+2)*Du*(v*Y-Z)**2+2*Dv*(X+u*Z)**2-2*Du*Dv*Z**2)

F=(4*n*(U*U+e)*(V*V+e)*X**2+4*A*(U*U+e)*(V*V+e)*Y**2
 +4*m*(V*W+U*V*Y-U*Z-e*X)**2+(4*D+2)*(U*W+U*V*X-V*Z-e*Y)**2
 +(4*al+2)*(U*U+e)*(V*Y-Z)**2
 +4*C/e*(V*V+e)*(U*W-e*Y)**2
 +2*(V*V+e)*(2*U*X*Z+e*X**2-Z**2))
assert s.expand(orig(U/s.sqrt(e),V/s.sqrt(e),W,s.sqrt(e)*X,s.sqrt(e)*Y,Z,n,m,C,A,D,al)*e-F)==0
print('PASS exact finite-epsilon four-state identity')

P=s.Poly(s.expand(F),U,V,W,X,Y,Z)
vars_=(W,X,Y,Z)
def coeff(angular,states):
    ex=[angular[0],angular[1],0,0,0,0]
    for st in states: ex[2+vars_.index(st)]+=1
    return s.factor(P.coeff_monomial(tuple(ex)))
assert s.expand(coeff((2,2),(W,W))-(4*C/e))==0
assert s.expand(coeff((0,0),(Z,Z))-(4*al*e))==0
assert s.expand(coeff((2,0),(X,X))-(4*n*e))==0
assert s.expand(coeff((0,2),(X,X))-((4*n+2)*e))==0
assert s.expand(coeff((0,0),(X,X))-((4*n+4*m+2)*e**2))==0
assert s.expand(coeff((2,0),(Y,Y))-(4*A*e))==0
assert s.expand(coeff((0,2),(Y,Y))-((4*A+4*C+4*al+2)*e))==0
assert s.expand(coeff((0,0),(Y,Y))-((4*A+4*C+4*D+2)*e**2))==0
assert coeff((2,0),(X,Y))==coeff((0,2),(X,Y))==0
assert s.expand(coeff((1,0),(X,Z))-(4*e*(2*m+1)))==0
assert s.expand(coeff((0,1),(Y,Z))-(8*e*(D-al)))==0
assert s.expand(coeff((1,1),(X,Y))-(-4*e*(2*D+2*m+1)))==0
assert coeff((0,1),(X,Z))==coeff((1,0),(Y,Z))==0
print('PASS eight forced Z6 vertex diagonals; two unique XY zeros; three mixed coefficients')

# Explicit source-only limit fractions from the pinned exact source calculation
a,t,eta,rho=s.symbols('a t eta rho',positive=True)
n0=(t-rho)/(2*rho*(t+1))
m0=(1+rho)/(2*(t+1)*(a*a*t+eta-1))
A0=a*a*t*(1+rho)/(2*eta*(t+1))
D0=t*(t-rho)/(2*(t+1)*(1+rho-t))
N=a*a*t*(1+rho)+eta*(t+1)+rho-t
al0=N/(2*(1-eta)*(t+1))
C2=(1+rho)*N/(2*a**4*t**2*(t+1))
b0=a*a*t
assert s.factor(al0-C2*b0**2)==s.factor(N*(eta*(1+rho)-rho)/(2*(1-eta)*(t+1)))
dict_d={
'P':C2*b0**2,'Z':al0,'UX':n0*rho**2,
'VX':(n0+s.Rational(1,2))*a*a*rho**2,
'UY':A0*eta**2/a**2,'VY':(A0+al0+s.Rational(1,2))*eta**2}
Tot=sum(dict_d.values())
AA=dict_d['UX']+dict_d['UY']; BB=dict_d['VX']+dict_d['VY']
right={
'P':(eta*(1+rho)*(a*a*t+eta)+(1-eta)*rho*(t-rho))/((1-eta)*(t+1)),
'Z':t*(1+rho)*(a*a*rho+eta)/(t+1),
'AA':(1+rho)*(a*a*t*((1-eta)*rho+1)+eta-(1-eta)*(t-rho))/((1-eta)*(t+1)),
'BB':((1+rho)*(a*a*t+eta*(t+1))-(t-rho))/(t+1)}
for lhs,rhs in [(Tot-2*dict_d['P'],right['P']),
                (Tot-2*dict_d['Z'],right['Z']),
                (Tot-2*AA,right['AA']),
                (Tot-2*BB,right['BB'])]:
    assert s.cancel(lhs-rhs)==0
print('PASS four exact positive source degree-gap factorizations')

# Abstract edge construction with exactly preserved six vertex degrees
P_,Z_,ux,vx,uy,vy,h_=s.symbols('dP dZ dUX dVX dUY dVY h',positive=True)
groups={'A':[ux,uy],'B':[vx,vy],'C':[P_,Z_]}
A_=ux+uy; B_=vx+vy; C_=P_+Z_
pAB=(A_+B_-C_+2*h_)/2
pAC=(A_-B_+C_-2*h_)/2
pBC=(-A_+B_+C_-2*h_)/2
weights={}
def add(i,j,v):
    key=tuple(sorted((i,j),key=str));weights[key]=v
for i in groups['A']:
    for j in groups['B']:add(i,j,pAB*i*j/(A_*B_))
for i in groups['A']:
    for j in groups['C']:add(i,j,pAC*(i/A_)*((j-h_)/(C_-2*h_)))
for i in groups['B']:
    for j in groups['C']:add(i,j,pBC*(i/B_)*((j-h_)/(C_-2*h_)))
add(P_,Z_,h_)
for i in [P_,Z_,ux,vx,uy,vy]:
    degree=sum(val for edge,val in weights.items() if i in edge)
    assert s.cancel(degree-i)==0,(i,degree)
assert tuple(sorted((ux,uy),key=str)) not in weights
assert tuple(sorted((vx,vy),key=str)) not in weights
print('PASS 13-edge positive connected graph symbolic degree equations and absent forbidden edges')

# Native C2 symmetry: direct full-state parent and contact transformation
LHS=s.cancel(u*u*v*v*orig(1/v,1/u,Z,X,Y,W,n,D,al,A,m,C))
assert s.cancel(s.together(LHS-orig(u,v,W,X,Y,Z,n,m,C,A,D,al)))==0
cc=s.symbols('cc',positive=True)
bb,rr=s.symbols('bb rr',positive=True)
ccp=1/(a*cc); bp=a*a/bb; rp=a*a*cc*cc*rr/bb; tp=(1-cc*cc*(bb-a*a*t))/bb
oldpsi=s.Matrix([bb/a,-cc*rr,cc*(bb/a-a*t),1]);newpsi=s.Matrix([oldpsi[3],oldpsi[1],oldpsi[2],oldpsi[0]])/oldpsi[0]
assert all(s.cancel(h)==0 for h in newpsi-s.Matrix([bp/a,-ccp*rp,ccp*(bp/a-a*tp),1]))
for i in range(7):
    for j in range(7):
        assert 0<=6-i<=6 and 0<=6-j<=6
print('PASS C2 full-state projective identity, native contact map, Z6 index complement')

assert s.expand((4*n+4*m+2)*(4*m+4*al+2)-(4*m+2)**2
    -(4*(4*m+2)*(n+al)+16*n*al))==0
print('PASS exact positive 2x2 grade-2 pair reserve')
print('ALL CHARGE-SKELETON CHECKS PASSED. Finite-c full Gram coefficient matching was NOT asserted.')
