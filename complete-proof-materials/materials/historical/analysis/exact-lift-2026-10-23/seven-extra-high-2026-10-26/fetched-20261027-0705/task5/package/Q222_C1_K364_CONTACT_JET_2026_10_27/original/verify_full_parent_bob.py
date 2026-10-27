#!/usr/bin/env python3
"""Independent exact symbolic checks of C1 Bob-harmonic block construction.
No SDP, random/sampled source, numerical sweep or source-rank recomputation.
"""
import itertools
import sympy as S

u,v=S.symbols('u v',real=True)
W,X,Y,Z=S.symbols('W X Y Z',real=True)
d,m,C,A,D,al,be=S.symbols('d m C A D alpha beta',positive=True)
psi=S.Matrix([W,X,Y,Z]); states=(W,X,Y,Z)
Du=1+u*u; Dv=1+v*v
F=S.expand(4*m*(v*W-X+u*v*Y-u*Z)**2+4*C*Dv*(u*W-Y)**2
 +4*A*Du*Dv*Y**2+(4*D+2*d)*(u*W+u*v*X-Y-v*Z)**2
 +(4*al+2*d)*Du*(v*Y-Z)**2+(4*be+2*d)*Dv*(X+u*Z)**2
 -2*d*Du*Dv*Z**2)
Q=S.hessian(F,states)/2
P0=Q.subs(u,0)
P1=Q.diff(u).subs(u,0)
P2=Q.diff(u,2).subs(u,0)/2
assert (Q-P0-u*P1-u*u*P2).applyfunc(S.expand)==S.zeros(4)
print('PASS unrestricted C1 full-state Alice degree-two extraction')

# Every expression is an identity in u,v,W,X,Y,Z and seven independent source weights.
F_u0=(4*m*(v*W-X)**2+4*(C+A)*Dv*Y**2+4*D*(Y+v*Z)**2
 +4*al*(v*Y-Z)**2+4*be*Dv*X**2+2*d*Dv*(X**2+Y**2))
F_u2=(4*m*(v*Y-Z)**2+4*C*Dv*W**2+4*A*Dv*Y**2
 +4*D*(W+v*X)**2+4*al*(v*Y-Z)**2+4*be*Dv*Z**2
 +2*d*((W+v*X)**2+(v*Y-Z)**2))
F_v0=((4*C+4*D+2*d)*(u*W-Y)**2
 +(4*be+2*d+4*m)*(X+u*Z)**2+4*A*Du*Y**2+4*al*Du*Z**2)
F_v2=(4*m*(W+u*Y)**2+4*C*(u*W-Y)**2+4*(A+al)*Du*Y**2
 +4*D*(u*X-Z)**2+4*be*(X+u*Z)**2+2*d*Du*(X**2+Y**2))
for name,res in [('u0',F.subs(u,0)-F_u0),('u2',S.Poly(F,u).coeff_monomial(u*u)-F_u2),
 ('v0',F.subs(v,0)-F_v0),('v2',S.Poly(F,v).coeff_monomial(v*v)-F_v2)]:
 assert S.expand(res)==0,name
print('PASS four full-state positive endpoint SOS factorizations')

# Exact ordinary, not contact-restricted, endpoint 8x8 matrix Grams.
# Each row is (a,b) with the square (a.psi + u*b.psi)^2.
I=S.eye(4)
def gram_of_endpoint(terms):
 return sum((weight*(aa.col_join(bb))*(aa.col_join(bb)).T for weight,aa,bb in terms),S.zeros(8))
O=S.zeros(4,1)
def vec(**args):
 m=S.zeros(4,1)
 for name,value in args.items():m['WXYZ'.index(name),0]=value
 return m
K=4*C+4*D+2*d;L=4*be+2*d+4*m
G0=gram_of_endpoint([
 (K,vec(Y=-1),vec(W=1)), (L,vec(X=1),vec(Z=1)),
 (4*A,vec(Y=1),O),(4*A,O,vec(Y=1)),
 (4*al,vec(Z=1),O),(4*al,O,vec(Z=1))])
Ginf=gram_of_endpoint([
 (4*m,vec(W=1),vec(Y=1)),(4*C,vec(Y=-1),vec(W=1)),
 (4*(A+al),vec(Y=1),O),(4*(A+al),O,vec(Y=1)),
 (4*D,vec(Z=-1),vec(X=1)),(4*be,vec(X=1),vec(Z=1)),
 (2*d,vec(X=1),O),(2*d,O,vec(X=1)),
 (2*d,vec(Y=1),O),(2*d,O,vec(Y=1))])
P0v2=P0.applyfunc(lambda x:S.expand(x).coeff(v,2))
P1v2=P1.applyfunc(lambda x:S.expand(x).coeff(v,2))
P2v2=P2.applyfunc(lambda x:S.expand(x).coeff(v,2))
Bzero=G0[:4,4:]; Binf=Ginf[:4,4:]
assert G0[:4,:4]==P0.subs(v,0)
assert G0[4:,4:]==P2.subs(v,0)
assert Bzero+Bzero.T==P1.subs(v,0)
assert Ginf[:4,:4]==P0v2
assert Ginf[4:,4:]==P2v2
assert Binf+Binf.T==P1v2
S0=(Bzero-Bzero.T)/2; S4=(Binf-Binf.T)/2
print('PASS full 8x8 endpoint Gram identities and positive source-weight square rows')

# Free all-orders source-only skew K(v) of degree at most 4.
ss=S.symbols('sk0:30')
blocks=[]
for j in range(5):
 Sj=S.zeros(4)
 for t,(i,k) in enumerate(itertools.combinations(range(4),2)):
  Sj[i,k]=ss[6*j+t];Sj[k,i]=-ss[6*j+t]
 blocks.append(Sj)
Lsk=sum((v**j*blocks[j] for j in range(5)),S.zeros(4))
x=psi.col_join(u*psi)
Gram=(Dv*P0).row_join(Dv*P1/2+Lsk).col_join(
 (Dv*P1/2-Lsk).row_join(Dv*P2))
assert S.expand((x.T*Gram*x)[0]-Dv*F)==0
assert Gram==Gram.T
assert S.degree(Gram[0,0],v)<=4
print('PASS source-only degree-four skew changes NO coefficient of Dv times full parent')

# The test-feature spaces of matrix half-line Markov-Lukacs are strict subsets of Z_(6,12).
features=[]
for aa in (0,1):
 for bb in (0,1,2):
  for ll in range(4): features.append((6*aa,6*bb,ll,'even'))
for aa in (0,1):
 for bb in (0,1):
  for ll in range(4): features.append((6*aa,6*bb+3,ll,'odd'))
assert len(features)==40
assert len(set((i,j,l) for i,j,l,s in features))==40
assert all(0<=i<=6 and 0<=j<=12 for i,j,l,s in features)
print('PASS forty fixed non-source-dependent features embed in Z_(6,12)')

# Contact-kernel repair in abstract coefficient algebra, using the exact kernel and u jet.
uc,vc=S.symbols('uc vc', positive=True)
z=S.Matrix(S.symbols('z0:4',real=True)); P0c=S.Matrix(4,4,lambda i,j:S.symbols('A%d%d'%(min(i,j),max(i,j))))
P1c=S.Matrix(4,4,lambda i,j:S.symbols('B%d%d'%(min(i,j),max(i,j))))
P2c=S.Matrix(4,4,lambda i,j:S.symbols('C%d%d'%(min(i,j),max(i,j))))
Dc=1+vc**2
bvec=-Dc*(P0c*z/uc+P1c*z/2)
Fjet=(z.T*(P0c+uc*P1c+uc**2*P2c)*z)[0]
Ujet=(z.T*(P1c+2*uc*P2c)*z)[0]
assert S.simplify((z.T*bvec)[0]+Dc*(Fjet/uc-Ujet/2))==0
Sk=(bvec*z.T-z*bvec.T)/(z.T*z)[0]
assert (Sk.T+Sk).applyfunc(S.cancel)==S.zeros(4)
Nnorm=(z.T*z)[0]
first=S.simplify(Dc*P0c*z+uc*(Dc*P1c/2+Sk)*z)
second=S.simplify((Dc*P1c/2-Sk)*z+uc*Dc*P2c*z)
# This polynomial identity exhibits all zero-jet conditions explicitly, and no inverse of a testing variable.
assert (first+uc*z*(z.T*bvec)[0]/Nnorm).applyfunc(S.cancel)==S.zeros(4,1)
assert (second-Dc*(P0c+uc*P1c+uc**2*P2c)*z/uc-z*(z.T*bvec)[0]/Nnorm).applyfunc(S.cancel)==S.zeros(4,1)
print('PASS rational contact-kernel skew formula modulo actual state and u-stationarity jets')

print('All exact checks passed. No universal polynomial-matrix PSD assertion has been made.')
