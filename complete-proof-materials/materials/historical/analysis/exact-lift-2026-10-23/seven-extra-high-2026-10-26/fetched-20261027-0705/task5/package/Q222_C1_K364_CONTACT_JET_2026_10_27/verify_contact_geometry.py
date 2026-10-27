#!/usr/bin/env python3
"""Pure exact algebra: full C1 Alice block, parity, contact skew, jet freedom.
No sources sampled, no SDP solve, no rank reconstruction of the K364 map.
Analytic positivity assumptions (native Hessian and Calderon) proved in PROOF.md.
"""
import sympy as s
u,v=s.symbols('u v',real=True)
w0,w1,y0,y1=s.symbols('w0 w1 y0 y1',real=True,nonzero=True)
p0,p1,q0=s.symbols('p0 p1 q0',real=True)
d,m,C,A,D,al,be=s.symbols('d m C A D alpha beta',real=True)
W,X,Y,Z=s.symbols('W X Y Z',real=True)
ps=s.Matrix([W,X,Y,Z]); x=ps.col_join(u*ps)
Du,Dv=1+u*u,1+v*v
F=(4*m*(v*W-X+u*v*Y-u*Z)**2+4*C*Dv*(u*W-Y)**2
 +4*A*Du*Dv*Y**2+(4*D+2*d)*(u*W+u*v*X-Y-v*Z)**2
 +(4*al+2*d)*Du*(v*Y-Z)**2+(4*be+2*d)*Dv*(X+u*Z)**2-2*d*Du*Dv*Z**2)
Q=s.hessian(F,(W,X,Y,Z))/2
P0=Q.subs(u,0); P1=Q.diff(u).subs(u,0); P2=Q.diff(u,2).subs(u,0)/2
assert (Q-P0-u*P1-u*u*P2).applyfunc(s.expand)==s.zeros(4)
T=s.Matrix(2,2,s.symbols('T00 T01 T10 T11'))
S=s.zeros(4)
S[:2,2:]=T;S[2:,:2]=-T.T
H=(Dv*P0).row_join(Dv*P1/2+S).col_join((Dv*P1/2-S).row_join(Dv*P2))
Ia=[0,1,6,7];Ib=[2,3,4,5]
assert H.extract(Ia,Ib)==s.zeros(4)
assert (x.T*H*x)[0].expand()==s.expand(Dv*F)
print('PASS exact unrestricted four-state parent, parity-separated 4+4 pencil')

ww=s.Matrix([w0,w1]); yy=s.Matrix([y0,y1]);pp=s.Matrix([p0,p1]);
# Choose q1 from exact source-only compatibility w.p=y.q; y1=Zc=1 in native chart.
qq=s.Matrix([q0,(ww.dot(pp)-y0*q0)/y1])
Sw=ww.dot(ww);Sy=yy.dot(yy)
Sstar=pp*yy.T/Sy+ww*qq.T/Sw-(ww.dot(pp))*ww*yy.T/(Sw*Sy)
wperp=s.Matrix([-w1,w0]);yperp=s.Matrix([-y1,y0]);Nfree=wperp*yperp.T
assert (Sstar*yy-pp).applyfunc(s.cancel)==s.zeros(2,1)
assert (Sstar.T*ww-qq).applyfunc(s.cancel)==s.zeros(2,1)
assert Nfree*yy==s.zeros(2,1) and Nfree.T*ww==s.zeros(2,1)
# Full affine contact fiber is Sstar+tau*w_perp*y_perp^T.
print('PASS closed source-rational affine line of all contact-killing 2x2 skews')

uc=s.symbols('uc',real=True,positive=True)
xA=s.Matrix([w0,w1,0,0,0,0,uc*y0,uc*y1]);xB=s.Matrix([0,0,y0,y1,uc*w0,uc*w1,0,0]);xp=xA+xB;xm=xA-xB
V=s.zeros(4);V[:2,2:]=T;V[2:,:2]=-T.T
dH=s.zeros(8);dH[:4,4:]=V;dH[4:,:4]=-V
# Skew changes vanish on every rank-one Alice-state tensor; all tangent first variations vanish.
assert s.expand((xp.T*dH*xp)[0])==0
pstate=s.Matrix([w0,w1,y0,y1]);
Etangent=[s.eye(4)[:,i].col_join(uc*s.eye(4)[:,i]) for i in range(4)]+[s.zeros(4,1).col_join(pstate)]
assert all(s.expand((t.T*dH*xp)[0])==0 for t in Etangent)
val=s.expand(2*uc*(ww.T*T*yy)[0])
assert s.expand((xA.T*dH*xA)[0]-val)==0
assert s.expand((xB.T*dH*xB)[0]+val)==0
assert s.expand((xm.T*dH*xp)[0]-2*val)==0
print('PASS skew derivatives annihilate all five physical tangent directions; split the two soft pivots oppositely')

# The normal-direction gauge image has dimension three; fixing its soft split leaves exactly two.
M=s.Matrix.vstack(s.Matrix.hstack(T*yy,T.T*ww)) if False else s.Matrix([
    [y0,y1,0,0], [0,0,y0,y1], [w0,0,w1,0], [0,w0,0,w1]])
assert M.det()==0
assert s.factor(M.extract([0,1,3],[0,1,3]).det())==-w0*y0*y1
J1=wperp*yy.T; J2=ww*yperp.T
assert s.expand((ww.T*J1*yy)[0])==0 and s.expand((ww.T*J2*yy)[0])==0
assert (J1*yy-Sy*wperp).applyfunc(s.expand)==s.zeros(2,1) and (J1.T*ww).applyfunc(s.expand)==s.zeros(2,1)
assert (J2*yy).applyfunc(s.expand)==s.zeros(2,1) and (J2.T*ww-Sw*yperp).applyfunc(s.expand)==s.zeros(2,1)
Tbase=ww*yy.T/(Sw*Sy)
assert s.cancel((ww.T*Tbase*yy)[0]-1)==0
print('PASS rank-three skew/tangent action; two independent source-only soft-neutral jet directions')
print('PASS source-only charge direction shifts the two Bob-contact quadratic Schur coefficients by opposite 2uc units')
print('SCOPE: exact algebra only. The Calderon anchor and native strict-Hessian estimates are proved separately.')
