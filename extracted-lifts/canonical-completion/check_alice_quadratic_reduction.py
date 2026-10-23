#!/usr/bin/env python3
"""New exact all-source Alice quadratic reduction; no old data recomputed."""
import sympy as s
u,v=s.symbols('u v', real=True)
lam,n,m,C,A,D,al,be,d=s.symbols('lambda n m C A D alpha beta d', nonnegative=True)
W,X,Y,Z=s.symbols('W X Y Z')
S,T,P,Q=s.symbols('S T P Q')
Du,Dv=1+u*u,1+v*v
F=4*lam*Du*(W+v*X)**2+4*n*Du*Dv*X**2+4*m*(v*W-X+u*v*Y-u*Z)**2+4*C*Dv*(u*W-Y)**2+4*A*Du*Dv*Y**2+(4*D+2*d)*(u*W+u*v*X-Y-v*Z)**2+(4*al+2*d)*Du*(v*Y-Z)**2+(4*be+2*d)*Dv*(X+u*Z)**2-2*d*Du*Dv*Z**2
qp=4*m*(v*S-T)**2+(4*be+2*d)*Dv*T*T
qm=4*C*Dv*P*P+(4*D+2*d)*(P+v*Q)**2
q0=4*A*Dv*(S-P)**2+(4*al+2*d)*(v*(S-P)-(T-Q))**2-2*d*Dv*(T-Q)**2
U=s.symbols('U')
Phi=(1+U)*qp+U*(1+U)*qm+U*q0+4*lam*((S+v*T)+U*(P+v*Q))**2+4*n*Dv*(T+U*Q)**2
trans={W:S+u*u*P,X:T+u*u*Q,Y:u*(S-P),Z:u*(T-Q)}
assert s.expand(F.subs(trans, simultaneous=True)-Du*Phi.subs(U,u*u))==0
N=s.Matrix([[1,0,u*u,0],[0,1,0,u*u],[u,0,-u,0],[0,u,0,-u]])
L=s.Matrix([[u,0,u*u,0],[0,u,0,u*u],[u,0,-1,0],[0,u,0,-1]])
assert s.simplify(N.det()-u*u*Du*Du)==0
assert all(s.expand(e)==0 for e in L*N-u*Du*s.eye(4))
vec=s.Matrix([W,X,Y,Z])
fixed=s.expand(Phi.subs(U,u*u).subs(dict(zip((S,T,P,Q),L*vec)), simultaneous=True)-u*u*Du*F)
assert fixed==0
coeffs=[s.hessian(s.expand(Phi).coeff(U,k),(S,T,P,Q))/2 for k in range(3)]
assert all(s.Poly(e,v).degree()<=2 for mat in coeffs for e in mat)
Ap=s.Matrix([[4*m*v*v,-4*m*v],[-4*m*v,4*m+(4*be+2*d)*Dv]])
Am=s.Matrix([[4*C*Dv+4*D+2*d,(4*D+2*d)*v],[(4*D+2*d)*v,(4*D+2*d)*v*v]])
C0=s.Matrix([[4*A*Dv+(4*al+2*d)*v*v,-(4*al+2*d)*v],[-(4*al+2*d)*v,4*al-2*d*v*v]])
K=4*lam*s.Matrix([[1,v],[v,v*v]])+4*n*Dv*s.diag(0,1)
y11,y12,y21,y22=s.symbols('y11 y12 y21 y22')
J=s.Matrix([[y11,y12],[y21,y22]])
G=s.BlockMatrix([[Ap+K,J],[J.T,Am+K]]).as_explicit()
H=s.BlockMatrix([[Ap+C0,K-C0-J],[K-C0-J.T,Am+C0]]).as_explicit()
av=s.Matrix([S,T]);bv=s.Matrix([P,Q]); ab=av.col_join(bv)
aub=av.col_join(U*bv)
assert s.expand((aub.T*G*aub)[0]+U*(ab.T*H*ab)[0]-Phi)==0
print('PASS: full nine-source polynomial identity; all four states independent')
print('PASS: pure-testing inverse and fixed multiplier u^2(1+u^2)')
print('PASS: quadratic U coefficients have Bob degree at most two')
print('PASS: two 4x4 LMI reconstruction for four independent selector entries')
print('No source weight was inverted; no previous matrices read or regenerated.')
