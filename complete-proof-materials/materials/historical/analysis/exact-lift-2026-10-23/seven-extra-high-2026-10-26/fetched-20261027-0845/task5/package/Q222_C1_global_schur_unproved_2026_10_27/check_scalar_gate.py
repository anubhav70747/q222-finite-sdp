#!/usr/bin/env python3
"""Exact generic 2x2 Schur-to-four-polynomial equivalence and endpoint skew.
No source sampling, SDP, approximate tests, or claims of universal feasibility.
"""
import sympy as s
v=s.symbols('v', real=True)
d,m,C,A,D,al,be=s.symbols('d m C A D alpha beta', positive=True)
h=1+v**2
P0w=s.Matrix([[4*m*v**2,-4*m*v],[-4*m*v,4*m+(4*be+2*d)*h]])
P0y=s.Matrix([[4*(A+C)*h+4*D+2*d*h+4*al*v**2,4*(D-al)*v],[4*(D-al)*v,4*D*v**2+4*al]])
P2w=s.Matrix([[4*C*h+4*D+2*d,(4*D+2*d)*v],[(4*D+2*d)*v,(4*D+2*d)*v**2]])
P2y=s.Matrix([[4*A*h+(4*al+2*d+4*m)*v**2,-(4*al+2*d+4*m)*v],[-(4*al+2*d+4*m)*v,4*al+4*be*h+2*d+4*m]])
E=s.Matrix([[2*m*v**2-2*C*h-2*D-d,-v*(2*D+d+2*m)],[-v*(2*D+d+2*m),2*be*h-2*D*v**2+d+2*m]])
assert s.factor(P0w.det())==8*m*(2*be+d)*v**2*h
assert s.factor(P2w.det())==8*C*(2*D+d)*v**2*h
assert s.expand(P0y.det()-8*h*(2*(A+C)*(D*v**2+al)+2*D*al*h+d*(D*v**2+al)))==0
assert s.expand(P2y.det()-8*h*(2*A*(al+be*h+m)+A*d+2*al*be*v**2+be*d*v**2+2*be*m*v**2))==0
print('PASS exact positive 2x2 endpoint-block determinants for all v>0')
A0=h*P0w;B0=h*P0y;A2=h*P2y;B2=h*P2w;L=h*E
T00=2*C+2*D+d;T11=2*be+d+2*m
T0=s.diag(T00,T11)
assert (L.subs(v,0)+T0)[0,:]==s.zeros(1,2)
assert (L.subs(v,0)-T0)[1,:]==s.zeros(1,2)
print('PASS forced endpoint skew T(0)=diag(2C+2D+d, 2beta+d+2m)')
# Independent 2x2 symmetric hard blocks and arbitrary source-only skew block
b0,b1,b2,c0,c1,c2,a0,a1,a2=s.symbols('b0 b1 b2 c0 c1 c2 a0 a1 a2')
B=s.Matrix([[b0,b1],[b1,b2]])
Z=s.Matrix([[c0,c1],[c1,c2]])
H=s.Matrix([[a0,a1],[a1,a2]])
R=s.simplify(Z-H.T*B.adjugate()*H/B.det())
M=B.row_join(H).col_join(H.T.row_join(Z))
assert s.factor(M.det()-B.det()*R.det())==0
trace_numer=s.expand(B.det()*s.trace(Z)-s.trace(H.T*B.adjugate()*H))
assert s.factor(B.det()*s.trace(R)-trace_numer)==0
print('PASS exact 4x4 determinant and 2x2 Schur-trace polynomial identities')
print('THEOREM: with B(v) positive definite, M(v)>=0 iff trace_numer(v)>=0 and det M(v)>=0.')
print('The source-wide 40-feature criterion requires these two inequalities in EACH parity block for all v>0.')
print('No universal source-only quartic T(v) is constructed or certified by this checker.')
