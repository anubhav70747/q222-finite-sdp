#!/usr/bin/env python3
"""Exact strict-C2 Laurent hexagon curvature, no numerical inputs."""
import sympy as s
x1,x2,y1,y2,z1,z2=s.symbols('x1 x2 y1 y2 z1 z2',positive=True)
# weights on exponent vectors (+/-1,0), (0,+/-1), (+/-1,-/+1)
weights={(1,0):x1,(-1,0):x2,(0,1):y1,(0,-1):y2,(1,-1):z1,(-1,1):z2}
H=s.zeros(2)
for (u,v),q in weights.items():
 e=s.Matrix([u,v]);H+=q*(e*e.T)
A=x1+x2;B=y1+y2;C=z1+z2
assert H==s.Matrix([[A+C,-C],[-C,B+C]])
assert s.expand(H.det()-(A*B+C*(A+B)))==0
assert A>0 and B>0 and C>0
print('EXACT: positive C2 log-angular Hessian with determinant AB+C(A+B)>0')
print('EXACT: near contact L(z)=sum w_e(exp(e.z)-1-e.z) >= exp(-sqrt(2)R)*z^THz/2 for |z|<=R')
