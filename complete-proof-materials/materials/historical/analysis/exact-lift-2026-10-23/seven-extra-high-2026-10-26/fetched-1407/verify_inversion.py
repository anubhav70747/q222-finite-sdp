#!/usr/bin/env python3
"""Exact native C2 contact/source inversion and full-parent covariance, symbolic."""
import sympy as s
z,a,b,t,r=s.symbols('z a b t r',nonzero=True)

def source(z,a,b,t,r):
 ac=1+z*r;h=b-a*a*t;g=(1+z)*b-(1-r+a*a*z*t);j=1+z*r-t
 dl=(1+z)*z*b*b+2*z*(1+z)*b*r+(1+z)*(1+a*a*z)*b+z*z*(1-a*a)*r*r+2*z*(1+a*a*z)*r+1+a*a*z
 nn=(1+z)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z)*(1-r)*t
 return [((b+r)*ac*(ac+a*a*z*t)-r*dl)/(2*r*dl),
 (1+a*a*z)*(1-r)*ac*(ac+a*a*z*t)/(2*g*dl),
 ac*nn/(2*a*a*t*dl),a*a*ac*(b+r)*(ac+a*a*z*t)/(2*h*dl),
 (1+z)*(b+r)*(b-a*a*z*r)*(ac+a*a*z*t)/(2*a*a*j*dl),
 z*z*(b+r)*nn/(2*(1-z*h)*dl)]
orig=source(z,a,b,t,r)
zp=1/(a*a*z); bp=a*a/b;rp=a*a*z*r/b;tp=(1-z*(b-a*a*t))/b
new=source(zp,a,bp,tp,rp)
permutation=(0,4,5,3,1,2)  # (n,m,C,A,D,alpha) -> (n,D,alpha,A,m,C)
for i,j in enumerate(permutation):
 assert s.cancel(new[i]-orig[j])==0,(i,j)
 print('VERIFIED native source weight',i,'maps to',j,flush=True)
# A parent-level homogeneous identity for arbitrary, independent test angles/state.
u,v,W,X,Y,Z=s.symbols('u v W X Y Z',nonzero=True)
n,m,C,A,D,alpha=s.symbols('n m C A D alpha')
def parent(u,v,W,X,Y,Z,n,m,C,A,D,alpha):
 du=1+u*u;dv=1+v*v
 return (4*n*du*dv*X**2+4*m*(v*W-X+u*v*Y-u*Z)**2+4*C*dv*(u*W-Y)**2+4*A*du*dv*Y**2
 +(4*D+2)*(u*W+u*v*X-Y-v*Z)**2+(4*alpha+2)*du*(v*Y-Z)**2+2*dv*(X+u*Z)**2-2*du*dv*Z**2)
assert s.cancel(u*u*v*v*parent(1/v,1/u,Z,X,Y,W,n,D,alpha,A,m,C)-parent(u,v,W,X,Y,Z,n,m,C,A,D,alpha))==0
print('VERIFIED full parent identity, source map and fixed Z6 index reversal',flush=True)
