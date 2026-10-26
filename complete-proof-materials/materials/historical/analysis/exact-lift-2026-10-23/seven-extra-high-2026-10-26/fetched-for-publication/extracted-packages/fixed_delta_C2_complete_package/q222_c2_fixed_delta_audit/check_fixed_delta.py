#!/usr/bin/env python3
"""Exact finite checks for the fixed Delta angular library.

Only one rational native C2 regression source is evaluated. There is NO
SDP optimization, numerical sweep, admission claim, or dual separator.
"""
from fractions import Fraction as F
from collections import defaultdict
import json
import sympy as sp

checks = {}

# (I) Sharp weighted-Motzkin lemma.
x,y,z,A0,B0=sp.symbols('x y z A B', positive=True)
M=x**4*y**2+x**2*y**4+z**6-3*x**2*y**2*z**2
r1=x*y*(2*z**2-x**2-y**2)
r2=x*z*(z**2-y**2)
r3=y*z*(z**2-x**2)
checks['isotropic_Motzkin_three_squares'] = sp.expand((x*x+y*y)*M-r1**2-r2**2-r3**2)==0
checks['anisotropy_principal_minor'] = sp.factor(A0*B0-(A0+B0)**2/sp.Integer(4)) == -(A0-B0)**2/sp.Integer(4)

# (II) Exact projective reciprocal identities and four-square decomposition.
u,v=sp.symbols('u v', nonzero=True)
D=(u*u+v*v)*(1+u*u*v*v)
checks['Delta_four_monomial_squares'] = sp.expand(D-u**2-v**2-u**4*v**2-u**2*v**4)==0
checks['Delta_reciprocity_u'] = sp.cancel(D.subs(u,1/u)-D/u**4)==0
checks['Delta_reciprocity_v'] = sp.cancel(D.subs(v,1/v)-D/v**4)==0
checks['Delta_exchange'] = sp.expand(D-D.xreplace({u:v,v:u}))==0

# (III) Source: exact rational strict nonmaximal C2 benchmark.
a,b,c,t,r=F(13,600),F(101,1000),F(47,15),F(27,35),F(13,154)
z2=c*c;ac=1+z2*r;h=b-a*a*t
g=(1+z2)*b-(1-r+a*a*z2*t);jr=1+z2*r-t
Delta=((1+z2)*z2*b*b+2*z2*(1+z2)*b*r+(1+z2)*(1+a*a*z2)*b
       +z2*z2*(1-a*a)*r*r+2*z2*(1+a*a*z2)*r+1+a*a*z2)
N=(1+z2)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z2)*(1-r)*t
Xi=(b+r)*ac*(ac+a*a*z2*t)
d=F(1)
n=(Xi-r*Delta)/(2*r*Delta)
m=(1+a*a*z2)*(1-r)*ac*(ac+a*a*z2*t)/(2*g*Delta)
C=ac*N/(2*a*a*t*Delta)
AA=a*a*ac*(b+r)*(ac+a*a*z2*t)/(2*h*Delta)
DD=(1+z2)*(b+r)*(b-a*a*z2*r)*(ac+a*a*z2*t)/(2*a*a*jr*Delta)
al=z2*z2*(b+r)*N/(2*(1-z2*h)*Delta)
checks['strict_native_source'] = all(s>0 for s in
   (a,b,c,t,r,1-r,h,g,jr,1-z2*h,b-a*a*z2*r,Delta,N,n,m,C,AA,DD,al))

# Full SG2 parent /4 (NO omitted state/normal/angle terms).
P=defaultdict(F)
def packet(weight, terms, shifts=((0,0),)):
    for su,sv in shifts:
        for i,j,l,v1 in terms:
            for k,hk,ll,v2 in terms:
                ell, mm=sorted((l,ll))
                P[(i+k+su,j+hk+sv,ell,mm)]+=weight*v1*v2
DU=((0,0),(2,0));DV=((0,0),(0,2));BOTH=((0,0),(2,0),(0,2),(2,2))
packet(n,[(0,0,1,1)],BOTH)
packet(m,[(0,1,0,1),(0,0,1,-1),(1,1,2,1),(1,0,3,-1)])
packet(C,[(1,0,0,1),(0,0,2,-1)],DV)
packet(AA,[(0,0,2,1)],BOTH)
packet(DD+F(1,2),[(1,0,0,1),(1,1,1,1),(0,0,2,-1),(0,1,3,-1)])
packet(al+F(1,2),[(0,1,2,1),(0,0,3,-1)],DU)
packet(F(1,2),[(0,0,1,1),(1,0,3,1)],DV)
packet(-F(1,2),[(0,0,3,1)],BOTH)
P={k:v for k,v in P.items() if v}
checks['parent_bidegree_2_2']=all(i<=2 and j<=2 for i,j,l,k in P)

uc,vc=c,a*c
psi=(b/a,-c*r,c*h/a,F(1))
def E(poly):
    return sum(w*uc**i*vc**j*psi[l]*psi[k] for (i,j,l,k),w in poly.items())
checks['contact_parent_zero']=E(P)==0
angle_grad=[]
for angle in (0,1):
    ang=F(0)
    for (i,j,l,k),coeff in P.items():
        exponent=(i,j)[angle]
        if exponent:
            ang+=exponent*coeff*uc**(i-(angle==0))*vc**(j-(angle==1))*psi[l]*psi[k]
    angle_grad.append(ang)
checks['contact_two_angular_gradients_zero']=angle_grad==[0,0]
state_grad=[]
for ell in range(4):
    row=F(0)
    for (i,j,l,k),w in P.items():
        if ell==l:row+=w*uc**i*vc**j*psi[k]
        if ell==k:row+=w*uc**i*vc**j*psi[l]
    state_grad.append(row)
checks['contact_four_state_gradients_zero']=state_grad==[0]*4

# Delta*F/4, as a polynomial in u,v and ALL four state variables.
mult=((2,0),(0,2),(4,2),(2,4))
T=defaultdict(F)
for (i,j,l,k),co in P.items():
    for di,dj in mult:T[(i+di,j+dj,l,k)]+=co
T={k:v for k,v in T.items() if v}
checks['multiplied_parent_bidegree_6_6'] = all(i<=6 and j<=6 for i,j,l,k in T)

# Exact diagonal supports: Newton vertices, universal at STRICT contacts.
vertices=[
   {(0,4),(4,0),(6,2),(6,4),(4,6),(2,6)},
   {(0,2),(2,0),(4,0),(6,2),(6,4),(4,6),(2,6),(0,4)},
   {(0,2),(2,0),(4,0),(6,2),(6,4),(4,6),(2,6),(0,4)},
   {(0,2),(2,6),(6,2),(4,0),(2,0),(0,4)}
]
# A half-Newton monomial r^i s^j psi_l occurs only inside
# these convex integer polygons: 4<=i+j<=10 (W), 2..10(X,Y),
# 2..8(Z), plus |i-j|<=4; 0<=i,j<=6.
allowed=[]
for l,(lo,hi) in enumerate(((4,10),(2,10),(2,10),(2,8))):
    candidates={(i,j) for i in range(7) for j in range(7)
                if lo<=i+j<=hi and abs(i-j)<=4}
    allowed.append(candidates)
checks['Newton_feature_counts']=[len(t0) for t0 in allowed]==[30,37,37,30]
parity=[sum(1 for l in range(4) for i,j in allowed[l]
            if (i%2,j%2)==p) for p in ((0,0),(0,1),(1,0),(1,1))]
checks['Newton_parity_orders']=parity==[44,28,28,34]
checks['Newton_total']=sum(parity)==134
checks['native_diagonal_support_inside_hulls']=all(
    (i,j) in allowed[ll] for i,j,ll,k in T if ll==k
)
# Every (i,j) exponent of target polynomial is an element of twice
# the feature Newton polygon. This checks the explicit polygon constraints.
checks['diagonal_Newton_vertices_present']=all(
    vertices[l].issubset({(i,j) for i,j,ll,k in T if ll==k==l})
    for l in range(4)
)

checks['regression_contact_all_components_nonzero']=all(x!=0 for x in psi)
checks['source_contact_positive_angles']=uc>0 and vc>0

result={'all_pass':all(checks.values()),'checks':checks,'parity_orders':parity,
        'compressed_orders':[x-1 for x in parity],
        'regression':{'native_parameters':list(map(str,[a,b,c,t,r])),
                      'contact_angles':list(map(str,[uc,vc])),
                      'contact_state':list(map(str,psi)),
                      'source_weights_quarter_parent':dict(zip(('d','n','m','C','A','D','alpha'),map(str,[d,n,m,C,AA,DD,al])))},
        'scope':'Exact source and scalar/corner identities; no numerical or exact SDP feasibility claim.'}
from pathlib import Path
Path(__file__).with_name('verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'all_pass':result['all_pass'],'checks':checks,'parity_orders':parity},indent=2))
if not result['all_pass']:raise SystemExit(1)
