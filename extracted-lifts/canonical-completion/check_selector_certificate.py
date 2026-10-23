#!/usr/bin/env python3
"""Independent standard-library-only verifier of the NEW selector data.

Uses Fraction throughout.  Does not import the builder or read any old pencil.
"""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
def zeros(n,m=None): return [[Q(0) for _ in range(n if m is None else m)] for _ in range(n)]
def add(A,B): return [[a+b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def scale(k,A): return [[k*a for a in ar] for ar in A]
def transpose(A): return [list(r) for r in zip(*A)]
def mv(A,x): return [sum(a*b for a,b in zip(ar,x)) for ar in A]
def block(A,B,C,D): return [a+b for a,b in zip(A,B)]+[c+d for c,d in zip(C,D)]
def determinant(A):
    if len(A)==1: return A[0][0]
    return sum((-1)**j*A[0][j]*determinant([r[:j]+r[j+1:] for r in A[1:]]) for j in range(len(A)))
def kernel_rank(A):
    R=[r[:] for r in A]; nr=len(R);nc=len(R[0]);pivot=0
    for j in range(nc):
        k=next((k for k in range(pivot,nr) if R[k][j]),None)
        if k is None: continue
        R[k],R[pivot]=R[pivot],R[k]
        val=R[pivot][j];R[pivot]=[x/val for x in R[pivot]]
        for k in range(nr):
            if k!=pivot:
                val=R[k][j];R[k]=[x-val*y for x,y in zip(R[k],R[pivot])]
        pivot+=1
        if pivot==nr:break
    return pivot

def elementary(c,v):
    d,la,n,m,C,A,D,al,be=c; dv=1+v*v
    ap=[[4*m*v*v,-4*m*v],[-4*m*v,4*m+(4*be+2*d)*dv]]
    am=[[4*C*dv+4*D+2*d,(4*D+2*d)*v],[(4*D+2*d)*v,(4*D+2*d)*v*v]]
    h=[[4*A*dv+(4*al+2*d)*v*v,-(4*al+2*d)*v],[-(4*al+2*d)*v,4*al-2*d*v*v]]
    k=[[4*la,4*la*v],[4*la*v,4*la*v*v+4*n*dv]]
    return ap,am,h,k

def matrices(c,v,J):
    ap,am,h,k=elementary(c,v)
    G=block(add(ap,k),J,transpose(J),add(am,k))
    off=add(add(k,scale(-1,h)),scale(-1,J))
    H=block(add(ap,h),off,transpose(off),add(am,h))
    L=block(G,zeros(4),zeros(4),H)
    return G,H,L

def eval_coo(data,parameters,v):
    L=zeros(8)
    for k,j,row,col,z in data['entries']:
        assert 0<=k<len(parameters) and 0<=j<=2 and 0<=row<8 and 0<=col<8
        L[row][col]+=parameters[k]*(v**j)*Q(z)
    return L

full=json.loads((HERE/'POINTWISE_SELECTOR_ORDER8.coo.json').read_text())
c1=json.loads((HERE/'C1_POLYNOMIAL_SELECTOR3_ORDER8.coo.json').read_text())
for data,kind in [(full,'full'),(c1,'c1')]:
    assert data['order']==8 and data['blocks']==[4,4] and data['max_testing_degree']==2
    count=len(data['parameter_order'])
    for k in range(count):
        param=[Q(int(j==k)) for j in range(count)]
        # Three distinct exact values determine each quadratic coefficient.
        for v in map(Q,[-1,0,1]):
            if kind=='full':
                c=param[:9]; J=[param[9:11],param[11:13]]
            else:
                d,m,C,A,D,al,be,p,r,t=param
                c=[d,Q(0),Q(0),m,C,A,D,al,be]
                J=[[p*v*v,Q(0)],[r*v,t*v*v]]
            expected=matrices(c,v,J)[2]
            assert eval_coo(data,param,v)==expected
            assert transpose(expected)==expected
print('PASS: every exported integer matrix coefficient, independently reconstructed')

cert=json.loads((HERE/'POLYNOMIAL_SELECTOR_COUNTEREXAMPLE.json').read_text())
z=dict(zip(cert['source_order'],map(Q,cert['source'])))
c=[z['d'],Q(0),Q(0),z['m'],z['C'],z['A'],z['D'],z['alpha'],z['beta']]
u=Q(cert['contact']['u']);v=Q(cert['contact']['v']);U=u*u;du=1+U
W,X,Y,Z=map(Q,cert['contact']['state'])
a=[(W+u*Y)/du,(X+u*Z)/du]
b=[(W-Y/u)/du,(X-Z/u)/du]
assert a==list(map(Q,cert['transformed_contact']['a']))
assert b==list(map(Q,cert['transformed_contact']['b']))
ap,am,h,k=elementary(c,v)
P=block(add(scale(1+U,ap),scale(U,h)),scale(-U,h),scale(-U,h),
        add(scale(U*(1+U),am),scale(U,h)))
assert mv(P,a+b)==[0]*4

forced=cert['forced_selector'];p,r,t=(Q(forced[k]) for k in ['p','r','s'])
J=[[p*v*v,Q(0)],[r*v,t*v*v]]
assert mv(J,b)==[-q/U for q in mv(ap,a)]
assert mv(transpose(J),a)==[-U*q for q in mv(am,b)]
coef=[]
for q in [(1,0,0),(0,1,0),(0,0,1)]:
    pp,rr,tt=map(Q,q); JJ=[[pp*v*v,Q(0)],[rr*v,tt*v*v]]
    coef.append(mv(JJ,b)+mv(transpose(JJ),a))
coef=transpose(coef)
rhs=[-q/U for q in mv(ap,a)]+[-U*q for q in mv(am,b)]
assert coef==[[Q(x) for x in row] for row in cert['kernel_equations']['coefficient_matrix']]
assert rhs==list(map(Q,cert['kernel_equations']['right_hand_side']))
assert kernel_rank(coef)==3
assert mv(coef,[p,r,t])==rhs
print('PASS: exact full-state contact kernel and uniquely forced three-parameter selector')

v=Q(1);G,H,L=matrices(c,v,[[p*v*v,Q(0)],[r*v,t*v*v]])
ind=[1,2,3]; minor=determinant([[G[i][j] for j in ind] for i in ind])
assert minor==Q(cert['negative_principal_minor']['value'])<0
aa,bb,cc,ee=4*z['m'],4*z['beta']+2*z['d'],4*z['C'],4*z['D']+2*z['d']
assert 3<ee<4 and 0<cc<Q(1,100) and 0<aa+2*bb<5 and r<Q(5,4) and t>2
assert minor==2*ee*cc*(aa+2*bb)-ee*(r-t)**2-2*cc*t*t
assert minor<Q(-103,80)
print('PASS: exact negative principal minor at v=1; hand bound < -103/80')

fac=Q(cert['determinant']['positive_prefactor']);assert fac>0
# Both sides are degree at most eight, so nine values verify the identity.
for v in map(Q,range(9)):
    G=matrices(c,v,[[p*v*v,Q(0)],[r*v,t*v*v]])[0]
    exact=fac*v**4*(v-16)*(v+16)*(1949306657149583773*v*v-44177987181070080)
    assert determinant(G)==exact
print('PASS: full determinant factorization by nine exact evaluations')
print('SCOPE: all polynomial Bob selectors are excluded for this established valid source.')
print('This does not exclude an additional fixed testing multiplier or algebraic testing family.')
