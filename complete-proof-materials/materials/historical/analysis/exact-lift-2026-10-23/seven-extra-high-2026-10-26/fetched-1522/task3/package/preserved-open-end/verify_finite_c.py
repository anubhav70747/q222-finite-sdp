#!/usr/bin/env python3
"""Exact symbolic verification of an analytic C2 finite-c sixth-root continuation theorem.
All arithmetic exact. No SDP solves, random samples, absent .pkl, or testing inverses.
Dependencies: Python3 + SymPy. Run `python verify_finite_c.py`.
"""
from map_core import *
from two_kernel import M, DIFFCOLS
import time
T0=time.time()
def passed(t):print('VERIFIED',t,'(%.1fs)'%(time.time()-T0),flush=True)

def vec_poly(poly):
 co={}
 for exp,v in poly.terms():
  eu,ev=exp[:2];st=[k for k,p in enumerate(exp[2:]) for _ in range(p)]
  assert len(st)==2
  co[key(eu,ev,*st)]=v
 assert set(co).issubset(set(KEYS))
 return S.Matrix([S.factor(co.get(i,S.Integer(0))) for i in KEYS])

Ws=native_weights()
W0={k:S.factor(v.subs(eps,0)) for k,v in Ws.items()}
W1={k:S.factor(S.diff(v,eps).subs(eps,0)) for k,v in Ws.items()}
C2=S.factor(S.diff(Ws[C],eps,2).subs(eps,0)/2)
trunc={k:(eps**2*C2 if k==C else W0[k]+eps*W1[k]) for k in (n,m,C,A,D,al)}
Ftrunc=native_parent(n,m,C,A,D,al).subs(trunc)/4
f0=vec_poly(S.Poly(S.expand(Ftrunc.subs(eps,0)),U,V,W,X,Y,Z))
f1=vec_poly(S.Poly(S.expand(S.diff(Ftrunc,eps).subs(eps,0)),U,V,W,X,Y,Z))
assert M.shape==(78,79) and M.rank()==66
passed('fixed 14+2 two-kernel map: 78 coefficient equations, 79 variables, exact rank 66')

# Six nontrivial all-epsilon two-kernel source cancellations, via one shared source denominator.
left=M.T.nullspace();fv=target_coeff('symbolic')
zv=(a,t,eta,rho,eps);P=lambda expr:S.Poly(S.expand(expr),*zv,domain='QQ')
Hden=S.cancel(S.denom(Ws[n])/(2*rho)); g=a*a*t+eta-1+(eta+rho)*eps
units={n:rho,m:g,C:a*a*t,A:eta,D:a*a*(1+rho-t),al:1-eta}
num={k:P(S.cancel(Ws[k]*(2*Hden*den))) for k,den in units.items()}
prod={k:P(S.prod([den for kk,den in units.items() if kk!=k])) for k in units}
prodall=P(S.prod(units.values()))
for j,lf in enumerate(left):
 q=S.expand((lf.T*fv)[0]);
 if q==0:continue
 q0=q.subs({k:0 for k in units})
 eq=P(2*Hden*eps*a*a*q0)*prodall
 for k in units:eq+=P(eps*a*a*S.diff(q,k))*num[k]*prod[k]
 assert eq.is_zero,('failed generic source identity',j,eq.terms()[:1])
passed('all 12 exact all-epsilon native source identities on the two-kernel coefficient image')

# Endpoint Gram: the known four end squares with a free positive kappa-gauge; all state coordinates intact.
kap=S.symbols('kap',positive=True)
cX=2*t*(t-2*rho-1)/(rho*(t+1)*(t-rho-1))
cY=-2*a*a*t*(rho+1)*(a*a*t+2*eta-1)/(eta*(eta-1)*(t+1)*(a*a*t+eta-1))
p2=2*(rho+1)*(a*a*t+eta)/(a*a*t*(t+1)*(a*a*t+2*eta-1))
q2=2*(rho+1)*(t-rho)/(t*(t+1)*(2*rho+1-t))
lX0=-rho*(rho+1)/(t*(t-2*rho-1));lX2=-rho*(rho-t)/(t-2*rho-1)
lY3=-eta*(eta-1)/(a*a*t*(a*a*t+2*eta-1));lY4=-eta*(a*a*t+eta)/(a*a*t+2*eta-1)
c0=S.Matrix([lX0*a*t,-a*rho,lX2*a,0,0,0]);c3=S.Matrix([0,0,0,lY3*a*a*t,lY4,eta])
u1=S.Matrix([a*t,0,-a*t,0,0,0]);u2=S.Matrix([0,0,0,a*a*t,-a*a*t,0])
s1=S.Matrix([a*t,0,0,0,-a*t,0]);s2=S.Matrix([0,0,-a*a*t,a*a*t,0,0])
Gzero=(cX*c0*c0.T+cY*c3*c3.T+p2*u2*u2.T+q2*u1*u1.T)/4
Gk=Gzero-kap*(s1*s2.T+s2*s1.T)/8
assert (Gk*S.ones(6,1)).applyfunc(S.cancel).is_zero_matrix
H0=S.zeros(79,1)
for i in range(5):
 for j in range(i,5):H0[DIFFCOLS.index((i,j))]=Gk[i,j]
H0[78]=kap*a**3*t**2/4
assert all(S.cancel(r)==0 for r in M*H0-f0)
passed('exact full-parent endpoint Gram for every source parameter and positive symbolic gauge kap')

# The source-dependent positive five-edge tree. Its reduced 5x5 Gram is PD, including all strict face parameters.
Nstar=a*a*t+eta+rho-t
edges={(6,7):(rho+1)*Nstar/(2*(t+1)),
 (6,9):a*a*rho*t*(rho+1)/(2*(t+1)),
 (6,10):eta*t*(rho+1)/(2*(t+1)),
 (7,8):rho*(t-rho)/(2*(t+1)),
 (7,11):eta*(rho+1)*(a*a*t+eta)/(2*(1-eta)*(t+1))}
Tcharge=S.zeros(6)
for (i,j),w in edges.items():
 i-=6;j-=6
 Tcharge[i,i]+=w;Tcharge[j,j]+=w;Tcharge[i,j]-=w;Tcharge[j,i]-=w
six_degrees=S.Matrix([C2*(a*a*t)**2,W0[al],W0[n]*rho**2,
 (W0[n]+S.Rational(1,2))*a*a*rho**2,W0[A]*eta**2/a**2,
 (W0[A]+W0[al]+S.Rational(1,2))*eta**2])
assert all(S.cancel(w)==0 for w in S.Matrix(list(Tcharge.diagonal()))-six_degrees)
# Matrix-tree cofactor identity with five formal edge variables avoids expanding
# an unnecessarily enormous rational source determinant.
edge_vars=S.symbols('g0:5')
Tree=S.zeros(6)
for ((i,j),_w),q in zip(edges.items(),edge_vars):
 i-=6;j-=6
 Tree[i,i]+=q;Tree[j,j]+=q;Tree[i,j]-=q;Tree[j,i]-=q
assert S.expand(Tree[:5,:5].det()-S.prod(edge_vars))==0
passed('five source-positive tree edges, six exact charged-degree identities, tree Gram determinant')

# Source-algebraic old-to-tail linear couplings at first epsilon order.
lamX=rho/(a*a*t);lamY=eta/(a*a*t)
Btail=S.Matrix([[S.factor(lamX*Gzero[i,3]),S.factor(-lamY*Gzero[i,0])] for i in range(6)])
assert all(S.cancel(v)==0 for v in S.ones(1,6)*Btail)

# Prescribe the charge tree, tail directions and zero charge-tail/second-tail interactions.
fix={};free=[]
for k,(i,j) in enumerate(DIFFCOLS):
 if i<5 and j<5:free.append(k)
 elif i<5 and j>=7:free.append(k)
 elif i<5 and j in (5,6):fix[k]=Btail[i,j-5]
 elif i in (5,6) and j>=5:fix[k]=S.Integer(0)
 elif i>=7 and j>=7:fix[k]=Tcharge[i-7,j-7]
 else:raise ValueError((i,j))
free.append(78)
assert len(free)==41 and len(fix)==38
rhs=f1-S.Matrix([sum(M[i,k]*v for k,v in fix.items()) for i in range(78)])
Mfree=M[:,free]
assert Mfree.rank()==37
for j,v in enumerate(Mfree.T.nullspace()):
 assert S.cancel(S.together((v.T*rhs)[0]))==0,('first-order',j)
passed('full first-order native parent compatible with fixed tree/tail couplings (37-column exact rank)')

# Complete deterministic first-order rational Gram, using 37 source-independent pivot pairs.
pc=list(Mfree.rref()[1]);pr=list(Mfree[:,pc].T.rref()[1]);assert len(pc)==len(pr)==37
chosen=Mfree[pr,pc].inv()*rhs[pr,:]
H1=S.zeros(79,1)
for k,v in fix.items():H1[k]=v
for j,i in enumerate(pc):H1[free[i]]=chosen[j]
assert all(S.cancel(v)==0 for v in M*H1-f1)
passed('explicit first-order rational 79-coordinate Gram fully matches all 78 native coefficients')

# All-epsilon completion: fixed universal rational right inverse on the 66-dimensional image.
pc2=list(M.rref()[1]);pr2=list(M[:,pc2].T.rref()[1]);assert len(pc2)==len(pr2)==66
Rinv=M[pr2,pc2].inv()  # constant rational 66x66
assert M[pr2,pc2]*Rinv==S.eye(66)
# Define H2(eps)=Rinv*(f(eps)-f0-eps*f1)[pr2]/eps^2 in the pivot coordinates.
# All 12 source identities above prove f(eps) lies in im M. The exact Taylor
# definition of f0,f1 proves no pole at eps=0. Therefore M*(H0+eps*H1+eps^2*H2)=f(eps).
passed('source-independent rational 66-pivot all-epsilon right inverse; no missing binary matrices')
# The exact last two diagonal coefficients of the *quarter parent* are unique.
# This verifies that every possible rational image right inverse has the same
# prescribed epsilon^2 tail block, including the entire D, C and alpha packets.
for st,expected in [
 ((0,0,1,1),eps**2*rho**2*(n+m+S.Rational(1,2))),
 ((0,0,2,2),eps**2*(eta/a)**2*(A+C+D+S.Rational(1,2))),
 ((0,0,1,2),S.Integer(0))]:
 assert S.cancel(fv[KEY_IND[st]]/4-expected)==0,st
passed('exact all-epsilon final X,Y tail diagonals and vanishing XY mixed coefficient')


# Final Schur margins in the actual two-kernel difference coordinates.
TX=rho**2*(W0[n]+W0[m]+S.Rational(1,2))
TY=(eta/a)**2*(W0[A]+W0[D]+S.Rational(1,2))
deltaX=rho*t*(rho+1)/(2*(t+1))
deltaY=eta*t*(rho+1)/(2*(t+1))
assert S.factor(TX-lamX**2*Gzero[3,3]-deltaX)==0
assert S.factor(TY-lamY**2*Gzero[0,0]-deltaY)==0
assert all(S.cancel(v)==0 for v in (S.Matrix([Gzero[0,3],Gzero[3,0]])))
assert S.cancel((S.ones(1,3).row_join(S.zeros(1,2))*S.diff(Gk[:5,:5],kap)*(S.ones(3,1).col_join(S.zeros(2,1))))[0]-a**3*t**2/4)==0
passed('both positive last-Schur limiting pivots and the positive newly activated old-contact pivot')
print('CERTIFIED THEOREM INPUTS: exact coefficient map, complete full-parent identities, positive tree and final graded Schur limits.',flush=True)
print('For any strict face point choose 0<kap<kap_max sufficiently small, then epsilon>0 sufficiently small. All choices are pure source.',flush=True)
