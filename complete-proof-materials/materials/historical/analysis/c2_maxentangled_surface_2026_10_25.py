"""One exact two-parameter native contact calculation; no solver or sweep."""
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix

c,r,u,v=s.symbols('c r u v',positive=True)
dc,nc,mc,Cc,Ac,Dc,alc=s.symbols('d n m C A D alpha')
active=[dc,nc,mc,Cc,Ac,Dc,alc]
I=s.eye(2); Z=s.diag(1,-1); X=s.Matrix([[0,1],[1,0]])
a0=s.kronecker_product(Z,I)
a1=s.kronecker_product(((1-u*u)*Z+2*u*X)/(1+u*u),I)
b0=s.kronecker_product(I,Z)
b1=s.kronecker_product(I,((1-v*v)*Z+2*v*X)/(1+v*v))
ops=[s.eye(4),a0,a1,b0,b1,a0*b0,a0*b1,a1*b0,a1*b1]
coeff=[dc+nc+mc+Cc+Ac+Dc+alc,nc-Ac-alc,mc-Cc-Dc,
       -nc+Cc+Ac,-mc+Dc-alc,-dc/2-nc-Ac,dc/2+alc,
       -dc/2-Cc,-dc/2-mc-Dc]
F=s.zeros(4)
for a,o in zip(coeff,ops): F+=a*o
psi=s.Matrix([1,-c*r,c*r,1]); sub={u:c,v:c}
equations=list(F.subs(sub)*psi)+[(psi.T*s.diff(F,t).subs(sub)*psi)[0] for t in (u,v)]
M,_=s.linear_eq_to_matrix([s.factor(a) for a in equations],active)
N=M[:,1:]; rhs=-M[:,0]
aug=N.row_join(rhs).applyfunc(s.cancel)
red,den,piv=DomainMatrix.from_Matrix(aug).to_field().rref_den()
R=red.to_Matrix()/den.as_expr()
assert piv==tuple(range(6))
answer=s.Matrix([1]+[s.factor(R[i,6]) for i in range(6)])
assert all(s.factor(a)==0 for a in M*answer)
z=s.symbols('z',positive=True)
source=[s.factor(a.subs(c*c,z)) for a in answer]
assert all(c not in a.free_symbols for a in source)
minor=s.factor(DomainMatrix.from_Matrix(N.applyfunc(s.cancel)).det().as_expr())
data={'active_source_order':list(map(str,active)),
      'source':list(map(str,source)),
      'contact_minor_delete_d':str(minor),
      'all_six_contact_equations_verified':True}
x=s.symbols('x',positive=True)
K0=r+x+2*r*x
scale=4*K0*(1-r)*(1-x)
polynomial=[scale,(1-r-x-3*r*x)*(1-r)*(1-x),
            (1-r)**2*(1-x*x),4*r*r*(1-x*x),
            (1-r*r)*(1-x*x),(1-r*r)*(1-x)**2,4*x*x*(1-r*r)]
assert all(s.factor(a.subs(z,x/r)*scale-b)==0 for a,b in zip(source,polynomial))
sd,sn,sm,sC,sA,sD,sal=polynomial
assert s.factor(sd-4*sA+2*sm+2*sD)==0
assert s.factor(sn-sm-sD+sA)==0
assert s.factor((sA-sm)**2-sm*sC)==0
assert s.factor((sA-sD)**2-sD*sal)==0
# Exact inverse/minimal-source decomposition for arbitrary cone coordinates.
aa,mm,dd=s.symbols('A m D',positive=True)
ri=(aa-mm)/(aa+mm); xi=(aa-dd)/(aa+dd)
lam=aa/((1-ri*ri)*(1-xi*xi))
minimal=[4*aa-2*mm-2*dd,mm+dd-aa,mm,(aa-mm)**2/mm,
         aa,dd,(aa-dd)**2/dd]
assert all(s.factor(a.subs({r:ri,x:xi},simultaneous=True)*lam-b)==0
           for a,b in zip(polynomial,minimal))
domain_inverse=s.factor(ri+xi+3*ri*xi-1)
assert s.factor(domain_inverse-4*aa*(aa-mm-dd)/((aa+mm)*(aa+dd)))==0
# The two independent excess event rays really are in the surface closure.
eps=s.symbols('eps',positive=True)
axisC=[s.limit(a.subs({r:1-eps,x:eps/8},simultaneous=True),eps,0,dir='+') for a in polynomial]
axisAlpha=[s.limit(a.subs({r:eps/8,x:1-eps},simultaneous=True),eps,0,dir='+') for a in polynomial]
assert axisC==[0,0,0,4,0,0,0]
assert axisAlpha==[0,0,0,0,0,0,4]
# Whole-parent physical validity: one exact four-variable determinant.
du=1+u*u; dv=1+v*v
packets=[(4*sn*du*dv,[0,1,0,0]),(4*sm,[v,-1,u*v,-u]),
         (4*sC*dv,[u,0,-1,0]),(4*sA*du*dv,[0,0,1,0]),
         (4*sD+2*sd,[u,u*v,-1,-v]),((4*sal+2*sd)*du,[0,0,v,-1]),
         (2*sd*dv,[0,1,0,u]),(-2*sd*du*dv,[0,0,0,1])]
P=s.zeros(4)
for coefficient,vector in packets:
    a=s.Matrix(vector); P+=coefficient*a*a.T
G=16*(1-r*r)**3*(1-x*x)**3
wa=8*r*x*(1-r); wb=8*r*x*(1-x)
wc=(1+r)*(1+x)*(1-r-x-3*r*x); wt=8*r*x*(r+x)
Xn,Yn=s.symbols('Xn Yn',positive=True)
Rnorm=(wa*(Xn+1/Xn-2)+wb*(Yn+1/Yn-2)
       +wc*(Xn/Yn+Yn/Xn-2)+wt*(Yn/Xn+Xn+1/Yn-3))
Rtest=s.cancel(Rnorm.subs({Xn:r*u*u/x,Yn:r*v*v/x}))
expected=s.cancel(16*G*du**2*dv**2*u*u*v*v*Rtest)
print('Checking one full four-variable polynomial determinant...',flush=True)
detP=DomainMatrix.from_Matrix(P.applyfunc(s.expand)).det().as_expr()
assert s.Poly(detP-expected,r,x,u,v).is_zero
data.update({
    'new_coordinate':'x=r*z',
    'domain':'0<r<1, 0<x, r+x+3*r*x<=1',
    'positive_polynomial_source_scale':str(scale),
    'polynomial_source':list(map(lambda a:str(s.expand(a)),polynomial)),
    'bidegrees':[(int(s.degree(a,r)),int(s.degree(a,x))) for a in polynomial],
    'linear_and_two_rank_one_PSD_identities_verified':True,
    'inverse_parameters':{'r':str(ri),'x':str(xi)},
    'inverse_domain_identity':str(domain_inverse),
    'all_minimal_source_decomposition_identities_verified':True,
    'event_C_limit':list(map(str,axisC)),
    'event_alpha_limit':list(map(str,axisAlpha)),
    'AMGM_weights':list(map(str,[wa,wb,wc,wt])),
    'determinant_common_factor':str(G),
    'full_four_variable_polynomial_determinant_verified':True,
    'exact_closed_source_cone_total_PSD_order':7,
})
Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
print(json.dumps(data,indent=2))
