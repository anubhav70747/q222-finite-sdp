"""Exact determinant and elementary AM--GM proof for one C2 curve.

No solver, sampling, or imported physical-validity lemma.
Requires SymPy; emits exact factorized data next to this file.
"""
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix

z,u,v=s.symbols('z u v',real=True)
K=8*(2-z)*(1+2*z)
d=K; la=0; be=0
n=(2-z)*(2-5*z); m=(2-z)*(2+z)
C=4*(2-z)*(2+z); A=3*(2-z)*(2+z)
D=3*(2-z)**2; al=12*z*z
du=1+u*u; dv=1+v*v
packets=[(4*n*du*dv,[0,1,0,0]),
         (4*m,[v,-1,u*v,-u]),
         (4*C*dv,[u,0,-1,0]),
         (4*A*du*dv,[0,0,1,0]),
         (4*D+2*d,[u,u*v,-1,-v]),
         ((4*al+2*d)*du,[0,0,v,-1]),
         (2*d*dv,[0,1,0,u]),
         (-2*d*du*dv,[0,0,0,1])]
P=s.zeros(4)
for coefficient,vector in packets:
    a=s.Matrix(vector); P+=coefficient*a*a.T
P=P.applyfunc(s.expand)
t=[None,d+2*(la+n+D),d+2*(al+m+A),d+2*(la+C+D),
   d+2*(al+be+m),d+2*(la+be+n),d+2*(al+C+A),
   d+2*(be+n+m),d+2*(C+A+D)]
q={(0,1):4*al*m*t[5]*t[8],(0,2):4*m*D*t[5]*t[6],
   (1,0):4*al*A*t[7]*t[3],(1,2):4*C*D*t[5]*t[2],
   (2,0):4*n*A*t[3]*t[4],(2,1):4*C*A*t[1]*t[4]}
q={e:s.factor(a) for e,a in q.items()}
h=s.factor(P.subs({u:1,v:1}).det()/256-sum(q.values()))
U,V=s.symbols('U V',positive=True)
L=h+sum(a*U**(i-1)*V**(j-1) for (i,j),a in q.items())
w={e:s.factor(a*z**(e[0]+e[1]-2)) for e,a in q.items()}
delta=s.factor(w[0,2]-w[2,0])
assert s.factor(w[2,1]-w[0,1]-delta)==0
assert s.factor(w[1,0]-w[1,2]-delta)==0
assert s.factor(h+sum(w.values()))==0
X,Y=s.symbols('X Y',positive=True)
Lnorm=(w[0,1]*(X+1/X-2)+w[1,2]*(Y+1/Y-2)
       +w[2,0]*(X/Y+Y/X-2)+delta*(Y/X+X+1/Y-3))
assert s.factor(L.subs({U:z*X,V:z*Y})-Lnorm)==0
print('Computing one exact polynomial determinant...',flush=True)
detP=DomainMatrix.from_Matrix(P).det().as_expr()
expected=s.cancel(16*du**2*dv**2*u**2*v**2*L.subs({U:u*u,V:v*v}))
assert s.Poly(detP-expected,u,v,z).is_zero
common=432*(2-z)**3*(z+2)**3
normalized_q={(0,1):16*z*z,(0,2):(6-z)*(2-z),(1,0):48*z*z,
              (1,2):16*(2-z),(2,0):3*(z+2)*(2-5*z),(2,1):16*(z+2)}
assert all(s.factor(q[e]-common*a)==0 for e,a in normalized_q.items())
assert s.factor(h-common*(14*z*z-96*z-24))==0
# Six contact equations, on all seven active source coordinates.
dc,nc,mc,Cc,Ac,Dc,alc=s.symbols('d n m C A D alpha')
active=[dc,nc,mc,Cc,Ac,Dc,alc]
c=s.symbols('c',positive=True)
I=s.eye(2); Z=s.diag(1,-1); Xp=s.Matrix([[0,1],[1,0]])
a0=s.kronecker_product(Z,I)
a1=s.kronecker_product(((1-u*u)*Z+2*u*Xp)/(1+u*u),I)
b0=s.kronecker_product(I,Z)
b1=s.kronecker_product(I,((1-v*v)*Z+2*v*Xp)/(1+v*v))
ops=[s.eye(4),a0,a1,b0,b1,a0*b0,a0*b1,a1*b0,a1*b1]
coeff=[dc+nc+mc+Cc+Ac+Dc+alc,nc-Ac-alc,mc-Cc-Dc,
       -nc+Cc+Ac,-mc+Dc-alc,-dc/2-nc-Ac,dc/2+alc,
       -dc/2-Cc,-dc/2-mc-Dc]
F=s.zeros(4)
for a,o in zip(coeff,ops): F+=a*o
psi=s.Matrix([1,-c/2,c/2,1]); sub={u:c,v:c}
equations=list(F.subs(sub)*psi)+[(psi.T*s.diff(F,t).subs(sub)*psi)[0] for t in (u,v)]
M,_=s.linear_eq_to_matrix([s.factor(a) for a in equations],active)
minor=s.factor(M[:,:6].det())
assert minor==-768*c**10/(1+c*c)**5
curve=s.Matrix([1,n/K,m/K,C/K,A/K,D/K,al/K]).subs(z,c*c)
assert all(s.factor(a)==0 for a in M*curve)
data={
    'positive_source_multiplier':str(s.factor(K)),
    'scaled_source_order_d_lambda_n_m_C_A_D_alpha_beta':list(map(str,[d,la,n,m,C,A,D,al,be])),
    'scaled_t':{str(i):str(s.factor(t[i])) for i in range(1,9)},
    'scaled_q':{str(e):str(a) for e,a in q.items()},
    'scaled_h':str(h),
    'scaled_contact_weights':{str(e):str(a) for e,a in w.items()},
    'triangle_weight_scaled':str(delta),
    'unscaled_pair_horizontal':str(s.factor(w[0,1]/K**4)),
    'unscaled_pair_vertical':str(s.factor(w[1,2]/K**4)),
    'unscaled_pair_diagonal':str(s.factor(w[2,0]/K**4)),
    'unscaled_triangle_weight':str(s.factor(delta/K**4)),
    'full_polynomial_determinant_identity_verified':True,
    'three_balance_and_value_identities_verified':True,
    'four_AMGM_terms_identity_verified':True,
    'normalized_q_by_common_factor':{str(e):str(s.factor(a)) for e,a in normalized_q.items()},
    'normalized_h':'14*z**2-96*z-24',
    'common_determinant_factor':str(common),
    'contact_rank_six_minor_delete_alpha':str(minor),
    'all_six_contact_equations_verified':True,
}
Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
print(json.dumps(data,indent=2))
