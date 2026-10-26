"""Exact obstruction and missing angular relation on m=D,C=alpha."""
import json
from pathlib import Path
import sympy as s

d,n,m,C,A,u,v=s.symbols('d n m C A u v',real=True)
Du=1+u*u; Dv=1+v*v
packets=[(4*n*Du*Dv,[0,1,0,0]),(4*m,[v,-1,u*v,-u]),
         (4*C*Dv,[u,0,-1,0]),(4*A*Du*Dv,[0,0,1,0]),
         (4*m+2*d,[u,u*v,-1,-v]),((4*C+2*d)*Du,[0,0,v,-1]),
         (2*d*Dv,[0,1,0,u]),(-2*d*Du*Dv,[0,0,0,1])]
P=s.zeros(4)
for c,L in packets:
    q=s.Matrix(L); P+=c*q*q.T
P=P.applyfunc(s.expand)
P1=P.subs({u:1,v:1})
T=s.Matrix([[1,0,0,1],[0,1,0,0],[0,0,1,0],[1,0,0,-1]])
B=s.Matrix([[2*C,d/2,-2*C-d/2],
            [d/2,3*d/4+m+2*n,-d/4-m],
            [-2*C-d/2,-d/4-m,2*A+2*C+3*d/4+m]])
assert (T.T*P1*T/8-s.diag(B,2*C+d+4*m)).applyfunc(s.expand)==s.zeros(4)
E=s.Matrix([[C+d/2,d/4+m,-d/4-m],
            [d/4+m,-d/4,0],[-d/4-m,0,-C-d/4]])
assert (T.T*(s.diff(P,u).subs({u:1,v:1})-P1)*T/8)[:3,:3]==E
assert (T.T*(s.diff(P,v).subs({u:1,v:1})-P1)*T/8)[:3,:3]==-E

r,y=s.symbols('r y',positive=True)
kernel=s.Matrix([1,-r,y,1])
kernel_source={C:d*(r+y)/(4*(1-y)),
 n:(d*(2-3*r-y)-4*m*(r+y))/(8*r),
 A:(d*(2+r-y)-4*m*(r+y))/(8*y)}
assert (P1.subs(kernel_source)*kernel).applyfunc(s.factor)==s.zeros(4,1)
gap=d*(1-r)*(2+r-y)-8*m*(r+y)
Fu=s.factor((kernel.T*s.diff(P,u)*kernel)[0].subs({u:1,v:1}).subs(kernel_source))
Fv=s.factor((kernel.T*s.diff(P,v)*kernel)[0].subs({u:1,v:1}).subs(kernel_source))
assert s.factor(Fu-2*gap)==0
assert s.factor(Fv+2*gap)==0
stationary_m=d*(1-r)*(2+r-y)/(8*(r+y))
stationary=[d,
 d*((1+r)*(2+r-y)-8*r)/(16*r),
 stationary_m,d*(r+y)/(4*(1-y)),
 d*(1+r)*(2+r-y)/(16*y)]
assert all(s.factor(kernel_source.get(q,q).subs(m,stationary_m)-val)==0
           for q,val in zip([d,n,m,C,A],stationary))
assert s.factor(gap.subs(m,stationary_m))==0

# No additional homogeneous linear source equation holds on this surface.
polys=[s.Poly(s.cancel(q.subs(d,1)*16*r*y*(1-y)*(r+y)),r,y)
       for q in stationary]
mons=sorted(set().union(*(set(q.monoms()) for q in polys)))
coefficient_matrix=s.Matrix([[q.coeff_monomial(e) for q in polys] for e in mons])
assert coefficient_matrix.rank()==5

example={d:1,n:s.Rational(1,4),m:s.Rational(1,4),C:s.Rational(1,6),A:s.Rational(3,4)}
psi=s.Matrix([1,-s.Rational(1,4),s.Rational(1,4),1])
L=s.Matrix([[1,0,0],[s.Rational(3,2),1,0],[-s.Rational(5,2),1,1]])
assert B.subs(example)==L*s.diag(s.Rational(1,3),s.Rational(3,4),0)*L.T
assert (2*C+d+4*m).subs(example)==s.Rational(7,3)
F=s.factor((psi.T*P.subs(example)*psi)[0])
assert s.factor(F.subs(v,1)-(u-1)*(67*u-55)/12)==0
assert F.subs({u:s.Rational(7,8),v:1})==-s.Rational(29,768)
assert s.diff(F,u).subs({u:1,v:1})==1
assert s.diff(F,v).subs({u:1,v:1})==-1

data={
 'slice':'lambda=beta=0, D=m, alpha=C',
 'one_point_block_congruence_exact':True,
 'normalized_angle_differentials_on_W_equals_Z_are_opposite':True,
 'kernel_equations_and_remaining_angular_equation_exact':True,
 'missing_relation':'8*m*(r+s)=d*(1-r)*(2+r-s)',
 'homogeneous_kernel_relation':'d*(w+X)*(2*w-X-Y)=8*m*w*(Y-X)',
 'stationary_surface_source_linear_span_rank':5,
 'no_additional_homogeneous_linear_source_equation':True,
 'counterexample_source_d_lambda_n_m_C_A_D_alpha_beta':[
     '1','0','1/4','1/4','1/6','3/4','1/4','1/6','0'],
 'counterexample_P1_PSD_LDL_pivots':['1/3','3/4','0','7/3'],
 'counterexample_kernel':['1','-1/4','1/4','1'],
 'angle_derivatives':['1','-1'],
 'negative_testing_pair':['7/8','1'],
 'negative_full_parent_numerator':'-29/768',
 'scope':'Disproves one-point PSD validity test on this linear slice; does not disprove another finite lift.',
 'no_solver_or_numerical_sweep':True}
Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
print(json.dumps(data,indent=2))
