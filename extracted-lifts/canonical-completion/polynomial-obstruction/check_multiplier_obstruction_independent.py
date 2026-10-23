"""Independent interfaces not tested by the supplied obstruction checker."""
import json
from pathlib import Path
import sympy as s

u,v,U=s.symbols('u v U')
d,la,n,m,C,A,D,al,be=s.symbols('d la n m C A D al be')
W,X,Y,Z=s.symbols('W X Y Z')
Du,Dv=1+u*u,1+v*v
F=(4*la*Du*(W+v*X)**2+4*n*Du*Dv*X**2
   +4*m*(v*W-X+u*v*Y-u*Z)**2+4*C*Dv*(u*W-Y)**2
   +4*A*Du*Dv*Y**2+(4*D+2*d)*(u*W+u*v*X-Y-v*Z)**2
   +(4*al+2*d)*Du*(v*Y-Z)**2+(4*be+2*d)*Dv*(X+u*Z)**2
   -2*d*Du*Dv*Z**2)
P=s.Matrix([[4*m*v*v,-4*m*v],[-4*m*v,4*m+(4*be+2*d)*Dv]])
N=s.Matrix([[4*C*Dv+4*D+2*d,(4*D+2*d)*v],[(4*D+2*d)*v,(4*D+2*d)*v*v]])
E=s.Matrix([[4*A*Dv+(4*al+2*d)*v*v,-(4*al+2*d)*v],[-(4*al+2*d)*v,4*al-2*d*v*v]])
K=4*la*s.Matrix([[1,v],[v,v*v]])+4*n*Dv*s.diag(0,1)
a=s.Matrix(s.symbols('a1 a2')); b=s.Matrix(s.symbols('b1 b2'))
q=s.symbols('q')
J=s.Matrix(2,2,s.symbols('j11 j12 j21 j22'))
G=(q*(P+K)).row_join(J).col_join(J.T.row_join(q*(N+K)))
H=(q*(P+E)).row_join(q*(K-E)-J).col_join((q*(K-E)-J.T).row_join(q*(N+E)))
Phi=((1+U)*(a.T*P*a)[0]+U*(1+U)*(b.T*N*b)[0]
     +U*((a-b).T*E*(a-b))[0]+((a+U*b).T*K*(a+U*b))[0])
apub=a.col_join(U*b); apb=a.col_join(b)
assert s.expand(q*Phi-(apub.T*G*apub)[0]-U*(apb.T*H*apb)[0])==0
L=s.Matrix([[u,0,u*u,0],[0,u,0,u*u],[u,0,-1,0],[0,u,0,-1]])
state=s.Matrix([W,X,Y,Z]); transformed=L*state
sub=dict(zip(list(a)+list(b),transformed)); sub[U]=u*u
assert s.expand(Phi.subs(sub,simultaneous=True)-u*u*Du*F)==0

# Derive the standard Bell matrix from Pauli tensors, rather than taking
# the displayed matrix's physical interpretation as an input.
c,e,r,z0,T=s.symbols('c e r z0 T')
sz=s.diag(1,-1); sx=s.Matrix([[0,1],[1,0]]); I=s.eye(2)
At=c*sz+r*sx; Bt=e*sz+z0*sx; kron=s.kronecker_product
bell=(kron(sz,sz)+kron(sz,Bt)+kron(At,sz)-kron(At,Bt)
      +(1-T)*(kron(sz,I)+kron(I,sz)))
Omega=1+c+e-c*e
printed=s.Matrix([[Omega+2*(1-T),(1-c)*z0,(1-e)*r,-r*z0],
 [(1-c)*z0,-Omega,-r*z0,(e-1)*r],
 [(1-e)*r,-r*z0,-Omega,(c-1)*z0],
 [-r*z0,(e-1)*r,(c-1)*z0,Omega-2*(1-T)]])
assert all(s.expand(a)==0 for a in bell-printed)

# Independently establish the uniform testing degree used in SOS closedness.
t,x,y=s.symbols('t x y')
psit=[1+s.Rational(2,3)*t**4*x*x*y*y,-s.Rational(1,3)*t**3*x*x*y,
      -t*x+s.Rational(2,3)*t**3*x*y*y,-t*t*x*y-s.Rational(2,3)*t**4*x**3*y]
family={d:2*T-4*T**3,la:T**3,n:1-T+T**3,m:T**3,C:(1-T)/2+T**3,
        A:T**3,D:1-T+T**3,al:1-T+T**3,be:(1-T)/2+T**3}
Ft=s.expand(F.subs(family).subs(T,t*t).subs(dict(zip([u,v,W,X,Y,Z],[-t*x,t*y]+psit)),simultaneous=True))
assert s.Poly(Ft,x,y).total_degree()==12
assert s.Poly(t*t*x*x*(1+t*t*x*x)*Ft,x,y).total_degree()==16
pt=s.Poly(Ft,t)
M=x**4*y*y+x*x*y**4-3*x*x*y*y+1
assert all(pt.coeff_monomial(t**j)==0 for j in range(6))
assert s.expand(pt.coeff_monomial(t**6)-s.Rational(4,3)*(M+2))==0
assert all(power[0]%2==0 for power,coef in pt.terms())

result={'generic_two_LMI_identity':True,'generic_full_parent_pullback':True,
        'physical_Bell_matrix_derived_from_Pauli_operators':True,
        'critical_full_parent_degree_xy':12,
        'critical_parent_times_u2_Du_degree_xy':16,
        'critical_leading_coefficient_exact':True,
        'bounded_degree_consequence':'q of degree 2r gives total degree at most 2r+16, hence SOS root degree at most r+8',
        'scope':'Exact identity and degree checks only; positivity and SOS closure arguments are reviewed analytically.'}
Path(__file__).with_name('multiplier_obstruction_independent_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
