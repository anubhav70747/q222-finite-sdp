"""Exact obstruction to replacing the fixed-angle-ratio source surface by its quadrics."""
import json
from pathlib import Path
import sympy as s

a=s.Integer(2)
d,n,m,C,A,D,h=s.symbols('d n m C A D alpha')
coords=s.Matrix([d,n,m,C,A,D,h])
direction=s.Matrix([1,-s.Rational(1,3),0,-s.Rational(1,2),s.Rational(1,3),0,-s.Rational(1,2)])
base=s.Matrix([1,s.Rational(43,299),s.Rational(77,1495),s.Rational(891,598),
               s.Rational(385,299),s.Rational(119,598),s.Rational(9,598)])
epsilon=s.Rational(1,1000)
shifted=base+epsilon*direction
assert all(t>0 for t in shifted)
b=A+n; P=C+d/2; Q=h+d/2; Gamma=a*a*(b+P+Q)-b
equations=[a*(d+2*n)-2*A,m*P-D*Q,m*Gamma-b*Q,D*Gamma-b*P]
sub=dict(zip(coords,shifted))
assert all(s.factor(f.subs(sub))==0 for f in equations)
eta=s.symbols('eta')
assert all(s.expand(f.subs(dict(zip(coords,coords+eta*direction)),simultaneous=True)-f)==0
           for f in equations)

u,v=s.symbols('u v',real=True)
I=s.eye(2); Z=s.diag(1,-1); X=s.Matrix([[0,1],[1,0]])
A0=s.kronecker_product(Z,I); B0=s.kronecker_product(I,Z)
A1=s.kronecker_product(((1-u*u)*Z+2*u*X)/(1+u*u),I)
B1=s.kronecker_product(I,((1-v*v)*Z+2*v*X)/(1+v*v))
ops=[s.eye(4),A0,A1,B0,B1,A0*B0,A0*B1,A1*B0,A1*B1]
coef=[d+n+m+C+A+D+h,n-A-h,m-C-D,-n+C+A,-m+D-h,
      -d/2-n-A,d/2+h,-d/2-C,-d/2-m-D]
bell=sum((k*op for k,op in zip(coef,ops)),s.zeros(4))
contact={u:s.Rational(1,4),v:s.Rational(1,2)}
psi=s.Matrix([1,-s.Rational(1,8),s.Rational(1,8),1])
base_bell=bell.subs(dict(zip(coords,base)))
assert base_bell.subs(contact)*psi==s.zeros(4,1)
assert all(s.factor((psi.T*s.diff(base_bell,t).subs(contact)*psi)[0])==0 for t in (u,v))
value=s.factor((psi.T*bell.subs(sub).subs(contact)*psi)[0]/psi.dot(psi))
assert value==-s.Rational(963,1105000)
out={'a':2,'active_order':list(map(str,coords)),
     'base':list(map(str,base)),'direction':list(map(str,direction)),
     'epsilon':str(epsilon),'all_shifted_weights_positive':True,
     'angle_recovery_and_three_quadrics_preserved_identically':True,
     'base_kernel_and_two_angle_equations_exact':True,
     'negative_normalized_slack':str(value),
     'scope':'Quadratic identities and source signs do not characterize the fixed-a native hull.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))
