"""Exact Sylvester/minor certificates for the R2O20 full X,Z,V state matrix."""
import sympy as sp
u,j,p,b=sp.symbols('u j p b', nonzero=True)
Du=1+u*u; k=(j+1-p)/(j+1)
mc,nc=sp.symbols('mc nc')
X,Z,V=sp.symbols('X Z V')
Q=mc*((Du/u)*V-X-u*Z)**2+nc*Du*X**2+sp.Rational(1,2)*Du*V**2-Du*V*Z+sp.Rational(1,2)*(X+u*Z)**2-2*b*V*(u*X-Z)
M=sp.hessian(Q,(X,Z,V))/2
m0=j*p/(2*(1-p)*(j+1));n0=p/(2*(j-p)*(j+1));b0=k/2
R=sp.simplify(M.subs({mc:m0,nc:n0,b:b0}));exps=[m0+n0*Du+sp.Rational(1,2),p*u*u*Du*(j+1-p)/(4*(j-p)*(j+1)**2*(1-p)),p*u*u*Du*(j+1-p)**2/(8*(j-p)*(j+1)**3)]
minors=[R[0,0],R[:2,:2].det(),R.det()]
for n,(lhs,rhs) in enumerate(zip(minors,exps),1):
  ok=sp.factor(lhs-rhs)==0
  print('R2O21 minor',n,ok,flush=True)
  if not ok:print(sp.factor(lhs-rhs));raise SystemExit(1)
# R2O22 j=p<1, restrict X=0; b=1/[2(j+1)]
Q22=sp.factor(Q.subs({mc:(j*j/(2*(1-j)*(j+1))),X:0,b:1/(2*(j+1))}))
M22=sp.hessian(Q22,(Z,V))/2
print('R2O22 diagonal',sp.factor(M22[0,0]),flush=True)
print('R2O22 determinant',sp.factor(M22.det()),flush=True)
assert sp.factor(M22[0,0]-u*u*(j*j/(2*(1-j)*(j+1))+sp.Rational(1,2)))==0
assert sp.factor(M22.det()-u*u/(4*(j+1)**2))==0
# R2O23 p=1<j; on kernel of m square X=Du/u V-uZ, use X,V as independent => Z=(Du/u V-X)/u
n0_23=1/(2*(j*j-1))
Q23=sp.factor(Q.subs({mc:0,nc:n0_23,p:1,Z:(Du/u*V-X)/u}))
M23=sp.hessian(Q23,(X,V))/2
print('R2O23 diagonal',sp.factor(M23[0,0]),flush=True)
print('R2O23 determinant',sp.factor(M23.det()),flush=True)
assert sp.factor(M23[0,0]-Du/(2*(j*j-1)))==0
expected=-Du**2*(2*b*(j-1)-j)*(2*b*(j+1)-j)/(4*u*u*(j*j-1))
assert sp.factor(M23.det()-expected)==0
# R2O24 p=j=1, X=0,Z=Du V/u^2, mc=nc=0
Q24=sp.factor(Q.subs({mc:0,nc:0,p:1,j:1,X:0,Z:Du*V/u**2}))
print('R2O24',Q24,flush=True)
assert sp.factor(Q24-Du*(4*b-1)*V*V/(2*u*u))==0
print('all R2O finite cap identities PASS',flush=True)
