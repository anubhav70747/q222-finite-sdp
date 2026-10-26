"""Fresh boundary checks, using only an independently installed symbolic engine.

An exact-identity diagnostic, not a global representability proof.
"""
import sys
import sympy as s
import json
from pathlib import Path

out = {}
def check(name, value):
    z = s.cancel(value)
    out[name] = z == 0
    if z != 0:
        print('FAIL', name, s.factor(z), flush=True)

# The shared reciprocal parent.  The symbols k,q,t,l remove all radicals:
# x=k^2, z=q^2, sigma=t^2, lambda=l^2, rho=l^2 k^2/t^2.
u,w,W,Z,k,q,t,l = s.symbols('u w W Z k q t l', nonzero=True)
x,z,sigma,lam,rho = k*k,q*q,t*t,l*l,l*l*k*k/(t*t)
sr = l*k/t
Gamma = sigma*Z-u*w*W
R = w*Z/q-u*q*W
A = q*(u*u-x)*W
T = k*q*(u-k)*W
Pb = u*Z-x/sigma*w*W
a1 = l*t*u*Z
a2 = sr*u*w*W
a3 = l*x*q*W
a4 = l*k*w*Z/q
D12, D23 = a1-a2, a2-a3
packet = lam*sigma*Pb**2+rho*(u*u-x)*(w*w-sigma*z)*W**2
check('SF15 first normal', u*Gamma-k/sr*D12+(u-k)/sr*D23+t*T)
check('SF15 second normal', q*w*Gamma-sigma*z*R+(q*w+z*t)/sr*D23-(sigma*z/k+t*q*w/k)*T)
check('SF13 angular domination', A*A-T*T-q*q*(u-k)**2*W**2*(u*u+2*u*k))
check('US8 four-monomial circuit', packet+lam*x*R*R-sum(a*a for a in (a1,a2,a3,a4))+4*lam*x*u*w*W*Z)
check('US13 fourth-monomial normal', q*w*Gamma-(q*w+z*t)/sr*(a4-a2)+t*q*w*R)
check('US fourth-monomial angular conductor', T-(a4-a3)/l+k*R)
J = q*(u-k)*W
Fa = u*(Z-k/sigma*w*W)
check('US6 complete scalar', packet-lam*sigma*Fa*Fa-2*lam*k*u*J*J-lam*k*k*J*J-2*lam*k*u*J*R)

# GA domination and reserve allocation, with a free positive real reserve C.
C,U,V=s.symbols('C U V')
check('GA5 reserve algebra', C*(U+V)**2-2*C*U*U-2*C*V*V+C*(U-V)**2)
E=q*u*W-t*Z
check('GA3 first', Gamma+t*R+(w+t*q)*E/q)
check('GA3 second', A-(a1-a3)/l-u*E)
check('GA4 exact SOS domination', (q*u/w)**2*(Gamma+t*R)**2-(u*E)**2-u*u*E*E*(2*t*q/w+(t*q/w)**2))

# JR6 and the exact real Schur-complement calculation.
V=u*w*(u*u-x)*W
check('JR6', V*V-sigma*u*u*A*A-u*sigma*(u*u-x)*A*R+(u*u-x)*V*Gamma)
tt=s.symbols('tt')
ll=2*(tt*tt+tt+1)/(1+tt)**2
check('JR8 Schur complement', ll*(1+tt*tt)+2*tt*tt/(1+tt)-s.Rational(1,2)*(1+tt+ll*tt/2)**2-(4*tt+3)/(2*(1+tt)**4))

# A6: derive the signed core matrix directly from its four factored rows.
u,v=s.symbols('u v', nonzero=True)
W,X,Y,Z=s.symbols('W X Y Z')
Du,Dv=1+u*u,1+v*v
rm=s.Matrix([v,-1,u*v,-u]); rd=s.Matrix([u,u*v,-1,-v])
ra=s.Matrix([0,0,v,-1]); rb=s.Matrix([0,1,0,u])
ez=s.Matrix([0,0,0,1])
H=2*(rd*rd.T+Du*ra*ra.T+Dv*rb*rb.T-Du*Dv*ez*ez.T)
swap=s.Matrix([[0,0,0,1],[0,0,-1,0],[0,-1,0,0],[1,0,0,0]])
square=(H-Du*Dv*s.eye(4))**2-(Du*Dv)**2*s.eye(4)-4*u*v*Du*Dv*swap
check('CF core matrix identity', sum(a*a for a in square))

# Direct two-hard constrained determinants. Positive monomial coefficients
# determine their active Laurent supports. Report supports separately.
n,m,c,A,D,alpha,beta=s.symbols('n m c A D alpha beta', nonnegative=True)
rn=s.Matrix([0,1,0,0]); rc=s.Matrix([u,0,-1,0]); rA=s.Matrix([0,0,1,0])
M=H+4*n*Du*Dv*rn*rn.T+4*m*rm*rm.T+4*c*Dv*rc*rc.T+4*A*Du*Dv*rA*rA.T+4*D*rd*rd.T+4*alpha*Du*ra*ra.T+4*beta*Dv*rb*rb.T
bases={
 'lambda,c':s.Matrix([[-v,0],[1,0],[-u*v,0],[0,1]]),
 'lambda,m':s.Matrix([[-v,0],[1,0],[0,1],[-Dv/u,v]]),
 'lambda,alpha':s.Matrix([[-v,0],[1,0],[0,1],[0,v]])}
support={}
for key,B in bases.items():
    det=s.factor((B.T*M*B).det())
    # Divide common positive factors, leaving the Newton-square polynomial.
    # The factor is identified by exact factorization, not an inequality.
    numerator, denominator=s.fraction(det)
    factors=s.factor_list(numerator, u,v)[1]
    nontrivial=[f for f,p in factors if s.Poly(f,u,v).total_degree()>2]
    if len(nontrivial) != 1:
        raise ValueError((key, s.factor(det)))
    p=s.Poly(nontrivial[0],u,v)
    ss={}
    for (i,j),co in p.terms():
        if (i,j)==(2,2):
            continue
        mon=s.Poly(co,n,m,c,A,D,alpha,beta)
        positive=all(cv>0 for cv in mon.coeffs())
        out['A6 nonnegative outer coefficient '+key+' '+str((i,j))]=positive
        ss[str((i//2-1,j//2-1))]=str(s.factor(co))
    support[key]=ss
    print('Constrained determinant completed', key, flush=True)

dest=Path(__file__).with_suffix('.json')
dest.write_text(json.dumps({'engine':'SymPy '+s.__version__,'checks':out,'constrained_outer_coefficients':support,'scope':'Exact identities and coefficient signs only; coverage and estimates require the written proof.'},indent=2)+'\n')
print(f'{sum(out.values())}/{len(out)} independent boundary checks passed', flush=True)
if not all(out.values()):
    raise SystemExit(1)
