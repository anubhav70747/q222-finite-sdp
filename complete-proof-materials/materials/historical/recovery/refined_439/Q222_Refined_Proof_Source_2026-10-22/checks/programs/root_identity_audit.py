"""Independent exact diagnostics for the C1 proof; no archive code is imported.

These calculations corroborate identities. They do not prove chart coverage,
ordered-field estimates, or semidefinite representability.
"""
import sys
import json
from pathlib import Path
import sympy as s

results = {}
def check(name, expression):
    value = s.cancel(expression)
    results[name] = value == 0
    if value != 0:
        print('FAILED', name, s.factor(value), flush=True)

# The cubic six-root identity and all six conductor identities (EZ15--16).
r,z,h,c = s.symbols('r z h c')
u,v,P,G,kappa = r**3,z**3,h**6*c**6,h**6,c**6
U,V = u*u,v*v
ls = [r**4*z**2,c**2*r**2*z**4,h**6*c**4]
fs = [r**2*z,c*r*z**2,h**3*c**2]
roots = [f*(ls[0]-ls[j]) for f in fs for j in (1,2)]
R1,R2,R3,R4,R5,R6 = roots
D,E = U-P,V-G
J1 = r**3*R1+c**2*r*z**2*R2
J2 = -h**6*R1+z**6*R2-h**3*z**3*R5
check('EZ15', U**2*V+kappa*U*V**2+P**2*G-3*P*U*V-
      sum(roots[i]**2-roots[i]*roots[i+1]+roots[i+1]**2 for i in (0,2,4)))
check('EZ16a', u*v*D-J1)
check('EZ16b', u*u*v*E-J2)
check('EZ16c', kappa*u*v*E-c**2*r*z**2*(R2-R1)+c**3*z**3*R3)
check('EZ16d', U*V*D**2-J1**2)
check('EZ16e', U*V*D*E-z**6*R2**2+h**6*R1**2+h**3*z**3*R1*(R5+R6))
check('EZ16f', kappa*U*V*E**2-z**6*(2*R1**2-3*R1*R2+R2**2+R3**2)-h**3*z**3*R1*(R5-R6))

# Exact native source identities and division-free elimination (NC2, NC6--7).
x,r,h,z,eta = s.symbols('x r h z eta')
theta=1-h; K=1+x; R=1+x*r; T=1+x*h; G=h+z
p0=1+2*x*h+r*x*x*G+x*z
p1=x*(K*z+x*h*theta)
Pd=p0+eta*p1
B5=K*(1+x*eta*z)*(r+eta)-(eta*T-1+r)*Pd
Bbeta=K*x*(r+eta)*G-theta*Pd
AA=theta*T+K*z
CC=R*G*(1+r*x*h)
BB=h*h*r*x*x+h*r*x*(2-h)+h*h*x+h*r*x*x*z+h*x*z+h*x+h+r*x*z+x*z+2*z+1
E=1-x*eta*h
check('NC6', h*B5-(-AA*E*E+BB*E-CC))
a2=K*x*z-T*p1
a1=K*(1+x*z*r)-T*p0+(1-r)*p1
a0=K*r+(1-r)*p0
b1=K*x*G-theta*p1
b0=K*x*r*G-theta*p0
check('NC7 first', b1*b1*B5-a2*Bbeta**2-(b1*a1-2*a2*b0)*Bbeta-
      (b1*b1*a0-b1*a1*b0+a2*b0*b0))
check('NC7 positive factor', b1*b1*a0-b1*a1*b0+a2*b0*b0-
      x*x*h*h*K*K*R*G*G*(1+2*x*h+r*x*x*h*h+(1-r)*x*z))
check('NC7 b1', b1-x*h*(K*(1+z)-x*theta**2))

# The common polynomial contact conductor (CC1--6).
d,k,b,u,w,v = s.symbols('d k b u w v')
A=w*v; B=d*v; C=d*b*b*u*w; D0=d*d*k**3*b*b
X=B-C; Y=A-k*C; Z=u*B-k*C; H=k*C-D0
star=Y+Z-H-2*k*X
check('CC-u', C*(u-k)-Z+u*X)
check('CC-w', C*(w-d*k)-d*Y+w*X)
check('CC-v', C*(v-d*k*k*b*b)-b*b*(w*Z+d*k*Y-k*w*X))
check('CC-star', star-(v-d*k*k*b*b)*((w-d*k)+d*(u-k))+d*k*b*b*(u-k)*(w-d*k))
# q^2 = (a1 a2 - a3 a4)^2. The product of all a_i is k^2 BC.
q2=A*(u*B)+(k*C)*D0-2*k*k*B*C
check('CC-repair', k*C*star-q2+Y*Z)
jac=s.Matrix([X,Y,Z]).jacobian([u,w,v]).subs({u:k,w:d*k,v:d*k*k*b*b})
check('CC gradient determinant', jac.det()+d**4*k**4*b**4)

# Two-pairing strict circuit identity.
a,b,c,d=s.symbols('a b c d')
H1=a*a-c*c; H2=b*b-c*c; H3=c*c-d*d
q1=a*b-c*d; q2=a*c-b*d
check('two-pairing circuit', a**4+b**4+c**4+d**4-4*a*b*c*d-
      (H1**2+H2**2+H3**2-H1*H2+H2*H3+q1**2+q2**2))

# Small-e conductor, including its extra generator (SE5--7).
eps,j,E,u,w,v=s.symbols('eps j E u w v')
k=eps*j; A=eps*w*v; B=v; C=E*u*w; D0=eps**2*j**3*E
X=B-C; Y=A-k*C; Z=u*B-k*C
check('SE5u', C*(u-k)-Z+u*X)
check('SE5w', eps*C*(w-j)-Y+eps*w*X)
check('SE5v', C*(v-eps*j*j*E)-E*(w*Z+j*Y-k*w*X))
check('SE6', C*v*(w-j)-E*w*(w-j)*Z-E*j*w*Y+E*k*w*w*X)
Tu,Tv,U,V=s.symbols('Tu Tv U V')
S=(u-k)*Tu+(w-j)*(eps*U+v*V)+(v-eps*j*j*E)*Tv
rhs=X*(-u*Tu-eps*w*U-E*k*w*w*V-E*k*w*Tv)+Y*(U+E*j*w*V+E*j*Tv)+Z*(Tu+E*w*(w-j)*V+E*w*Tv)
check('SE7', C*S-rhs)

# Full signed core and critical-parent identity (SEO11), not only scalar jets.
u,v,W,X,Y,Z,c,A,D,d=s.symbols('u v W X Y Z c A D d')
Du=1+u*u; Dv=1+v*v
core=2*(u*W+u*v*X-Y-v*Z)**2+2*Du*(v*Y-Z)**2+2*Dv*(X+u*Z)**2-2*Du*Dv*Z**2
parent=d*core+4*c*Dv*(u*W-Y)**2+4*A*Du*Dv*Y**2+4*D*(u*W+u*v*X-Y-v*Z)**2
B=X+u*Z; J=u*W-Y; V=v*Z; CC=4*D+2*d
rhs=4*c*Dv*J**2+4*A*Du*Dv*Y**2+CC*(J+u*v*B-Du*V)**2+2*d*Dv*B**2+2*d*Du*(v*v*Y**2-2*Y*V-V*V)
check('SEO11 whole parent', parent-rhs)

# Critical zero-section scalars in both hard-normal choices.
cstar=d*(d+2*A)*(d+2*D)/(2*(4*A*D-d*d))
b=(4*A*D-d*d)/(A*CC); y=d/(2*A)
B,J,V=s.symbols('B J V')
Fc=4*cstar*Dv*J**2+4*A*Du*Dv*Y**2+CC*(J+u*v*B-Du*V)**2+2*d*Dv*B**2+2*d*Du*(v*v*Y**2-2*Y*V-V*V)
g=d*(2*A+d)/A
gamma=D*d*(2*A+d)**2/(A*A*(2*D+d))
zero1=Fc.subs({B:0,Y:y*V,J:(b+u*u)*V})/V**2
check('SEO13', zero1-(g*u*u+gamma*v*v+4*cstar*u**4+g*(y+2)*u*u*v*v+4*cstar*u**4*v*v))
zero2=Fc.subs({B:0,Y:y*V,J:b*V})/V**2
check('SEO15', zero2-(g*u*u+gamma*v*v+CC*u**4+g*y*u*u*v*v))

# Native-independent toric identities: use positive monomial square roots.
# u=U^4, v=V^4, w=W^4, b=B^4, g=G^4, c=C^4 removes radicals.
U,V,W,B,G,C=s.symbols('U V W B G C')
u,v,w,b,g,c=U**4,V**4,W**4,B**4,G**4,C**4
A=c*u*w/b**2; BB=c*v*w/(b*g); CC=u*v/b; DD=c*c*g
Gamma=C**4*U**2*V**2*W**2/B**4
sqAC=C**2*U**4*V**2*W**2/B**6
sqBD=C**6*V**2*W**2/B**2
sqBC=C**2*U**2*V**4*W**2/(B**4*G**2)
sqAD=C**6*G**2*U**2*W**2/B**4
sqAB=C**4*U**2*V**2*W**4/(B**6*G**2)
sqCD=C**4*G**2*U**2*V**2/B**2
check('SC18u', Gamma*(u-b*c)-U**2*B**2*C**2*(sqAC-sqBD))
check('SC18v', Gamma*(v-c*g)-V**2*C**2*G**2*(sqBC-sqAD))
check('SC18w', Gamma*(w-b*g)-W**2*B**2*G**2*(sqAB-sqCD))
# SC20 follows by factoring a difference of squares in SC18; verify directly.
sqrtGamma=C**2*U*V*W/B**2
check('SC20u', Gamma*(U**2-B**2*C**2)-U*B*C*sqrtGamma*
      (C*U**2*V*W/B**3-C**3*V*W/B))
check('SC20v', Gamma*(V**2-C**2*G**2)-V*C*G*sqrtGamma*
      (C*U*V**2*W/(B**2*G)-C**3*G*U*W/B**2))

out=Path(__file__).with_suffix('.json')
out.write_text(json.dumps({'engine':'SymPy '+s.__version__, 'checks':results,
                          'scope':'Exact identities only; not a proof of the global theorem.'},indent=2)+'\n')
print(f'{sum(results.values())}/{len(results)} independent exact identities passed',flush=True)
if not all(results.values()):
    raise SystemExit(1)
