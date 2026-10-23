"""Exact general GE1 contact identities; cancel source denominators termwise."""
from pathlib import Path
import json
import sympy as s

z,y,H,L,r,x,w=s.symbols('z y H L r x w')
n,m,C,A,D,al=s.symbols('n m C A D al')
half=s.Rational(1,2)
B=H+L; a=1+z*r; T=B+z*r; U=B-z*y*r; SS=a+L
yj=y*a-L; G=B+z*(H+r-1)
de=z+(1+z)*B+(1+z)*B**2+2*z*(1+z)*B*r+z**3*r**2+2*z**2*r+y*((1+z)*B+z+z*z*r*(2-r))
M=(1+z)*T*H+y*G+z*yj+z*r*L*(1+y)
# These entries are exactly 2 Delta times the corresponding source weight.
scaled={n:(T*a*SS-r*de)/r,m:z*z*(1+y)*(1-r)*a*SS/G,
        C:a*M/L,A:y*a*T*SS/H,D:(1+z)*T*U*SS/yj,al:T*M/(1-H)}
K=B+x*(H+r-1); J=L-w*(1+r*x)
normals=[m*w*K+C*(1+w)*x*L+(D+half)*x*J,
 m*w*(1+x)*K+A*(1+x)*(1+w)*x*H+(al+half)*(1+x)*x*w*(H-1),
 -n*r*(1+x)*(1+w)*x*x*w-m*w*x*K+(D+half)*x*x*w*J+half*(1+w)*x*x*w*(1-r)]
p=(n*(1+x)*(1+w)*r*r*x*x*w+m*w*K*K+C*(1+w)*x*L*L
   +A*(1+x)*(1+w)*H*H*x+(D+half)*x*J*J
   +(al+half)*(1+x)*x*w*(1-H)**2+half*(1+w)*x*x*w*(1-r)**2
   -half*(1+x)*(1+w)*x*w)
exprs=dict(zip(['normal1','normal2','normal3'],normals))
exprs.update(value=p,partial_x=s.diff(p,x),partial_w=s.diff(p,w))
results={}
zero={v:0 for v in scaled}
for name,expr in exprs.items():
    at=expr.subs({x:z,w:y})
    pieces=[2*de*at.subs(zero)]
    for variable,weight in scaled.items():
        coefficient=s.factor(s.diff(at,variable))
        piece=s.cancel(coefficient*weight)
        assert s.denom(piece)==1,(name,variable,s.denom(piece))
        pieces.append(piece)
    numerator=s.Poly(s.expand(sum(pieces)),z,y,H,L,r)
    assert numerator.is_zero,name
    results[name]='zero polynomial after multiplication by 2 Delta'
    print(name+': exact zero',flush=True)
out={'identities':results,'assumptions':'Positive GE1 denominators; 0<H<1, z,y,L>0, 0<r<1, j>0. Global validity additionally uses nonnegative source weights and SG7.',
     'scope':'All general GE1 contact parameters, not only two special families. No finite testing span or lift completeness claim.'}
# Verify the full state change independently of the contact equations.
u,v,N1,N2,N3,T0=s.symbols('u v N1 N2 N3 T0')
W=N1+N2+B*T0; X=u*u*v*(N3-r*T0); Y=u*(N2+H*T0); Z=u*v*T0
Du=1+u*u; Dv=1+v*v
parent=(n*Du*Dv*X*X+m*(v*W-X+u*v*Y-u*Z)**2
        +C*Dv*(u*W-Y)**2+A*Du*Dv*Y*Y
        +(D+half)*(u*W+u*v*X-Y-v*Z)**2
        +(al+half)*Du*(v*Y-Z)**2+half*Dv*(X+u*Z)**2-half*Du*Dv*Z*Z)
contact_form=(n*(1+x)*(1+w)*x*x*w*(N3-r*T0)**2
 +m*w*(N1+(1+x)*N2-x*N3+K*T0)**2+C*(1+w)*x*(N1+L*T0)**2
 +A*(1+x)*(1+w)*x*(N2+H*T0)**2
 +(D+half)*x*(N1+x*w*N3+J*T0)**2
 +(al+half)*(1+x)*x*w*(N2+(H-1)*T0)**2
 +half*(1+w)*x*x*w*(N3+(1-r)*T0)**2-half*(1+x)*(1+w)*x*w*T0*T0)
assert s.Poly(s.expand(parent-contact_form.subs({x:u*u,w:v*v})),u,v,N1,N2,N3,T0,n,m,C,A,D,al,H,L,r).is_zero
state_map=s.Matrix([W,X,Y,Z]).jacobian([N1,N2,N3,T0])
assert s.factor(state_map.det())==-u**4*v*v
out['full_state_identity']='Exact polynomial identity; state-map determinant -u^4 v^2.'
print('Full four-state congruence: exact polynomial identity',flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
