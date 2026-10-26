"""Reconstruct the C2 stationary Jacobian directly from the six event squares."""
from pathlib import Path
import json, time
import sympy as s
start=time.monotonic()
W,X,Y,Z,u,v,a,b,c,r,t=s.symbols('W X Y Z u v a b c r t')
du,dv=1+u*u,1+v*v
events=[du*dv*X**2,(v*W-X+u*v*Y-u*Z)**2,dv*(u*W-Y)**2,
        du*dv*Y**2,(u*W+u*v*X-Y-v*Z)**2,du*(v*Y-Z)**2]
subs={u:c,v:a*c,W:b/a,X:-c*r,Y:c*(b/a-a*t),Z:1}
M=s.Matrix([[s.diff(e,x).subs(subs) for e in events] for x in [W,X,Y,Z,u,v]])
z=c*c; D=1+z; E=1+a*a*z; h=b-a*a*t
g=D*b-(1-r+a*a*z*t); j=1+z*r-t
Delta=1+(1+z)*b+z*(1+z)*b*b+2*z*(1+z)*b*r+z*z*r*r+2*z*r
Delta+=a*a*z*((1+z)*b+1+z*r*(2-r))
factors=[2*c*r,2*c*g,2*a*c*t,-2*c*h/a,2*a*c*j,2*(1-z*h)]
B=M.copy()
for k, f in enumerate(factors):
    B[:,k]=B[:,k].applyfunc(lambda x:s.cancel(x/f))
B[0,:]=B[0,:]/c
B_expected=s.Matrix([[0,a,E,0,-1,0],[-D*E,-1,0,0,-a*z,0],
    [0,a*z,-E,-D*E,1,-a*c*D],[0,-c,0,0,a*c,D],
    [z*r*E,z*h-1,b*E/a,-z*h*E/a,a*z*r-b/a,c*(1-z*h)],
    [a*z*r*D,b/a+z*h/a,a*a*z*t,-z*D*h,1+z*r,-c*D*h/a]])
assert all(s.cancel(x)==0 for x in B-B_expected)
det=s.factor(B.det(method='domain-ge'))
assert s.cancel(det-D**2*E**2*Delta)==0
expected=-64*a*c**6*r*g*t*h*j*(1-z*h)*D**2*E**2*Delta
assert s.cancel(det*c*s.prod(factors)-expected)==0
report={'event_gradient_matches_all_36_printed_entries':True,
        'reduced_determinant_matches_Delta':True,'EN16_full_stationary_determinant':True,
        'seconds':time.monotonic()-start,
        'scope':'Exact identities. Positivity follows from source domain, checked separately.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
