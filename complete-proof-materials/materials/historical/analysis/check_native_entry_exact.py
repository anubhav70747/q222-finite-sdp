"""Exact source-field contact identities used in EN, independent of tensor SOS."""
from pathlib import Path
import json, time
import sympy as s
started=time.monotonic()
u,v,S,Y0,theta=s.symbols('u v S Y0 theta')
d,m,C,A,D,al,be=s.symbols('d m C A D al be')
du,dv=1+u*u,1+v*v
W,X,Y,Z=s.symbols('W X Y Z')
psi=s.Matrix([W,X,Y,Z])
F=4*m*(v*W-X+u*v*Y-u*Z)**2+4*C*dv*(u*W-Y)**2+4*A*du*dv*Y**2
F+=(4*D+2*d)*(u*W+u*v*X-Y-v*Z)**2+(4*al+2*d)*du*(v*Y-Z)**2
F+=(4*be+2*d)*dv*(X+u*Z)**2-2*d*du*dv*Z**2
contact={W:S,X:u*(theta-1),Y:u*Y0,Z:1}
eqs=[s.diff(F,x)/2 for x in psi]+[s.diff(F,u),s.diff(F,v)]
eqs=[s.expand(e.subs(contact).subs(C,2)) for e in eqs]
unknowns=[d,m,A,D,al,be]
M,rhs=s.linear_eq_to_matrix(eqs,unknowns)
det=s.factor(M.det(method='domain-ge'))
Rm=v*S-u*(theta-1)+u*u*v*Y0-u
RD=u*S+u*u*v*(theta-1)-u*Y0-v
calD=Y0*du*(u*Y0+v)*(1-u*v*Y0)*(S+u*v*(theta-1))*Rm*(-RD)
expected=-8192*theta*u*u*v*du**2*dv**2*calD
ratio=s.cancel(det/expected)
assert ratio==1,('EN8 normalization',ratio)
M_m=M.copy(); M_m[:,1]=rhs
det_m=s.factor(M_m.det(method='domain-ge'))
res=s.factor(det_m*(S+u*v*(theta-1))*Rm-det*2*u**3*(S-Y0)*(1-theta)*dv)
assert res==0,('EN7 Cramer identity',res)
seed={u:1,v:1,S:1,Y0:s.Rational(1,4),theta:s.Rational(1,2),d:s.Rational(27,2),m:4,A:9,D:s.Rational(5,4),al:s.Rational(13,4),be:s.Rational(9,4)}
assert all(e.subs(seed)==0 for e in eqs)
assert det.subs(seed)!=0
report={'EN8_full_stationary_minor':True,'EN7_cross_multiplied_identity':True,
        'EN9_stationary_seed_and_nonzero_minor':True,'seconds':time.monotonic()-started,
        'scope':'Exact polynomial identities; connected-graph and boundary arguments are checked analytically.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
