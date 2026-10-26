"""Independent native-entry cofactor and contact-minor checks."""
import sys
import sympy as s
from pathlib import Path
import json
u,v,W,X,Y,Z=s.symbols('u v W X Y Z')
d,m,C,A,D,alpha,beta=s.symbols('d m C A D alpha beta')
Du,Dv=1+u*u,1+v*v
rm=s.Matrix([v,-1,u*v,-u]);rc=s.Matrix([u,0,-1,0])
rA=s.Matrix([0,0,1,0]);rd=s.Matrix([u,u*v,-1,-v])
ra=s.Matrix([0,0,v,-1]);rb=s.Matrix([0,1,0,u]);ez=s.Matrix([0,0,0,1])
H=2*(rd*rd.T+Du*ra*ra.T+Dv*rb*rb.T-Du*Dv*ez*ez.T)
P=d*H+4*m*rm*rm.T+4*C*Dv*rc*rc.T+4*A*Du*Dv*rA*rA.T+4*D*rd*rd.T+4*alpha*Du*ra*ra.T+4*beta*Dv*rb*rb.T
results={}
def check(name,expression):
    f=s.factor(expression)
    results[name]=f==0
    if f!=0:print('FAIL',name,f,flush=True)
cof44=P.minor_submatrix(3,3).det(method='domain-ge')
cof24=P.minor_submatrix(3,1).det(method='domain-ge')
f=(d+2*D)*((d+2*alpha)*(m*Dv+C*u*u*v*v)+2*A*Dv*(m+C*u*u))+2*m*C*(d*(u*u*v*v+Dv)+2*alpha+2*D*v*v*Du)
check('EN6 exact omitted-row cofactor',u*cof44+cof24-16*u*v*v*Du*Du*Dv*f)
print('Cofactor completed',flush=True)
S,Y0,theta=s.symbols('S Y0 theta')
psi=s.Matrix([W,X,Y,Z]);F=(psi.T*P*psi)[0]
contact={W:S,X:u*(theta-1),Y:u*Y0,Z:1}
eqs=list((P*psi).subs(contact))+[s.diff(F,u).subs(contact),s.diff(F,v).subs(contact)]
cols=[d,m,A,D,alpha,beta]
M=s.Matrix([[s.diff(eq,a) for a in cols] for eq in eqs])
det=s.factor(M.det(method='domain-ge'))
Rm=v*(S+u*u*Y0)-u*theta
negRd=v-u*(S-Y0)-u*u*v*(theta-1)
calD=Y0*Du*(u*Y0+v)*(1-u*v*Y0)*(S+u*v*(theta-1))*Rm*negRd
expected=-8192*theta*u*u*v*Du*Du*Dv*Dv*calD
ratio=s.cancel(det/expected)
print('Native selected-minor ratio to EN8:',ratio,flush=True)
check('EN8 selected contact minor',det-expected)
# EN7 after solving by Cramer's rule at C=2: its cross-multiplied form
# is checked without dividing by the minor or by any contact amplitude.
rhs=-s.Matrix([eq.subs({**{a:0 for a in cols},C:2}) for eq in eqs])
Mm=M.copy();Mm[:,1]=rhs
detm=s.factor(Mm.det(method='domain-ge'))
check('EN7 Cramer numerator identity',detm*(S+u*v*(theta-1))*Rm-det*2*u**3*(S-Y0)*(1-theta)*Dv)
Path(__file__).with_suffix('.json').write_text(json.dumps({'engine':'SymPy '+s.__version__,'checks':results,'minor_ratio':str(ratio)},indent=2)+'\n')
print(f'{sum(results.values())}/{len(results)} native-entry identities passed',flush=True)
if not all(results.values()):raise SystemExit(1)
