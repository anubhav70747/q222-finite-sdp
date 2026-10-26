"""Actual imported UI Laurent model: full determinant, D-kernel and contact jets."""
import json
from pathlib import Path
import sympy as s
u,v,W,Y,Z=s.symbols('u v W Y Z',nonzero=True)
h,C,A,d,D,al=s.symbols('h C A d D al')
U,T=s.symbols('U T',nonzero=True)
du=1+u*u
f=4*h*v*v*(W+u*Y)**2+4*C*(u*W-Y)**2+4*A*du*Y*Y
f+=(4*D+2*d)*(u*W-Y-v*du*Z)**2+4*al*du*Z*Z-4*d*du*v*Y*Z-2*d*du*v*v*Z*Z
L0=A*C*D-A*D*d/2+h*al*(C+D+d/2)-d*d*(A+C+D)/4-d**3/8
L=A*C*(D+d/2)*U+h*al*(A+C+D+d/2)/U+A*al*(C+D+d/2)/T
L+=h*D*(A+C+d/2)*T/U+h*C*(D+d/2)*U*T
L+=h*(A*D+2*C*D+A*d/2+C*d/2+D*d/2+d*d/4)*T+L0
M=s.hessian(f,(W,Y,Z))/2
out={'CL3_full_determinant':s.expand(M.det()-64*du**2*u*u*v*v*L.subs({U:u*u,T:v*v}))==0}
restricted=s.expand(f.subs(W,(Y+v*du*Z)/u))
ML=s.hessian(restricted,(Y,Z))/2
LD=h*T*((C+A+d/2)/U+(2*C+A+d/2)+C*U)+h*al*(1+1/U)+A*C*(U+1)+A*al/T-A*d/2-d*d/4
out['CL6_complete_D_kernel_determinant']=s.factor(ML.det()-16*du**2*v*v*LD.subs({U:u*u,T:v*v}))==0
x,eta,a=s.symbols('x eta a',nonzero=True)
r=1-a*eta; K=1+x; R=1+x*r; E=1-x*eta; J=K-a
P4=K*E+x*(1-r); Pd=1+2*x+r*x*x
src={h:2*x*eta/(R*E),C:2/a,A:2*eta*eta*P4/(K*r*(r+eta)*E),d:4*eta*Pd/(K*(r+eta)*E),D:2/(J*E)-2*eta*Pd/(K*(r+eta)*E),al:2*x*x*eta*R/(K*E)}
for name,expr in [('value',L),('U_jet',s.diff(L,U)),('T_jet',s.diff(L,T))]:
    out['CL5_'+name]=s.factor(expr.subs({U:x,T:x*eta}).subs(src,simultaneous=True))==0
out['all_zero']=all(out.values())
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
