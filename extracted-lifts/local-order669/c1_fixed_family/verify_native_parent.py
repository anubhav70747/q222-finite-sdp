"""Check P4=SP6 for all native stationary-source parameters, exactly."""
import json,time
from pathlib import Path
import sympy as s
here=Path(__file__).resolve().parent
x,sig,z,b,h,u,w=s.symbols('x sigma z b h u w')
q=s.symbols('q',positive=True)
W,X,Y,Z=s.symbols('W X Y Z')
m,C,A,D,al,be,p,l,kap=s.symbols('m C A D alpha beta p l kappa')
Q=1-b;th=1-h;r=b*z;H=h/sig;e=x*h/sig
Du=1+u*u;Dw=1+w*w
Nx=X+u*h*Z;Ny=(Y-b*w*Z)/q;R=w*Z/q-u*q*W
Gam=sig*Z-u*w*W;At=q*(u*u-x)*W
Pbeta=u*Z-x/sig*w*W;Xe=X+e*w*W
JC=Ny+b*R;JL=Ny-Q*R
nu=(2*D-2*be*th+h)/sig
lam=(8*D*h+4*be*th**2+2*h*h)/sig-4*b*p*sig*z
KA=4*kap*H-lam
Bx=8*kap*e
g=h-Q
Ay=8*(p-l)+4*z*((1+2*m)*g+2*D*h+2*p*sig)
Rx=-4*(2*D*b+(1+2*m)*g-2*be*th)
Rgam=-4*(2*(C+D)*(r-1)+2*p+r-1+2*p*u*u)
chi=8*(p*b+l*Q)-4*h+4*e*h*(2*D+1)-8*D*h+4*z*(4*D*b*h+4*p*b*sig+(1+2*m)*g*g+2*be*th*th)
k1=2*(1+2*m)*g*g+8*D*h*b+8*p*b*sig+4*be*th*th
Eyy=2*z*((1+2*m)*u*u+(1+2*D)*w*w)
ERR=z*(8*b*D+8*b*p*sig+4*b*sig*(2*C+2*D+1)*(r-1)+k1*u*u+2*b*b*(2*D+1)*w*w)
EyR=4*z*(u*u*(2*D*h+2*p*sig+(1+2*m)*g)+b*(2*D+1)*w*w+sig*(2*C+1)*(r-1)+2*D*(1+sig*(r-1))+2*p*sig)
Eplus=4*be*Dw*Nx*Nx+2*(2*D+1)*(u*Xe)**2+4*b*p*z*u*u*Gam*Gam
Eplus+=4*(p/b)*(Du*Dw-1)*JC*JC+4*(l/Q)*u*u*JL*JL+4*C*z*w*w*JC*JC
SP=Eplus+(2*(1+w*w)+4*m*w*w)*Nx*Nx+(4*p/b+2*z*(2*C+1))*JC*JC+4*l/Q*JL*JL
SP+=Eyy*Ny*Ny+ERR*R*R+EyR*Ny*R+4*D/sig**2*Gam*Gam
SP+=-4*nu*Nx*u*Gam+Rgam*JC*q*w*Gam
SP+=-4*u*q*w*(2*m+2*D+1)*Nx*Ny+u*q*w*Rx*Nx*R+q*w*Bx*Nx*At+u*Ay*Ny*At
SP+=KA*At*At+chi*u*At*R+lam*sig*Pbeta*Pbeta+lam*x/sig*(u*u-x)*(w*w-sig*z)*W*W
P4=4*m*(W-w*X+u*Y-u*w*Z)**2+4*C*Dw*(u*W-Y)**2+4*A*Du*Dw*Y*Y
P4+=(4*D+2)*(u*w*W+u*X-w*Y-Z)**2+(4*al+2)*Du*(Y-w*Z)**2
P4+=(4*be+2)*Dw*(X+u*Z)**2-2*Du*Dw*Z*Z
# q^2=z cancels algebraically before source substitution.
raw=s.Poly(s.expand((SP-P4)*q*q),W,X,Y,Z,u,w)
rawterms=[]
for mon,cf in raw.terms():
 cf=s.cancel(cf.subs(q*q,z))
 if cf!=0: rawterms.append((mon,cf))
print('Pre-source nonzero polynomial coefficients:',len(rawterms),flush=True)
K=1+x;Rs=1+x*r;M=Rs-x*z*th;J=1+x*h-sig*(1-r);E=1-e
P=sig+K+sig*z*(1+b*x+2*x*h)+sig*b*x*x*h*z*z+x*x*h*th*z
pn=(z*(1+x*h)*E+(1-r)*(1+x*h*z))/(2*P)
ln=(1+sig*r)*(e*Rs*(1+sig*z)+M*E)/(2*P)
kn=K*(1+sig*r)*(1+sig*z)/(2*M*P)
native={p:pn,l:ln,kap:kn,
        m:kn*x*x*z*H,C:K*z*(1+sig*r)*E/(2*(1-r)*P),A:pn/r,
        al:ln/(z*Q),D:(K*(1+sig*z)*(1+sig*r)-J*P)/(2*J*P),
        be:(K*(1+sig*r)*(1+x*h*z)-th*P)/(2*th*P)}
results=[]
for i,(mon,cf) in enumerate(rawterms):
 val=s.cancel(cf.subs(native))
 ok=val==0
 results.append({'monomial':list(mon),'passed':ok})
 print(i+1,'/',len(rawterms),ok,flush=True)
 if not ok:
  (here/'native_failure.txt').write_text(s.sstr(s.factor(val)))
  raise AssertionError(mon)
report={'identity':'P4 = SP6','method':'All coefficients in the independent polynomial variables W,X,Y,Z,u,w vanish after exact rational P1 source substitution and q^2=z.',
 'pre_source_nonzero_coefficients':len(rawterms),'all_passed':True,'checks':results}
(here/'NATIVE_PARENT_VERIFIED.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS full rational native parent',flush=True)
