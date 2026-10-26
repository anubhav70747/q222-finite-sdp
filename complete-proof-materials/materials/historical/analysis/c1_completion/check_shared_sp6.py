import json
from pathlib import Path
import sympy as s

u,w,W,X,Y,Z=s.symbols('u w W X Y Z')
x,sg,z,b,h=s.symbols('x sg z b h',nonzero=True)
q=s.symbols('q',nonzero=True)
m,C,A,D,al,be,kap=s.symbols('m C A D al be kap')
r=b*z; Q=1-b; theta=1-h; e=x*h/sg
p=A*r; ell=al*z; l=ell*Q
du=1+u*u; dw=1+w*w; t=u*u; xi=w*w; g=h-Q
MW=W-w*X+u*Y-u*w*Z
RW=u*W-Y; JW=u*w*W+u*X-w*Y-Z; BW=X+u*Z
FW=4*m*MW**2+4*C*dw*RW**2+4*A*du*dw*Y**2
FW+=4*D*JW**2+4*al*du*(Y-w*Z)**2+4*be*dw*BW**2
FW+=2*(JW**2+du*(Y-w*Z)**2+dw*BW**2-du*dw*Z**2)
NX=X+u*h*Z; NY=(Y-b*w*Z)/q
GAM=sg*Z-u*w*W; RR=w*Z/q-u*q*W; AA=q*(t-x)*W
PB=u*Z-(x/sg)*w*W; XE=X+e*w*W
nu=(2*D-2*be*theta+h)/sg
lam=(8*D*h+4*be*theta**2+2*h*h)/sg-4*b*p*sg*z
chi=8*(p*b+l*Q)-4*h+4*e*h*(2*D+1)-8*D*h
chi+=4*z*(4*D*b*h+4*p*b*sg+(1+2*m)*g*g+2*be*theta**2)
BX=8*kap*e
AY=8*(p-l)+4*z*((1+2*m)*g+2*D*h+2*p*sg)
RX=-4*(2*D*b+(1+2*m)*g-2*be*theta)
RG=-4*(2*(C+D)*(r-1)+2*p+r-1+2*p*t)
RY=4*(t*(2*D*h+2*p*sg+(1+2*m)*g)+b*(2*C+2*D+1)*(xi+sg*z)
         +(2*C+1)*(b-sg)+2*D*(1-sg)+2*p*sg)
k0=2*b*((2*C+1)*b+4*D+4*p*sg+2*(2*C+2*D+1)*sg*(r-1))
k1=2*(1+2*m)*g*g+8*D*h*b+8*p*b*sg+4*be*theta**2
k2=2*b*b*(2*C+2*D+1); KR=k0+k1*t+k2*xi
HY=4*du*(A*z*dw+ell)+4*m*t*z+4*C*z*dw+(4*D+2)*xi*z+2*du*z
CG=4*du*(p*b*dw+l*Q); CN=8*du*(p*dw-l)
GG=2*(2*D+1)*(u*XE)**2+(2*(1+xi)+4*m*xi)*NX**2+HY*NY**2
GG+=(4*D/sg**2+4*b*p*z*t)*GAM**2+(CG+z*KR)*RR**2-4*u*nu*NX*GAM
GG+=-4*u*q*w*(2*m+2*D+1)*NX*NY+q*w*RG*(NY+b*RR)*GAM
GG+=u*q*w*RX*NX*RR+(CN+z*RY)*NY*RR+q*w*BX*NX*AA+u*AY*NY*AA
GG+=4*kap*h/sg*AA**2+lam*u*(GAM*PB+AA*RR)+chi*u*AA*RR
GG+=4*be*dw*NX**2
diff=s.Poly(s.expand(FW-GG),W,X,Y,Z)
K=1+x; J=1+x*h-sg*(1-r); M=1+x*r-x*z*theta
P=sg+K+sg*z*(1+b*x+2*x*h)+sg*b*x*x*h*z*z+x*x*h*theta*z
E=1-e
pp=(z*(1+x*h)*E+(1-r)*(1+x*h*z))/(2*P)
ll=(1+sg*r)*(e*(1+x*r)*(1+sg*z)+M*E)/(2*P)
kk=K*(1+sg*r)*(1+sg*z)/(2*M*P)
source={A:pp/r,al:ll/(z*Q),
 C:K*z*(1+sg*r)*E/(2*(1-r)*P),
 D:(K*(1+sg*z)*(1+sg*r)-J*P)/(2*J*P),
 be:(K*(1+sg*r)*(1+x*h*z)-theta*P)/(2*theta*P),
 kap:kk,m:kk*x*x*z*h/sg}
out={}
for mon,coeff in diff.terms():
    # q occurs only with even powers after expansion.
    coeff=s.cancel(coeff).subs(q*q,z)
    value=s.factor(coeff.subs(source,simultaneous=True))
    out[str(mon)]={'zero':value==0}
    if value!=0: out[str(mon)]['residual']=str(value)
out['all_zero']=all(v['zero'] for v in out.values())
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
