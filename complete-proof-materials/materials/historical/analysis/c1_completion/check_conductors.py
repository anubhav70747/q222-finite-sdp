"""Independent small, exact replays of operative conductor identities; no solver."""
import json
from pathlib import Path
import sympy as s
out={}
def ck(name,expr):
    val=s.factor(expr)
    out[name]={'zero':val==0}
    if val!=0: out[name]['residual']=str(val)

# EZ15--16: polynomialize all separated sixth/cube roots.
r,t,h,c=s.symbols('r t h c')
u=r**3; z=t**3; G=h**6; kap=c**6; P=kap*G
L=[r**4*t**2,c**2*r**2*t**4,h**6*c**4]
f=[r**2*t,c*r*t**2,h**3*c**2]
R=sum(([v*(L[0]-L[1]),v*(L[0]-L[2])] for v in f),[])
delta=u*u-P; gam=z*z-G; UV=u*u*z*z
M=u**4*z*z+kap*u*u*z**4+P**2*G-3*P*u*u*z*z
ck('EZ15_cubic',M-sum(R[i]**2-R[i]*R[i+1]+R[i+1]**2 for i in (0,2,4)))
ck('EZ16_uz_delta',u*z*delta-(r**3*R[0]+c**2*r*t**2*R[1]))
ck('EZ16_u2z_gamma',u*u*z*gam-(-h**6*R[0]+t**6*R[1]-h**3*t**3*R[4]))
ck('EZ16_kapuz_gamma',kap*u*z*gam-(c**2*r*t**2*(R[1]-R[0])-c**3*t**3*R[2]))
ck('EZ16_cross',UV*delta*gam-(t**6*R[1]**2-h**6*R[0]**2-h**3*t**3*R[0]*(R[4]+R[5])))
ck('EZ16_gamma_square',kap*UV*gam**2-(t**6*(2*R[0]**2-3*R[0]*R[1]+R[1]**2+R[2]**2)+h**3*t**3*R[0]*(R[4]-R[5])))

# IC12--13 are exact finite identities, including a unit second coordinate.
u,t,k,ss,a,b=s.symbols('u t k ss a b')
du=u-k; dt=t-ss
hc=a*(u*u-k*k)+b*(u*u*t*t-k*k*ss*ss)
ck('IC12',hc-(2*k*(a+b*ss*ss)*du+2*k*k*ss*b*dt+(a+b*t*t)*du**2+2*k*b*(t+ss)*du*dt+k*k*b*dt**2))
ck('IC13',hc-(du*(u+k)*(a+b*t*t)+dt*k*k*(t+ss)*b))

# JR8 and the exact real finite-penalty Schur identity.
t=s.symbols('t'); lam=2*(t*t+t+1)/(1+t)**2
ck('JR8',lam*(1+t*t)+2*t*t/(1+t)-s.Rational(1,2)*(1+t+lam*t/2)**2-(4*t+3)/(2*(1+t)**4))

# Endpoint E1/E2: retain every testing normal, the product root and D_t.
u,v,X,Y,Z,e,k,b,a,t=s.symbols('u v X Y Z e k b a t')
Rc=1+k*k*(1-a*e*e*b*b); T0=1+k*k-t*k**4*b*b
c0=1-t*k*k*b*b; U0=1+u*u; P0=T0+u*u
vg=e*k*b*b*u; Xg=-u*Z+t*vg; Yg=(1-a*e*e*b*b)*u
N1=X-Xg; N2=Yg-Y; dv=v-vg; Dt=Z-u*X
S=u+u*v*X-Y-v*Z+(T0-a)*e*e*b*b*u
soft=e*b*b*u*(e-k*Z); qq=e*b*b*(e*u*u-k**3*Z)
L=P0*soft-S+N2-Dt*dv+u*vg*N1-u*c0*qq
ck('E1',L)
xx=e*N1; yy=e*dv; ww2=e*k*Dt*dv**2
hh=k*k*soft-u*qq; rr=k*k*Dt*dv
rhs=k*u*u*(P0-k*k)*yy*(S-N2)-k*u*(P0+c0*k*k*u*u)*yy*qq
rhs+=u*u*(P0-k*k)*ww2-k**4*b*b*u*u*T0*xx*yy
ck('E2',P0*hh*rr-rhs)
out['all_zero']=all(v['zero'] for v in out.values())
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
