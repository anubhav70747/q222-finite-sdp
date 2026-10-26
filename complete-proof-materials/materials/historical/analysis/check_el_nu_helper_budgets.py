import sympy as s
u,v,w,W,X,Y,Z,S,k,ell,nu,b,zeta,x,q=s.symbols('u v w W X Y Z S k ell nu b zeta x q', nonzero=True)
checks=[]
def zero(name,expr):
    r=s.factor(expr)
    assert r==0,(name,r)
    checks.append(name)
# Native core: no source specialization, full four-state identity.
Du,Dv=1+u*u,1+v*v
J=u*W+u*v*X-Y-v*Z
K=J**2+Du*(v*Y-Z)**2+Dv*(X+u*Z)**2-Du*Dv*Z**2
Bt,Rt,Vt=X+u*Z,u*W-Y,v*Z
zero('CE2 entire four-state signed core',K-((Rt+u*v*Bt)**2+Dv*Bt**2+Du*v*v*Y*Y+Du*u*u*Vt*Vt-2*Du*u*Vt*(W+v*Bt)))
# NU exact root identities with S=s0^2 and ell=k*zeta.
s0=s.symbols('s0',nonzero=True)
delta=w*w-q
N=W-zeta*u*w*Z
G=u*W/s0-ell*s0*w*Z
B=s0*delta*Z
Bu=u*delta*Z
subs={ell:k*zeta,x:k*s0*s0}
zero('N10 complete mixed radial identity',(w*w*(u*u-x)*delta*Z*Z-(G*w*B-N*w*Bu)/zeta).subs(subs))
zero('N13 first helper budget',(delta*delta*N*N+2*zeta*delta*w*G*B+2*zeta*ell*w*w*B*B-zeta*zeta*(w*Bu)**2-(delta*W)**2).subs(subs))
zero('N13 second helper budget',(delta*delta*G*G+2*ell*delta*w*N*Bu+2*ell*zeta*w*w*Bu*Bu-ell*ell*(w*B)**2-(delta*u*W/s0)**2).subs(subs))
# EL exact radial identities: sqrt(x)=a0, k=a0^2/s0^2, ell=nu*b^2.
a0=s.symbols('a0',nonzero=True)
ke=a0*a0/(s0*s0)
j=nu*b
P,Q=u*W/s0,s0*Z
Gamma=P-nu*v*Q
Lambda=v*P-nu*b*b*Q
d=v-b
Bo=d*Q
G0=P-j*Q
E12=(v/b)*(ke*W-j*u*Z)
E23=nu*(u*v-b*a0)*Z
Be=s0*(v*v-b*b)*Z
Vh=v*d*u*Z
Rv=v*(u*u-a0*a0)*Z
R0=nu*b*b*s0*(u*u-a0*a0)*Z
eta=ke*(v+b)/j
zero('E13 complete radial allocation',v*v*(u*u-a0*a0)*(v*v-b*b)*Z*Z-((v*v-b*b)*E23*E23/(nu*nu)+2*a0*b*E23*Be/(nu*s0)-ke*Be*Be))
cost=d*d*E12*E12/(nu*nu)+eta*(Lambda+b*Gamma)**2+(eta+2*ke)*(Lambda-v*Gamma)**2/(nu*nu)
rem=(v*d*ke*W/j)**2+eta*((v+b)**2-v*v)*(G0*G0+Bo*Bo)+(ke*v/j)*(v*G0-v*Bo)**2+(ke*b/j)*(v*G0+v*Bo)**2+2*ke*((v+b)**2-v*v)*Bo*Bo
zero('E15 full separated SOS remainder',cost-Vh*Vh-rem)
zero('E16 Rv root readback',Rv-(ke*v*s0*Gamma/(nu*b)-u*E12/nu+a0*a0*v*(v-b)*Z/b))
zero('E16 B budget exact remainder',a0**4*Be*Be/(b*b*s0*s0)-(a0*a0*v*(v-b)*Z/b)**2-a0**4*(2*v*b+b*b)*(v-b)**2*Z*Z/(b*b))
zero('E16 R0 budget exact remainder',2*nu*nu*s0*s0*v*v*Rv*Rv+2*nu*nu*(u*u-a0*a0)**2*Be*Be-R0*R0-(nu*s0*(2*v*v-b*b)*(u*u-a0*a0)*Z)**2)
for name in checks:print('PASS',name)
print('All',len(checks),'exact identities passed; native contact bounds and full-parent admissions remain separate proof obligations.')

