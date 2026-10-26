import sympy as sp
u,v,w,W,Z=sp.symbols('u v w W Z')
t,kappa,nu0,b,ell,zeta=sp.symbols('t kappa nu0 b ell zeta',nonzero=True)
S=t*t;q=b*b;j=nu0*b;x=kappa*S
# NU budgets
Dlt=w*w-q
N=W-zeta*u*w*Z;G=u*W/t-ell*t*w*Z;B=t*Dlt*Z;Bu=u*Dlt*Z;P=u*W/t
N1=Dlt**2*N**2+2*zeta*Dlt*w*G*B+2*zeta*ell*w*w*B**2-zeta**2*(w*Bu)**2-(Dlt*W)**2
N2=Dlt**2*G**2+2*ell*Dlt*w*N*Bu+2*ell*zeta*w*w*Bu**2-ell**2*(w*B)**2-(Dlt*P)**2
print('NU10.1',sp.factor(N1)==0,flush=True)
print('NU10.2',sp.factor(N2)==0,flush=True)
# NU15, x=ell*S/zeta
mixed=w*w*(u*u-ell*S/zeta)*(w*w-q)*Z**2-(G*(w*B)-N*(w*Bu))/zeta
print('NU15',sp.factor(mixed)==0,flush=True)
# EL helper identities
# relations ell=nu0*b**2, zeta=ell/kappa.
P=u*W/t;Q=t*Z;T=kappa*W;U=u*Z
Gamma=P-nu0*v*Q;Lambda=v*P-ell*Q;delta=v-b
G0=P-j*Q;Bo=delta*Q;E12=(v/b)*(T-j*U);Vh=v*delta*U
eta=kappa*(v+b)/j

def red(z):
 z=sp.cancel(sp.together(z))
 num,den=sp.fraction(z)
 num=sp.expand(num).subs(ell,nu0*b*b)
 return sp.factor(num)
print('EL7.1',red((v+b)*G0-(Lambda+b*Gamma))==0,flush=True)
print('EL7.2',red(nu0*(v+b)*Bo-(Lambda-v*Gamma))==0,flush=True)
C=delta**2/nu0**2*E12**2+eta*(Lambda+b*Gamma)**2+(eta+2*kappa)/nu0**2*(Lambda-v*Gamma)**2
rhs=(v*delta*kappa*W/j)**2+eta*((v+b)**2-v*v)*(G0**2+Bo**2)+(kappa*v/j)*(v*G0-v*Bo)**2+(kappa*b/j)*(v*G0+v*Bo)**2+2*kappa*((v+b)**2-v*v)*Bo**2
print('EL8 budget exact',red(C-Vh**2-rhs)==0,flush=True)
# EL9, EL10 and EL11
Rv=v*(u*u-x)*Z
B=t*(v*v-b*b)*Z
E9=kappa*v*t/(nu0*b)*Gamma-u/nu0*E12+(x*v*(v-b)/b)*Z
print('EL9',red(Rv-E9)==0,flush=True)
E10=x*x/(q*S)*B**2-(x*v*(v-b)*Z/b)**2-(x*x/q)*(2*v*b+b*b)*(v-b)**2*Z**2
print('EL10',red(E10)==0,flush=True)
R0=ell*t*(u*u-x)*Z
E11=2*nu0**2*S*v*v*Rv**2+2*nu0**2*(u*u-x)**2*B**2-R0**2-(nu0*t*(2*v*v-q)*(u*u-x)*Z)**2
print('EL11',red(E11)==0,flush=True)
# EL18 mixed radial identity
E23=nu0*(u*v-sp.sqrt(x*q))*Z
# substituting sqrt(xq)=sqrt(kS q) via separate positive symbols requires formal b,t,sqrt(kappa).
kk=sp.symbols('kk',nonzero=True)
E23=nu0*(u*v-kk*t*b)*Z
expr=v*v*(u*u-kappa*S)*(v*v-q)*Z*Z - (v*v-q)/nu0**2*E23**2 - 2*kk*b/nu0*E23*B +kappa*B**2
expr=expr.subs(kappa,kk*kk)
print('EL18',red(expr)==0,flush=True)
