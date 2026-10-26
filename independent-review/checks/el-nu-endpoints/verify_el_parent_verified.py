import sympy as sp
u,v,W,Z,N,Ny=sp.symbols('u v W Z N Ny')
r,p,S,t,gam,c0,a0,m,D,alpha,beta,kappa,nu0,ell=sp.symbols('r p S t gam c0 a0 m D alpha beta kappa nu0 ell')
a=u*u;c=v*v; Du=1+a;Dv=1+c
C=c0/(S*(1-r));A=a0/(S*r);d=D+sp.Rational(1,2);mu=m*(1-p)
Y=r*u*W+t*gam*Ny;V=v*(W+u*Y);X=N-u*Z+p*V
P=m*(v*W-X+u*v*Y-u*Z)**2+C*Dv*(u*W-Y)**2+A*Du*Dv*Y**2+d*(u*W+u*v*X-Y-v*Z)**2+(alpha+sp.Rational(1,2))*Du*(v*Y-Z)**2+(beta+sp.Rational(1,2))*Dv*(X+u*Z)**2-sp.Rational(1,2)*Du*Dv*Z**2

def norm(expr):return sp.factor(sp.cancel(sp.expand(expr).subs(S,t**2).subs(gam**2,r*(1-r))))
def chk(tag,expr):
 x=norm(expr)
 print(tag, 'OK' if x==0 else 'FAIL '+str(x)[:600])
L0=mu*(1-p)+(beta+sp.Rational(1,2))*p**2;L1=(beta+sp.Rational(1,2))*p**2
Bc=1-r+p*c*(1+r*a)
AP=(1+c)*(c0*(1-r)+a0*r*(1+a))+S*(d*Bc**2+(alpha+sp.Rational(1,2))*(1+a)*r**2*c)
K=(1+a)*(d*Bc+(alpha+sp.Rational(1,2))*r)
FZ=(1+a)*(alpha+D*c+d*a*c)
P0=u*W/t;Q=t*Z
F0=c*(1+r*a)**2*(L0+L1*c)*W**2+AP*P0**2-2*v*K*P0*Q+FZ*Z**2
chk('EL3 full scalar graph',P.subs({N:0,Ny:0})-F0)
B0=(-mu+(beta+sp.Rational(1,2))*p*(1+c))/kappa
Fa=r*kappa*B0+d*(1-r)+d*p*c*(1+r*a)
Dt=d*(1+a)-nu0*S*Fa
L_N=v*B0*(kappa*W) +u*v*0  # use algebraic expression without conductors: L_N = ? via EA normal check E13 requires source jets.
YP=(1+c)*(-c0+a0*(1+a))+S*(c*(1+r*a)*(L0+L1*c)+d*(-1+p*a*c)*Bc+(alpha+sp.Rational(1,2))*(1+a)*r*c)
YK=nu0*YP+(1+a)*(D-alpha-d*p*a*c)
Gamma=P0-nu0*v*Q
Ly=gam*(YP*Gamma+v*YK*Q)
H_N=m+d*a*c+(beta+sp.Rational(1,2))*(1+c)
H_y=r*c0*(1+c)+(1-r)*a0*(1+a)*(1+c)+S*r*(1-r)*(d*(-1+p*a*c)**2+(alpha+sp.Rational(1,2))*(1+a)*c+a*c*(L0+L1*c))
H_Ny=t*gam*u*v*(kappa*B0-d+d*p*a*c)
pol=sp.expand(P)
chk('EL11 Y linear',pol.coeff(Ny,1).subs(N,0)/2-Ly)
chk('EL12 N hard',pol.coeff(N,2)-H_N)
chk('EL12 Ny hard',pol.coeff(Ny,2)-H_y)
chk('EL12 mixed hard',pol.coeff(N*Ny)-2*H_Ny)
# Raw N coefficient based on native source p, kappa relation as needed
LN=pol.coeff(N,1).subs(Ny,0)/2
# derive raw LN target using graph expressed in auxiliary rows E12; E13 eq depends on contact and divisibility
# Could express elementary LN first via N-dependent rows E13 with note missing contact identities
print('raw L_N width',len(str(sp.factor(LN))))

b,aa,cc=sp.symbols('b aa cc', nonzero=True)
q=b**2
x=kappa*S
j=nu0*b
BB0=lambda C:(-mu+(beta+sp.Rational(1,2))*p*(1+C))/kappa
FFa=lambda AA,CC:r*kappa*BB0(CC)+d*(1-r)+d*p*CC*(1+r*AA)
DTT=lambda AA,CC:d*(1+AA)-nu0*S*FFa(AA,CC)
Da=sp.cancel((DTT(a,c)-DTT(x,c))/(a-x))
Dc=sp.cancel((DTT(x,c)-DTT(x,q))/(c-q))
B01=sp.cancel((BB0(c)-BB0(q))/(c-q))
E12=(v/b)*(kappa*W-j*u*Z)
Rv=v*(a-x)*Z
Vh=v*(v-b)*u*Z
P0=u*W/t;Q=t*Z
Gamma=P0-nu0*v*Q
LnEL=b*BB0(c)*E12+u*v*t*FFa(a,c)*Gamma-u*v*Da*Rv +(-DTT(x,q)+nu0*b*B01*(v+b)-v*(v+b)*Dc)*Vh
# assume S=t**2 but generic b relation; E13 may require source stationarity.
chk('EL13 modulo native contact jet',LN-LnEL - u*v*Z*b*(nu0*BB0(q)-DTT(x,q)))
