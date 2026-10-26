import sympy as S
u,v,w,N,MB=S.symbols('u v w N MB',real=True)
p,c,E,mu1,mu2,n1,n2=S.symbols('p c E mu1 mu2 n1 n2',real=True)
C,A,D,Ah,bd,H,Km,h,al=S.symbols('C A D Ah bd H Km h al',real=True)
q=1-c; du=1+u*u;dv=1+v*v
A0=v*w; B0=v; C0=E*u*w
F0=p*p/2-n2-2*n1*p+mu2; AN=p-2*mu1-2*n1
vzero=1-p*u*v*w
X=N+p*du*A0-u
qb=v*(du-u*N-p*u*du*v*w)
T=MB+c*qb;J=MB-q*qb
Y=u*w-T/du; uW=u*w+u*u*T/du
L=A*dv+H*v*v
G_RI=(Km*N*N-2*mu1*du*A0*N+mu2*du**2*A0**2-2*n1*du*dv*X*A0-n2*du*dv*A0**2+C*dv*T**2+A*du*dv*Y**2+(D+S.Rational(1,2))*J**2+H*du*(v*Y-1)**2+S.Rational(1,2)*dv*(X+u)**2-S.Rational(1,2)*du*dv)
G_RU=(-2*mu1*du*A0*N+mu2*du**2*A0**2-2*n1*du*dv*X*A0-n2*du*dv*A0**2+C*dv*T**2+Ah*du*dv*(2*uW*Y-h*(uW)**2)-2*bd*MB*qb+bd*q*qb**2+J**2/2+H*du*(v*Y-1)**2+S.Rational(1,2)*dv*(X+u)**2-S.Rational(1,2)*du*dv)
# scalar raw expansions, using physically alpha=H-1/2, and Dcal^2 = s^3 b^2 t
fAA=du**2*dv*F0
fBB_RI=(C*c*c*du**2*dv+(D+S.Rational(1,2))*q*q*du**2+du*L*c*c)*vzero**2+2*H*du*c*vzero-du/2-du**2*mu2*v*v*w*w
fBC_RI=(-2*du*L*c*vzero-2*H*du+2*n1*du*dv)/E
fCC_RI=(du*L+du*dv*n2*v*v)/E**2
RI_model=fAA*A0*A0+fBB_RI*B0*B0+fBC_RI*B0*C0+fCC_RI*C0*C0+(H-S.Rational(1,2))*du
fBB_RU=(C*c*c*du**2*dv+(bd*q+q*q/2)*du**2+H*du*v*v*c*c -Ah*du*dv*(2*u*u+h*u**4)*c*c)*vzero**2+2*H*du*c*vzero-du/2-du**2*mu2*v*v*w*w
fBC_RU=(2*Ah*du*dv*( (1-h)*u*u-1)*c*vzero-2*H*du*v*v*c*vzero-2*H*du+2*n1*du*dv)/E
fCC_RU=(Ah*(2-h)*du*dv+H*du*v*v+du*dv*n2*v*v)/E**2
RU_model=fAA*A0*A0+fBB_RU*B0*B0+fBC_RU*B0*C0+fCC_RU*C0*C0+(H-S.Rational(1,2))*du

R_A=du*AN # multiply A0, not calA
RI_RC=2*c*u*v*L/E
RU_RC=2*c*u*v*(Ah*dv*(1-(1-h)*u*u)+H*v*v)/E
RI_RB=du*(p-2*n1)*v*v*w -2*u*v*(C*c*c*du*dv+A*c*c*dv+(D+S.Rational(1,2))*q*q*du+H*c*c*v*v)*vzero -2*H*c*u*v
RU_RB=du*(p-2*n1)*v*v*w -2*u*v*(C*c*c*du*dv+(bd*q+q*q/2)*du+H*c*c*v*v-Ah*c*c*dv*(2*u*u+h*u**4))*vzero-2*H*c*u*v
RI_SC=-2*L/E
RI_SB=(2*C*c*du*dv-2*(D+S.Rational(1,2))*q*du+2*A*c*dv+2*H*c*v*v)*vzero+2*H
RU_SC=2*(Ah*dv*((1-h)*u*u-1)-H*v*v)/E
RU_SB=(2*C*c*du*dv-(2*bd+q)*du-2*Ah*c*dv*(2*u*u+h*u**4)+2*H*c*v*v)*vzero+2*H

def coef(poly,i,j):return S.expand(poly).coeff(N,i).coeff(MB,j)
checks={}
def check(label, expr):
    zz=S.factor(S.together(expr));checks[label]=zz==0
    if zz!=0:print('FAIL',label, str(zz)[:1000])

check('RI scalar',coef(G_RI,0,0)-RI_model)
check('RU scalar',coef(G_RU,0,0)-RU_model)
check('RI N',coef(G_RI,1,0)-(R_A*A0+RI_RB*B0+RI_RC*C0))
check('RI MB',coef(G_RI,0,1)-(RI_SB*B0+RI_SC*C0))
check('RU N',coef(G_RU,1,0)-(R_A*A0+RU_RB*B0+RU_RC*C0))
check('RU MB',coef(G_RU,0,1)-(RU_SB*B0+RU_SC*C0))

RiNN=Km+dv/S.Integer(2)+C*dv*c*c*u*u*v*v+A*dv/du*c*c*u*u*v*v+(D+S.Rational(1,2))*q*q*u*u*v*v+H/du*c*c*u*u*v**4
RiNM=-C*dv*c*u*v-A*dv/du*c*u*v+(D+S.Rational(1,2))*q*u*v-H/du*c*u*v**3
RiMM=C*dv+A*dv/du+(D+S.Rational(1,2))+H*v*v/du
RuNN=dv/2+C*dv*c*c*u*u*v*v-Ah*dv/du*(2*u*u+h*u**4)*c*c*u*u*v*v+(bd*q+q*q/2)*u*u*v*v+H/du*c*c*u*u*v**4
RuNM=-C*dv*c*u*v+Ah*dv/du*(2*u*u+h*u**4)*c*u*v+(bd+q/2)*u*v-H/du*c*u*v**3
RuMM=C*dv-Ah*dv/du*(2*u*u+h*u**4)+S.Rational(1,2)+H*v*v/du
for name,G,hNN,hNM,hMM in [('RI',G_RI,RiNN,RiNM,RiMM),('RU',G_RU,RuNN,RuNM,RuMM)]:
    check(name+' hard NN',coef(G,2,0)-hNN)
    check(name+' hard NM',coef(G,1,1)-2*hNM)
    check(name+' hard MM',coef(G,0,2)-hMM)
print('PASS',sum(checks.values()),'/',len(checks))
if not all(checks.values()):raise SystemExit(1)
