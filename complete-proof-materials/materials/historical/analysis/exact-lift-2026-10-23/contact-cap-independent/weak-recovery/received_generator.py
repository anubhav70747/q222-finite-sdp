import sympy as s
from fractions import Fraction as F
h=s.Rational(1,2)
B,T,j,p,H=s.symbols('B T j p H', positive=True)
n,mu,m,C,A,D,alpha,fc,qc,rc,gc,kc,dc=s.symbols(
    'n mu m C A D alpha fc qc rc gc kc dc', real=True)
x=j*B; y=T*B; qz=j*T
K=j+1-p; E=1+B*(H*T+p-1); W0=(1-H)*j+1-p
P=(j+1+B*((j+1)**2-p)
   +T*B*(H*j+j+1+B*(H*j*(2*j+2-p)+p*(2*j+1-p)))
   +H*j*T*T*B*B*(1+B*p+B*H*j))
L=K+T*(W0-H)+B*T*(p*K+H*j*(1-p))+B*H*T*T*(W0-1)
NN=(p+B*(-j*j+3*j*p-j-p*p+p)
    +B*T*(H*j*j+H*j-j*j+j*p-j+p)
    +B*B*T*(2*H*j*j*p-H*j*j-H*j*p*p+2*H*j*p
             -2*j*j*p+3*j*p*p-j*p-p**3+p*p)
    +B*B*H*T*T*j*(H*j-j+p)*(1+B*p))
F1=(-2*B*H*j*p+B*H*j+2*B*j*p-2*B*p*p+2*B*p
    -H*j-H*p-H+j+1-p*p/j+p/j)
F2=(-B*B*H*H*j*p+B*B*H*j*p-B*B*H*p*p+B*B*H*p
    -B*H*H*j-B*H*H*p+B*H*j-B*H*p*p/j+B*H*p/j)
FF=K+T*F1+T*T*F2
KK=(p*W0+B*(1-p)*(j*(H*(p+1)-p)-(1-p)*(p-H))
    +T*B*H*((p+1)*(1-H)*j+(1-p)*(B*H*j+p+1-H))
    +T*T*B*B*H*H*j*(1-H))
src={
 mu:j*p*(1+y)*(1+x)*(1+H*y)/(2*P),
 n:NN/(2*(j-p)*P),
 C:L/(2*H*T*P),
 A:T*(1+x)*(1+H*y)*(K+H*j*y)/(2*(1-H*T)*P),
 D:(1+H*y)*(K+H*j*y)*(1-y*((1-H)*j-p))/(2*T*(1-H)*P),
 alpha:B**2*(K+H*j*y)*L/(2*E*P),
 fc:T*x*x*FF/(2*P),
 rc:B*B*j*T*(1-p)*(1+x)*(1+H*y)*L/(2*E*P),
 gc:H*T*T*j*j*B**3*L/(2*P),
 kc:-T*T*B*B*j*(1+x)*KK/(2*E*P)}
src[m]=src[mu]/(1-p)
src[qc]=src[A]*x-src[fc]
src[dc]=H**2*(src[A]+src[C]+src[alpha]+h)*y*y*qz-src[kc]
# Native admissibility includes src[n]>=0, equivalently NN>=0.

h0=mu+h
cb=n*j/(n+h0); Me=n*h0*j*j/(n+h0)
a=1/s.sqrt(n+h0); b=s.sqrt(n/(h0*(n+h0)))
tau=s.sqrt(T*H*(1-H))
rhoQ=T*s.sqrt(x/qc); rhoR=s.sqrt(qz/rc)
Sigma=Me+fc/y; gamma=p*p/Sigma; kap=-kc
Kc=p/2+j*(H*(alpha+h)-h)+n*(p-j)
KR=-mu-(alpha+h)*(1-p)+(T/j)*Kc
kR=KR+alpha/B

# Exactly the saved 44 actual upper-triangular coefficient records.
# Each tuple is (row,column,u_exponent,v_exponent,coefficient).
rec=[]
def add(ii,jj,eu,ev,c):
    rec.append((min(ii,jj),max(ii,jj),eu,ev,s.sympify(c)))
diag=[1,tau**2*(C+D+h+A),1+gc/qc,1,1+(gc-2*kap)/rc,
      1-gc/fc-64*x*x*y/fc,1-gc/fc-64*x*x*y/fc,
      1-gc/fc-64*x*x*y/fc,16,8,1]
for ii,c in enumerate(diag,1): add(ii,ii,0,0,c)
for ii,jj,c in [(1,9,a*kR),(4,9,-b*kR),(1,10,a/2),(4,10,-b/2),
               (2,3,tau*rhoQ*A/T),(2,11,(mu/p)*s.sqrt(H*(1-H)))]:
    add(ii,jj,0,0,c)
add(1,5,1,0,a*alpha*rhoR); add(4,5,1,0,-b*alpha*rhoR)
add(2,5,0,1,tau*alpha*rhoR)
add(3,3,2,0,-64*x*y/qc); add(5,5,2,0,-64*x*y/rc)
for ii in [6,7,8]: add(ii,ii,2,0,-gamma*T*y)
add(3,3,0,2,-2*kap/(y*qc))
add(2,3,0,2,tau*rhoQ*H*(C+A+alpha+h))
for ii in [6,7,8]: add(ii,ii,0,2,-8*gamma*y)
for ii,jj,c in [(1,3,a*(Kc/j)*rhoQ),(3,4,-b*(Kc/j)*rhoQ),
               (1,2,(alpha+h)*tau*a),(2,4,-(alpha+h)*tau*b)]:
    add(ii,jj,1,1,c)
add(2,5,2,1,-tau*B*T*C*H*rhoR)
add(2,3,2,2,tau*rhoQ*H*C)
for ii in [4,6,7,8]: add(ii,ii,2,2,-gamma*T)
for ii in [4,6,7,8]: add(ii,ii,0,4,-8*gamma)
assert len(rec)==44

u,v=s.symbols('u v', nonnegative=True)
t,r,w,zeta=s.symbols('t r w zeta', real=True)
Qab={e:s.zeros(11) for e in {(eu,ev) for _,_,eu,ev,_ in rec}}
for ii,jj,eu,ev,c in rec:
    Qab[eu,ev][ii-1,jj-1]+=c
    if ii!=jj: Qab[eu,ev][jj-1,ii-1]+=c
Q=s.zeros(11)
for (eu,ev),mat in Qab.items(): Q+=u**eu*v**ev*mat

# The eleven rows, before their common division by sqrt(1+v^2).
N0=t-p*v*w; M0=-r-H*v*zeta
bcur=v*w-u*zeta/j; q=v*zeta-T*u*w
R=zeta-u*v*w/B; Ru=v*(u*u-x)*w
U=u/s.sqrt(x); V=v/s.sqrt(y); Zs=zeta/s.sqrt(qz)
Qh=U*w-V*Zs; Bh=V*w-U*Zs; Rh=Zs-U*V*w
Psi=(U-V)**2*w*w+Zs**2*(1-U*V)**2+2*U*V*(w-Zs)**2
f=s.Matrix([
 s.sqrt(n+h0)*(N0+cb*bcur), M0/tau, q/rhoQ, s.sqrt(Me)*bcur,
 R/rhoR, s.sqrt(fc)*(U-V)*w,
 s.sqrt(fc/qz)*(1-U*V)*zeta,
 s.sqrt(2*fc)/(x*y)**s.Rational(1,4)*s.sqrt(u*v)*(w-Zs),
 Ru, p*v*v*bcur, s.sqrt(T)*p*u*v*bcur])

# Original RFN18b and RFN19, homogenized without a state denominator.
E18=Qh**2+Rh**2+V**2*(w-U*V*Zs)**2
W19=4*Psi+4*U**2*Qh**2+4*U**2*Rh**2-V**2*(U**2-1)**2*w*w
# RFN18b: E18=(1+V^2)*Psi-V^2*Bh^2.
# RFN19 supplies the unchanged positive rational 7-row/3-row SOS for W19.

# ALL retained whole squares, including the original mp*(1+u^2)*N0^2.
mp=mu*p/(1-p)
RP=(mp*N0**2+(mp+h0+alpha)*(u*N0)**2+h*(v*N0)**2
    +n*(v*(N0+p*bcur))**2
    +tau**2*(C+D+h)*(u*f[1])**2
    +tau**2*(C+A+alpha+h)*(v*f[1])**2
    +tau**2*C*(u*v*f[1])**2
    +dc*(V*Qh)**2+h*p*p*y*y*(V*Bh)**2
    +gc*(V*(w-U*V*Zs))**2+kap*(V*Qh-Rh)**2)
certificate_numerator=((f.T*Q*f)[0]+RP+16*x*x*y*W19
                       +8*gamma*fc*y*v*v*E18
                       +gamma*T*fc*y*u*u*E18)
# F_h/4 = certificate_numerator/(1+v^2).

# Complete saved RFN8/RFN12 parent for independent comparison.
Fs=(qc*Qh**2+n*j*j*y*Bh**2+rc*Rh**2+dc*V**2*Qh**2
    +(n*p*p*y*y+h*p*p*y*y-gc)*V**2*Bh**2
    +fc*Psi+gc*V**2*Psi+kc*Zs**2*(V**2-1)**2)
LN=(n*j*bcur+n*p*v*v*bcur-j*alpha*bcur+h*p*v*v*bcur
    +KR*Ru+(Kc/j)*u*v*q)
Cq=A/T+v*v*(C*H*(1+u*u)+A*H+(alpha+h)*H)
LM=Cq*q+mu*u*v*bcur+v*(alpha-B*T*C*H*u*u)*R
HN=h0*(1+u*u)+h*v*v+alpha*u*u+n*(1+v*v)
HM=(1+u*u)*(C*(1+v*v)+D+h)+A*(1+v*v)+(alpha+h)*v*v
parent_numerator=(Fs+2*LN*N0+2*LM*M0+HN*N0**2+HM*M0**2
                  +2*(alpha+h)*u*v*N0*M0+mp*(1+u*u)*N0**2)

# Literal fixed gate prescription: no free thresholds.
rb={(1,0):1,(0,1):1,(2,0):512,(0,2):8,
    (1,1):16,(2,1):4,(2,2):16,(0,4):128}
I=s.eye(11); reserve=Qab[0,0].copy(); gates={}
for (eu,ev),bound in rb.items():
    scale=s.Rational(1,2**(12*eu+25*ev))
    budget=bound*scale; Ac=scale*Qab[eu,ev]
    reserve-=budget*I
    gates[eu,ev]=((budget*I+Ac)/2,(budget*I-Ac)/2)
# Q=reserve + sum Gplus*(1+p_ab)+Gminus*(1-p_ab),
# p_ab=(4096*u)^a*(33554432*v)^b.
# Use expr.subs(src,simultaneous=True) for the native source witness.
for row in rec: print(row)
