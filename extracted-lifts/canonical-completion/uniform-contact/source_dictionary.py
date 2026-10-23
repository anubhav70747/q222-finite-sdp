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

u,v=s.symbols('u v', nonnegative=True)
t,r,w,zeta=s.symbols('t r w zeta', real=True)
