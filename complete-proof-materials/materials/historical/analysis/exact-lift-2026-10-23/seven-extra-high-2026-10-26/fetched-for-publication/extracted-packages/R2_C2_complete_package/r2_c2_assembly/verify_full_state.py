"""Generic exact four-state source decomposition and full hard coefficient check."""
import sympy as sp
u,v,z,N,M=sp.symbols('u v z N M')
p,H,mu,n,C,A,D,al=sp.symbols('p H mu n C A D alpha')
DU=1+u*u; DV=1+v*v
X=N+p*v-u*z; Z=z+u*N+p*u*v
W=1+u*(M+H*v*z);Y=u-M-H*v*z
m=mu/(1-p)
R_m=v*W-X+u*v*Y-u*Z
R_D=u*W+u*v*X-Y-v*Z
P=m*R_m**2+n*DU*DV*X**2+C*DV*(u*W-Y)**2+A*DU*DV*Y**2+(D+sp.Rational(1,2))*R_D**2+(al+sp.Rational(1,2))*DU*(v*Y-Z)**2+sp.Rational(1,2)*DV*(X+u*Z)**2-sp.Rational(1,2)*DU*DV*Z**2
G=sp.cancel((P-m*p*DU**2*N**2)/DU)
G0=mu*DU*(N**2-2*v*N+(1-p)*v**2)+n*DV*(N+p*v-u*z)**2+C*DU*DV*(M+H*v*z)**2+A*DV*(u-M-H*v*z)**2+(D+sp.Rational(1,2))*DU*(M-(1-H)*v*z)**2+(al+sp.Rational(1,2))*(u*(1-p)*v-(1+H*v*v)*z-u*N-v*M)**2+sp.Rational(1,2)*DU*DV*(N+p*v)**2-sp.Rational(1,2)*DV*(z+u*N+p*u*v)**2
assert sp.expand(G-G0)==0
print('RFN7 generic full-state equality PASS',flush=True)
HN=sp.expand(G0).coeff(N,2).subs(M,0)
HM=sp.expand(G0).coeff(M,2).subs(N,0)
HNm= (mu+sp.Rational(1,2))*DU+v*v/sp.Integer(2)+al*u*u+n*DV
HMm=DU*(C*DV+D+sp.Rational(1,2))+A*DV+(al+sp.Rational(1,2))*v*v
assert sp.expand(HN-HNm)==0 and sp.expand(HM-HMm)==0
assert sp.expand(sp.diff(G0,N,M)/2-(al+sp.Rational(1,2))*u*v)==0
print('RFN9 2x2 hard matrix PASS',flush=True)
# Full N first derivative direct and shifted identities: general scalar equation RFN22 uses source-specific stationarity
LN=sp.diff(G0,N).subs({N:0,M:0})/2
LM=sp.diff(G0,M).subs({N:0,M:0})/2
# derived unconditional derivative and shifted second-normal coefficient identity:
Kd=C*H-(D+sp.Rational(1,2))*(1-H)
Cq=Kd+A*H+sp.Rational(1,2)+v*v*(C*H*DU+A*H+(al+sp.Rational(1,2))*H)
assert sp.expand(Cq-H*HMm+D+H*(C+D+sp.Rational(1,2))*u*u)==0
print('RL23 second-normal shift cancellation independent of stationary source PASS',flush=True)
print('LN expression:',str(sp.factor(LN))[:3500])
print('LM expression:',str(sp.factor(LM))[:3500])
