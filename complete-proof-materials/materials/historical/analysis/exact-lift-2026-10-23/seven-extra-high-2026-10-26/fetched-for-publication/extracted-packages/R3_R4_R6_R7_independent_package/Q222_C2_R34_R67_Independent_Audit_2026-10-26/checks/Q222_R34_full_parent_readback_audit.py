import sympy as s
u,v,z,N,M=s.symbols('u v z N M')
p,H,mu,n,C,A,D,al,m=s.symbols('p H mu n C A D al m')
Du=1+u*u; Dv=1+v*v
X=N+p*v-u*z; Z=z+u*N+p*u*v; W=1+u*(M+H*v*z);Y=u-M-H*v*z
raw=m*(v*W-X+u*v*Y-u*Z)**2+n*Du*Dv*X**2+C*Dv*(u*W-Y)**2+A*Du*Dv*Y**2+(D+s.Rational(1,2))*(u*W+u*v*X-Y-v*Z)**2+(al+s.Rational(1,2))*Du*(v*Y-Z)**2+s.Rational(1,2)*Dv*(X+u*Z)**2-s.Rational(1,2)*Du*Dv*Z**2
expr=(raw-m*p*Du**2*N*N)/Du
expected=mu*Du*(N**2-2*v*N+(1-p)*v**2)+n*Dv*(N+p*v-u*z)**2+C*Du*Dv*(M+H*v*z)**2+A*Dv*(u-M-H*v*z)**2+(D+s.Rational(1,2))*Du*(M-(1-H)*v*z)**2+(al+s.Rational(1,2))*(u*(1-p)*v-(1+H*v*v)*z-u*N-v*M)**2+s.Rational(1,2)*Du*Dv*(N+p*v)**2-s.Rational(1,2)*Dv*(z+u*N+p*u*v)**2
# m(1-p)=mu, symbolic no division by opening in grammar
assert s.cancel(expr-expected.subs(mu,m*(1-p)))==0
print('PASS full R3/R4 original parent identity (9)')
HN=(mu+s.Rational(1,2))*Du+s.Rational(1,2)*v*v+al*u*u+n*Dv
HM=Du*(C*Dv+D+s.Rational(1,2))+A*Dv+(al+s.Rational(1,2))*v*v
HNM=(al+s.Rational(1,2))*u*v
for name,test in [('HN',s.diff(expected,N,2)/2-HN),('HM',s.diff(expected,M,2)/2-HM),('HNM',s.diff(expected,N,M)/2-HNM)]:
 assert s.cancel(test)==0,(name,s.cancel(test))
 print('PASS',name)
# Full RSR21/22 derivatives as polynomial identities up to exact source balances.
j,B,T=s.symbols('j B T'); x=j*B; y=T*B
b0=v-u*z/j; q0=v*z-T*u; R0=z-u*v/B; Ru=v*(u*u-x)
Kc=p/s.Integer(2)+j*(H*(al+s.Rational(1,2))-s.Rational(1,2))+n*(p-j)
KR=-mu-(al+s.Rational(1,2))*(1-p)+(T/j)*Kc
Kd=C*H-(D+s.Rational(1,2))*(1-H)
Cq=Kd+A*H+s.Rational(1,2)+v*v*(C*H*(1+u*u)+A*H+(al+s.Rational(1,2))*H)
Ln=(j*(n-al)+p*(n+s.Rational(1,2))*v*v)*b0+KR*Ru+(Kc/j)*u*v*q0
Lm=Cq*q0+mu*u*v*b0+v*(al-B*T*C*H*u*u)*R0
SN=n*p-mu+p/s.Integer(2)+(-mu-(al+s.Rational(1,2))*(1-p))*x+j*(al-n)+y*Kc
S0=T*(Kd+A*H+s.Rational(1,2))-A
S1=j*Kd+mu+x*T*C*H
S2=j*Kd+x*T*C*H+T*(C*H+A*H+(al+s.Rational(1,2))*H)-A-(al+s.Rational(1,2))*(1-p)+al/B
Dln=s.expand(s.diff(expected,N).subs({N:0,M:0})/2-Ln)
Dlm=s.expand(s.diff(expected,M).subs({N:0,M:0})/2-Lm)
# Original states difference equals v*S_N, u*S0+u²v z*S1/j+uv²(S2-S1), given source balances.
for fac in [1,2]:
 dr=s.cancel(Dln-fac*v*SN)
 print('Dln with factor',fac,'zero',dr==0)
 if dr==0:break
for fac in [1,2]:
 dr=s.cancel(Dlm-fac*(u*S0+u*u*v*z*S1/j+u*v*v*(S2-S1)))
 print('Dlm with factor',fac,'zero',dr==0)
 if dr==0:break
