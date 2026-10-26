"""Exact generic strong-II full-parent and packet readback, no physical-source approximations."""
import sympy as s
u,v,W,xi,eta,ze,z,r,B,k,t,C,D,n,m,A,alpha,ellstar,alphastar=s.symbols('u v W xi eta ze z r B k t C D n m A alpha ellstar alphastar')
y=k*k;U=u*u;V=v*v;Du=1+U;Dv=1+V
H=B-y*t;j=1+z*r-t;cc=C*y*t;dd=D*y*j;ss=(1+z*r)/B
ell=y*(C+D);q=D/(C+D);g=y*C*D/(C+D)
XX= -z*r*v/B*W+xi
YY=H*u/B*W+eta
ZZ=u*v/B*W+ze
Mm=v*W-XX+u*v*YY-u*ZZ
JJ=u*W+u*v*XX-YY-v*ZZ
QQ=u*W-YY
F=(4*n*Du*Dv*XX**2+4*m*Mm**2+4*C*Dv*QQ**2+4*A*Du*Dv*YY**2+4*D*JJ**2+4*alpha*Du*(v*YY-ZZ)**2+2*(JJ**2+Du*(v*YY-ZZ)**2+Dv*(XX+u*ZZ)**2-Du*Dv*ZZ**2))
Qn=u*xi-ze
Rh=(v*Qn-ss*u*(V-y)*W)/k
NN=eta/k-q*Rh
K=u*k*W
first=4*C*QQ**2+4*D*JJ**2
first_r=4*ell*NN**2+4*g*Rh**2+8*(dd-cc)/B*K*NN-8*g*ss*K*Rh+4*(cc*t+dd*j)/B**2*K**2
third=4*C*V*QQ**2
third_r=4*C*V*eta**2+4*cc*t/B**2*V*K**2-8*cc/B*V*K*(NN+q*Rh)
for name,lhs in [('LII4-main',first-first_r),('LII4-CV',third-third_r)]:
    print(name,flush=True)
    assert s.cancel(lhs)==0,(name,s.factor(lhs))
    print('PASS',flush=True)
Pi=4*(ell-ellstar)*NN**2+4*C*V*eta**2+4*n*Du*Dv*xi**2+4*(alpha-alphastar)*Du*(v*eta-ze)**2
G=s.expand(F-Pi)
Gr=4*ellstar*NN**2+4*g*Rh**2+8*(dd-cc)/B*K*NN-8*g*ss*K*Rh+4*(cc*t+dd*j)/B**2*K**2 +4*cc*t/B**2*V*K**2-8*cc/B*V*K*(NN+q*Rh) + 4*n*Du*Dv*(z*z*r*r*V/B**2*W**2-2*z*r*v/B*W*xi) +4*alphastar*Du*(v*eta-ze)**2+4*alpha*Du*((1-H)**2*U*V/B**2*W**2-2*(1-H)*u*v/B*W*(v*eta-ze)) +4*m*Mm**2+4*A*Du*Dv*YY**2+2*(JJ**2+Du*(v*YY-ZZ)**2+Dv*(XX+u*ZZ)**2-Du*Dv*ZZ**2)
print('LII6 residual parent readback',flush=True)
assert s.cancel(G-Gr)==0
print('PASS: original full parent = exact residual + four positive packets')
