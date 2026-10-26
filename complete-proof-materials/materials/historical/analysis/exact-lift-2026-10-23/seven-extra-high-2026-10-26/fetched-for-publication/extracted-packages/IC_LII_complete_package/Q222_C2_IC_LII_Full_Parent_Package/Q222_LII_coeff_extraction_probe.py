import sympy as s
u,v,W,xi,eta,ze,z,r,B,k,t,C,D,n,m,A,alpha,ellstar,alphastar=s.symbols('u v W xi eta ze z r B k t C D n m A alpha ellstar alphastar')
y=k*k;U=u*u;V=v*v;Du=1+U;Dv=1+V
H=B-y*t;j=1+z*r-t;cc=C*y*t;dd=D*y*j;ss=(1+z*r)/B
ell=y*(C+D);q=D/(C+D);g=y*C*D/(C+D)
XX=-z*r*v/B*W+xi;YY=H*u/B*W+eta;ZZ=u*v/B*W+ze
Mm=v*W-XX+u*v*YY-u*ZZ; JJ=u*W+u*v*XX-YY-v*ZZ;QQ=u*W-YY
F=4*n*Du*Dv*XX**2+4*m*Mm**2+4*C*Dv*QQ**2+4*A*Du*Dv*YY**2+4*D*JJ**2+4*alpha*Du*(v*YY-ZZ)**2+2*(JJ**2+Du*(v*YY-ZZ)**2+Dv*(XX+u*ZZ)**2-Du*Dv*ZZ**2)
N,rh=s.symbols('N rh')
eta_form=k*(N+q*rh)
XX=-z*r*v/B*W+xi;YY=H*u/B*W+eta_form;ZZ=u*v/B*W+ze
Mm=v*W-XX+u*v*YY-u*ZZ
K=u*k*W
JJ=k*((1-q)*rh-N-j/B*K)
Qn=u*xi-ze
Gr=(4*ellstar*N**2+4*g*rh**2+8*(dd-cc)/B*K*N-8*g*ss*K*rh+4*(cc*t+dd*j)/B**2*K**2
 +4*cc*t/B**2*V*K**2-8*cc/B*V*K*(N+q*rh)
 +4*n*Du*Dv*(z*z*r*r*V/B**2*W**2-2*z*r*v/B*W*xi)
 +4*alphastar*Du*(v*eta_form-ze)**2
 +4*alpha*Du*((1-H)**2*U*V/B**2*W**2-2*(1-H)*u*v/B*W*(v*eta_form-ze))
 +4*m*Mm**2+4*A*Du*Dv*YY**2+2*(JJ**2+Du*(v*YY-ZZ)**2+Dv*(XX+u*ZZ)**2-Du*Dv*ZZ**2))
G=s.expand(Gr)
print('G terms',len(s.Poly(G,[W,N,xi,ze,rh]).terms()),flush=True)
poly=s.Poly(G,[W,N,xi,ze,rh]); old=poly.coeff_monomial(W*rh)
rho_bounded=(-8*g*ss-8*cc*q*V/B-8*alpha*(1-H)*q*Du*V/B+8*m*q*V*(B+z*r+U*(H-1))/B+8*A*H*q*Du*Dv/B-4*y*j*(1-q)/B+4*Du*V*(H-1)*q/B)
assert s.cancel(old/(u*k)-rho_bounded)==0
print('LII7 formal bounded gauge coefficient PASS',flush=True)
gauged=s.expand(G-old*W*rh + old/(u*k)*(u*v*W*(u*xi-ze)-ss*U*(V-y)*W**2))
poly=s.Poly(gauged,[W,N,xi,ze,rh])
print('gauge leftover rhW',s.factor(poly.coeff_monomial(W*rh)))
px=s.cancel(poly.coeff_monomial(W*xi)/v)
pz=s.cancel(poly.coeff_monomial(W*ze)/(u*v))
pn=s.cancel(poly.coeff_monomial(W*N)/(u*k))
phi=s.cancel(poly.coeff_monomial(W*W))
# Exact finite support and degree checks before native source substitution.
assert s.cancel(px.subs({u:0,z:0})+8*m)==0
assert s.cancel(px.subs(u,-u)-px)==0
assert s.cancel(px.subs(v,-v)-px)==0
assert s.cancel(pz.subs(u,-u)-pz)==0
assert s.cancel(pz.subs(v,-v)-pz)==0
assert s.cancel(pn.subs(u,-u)-pn)==0
assert s.cancel(pn.subs(v,-v)-pn)==0
assert s.degree(phi,u)<=4
assert s.cancel(phi.subs(u,-u)-phi)==0
ds=4/B**2*(m*(B+z*r)**2+(n+s.Rational(1,2))*z*z*r*r)
es=4/B**2*(n+s.Rational(1,2))*z*z*r*r
assert s.cancel(phi.subs(u,0)-V*(ds+es*V))==0
print('PASS LII8: bounded-gauge formula, coefficient parities, Pxi ideal support, deg_U(phi)<=2 and exact phi(0,V)')
