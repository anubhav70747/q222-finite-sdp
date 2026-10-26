import sympy as S
s,b,h,t=S.symbols('s b h t');l=1-h;x=s*b;B=s*s*b;a=1/s;R=1+h*x;K=1+x;w=R;r=1-s*w;T=1+x*r;Sc=1+B*(l+a*r)
Del=S.cancel(s*(a+K+K*B+2*x*K*r+a**3*B**2*r**2+2*a**2*B*r+t*B*(1+a+x+a*a*B*r*(2-r))))
Mstar=K*(s+r)*(h+t)-s*w*(1+t*B)*l
C=S.cancel(T*Mstar/(2*l*Del));Ah=S.cancel(t*T*(s+r)*Sc/(2*Del));A=S.cancel(Ah/h)
betaD=S.cancel(K*(s+r)*(1-t*x*r)*Sc/(2*t*T*Del));c=S.cancel(l/(t*T));q=1-c
alpha=S.cancel(s**3*b*b*(s+r)*Mstar/(2*(1-h*B)*Del));H=alpha+S.Rational(1,2);E=t*K/R
Q=1+B*(1+t+a*h+2*a*r+a*t)+B*B*(a*(h+t)+a*a*(h*(r+1)+t*(1+r-r*r)))
Kstar=(1+s)*(h+t)*K;num=S.cancel(w*Q-b*Kstar)
n1=S.cancel(s*num*K/(2*Del*R));n2=S.cancel(num*r*K*K/(2*Del*R*R));mu1=S.cancel((1+t*B)*w*T*Sc/(2*Del*R))
# Use dimensionless RI fBB at v=0, u=0, and u^2 terms
BB0=C*c*c+betaD*q+q*q/2+A*c*c+2*H*c-S.Rational(1,2)
BC0=(-2*A*c-2*H+2*n1)/E
CC0=A/E**2
BBu=2*C*c*c+2*(betaD*q+q*q/2)+A*c*c+2*H*c-S.Rational(1,2)
BCu=BC0;CCu=CC0
# epsilon^2 = bt(1+2l*k^2)+... ; s=k^2/b
k2=s*b
eps2=S.cancel(b*t*K*K/(R*R))
base=b/(2*h); qb=b*(h+t)
models=[('BB0',BB0,base, -(1+b*h*(h+t))/(2*h)),('BC0',BC0,-2*base,(1-b*h*(h+t))/h),('CC0',CC0,base,-1/(2*h))]
for nm,ex,c0,kcoeff in models:
 actual0=S.factor((eps2*ex).subs(s,0))
 actuals=S.factor(S.diff(eps2*ex,s).subs(s,0))
 print(nm,'res difference',S.factor(actual0-c0),'k² coefficient diff',S.factor(actuals/b-kcoeff))
 assert S.factor(actual0-c0)==0
 assert S.factor(actuals/b-kcoeff)==0
for nm,ex,expect in [('BBu',BBu,base+qb/2),('BCu',BCu,-2*base),('CCu',CCu,base)]:
 print(nm,'diff',S.factor((eps2*ex).subs(s,0)-expect))
 assert S.factor((eps2*ex).subs(s,0)-expect)==0
print('RB17 9 coefficient equations EXACT')
