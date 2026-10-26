import sympy as S
s,b,h,t=S.symbols('s b h t');g=S.Integer(0);l=1-h;x=s*b;B=s*s*b;a=1/s;R=1+h*x;K=1+x;w=R;r=1-s*w;e=s*w;T=1+x*r
Sc=1+B*(l+a*r);j=T-l/t
DeltaStar=S.cancel(s*(a+K+K*B+2*x*K*r+a**3*B**2*r**2+2*a**2*B*r+t*B*(1+a+x+a*a*B*r*(2-r))))
Mstar=K*(s+r)*(h+t)-s*w*(1+t*B)*l
C=S.cancel(T*Mstar/(2*l*DeltaStar));Ah=S.cancel(t*T*(s+r)*Sc/(2*DeltaStar));betaD=S.cancel(K*(s+r)*(1-t*x*r)*Sc/(2*t*T*DeltaStar))
D=S.cancel(betaD*T/j)
mg=S.cancel((1+t*B)*w*T*Sc/(2*DeltaStar));alpha0=S.cancel((s+r)*Mstar/(2*(1-h*B)*DeltaStar))
Q=1+B*(1+t+a*h+2*a*r+a*t)+B*B*(a*(h+t)+a*a*(h*(r+1)+t*(1+r-r*r)))
Kstar=(1+s)*(h+t)*K;ns2=S.cancel((w*Q-b*Kstar)/(2*r*DeltaStar))
table=[('n/s²',ns2, (1-b*(h+t))/2),('mg',mg,S.Rational(1,2)),('C',C,(h+t)/(2*l)),('Ah',Ah,t/2),('D',D,1/(2*(t-l))),('alpha/s³b²',alpha0,(h+t)/2)]
for name,ex,target in table:
    v0=S.factor(S.cancel(ex.subs(s,0)));assert S.factor(v0-target)==0,(name,v0)
    print(name,'residue',v0)
for name,ex,target in [('mg',mg,(b*h-1)/2),('C',C,-(t+1)/(2*l)),('Ah',Ah,-t/2),('D',D,-(b*t*(h+t)+t-l)/(2*(t-l)**2))]:
    deriv=S.factor(S.diff(ex,s).subs(s,0))
    print('first s',name,deriv,'expected difference',S.factor(deriv-target))
    assert S.factor(deriv-target)==0
