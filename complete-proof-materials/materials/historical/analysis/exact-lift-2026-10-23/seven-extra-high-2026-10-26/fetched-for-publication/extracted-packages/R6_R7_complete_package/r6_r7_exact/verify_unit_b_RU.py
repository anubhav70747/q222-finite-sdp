import sympy as S
s,b,t=S.symbols('s b t');h=g=S.Integer(0);l=1;x=s*b;B=s*s*b;a=1/s;R=1;K=1+x;ws=1;e=s;r=1-s;T=1+x*r
Sc=1+B*(l+a*r);DeltaStar=S.cancel(s*(a+K+K*B+2*x*K*r+a**3*B**2*r**2+2*a**2*B*r+t*B*(1+a+x+a*a*B*r*(2-r))))
Mstar=K*(s+r)*(h+t)-s*ws*(1+t*B)*l
C=S.cancel(T*Mstar/(2*l*DeltaStar));Ah=S.cancel(t*T*(s+r)*Sc/(2*DeltaStar));bd=S.cancel(K*(s+r)*(1-t*x*r)*Sc/(2*t*T*DeltaStar));c=1/(t*T);q=1-c;E=t*K;H=S.Rational(1,2)
Q=1+B*(1+t+a*h+2*a*r+a*t)+B*B*(a*(h+t)+a*a*(h*(r+1)+t*(1+r-r*r)))
Kstar=(1+s)*(h+t)*K;ng=S.cancel(ws*Q-b*Kstar)
n1=S.cancel(s*ng*K/(2*DeltaStar*R))
BB0=C*c*c+bd*q+q*q/2+2*H*c-S.Rational(1,2)
BC0=(-2*Ah*c-2*H+2*n1)/E
CC0=2*Ah/E**2
BBu2=2*C*c*c+2*(bd*q+q*q/2)-2*Ah*c*c+2*H*c-S.Rational(1,2)
BCu2=(-2*H+2*n1)/E
CCu2=2*Ah/E**2
pairs=[('BB0',BB0,1/t, -(1/t+(1/t+S.Rational(1,2))*b)),('BC0',BC0,-2/t, 2/t+(3/t-1)*b),('CC0',CC0,1/t,-1/t-2*b/t),('BBu2',BBu2,S.Rational(1,2),None),('BCu2',BCu2,-1/t,None),('CCu2',CCu2,1/t,None)]
for name,expr,v0,d0 in pairs:
    a0=S.factor(expr.subs(s,0))
    print(name,'res',a0,'diff expected',S.factor(a0-v0))
    assert S.factor(a0-v0)==0
    if d0 is not None:
        dd=S.factor(S.diff(expr,s).subs(s,0))
        print('first s',dd,'diff expected',S.factor(dd-d0))
        assert S.factor(dd-d0)==0
print('unit-b R7 9 tests done')
