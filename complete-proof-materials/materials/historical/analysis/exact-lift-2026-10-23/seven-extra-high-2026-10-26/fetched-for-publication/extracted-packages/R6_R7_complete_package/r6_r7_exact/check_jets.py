import sympy as S
s,b,h,t,g,lam=S.symbols('s b h t g lam')
l=1-h;x=s*b;B=s*s*b;a=1/s;R=1+h*x;K=1+x;w=R-g;e=s*w;r=1-e;T=1+x*r
Sc=1+B*(l+a*r);j=T-l/t
DeltaStar=S.cancel(s*(a+K+K*B+2*x*K*r+a**3*B**2*r**2+2*a**2*B*r+t*B*(1+a+x+a*a*B*r*(2-r))))
Mstar=K*(s+r)*(h+t)-s*w*(1+t*B)*l
C=S.cancel(T*Mstar/(2*l*DeltaStar));Ah=S.cancel(t*T*(s+r)*Sc/(2*DeltaStar));betaD=S.cancel(K*(s+r)*(1-t*x*r)*Sc/(2*t*T*DeltaStar))
Q=1+B*(1+t+a*h+2*a*r+a*t)+B*B*(a*(h+t)+a*a*(h*(r+1)+t*(1+r-r*r)))
Kstar=(1+s)*(h+t)*K
num=S.cancel(w*Q-b*Kstar)
n1=S.cancel(s*num*K/(2*DeltaStar*R))

# obtain exact polynomial rational Taylor in lam up to degree 2 by truncated polynomials. no Series large.
def jet(expr):
    xx=S.cancel(expr.subs(g,0));num,den=S.fraction(xx)
    num=S.Poly(S.expand(num.subs({s:lam*s,b:lam*b})),lam)
    den=S.Poly(S.expand(den.subs({s:lam*s,b:lam*b})),lam)
    nv=[num.nth(i) for i in range(3)];dv=[den.nth(i) for i in range(3)]
    out=[]
    for k in range(3):out.append(S.cancel((nv[k]-sum(dv[j]*out[k-j] for j in range(1,k+1)))/dv[0]))
    return [S.factor(x) for x in out]
for name,ex in [('C',C),('Ah',Ah),('betaD',betaD),('n1',n1)]:
    print(name,'coeffs:',jet(ex))
