"""Exact focused R2 source and complete normal identities. No numerical solver."""
import sympy as sp
S,z,q,w,h=sp.symbols('S z q w h', nonzero=True)
e=z*w;r=1-e;l=1-h
A=1/z;t=1/S;B=S*z*q;x=A*B;y=t*B
ca=1+x*r; sc=1+B*(l+A*r); R=1+h*x; g=R-A*e
jD=ca-l/t
Delta=A+(1+x)+(1+x)*B+2*x*(1+x)*r+A**3*B**2*r**2+2*A**2*B*r+t*B*(1+A+x+A**2*B*r*(2-r))
MM=(1+x)*(1+A*r)*(h+t)-A*(1+y)*e*l
n=((1+A*r)*ca*sc-r*Delta)/(2*r*Delta)
m=A**2*(1+y)*e*ca*sc/(2*g*Delta)
C=ca*MM/(2*l*Delta)
AA=t*ca*(1+A*r)*sc/(2*h*Delta)
D=(1+x)*(1+A*r)*(1-t*x*r)*sc/(2*t*jD*Delta)
al=B**2*(1+A*r)*MM/(2*(1-h*B)*Delta)
K=1+x; TC=1+x*r; GC=r+z; EC=1-q*r; QC=1-S*q*z*h
hC=r; η=1/S;ζ=z
P4=TC*EC+x*(1-h)*GC
Palpha=R*hC*(1+y)+ζ*g*EC
Pd=1+2*x*hC+h*x**2*GC+x*ζ+η*x*(K*ζ+x*hC*e)
B5=K*(1+y)*(h+η)-(η*TC-1+h)*Pd
Bβ=K*x*(h+η)*GC-e*Pd
c1=K*(h+η)*EC/(2*(1-h)*Pd)
a1=η*P4/(2*h*Pd)
m1=K*x*(h+η)*hC*(1+y)/(2*g*ζ*Pd)
d1=B5/(2*jD*η*Pd)
al1=x**2*ζ*(h+η)*Palpha/(2*QC*Pd)
β1=Bβ/(2*e*Pd)
P=1+q*z+S*q*(2+z-2*e)+S*q**2*(z+e-e**2)+S**2*h*q**2*(1+z-e)
P_C=1+q*z+S*(1+q*(2+2*z-e)+q**2*z)+S**2*(q*(h*(1+z)+1-e)+q**2*(h*e*z-e*(1+z)+2*z+1))+S**3*h*q**2*(1+z-e)
P_A=1+q*z+S*(q*(2+z-2*e)+q**2*z*(1-e))+S**2*q**2*(1+z-e)*(1+z-e-h*z)
P_D=1+q*(z+e)+S*(q*(h*(1+z)+1-e)+q**2*(z+e-e**2))+S**2*h*q**2*(1+z-e)
P_m=1+q*z+S*(q*(2+z-2*e)+q**2*z)+S**2*q**2*((1-e)**2+z*(1-e)+h*e*z)
P_b=(1+z)*(1+q*z)+S*q*(1+z)*(2+z-2*e)+S*q**2*z*(1+z-e**2)+S**2*q**2*(1+z-e)**2
lam=n/(q*z)
checks={
'TU6 Pd equals P': (Pd,P),
'TU9 C': ((1-h)*S*(C-c1),-n*(1-e)*P_C/(z*P)),
'TU9 A': (h*S*(AA-a1),-n*(1-e)*P_A/(z*P)),
'TU9 D': ((D-d1)/S,-n*(1-e)*(1+q*z)*P_D/(z*jD*P)),
'TU9 m': (z*g*(m-m1),n*(1-e)*(1+q*z)*P_m/(z*P)),
'TU9 hb': (z*w*(sp.Rational(1,2)-(β1+sp.Rational(1,2))),n*(1-e)*P_b/(z*P)),
'TU9 alpha': ((al-al1)/(S*z),n*(1-e)*q**2*(1+q*z)*((1+S)*e-S*(1-h)*(1+z))/(z*QC*P)),
'TU10 poly':((P_b-(1+q*z)*P_m).subs(z,0),sp.Integer(0)),
}
for key,(lhs,rhs) in checks.items():
    xexpr=sp.together(lhs-rhs, deep=False)
    num,den=xexpr.as_numer_denom()
    ok=(sp.Poly(num,S,z,q,w,h).is_zero if num!=0 else True)
    if not ok:
        # try cancel first; can be large
        num=sp.cancel(xexpr).as_numer_denom()[0]
        ok=sp.Poly(num,S,z,q,w,h).is_zero
    print(key,ok,flush=True)
    if not ok:
        print('FAILED NUMERATOR: ', sp.factor(num))
        raise SystemExit(1)
print('TU7 contact normal dictionary independent above; all TU source polynomial rows PASS')
