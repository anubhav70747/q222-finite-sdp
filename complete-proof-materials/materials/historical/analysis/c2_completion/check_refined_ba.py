"""Independent exact verification of the new refined BA parent and joins."""
import json,time
from pathlib import Path
import sympy as S
start=time.monotonic(); out={}
def check(name,expr):
    ok=S.expand(S.together(expr).as_numer_denom()[0])==0
    out[name]=bool(ok)
    if not ok: raise AssertionError(name)
a,B,t,h,r=S.symbols('a B t h r'); u,v,W,X,Y,Z,p,qv=S.symbols('u v W X Y Z p qv')
n,m,C,A,D,alpha=S.symbols('n m C A D alpha')
e=1-r; l=1-h; x=a*B; y=t*B
ca=1+x*r; Sc=1+B*(l+a*r); g=1-a*e+h*x; j=ca-l/t
Delta=a+(1+x)+(1+x)*B+2*x*(1+x)*r+a**3*B**2*r**2+2*a*a*B*r+t*B*(1+a+x+a*a*B*r*(2-r))
Mass=(1+x)*(1+a*r)*(h+t)-a*(1+t*B)*e*l
weights={n:((1+a*r)*ca*Sc-r*Delta)/(2*r*Delta),m:a*a*(1+t*B)*e*ca*Sc/(2*g*Delta),C:ca*Mass/(2*l*Delta),A:t*ca*(1+a*r)*Sc/(2*h*Delta),D:(1+x)*(1+a*r)*(1-t*x*r)*Sc/(2*t*j*Delta),alpha:B*B*(1+a*r)*Mass/(2*(1-h*B)*Delta)}
Du=1+u*u; Dv=1+v*v; U=u*u; V=v*v
P=m*(v*W-X+u*v*Y-u*Z)**2+n*Du*Dv*X*X+C*Dv*(u*W-Y)**2+A*Du*Dv*Y*Y+(D+S.Rational(1,2))*(u*W+u*v*X-Y-v*Z)**2+(alpha+S.Rational(1,2))*Du*(v*Y-Z)**2+Dv*(X+u*Z)**2/2-Du*Dv*Z*Z/2
q=n*r*r+(m+S.Rational(1,2))*e*e+alpha; mu=m-a*a*q; lam=mu/(t*a*a*B)
E=t*u*W-v*Z; eta=Y-h*u*W; xi=X+r*u*Z; Fv=(V-y)*Z; RD=v*(U-x)*Z
Ystar=h*E/t+eta; LQ=a*v*W-u*Z-a*xi
Qsq=q*LQ**2
Nsq=n*(r*((a/t)*V-U)*Z+u*xi)**2
circuit=lam*(t*t*x*(U+(a/t)*V)*W*W-4*t*x*u*v*W*Z+(U*V+x*y)*Z*Z)
b2=2*m*h+C*l*l+(A+alpha+S.Rational(1,2))*h*h
b4=h*h*(m+A+alpha+S.Rational(1,2)); d4=(n+D+S.Rational(1,2))*r*r
e1=m*h*e+r*(D+S.Rational(1,2))*l+(alpha+S.Rational(1,2))*h
cF=(b2-a*a*n*r*r)/(t*t)
aE=2*A*h*h/(t*t)-2*e1/t; bE=2*b2/(t*t); cE=2*b4/(t*t)
aY=2*(A*h/t-m*e+r*(D+S.Rational(1,2))-alpha-S.Rational(1,2))
bY=2*(-C*l+A*h+m+(alpha+S.Rational(1,2))*h)/t
cY=2*(A*h+h*(m+alpha+S.Rational(1,2)))/t
x3=2*((D+S.Rational(1,2))*(l/t-1)+e/2-n*r*(1+a/t)-m*h/t)
x5=-2*r*(n+D+S.Rational(1,2))
H0=C*Dv*(l*E/t-eta)**2+A*Du*Dv*Ystar**2+(D+S.Rational(1,2))*(l*E/t-eta+u*v*xi)**2+(alpha+S.Rational(1,2))*Du*V*Ystar**2+(n*(1+V+U*V)+Dv/2)*xi**2+2*m*(V*E/t-u*v*xi)*Ystar+m*U*V*Ystar**2
rhs=Qsq+Nsq+circuit+d4*RD**2+(cF+b4*x/(t*t))*Fv**2+b4*v*RD*Fv/(t*t)+H0-lam*x*E**2-mu*xi**2-2*t*lam*x*LQ*xi+E*((aE+cE*V)*RD+(bE+cE*x)*v*Fv)+eta*((aY+cY*V)*RD+(bY+cY*x)*v*Fv)+xi*(x5*u*v*RD+u*(x3+x5*x)*Fv)
poly=S.Poly(S.expand(P-rhs),W,X,Y,Z,u,v)
coefficients=set(poly.coeffs())
for i,coefficient in enumerate(coefficients): check('BA7 coefficient '+str(i),coefficient.subs(weights,simultaneous=True))
out['BA7_number_of_distinct_coefficients']=len(coefficients)

# Opening regularity is checked on canceled rational functions with a free.
lam_source=S.cancel(lam.subs(weights,simultaneous=True))
num,den=lam_source.as_numer_denom()
corner={B:0,h:0,r:1}
out['lambda_canceled_denominator_corner']=str(S.factor(den.subs(corner)))
check('BA3 lambda joint residue',lam_source.subs(corner)-S.Rational(1,2))
if S.factor(den.subs(corner).subs(a,0))==0: raise AssertionError('lambda still has an a pole')
for name,expr,target in [('C',weights[C],t/2),('Ah',h*weights[A],t/2),('D',weights[D],1/(2*(t-1))),('alpha/B2',weights[alpha]/B**2,t*(1+a)/2),('n',weights[n],0),('m',weights[m],0)]: check('BA3 '+name,S.cancel(expr).subs(corner)-target)
Qs=1+B*(1+t+a*h+2*a*r+a*t)+B*B*(a*(h+t)+a*a*(h*(r+1)+t*(1+r-r*r)))
Ks=(1+a)*(h+t)*(1+a*B)
check('BA9 original physical numerator',2*r*Delta*weights[n]-e*Qs+B*Ks)
sigma_bound=t*a*a/(2*(1+a)*(1+t*(1+a)))
check('BA10 source Schur expression',S.Rational(1,2)-sigma_bound-(1+a+t*(1+2*a))/(2*(1+a)*(1+t*(1+a))))

# Complete homogeneous remote identities: no stationary/source specialization.
Dp=1+p*p; K=p*X+Z; VX=v*X; JI=W+Dp*VX-p*Y-p*v*K
rhs11=Dp*Dv*n*X*X+m*(p*v*W+v*Y-K)**2+C*Dv*(W-p*Y)**2+Dp*Dv*A*Y*Y+D*JI**2+Dp*alpha*(p*X+v*Y-K)**2+(JI**2+(1-p*p*v*v)*K*K+Dp*(v*v*Y*Y-2*v*Y*K+2*p*Y*VX+2*p*v*K*VX-p*p*VX*VX))/2
check('BA11 full homogeneous parent',p*p*P.subs(u,1/p)-rhs11)
PI=p*Z; QI=qv*Z
rhs15=n*(1+p*p)*(1+qv*qv)*X*X+m*(p*W-p*qv*X+Y-QI)**2+C*(1+qv*qv)*(W-p*Y)**2+A*(1+p*p)*(1+qv*qv)*Y*Y+D*(qv*W+X-p*qv*Y-PI)**2+alpha*(1+p*p)*(Y-QI)**2+((qv*W+X-p*qv*Y-PI)**2+(1+p*p)*(Y-QI)**2+(1+qv*qv)*(p*p*X*X+2*X*PI-PI*PI))/2
check('BA15 full double homogeneous parent',p*p*qv*qv*P.subs({u:1/p,v:1/qv})-rhs15)
tau=S.symbols('tau'); c0=tau/2; d0=1/(2*(tau-1)); VV=S.symbols('V')
check('BA12 determinant',(c0*(1+VV)+d0+S.Rational(1,2))*d0*VV-(d0+S.Rational(1,2))**2*VV-c0*d0*VV**2)
out['elapsed_seconds']=round(time.monotonic()-start,3)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({'all_equalities_passed':True,'distinct_BA7_coefficients':len(coefficients),'seconds':out['elapsed_seconds'],'lambda_denominator_corner':out['lambda_canceled_denominator_corner']},indent=2))
