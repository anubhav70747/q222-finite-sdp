"""Exact independent BN/UN symbolic identities, pinned source 171a7e94; no samples.
Call python verify_bn_un_independent.py ."""
import sympy as s

def check(label, expr):
    val = s.cancel(s.together(expr)); assert val == 0, (label,s.factor(val)); print('PASS',label)
# Original source data: ES1-ES6 and FM2
k,nu,S,h,r=s.symbols('k nu S h r', real=True)
ell=1-nu; q=nu/ell; x=k*S; sig=nu*S; L=ell*S; th=1-h
R=1+x*r; J=1+x*h-sig*(1-r); M=R-k*th/ell
P0=1+S*(nu*ell+2*nu*k*h+k*ell+k**2*h*th)+S**2*(nu*r*k*ell+nu*r*k**2*h)
H0=1+L*(1-r); Bs=P0+k*(J-sig); K=1+x; e=k*h/nu
mathH=R*H0+(1-r)*x*h*(2+x*r)-e*(2+x*(1+h+r*th))+K*e**2
beta=(h*Bs-nu*H0)/(2*th*P0)
D=s.cancel((nu**2*S*mathH/(2*P0)-k*S*th*beta)/J)
d_s=s.cancel(D/(S*q**2)); mu=K*(1+sig*r)*h/(2*M*ell*nu*P0)
check('source BN beta/q', beta/q-(h/q*Bs-ell*H0)/(2*th*P0))
check('source sphere', J*d_s+k*th*beta/q**2 - ell**2*mathH/(2*P0))
check('source SN mu', mu*ell**2-(h/q)*K*(1+sig*r)/(2*M*P0))
# Physical sphere coordinates (no positivity computation, pure algebra).
b0sq=s.cancel(2*k*th*beta*P0/(nu**2*mathH))
d0sq=s.cancel(2*P0*J*D/(nu**2*S*mathH))
check('ES5 d0²+b0²=1',d0sq+b0sq-1)
check('ES7 source sphere',e**2*Bs-k*(e*H0+h*mathH*b0sq))
zn=s.cancel(k*(beta*th**2+2*D*h+h**2/2)/q**2)
znexpected=ell**2*mathH*th*b0sq/(2*P0)+ell**2*(e*H0+h*mathH*b0sq)/(2*Bs)+2*k*S*d_s*h
check('UN3 bounded zn identity',zn-znexpected)

# exact UN4 cancellation, T is polynomial in k,nu,S,h,r
U=K*(1+sig*r); Mstar=ell*M
T=s.expand(ell*J*U+J*Mstar*P0-J*Mstar*k*th*(J-sig)-2*Mstar*U)
Tnu=r-1-nu-nu*S*(1-r)*(r+ell)
Tkap=s.div(s.Poly(s.expand(T-ell*nu*S*Tnu),k), s.Poly(k,k))[0].as_expr()
check('T kappa monic remainder', T - k*Tkap-ell*nu*S*Tnu)
n=beta*th**2+2*D*h+h**2/2
A00=s.cancel(mu*ell**2-n/q)
A00expected=th*ell*H0/(2*P0)+(e*Tkap+h*ell*S*Tnu)/(2*J*M*P0)
check('UN4 bounded radial A00', A00-A00expected)
t0=(D+h/s.Integer(2)-beta*th)/q
t0expected=D/q + (ell*H0-k*(h/q)*(J-sig))/(2*P0)
check('UN15 X pole',t0-t0expected)
# Contact and conductor (formal). Independent testing/source symbols with only x=k*S, zeta=ell/k.
u,w,W,Z=s.symbols('u w W Z',real=True)
a=u**2;c=w**2;delta=c-q;zeta=ell/k
N=k*W-ell*u*w*Z;N0=N/ell
G=u*W/s.sqrt(S)-ell*s.sqrt(S)*w*Z
B=s.sqrt(S)*delta*Z;Bu=u*delta*Z
Rad=w*(a-x)*Z
check('correct N0 = N/ell', N0-(W/zeta-u*w*Z))
check('radial conductor SN7', Rad-(k*s.sqrt(S)*G-u*N)/ell)
check('cross product SN14 conductor',c*(a-x)*delta*Z**2-(k*G*(w*B)-N*(w*Bu))/ell)
check('SN8 normalized helper identity',delta**2*N0**2+2*k*delta/ell*G*(w*B)+2*k*(w*B)**2-(w*Bu)**2-(delta*W/zeta)**2)
# Normalize four circuit rows in rational source roots to avoid nested radicals.
t=s.symbols('t',positive=True) # represents sqrt(q)
j=s.symbols('j',positive=True) # represents sqrt(S)
N2=n # UN uses n>0; BN use beta. Generic circuit coefficient C>0
C=s.symbols('coef',positive=True)
# represent sqrt(coef), sqrt(k), sqrt(S), sqrt(q) via abstract positive roots and relation
z,kk=s.symbols('z kk',positive=True)
# put k=kk^2, S=j^2, q=t^2, ell remains a symbol for circuit because q = nu/ell not used
E1=z*(u*Z-W/( (ell/kk**2)*t));E2=z*(kk*j*w/t*Z-u*w*W/((ell/kk**2)*t**2*kk*j));E3=z*(u*Z-kk*j*w/t*Z)
# substitutions sqrt(C)=z, sqrt(k)=kk, sqrt(S)=j, sqrt(q)=t; general symbolic relation zeta=ell/kk^2
E1=s.cancel(E1.subs(zeta,ell/kk**2));E2=s.cancel(E2.subs(zeta,ell/kk**2));E3=s.cancel(E3.subs(zeta,ell/kk**2))
N0c=kk**2*W/ell-u*w*Z;Gc=u*W/j-ell*j*w*Z;Bc=j*(w**2-t**2)*Z;Buc=u*(w**2-t**2)*Z
bn_b=z**2/t**2;bn_z=kk**2*z**2/t**4
# BN14 helper w B: sqrt(zeta_BN) = sqrt(k*coef/q²) = kk*z/t²
check('BN14 wB helper',w*Bc+(w+t)/(kk*z/t**2)*E2+w*(w+t)/ell*Gc)
check('BN17 bBu helper',bn_b*Buc+(w+t)*(z/t)*E1+bn_b*(w+t)*N0c)
check('UN20 Bu helper',Buc+(w+t)*(t/z)*E1+(w+t)*N0c)
# UN21: source sqrt(k/q)=kk/t, sqrt(q/n)=t/z
check('UN21 N0 pivot',N0c-(t/z)*(E2-E1)+(w-t)/z*E3-(kk/t)/ell*w*Gc)

# UN16 complete hard-square exchange (all full-state values, no contact substitution).
Hsym=s.symbols('Hsym',real=True)
tau0,ee,dd=s.symbols('tau0 ee dd',real=True)
Nx=s.symbols('Nx',real=True)
QQ=w*N0c+Buc
Gh=-tau0*u*QQ+ee*w*j*Gc
Xe_new=Nx+tau0*QQ
check('UN16 exact square exchange', dd*(Gh**2+u**2*Nx**2+2*tau0*u**2*QQ*Nx)- (dd*u**2*Xe_new**2-2*dd*ee*tau0*u*w*j*QQ*Gc+dd*ee**2*j**2*w**2*Gc**2))

print('All source/circuit identities PASSED')
