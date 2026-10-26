"""Exact, nonsampled UL1 complete-parent interface checks.
These identities are field-polynomial and hold before independent source and test specialization.
They verify the full SC3 parent, the six normalized source-difference
coefficient channels (including two first normals), and UL8/9/15 root algebra.
They do NOT assert source-uniform strict residue by themselves.
"""
import sympy as s

def check(label, expr):
    ans=s.cancel(expr)
    if ans!=0:
        raise AssertionError(f'{label}: {s.factor(ans)}')
    print('PASS',label)

u,v,W,X,Y,Z,b,g,w,R,h2,hb,cb,ab,Db,alphab=s.symbols('u v W X Y Z b g w R h2 hb cb ab Db alphab',nonzero=True)
Du=1+u*u;Dv=1+v*v;S=b*b; z=g*g; p=w/R; M0=R-w
Qb=(W+u*Y)/b;Yb=Y/S
N=((X+u*Z)/b-p*v*Qb)/g
H=u*Qb/b;L=v*Qb/g;V=u*v*Z/b
Bb=g*N+p*v*Qb
Jb=S*(H-Du*Yb)+u*b*v*Bb-Du*v*Z
Kb=S*(H-Du*Yb)**2+2*b*u*v*(H-Du*Yb)*Bb+u*u*v*v*Bb**2+Du*v*v*S*Yb**2+Du*V**2-2*Du*V*(Qb+v*Bb)+2*Du*u*u*v*Z*Yb
G= 4/R*((h2+hb*Dv)*N*N+2*(hb*Dv-h2)*L*N+(h2*(1-p)+hb*p*Dv)*L*L)+4*Dv*(cb*(H-Du*Yb)**2+ab*Du*Yb**2)+4*Db*Jb**2+4*alphab*z*Du*(v*S*Yb-Z)**2+2*Kb
pair=4*h2*w/(M0*R)*N*N+4*hb*M0/(w*R)*Dv*N*N
M_ev=v*W-X+u*v*Y-u*Z;J_ev=u*W+u*v*X-Y-v*Z;B_ev=X+u*Z
# Beta-paired coefficient is 4hb/(zw) and alpha weight=alphab*S*z
FdivS=(4*h2/(z*M0)*M_ev**2+4*cb/S*Dv*(u*W-Y)**2+4*ab/S*Du*Dv*Y**2+(4*Db*S+2)*J_ev**2+(4*alphab*S*z+2)*Du*(v*Y-Z)**2+4*hb/(z*w)*Dv*B_ev**2-2*Du*Dv*Z**2)/S
check('SC3 entire 4-state parent, arbitrary normalized source weights', FdivS-G-pair)

# Source-difference readback from SC3. Hard variables Nn, Ny remain unrestricted.
Nn,Ny,HH,LL,sigma_rho,rho0=s.symbols('Nn Ny HH LL sigma_rho rho0',nonzero=True)
tc,ta,td,tm,tbeta,talpha,lam,c=s.symbols('tc ta td tm tbeta talpha lam c',nonzero=True)
# Differences in normalized weights are lam*c^2 times their indicated bounded quotients.
Lmass=lam*c*c
at=1-rho0+p*v*v
YY=(rho0*HH+Ny/sigma_rho)/Du
JJ=S*(HH-Du*YY)+u*b*v*(g*Nn+p*v*(b*HH/u))-Du*v
Delt=4/R*((Lmass*tm+Lmass*tbeta*Dv)*Nn**2+2*(Lmass*tbeta*Dv-Lmass*tm)*LL*Nn+(Lmass*tm*(1-p)+Lmass*tbeta*p*Dv)*LL**2)+4*Dv*(Lmass*tc*(HH-Du*YY)**2+Lmass*ta*Du*YY**2)+4*Lmass*td*JJ**2+4*Lmass*talpha*z*Du*(v*S*YY-1)**2
# UL6 definitions, as absolute normalized differences divided by lam*q.
ma=4*Dv*R**2*(tc*(1-rho0)**2+ta*rho0**2/Du)+4*talpha*z*S**2*v**2*R**2*rho0**2/Du+4*td*S**2*R**2*at**2
mb=-8*talpha*z*S*R*rho0-8*td*S*R*at*Du
mc=4*td*Du**2
ml=4*R*(tm*(1-p)+tbeta*p*Dv)
ht,br=s.symbols('ht br')
raw_scalar=Lmass*(ma*ht**2+mb*ht*v+mc*v**2+ml*br**2/(c*c)+4*talpha*z*Du)
check('Full scalar TU source difference; all six weights',Delt.subs({Nn:0,Ny:0,HH:R*ht,LL:R*br/c})-raw_scalar)
# QI scalar normal form, using Gamma^2=q ht v, A=c ht, E=ht-v, D=q g.
Astar=c*ht; Bstar=br; Dstar=c*c*g; Estar=ht-v; Gam2=c*c*ht*v
scalar_toric=lam*((ma-mc)*Astar**2+(mb+2*mc)*Gam2+ml*Bstar**2)+Lmass*mc*Estar**2+4*Lmass*talpha*Du/(c**4)*Dstar**2
check('Entire scalar QI6 correction, before actual current',raw_scalar-scalar_toric)
# SC3 first-normal derivatives exactly give UL10 and UL12, including every D/alpha contribution.
DL_N=s.diff(Delt,Nn).subs({Nn:0,Ny:0,HH:HH,LL:R*Bstar/c})
FN_base=(tbeta*Dv-tm)
check('SC3 source-only N correction', DL_N-(8*lam*c*FN_base*Bstar+8*Lmass*td*u*b*v*g*(S*at*HH-Du*v)))
# Add the complete TU17 *signed* first-normal remainder (-8 lambda D_uD_v (1-e)c z B_* N); the retained current square is separate.
e=s.symbols('e')
FN_full=FN_base/(1-e) - z*Du*Dv
DL_N_full=DL_N-8*lam*Du*Dv*(1-e)*c*z*Bstar
ULN=8*lam*c*(1-e)*FN_full*Bstar+8*Lmass*td*u*b*v*g*(S*at*HH-Du*v)
check('Entire first hard normal UL10 INCLUDING TU17 signed current', DL_N_full-ULN)
DL_y=s.diff(Delt,Ny).subs({Nn:0,Ny:0})
FH=Dv*(rho0*ta-(1-rho0)*tc)-Dv*ta*u*u*rho0/Du+talpha*z*S*S*v*v*rho0/Du-td*S*S*at
Fv=td*S*Du-talpha*z*S
ULy=8*Lmass/sigma_rho*(FH*HH+Fv*v)
check('Entire harmonic normal UL12 algebra',DL_y-ULy)

# Root identities UL8/UL9, every independent testing variable kept symbolic.
Qbb=s.symbols('Qbb',nonzero=True)
HHt=u*Qbb/(b*R)
AA=c*HHt; BB=c*v*Qbb/(g*R); CC=u*v/b; DD=c*c*g; UU=c*g*u/b
EE=HHt-v; FF=v*HHt-g*DD
Gam2=c*c*HHt*v
nr=lam*c*c*g*g
mone=nr*S/(z*z)
check('TU16 product U* B*=Gamma^2',UU*BB-Gam2)
check('UL8 U A^2',lam*u*u*AA**2-(lam*c*c)*(u*EE+b*CC)**2)
check('UL8 V A^2',lam*v*v*AA**2-(lam*c*c)*(FF+g*DD)**2)
check('UL8 U Gamma^2',lam*u*u*Gam2-(lam*c*c)*b*CC*(u*EE+b*CC))
check('UL8 V Gamma^2',lam*v*v*Gam2-(lam*c*c)*(FF+g*DD)*(FF+g*DD-v*EE))
check('UL8 U B^2',lam*u*u*BB**2-mone*(FF+g*DD)**2)
rootlam=s.symbols('rootlam',nonzero=True)
J3=rootlam*(UU-BB)
check('UL8 V B, normal including J3',rootlam*v*BB-(rootlam*c*g*CC-v*J3))
check('UL9 dangerous BC coefficient', (w*g/c)*v*BB*CC-(w*z*CC**2-(w*g/(c*rootlam))*v*CC*J3))

sqrtU,sqrtB=s.symbols('sqrtU sqrtB')
GM=sqrtU*sqrtB; JJ=rootlam*(sqrtU*sqrtU-sqrtB*sqrtB)
VV2=rootlam**2*GM*(sqrtU-sqrtB)**2
check('UL15 J3^2 split',JJ**2-4*VV2-rootlam**2*(sqrtU-sqrtB)**4)
check('UL15 Gamma J3 cross',GM*JJ-VV2/rootlam+2*rootlam*GM*(sqrtB**2-GM))
d,hard=s.symbols('d hard')
check('UL15 retained half packet and hard block split',4*(JJ-d*hard)**2-(2*(JJ-d*hard)**2+JJ**2+4*VV2+rootlam**2*(sqrtU-sqrtB)**4-4*d*JJ*hard+2*d*d*hard**2))
print('Every displayed identity passed exactly.')

# QI7/QI8: all positive square corrections at the contact, exact.
ri,rj=s.symbols('ri rj',nonzero=True)
# ri = sqrt(A_*), rj=sqrt(c v), Gamma=ri*rj, A*=ri^2, c v=rj^2.
Gamma=ri*rj; As=ri**2; cv=rj**2; Ds=c*c*g
EE1=(As-cv)/c
RR_E2=Gamma*(ri-rj)**2
check('QI7 c Gamma E*',c*Gamma*EE1-(2*Gamma*(As-Gamma)-RR_E2))
RR_F2=(g*Gamma/(c*c))*(Gamma-Ds)**2
FF1=(Gamma**2/(c*c)-g*Ds)
check('QI7 Gamma F*',Gamma*FF1-(-2*g*Gamma*(Ds-Gamma)+RR_F2/g))
check('QI7 full (v h)^2 split', (Gamma**2/(c*c))**2-(FF1**2+2*RR_F2+2*g*FF1*(Ds-Gamma)+g*g*Ds**2-4*g*g*Gamma*(Ds-Gamma)))
check('QI8 q E^2 roots', c*c*EE1**2-(4*RR_E2+(ri-rj)**4))

# SC16: exact strict positive toric average, using all 16 roots.
x1,x2,x3,x4=s.symbols('x1 x2 x3 x4')
xs=[x1,x2,x3,x4]; Ai=[x**4 for x in xs]; GM=x1*x2*x3*x4
T4=sum(ai**2 for ai in Ai)-4*GM**2
pairs=[(0,1,2,3),(0,2,1,3),(0,3,1,2)]
for idx,(i,j,k,l) in enumerate(pairs):
    alpha=[ai-GM for ai in Ai]
    Q=(xs[i]**2*xs[j]**2-xs[k]**2*xs[l]**2)
    D=(xs[i]**2-xs[j]**2)
    D2=(xs[k]**2-xs[l]**2)
    P=(xs[i]*xs[j]-xs[k]*xs[l])
    Tdiff=(alpha[i]-alpha[j])**2+(alpha[k]-alpha[l])**2+2*Q**2
    Tvar=sum(t**2 for t in alpha)+2*GM*D**2+2*GM*D2**2+4*GM*P**2
    check(f'SC16 difference identity pairing {idx}',T4-Tdiff)
    check(f'SC16 variance identity pairing {idx}',T4-Tvar)
print('All 24-root identities and complete-parent interfaces passed.')

# Actual native UL7/UL13 coefficient specialization (including the h and 1-h opening cancellations).
Ss,z0,q0,ww,hh,U0,V0=s.symbols('Ss z0 q0 ww hh U0 V0')
ee=z0*ww; RR=1+Ss*q0*hh; KK=1+Ss*q0
JJ=KK-Ss*q0*ee-Ss*(1-hh)
PP=1+q0*z0+Ss*q0*(2+z0-2*ee)+Ss*q0*q0*(z0+ee-ee**2)+Ss**2*hh*q0*q0*(1+z0-ee)
GG=1+z0-ee
PC=1+q0*z0+Ss*(1+q0*(2+2*z0-ee)+q0*q0*z0)+Ss**2*(q0*(hh*(1+z0)+1-ee)+q0*q0*(hh*ee*z0-ee*(1+z0)+2*z0+1))+Ss**3*hh*q0*q0*GG
PA=1+q0*z0+Ss*(q0*(2+z0-2*ee)+q0*q0*z0*(1-ee))+Ss**2*q0*q0*GG*(GG-hh*z0)
PD=1+q0*(z0+ee)+Ss*(q0*(hh*(1+z0)+1-ee)+q0*q0*(z0+ee-ee**2))+Ss**2*hh*q0*q0*GG
PM=1+q0*z0+Ss*(q0*(2+z0-2*ee)+q0*q0*z0)+Ss**2*q0*q0*((1-ee)**2+z0*(1-ee)+hh*ee*z0)
PB=(1+z0)*(1+q0*z0)+Ss*q0*(1+z0)*(2+z0-2*ee)+Ss*q0*q0*z0*(1+z0-ee**2)+Ss**2*q0*q0*GG**2
QQ=1-Ss*q0*z0*hh
TC=-(1-ee)*PC/((1-hh)*PP); TA=-(1-ee)*PA/(hh*PP)
TD=-(1-ee)*(1+q0*z0)*PD/(JJ*PP); TM=(1-ee)*(1+q0*z0)*PM/PP
TB=(1-ee)*PB/PP; TT=(1-ee)*q0*q0*(1+q0*z0)*((1+Ss)*ee-Ss*(1-hh)*(1+z0))/(QQ*PP)
Duu=1+U0; Dvv=1+V0; rrho=KK*hh/RR; a_tt=1-rrho+(ww/RR)*V0
MA=4*Dvv*RR**2*(TC*(1-rrho)**2+TA*rrho**2/Duu)+4*TT*z0*Ss**2*V0*RR**2*rrho**2/Duu+4*TD*Ss**2*RR**2*a_tt**2
MB=-8*TT*z0*Ss*RR*rrho-8*TD*Ss*RR*a_tt*Duu
MC=4*TD*Duu**2
ML=4*RR*(TM*(1-ww/RR)+TB*(ww/RR)*Dvv)
for nm,F in [('UL7 X',MA-MC),('UL7 Y',MB+2*MC+8*Duu*Dvv*(1-ee)),('UL7 Z',ML-4*Duu*Dvv*(1-ee**2))]:
    check(nm+' exact constant',F.subs({q0:0,U0:0,V0:0}))
FH0=Dvv*(rrho*TA-(1-rrho)*TC)-Dvv*TA*U0*rrho/Duu+TT*z0*Ss**2*V0*rrho/Duu-TD*Ss**2*a_tt
FV0=TD*Ss*Duu-TT*z0*Ss
check('UL13 F_H exact constant',FH0.subs({q0:0,U0:0,V0:0})-(1-ee)*Ss/(1-Ss*(1-hh)))
check('UL13 F_v exact constant',FV0.subs({q0:0,U0:0,V0:0})+(1-ee)*Ss/(1-Ss*(1-hh)))
check('TU10 exact q subtraction source polynomial',s.cancel((PB-(1+q0*z0)*PM).subs({ee:z0*ww,z0:0})))
check('UL10 F_N q=0 exact',s.cancel(((PB*Dvv-(1+q0*z0)*PM)/PP-z0*Duu*Dvv).subs(q0,0)-(V0-z0*U0*Dvv)))
print('Native UL coefficient specializations passed.')
