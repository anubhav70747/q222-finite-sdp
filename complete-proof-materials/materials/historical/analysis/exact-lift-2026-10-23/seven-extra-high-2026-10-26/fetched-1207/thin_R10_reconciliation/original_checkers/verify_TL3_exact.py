#!/usr/bin/env python3
"""Exact source-only TL3 verification from native C2 BA1 and native C1 NC2/NC4.
No SDP, numerical tests, external source files or hidden assumptions.
Run: python verify_TL3_exact.py
"""
import sympy as S

a, B, r, h, t = S.symbols('a B r h t', nonzero=True)
s = 1/a
e = 1-r
ell = 1-h
x = a*B
y = t*B
K = 1+x
ca = 1+x*r
Sc = 1+B*ell+x*r
Qc = 1-h*B
R = 1+x*h
E = 1-t*x*r
g = R-a*e
jR = ca-ell/t

# Native C2: exact BA1, with positive denominator Delta and numerator M.
Delta = (a+(1+x)+(1+x)*B+2*x*(1+x)*r+a**3*B**2*r**2
         +2*a**2*B*r+t*B*(1+a+x+a**2*B*r*(2-r)))
M = K*(1+a*r)*(h+t)-a*(1+y)*e*ell
N = (1+a*r)*ca*Sc-r*Delta
n = N/(2*r*Delta)
c2 = {'m': a**2*(1+y)*e*ca*Sc/(2*g*Delta),
      'C': ca*M/(2*ell*Delta),
      'Ah': t*ca*(1+a*r)*Sc/(2*Delta),
      'D': K*(1+a*r)*E*Sc/(2*t*jR*Delta),
      'alpha': B**2*(1+a*r)*M/(2*Qc*Delta)}

# Formal same-contact C1: substitute OP4 into NC2/NC4, then d1=1.
# NC1 r1=h, theta1=e, h_C1=r, a1=ell/t, t1=a*e, eta=t,zeta=s.
Pd = 1+B*(1+t)+2*x*r+x*B*(h+t)+x**2*r*(h+t*e)
P4 = ca*E+x*ell*(r+s)
Palpha1 = R*r*(1+y)+s*g*E
c1 = {'m': K*x*(h+t)*r*(1+y)/(2*g*s*Pd),
      'C': K*(h+t)*E/(2*ell*Pd),
      'Ah': t*P4/(2*Pd),
      'D': K*(1+y)*(h+t)/(2*t*jR*Pd)-S.Rational(1,2),
      'alpha': x**2*s*(h+t)*Palpha1/(2*Qc*Pd),
      'beta': K*x*(h+t)*(r+s)/(2*e*Pd)-S.Rational(1,2)}

Pkap = ca**2+s*x*(1+t)+s*x**2*(r+h*e+t)
PD = t*x**2*r**2+x*(t-1+e*(ell-2*t))-1
PA = 1+s-2*e+x*(r**2+s*(r+h*e))
pi = K*B*(h+t)*(1+a*r)/(e*Pd)
kappa = (1+2*n)/pi
kappa_closed = e*ca*Sc*Pd/(r*Delta*K*B*(h+t))
Theta = n/B

checks=[]
def check(name,expr):
    assert S.cancel(expr)==0, f'FAILED: {name}: {S.factor(expr)}'
    checks.append(name)

# All master factorizations, including the printed N polynomial.
Npoly=e+B*(e-h-r*t+a*r*(2*e-h-t))-a*r*B**2*((h+t)+a*r*(h+t*(2-r)))
for name,expr in [('BA1 numerator in B',N-Npoly),
                  ('OP master',N-(e*Pd-K*B*(h+t)*(1+a*r))),
                  ('P4-rPd',N-(P4-r*Pd)),
                  ('C bridge',N-(e*Sc*E-r*B*M))]:check(name,expr)

# Small algebraic auxiliaries, not fitted to TL3.
for name,expr in [
 ('Pk',Pkap-(ca*Sc+K*B*(h+t))),
 ('PD',PD-(ca*t*jR-K*(h+t))),
 ('PA',PA-((r+s)*ca-e*Qc)),
 ('Sc composition',Sc-(Qc+x*(r+s))),
 ('E+tjR',E+t*jR-(h+t)),
 ('K-ca',K-ca-x*e),
 ('D aux 1',K*E+x*t*jR-(1+x*(h+t*e))),
 ('D aux 2',Pd-ca*(1+y)-B*(1+a*r)*(1+x*(h+t*e))),
 ('alpha aux',ca*Palpha1-E*PA-r*B*M),
 ('M aux',M-(ca*(a*h+t+a*t*r+a*t*B*h)+E*(h*K-a*e))),
]:check(name,expr)

# Original native C2 -> BA1 weak-source dictionary.  a_orig^2=t/a,
# z_orig=x, b_orig=1/a, t_orig=(1-h)/t; all coefficients are rational
# in these combinations, so no square-root choice is needed here.
aosq=t/a; zo=x; bo=1/a; to=ell/t
Delta_o=((1+zo)*zo*bo**2+2*zo*(1+zo)*bo*r
         +(1+zo)*(1+aosq*zo)*bo+zo**2*(1-aosq)*r**2
         +2*zo*(1+aosq*zo)*r+1+aosq*zo)
No=((1+zo)*(bo+r)*(bo+aosq*(1-to))
    -aosq*(1+aosq*zo)*(1-r)*to)
check('native old Delta to BA1',a*Delta_o-Delta)
check('native old positive N to BA1 M',a**2*No-M)
check('native old contact Sc',ca+aosq*zo*to-Sc)
check('native old gap g',((1+zo)*bo-(e+aosq*zo*to))-g/a)
check('native old gap h',bo-aosq*to-h/a)
check('native old gap j',ca-to-jR)
check('native old gap positive E',bo-aosq*zo*r-E/a)
check('native old gap positive Qc',1-zo*(bo-aosq*to)-Qc)
old2={
 'm': (1+aosq*zo)*e*ca*(ca+aosq*zo*to)/(2*((1+zo)*bo-(e+aosq*zo*to))*Delta_o),
 'C': ca*No/(2*aosq*to*Delta_o),
 'Ah': h*aosq*ca*(bo+r)*(ca+aosq*zo*to)/(2*(bo-aosq*to)*Delta_o),
 'D': (1+zo)*(bo+r)*(bo-aosq*zo*r)*(ca+aosq*zo*to)/(2*aosq*jR*Delta_o),
 'alpha':zo**2*(bo+r)*No/(2*(1-zo*(bo-aosq*to))*Delta_o)
}
oldn=((bo+r)*ca*(ca+aosq*zo*to)-r*Delta_o)/(2*r*Delta_o)
check('native old n to BA1',oldn-n)
for key in old2:check('native old '+key+' to BA1',old2[key]-c2[key])
# Original native C1 NC1, NC2 under OP4.
C1_p0=1+2*x*r+h*x**2*(r+s)+x*s
C1_p1=x*(K*s+x*r*e)
check('NC2 Pd to compact',C1_p0+t*C1_p1-Pd)
check('NC2 P4 to compact',ca*E+x*(1-h)*(r+s)-P4)
check('NC2 Palpha to compact',R*r*(1+y)+s*g*E-Palpha1)

# Parent half-paired rows and scale.
check('pi formal beta',c1['beta']+S.Rational(1,2)-pi/2)
check('kappa exact',kappa-kappa_closed)
check('kappa-1 factored',kappa-1-N*Pkap/(r*Delta*K*B*(h+t)))
check('native m exact',c2['m']-kappa_closed*c1['m'])
check('C',c2['C']-kappa_closed*c1['C']+Theta*ca/ell)
check('centered Ah',c2['Ah']-kappa_closed*c1['Ah']+Theta*t*ca*Sc/(K*(h+t)))
check('D',c2['D']-kappa_closed*c1['D']-Theta*Sc*PD/(K*(h+t)*t*jR))
check('alpha',c2['alpha']-kappa_closed*c1['alpha']+n*x*PA/(K*Qc))
check('d',1-kappa_closed+2*Theta*Pkap/(K*(h+t)))
check('beta+d/2',-kappa_closed*c1['beta']+(1-kappa_closed)/2+n)
check('D paired simpler',c2['D']+S.Rational(1,2)-kappa_closed*(c1['D']+S.Rational(1,2))+Theta*(B+Sc/(t*jR)))
check('alpha raw factor',e*ca*Sc*Palpha1-r*K*B*(r+s)*M-N*PA)

# Recover R10.2 identically from BA1, with no prior b bound.
sR,bR,wR,hR,tR=S.symbols('s_R b_R w_R h_R t_R',nonzero=True)
aR=1/sR;BR=sR**2*bR;xR=sR*bR;eR=sR*wR;rR=1-eR
subs={a:aR,B:BR,r:rR,h:hR,t:tR}
NR=S.factor(N.subs(subs, simultaneous=True))
Qstar=(1+BR*(1+tR+aR*hR+2*aR*rR+aR*tR)
      +BR**2*(aR*(hR+tR)+aR**2*(hR*(rR+1)+tR*(1+rR-rR**2))))
Kstar=(1+sR)*(hR+tR)*(1+xR)
check('R10.2 native reduction',NR-sR*(wR*Qstar-bR*Kstar))
Qstar_bound=(1+sR*xR*(1+tR)+xR*(hR+2*rR+tR)+sR*xR**2*(hR+tR)
             +xR**2*(hR*(rR+1)+tR*(1+rR-rR**2)))
check('Qstar no b inverse',Qstar-Qstar_bound)
DeltaStarR=S.cancel(sR*Delta.subs(subs, simultaneous=True))
DeltaStarBound=(1+sR*(1+xR)+sR**2*xR*(1+xR)
    +2*sR*xR*(1+xR)*rR+xR**2*rR**2+2*xR*rR
    +tR*(sR**2*xR*(1+xR)+sR*xR+sR*xR**2*rR*(2-rR)))
check('DeltaStar no b inverse',DeltaStarR-DeltaStarBound)

print('PASS:',len(checks),'exact symbolic identities')
for ix,c in enumerate(checks,1): print(f'{ix:02d}. {c}')
