#!/usr/bin/env python3
"""Native strict C1 exact 16-feature infeasibility certificate, no SDP and no sampling.

All symbolic computations are rational. The two necessary Gram inequalities
come from an exact 4-contact Fourier argument and a forced 4x4 principal block.
"""
import sys
from pathlib import Path
import sympy as S
BASE=Path(__file__).resolve().parent/'original_baseline'
if not (BASE/'c1_native.py').exists():
    BASE=Path(__file__).resolve().parents[1]/'Q222_C1_fixed_admission_2026_10_26'
sys.path.insert(0,str(BASE))
from c1_native import k,f,g,a,t,Lc,Hc,Qc,R,T,M,J,P4,Pa,Pd,c2,cb,c3,c5,c4,cal,cd,D,be,NATIVE,CONTACT
from c1_small_map import u,v,W,X,Y,Z,FEATS,sym_parent

vals={k:S.Integer(4),f:S.Rational(1,1000),g:S.Rational(1,8),a:S.Integer(1),t:S.Integer(1)}
checks=[k,f,g,a,t,Lc,Hc,Qc,R,T,M,J,1-t*g*g,1-k*k*f*f+t*k*k*f*f*g*g,
        P4,Pa,Pd,c2,cb,c3,c5,c4,cal,cd,D,be]
for item in checks:
    z=S.cancel(item.subs(vals))
    assert z>0,('strict native condition',item,z)
source=tuple(S.cancel(t_.subs(vals)) for t_ in NATIVE)
contact=tuple(S.cancel(t_.subs(vals)) for t_ in CONTACT)
d,la,nn,mm,C,A,De,alpha,beta=source
uc,vc,wc,xc,yc,zc=contact
assert la==nn==0 and all(q>0 for q in (d,mm,C,A,De,alpha,beta))
assert all(q!=0 for q in contact)
assert xc<0 and all(q>0 for q in (uc,vc,wc,yc,zc))
print('PASS: every original strict A1 native denominator/coefficient and contact sign.')

# Exact SG2 sign gauges hold as polynomial identities with arbitrary source weights.
from c1_small_map import u,v,W,X,Y,Z
DD,MM,CC,AA,EE,AL,BE=S.symbols('dd mm cc aa ee al be')
Fgeneric=sym_parent((DD,0,0,MM,CC,AA,EE,AL,BE))
assert S.Poly(S.expand(Fgeneric.subs({u:-u,Y:-Y,Z:-Z})-Fgeneric),u,v,W,X,Y,Z).is_zero
assert S.Poly(S.expand(Fgeneric.subs({v:-v,X:-X,Z:-Z})-Fgeneric),u,v,W,X,Y,Z).is_zero
print('PASS: independent original SG2 sign congruences for every C1 coefficient tuple.')

# Six exact stationary identities on the whole physical parent. No state restriction.
F=sym_parent(source)
stationary={u:uc,v:vc,W:wc,X:xc,Y:yc,Z:zc}
assert S.cancel(F.subs(stationary))==0
for x in (u,v,W,X,Y,Z):
    assert S.cancel(S.diff(F,x).subs(stationary))==0, ('stationarity',x)
print('PASS: full four-state physical contact and all six exact stationarity equations.')

# Fourier necessary lower bound on odd half-power coefficient gamma.
sc=S.factor((4*beta+2*d)*vc**2*xc**2)
sb=S.factor(4*A*uc**2*yc**2)
se=S.factor(4*beta*uc**2*vc**2*zc**2)
assert sb>S.Rational(1,200)**2
assert sc<S.Rational(1,4000)**2
assert se<S.Rational(1,4000)**2
assert sb>sc and sb>se
# From sqrt(gamma) >= sqrt(sb)-sqrt(sc)-sqrt(se)
gamma_lower=S.Rational(9,2000)**2
assert gamma_lower==S.Rational(81,4_000_000)
print('PASS: exact Fourier inequalities; every 16-Gram would force gamma >',gamma_lower)

# Forced positive H7[1,2] pivot, forced H7[0,4] coupling,
# and the resulting upper bound on gamma^2 (all contact factors included).
B=2*C+2*De+d
h00=2*B*uc**2*wc**2
h04=-2*B*uc*wc*yc
h44=2*(2*A+B)*yc**2
h11=4*mm*vc**2*wc**2
h12=-4*mm*vc*wc*xc
h22=2*(2*beta+d+2*mm)*xc**2
assert S.factor(h11*h22-h12*h12)==8*mm*(2*beta+d)*vc**2*wc**2*xc**2
assert S.factor(h00*h44-h04*h04)==8*A*B*uc**2*wc**2*yc**2
delta=h11*h22-h12*h12
u_gamma2=S.factor(4*delta/h22*(h00-h04*h04/h44))
expect=S.factor(64*A*mm*(2*beta+d)*B*uc**2*vc**2*wc**4/((2*beta+d+2*mm)*(2*A+B)))
assert S.cancel(u_gamma2-expect)==0
assert u_gamma2 < S.Rational(1,100_000)**2
assert gamma_lower>S.Rational(1,100_000)
print('PASS: forced positive Schur pivots; every 16-Gram would require 0 <= gamma < 1/100000.')
print('PASS: incompatible lower and upper bounds, exact strict C1 source outside 16-feature cone.')
# Stronger entirely linear, face-local infeasibility certificate: Cauchy gives
# gamma >= S_b/3 - S_c - S_e; a fixed rational H7 test vector gives
# Q_H7(gamma) >=0. Multiply the first by a positive rational and add to
# eliminate gamma, leaving a strictly negative rational coefficient functional.
gamma0=S.Rational(1,100_000)
lower_lin=S.factor(sb/3-sc-se)
assert lower_lin>S.Rational(1,100_000)
rat_a=S.factor(gamma0*h22/(2*delta))
rat_b=S.factor(-h12/h22*rat_a)
rat_c=S.factor(-h04/h44)
assert rat_a>0
# (H7[0,0],H7[0,1]=-gamma/2,H7[0,4],H7[{1,2}],H7[4,4]).
# Scalar PSD on q=(e0 + rat_a e1 + rat_b e2 + rat_c e4):
Q_at_gamma0=S.factor(h00+2*h04*rat_c+h44*rat_c**2+
           h11*rat_a**2+2*h12*rat_a*rat_b+h22*rat_b**2-gamma0*rat_a)
assert Q_at_gamma0<0
L_sep=S.factor(Q_at_gamma0+rat_a*(gamma0-lower_lin))
assert L_sep<0
assert rat_a==S.Rational(276058888004033,102400025600000)
assert rat_b==-S.Rational(3999937,1612800000)
assert rat_c==S.Rational(68000000000000,67999932000063)
print('PASS: fixed rational vector + contact Fourier/Cauchy give exact LINEAR separating functional.')
print('GAMMA_LINEAR_LOWER',lower_lin)
print('FIXED_H7_TEST_VECTOR', (S.Integer(1),rat_a,rat_b,rat_c))
print('SEPARATION_VALUE',L_sep)

print('C1_PARAMETERS', [vals[x] for x in (k,f,g,a,t)])
print('NATIVE_WEIGHTS',source)
print('CONTACT',contact)
print('SQUARE_LENGTHS',sc,sb,se)
print('GAMMA_UPPER_SQUARED',u_gamma2)

# The four sign-related contact data, and the character-average property of
# integer feature polynomials. No sign argument is applied to sqrt(uv).
ps=(W,X,Y,Z)
chars=[]
for label,iu,iv,idx in FEATS:
    xchar=(iu+(idx in (2,3)))%2
    ychar=(iv+(idx in (1,3)))%2
    chars.append((label,xchar,ychar))
trivial=[label for label,xx,yy in chars if xx==yy==0]
assert trivial==['vX','uY','uvZ'],trivial
for sign_u in (+1,-1):
 for sign_v in (+1,-1):
    xx={u:sign_u,v:sign_v,W:1,X:sign_v,Y:sign_u,Z:sign_u*sign_v}
    odd=sign_u*sign_v*(1-sign_u*sign_v)**2
    assert odd==(0 if sign_u*sign_v==1 else -4)
    assert all(S.cancel(S.expand(F.subs({u:uc*u,v:vc*v,W:wc*W,X:xc*X,Y:yc*Y,Z:zc*Z})).subs(xx))==0 for _ in [0])
print('PASS: four sign contacts, fractional square values 0,0,-4,-4 and trivial Fourier sector.')
