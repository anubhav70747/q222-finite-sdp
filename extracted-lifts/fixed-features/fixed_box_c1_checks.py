"""Exact polynomial checks for fixed_box_c1.tex; no numerical fitting."""
import json
from fractions import Fraction
from pathlib import Path
import sympy as sp

out = {}
x, sig, z, b, h, r = sp.symbols('x sig z b h r')
K = 1+x
M = 1+x*r-x*z*(1-h)
P = sig+K+sig*z*(1+b*x+2*x*h)+sig*b*x*x*h*z*z+x*x*h*(1-h)*z
Nk = sp.expand(K*(1+sig*r)*(1+sig*z)-M*P)
Fx = sp.cancel((Nk-Nk.subs(x,0))/x)
Fs = sp.cancel(Nk.subs(x,0)/sig)
assert sp.expand(Fs-(r-1+sig*r*z)) == 0
fx_l1 = sum(abs(a) for a in sp.Poly(Fx,x,sig,z,b,h,r).coeffs())
assert fx_l1 == 37
out['source_cancellation'] = {'Fx_coefficient_l1': int(fx_l1),
    'Fs': str(Fs), 'identity_exact': True}

J = 1+x*h-sig*(1-r)
H0 = 1+z*(1-b)
Bs = P+x*z*(J-sig)
dn_raw = K*(1+sig*z)*(1+sig*r)-J*P
bn_raw = K*(1+sig*r)*(1+x*h*z)-(1-h)*P
ee = x*h/sig
hcal = (1+x*r)*H0+z*((1-r)*x*h*(2+x*r)
    -ee*(2+x*(1+h+r*(1-h)))+K*ee*ee)
be3 = sp.expand((bn_raw-h*Bs+sig*H0).subs(r,b*z))
be4 = sp.expand((dn_raw+x*bn_raw-sig**2*hcal).subs(r,b*z))
assert be3==0 and be4==0
out['source_identities_BE3_BE4_exact'] = True

# An absolute-coefficient majorant of the correction to the twelve-root
# diagonal Gram, after Ny=Q*JC+b*Jl and R=JC-Jl. The source parameter s
# below stands for sqrt(epsilon); v=sqrt(u). Every displayed term has
# nonnegative integer coefficients, so evaluation at one gives its l1 norm.
s, v, w = sp.symbols('s v w', nonnegative=True)
u = v*v
N, C, L, Ga, A, T, d1, d2, d3, d4, U, V = roots = sp.symbols(
    'N C L Ga A T d1 d2 d3 d4 U V', nonnegative=True)
ny = rr = C+L
f1 = u*Ga+s*d1+u*d2+(u+s)*d3+s*T
f2 = (512*s*s*w*Ga+512*s**4*rr
      +512*s*s*(w+s)*(d2+d3)+(2048*s*s+1024*s*w)*T
      +16*u*w*f1)
Q = ((4*u*u+4*w*w)*ny*ny
     +(64*s*s+16*u*u+4*w*w)*rr*rr
     +(32*s*s+12*u*u+8*w*w)*ny*rr
     +12*N*f1+C*f2+16*u*w*N*ny+32*u*w*N*rr
     +48*s*s*w*N*A+32*u*ny*A+64*u*A*rr
     +sp.Rational(5,2)*s*s*rr*rr+2*v*U*rr
     +6*(16*s*s+u+s)*(d1*d1+d2*d2+d3*d3+d4*d4)
     +6*u*u*rr*rr+(144*u+48*s)*Ga*Ga
     +48*u*U*U+96*v*U*V+48*s*V*V)
Q = sp.expand(Q)
poly = sp.Poly(Q, s,v,w,*roots)
assert all(c>=0 for c in poly.coeffs())
norm = Q.subs({a:1 for a in (s,v,w,*roots)})
Q0 = sp.expand(Q.subs({v:0,w:0}))
norm0 = Q0.subs({a:1 for a in (s,*roots)})
assert Q0.subs(s,0) == 0
assert norm < 2**20 and norm0 < 2**20
support = sorted({(mon[1],mon[2]) for mon,c in poly.terms()})
assert all(0<=a<=4 and 0<=bb<=2 for a,bb in support)
out['gram_majorant'] = {'coefficient_l1':str(norm),
    'constant_testing_part_l1':str(norm0),
    'testing_exponents_sqrt_u_w':support,
    'all_coefficients_nonnegative':True,
    'constant_testing_part_divisible_by_s':True}

epsilon = Fraction(1,2**160)
sqrt_epsilon = Fraction(1,2**80)
delta = Fraction(1,2**100)
sqrt_delta = Fraction(1,2**50)
c0 = Fraction(1,256)
channel_A = 2**30
gamma = c0/(14*channel_A)
assert 2**20*sqrt_epsilon < Fraction(1,192)
assert Fraction(1,2)-Fraction(11,48)-2**20*sqrt_epsilon > c0
assert Fraction(1,96)-2**20*sqrt_epsilon/2 > c0
assert sqrt_delta < gamma
assert sqrt_delta < c0/(7*channel_A)
assert Fraction(1,2)-6*(16*epsilon+sqrt_epsilon)>c0
out['evaluated_constants'] = {'epsilon':str(epsilon),'delta':str(delta),
    'scalar_c0':str(c0),'scalar_A':channel_A,'gamma':str(gamma),
    'main_features':12*(1+2*14),'reduced_main_features':12*(1+2*7),
    'reduced_gamma':str(c0/(7*channel_A)),
    'all_margin_comparisons_exact':True}

# Nonempty source and both zero faces: t>0, x=sigma=t, z=1, b=1/2,
# h=H*t. The numerator polynomials are exact; signs use coefficient
# bounds on H in [1,4] and t<=1/1024, a much larger interval than needed.
t,H = sp.symbols('t H')
sub = {x:t,sig:t,z:1,b:sp.Rational(1,2),h:H*t,r:sp.Rational(1,2)}
dn = sp.expand((K*(1+sig*z)*(1+sig*r)-J*P).subs(sub))
bn = sp.expand((K*(1+sig*r)*(1+x*h*z)-(1-h)*P).subs(sub))
dn_expected = ((3-H)*t*t+(sp.Rational(3,4)-6*H)*t**3
    +(H*H+H/2)*t**4+(-sp.Rational(7,2)*H*H+H/4)*t**5
    +(H**3-H*H/2)*t**6)
bn_expected = ((H-sp.Rational(3,2))*t+4*H*t*t-H*t**3
    +4*H*H*t**4+(-H**3+H*H/2)*t**5)
assert sp.expand(dn-dn_expected)==0
assert sp.expand(bn-bn_expected)==0
def t_remainder_l1(expr,degree):
    q=sp.Poly(sp.expand(expr/t**degree),t,H)
    # Absolute coefficient bound on H in [0,4], t in [0,1].
    return sum(abs(c)*4**mon[1] for mon,c in q.terms())
dn_rest=sp.expand(dn-(3-H)*t*t)
bn_rest=sp.expand(bn-(H-sp.Rational(3,2))*t)
dn_bound=t_remainder_l1(dn_rest,3)
bn_bound=t_remainder_l1(bn_rest,2)
assert dn_bound*sp.Rational(1,1024)<1
assert bn_bound*sp.Rational(1,1024)<sp.Rational(1,4)
out['nonempty_faces'] = {'D_numerator':str(dn),
    'beta_numerator':str(bn),
    'D_remainder_bound_t3_H_0_4':str(dn_bound),
    'beta_remainder_bound_t2_H_0_4':str(bn_bound),
    'D_zero_bracket_H':[2,4],'beta_zero_bracket_H':[1,2],
    'IVT_brackets_valid_on_t_le_1_1024':True}

sub_inf = {x:t*t/16,sig:t,z:1,b:sp.Rational(1,2),
    h:sp.Rational(3,4),r:sp.Rational(1,2)}
dn_inf = sp.factor(dn_raw.subs(sub_inf))
bn_inf = sp.factor(bn_raw.subs(sub_inf))
dn_inf_expected = -3*t*t*(6*t**5-61*t**4+608*t**3-7872*t*t
                         +8192*t-126976)/262144
bn_inf_expected = 3*(3*t*t+64)*(2*t**3+5*t*t+64)/16384
assert sp.expand(dn_inf-dn_inf_expected)==0
assert sp.expand(bn_inf-bn_inf_expected)==0
out['unbounded_H_path'] = {'D_numerator':str(dn_inf),
    'beta_numerator':str(bn_inf),'H':'3/(4t)',
    'KA_comparison_term_exactly_zero':True}

Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
