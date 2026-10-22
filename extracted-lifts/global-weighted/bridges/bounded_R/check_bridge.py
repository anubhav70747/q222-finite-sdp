#!/usr/bin/env python3
"""Exact checks for the bounded-R bridge, not a regeneration of any pencil.

Requires Python >=3.10 and SymPy.  All arithmetic is exact.  The source
formulas are the original P1 definitions with r=b*z; native admissibility
is used in the pen-and-paper inequalities, not replaced by independent
intervals for the derived variables.  This script checks the new finite
polynomial majorant, denominator removals, and rational reserve constants.
"""
import json
from pathlib import Path
from fractions import Fraction as F
import sys
import sympy as S

B, tau, s, w = S.symbols("B tau s w", nonnegative=True)
u = s**2
Nx, Ny, R, G, At, Tt, d1, d2, d3, d4, U, Kt = roots = S.symbols(
    "Nx Ny R G At Tt d1 d2 d3 d4 U Kt", nonnegative=True
)
f1 = u*G + tau*d1 + u*d2 + (u+tau)*d3 + tau*Tt
f2 = (
    64*B**2*tau**2*w*G + 64*B**2*tau**4*R
    + 64*B**2*tau**2*(w+tau)*(d2+d3)
    + (128*B**2*tau**2 + 128*B**2*tau*w)*Tt
    + 8*B**2*u*w*f1
)
majorant = S.expand(
    (64*B*tau**2 + S.Rational(3,2)*B*tau**2 + 4*u**2)*R**2
    + (64*B**2*tau**2 + 16*B*u**2 + 8*w**2)*Ny*R
    + 8*B*Nx*f1 + (Ny+R)*f2
    + 8*B*u*w*Nx*Ny + 16*B*u*w*Nx*R
    + 40*B**2*tau**2*w*Nx*At + 32*B*u*Ny*At
    + 64*B*u*At*R + 2*B*s*U*R
    + 32*B*u*U**2 + 64*B*s*U*Kt + 32*B*tau*Kt**2
    + ((64*B+32)*u+32*tau)*G**2
    + (64*tau**2+4*u+4*tau)*(d1**2+d2**2+d3**2+d4**2)
)
poly = S.Poly(majorant, B, tau, s, w, *roots)
assert all(c >= 0 for c in poly.coeffs())
support = sorted(set((m[2], m[3]) for m in poly.monoms()))
assert support == [(0,0),(0,1),(0,2),(1,0),(2,0),(2,1),(4,0),(4,1)]
mass = S.expand(majorant.subs({v:1 for v in (tau,s,w,*roots)}))
expected_mass = 1480*B**2 + S.Rational(887,2)*B + 364
assert mass == expected_mass
# For B>=2, write B=t+2 and check coefficientwise positivity.
t = S.symbols("t", nonnegative=True)
gap = S.Poly(S.expand((2048*B**2-mass).subs(B,t+2)),t)
assert all(c > 0 for c in gap.all_coeffs())
constant_part = S.expand(majorant.subs({s:0,w:0}))
assert constant_part != 0
assert all(m[1] >= 1 for m in S.Poly(constant_part,B,tau,*roots).monoms())
# The 220-function organization puts four already retained positive terms
# into the main Gram.  They require only this extra positive majorant.
majorant_220 = S.expand(majorant + 4*B*(u**2+w**2)*Ny**2
                       +(16*B*u**2+4*B*w**2)*R**2)
mass_220 = S.expand(majorant_220.subs({v:1 for v in (tau,s,w,*roots)}))
assert mass_220 == 1480*B**2 + S.Rational(943,2)*B + 364
assert all(c > 0 for c in S.Poly(S.expand(
    (2048*B**2-mass_220).subs(B,t+2)),t).all_coeffs())

# Pole removals and the z-divisibility statement: no PSD matrices rebuilt.
x, sig, z, b, h, e, D = S.symbols("x sig z b h e D")
r = b*z
theta, Rs = 1-h, 1+x*r
P = sig+1+x+sig*z*(1+b*x+2*x*h)+sig*b*x*x*h*z*z+x*x*h*theta*z
M = Rs-x*z*theta
Nk = (1+x)*(1+sig*r)*(1+sig*z)-M*P
Fx = S.cancel((Nk-Nk.subs(x,0))/x)
Fs = S.cancel(Nk.subs(x,0)/sig)
assert S.expand(Fs-(r-1+sig*r*z)) == 0
fxpoly = S.Poly(Fx,x,sig,z,b,h)
fx_mass = sum(abs(c) for c in fxpoly.coeffs())
assert fx_mass == 37 and fxpoly.degree(z) == 3
Np = z*(1+x*h)*(1-e)+(1-r)*(1+x*h*z)
Nl = (1+sig*r)*(e*Rs*(1+sig*z)+M*(1-e))
zquotient = S.cancel((Np-Nl)/z)
assert S.denom(zquotient) == 1
zp = S.Poly(zquotient,x,sig,z,b,h,e)
z_mass, z_degree = sum(abs(c) for c in zp.coeffs()), zp.degree(z)
minus_PR0_over4 = (
    -(1+x)*z*(1+sig*r)*(1-e)
    + 2*D*P*(r-1) + Np + P*(r-1)
)
removed_R0 = (
    z*(1-e)*(1+x*h-(1+x)*(1+sig*r))
    + (1-r)*(1+x*h*z-P)-2*D*P*(1-r)
)
assert S.expand(minus_PR0_over4-removed_R0) == 0

# All scalar conclusions use B>=2; maxima occur at B=2.
epsilon = S.Rational(1,2**40)/B**8
delta = S.Rational(1,2**20)/B**2
envelope = 2048*B**2
constant_error = 2048*B**2*S.Rational(1,2**20)/B**4
reserve = envelope*(2*delta+2*delta**2+delta**3+delta**4+delta**5)
loose_reserve = 7*envelope*delta
assert loose_reserve == S.Rational(7,512)
assert loose_reserve < S.Rational(1,48)
assert constant_error.subs(B,2) == S.Rational(1,2048)
assert S.Rational(1,3)-S.Rational(1,2048)-loose_reserve > S.Rational(1,6)
assert (128*B**2*epsilon).subs(B,2) < S.Rational(1,12)
assert 128*B**2*delta == S.Rational(1,8192) < S.Rational(1,56)

# d_G lower bound in the lambda<=0 branch. Clearing its positive denominators
# gives the polynomial below, >=1-72 B^2 epsilon, itself positive.
eps = S.symbols("eps", positive=True)
dg_numerator = S.expand(
    2*(1-32*B**2*eps)*(1-eps)
    -(1+6*B*eps)*(1-eps)-8*B*eps*(1+6*B*eps)
)
assert S.expand(dg_numerator-(1-(64*B**2+14*B+1)*eps+(16*B**2+6*B)*eps**2)) == 0
assert all(c >= 0 for c in S.Poly(
    S.expand((72*B**2-(64*B**2+14*B+1)).subs(B,t+2)),t).all_coeffs())
assert (72*B**2*epsilon).subs(B,2) < 1

# Exact admissible path exposing the unchanged-principal-Gram obstruction.
# z>=1, sigma=t, x=t^2, r=1/2, h=(z+1/2)t.
tt = S.symbols("tt", positive=True)
H0path = z+S.Rational(1,2)
path = {x:tt**2,sig:tt,b:S.Rational(1,2)/z,h:H0path*tt,e:H0path*tt**2}
Ppath = S.expand(P.subs(path, simultaneous=True))
Jpath = 1+tt**3*H0path-tt/2
NDpath = S.expand((1+tt**2)*(1+tt*z)*(1+tt/2)-Jpath*Ppath)
Dquot = S.cancel(NDpath/tt**2)
assert S.denom(Dquot) == 1
assert S.expand(Dquot.subs(tt,0)-H0path) == 0
Drest = S.Poly(S.expand(Dquot-H0path),tt,z)
assert all(m[0]>=2 for m in Drest.monoms())
assert sum(abs(c) for c in Drest.coeffs()) == 19
assert Drest.degree(z) == 4
Bspath = S.expand(Ppath+tt**2*z*(Jpath-tt))
assert S.expand((Bspath-1).coeff(tt,1)-(z+1)) == 0
betapath = H0path*tt*(Bspath-1)/(2*(1-H0path*tt)*Ppath)
Dpath = NDpath/(2*Jpath*Ppath)
ppath = Np.subs(path,simultaneous=True)/(2*Ppath)
lampath = (8*Dpath*H0path*tt+4*betapath*(1-H0path*tt)**2
           +2*(H0path*tt)**2)/tt-2*ppath*tt
assert S.cancel(S.limit(betapath/tt**2,tt,0)-H0path*(z+1)/2) == 0
assert S.cancel(S.limit(Dpath/tt**2,tt,0)-H0path/2) == 0
assert S.cancel(S.limit(lampath/tt,tt,0)-4*H0path**2) == 0
uu = S.symbols("uu")
minor_limit = 2-2*H0path*uu**2

# An evaluated rational source and a rationally congruent 2x2 matrix,
# for the fixed testing value sqrt(u)=2^-20, w=0.
# Only this new obstruction matrix is evaluated, not an old pencil.
if hasattr(sys,"set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)
uf = F(1,2**40)
zf = 32/uf**2
Bf = zf+1
tf = F(1,2**40)/Bf**8
xf, sf, bf, rf = tf**2, tf, F(1,2)/zf, F(1,2)
Hf, hf = zf+F(1,2), (zf+F(1,2))*tf
thf, ef = 1-hf, xf*hf/sf
Kf, Rsf = 1+xf, 1+xf*rf
Mf, Jf = Rsf-xf*zf*thf, 1+xf*hf-sf*(1-rf)
Pf = sf+Kf+sf*zf*(1+bf*xf+2*xf*hf)+sf*bf*xf**2*hf*zf**2+xf**2*hf*thf*zf
pf = (zf*(1+xf*hf)*(1-ef)+(1-rf)*(1+xf*hf*zf))/(2*Pf)
Df = (Kf*(1+sf*zf)*(1+sf*rf)-Jf*Pf)/(2*Jf*Pf)
betaf = (Kf*(1+sf*rf)*(1+xf*hf*zf)-thf*Pf)/(2*thf*Pf)
kappaf = Kf*(1+sf*rf)*(1+sf*zf)/(2*Mf*Pf)
nuf = (2*Df-2*betaf*thf+hf)/sf
lamf = (8*Df*hf+4*betaf*thf**2+2*hf**2)/sf-4*bf*pf*sf*zf
KAf = 4*kappaf*Hf-lamf
dGf, rhof = 4*Df/sf**2, lamf*xf/sf
Sf = dGf+rhof
assert 0<bf<1 and 0<hf<1 and 0<rf<1
assert 0<ef<1 and Mf>0 and Jf>0 and Df>0 and betaf>0 and lamf>0 and KAf>1
af = -2*nuf*dGf*uf/Sf
minor_det = 2*dGf-af**2
assert minor_det < 0
fixture = {
    "description":"Native rational source; c=0; unchanged principal Gram is indefinite",
    "test":{"u":str(uf),"w":"0","sqrt_u":"1/1048576"},
    "source":{name:str(value) for name,value in
              [("x",xf),("sigma",sf),("z",zf),("b",bf),("h",hf)]},
    "native_guards":{"E":"positive","M":"positive","J":"positive",
                     "D":"positive","beta":"positive","lambda":"positive","KA_minus_1":"positive"},
    "matrix_coordinates":["Nx","Gamma"],
    "congruence_from_Nx_Gamma_star":"diag(1,sqrt(dG))",
    "rational_matrix":[["2",str(af)],[str(af),str(dGf)]],
    "determinant":str(minor_det),
    "determinant_sign":-1,
    "scope":"Obstructs the unchanged principal Gram and radial-only reserve repair; not the full slack",
}
Path(__file__).with_name("fixed_delta_obstruction.json").write_text(
    json.dumps(fixture,indent=2)+"\n")

result = {
    "status": "all exact checks passed",
    "scope": "new bounded-R estimates only; no pencil reassembly or global adequacy claim",
    "constants": {"B":"1+R", "epsilon_R":str(epsilon), "delta_R":str(delta)},
    "positive_testing_support": support,
    "positive_majorant_mass": str(mass),
    "positive_220_majorant_mass": str(mass_220),
    "positive_matrix_envelope": str(envelope),
    "constant_error_max": "1/2048",
    "reserve_upper_bound": str(loose_reserve),
    "nonpositive_matrix_envelope": "128*B**2",
    "nonpositive_scaled_coefficient_bound": "1/8192",
    "Fx_coefficient_mass": int(fx_mass),
    "Fx_z_degree": int(fxpoly.degree(z)),
    "p_minus_l_quotient_coefficient_mass": int(z_mass),
    "p_minus_l_quotient_z_degree": int(z_degree),
    "dG_guard_cleared_numerator": str(dg_numerator),
    "principal_minor_path": {"r":"1/2", "b":"1/(2*z)", "sigma":"t",
                              "x":"t**2", "h":"(z+1/2)*t"},
    "path_D_numerator_remainder_mass": 19,
    "path_beta_over_t2_limit": str(H0path*(z+1)/2),
    "path_D_over_t2_limit": str(H0path/2),
    "path_lambda_over_t_limit": str(4*H0path**2),
    "path_Nx_Gamma_minor_limit": str(minor_limit),
    "evaluated_rational_obstruction": "fixed_delta_obstruction.json",
    "evaluated_obstruction_sqrt_u": "2^-20",
    "evaluated_obstruction_z": str(zf),
    "evaluated_obstruction_determinant_sign": -1,
}
out = Path(__file__).with_name("check_results.json")
out.write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps(result,indent=2))
