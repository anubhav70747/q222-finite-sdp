"""Small exact checks accompanying bounded_z_R_extension.tex.

Checks the source path, the retained w=0 repair, and the polynomial
coefficient bound used for the unbounded-H source estimates.  It is not
a global Q222 inclusion test or a substitute for the displayed inequalities.
"""
from pathlib import Path
import json
import sympy as S

x, sig, z, b, h, t, R = S.symbols('x sig z b h t R', positive=True)
u, w, W, X, Y, Z = S.symbols('u w W X Y Z', real=True)
m, C, A, D, al, be = S.symbols('m C A D al be', nonnegative=True)
r=b*z; th=1-h; K=1+x
J=1+x*h-sig*(1-r); M=1+x*r-x*z*th
P=sig+K+sig*z*(1+b*x+2*x*h)+sig*b*x*x*h*z*z+x*x*h*th*z
Nk=S.expand(K*(1+sig*r)*(1+sig*z)-M*P)
Fs=r-1+sig*r*z
Fx=S.cancel((Nk-sig*Fs)/x)
assert S.denom(Fx)==1
pol=S.Poly(Fx,x,sig,z,b,h)
fxnorm=sum(abs(c) for c in pol.coeffs())
fxzdegree=S.degree(Fx,z)
assert fxnorm==37 and fxzdegree<=3

F=(4*m*(W-w*X+u*Y-u*w*Z)**2
  +4*C*(1+w*w)*(u*W-Y)**2+4*A*(1+u*u)*(1+w*w)*Y*Y
  +(4*D+2)*(u*w*W+u*X-w*Y-Z)**2
  +(4*al+2)*(1+u*u)*(Y-w*Z)**2
  +(4*be+2)*(1+w*w)*(X+u*Z)**2
  -2*(1+u*u)*(1+w*w)*Z*Z)
F0=(4*m*(W+u*Y)**2+4*C*(u*W-Y)**2
  +4*(A+al)*(1+u*u)*Y*Y+2*(1+u*u)*(X*X+Y*Y)
  +4*D*(u*X-Z)**2+4*be*(X+u*Z)**2)
assert S.expand(F.subs(w,0)-F0)==0
F2=(4*m*(X+u*Z)**2+(4*C+4*D+2)*(u*W-Y)**2
  +4*A*(1+u*u)*Y*Y+4*al*(1+u*u)*Z*Z
  +(4*be+2)*(X+u*Z)**2)
F1=(-8*m*(W+u*Y)*(X+u*Z)
  +2*(4*D+2)*(u*X-Z)*(u*W-Y)
  -2*(4*al+2)*(1+u*u)*Y*Z)
assert S.expand(F-F0-w*F1-w*w*F2)==0

e=x*h/sig
p=(z*(1+x*h)*(1-e)+(1-r)*(1+x*h*z))/(2*P)
d=(K*(1+sig*z)*(1+sig*r)-J*P)/(2*J*P)
beta=(K*(1+sig*r)*(1+x*h*z)-th*P)/(2*th*P)
kap=K*(1+sig*r)*(1+sig*z)/(2*M*P)
lam=(8*d*h+4*beta*th*th+2*h*h)/sig-4*b*p*sig*z
nu=(2*d-2*beta*th+h)/sig
KA=4*kap*h/sig-lam
path={x:t*t,sig:t,z:R,b:1/(2*R),h:S.Rational(1,2)}
path_limit=lambda expr:S.simplify(S.limit(S.cancel(expr.subs(path)),t,0,dir='+'))
limits={
 'D_over_t2':path_limit(d/t**2),'beta':path_limit(beta),
 't_lambda':path_limit(t*lam),'dG':path_limit(4*d/sig**2),
 'rho':path_limit(lam*x/sig),'nu':path_limit(nu),'KA':path_limit(KA)}
expected={'D_over_t2':R/2,'beta':S.Rational(1,2),'t_lambda':S.Integer(1),
 'dG':2*R,'rho':S.Integer(1),'nu':R+S.Rational(1,2),'KA':R}
assert all(S.simplify(limits[k]-v)==0 for k,v in expected.items())
det_limit=S.simplify(2-4*limits['nu']**2*u**2*limits['dG']/(limits['dG']+limits['rho'])**2)
assert S.simplify(det_limit-2*(1-R*u*u))==0

out={'scope':'Exact identities and path limits only; the R-box inequalities are proved in the TeX.',
 'Fx_coefficient_l1':int(fxnorm),'Fx_z_degree':int(fxzdegree),
 'w0_full_parent_SOS_identity':True,'w_quadratic_expansion':True,
 'source_path_limits':{k:str(v) for k,v in limits.items()},
 'unchanged_principal_minor_limit':str(det_limit),
 'testing_box_independent_of_R_obstructed_for_unchanged_principal_Gram_only':True}
(Path(__file__).parent/'bounded_z_R_extension_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
