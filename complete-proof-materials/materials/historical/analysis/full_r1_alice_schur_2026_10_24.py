"""Exact two-normal Schur elimination and its finite-feature obstruction.

No solver, sampling sweep, source search, or public edit. The two rational
source evaluations establish nonconstancy of one actual native family;
they do not estimate positivity on the missing testing strip.
"""
from pathlib import Path
from fractions import Fraction as F
import json,runpy
import sympy as s
HERE=Path(__file__).resolve().parent
project=HERE.parent
u,v,N,M,w,z=s.symbols('u v N M w zeta')
n,m,C,A,D,al,p,H=s.symbols('n m C A D alpha p H')
Du=1+u*u;Dv=1+v*v;q=al+s.Rational(1,2)
t=N+p*v*w;r=-M-H*v*z
W=w-u*r;X=t-u*z;Y=r+u*w;Z=z+u*t
TD=Du*(r+v*z)**2;TA=(v*Y-Z)**2;TB=Du*Dv*t*t
T=n*Dv*X*X+m*Du*(v*w-t)**2+C*Du*Dv*r*r+A*Dv*Y*Y+D*TD+al*TA+(TD+TA+TB-Dv*Z*Z)/2
RN=m*Du+(n+s.Rational(1,2))*Dv
RM=Du*(C*Dv+D+s.Rational(1,2))+A*Dv
hN=q*u*u+RN;hM=q*v*v+RM;h=q*u*v
LN=(v*w*(n*p*Dv-m*(1-p)*Du-q*u*u*(1-p)+p*Dv/2)
    +u*z*(al*(1+H*v*v)-n*Dv+(H-1)*v*v/2))
LM=(v*z*(Du*(H*C*Dv-(1-H)*(D+s.Rational(1,2)))+H*A*Dv+q*(1+H*v*v))
    -u*w*(A*Dv+q*(1-p)*v*v))
F0=T.subs({N:0,M:0})
assert s.expand(T-hN*N*N-2*h*N*M-hM*M*M-2*LN*N-2*LM*M-F0)==0
Delta=hN*hM-h*h
assert s.expand(Delta-RN*RM-q*u*u*RM-q*v*v*RN)==0

# Generic fraction-free block identity: this operation does not claim
# an inverse of a mixed source/testing expression in the tensor ring.
a,b,c,L1,L2,Fg=s.symbols('a b c L1 L2 Fg')
G=s.Matrix([[a,b],[b,c]]);adj=s.Matrix([[c,-b],[-b,a]])
zz=s.Matrix([N,M]);ll=s.Matrix([L1,L2]);det=a*c-b*b
assert s.expand(det*((zz.T*G*zz)[0]+2*(zz.T*ll)[0]+Fg)
                -((G*zz+ll).T*adj*(G*zz+ll))[0]-det*Fg+(ll.T*adj*ll)[0])==0

# Exact v=0 form and the rational entry of the separated Schur square.
AN=n+m+s.Rational(1,2);BN=m+al+s.Rational(1,2)
AM=C+D+A+s.Rational(1,2);BM=C+D+s.Rational(1,2)
vzero=(AN+BN*u*u)*N*N+2*(al-n)*u*N*z+(n*u*u+al)*z*z
vzero+=(AM+BM*u*u)*M*M-2*A*u*M*w+A*u*u*w*w
assert s.expand(T.subs(v,0)-vzero)==0
assert s.expand((AN+BN*u*u)*(n*u*u+al)-(n-al)**2*u*u
     -(al*AN+((n+al)*(m+s.Rational(1,2))+2*n*al)*u*u+n*BN*u**4))==0

# Actual native family j=1,H=1/2,B=T=eps, p in[1/4,1/2].
dictionary=project/'github-q222-finite-sdp/extracted-lifts/canonical-completion/uniform-contact/source_dictionary.py'
ctx=runpy.run_path(str(dictionary))
eps=s.symbols('eps',positive=True)
family={ctx['j']:1,ctx['H']:s.Rational(1,2),ctx['B']:eps,ctx['T']:eps}
pp=ctx['p'];NN=s.expand(ctx['NN'].subs(family))
error=sum(abs(co)*eps**power_eps*s.Rational(1,2)**power_p
          for (power_eps,power_p),co in s.Poly(NN-pp,eps,pp).terms())
eps0=s.Rational(1,2**20)
assert error.subs(eps,eps0)<8*eps0<s.Rational(1,4)
ratios=[];positive=[]
for pv in [s.Rational(1,4),s.Rational(1,2)]:
    sub={**family,eps:eps0,pp:pv}
    # Substitution is sequential because family introduces eps.
    nn=s.cancel(ctx['src'][ctx['n']].subs(family).subs({eps:eps0,pp:pv}))
    mm=s.cancel(ctx['src'][ctx['m']].subs(family).subs({eps:eps0,pp:pv}))
    aa=s.cancel(ctx['src'][ctx['alpha']].subs(family).subs({eps:eps0,pp:pv}))
    assert nn>aa>0
    ratios.append(s.cancel((nn+mm+s.Rational(1,2))/(mm+aa+s.Rational(1,2))))
    positive.append({'p':str(pv),'n_minus_alpha_positive':True,'ratio':str(ratios[-1])})
assert ratios[0]!=ratios[1]

# Outer contact separation alone is source/testing polynomial arithmetic.
delta=F(1,1024);s0=F(81,65536)
assert delta/s0==F(64,81) and 1-delta/s0==F(17,81)
out={'full_parent_two_normal_block_and_all_linear_terms_exact':True,
 'normal_determinant_positive_sum_exact':True,'generic_fraction_free_Schur_identity_exact':True,
 'v_zero_parent_exact':True,'v_zero_zeta_Schur_numerator_positive_coefficients_exact':True,
 'native_family':{'B':str(eps0),'T':str(eps0),'j':'1','H':'1/2','p_interval':['1/4','1/2'],
    'physicality_N_error_upper':str(error.subs(eps,eps0)),
    'physicality_N_gt_p_minus_8eps':True,'two_exact_ratio_values':positive,
    'normal_pivot_ratio_nonconstant':True},
 'outer_gap_factor':'17/81',
 'conclusion':'Literal separated Schur-positive component has no one finite source-independent testing family, by its infinitely many distinct rational poles. This is not an obstruction to a different full-parent fixed Gram.',
 'no_missing_strip_admission_claim':True}
outfile=Path(__file__).with_suffix('.json')
outfile.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps(out,indent=2))
