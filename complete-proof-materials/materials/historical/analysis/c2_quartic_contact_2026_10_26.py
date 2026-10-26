"""One exact quartic contact check; no hull, converse, or SDP claim."""
import json
from pathlib import Path
import sympy as S

folder=Path(__file__).parent
prior=json.loads((folder/'c2_maxentangled_general_angles_2026_10_25.json').read_text())
a,r,z=S.symbols('a r z',positive=True)
local={'a':a,'r':r,'z':z}
clean=[S.sympify(t,locals=local) for t in prior['clean_source']]
scale=S.sympify(prior['positive_common_polynomial_scale'],locals=local)
p=[S.Poly(S.cancel(t*scale),a,r,z,domain=S.QQ) for t in clean]
dd,nn,mm,CC,AA,DD,al=p
half=dd.mul_ground(S.Rational(1,2))
bb=AA+nn;UU=CC+half;VV=al+half;SS=UU+VV;MM=mm+DD
angle=a*(dd+2*nn)-2*AA
old= a*a*MM*(bb+SS)-bb*(MM+SS)
quartic=SS*SS*(bb+half)*(MM+half)-(a+1)**2*UU*VV*MM*(bb+SS)
assert angle.is_zero and old.is_zero and quartic.is_zero

k,t,y,theta,x=S.symbols('k t y theta x')
D0=2*k*t+k+t
W=(k+t+1)**2*x*theta-(x+(k-t)**2*theta)*(x+theta)
quadratic=(x-D0*theta)**2-4*k*t*(k+1)*(t+1)*theta**2
assert S.expand(W+quadratic)==0
p0=t*(k+1); q0=k*(t+1)
assert S.expand(D0-p0-q0)==0
factored=(x-(a+1)**2*p0*theta)*(x-(a-1)**2*p0*theta)
angle_relation=a*a*p0-q0
assert S.rem(S.Poly(S.expand(factored-quadratic),a),
             S.Poly(angle_relation,a)).is_zero

counter=json.loads((folder/'c2_quadratic_constraints_counterexample_2026_10_26.json').read_text())
base=S.Matrix([S.Rational(t) for t in counter['base']])
direction=S.Matrix([S.Rational(t) for t in counter['direction']])
eps=S.Rational(counter['epsilon']); av=S.Rational(counter['a'])
delta=S.symbols('epsilon',real=True)
def quantities(c):
    d,n,m,C,A,D,al=c
    B=A+n;U=C+d/2;V=al+d/2;SS=U+V;M=m+D
    F=SS**2*(B+d/2)*(M+d/2)-(av+1)**2*U*V*M*(B+SS)
    return d,B,U,V,SS,M,S.factor(F)
d,B,U,V,SS,M,Fbase=quantities(base)
assert Fbase==0
shifted=base+delta*direction
_,Bs,Us,Vs,Ss,Ms,Fshift=quantities(shifted)
assert [Bs,Us,Vs,Ss,Ms]==[B,U,V,SS,M]
formula=SS**2*(delta*(B+M+d)/2+delta**2/4)
assert S.factor(Fshift-formula)==0
counter_residual=S.factor(Fshift.subs(delta,eps))
assert counter_residual>0
kn=S.factor(B/SS);tn=S.factor(M/SS);thn=S.factor(U*V/SS**2)
yn=S.factor((d+eps)/(2*SS));xn=S.factor((kn+yn)*(tn+yn))
weighted=S.factor(W.subs({k:kn,t:tn,theta:thn,x:xn}))
assert weighted<0
out={
    'full_three_parameter_quartic_verified_as_polynomial':True,
    'positive_source_clearing_factor':str(scale),
    'angle_and_summed_minor_identity_verified':True,
    'weighted_CHSH_factorization_verified':True,
    'normalized_definitions':{'k':'B/S','t':'M/S','y':'d/(2*S)',
                             'theta':'U*V/S^2','x':'(k+y)*(t+y)'},
    'upper_root':'(a+1)^2*t*(k+1)*theta',
    'lower_root':'(a-1)^2*t*(k+1)*theta',
    'root_factorization_modulo_angle_relation_verified':True,
    'branch_reason':'p=t(k+1)>0, q=k(t+1)=a^2*p, sqrt(p*q)=a*p since a>0',
    'counterexample_a':str(av),
    'counterexample_epsilon':str(eps),
    'counterexample_invariants':{'B':str(B),'U':str(U),'V':str(V),'S':str(SS),'M':str(M)},
    'quartic_residual_along_invisible_direction':str(S.factor(Fshift)),
    'counterexample_quartic_residual':str(counter_residual),
    'counterexample_normalized_weighted_CHSH_polynomial':str(weighted),
    'scope':'Necessary identity for individual native contacts; no converse or full-hull condition proved.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))
