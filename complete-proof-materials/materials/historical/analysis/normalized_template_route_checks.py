"""Exact native subfamily and continuously rotating event-square limit.

No optimizer, numerical sweep or public edit. The complete parent is used
to isolate its m-event coefficient; all remaining native coefficients have
finite p->1 limits by their displayed positive denominators.
"""
from pathlib import Path
import json,runpy
import sympy as s
HERE=Path(__file__).resolve().parent
project=HERE.parent
dct=runpy.run_path(str(project/'github-q222-finite-sdp/extracted-lifts/canonical-completion/uniform-contact/source_dictionary.py'))
j,p=dct['j'],dct['p'];eps=s.Rational(1,4096);H=s.Rational(1,2)
sub={dct['B']:eps,dct['T']:eps,dct['H']:H}
ee=s.symbols('eps',positive=True)
NN=s.expand(dct['NN'].subs({dct['B']:ee,dct['T']:ee,dct['H']:H}))
error=sum(abs(co)*3**pj*ee**pe for (pe,pj,pp),co in s.Poly(NN-p,ee,j,p).terms())
error=s.cancel(error.subs(ee,eps))
assert error<24*eps<s.Rational(3,4)
assert 3*eps<s.Rational(1,1024)
assert s.cancel(dct['src'][dct['m']]*(1-p)-dct['src'][dct['mu']])==0

S,T,ew,ex,ey,ez=s.symbols('s t etaW etaX etaY etaZ')
uc=s.sqrt(j*eps);vc=eps;y=eps**2
Wc=1+H*j*y;Xc=-(j-p)*eps
Yc=(1-H*eps)*uc;Zc=(1+p*eps)*uc
u=uc*S;v=vc*T
W=Wc*ew;X=Xc*ex;Y=Yc*ey;Z=Zc*ez
L=(1+H*j*y)*T*ew+(j-p)*ex+j*eps*(1-H*eps)*S*T*ey-j*(1+p*eps)*S*ez
event=v*W-X+u*v*Y-u*Z
assert s.expand(event-eps*L)==0
n,m,C,A,D,al=s.symbols('n m C A D alpha')
Du=1+u*u;Dv=1+v*v
F=(4*n*Du*Dv*X*X+4*m*event**2+4*C*Dv*(u*W-Y)**2+4*A*Du*Dv*Y*Y
   +(4*D+2)*(u*W+u*v*X-Y-v*Z)**2+(4*al+2)*Du*(v*Y-Z)**2
   +2*Dv*(X+u*Z)**2-2*Du*Dv*Z*Z)
assert s.expand(s.diff(F,m)-4*y*L*L)==0
Llim=L.subs(p,1)
L0=T*ew-ex
L1=H*y*T*ew+ex+eps*(1-H*eps)*S*T*ey-(1+eps)*S*ez
assert s.expand(Llim-L0-j*L1)==0
assert s.expand(L0.subs({S:1,T:1,ew:1,ex:1,ey:1,ez:1}))==0
assert s.expand(L1.subs({S:1,T:1,ew:1,ex:1,ey:1,ez:1}))==0
k=s.symbols('k')
vec=s.Matrix([1+H*j*y,j-1,j*eps*(1-H*eps),-j*(1+eps)])
det=s.det(s.Matrix([[vec[0],vec[0].subs(j,k)],[vec[3],vec[3].subs(j,k)]]))
assert s.expand(det-(j-k)*(1+eps))==0
out={'native_family':{'B':str(eps),'T':str(eps),'H':'1/2','j':'[2,3]','p':'[3/4,1)'},
 'jB_lt_delta':True,'physicality_N_error_upper':str(error),
 'physicality_N_gt_p_minus_24B':True,'m_times_1_minus_p_equals_positive_mu':True,
 'complete_parent_m_coefficient_exact':True,'normalized_m_root':str(L),
 'limit_root':str(Llim),'fixed_two_root_decomposition_exact':True,
 'fixed_roots':[str(L0),str(L1)],'rank_one_direction_minor':'(j-k)*(1+B)',
 'finite_fixed_PSD_template_conic_hull_impossible':True,
 'scope':'Obstruction to a finite polyhedral positive-template cover of all normalized native targets; not an obstruction to fixed functions or L6.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps(out,indent=2))
