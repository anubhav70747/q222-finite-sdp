"""Independent exact algebra for inverse_alice_all_z_completion.tex.

This checks identities and rational margins, not global Q222 adequacy.
The elementary source inequalities are proved in the accompanying note.
"""
from fractions import Fraction
from pathlib import Path
import json
import sympy as s

p,w,W,X,Y,Z=s.symbols('p w W X Y Z',real=True)
q,b,Q=s.symbols('q b Q',positive=True)
z=q*q; r=b*z; H=1+z*(1-b)
C0=z/(2*(1-r)); A0=H/(2*b*z); al0=1/(2*(1-b)*z)
m,C,A,D,al,be,f=s.symbols('m C A D al be f',real=True)
E=q*W; F=Y/q; G=w*Z/q; R=X+w*W; S=p*X
K=(R*R+(1+w*w)*S*S+z*(G-F)**2+p*p*(1+w*w)*z*F*F
   -2*p*E*(G+w*w*F)+2*q*w*S*(G-F))
K0=R*R+S*S+z*(G-F)**2+p*p*z*F*F-2*p*E*G
raw=(4*m*(p*W-p*w*X+Y-w*Z)**2
     +4*C*(1+w*w)*(W-p*Y)**2
     +4*A*(1+p*p)*(1+w*w)*Y*Y
     +(4*D+2)*(w*W+X-p*w*Y-p*Z)**2
     +(4*al+2)*(1+p*p)*(Y-w*Z)**2
     +(4*be+2)*(1+w*w)*(p*X+Z)**2
     -2*(1+p*p)*(1+w*w)*Z*Z)
core=raw.subs({m:0,C:0,A:0,D:0,al:0,be:0})
assert s.expand(core-2*K)==0
Q0=2*K0+4*C0/z*(E-p*z*F)**2+4*A0*z*(1+p*p)*F*F+4*al0*z*(1+p*p)*(F-G)**2
J0=F-b*G; B=H/(b*(1-b))
leading=2*R*R+2*S*S+2/(1-r)*(E-p*(z*J0+G))**2+2*B*(1+p*p)*J0*J0+2*H*G*G
assert s.cancel(Q0-leading)==0
retained=(4*m*(p*W-p*w*X+Y-w*Z)**2
   +4*D*(w*W+X-p*w*Y-p*Z)**2
   +4*be*(1+w*w)*(p*X+Z)**2
   +4*(C-f*C0)/z*(1+w*w)*(E-p*z*F)**2
   +4*(A-f*A0)*z*(1+p*p)*(1+w*w)*F*F
   +4*(al-f*al0)*z*(1+p*p)*(F-G)**2)
positive=(2*(1-f)*(R*R+S*S+z*(G-F)**2+p*p*z*F*F)
   +2*w*w*S*S+2*p*p*w*w*z*F*F
   +4*f*C0/z*w*w*(E-p*z*F)**2+4*f*A0*z*(1+p*p)*w*w*F*F)
error=-4*(1-f)*p*E*G-4*p*w*w*E*J0-4*p*b*w*w*E*G-4*q*w*S*J0+4*q*(1-b)*w*S*G
assert s.cancel(raw-f*Q0-error-positive-retained)==0

# The normalized baseline identity is polynomial and requires no sign assumptions.
eb,js,gs,rr,ss,t,a,c=s.symbols('eb js gs rr ss t a c',real=True)
L=a*js+c*gs
baseline=2*rr**2+2*ss**2+2*(eb-t*L)**2+2*(1+t*t)*js*js+2*gs*gs
baseline_sos=(s.Rational(3,2)*(rr*rr+ss*ss)
  +s.Rational(3,2)*(eb-s.Rational(4,3)*t*L)**2
  +s.Rational(3,2)*js*js+s.Rational(5,6)*gs*gs
  +s.Rational(2,3)*(1-t*t)*gs*gs+s.Rational(4,3)*t*t*js*js
  +s.Rational(2,3)*t*t*((1-a*a-c*c)*(js*js+gs*gs)+(c*js-a*gs)**2))
assert s.expand(baseline-s.Rational(1,2)*(eb*eb+js*js+gs*gs+rr*rr+ss*ss)-baseline_sos)==0
assert s.expand(H-1-b*(1-b)*z*z-z*(1-b)*(1-r))==0

# Three source matrices encode exactly the full remaining correction.
cut,hr,br=s.symbols('cut hr br',positive=True)
y=s.Matrix([eb,js,gs,rr,ss])
C1=s.zeros(5); C2=s.zeros(5); C3=s.zeros(5)
C1[0,2]=C1[2,0]=-2*(1-f)*cut**2/hr
C2[0,1]=C2[1,0]=-2*cut**2/br
C2[0,2]=C2[2,0]=-2*b*cut**2/hr
C3[1,4]=C3[4,1]=-2*q/br
C3[2,4]=C3[4,2]=2*q*(1-b)/hr
formalerror=s.expand((y.T*(t*C1+t*w*w*C2+w*C3)*y)[0])
expected=(-4*(1-f)*cut**2*t*eb*gs/hr
 -4*cut**2*t*w*w*eb*js/br-4*b*cut**2*t*w*w*eb*gs/hr
 -4*q*w*ss*js/br+4*q*(1-b)*w*ss*gs/hr)
assert s.expand(formalerror-expected)==0

eta=Fraction(1,2**30)
delta=Fraction(1,2**12)
reserve=50*eta+4*delta*delta+4*delta
assert reserve<Fraction(1,512)
assert (1-25*eta)/2-reserve>Fraction(1,4)
out={'scope':'Exact local inverse-Alice identities; no full-Q222 inclusion claim.',
 'inverse_parent_equals_twice_CE20_core_plus_events':True,
 'leading_completion_all_z':True,'strict_baseline_SOS':True,
 'source_normalization_gap_zQ_times_one_minus_r':True,
 'complete_retained_square_identity':True,'three_coefficient_matrices_exact':True,
 'worst_case_reserve_at_cutoff_one':str(reserve),
 'reserve_below_one_over_512':True,'final_margin_above_one_over_4':True,
 'general_cutoff_source_eta':'2^-30 / P_*^2',
 'general_cutoff_testing_delta':'2^-12 / P_*',
 'native_source_condition':'0 < x,sigma <= eta/(1+z)',
 'testing_family_independent_of_z_and_R':True}
(Path(__file__).parent/'inverse_alice_all_z_completion_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
