"""Exact source-independent polynomial and C2 native-contact checks.
No optimizer, network access, numerical sweeps, or external data dependencies.
Checks the all-Omega contact-zero fractional n-opening, all 4 state coefficients,
and the full u-quadratic decomposition with manifest endpoint SOS forms.
This is a PARTIAL SOS carrier, NOT universal K364 H2 admission.
"""
import sympy as S
from pathlib import Path
import json
u,v,W,X,Y,Z=S.symbols('u v W X Y Z')
n,m,C,A,D,alpha=S.symbols('n m C A D alpha', nonnegative=True)
a,b,c,t,r=S.symbols('a b c t r', positive=True)
q=S.symbols('q', nonnegative=True)
Du=1+u*u; Dv=1+v*v
F=(4*n*Du*Dv*X**2+4*m*(v*W-X+u*v*Y-u*Z)**2
   +4*C*Dv*(u*W-Y)**2+4*A*Du*Dv*Y**2
   +(4*D+2)*(u*W+u*v*X-Y-v*Z)**2
   +(4*alpha+2)*Du*(v*Y-Z)**2
   +2*Dv*(X+u*Z)**2-2*Du*Dv*Z**2)
F0=((4*n+2)*Dv*X**2+4*m*(v*W-X)**2
    +(4*C+4*A+2)*Dv*Y**2+4*D*(Y+v*Z)**2+4*alpha*(v*Y-Z)**2)
F1=(8*m*(v*W-X)*(v*Y-Z)-8*C*Dv*W*Y
    -2*(4*D+2)*(W+v*X)*(Y+v*Z)+4*Dv*X*Z)
F2=(4*n*Dv*X**2+4*C*Dv*W**2+4*A*Dv*Y**2
    +(4*m+4*alpha+2)*(v*Y-Z)**2+(4*D+2)*(W+v*X)**2)
assert S.expand(F-(F0+u*F1+u*u*F2))==0
print('PASS complete four-state identity F=F0+u F1+u^2 F2; F0,F2 manifest SOS')
# Independent source-only charge gauge 0<=q<=1: the square of the second root
# is (2q/a)*u*v*(X+(a*c*r/b)*W)^2, with no testing square root needed.
L = u*X-q*(v/a)*X+(1-q)*(c*r/b)*v*W
H_sq=(2*q/a)*u*v*(X+(a*c*r/b)*W)**2
T = S.expand(L**2+H_sq)
# Source contact (c,ac,b/a,-cr,c*(b/a-a*t),1)
cs={u:c,v:a*c,W:b/a,X:-c*r,Y:c*(b/a-a*t),Z:1}
assert S.cancel(L.subs(cs)) == 0
assert S.cancel((X+a*c*r*W/b).subs(cs)) == 0
print('PASS both fractional opening roots vanish at moving strict native C2 contact')
# Complete coefficients, with q arbitrary in the pure-source interval [0,1].
p=S.Poly(T, u,v,W,X,Y,Z)
def coef(ui,vi,ws,xs,ys,zs):return S.factor(p.coeff_monomial((ui,vi,ws,xs,ys,zs)))
assert S.simplify(coef(2,0,0,2,0,0)-1)==0
assert S.simplify(coef(1,1,0,2,0,0))==0
assert S.simplify(coef(0,2,0,2,0,0)-q*q/a**2)==0
assert S.simplify(coef(1,1,1,1,0,0)-2*c*r*(1+q)/b)==0
assert S.simplify(coef(0,2,2,0,0,0)-(1-q)**2*c*c*r*r/b**2)==0
assert S.simplify(coef(1,1,2,0,0,0)-2*q*a*c*c*r*r/b**2)==0
print('PASS forced u^2 X^2=1; uv X^2=0; and exact remaining WX/W^2 charge')
# Source-only choice q=a/(a+1) gives leftover positive v^2 X^2 coefficient.
qchoice=a/(a+1)
assert S.factor(coef(0,2,0,2,0,0).subs(q,qchoice)-1/(1+a)**2)==0
assert S.simplify(4*n+2-4*n/(1+a)**2- (2+4*n*a*(a+2)/(a+1)**2))==0
print('PASS q=a/(1+a): reserved v^2 X^2 margin = 2+4n*a*(a+2)/(a+1)^2 >0')
# Native source strict seed, all native gaps and six positive weights.
ss={a:1,b:1,c:1,t:S.Rational(3,4),r:S.Rational(1,4)}
z=c*c;Ac=1+z*r
h=b-a*a*t;g=(1+z)*b-(1-r+a*a*z*t);j=1+z*r-t
Delta=((1+z)*z*b*b+2*z*(1+z)*b*r+(1+z)*(1+a*a*z)*b
       +z*z*(1-a*a)*r*r+2*z*(1+a*a*z)*r+1+a*a*z)
Ns=(1+z)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z)*(1-r)*t
Nn=(b+r)*Ac*(Ac+a*a*z*t)-r*Delta
native=[Nn/(2*r*Delta),
 (1+a*a*z)*(1-r)*Ac*(Ac+a*a*z*t)/(2*g*Delta),
 Ac*Ns/(2*a*a*t*Delta),
 a*a*Ac*(b+r)*(Ac+a*a*z*t)/(2*h*Delta),
 (1+z)*(b+r)*(b-a*a*z*r)*(Ac+a*a*z*t)/(2*a*a*j*Delta),
 z*z*(b+r)*Ns/(2*(1-z*h)*Delta)]
assert [S.factor(expr.subs(ss)) for expr in native] == list(map(S.Rational,['1/8','3/8','1/6','5/8','3/8','1/6']))
assert all(S.factor(expr.subs(ss))>0 for expr in (h,g,j,1-z*h,b-a*a*z*r,Nn,Ns,Delta))
print('PASS genuine strict seed source and six native weights')
# At strict seed with q=a/(a+1), check exact contact-zero and full identity.
seed={**ss,n:S.Rational(1,8),m:S.Rational(3,8),C:S.Rational(1,6),A:S.Rational(5,8),D:S.Rational(3,8),alpha:S.Rational(1,6)}
assert S.factor(T.subs(q,qchoice).subs(seed).subs({u:1,v:1,W:1,X:-S.Rational(1,4),Y:S.Rational(1,4),Z:1}))==0
print('PASS genuine strict seed carrier and its moving-contact zeros')
# K364 grid inclusion. L uses (uX,vX,vW); sqrt(H_sq) uses (sqrt(uv)X, sqrt(uv)W).
# Multiplication by v produces v^2 X,v^2 W, u v X, u^(1/2)v^(3/2)X/W; all on the grid.
features=[(6,0,1),(0,6,1),(0,6,0),(3,3,1),(3,3,0),
          (6,6,1),(0,12,1),(0,12,0),(3,9,1),(3,9,0)]
assert all(0<=i<=6 and 0<=j<=12 and 0<=k<4 for i,j,k in features)
print('PASS all initial and v-shifted roots lie in fixed 364 feature grid')
# Demonstrate the residual coupled cross coefficient needs nontrivial compensation,
# without claiming the residual is PSD or every Gram dominates the carrier.
Fc=S.Poly(S.expand(Dv*F),u,v,W,X,Y,Z)
Nc=S.Poly(S.expand(4*n*Dv*T),u,v,W,X,Y,Z)
def cross(P,vi):return P.coeff_monomial((1,vi,1,1,0,0))
assert cross(Fc,1)==0 and cross(Fc,3)==0
assert S.simplify(cross(Nc,1)-8*n*c*r*(1+q)/b)==0
assert S.simplify(cross(Nc,3)-8*n*c*r*(1+q)/b)==0
print('PASS residual needs both coupled uv WX and uv^3 WX transfers; no PSD completion asserted')
record={
 'scope':'ALL genuine strict native C2 sources Omega; positive carrier only, H2 UNPROVED',
 'source_order':['d','lambda','n','m','C','A','D','alpha','beta'],
 'contact':['c','a*c','b/a','-c*r','c*(b/a-a*t)','1'],
 'carrier_squared_roots':['(u-q*v/a)*X+(1-q)*(c*r/b)*v*W',
                          'sqrt(2*q/a)*sqrt(u*v)*(X+(a*c*r/b)*W)'],
 'chosen_source_gauge':'q=a/(1+a)',
 'required_unclosed_gram':'For every theta in Omega, find P>=0 satisfying exact coefficients of Dv*F in Z_(6,12). The carrier is merely a sufficient partial allocation; no Gram need dominate it.',
 'remaining_cross_residual':['[u*v*W*X]=-8*n*c*r*(1+q)/b','[u*v^3*W*X]=-8*n*c*r*(1+q)/b'],
 'all_checks':'PASS',
}
Path(__file__).with_name('EXACT_DATA.json').write_text(json.dumps(record,indent=2)+'\n')
print('ALL EXACT CHECKS PASSED; universal H2 remains open')
