"""Independent symbolic expansions of selected 94/95/97 and 38 identities.

This script contains no calls to included manuscript checkers.  It checks exact
algebra, not positivity, tensor-ring membership, or atlas exhaustion.
"""
import json
from pathlib import Path
import sympy as s

checks = []

def zero(name, expression):
    value = s.cancel(s.expand(expression))
    ok = value == 0
    checks.append({"name": name, "passed": ok})
    if not ok:
        raise AssertionError((name, s.factor(value)))

u,v,W,X,Y,Z,n,m,C,A,D,al,be,r,h,et = s.symbols(
    'u v W X Y Z n m C A D al be r h et')
Du,Dv = 1+u*u,1+v*v
M = v*W-X+u*v*Y-u*Z
J = u*W+u*v*X-Y-v*Z
core = 2*(J**2+Du*(v*Y-Z)**2+Dv*(X+u*Z)**2-Du*Dv*Z**2)
F = 4*n*Du*Dv*X**2+4*m*M**2+4*C*Dv*(u*W-Y)**2+4*A*Du*Dv*Y**2+4*D*J**2+4*al*Du*(v*Y-Z)**2+core

# Native C1 dictionary SN1--SN2, with the original beta event retained.
Fc1 = F.subs(n,0)+4*be*Dv*(X+u*Z)**2
gr1 = {X:-h*u*Z, Y:r*u*W}
theta = 1-h
fw = -m*(1+r*u*u)+(D+s.Rational(1,2))*(1-r)*u*u
fz = m*theta-(D+s.Rational(1,2))*v*v*(1+h*u*u)+(be+s.Rational(1,2))*Dv*theta
kw = m*v*v*(1+r*u*u)-C*Dv*(1-r)+A*Du*Dv*r-(D+s.Rational(1,2))*(1-r)+(al+s.Rational(1,2))*Du*r*v*v
kz = -m*theta*u*u+(D+s.Rational(1,2))*(1+h*u*u)-(al+s.Rational(1,2))*Du
Hxx=m+be*Dv+(D+s.Rational(1,2))*u*u*v*v+Dv/2
Hyy=m*u*u*v*v+C*Dv+A*Du*Dv+D+s.Rational(1,2)+(al+s.Rational(1,2))*Du*v*v
Hxy=-(m+D+s.Rational(1,2))*u*v
zero('SN2 full X first normal', s.diff(Fc1,X).subs(gr1)-8*(v*W*fw+u*Z*fz))
zero('SN2 full Y first normal', s.diff(Fc1,Y).subs(gr1)-8*(u*W*kw+v*Z*kz))
zero('SN2 XX', s.diff(Fc1,X,2)/8-Hxx)
zero('SN2 YY', s.diff(Fc1,Y,2)/8-Hyy)
zero('SN2 XY', s.diff(Fc1,X,Y)/8-Hxy)

# The homogeneous scalar lemma.  t is sqrt(B), so radicals do not obscure
# the universal polynomial identities.
t,sg,k = s.symbols('t sg k', nonzero=True)
B=t*t; x=sg*B; xi=u*u; tau=v-k
R=(u*W-B*Z/k)/t
S=u*Z-sg*k*W
U=t*tau*Z
T=u*tau*Z
zero('HL4 W conductor', (xi-x)*W-u*t*R-B*S/k)
zero('HL4 Z conductor', (xi-x)*Z-u*S-sg*k*t*R)
zero('HL11 exact angular zero', T*T-sg*U*U-S*T*tau-sg*k*R*U*tau)
d,ss,hh=s.symbols('d ss hh')
cc=ss*(xi-x)-B*d/k
aa=(hh*(xi-x)**2-k*xi*d-sg*k*k*cc-sg*B*k*d)/B
qq=-k*d-sg*k*k*ss
zero('HL7 entire scalar identity', aa*W*W+2*u*d*W*Z+cc*Z*Z-qq*R*R-ss*S*S-hh*(u*R+t*S/k)**2)
Q1,a01,d1,b1,Jxi,Jv=s.symbols('Q1 a01 d1 b1 Jxi Jv')
L1=a01/k**2+b1
e1=(xi-x)*Jxi+tau*Jv-sg*L1-Q1/k**2-2*d1/k
D1=L1+B*Jxi
df=tau*((xi*Q1/B+sg**2*a01)*W**2+2*u*d1*W*Z+(B*e1+xi*b1)*Z**2)
lift=tau*Q1*R**2+tau*a01*S**2/k**2+(2*(Q1/k+d1)+sg*k*D1)*R*U+(D1-2*a01/k**2)*S*T+Jv*U**2
zero('HL10 entire angular coquotient',df-lift)

# Original rational source formulas; independently substitute the two source
# scales to compare with 97's polynomials and 95's cancelled formulas.
xs,bs,hs,ks,rs=s.symbols('xs bs hs ks rs', nonzero=True)
a2=ks**2/xs; ts=xs*(bs-hs)/ks**2
Ac=1+xs*rs; Sc=Ac+a2*xs*ts
gs=bs-1+rs+xs*hs; js=1+xs*rs-ts
Delta=(1+xs)*xs*bs**2+2*xs*(1+xs)*bs*rs+(1+xs)*(1+a2*xs)*bs+xs**2*(1-a2)*rs**2+2*xs*(1+a2*xs)*rs+1+a2*xs
Ns=(1+xs)*(bs+rs)*(bs+a2*(1-ts))-a2*(1+a2*xs)*(1-rs)*ts
old={
 'nr':((bs+rs)*Ac*Sc-rs*Delta)/(2*Delta),
 'm':(1+a2*xs)*(1-rs)*Ac*Sc/(2*gs*Delta),
 'C':Ac*Ns/(2*a2*ts*Delta),
 'A':a2*Ac*(bs+rs)*Sc/(2*hs*Delta),
 'D':(1+xs)*(bs+rs)*(bs-a2*xs*rs)*Sc/(2*a2*js*Delta),
 'alpha':xs**2*(bs+rs)*Ns/(2*(1-xs*hs)*Delta)}
BB,HH,LL=s.symbols('BB HH LL')
sub97={bs:(HH+LL)/xs,hs:HH/xs}
zero('SR x squared N polynomial',xs**2*Ns.subs(sub97)-((1+xs)*(HH+LL+xs*rs)*(HH+ks**2)-xs*(1+ks**2)*(1-rs)*LL))
delta97=(1+xs)*xs*bs+2*xs*(1+xs)*rs+(1+xs)*(1+ks**2)+(xs**2-ks**2*xs)*rs**2/bs+2*xs*(1+ks**2)*rs/bs+(1+ks**2)/bs
zero('SR Delta over b polynomial', Delta/bs-delta97)
sig,Bs,eta=s.symbols('sig Bs eta',nonzero=True)
sub95={xs:sig*Bs,bs:1/sig,hs:eta/sig}
dh=s.cancel(sig*Delta.subs(sub95)); nh=s.cancel(sig**2*Bs*Ns.subs(sub95))
ac=1+sig*Bs*rs; z=Bs*(1-eta); y=Bs*eta; sc=ac+z
gh=1+sig*(rs-1+Bs*eta); jj=ac-z/ks**2
new={
 'nr':((1+sig*rs)*ac*sc-rs*dh)/(2*dh),
 'm':sig**2*(1+ks**2)*(1-rs)*ac*sc/(2*gh*dh),
 'C':ac*nh/(2*z*dh),
 'A':ks**2*ac*(1+sig*rs)*sc/(2*y*dh),
 'D':Bs*(1+sig*Bs)*(1+sig*rs)*(1-sig*ks**2*rs)*sc/(2*ks**2*jj*dh),
 'alpha':Bs*(1+sig*rs)*nh/(2*(1-y)*dh)}
for key in old:
    zero('UB1 '+key,s.cancel(old[key].subs(sub95))-new[key])
residue_exprs={
 'nr':(new['nr'], (1/(1+ks**2)-rs)/2),
 'm over sigma squared':(new['m']/sig**2,(1-rs)/2),
 'Cz':(new['C']*z,ks**2/(2*(1+ks**2))),
 'Ay':(new['A']*y,ks**2/(2*(1+ks**2))),
 'D over B':(new['D']/Bs,1/(2*ks**2*(1+ks**2))),
 'alpha over B':(new['alpha']/Bs,ks**2/(2*(1+ks**2)))}
for key,(expr,want) in residue_exprs.items():
    zero('UB2 '+key,s.cancel(expr).subs({sig:0,Bs:0})-want)
zero('UB3 exact Nr identity',4*new['nr']+2*rs-2*(1+sig*rs)*ac*sc/dh)

# Entire X row, not only the scalar source graph.
graph={X:-r*u*Z,Y:et*u*W}
AX=-8*n*r*Du*Dv+8*m*(1-r)-(8*D+4)*v*v*(1+r*u*u)+4*Dv*(1-r)
BB1=v*(8*D+4)*(1-et)
zero('UB11 full X derivative',s.diff(F,X).subs(graph)-u*Z*AX-W*(u*u*BB1-8*m*v*(1+et*u*u)))
NN,rv=s.symbols('NN rv',nonzero=True)
zero('UB12 exact residue Schur expression',NN*rv**2*Dv-4*rv+2-(4-2*NN*rv*Dv)**2/(4*NN*Dv)-(2-4/(NN*Dv)))
K=s.symbols('K',nonzero=True)
qv=2*K*K*(1+v*v)/(1+K*K)
ev=2*(v*v/(K*K)+K*K)/(1+K*K)
zero('UB9 global remote determinant',qv*ev-4*v*v-4*(v*v-K*K)**2/(1+K*K)**2)

# Literal n=0 double-zero corner: derive the unnormalized original parent.
Hc=4*n+4*m+2; zeta=(4*m+2)/Hc; xc=4*m/Hc; ec=C/(C+A)
cg={X:xc*v*W-zeta*u*Z,Y:ec*u*W}
fg=s.Poly(s.expand(F.subs(cg)),W,Z)
fzz=s.cancel(fg.coeff_monomial(Z**2))
eps=4*al+4*n*(4*m+2)/Hc
bb=s.cancel((fzz-4*al-eps*u*u-4*n*zeta**2*u**4-4*D*v*v)/(u*u*v*v))
zero('UB double-zero fzz B_b residue',bb.subs({u:0,v:0,n:0,m:0,D:0,al:0})-2)
if not s.denom(bb).has(u,v):
    checks.append({'name':'UB double-zero fzz remainder is polynomial in testing angles','passed':True})
else: raise AssertionError(s.factor(bb))
fww0=s.cancel(fg.coeff_monomial(W**2).subs(u,0))
zero('UB double-zero W squared at u zero',fww0-4*m*(2*n+1)*v*v/(2*n+2*m+1)-(4*n+2)*xc**2*v**4)
fx=s.expand(s.diff(F,X).subs(cg))
zero('UB double-zero full X derivative at v zero',fx.subs(v,0)+8*n*zeta*u**3*Z)
fxz=s.cancel(s.diff(fx,Z)+8*n*zeta*u**3)
fxw0=s.cancel(s.diff(fx,W).subs(u,0))
if not s.denom(s.cancel(fxz/(u*v*v))).has(u,v):
    checks.append({'name':'UB double-zero remaining X-Z derivative has factor u v squared','passed':True})
else: raise AssertionError(s.factor(fxz))
if not s.denom(s.cancel(fxw0/(m*v**3))).has(u,v):
    checks.append({'name':'UB double-zero X-W derivative at u zero has factor m v cubed','passed':True})
else: raise AssertionError(s.factor(fxw0))

report={'scope':'Selected independent exact algebra; no automatic proof of inequalities or atlas coverage',
        'passed':sum(c['passed'] for c in checks),'total':len(checks),'checks':checks}
out=Path(__file__).with_suffix('.json')
out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
