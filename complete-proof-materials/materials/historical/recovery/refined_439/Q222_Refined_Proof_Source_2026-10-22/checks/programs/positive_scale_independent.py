"""Fresh targeted algebra for the positive-scale audit; no archive engine.

This corroborates the accompanying hand proof.  In particular, exact identities
do not certify field bounds, positivity, contact entry, or chart coverage.
"""
import sys
import sympy as s
import json
from pathlib import Path

records=[]
def check(name, expr):
    value=s.cancel(s.expand(expr))
    if value != 0:
        raise AssertionError((name, s.factor(value)))
    records.append(name)

u,v,W,X,Y,Z=s.symbols('u v W X Y Z', real=True)
n,m,C,A,D,al,be=s.symbols('n m C A D al be', real=True)
r,x,y,P,Q,V,T=s.symbols('r x y P Q V T', real=True)
b,z,a,w=s.symbols('b z a w', real=True)
Du=1+u*u; Dv=1+v*v
F=(n*Du*Dv*X**2 +m*(v*W-X+u*v*Y-u*Z)**2+C*Dv*(u*W-Y)**2
   +A*Du*Dv*Y**2+(D+s.Rational(1,2))*(u*W+u*v*X-Y-v*Z)**2
   +(al+s.Rational(1,2))*Du*(v*Y-Z)**2
   +(be+s.Rational(1,2))*Dv*(X+u*Z)**2-s.Rational(1,2)*Du*Dv*Z**2)
F=s.expand(F)
Fr=s.expand(b*b*F.subs(v,1/b))
assert s.denom(s.cancel(Fr))==1

# Weak-alpha full scalar, including every high-order monomial.
fg=s.expand(F.subs({W:1,Z:z,u:r*v*z,X:x*v*z*z,Y:y*v*z}, simultaneous=True))
cc=A*y*y+C*(r-y)**2+(al+s.Rational(1,2))*y*y+2*m*r*y
ee=al*r*r+(m+be+s.Rational(1,2))*(r+x)**2+n*x*x
gg=r*r*y*y*(A+al+m+s.Rational(1,2))
jj=n*r*r*x*x
kk=(D+n+s.Rational(1,2))*r*r*x*x
ff=(A*r*r*y*y+(2*D+1)*r*x*(r-y-1)-(2*al+1)*r*r*y
    -2*m*r*y*(r+x)+be*(r+x)**2+n*x*x+r*x+x*x/2)
hh=s.Poly(fg,v,z).coeff_monomial(v*v*z*z)
expected=m*v*v+al*z*z+hh*v*v*z*z+cc*v**4*z*z+ee*v*v*z**4+ff*v**4*z**4+gg*v**6*z**4+jj*v**4*z**6+kk*v**6*z**6
check('WA14: entire graph scalar', fg-expected)
normal={j:s.expand(s.diff(F,j).subs({W:1,Z:z,u:r*v*z,X:x*v*z*z,Y:y*v*z},simultaneous=True)/2) for j in [u,X,Y]}
au=s.cancel(normal[u]/(v*z)); ay=s.cancel(normal[Y]/(v*z)); bx=s.cancel((normal[X]/v+m)/(z*z))
check('WA19: X-normal c2',s.Poly(bx,v,z).coeff_monomial(v*v*z*z)-n*r*r*x)
check('WA19: X-normal c3',s.Poly(bx,v,z).coeff_monomial(v**4*z*z)-(D+n+s.Rational(1,2))*r*r*x)
check('WN11: u-normal a01',s.Poly(au,v,z).coeff_monomial(z*z)-(al*r+(m+be+s.Rational(1,2))*(r+x)))
check('WN11: u-normal a12',s.Poly(au,v,z).coeff_monomial(v*v*z**4)-n*r*x*x)
c0=s.Poly(bx,v,z).coeff_monomial(1)
check('WN4: exact positive E completion',ee-(al*r*r+c0*c0/(n+m+be+s.Rational(1,2))+n*(m+be+s.Rational(1,2))*r*r/(n+m+be+s.Rational(1,2))))

# C1 reciprocal full graph and all special first-normal coefficients.
fc=Fr.subs(n,0)
cg=s.expand(fc.subs({W:1,Z:z,u:r*b*z,X:x*b*z*z,Y:y*b*z},simultaneous=True))
CC=A*y*y+(C+D+s.Rational(1,2))*(r-y)**2
EE=be*(r+x)**2+x*x/2-2*D*r*x
GG=A*r*r*y*y; JJ=(D+s.Rational(1,2))*r*r*x*x
FF=(A*r*r*y*y+(al+s.Rational(1,2))*r*r*(y-1)**2+be*(r+x)**2
    +m*(r*(y-1)-x)**2+(2*D+1)*r*x*(r-y)+r*x+x*x/2)
HH=A*y*y+C*(r-y)**2+al*(y-1)**2+y*y/2-r-2*D*(r-y)+2*m*(r*(y-1)-x)
check('CR11--12: entire graph scalar',cg-(m+D*z*z+HH*b*b*z*z+CC*b**4*z*z+EE*b*b*z**4+FF*b**4*z**4+GG*b**6*z**4+JJ*b**4*z**6))
cnormal={j:s.expand(s.diff(fc,j).subs({W:1,Z:z,u:r*b*z,X:x*b*z*z,Y:y*b*z},simultaneous=True)/2) for j in [u,X,Y]}
cbx=s.cancel((cnormal[X]/b+m)/(z*z)); cau=s.cancel(cnormal[u]/(b*z))
check('CR17: entire X bracket',cbx-(-D*r+be*(r+x)+x/2+((D+s.Rational(1,2))*r*(r-y)+(be+s.Rational(1,2))*(r+x)+m*(x-r*(y-1)))*b*b+(D+s.Rational(1,2))*r*r*x*b*b*z*z))
check('CR18: u-normal a01',s.Poly(cau,b,z).coeff_monomial(z*z)-(be*(r+x)-D*x))
check('CR18: u-normal a12',s.Poly(cau,b,z).coeff_monomial(b*b*z**4)-(D+s.Rational(1,2))*r*x*x)

# Weighted finite-edge expansion: no test infinitesimal is divided.
d=s.symbols('d',real=True)
Fraw=4*F # d=1 suffices for polynomial divisibility and structural jets.
edge=s.expand(Fraw.subs(W,(1+X-u*v*Y+u*Z)/v))
ax,ay0,az=s.symbols('ax ay az',real=True)
edge_graph=s.expand(edge.subs({X:u*u*ax,Y:u*ay0,Z:u*az},simultaneous=True))
Kedge=s.cancel((edge_graph-4*m)/(u*u))
assert s.denom(Kedge).free_symbols <= {v}
check('finite edge: K_x divisible by U',s.diff(Kedge,ax).subs(u,0))
Ex=s.cancel(s.diff(Kedge,ax)/(u*u))
check('finite edge: divided x pivot',s.diff(Ex,ax).subs(u,0)-4*(1+2*n+2*be)*(1+v*v))
for j,mult,center in [(X,u*u,ax),(Y,u,ay0),(Z,u,az)]:
    check('finite edge: full linear coefficient '+str(j),s.diff(edge,j).subs({X:u*u*ax,Y:u*ay0,Z:u*az},simultaneous=True)-s.diff(Kedge,center)*u*u/mult)

# Finite-ratio full tail in the undivided ideal (N1,N2,N3,v*zeta)^2.
corner=s.expand(Fraw.subs(Z,1+v*Y))
xc,yc,tc,wc,zeta,N1,N2,N3=s.symbols('xc yc tc wc zeta N1 N2 N3',real=True)
cgraph=s.expand(corner.subs({X:v*xc,Y:v*yc,u:v*tc,W:wc},simultaneous=True))
Kcorner=s.cancel((cgraph-4*al)/(v*v))
assert s.denom(Kcorner)==1
expanded=s.expand(corner.subs({X:v*xc+N1,Y:v*yc+N2,u:v*tc+N3,W:wc+zeta},simultaneous=True))
lin=(v*s.diff(Kcorner,xc)*N1+v*s.diff(Kcorner,yc)*N2+v*s.diff(Kcorner,tc)*N3+v*v*s.diff(Kcorner,wc)*zeta)
tail=s.Poly(s.expand(expanded-cgraph-lin),N1,N2,N3,zeta,v)
for powers,coeff in tail.terms():
    hdeg=sum(powers[:3]); zd=powers[3]; vd=powers[4]
    # Allocate at most two zeta's, each requiring one v. The rest is coefficient.
    if hdeg>=2 or (hdeg>=1 and zd>=1 and vd>=1) or (zd>=2 and vd>=2): continue
    raise AssertionError(('finite ratio: unallocated full tail',powers,coeff))
records.append('finite ratio: all tail monomials in the undivided square ideal')

# Source-only radical substitutions make the displayed conductor identities
# ordinary polynomial equalities. These tests include literal zero variables.
t,ss,h,c=s.symbols('t ss h c',real=True)
aa=t**3; ww=ss**3; ep=h**6; tau=c**6; La=tau*ep
UU=aa**2; TT=ww**2
L1=t**4*ss**2; L2=c*c*t*t*ss**4; L3=h**6*c**4
f=[t*t*ss,c*t*ss*ss,h**3*c*c]
rr=[]
for q0 in f: rr.extend([q0*(L1-L2),q0*(L1-L3)])
A1=aa*(TT-ep); K1=UU*TT-La*ep
MM=UU*UU*TT+tau*UU*TT*TT+La*La*ep-3*La*UU*TT
check('WN24: complete cubic Gram',MM-sum(rr[2*i]**2-rr[2*i]*rr[2*i+1]+rr[2*i+1]**2 for i in range(3)))
check('WN25.1',aa*ww*(UU-La)-(t**3*rr[0]+c*c*t*ss*ss*rr[1]))
check('WN25.2',aa*aa*ww*(TT-ep)-(-h**6*rr[0]+ss**6*rr[1]-h**3*ss**3*rr[4]))
check('WN25.3',tau*aa*ww*(TT-ep)-(c*c*t*ss*ss*(rr[1]-rr[0])-c**3*ss**3*rr[2]))
check('WN25.4',tau*h**3*A1-(t**3*rr[5]-t*(t*t+c*c*ss*ss)*rr[4]))
check('WN26',La*(UU-La)*K1-(La*MM-(tau*h**3*A1)**2))
check('WN27',c*c*K1-(f[0]*(rr[1]-rr[0])+f[2]*rr[5]))

# Whole X-normal quadratic coefficient retained in the thin branch.
Hx=n+m+be+s.Rational(1,2)
DX=(n+be+s.Rational(1,2))*v*v+n*u*u+(n+D+s.Rational(1,2))*u*u*v*v
check('WN37: complete X mass',s.diff(F,X,2)/2-Hx-DX)

# Exact projective core formulas: these keep arbitrary states.
core=F.subs({n:0,m:0,C:0,A:0,D:0,al:0,be:0})*2
sv=s.symbols('sv',real=True)
Vx,Rx,rho,T0=s.symbols('Vx Rx rho T0',real=True)
alice=s.cancel(sv*sv*core.subs({u:1/sv,X:Vx/v,Z:Rx-sv*Vx/v},simultaneous=True))
rc20=(W-sv*Y+(1+sv*sv)*Vx-sv*v*Rx)**2+(1+v*v)*Rx**2+(1+sv*sv)*(v*v*Y*Y-2*Y*(v*Rx-sv*Vx)-(v*Rx-sv*Vx)**2)
check('RC20: entire reciprocal Alice core',alice-rc20)
bob=s.cancel(b*b*core.subs({v:1/b,W:(rho+Y)/u,Z:T0/b},simultaneous=True))
rc22=(1+u*u+b*b)*X*X+2*u*b*X*T0+b*b*rho*rho+2*u*b*rho*X-2*rho*T0-T0*T0+(1+u*u)*(Y-T0)**2
check('RC22: entire reciprocal Bob core',bob-rc22)
both=s.cancel(sv*sv*b*b*core.subs({u:1/sv,v:1/b,Z:T0/b},simultaneous=True))
rc24=(b*W+X-sv*b*Y)**2+(1+sv*sv)*(Y-T0)**2+sv*sv*(1+b*b)*X*X-2*sv*W*T0+2*sv*sv*Y*T0+2*sv*b*X*T0-sv*sv*T0*T0
check('RC24: entire doubly reciprocal core',both-rc24)

# Independently derive RC4's two soft pivots directly from the actual contact
# derivatives; do not substitute the archived native coefficient formulas.
k,p,S,R,mu,delta=s.symbols('k p S R mu delta',real=True)
rcgraph=s.expand(Fr.subs(be,0).subs({u:k,b:p,W:S*p/k,X:-R*k*p*p,Y:y*p,Z:1},simultaneous=True))
opens={m:mu*k**4,D:delta*p**4}
mu_eq=s.cancel((s.diff(rcgraph,k)/(2*k*p*p)).subs(opens))
delta_eq=s.cancel(((p*s.diff(rcgraph,p)/2-rcgraph)/(p**4)).subs(opens))
assert not s.denom(mu_eq).has(k,p)
assert not s.denom(delta_eq).has(k,p)
mu_res=mu_eq.subs({k:0,p:0})
delta_res=delta_eq.subs({k:0,p:0})
critical={C:1/(2*(S-y)),A:(1+S-y)/(2*y),al:S/(2*(1-y)),n:(1+S-y-R)/(2*R)}
check('RC4: mu soft pivot from full angle derivative',mu_res.subs(critical)-(-mu*S*S+(1+S-y)/2))
check('RC4: delta soft pivot from full contact/angle derivative',delta_res.subs(critical)-(-delta+S*(1+S-y)/2))
check('RC5: critical residue identity',(C*(al*(A-s.Rational(1,2))-s.Rational(1,4))-(A+al+s.Rational(1,2))/4).subs(critical))

# These are the actual normalized full state jets, rather than a scalar test.
rcstate={j:s.cancel(s.diff(Fr.subs(be,0),j).subs({u:k,b:p,W:S*p/k,X:-R*k*p*p,Y:y*p,Z:1},simultaneous=True)/2) for j in [W,X,Y,Z]}
scales={W:k*p,X:k*p*p,Y:p,Z:p*p}
for j in [W,X,Y,Z]:
    normalized=s.cancel((rcstate[j]/scales[j]).subs(opens))
    assert not s.denom(normalized).has(k,p)
    check('RC4: normalized full '+str(j)+'-state equation at residue',normalized.subs({k:0,p:0}).subs(critical))

out={'status':'all exact targeted identities vanish','count':len(records),'identities':records,
     'limit':'These checks corroborate algebra only; the hand audit establishes boundedness, positivity, source entry, critical centers, and coverage.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
