"""Factored exact assembly of all purchases; no large parent expansion.

The parent/core interface is independently checked in
../full_r1_enlarged_gram.py. Each payment here is reduced in the small
normalized circuit and its coefficient vector is compared with the
exported loss polynomials. All computations are symbolic.
"""
from pathlib import Path
import json
import sympy as s
HERE=Path(__file__).resolve().parent
U,V,w,Z=s.symbols('U V w Z')
x,y,fc,qc,rc,Me,Sigma,theta,phi,j,T,alpha,c,kappa,gc,gamma,B=s.symbols(
    'x y fc qc rc Me Sigma theta phi j T alpha c kappa gc gamma B')
ss,vv=s.symbols('s vv')
P=(U*U+V*V)*w*w+(1+U*U*V*V)*Z*Z-4*U*V*w*Z
Q=(U*w-V*Z)**2;R=(Z-U*V*w)**2;Bs=(V*w-U*Z)**2
closed=(w-U*V*Z)**2
A27=(1+U*U)*P-R
B29=2*(P+V*V*Q+U*U*Bs)-R
E18=(1+V*V)*P-V*V*Bs
W=4*P+4*U*U*(Q+R)-V*V*(U*U-1)**2*w*w
RA=W+4*U*U*A27
RB=8*U*U*P+8*U**4*Q+(8*U**4+2*(U*U-1)**2)*Bs-V*V*(U*U-1)**2*w*w
checks={}
def check(name,expr):
    num=s.cancel(expr).as_numer_denom()[0]
    assert s.expand(num)==0,name
    checks[name]=True
def harmonic(expr):
    return expr.subs(theta,Me/Sigma).subs(Sigma,Me+fc/y)

rPsi=(1-phi)*4*x*x*y/fc+phi*(8*theta*x*y*ss/fc+4*(x*x+x*ss+ss*ss)/Sigma)
rQ=(1-phi)*4*x*y*ss/qc+phi*(8*theta*y*ss*ss+4*(1-theta)*x*y*ss)/qc
rE=phi*(8*ss*ss+2*(ss-x)**2)/Sigma
rR=(1-phi)*4*x*y*ss/rc
sub={ss:x*U*U,vv:y*V*V}
radpay=(fc*rPsi*P+qc*rQ*Q+Me*rE*y*Bs+rc*rR*R).subs(sub)
radrem=x*x*y*((1-phi)*W+phi*(1-theta)*RA+phi*theta*RB)
check('smooth_radial_purchase',harmonic(radpay-x*x*y*V*V*(U*U-1)**2*w*w-radrem))

Abar=V*(U*U-1)*w;Bbar=U*(Z-U*V*w)
check('whole_Jalpha_linear_relation',Abar+Bbar+(V*w-U*Z))
check('whole_Jalpha_purchase',
      2*(1-phi)*(Abar*Abar+Bbar*Bbar)+phi*(1+U*U)*P-Bs
      -(1-phi)*(Abar-Bbar)**2-phi*(R+U*U*Q+U*U*closed))

h=[phi*((1+theta)/fc+ss/(Sigma*y*x)),
   2*phi*theta*vv/(qc*y),2*phi*ss/(Sigma*y*x),(1-phi)/rc]
hpay=(fc*h[0]*P+qc*h[1]*Q+Me*h[2]*y*Bs+rc*h[3]*R).subs(sub)
check('hybrid_Rstar_purchase',harmonic(hpay-R-phi*((1-theta)*A27+theta*B29)))
check('O_majorization_positive_square',2*(alpha*alpha+c*c*ss*ss)-(alpha-c*ss)**2-(alpha+c*ss)**2)

# Here gamma*Sigma=p^2. No current coefficient is inverted.
hp1=gamma*fc*(y*vv+vv*vv)*P+gamma*Me*vv*vv*y*Bs
hp2=gamma*fc*T*ss*(y+vv)*P+gamma*Me*T*ss*vv*y*Bs
check('helper_H1',((hp1-gamma*Sigma*vv*vv*y*Bs).subs(sub)
      -gamma*fc*y*y*V*V*E18).subs(Sigma,Me+fc/y))
check('helper_H2',((hp2-gamma*Sigma*T*ss*vv*y*Bs).subs(sub)
      -gamma*T*fc*y*x*U*U*E18).subs(Sigma,Me+fc/y))

# Assemble coefficient vectors with abstract pivot squares, avoiding a
# simultaneous expansion of source denominators and state monomials.
PP,QQ,EE,RR,Ru2,J2,O2,H12,H22=s.symbols('PP QQ EE RR Ru2 J2 O2 H12 H22')
radvec=[rPsi,rQ,rE,rR]
beta=16+32*(1-phi)*(alpha/B)**2
Cbar=2*j*T*T*vv*(alpha*alpha+c*c*ss*ss)
expected=[beta*q+16*Cbar*z+2*kappa*z for q,z in zip(radvec,h)]
expected[0]+=16*phi*j*j*alpha*alpha*y*(1+ss/x)/fc+gc/fc
expected[3]+=32*(1-phi)*alpha*alpha*j*T*ss/rc
expected[1]+=2*kappa*vv/(qc*y)
expected[0]+=8*gamma*(y*vv+vv*vv+T*ss*(y+vv))
expected[2]+=8*gamma*(vv*vv+T*ss*vv)
catalog=json.loads((HERE/'root_catalog.json').read_text())
local={str(a):a for a in [x,y,fc,qc,rc,Me,Sigma,theta,phi,j,T,alpha,c,kappa,gc,gamma,B]}
uu,v=s.symbols('u v');local.update(u=uu,v=v)
for name,ex in zip(['Psi','Q','E','R'],expected):
    actual=s.sympify(catalog['loss_polynomials'][name],locals=local).subs({uu*uu:ss,v*v:vv})
    check('exported_loss_'+name,actual-ex)

L=lambda a:fc*a[0]*PP+qc*a[1]*QQ+Me*a[2]*EE+rc*a[3]*RR
Jother=phi*j*j*alpha*alpha*y*(1+ss/x)*PP+2*(1-phi)*alpha*alpha*j*T*ss*RR
Jpay=2*(1-phi)*(alpha/B)**2*Ru2+Jother
Hpay=gamma*fc*(y*vv+vv*vv+T*ss*(y+vv))*PP+gamma*Me*(vv*vv+T*ss*vv)*EE
paid=beta*(L(radvec)-Ru2)+16*(Jpay-J2)+16*(Cbar*L(h)-O2)
paid+=2*kappa*(L(h)-RR)+8*(Hpay-H12-H22)
residual=s.Rational(3,4)*fc*PP+s.Rational(1,2)*(qc*QQ+Me*EE+rc*RR)
residual-=gc*PP+2*kappa*vv/y*QQ+2*kappa*RR+16*Ru2+16*J2+16*O2+8*H12+8*H22
remaining=(s.Rational(3,4)-expected[0])*fc*PP
remaining+=(s.Rational(1,2)-expected[1])*qc*QQ+(s.Rational(1,2)-expected[2])*Me*EE+(s.Rational(1,2)-expected[3])*rc*RR
check('entire_scalar_residual_assembly',residual-remaining-paid)

(HERE/'certificate_identity_checks.json').write_text(json.dumps({
 'checks':checks,'method':'Factored exact polynomial identities and independent payment coefficient aggregation.',
 'parent_interface':'../normal_core_checks.json: entire_RFN_parent_equals_core_plus_retained_plus_scalar_residual',
 'no_optimizer_or_numerical_sampling':True},indent=2)+'\n',encoding='utf8')
print('All',len(checks),'factored payment, export and full-residual identities pass exactly.')
