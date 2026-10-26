"""All six general native coquotient contact balances, before scale assumptions."""
import json
from pathlib import Path
import sympy as s
x,eta,zeta,a,t=s.symbols('x eta zeta a t')
K=1+x; V=x*eta*zeta; L=1-a*eta; H=1+(1-a)*eta
Q=1-x*zeta+a*V; R=1+x*L; T=1+x-t*x*zeta
P4=1+x+eta*x*(a-x-1)-zeta*x*t+V*(-zeta*x*t*t+2*x*t+a*(1-t)+t)
Pal=1+x-a*eta*x+zeta*((1+x)-(2+x)*t)+V*(-zeta*t*t+a*(t-1)+t)
Pd=(1+x)**2-a*eta*x*x+zeta*x*((x+1)-(x+2)*t)+V*(-zeta*x*t*t+x*a*(t-1)+x*t+x+1)
c2=2*K*(1-t*zeta)*L*(1+V)*H*Q
cb=2*K*L*(1+(1-t)*zeta)*H*Q
c3=2*K*L*H*Q*(1-x*eta+t*V)
c5=2*K*L*H*Q*(1+V)
c4=2*Q*P4; ca=2*L*H*Pal; cd=4*L*Q*Pd
balances={
 'cb-c2-zeta_c3':cb-c2-zeta*c3,
 'c3-c5+x_eta_cb':c3-c5+x*eta*cb,
 'c2-(1-t_zeta)c5':c2-(1-t*zeta)*c5,
 'fourth_balance':2*x*c2+2*(1+V)*L*c4-2*x*x*zeta*Q*ca-Q*cd,
 'fifth_balance':x*L*c2+c5+(1+V)*L*L*c4+x*Q*Q*ca+s.Rational(1,2)*(x*zeta*L*L-2*L-eta)*cd,
 'sixth_balance':x*R*c2+t*x*V*cb+a*V*c3+T*c5+K*V*L*L*c4-x*x*zeta*K*Q*L*ca-s.Rational(1,2)*K*(Q*L+eta)*cd,
}
out={key:s.expand(val)==0 for key,val in balances.items()}
q,r=s.symbols('q r',nonzero=True)
mapped={eta:q/x,a:x*(1-r)/q}
pdq=s.cancel(Pd.subs(mapped,simultaneous=True))
P4q=s.cancel(P4.subs(mapped,simultaneous=True))
Hq=2+x*r+zeta*(K-(x+2)*t)+q*zeta*(-zeta*t*t+t+1)+x*zeta*(1-r)*(t-1)
B5=K*(1+q*zeta)*(r+q/x)-(K-x*zeta*t-x*(1-r)/q)*(q/x)*pdq
rhs=K*r*(1+q*zeta)+(1-r)*pdq+q*(zeta*t*(1+q*zeta)-(K-x*zeta*t)*Hq)
out['CE5_bounded_Pd']=s.factor(pdq-(1+q*zeta+x*Hq))==0
out['CE5_B5']=s.factor(B5-rhs)==0
rhs4=K*(1-q)+x*(1-r)-x*zeta*t+q*zeta*(-x*zeta*t*t+2*x*t+t)+x*zeta*(1-r)*(1-t)
out['CE7_P4']=s.factor(P4q-rhs4)==0
out['all_zero']=all(out.values());out['scope']='Polynomial identities in unconstrained x,eta,zeta,a,t; no asymptotic specialization'
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
