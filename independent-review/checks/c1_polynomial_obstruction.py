"""Exact C1 signed-contact support obstruction; one prescribed rational source."""
from pathlib import Path
import json
import sympy as s

R=s.Rational
# Current SD1, specialized only after defining every coquotient.
k=f=g=s.S.One;a=R(7,8);t=R(3,5)
eta=f*f;zeta=g*g;kap=k*k;e=k*f/g;inv_e=1/e
Rc=1+kap*(1-a*eta);Tc=1+kap-kap*zeta*t;M=Rc-t;J=Tc-a
Vc=kap*eta*zeta;Lc=1-a*eta;Hc=1+(1-a)*eta;Qc=1-kap*zeta+a*Vc
P4=1+kap+eta*kap*(a-kap-1)-zeta*kap*t+Vc*(-zeta*kap*t*t+2*kap*t+a*(1-t)+t)
Pa=1+kap-a*eta*kap+zeta*((1+kap)-(2+kap)*t)+Vc*(-zeta*t*t+a*(t-1)+t)
Pd=(1+kap)**2-a*eta*kap**2+zeta*kap*((kap+1)-(kap+2)*t)+Vc*(-zeta*kap*t*t+kap*a*(t-1)+kap*t+kap+1)
c2=2*(1+kap)*(1-t*zeta)*Lc*(1+Vc)*Hc*Qc
cb=2*(1+kap)*Lc*(1+(1-t)*zeta)*Hc*Qc
c3=2*(1+kap)*Lc*Hc*Qc*(1-kap*eta+t*Vc)
c5=2*(1+kap)*Lc*Hc*Qc*(1+Vc)
c4=2*Qc*P4;cal=2*Lc*Hc*Pa;cd=4*Lc*Qc*Pd
dd=2*a*eta*cd/c3;mm=2*a*e*e*c2/(M*c3);CC=s.Integer(2)
AA=2*a*eta**2*c4/c3;al=2*a*inv_e**2*k**6*eta**2*cal/c3
DD=2*a*c5/(J*c3)-dd/2;be=2*a*e*e*cb/(t*c3)-dd/2
src=dict(zip(['d','lambda','n','m','C','A','D','alpha','beta'],[dd,0,0,mm,CC,AA,DD,al,be]))
src={key:s.sympify(value) for key,value in src.items()}
S=k/e;rc=1-a*eta;theta=t*zeta;ell=k*f*g
uc,vc=k,ell
psi=s.Matrix([S,k*(theta-1),k*S*rc,1])

u,v,W,X,Y,Z=s.symbols('u v W X Y Z')
d,n,m,C,A,D,alpha,beta=s.symbols('d n m C A D alpha beta')
Du,Dv=1+u*u,1+v*v
F=(4*n*Du*Dv*X*X+4*m*(v*W-X+u*v*Y-u*Z)**2
 +4*C*Dv*(u*W-Y)**2+4*A*Du*Dv*Y*Y
 +(4*D+2*d)*(u*W+u*v*X-Y-v*Z)**2
 +(4*alpha+2*d)*Du*(v*Y-Z)**2
 +(4*beta+2*d)*Dv*(X+u*Z)**2-2*d*Du*Dv*Z*Z)
poly=s.Poly(F,u,v,W,X,Y,Z)
checks={}
checks['forced_zero_coefficients']=poly.coeff_monomial(W*W)==0 and poly.coeff_monomial(u*u*X*X)==4*n
checks['invariant_character_diagonals']=[poly.coeff_monomial(v*v*X*X),poly.coeff_monomial(u*u*Y*Y),poly.coeff_monomial(u*u*v*v*Z*Z)]==[4*n+4*beta+2*d,4*A,4*beta]
checks['uv_character_diagonals']=[poly.coeff_monomial(u*u*v*v*W*W),poly.coeff_monomial(v*v*Y*Y),poly.coeff_monomial(Z*Z)]==[4*C,4*C+4*A+4*alpha+2*d,4*alpha]
ss={d:dd,n:0,m:mm,C:CC,A:AA,D:DD,alpha:al,beta:be}
ct={**ss,u:uc,v:vc,W:psi[0],X:psi[1],Y:psi[2],Z:1}
checks['all_original_retained_events_positive']=all(src[key]>0 for key in ['d','m','C','A','D','alpha','beta'])
checks['native_gaps_and_contact_amplitudes']=all(q>0 for q in [a,t,k,f,g,M,J,Qc,1-kap*eta*(1-theta),Rc,Tc,rc,1-rc,theta,1-theta]) and all(q!=0 for q in psi)
checks['actual_six_stationary_equations']=all(s.diff(F,q).subs(ct)==0 for q in [W,X,Y,Z,u,v])
for su in [-1,1]:
    for sv in [-1,1]:
        test={**ss,u:su*uc,v:sv*vc,W:psi[0],X:sv*psi[1],Y:su*psi[2],Z:su*sv}
        checks['signed_contact_'+str((su,sv))]=F.subs(test)==0

H=lambda ds:s.factor(2*sum(ds[i]*ds[j] for i in range(3) for j in range(i+1,3))-sum(q*q for q in ds))
inv_diag=[(4*be+2*dd)*(vc*psi[1])**2,4*AA*(uc*psi[2])**2,4*be*(uc*vc)**2]
uv_diag=[4*CC*(uc*vc*psi[0])**2,(4*CC+4*AA+4*al+2*dd)*(vc*psi[2])**2,4*al]
hinv=H(inv_diag);huv=H(uv_diag)
checks['strict_invariant_triangle_violation']=hinv<0
checks={key:bool(value) for key,value in checks.items()}
res={'checks':checks,'all_pass':all(checks.values()),'source_order':list(src),
 'source_C_equals_2':{key:str(value) for key,value in src.items()},
 'source_d_equals_1':{key:str(value/dd) for key,value in src.items()},
 'native_parameters':{'k':'1','f':'1','g':'1','a':str(a),'t':str(t)},
 'contact':[str(uc),str(vc)]+[str(q) for q in psi],
 'invariant_scaled_diagonal':[str(q) for q in inv_diag],
 'uv_scaled_diagonal':[str(q) for q in uv_diag],
 'invariant_Heron':str(hinv),'uv_Heron':str(huv),
 'scope':'One strict native C1 source is not an ordinary polynomial SOS. General necessary triangle inequalities, not a dense generic exclusion theorem.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(res,indent=2)+'\n',encoding='utf-8')
print(json.dumps(res,indent=2))
if not res['all_pass']:raise SystemExit(1)
