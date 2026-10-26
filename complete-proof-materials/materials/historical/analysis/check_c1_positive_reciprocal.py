"""Independent exact expansion of the positive-total C1 reciprocal scalar.

Rebuilds the parent from event amplitudes. No original checker is imported.
"""
from pathlib import Path
import sympy as s
import json,time
start=time.monotonic(); report={}
W,X,Y,Z,u,b,z,r,x,y,m,C,A,D,al,be=s.symbols('W X Y Z u b z r x y m C A D al be')
P,Q,V,T,h,cc,E,Fc,G,J=s.symbols('P Q V T h cc E Fc G J')
du,db=1+u*u,1+b*b
F=m*(W-b*X+u*Y-u*b*Z)**2+C*db*(u*W-Y)**2+A*du*db*Y**2
F+=(D+s.Rational(1,2))*(u*b*W+u*X-b*Y-Z)**2
F+=(al+s.Rational(1,2))*du*(Y-b*Z)**2+(be+s.Rational(1,2))*db*(X+u*Z)**2-du*db*Z**2/2
graph={W:1,Z:z,u:r*b*z,X:x*b*z*z,Y:y*b*z}
coefs={cc:A*y*y+(C+D+s.Rational(1,2))*(r-y)**2,
 E:be*(r+x)**2+x*x/2-2*D*r*x,G:A*r*r*y*y,J:(D+s.Rational(1,2))*r*r*x*x,
 Fc:A*r*r*y*y+(al+s.Rational(1,2))*r*r*(y-1)**2+be*(r+x)**2+m*(r*(y-1)-x)**2+(2*D+1)*r*x*(r-y)+r*x+x*x/2,
 h:A*y*y+C*(r-y)**2+al*(y-1)**2+y*y/2-r-2*D*(r-y)+2*m*(r*(y-1)-x)}
phi=m+D*T+h*V*T+cc*V*V*T+E*V*T*T+Fc*V*V*T*T+G*V**3*T*T+J*V*V*T**3
assert s.expand(F.subs(graph)-phi.subs(coefs).subs({V:b*b,T:z*z}))==0
report['CR9_to_CR11_all_scalar_terms']=True
balance={D:P*P*cc-P*Q*E+P**3*Q*G-P*P*Q*Q*J,
 m:P*Q*Q*E+P*P*Q*Q*Fc+P**3*Q*Q*G+2*P*P*Q**3*J,
 h:-2*P*cc-Q*E-2*P*Q*Fc-3*P*P*Q*G-2*P*Q*Q*J}
assert all(s.expand(g.subs({V:P,T:Q}).subs(balance))==0 for g in [phi,s.diff(phi,V),s.diff(phi,T)])
report['CR13_all_three_contact_identities']=True
R=z*(b*b-P); K=b*b*z*z-P*Q
M=b**4*z*z+(P/Q)*b*b*z**4+P*P*Q-3*P*b*b*z*z
balanced=(D/P**2)*R**2+(E*Q/P)*M+Fc*K*K+G*(P*K*K+b*b*z*R*K)+J*((Q*R-z*K)**2+2*Q*K*K)
assert s.cancel((phi.subs({V:b*b,T:z*z})-balanced).subs(balance))==0
report['CR15_full_balanced_scalar']=True
assert s.expand(2*Q*M+z*z*R*R-K*K-(K-z*R)**2-2*Q*R*R)==0
report['CR23_exact_zero_for_real_pivot_transfer']=True
assert s.cancel((E*Q/P)*M-(E*Q/P)*R*R-(E/P)*(K*K-z*R*K))==0
report['CR27_unbounded_source_ratio_identity']=True
bx,bz,pp,qq,ee=s.symbols('bx bz pp qq ee',positive=True)
six={b:bx**6,z:bz**6,P:pp**6,Q:qq**6,E:ee**6}
AA=bx**8*bz**4; BB=pp**2/qq**2*bx**4*bz**8; CC=pp**4*qq**2
sqrtA=bx**4*bz**2; sqrtB=pp/qq*bx**2*bz**4; sqrtC=pp**2*qq
theta=ee**6*qq**6/pp**6; sqrtTheta=ee**3*qq**3/pp**3; rho=pp**2/qq**2
assert s.cancel(M.subs(six)-(AA+BB+CC)*((AA-BB)**2+(AA-CC)**2+(BB-CC)**2)/2)==0
report['CR21_cubic_SOS']=True
first=bx**6*sqrtA*(AA-BB)+rho*bx**2*bz**4*sqrtA*(AA-CC)
second=ee**3*bz**6*(sqrtTheta*sqrtB*(BB-AA))+(ee**3*qq**2/pp**2)*bx**4*bz**2*(sqrtTheta*sqrtB*(BB-CC))
assert s.cancel((b*R).subs(six)-first)==0
assert s.cancel((E*b*z*(z*z-Q)).subs(six)-second)==0
report['CR22_both_complete_conductors']=True
assert s.cancel(K.subs(six)-rho**(-1)*(sqrtA**2*(BB-CC)+sqrtC**2*(AA-CC)))==0
report['CR25_infinitesimal_K_representation']=True
report['seconds']=time.monotonic()-start
report['scope']='Exact complete graph scalar and conductor identities; native inequalities and full normal estimates reviewed analytically.'
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
