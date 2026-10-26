"""Bounded independent single-seed source/parent-budget validation."""
from pathlib import Path
import json,time
import sympy as s
START=time.perf_counter();Q=s.Rational
u,v,W,X,Y,Z=s.symbols('u v W X Y Z')
a=b=c=Q(1);t=Q(3,4);r=Q(1,4);z=c*c;Ac=1+z*r
h=b-a*a*t;g=(1+z)*b-(1-r+a*a*z*t);j=1+z*r-t
Delta=1+(1+z)*b+z*(1+z)*b*b+2*z*(1+z)*b*r+z*z*r*r+2*z*r+a*a*z*((1+z)*b+1+z*r*(2-r))
N=(1+z)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z)*(1-r)*t
Nn=(b+r)*Ac*(Ac+a*a*z*t)-r*Delta
gaps={'h':h,'g':g,'j':j,'1-c2h':1-z*h,'b-a2c2r':b-a*a*z*r,'N_star':N,'N_n':Nn,'Delta':Delta}
assert all(q>0 for q in gaps.values()) and all(q>0 for q in (a,b,c,t,r)) and r<1
n=Nn/(2*r*Delta)
m=(1+a*a*z)*(1-r)*Ac*(Ac+a*a*z*t)/(2*g*Delta)
C=Ac*N/(2*a*a*t*Delta)
A=a*a*Ac*(b+r)*(Ac+a*a*z*t)/(2*h*Delta)
D=(1+z)*(b+r)*(b-a*a*z*r)*(Ac+a*a*z*t)/(2*a*a*j*Delta)
al=z*z*(b+r)*N/(2*(1-z*h)*Delta)
weights=[n,m,C,A,D,al]
assert weights==[Q(1,8),Q(3,8),Q(1,6),Q(5,8),Q(3,8),Q(1,6)]
contact=[c,a*c,b/a,-c*r,c*(b/a-a*t),Q(1)]
assert contact==[1,1,1,-Q(1,4),Q(1,4),1]
Du=1+u*u;Dv=1+v*v
# Direct SG2 full-state parent: independently expanded from nine source/event
# terms, not by importing the packet rows or projection checker.
F=4*n*Du*Dv*X**2+4*m*(v*W-X+u*v*Y-u*Z)**2+4*C*Dv*(u*W-Y)**2+4*A*Du*Dv*Y**2+(4*D+2)*(u*W+u*v*X-Y-v*Z)**2+(4*al+2)*Du*(v*Y-Z)**2+2*Dv*(X+u*Z)**2-2*Du*Dv*Z**2
P=s.expand(Dv*F)
physical={W:1,X:-Q(1,4),Y:Q(1,4),Z:1}
Pu=s.Poly(P.subs(physical).subs(v,0),u)
u_capacity=Pu.coeff_monomial(u**2)
origin_capacity=P.subs({u:0,v:0,W:0,X:0,Y:0,Z:1})
assert u_capacity==Q(409,48) and origin_capacity==Q(2,3)
assert u_capacity+origin_capacity==Q(147,16)
assert Q(147,16)-Q(57818093029,6288328800)==-Q(44072179,6288328800)
report={'result':'BOUNDED SEED SOURCE AND TWO DIRECT-PARENT CAPACITIES PASSED','gaps':{k:str(v) for k,v in gaps.items()},'native_weights':[str(x) for x in weights],'native_contact':[str(x) for x in contact],'u_squared_capacity_at_v0_physical_state':str(u_capacity),'origin_capacity_Z_only_state':str(origin_capacity),'combined_capacity':str(u_capacity+origin_capacity),'elapsed_seconds':round(time.perf_counter()-START,6),'scope':'One genuine native source and two necessary budgets only; no module/projection replay, sweep, PSD or K364 admission claim.'}
out=Path(__file__).resolve().parent
(out/'bounded_parent_budget_results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
