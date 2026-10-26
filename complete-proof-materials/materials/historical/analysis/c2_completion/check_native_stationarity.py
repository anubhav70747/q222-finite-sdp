"""Rebuild the entire parent and verify all six WE/R9 candidate equations."""
import json
import time
from pathlib import Path
import sympy as S

started=time.monotonic()
a,B,t,h,r=S.symbols('a B t h r')
u,v,W,X,Y,Z,k,l=S.symbols('u v W X Y Z k l')
n,m,C,A,D,alpha=S.symbols('n m C A D alpha')
ca=1+a*B*r
Sc=1+B*(1-h+a*r)
g=1+a*(h*B+r-1)
j=ca-(1-h)/t
Delta=a+(1+a*B)+(1+a*B)*B+2*a*B*(1+a*B)*r+a**3*B**2*r**2+2*a**2*B*r+t*B*(1+a+a*B+a*a*B*r*(2-r))
Mass=(1+a*B)*(1+a*r)*(h+t)-a*(1+t*B)*(1-r)*(1-h)
weights={
 n:((1+a*r)*ca*Sc-r*Delta)/(2*r*Delta),
 m:a*a*(1+t*B)*(1-r)*ca*Sc/(2*g*Delta),
 C:ca*Mass/(2*(1-h)*Delta),
 A:t*ca*(1+a*r)*Sc/(2*h*Delta),
 D:(1+a*B)*(1+a*r)*(1-a*t*B*r)*Sc/(2*t*j*Delta),
 alpha:B*B*(1+a*r)*Mass/(2*(1-h*B)*Delta)}
Du=1+u*u; Dv=1+v*v
M=v*W-X+u*v*Y-u*Z
J=u*W+u*v*X-Y-v*Z
P=n*Du*Dv*X*X+m*M*M+C*Dv*(u*W-Y)**2+A*Du*Dv*Y*Y+(D+S.Rational(1,2))*J*J+(alpha+S.Rational(1,2))*Du*(v*Y-Z)**2+S.Rational(1,2)*Dv*(X+u*Z)**2-S.Rational(1,2)*Du*Dv*Z*Z

# Multiply the actual Z=1 contact state by kl, keeping k^2=aB,l^2=tB.
# This avoids introducing radicals into any identity under inspection.
contact={u:k,v:l,W:B,X:-r*a*B*l,Y:h*B*k,Z:k*l}
results={}
for variable in (W,X,Y,Z,u,v):
    expression=S.diff(P,variable).xreplace(contact).subs(weights,simultaneous=True)
    num=S.together(expression).as_numer_denom()[0]
    poly=S.Poly(S.expand(num),k,l)
    remainder=0
    for (ik,il),coefficient in poly.terms():
        remainder+=coefficient*(a*B)**(ik//2)*(t*B)**(il//2)*k**(ik%2)*l**(il%2)
    ok=S.expand(remainder)==0
    results['dP_d'+str(variable)]=ok
    if not ok:
        raise AssertionError(str(variable)+': '+str(S.factor(remainder)))
results['statement']='All six complete-parent stationary equations vanish identically for the exact six printed coefficients, modulo k^2=aB,l^2=tB; no positivity or numerical specialization is used.'
results['elapsed_seconds']=round(time.monotonic()-started,3)
Path(__file__).with_suffix('.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
print(json.dumps(results,indent=2))
