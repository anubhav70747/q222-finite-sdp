"""Exact counterexample to the illustrative w^2 row <=4 bound only.

No counterexample to the128 coefficient bound or local certificate follows.
"""
from fractions import Fraction as F
import json
from pathlib import Path

t=F(1,2**30)
x=sigma=h=t
b=z=1-t**4
r=b*z
theta=1-h
e=x*h/sigma
E=1-e
K=1+x
J=1+x*h-sigma*(1-r)
M=1+x*r-x*z*theta
P=(sigma+K+sigma*z*(1+b*x+2*x*h)
   +sigma*b*x*x*h*z*z+x*x*h*theta*z)
p=(z*(1+x*h)*E+(1-r)*(1+x*h*z))/(2*P)
D=(K*(1+sigma*z)*(1+sigma*r)-J*P)/(2*J*P)
beta=(K*(1+sigma*r)*(1+x*h*z)-theta*P)/(2*theta*P)
lam=(8*D*h+4*beta*theta*theta+2*h*h)/sigma-4*b*p*sigma*z
row=2*z*(1+2*D)*(1+b)
assert 0<x<=F(1,2**30) and 0<sigma<=F(1,2**30)
assert 0<b<1 and 0<z<=1 and 0<r<1 and 0<h<1
assert min(E,J,M,D,beta,lam)>0
assert 4<row<5
result={
 'scope':'Counterexample to uniform_bounds.tex illustrative w^2 row <=4; bound5 suffices and128 is unaffected.',
 'source_parameters':{k:str(v) for k,v in {'x':x,'sigma':sigma,'h':h,'b':b,'z':z}.items()},
 'D':str(D),'beta':str(beta),'lambda':str(lam),
 'w2_row':str(row),'w2_row_minus4':str(row-4),
 'all_native_gaps_and_positive_lambda':True,
 'row_strictly_between4_and5':True}
(Path(__file__).parent/'pro669_w2_bound_counterexample.json').write_text(json.dumps(result,indent=2)+'\n')
print('Exact admissible source has4 < w^2 coefficient-row norm <5;128 bound is unaffected.')
