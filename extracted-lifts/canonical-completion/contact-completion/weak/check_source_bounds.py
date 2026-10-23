"""Exact coefficient-envelope/cancellation checks; no point sampling or solver."""
from pathlib import Path
from fractions import Fraction
import json
import sympy as sp
root=Path(__file__).parent
import runpy, contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    ns=runpy.run_path(str(root/'received_generator.py'))
B,T,j,p,H=[ns[k] for k in ["B","T","j","p","H"]]
ss=sp.symbols("ss")
delta=sp.Rational(1,1024)
def envelope(expr):
    poly=sp.Poly(sp.cancel(expr.subs(p,j*ss)),B,T,j,H,ss)
    return sum(abs(c)*delta**(eb+et+ej)
               for (eb,et,ej,eh,es),c in poly.terms())
items={"P_minus_1":ns["P"]-1,"L_minus_1":ns["L"]-1,
       "FF_minus_1":ns["FF"]-1,"F1":ns["F1"],"F2_div_B":ns["F2"]/B}
env={name:envelope(expr) for name,expr in items.items()}
thresholds={"P_minus_1":sp.Rational(1,128),"L_minus_1":sp.Rational(1,128),
            "FF_minus_1":sp.Rational(1,128),"F1":4,"F2_div_B":2}
assert all(v<thresholds[k] for k,v in env.items())
n,mu,C,A,D,alpha=[ns[k] for k in ["n","mu","C","A","D","alpha"]]
half=sp.Rational(1,2); x=j*B; y=T*B
Kd=C*H-(D+half)*(1-H)
Kc=ns["Kc"]; src=ns["src"]
identities=[
 T*(Kd+A*H+half)-A,
 j*Kd+mu+x*T*C*H,
 j*Kd+x*T*C*H+T*(C*H+A*H+(alpha+half)*H)-A-(alpha+half)*(1-p)+alpha/B,
 n*p-mu+p/2+(-mu-(alpha+half)*(1-p))*x+j*(alpha-n)+y*Kc]
zeros=[]
for expr in identities:
    # Cancel each term's known source factors before combining numerators.
    terms=[sp.cancel(term.subs(src,simultaneous=True))
           for term in sp.Add.make_args(sp.expand(expr))]
    numerator=sp.together(sum(terms)).as_numer_denom()[0]
    answer=sp.expand(numerator)==0
    zeros.append(answer)
assert all(zeros),zeros
out={"scope":"Exact envelopes at delta=1/1024 and four rational source identities",
     "record_count":len(ns["rec"]),
     "source_envelopes":{k:str(v) for k,v in env.items()},
     "all_envelopes_below_stated_thresholds":True,
     "four_source_cancellation_numerators_zero":zeros}
(root/"weak_recovery_source_envelopes_check.json").write_text(
    json.dumps(out,indent=2)+"\n",encoding="utf-8")
print(json.dumps(out,indent=2))
