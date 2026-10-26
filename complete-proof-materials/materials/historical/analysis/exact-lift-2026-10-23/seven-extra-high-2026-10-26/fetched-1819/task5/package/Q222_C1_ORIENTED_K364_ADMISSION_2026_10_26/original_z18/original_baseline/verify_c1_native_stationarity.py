#!/usr/bin/env python3
"""Reproduce 6 native C1 contact equations from the original SD2, entirely exact."""
from c1_native import *
from time import monotonic
import sys
s0=monotonic()
subs=dict(zip((u,v,W,X,Y,Z), CONTACT));
phis=phi_state();print('Built exact original SD2 source parent after cancellations, time',round(monotonic()-s0,2),flush=True)
# Evaluate each individual derivative contribution first; then collect in seven coquotients.
co=s.symbols('C2 CB C3 C5 C4 CA CD')
dictionary=dict(zip(COQ,co))
# avoid replace full coquotients in expanded objects: rebuild linear coefficient form from symbolic c's
from c1_native import c2,cb,c3,c5,c4,cal,cd
Du=1+u*u;Dv=1+v*v;bb=X+u*Z;qa=W+u*Y;ra=u*W-Y;r2=v*qa-bb;r5=ra-v*(Z-u*X)
cs=s.symbols('c2 cb c3 c5 c4 cal cd')
H2=4*T*cs[0];Hb=4*T*cs[1];H3=4*R*cs[2];H5=4*R*cs[3]
terms=[H2*(R/M)*r2**2,Hb*(R/t)*Dv*bb**2,Delta**2*H3*(T/a)*Dv*ra**2,
       Delta**2*H5*(T/J)*r5**2,4*R*T*Delta**2*eta**2*cs[4]*Du*Dv*Y**2,
       4*R*T*Delta**4*k**6*eta**2*cs[5]*Du*(v*Y-Z)**2,
       2*R*T*Delta**2*eta*cs[6]*Du*((v*Y-Z)**2-Dv*Z**2)]
# Source rational denominators M,J vanish by contact: keep factorizations in term gradients.
for variable in (W,X,Y,Z,u,v):
   items=[s.factor(s.diff(term,variable).subs(subs)) for term in terms]
   expr=s.Add(*items); expr=s.factor(expr)
   # insert coquotient formulas in a flat *linear* expression
   out=s.cancel(expr.subs(dict(zip(cs,COQ))))
   ok=out==0
   print(variable, 'EXACT_STATIONARITY',ok, 'time',round(monotonic()-s0,1),flush=True)
   if not ok:print('unresolved numerator',str(s.factor(s.together(out)))[:1500],flush=True);sys.exit(1)
print('PASS: all four full state derivatives and both angle derivatives exactly zero')
