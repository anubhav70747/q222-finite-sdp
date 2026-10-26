#!/usr/bin/env python3
"""Exact order-at-f=0 inequalities for a 1-parameter strict native C1 obstruction.
No numeric samples or solvers. The f->0 limit is NOT itself a native contact.
"""
import sys
from pathlib import Path
import sympy as S
BASE=Path(__file__).resolve().parent/'original_baseline'
if not (BASE/'c1_native.py').exists():
    BASE=Path(__file__).resolve().parents[1]/'Q222_C1_fixed_admission_2026_10_26'
sys.path.insert(0,str(BASE))
from c1_native import k,f,g,a,t,dB,m,A,be,D
fix={k:S.Integer(4),g:S.Rational(1,8),a:S.Integer(1),t:S.Integer(1)}
limits={'d/f2':(dB,f**2,S.Rational(1155,17)),
'm/f2':(m,f**2,S.Integer(126)),
'A/f4':(A,f**4,S.Rational(67,34)),
'beta/f2':(be,f**2,S.Rational(68477,34)),
'D':(D,S.Integer(1),S.Rational(8,63))}
for name,(x,scale,desired) in limits.items():
    val=S.limit(S.cancel(x.subs(fix)/scale),f,0)
    assert S.cancel(val-desired)==0,(name,val,desired)
    print(name,desired)
# Source contact: u=4, v=f/2, w=1/(8f), x=-63/16, y=(1-f^2)/(2f), z=1.
# S_b/f^2 -> 536/17, S_c,S_e = O(f^4), gamma_upper^2/f^4 -> 270144/18479.
k0=S.Integer(4);v0=f/2;w0=1/(8*f);x0=-S.Rational(63,16);y0=(1-f*f)/(2*f)
sc=S.cancel((4*be+2*dB)*v0**2*x0**2).subs(fix)
sb=S.cancel(4*A*k0**2*y0**2).subs(fix)
se=S.cancel(4*be*k0**2*v0**2).subs(fix)
upper=S.cancel(64*A*m*(2*be+dB)*(4+2*D+dB)*k0**2*v0**2*w0**4/((2*be+dB+2*m)*(2*A+4+2*D+dB))).subs(fix)
assert S.limit(S.cancel(sb/f**2),f,0)==S.Rational(536,17)
assert S.limit(S.cancel(sc/f**4),f,0)>0
assert S.limit(S.cancel(se/f**4),f,0)>0
coeff=S.limit(S.cancel(upper/f**4),f,0)
assert coeff==S.Rational(270144,18479),coeff
assert (S.Rational(536,17))**2>coeff
print('PASS: gamma lower/f^2 -> 536/17 and upper/f^2 -> sqrt(270144/18479); disjoint by exact rational square comparison.')
print('PASS: an entire sufficiently small positive-f strict native C1 branch lies outside 16 features.')
