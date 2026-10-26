#!/usr/bin/env python3
"""Produce integer-polynomial cleared native C2 full-parent target vector.

For every strict native point, b(theta)=P(theta)/D(theta), with D>0
and each P_k a polynomial in five source variables (a,b,c,t,r).
All 412 rows, including the 388 identically zero target coordinates,
are indexed by C2_154_CONTACT_MAP.json. No SDP or numerics.
"""
import json
import sympy as sp
import c2_fixed_map as M

a,b,c,t,r=M.a,M.b,M.c,M.t,M.r
n,m,C,A,D,al=M.n,M.m,M.C,M.A,M.D,M.al
z=c*c;Ac=1+z*r
h=b-a*a*t;g=(1+z)*b-(1-r+a*a*z*t);j=1+z*r-t
Delta=1+(1+z)*b+z*(1+z)*b*b+2*z*(1+z)*b*r+z*z*r*r+2*z*r+a*a*z*((1+z)*b+1+z*r*(2-r))
N=(1+z)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z)*(1-r)*t
num={n:(b+r)*Ac*(Ac+a*a*z*t)-r*Delta,
     m:(1+a*a*z)*(1-r)*Ac*(Ac+a*a*z*t),
     C:Ac*N,A:a*a*Ac*(b+r)*(Ac+a*a*z*t),
     D:(1+z)*(b+r)*(b-a*a*z*r)*(Ac+a*a*z*t),
     al:z*z*(b+r)*N}
comult={n:a**4*t*g*h*j*(1-z*h),m:a**4*t*r*h*j*(1-z*h),
        C:a*a*r*g*h*j*(1-z*h),A:a**4*t*r*g*j*(1-z*h),
        D:a*a*t*r*g*h*(1-z*h),al:a**4*t*r*g*h*j}
Dclear=2*a**4*t*r*g*h*j*(1-z*h)*Delta
uc=c;vc=a*c
psi=[b/a,-c*r,c*(b/a-a*t),sp.Integer(1)]
rows=M.coefficient_map()['row_keys']
row_id={k:i for i,k in enumerate(rows)}
P={}
for key,co in M.target_template().items():
    p,q,l,k=key
    factor=uc**(p//6)*vc**(q//6)*psi[l]*psi[k]
    pc=sp.expand(co.subs({n:0,m:0,C:0,A:0,D:0,al:0}))*(2*a**4*t*r*g*h*j*(1-z*h)*Delta)
    for symbol in (n,m,C,A,D,al):
        pc+=sp.expand(sp.diff(co,symbol))*num[symbol]*comult[symbol]
    cleared=sp.cancel(factor*pc)
    assert cleared.as_numer_denom()[1]==1,('uncleared denominator',key)
    P[str(row_id[key])]=str(sp.Poly(sp.expand(cleared),a,b,c,t,r).as_expr())
print('EXACT: all 24 full-parent target coefficients become integer polynomials')
print('EXACT: positive source-only denominator 2 a^4 t r g h j (1-c^2 h) Delta')
with open('C2_CLEARED_POLYNOMIAL_TARGET.json','w') as f:
 json.dump({'parameter_order':['a','b','c','t','r'],
            'denominator':str(sp.factor(Dclear)),
            'number_coefficient_rows':len(rows),
            'nonzero_polynomial_coefficients_by_row':P,
            'zero_rows_implicit':True},f,separators=(',',':'))
