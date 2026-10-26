#!/usr/bin/env python3
"""Exact native C2 full-parent stationarity at moving contact.

Runs six exact polynomial-numerator identities, with denominator clearing
performed as a monomial product of native source gaps. No sample points.
Requires local c2_fixed_map.py and SymPy. Runtime depends on SymPy version.
"""
import sympy as sp
import c2_fixed_map as core

a,b,c,t,r=core.a,core.b,core.c,core.t,core.r
n,m,C,A,D,al=core.n,core.m,core.C,core.A,core.D,core.al
z=sp.Symbol('z')
Ac=1+z*r;h=b-a*a*t;g=(1+z)*b-(1-r+a*a*z*t);j=1+z*r-t
Del=1+(1+z)*b+z*(1+z)*b*b+2*z*(1+z)*b*r+z*z*r*r+2*z*r+a*a*z*((1+z)*b+1+z*r*(2-r))
N=(1+z)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z)*(1-r)*t
P=r*g*a*a*t*h*j*(1-z*h)
num={n:(b+r)*Ac*(Ac+a*a*z*t)-r*Del,
     m:(1+a*a*z)*(1-r)*Ac*(Ac+a*a*z*t), C:Ac*N,
     A:a*a*Ac*(b+r)*(Ac+a*a*z*t),
     D:(1+z)*(b+r)*(b-a*a*z*r)*(Ac+a*a*z*t),
     al:z*z*(b+r)*N}
comult={n:g*a*a*t*h*j*(1-z*h),m:r*a*a*t*h*j*(1-z*h),
        C:r*g*h*j*(1-z*h),A:r*g*a*a*t*j*(1-z*h),
        D:r*g*t*h*(1-z*h),al:r*g*a*a*t*h*j}

contact={core.u:c,core.v:a*c,core.W:b/a,core.X:-c*r,
         core.Y:c*(b/a-a*t),core.Z:sp.Integer(1)}
F=core.full_parent()
for var,cpow in [(core.W,2),(core.X,1),(core.Y,1),(core.Z,0),(core.u,1),(core.v,1)]:
    raw=sp.expand(sp.diff(F,var).subs(contact))
    raw=sp.Poly(sp.cancel(raw/c**cpow),c)
    if any(k[0]%2 for k in raw.monoms()): raise ValueError('odd c power')
    reduced=sp.expand(sum(q*c**0*z**(power[0]//2) for power,q in raw.terms()))
    scalar=reduced.subs({n:0,m:0,C:0,A:0,D:0,al:0})
    polynomial=2*Del*P*scalar
    for symbol in (n,m,C,A,D,al):
        polynomial+=sp.diff(reduced,symbol)*num[symbol]*comult[symbol]
    polynomial=sp.Poly(sp.expand(polynomial),a,b,z,t,r)
    assert polynomial.is_zero, (var,len(polynomial.terms()))
    print('PASS native full-state/angle contact derivative:',str(var),flush=True)
print('ALL SIX EXACT NATIVE STATIONARITY POLYNOMIALS ZERO')
