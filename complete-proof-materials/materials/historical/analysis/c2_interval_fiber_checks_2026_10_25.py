"""Exact low-degree native C2 identities; no optimization or sampling."""
import json
from pathlib import Path
import sympy as s

z,k,b,r,t=s.symbols('z k b r t')
h=b-k*t
Ac=1+r*z
g=b+r-1+h*z
j=1-t+r*z
E=1-h*z
Delta=((1+z)*z*b*b+2*z*(1+z)*b*r+(1+z)*(1+k*z)*b
       +z*z*(1-k)*r*r+2*z*(1+k*z)*r+1+k*z)
N=(1+z)*(b+r)*(b+k*(1-t))-k*(1+k*z)*(1-r)*t
q=(b+r)*Ac*(Ac+k*z*t)-r*Delta
q0=b*(1-r)
q1=b*(b+2*r)*(1-r)-(b+r)*h-k*r*(b+1)
K=r*((b+r)*h+k*b+k*r*(2-r))
assert s.expand(q-(q0+q1*z-K*z*z))==0
assert s.expand(Delta-((b+1)+(b+1)*(b+k+2*r)*z
                   +((b+r)**2+k*(b+r*(2-r)))*z*z))==0
assert s.expand(N-((1+z)*(b+r)*h+k*(g+j+r*t*(1+k*z))))==0

source=[1,0,q/(2*r*Delta),
        (1+k*z)*(1-r)*Ac*(Ac+k*z*t)/(2*g*Delta),
        Ac*N/(2*k*t*Delta),
        k*Ac*(b+r)*(Ac+k*z*t)/(2*h*Delta),
        (1+z)*(b+r)*(b-k*z*r)*(Ac+k*z*t)/(2*k*j*Delta),
        z*z*(b+r)*N/(2*E*Delta),0]
scale=2*Delta*g*j*E
poly=[scale,0,q*g*j*E/r,
      (1+k*z)*(1-r)*Ac*(Ac+k*z*t)*j*E,
      Ac*N*g*j*E/(k*t),
      k*Ac*(b+r)*(Ac+k*z*t)*g*j*E/h,
      (1+z)*(b+r)*(b-k*z*r)*(Ac+k*z*t)*g*E/k,
      z*z*(b+r)*N*g*j,0]
degrees=[]
for c,p in zip(source,poly):
    assert s.cancel(scale*c-p)==0
    degrees.append(-1 if p==0 else s.Poly(s.cancel(p),z).degree())
assert max(degrees)<=5

# A simultaneous lower-end gap zero really needs common-factor removal.
theta={k:s.Rational(1,4),b:s.Rational(19,32),r:s.Rational(1,4),t:s.Rational(9,8)}
sp=[s.Poly(s.cancel(p.subs(theta)),z) for p in poly if p!=0]
common=sp[0]
for p in sp[1:]: common=s.gcd(common,p)
assert common.monic().as_expr()==z-s.Rational(1,2)
assert q.subs(theta).subs(z,s.Rational(1,2))==s.Rational(243,512)

data={
 'source_order':['d','lambda','n','m','C','A','D','alpha','beta'],
 'n_numerator_identity':True,'Delta_positive_coefficient_identity':True,
 'N_positive_identity':True,'all_scaled_source_identities':True,
 'source_degrees_in_z':degrees,'max_degree':max(degrees),
 'source_multiplier':'2*Delta*g*j*(1-z*h)',
 'fixed_shape_moment_PSD_orders':[3,3],
 'common_factor_example':{str(v):str(w) for v,w in theta.items()},
 'example_common_monic_factor':str(common.monic().as_expr()),
 'example_n_numerator_at_lower_endpoint':'243/512',
 'scope':'Exact identities for fixed-shape interval reduction; no full-C2 finite lift asserted.',
 'no_solver_or_numerical_sweep':True}
Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
print(json.dumps(data,indent=2))
