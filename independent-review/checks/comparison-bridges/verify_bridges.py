"""Exact arithmetic for the comparison bridges; no nonlinear search or H1/H2 claim."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import sys
import sympy as S

checks = []
signs = list(product([-1, 1], repeat=4))
chsh = [e for e in signs if e[0]*e[1]*e[2]*e[3] == -1]


def facets(p):
    p = list(map(F,p))
    aa, bb, ee = p[:2], p[2:4], p[4:]
    return [(1+a*aa[x]+b*bb[y]+a*b*ee[2*x+y])/4
            for x,y,a,b in product(range(2),range(2),[-1,1],[-1,1])] + [
            2-sum(e[k]*ee[k] for k in range(4)) for e in chsh]


def interval(p,q):
    gp,gq = facets(p),facets(q)
    if any(a < 0 for a,b in zip(gp,gq) if b == 0):
        return None
    lo = max([F(0)]+[a/b for a,b in zip(gp,gq) if b < 0])
    hi = min([F(1)]+[a/b for a,b in zip(gp,gq) if b > 0])
    return (lo,hi) if lo <= hi else None


vertices = [(a0,a1,b0,b1,a0*b0,a0*b1,a1*b0,a1*b1) for a0,a1,b0,b1 in signs]
assert len(set(vertices)) == 16 and len(facets(vertices[0])) == 24
assert all(g >= 0 for v in vertices for g in facets(v))
checks.append('24 stated facets valid on all 16 deterministic vertices')

# An exact physical nonlocal q: Phi+, Alice Z,X, Bob (3Z+4X)/5,(3Z-4X)/5.
q = (0,0,0,0,F(3,5),F(3,5),F(4,5),-F(4,5))
Id,pz,px = S.eye(2),S.diag(1,-1),S.Matrix([[0,1],[1,0]])
psi = S.Matrix([1,0,0,1])/S.sqrt(2)
Alice = [pz,px]; Bob = [(3*pz+4*px)/5,(3*pz-4*px)/5]
ops = [S.kronecker_product(a,Id) for a in Alice]+[S.kronecker_product(Id,b) for b in Bob]
ops += [S.kronecker_product(a,b) for a in Alice for b in Bob]
assert all(S.simplify((psi.T*op*psi)[0]-S.Rational(value.numerator,value.denominator)) == 0
           for op,value in zip(ops,map(F,q)))
checks.append('rational nonlocal test behavior has an exact two-qubit realization')

zero = (F(0),)*8
p = tuple(F(3,4)*x for x in q)
assert interval(q,q) == (F(1),F(1))
checks.append('nonlocal pure endpoint has feasible interval {1}')
assert interval(zero,q)[0] == 0
checks.append('local input includes tau=0')
lo,hi = interval(p,q)
assert lo == F(1,8) and hi == F(41,48)
for tau in (lo,F(3,4),hi):
    local = tuple((x-tau*y)/(1-tau) for x,y in zip(p,q))
    assert all(g >= 0 for g in facets(local))
checks.append('mixed rational input: interval [1/8,41/48] and three exact local residuals')
assert interval((0,0,0,0,1,1,1,-1),q) is None
checks.append('PR input rejected for the fixed rational q')
assert interval((2,0,0,0,0,0,0,0),vertices[-1]) is None
checks.append('zero-denominator facet condition is enforced')

# Exact rational coordinate transformation and generic Bell bidegree.
u,v,z1,z2 = S.symbols('u v z1 z2',real=True)
zu=(1+S.I*u)/(1-S.I*u); zv=(1+S.I*v)/(1-S.I*v)
circle=lambda z: ((z+1/z)/2)*pz+((z-1/z)/(2*S.I))*px
assert all(S.cancel(x)==0 for x in circle(zu)-((1-u*u)*pz+2*u*px)/(1+u*u))
assert all(S.cancel(x)==0 for x in circle(zv)-((1-v*v)*pz+2*v*px)/(1+v*v))
checks.append('Cayley substitution gives both exact SG2 measurement gauges')
mu=S.Symbol('mu'); a=S.symbols('a0 a1'); b=S.symbols('b0 b1'); c=S.symbols('c00 c01 c10 c11')
obsA=[pz,circle(z1)];obsB=[pz,circle(z2)]
Q=mu*S.eye(4)
for x in range(2):Q-=a[x]*S.kronecker_product(obsA[x],Id)
for y in range(2):Q-=b[y]*S.kronecker_product(Id,obsB[y])
for x,y in product(range(2),repeat=2):Q-=c[2*x+y]*S.kronecker_product(obsA[x],obsB[y])
for entry in Q:
    poly=S.Poly(S.expand(z1*z2*entry),z1,z2)
    assert all(0 <= i <= 2 and 0 <= j <= 2 for i,j in poly.monoms())
checks.append('generic 4x4 Bell slack has Laurent bidegree at most (1,1)')

result={'scope':'targeted exact algebra and endpoint examples; analytic proofs and cited hypotheses remain necessary',
        'checks':checks,'count':len(checks),'ok':True,
        'mixed_example_interval':[str(lo),str(hi)],'H1':'open','H2':'open'}
if '--write-receipt' in sys.argv:
    Path(__file__).with_name('result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2))
