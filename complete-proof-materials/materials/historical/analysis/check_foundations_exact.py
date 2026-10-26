"""Independent finite checks of the R21 parent and source geometry.

Exact polynomial arithmetic only. These checks support named identities;
the tensor-SOS coverage requires the separate analytical arguments.
"""
from pathlib import Path
from itertools import combinations, product
import json
import time
import sympy as s

start = time.monotonic()
u, v = s.symbols('u v', real=True)
d, la, n, m, C, A, D, al, be = s.symbols('d la n m C A D al be')
weights = [la, n, m, C, A, D, al, be]
du, dv = 1+u*u, 1+v*v
rows = [s.Matrix(x) for x in [
    [1,v,0,0], [0,1,0,0], [v,-1,u*v,-u], [u,0,-1,0],
    [0,0,1,0], [u,u*v,-1,-v], [0,0,v,-1], [0,1,0,u], [0,0,0,1]]]
scalars = [4*la*du, 4*n*du*dv, 4*m, 4*C*dv, 4*A*du*dv,
           4*D+2*d, (4*al+2*d)*du, (4*be+2*d)*dv, -2*d*du*dv]
P = sum((a*r*r.T for a,r in zip(scalars,rows)),s.zeros(4))
checks = {}
def check(name, expr):
    checks[name] = bool(s.cancel(expr) == 0)
    if not checks[name]:
        raise AssertionError((name,s.factor(expr)))

Z=s.diag(1,-1); X=s.Matrix([[0,1],[1,0]])
a1=((1-u*u)*Z+2*u*X)/du
b1=((1-v*v)*Z+2*v*X)/dv
H=s.eye(4)-(s.kronecker_product(Z,Z)-s.kronecker_product(Z,b1)
             +s.kronecker_product(a1,Z)+s.kronecker_product(a1,b1))/2
core=P.subs(dict.fromkeys(weights,0))/d
for i,j in product(range(4),repeat=2):
    check(f'Pauli_parent_{i}{j}',core[i,j]-du*dv*H[i,j])

t=[d+2*(la+n+D),d+2*(al+m+A),d+2*(la+C+D),d+2*(al+be+m),
   d+2*(la+be+n),d+2*(al+C+A),d+2*(be+n+m),d+2*(C+A+D)]
t1,t2,t3,t4,t5,t6,t7,t8=t
q={(0,0):4*la*al*t7*t8,(0,1):4*(la*D*t7*t6+al*m*t5*t8),
   (0,2):4*m*D*t5*t6,(1,0):4*(al*A*t7*t3+la*n*t8*t4),
   (1,2):4*(be*m*t1*t6+C*D*t5*t2),(2,0):4*n*A*t3*t4,
   (2,1):4*(be*n*t3*t2+C*A*t1*t4),(2,2):4*be*C*t1*t2}
q[1,1]=s.expand(P.subs({u:1,v:1}).det(method='domain-ge')/256-sum(q.values()))
# Cauchy--Binet is valid also for the one negative scalar packet.
det_terms=[]
for inds in combinations(range(9),4):
    minor=s.factor(s.Matrix.hstack(*(rows[i] for i in inds)).det())
    if minor!=0:
        det_terms.append(minor*minor*s.prod(scalars[i] for i in inds))
det_exact=s.Poly(s.expand(sum(det_terms)),u,v)
target=s.Poly(s.expand(16*du**2*dv**2*sum(a*u**(2*i)*v**(2*j) for (i,j),a in q.items())),u,v)
for mon,coef in (det_exact-target).terms():
    check('det_coefficient_'+str(mon),coef)
checks['full_determinant_dictionary']=det_exact==target

# Complete finite support implication, using actual event definitions.
events=[(1,1,0,1),(1,-1,0,0),(1,-1,1,1),(-1,1,1,0),
        (-1,1,0,0),(-1,1,1,1),(-1,-1,0,1),(1,-1,1,0)]
directions=[(1,0),(-1,0),(0,1),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1)]
for mask in product([0,1],repeat=8):
    sub={d:1,**dict(zip(weights,mask))}
    active=[(i-1,j-1) for (i,j),a in q.items() if (i,j)!=(1,1) and a.subs(sub)>0]
    coercive=all(any(x*dx+y*dy<0 for x,y in active) for dx,dy in directions)
    omitted=[events[i] for i,b in enumerate(mask) if b==0]
    compatible=all(not ((x==xx and a!=aa) or (y==yy and b!=bb))
                   for (a,b,x,y),(aa,bb,xx,yy) in combinations(omitted,2))
    if coercive!=compatible:
        raise AssertionError(('support_equivalence',mask))
checks['all_256_supports_compatibility']=True

# Exact witnesses used to exclude empty and one-dimensional Newton support.
for name,psi in [('Hardy',s.Matrix([u*v,-u,u*u*v,1+v*v+u*u*v*v]))]:
    for i in [0,2,3]:
        check(name+'_event_'+str(i),(rows[i].T*psi)[0])
    check(name+'_core',(psi.T*core*psi)[0]+2*u*u*v*v*du*dv)
eps=s.symbols('eps')
psi=s.Matrix([1,-eps,eps,1])
check('alternating_support_negative_witness',
      (psi.T*P.subs({u:1,v:1,la:0,C:0,al:0,be:0})*psi)[0]
      -16*eps*((d+n+m+A+D)*eps-d))

# C1 cofactor identity EN6, on the exact omitted face.
P1=P.subs({la:0,n:0})
adj=P1.adjugate(method='domain-ge')
Fcal=(d+2*D)*((d+2*al)*(m*dv+C*u*u*v*v)+2*A*dv*(m+C*u*u))
Fcal+=2*m*C*(d*(u*u*v*v+dv)+2*al+2*D*v*v*du)
check('C1_omitted_row_positive_cofactor',u*adj[3,3]+adj[1,3]-16*u*v*v*du**2*dv*Fcal)

report={'scope':__doc__.strip(),'checks':checks,'all_exact_checks':all(checks.values()),
        'seconds':round(time.monotonic()-start,3),'sympy':s.__version__}
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
