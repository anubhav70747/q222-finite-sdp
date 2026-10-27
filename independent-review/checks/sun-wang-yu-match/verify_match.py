"""Exact SG1/SG2-to-Sun-Wang-Yu matching; no numerical search or H1/H2 claim."""
from pathlib import Path
import json
import sys
import sympy as S

checks = []


def zero(name, value):
    entries = list(value) if isinstance(value, S.MatrixBase) else [value]
    assert all(S.cancel(S.expand(x)) == 0 for x in entries), name
    checks.append(name)


d, la, n, m, C, A, D, al, be = S.symbols('d lambda n m C A D alpha beta')
weights = [la, n, m, C, A, D, al, be]
u, v, W, X, Y, Z = S.symbols('u v W X Y Z', real=True)
av = S.symbols('A0 A1'); bv = S.symbols('B0 B1')
ev = S.Matrix(2, 2, S.symbols('E00 E01 E10 E11'))
events = [(0,1,1,1), (0,0,1,-1), (1,1,1,-1), (1,0,-1,1),
          (0,0,-1,1), (1,1,-1,1), (0,1,-1,-1), (1,0,1,-1)]
f = d*(1-(ev[0,0]-ev[0,1]+ev[1,0]+ev[1,1])/2)
for weight, (x,y,a,b) in zip(weights, events):
    f += weight*(1+a*av[x]+b*bv[y]+a*b*ev[x,y])
kappa = d+sum(weights)
avec = [A+al-la-n, C+D-m-be]
bvec = [n+be-C-A, m+al-la-D]
corr = S.Matrix([[d/2+n+A, -d/2-la-al], [d/2+C+be, d/2+m+D]])
bell = sum(avec[x]*av[x] for x in range(2))+sum(bvec[y]*bv[y] for y in range(2))
bell += sum(corr[x,y]*ev[x,y] for x in range(2) for y in range(2))
zero('SG1 affine constant, biases and four correlators', f-(kappa-bell))
zero('all-positive deterministic assignment has slack 4 lambda',
     f.subs({q:1 for q in [*av,*bv,*ev]})-4*la)
zero('frustration product factorization', S.prod(corr)+
     (d/2+n+A)*(d/2+la+al)*(d/2+C+be)*(d/2+m+D))

Du, Dv = 1+u*u, 1+v*v
Id, pz, px = S.eye(2), S.diag(1,-1), S.Matrix([[0,1],[1,0]])
obsA = [pz, ((1-u*u)*pz+2*u*px)/Du]
obsB = [pz, ((1-v*v)*pz+2*v*px)/Dv]
op = S.zeros(4)
for x in range(2):
    op += avec[x]*S.kronecker_product(obsA[x],Id)
for y in range(2):
    op += bvec[y]*S.kronecker_product(Id,obsB[y])
for x in range(2):
    for y in range(2):
        op += corr[x,y]*S.kronecker_product(obsA[x],obsB[y])
parent = (4*la*Du*(W+v*X)**2+4*n*Du*Dv*X**2
          +4*m*(v*W-X+u*v*Y-u*Z)**2+4*C*Dv*(u*W-Y)**2
          +4*A*Du*Dv*Y**2+(4*D+2*d)*(u*W+u*v*X-Y-v*Z)**2
          +(4*al+2*d)*Du*(v*Y-Z)**2+(4*be+2*d)*Dv*(X+u*Z)**2
          -2*d*Du*Dv*Z**2)
P = S.hessian(parent, (W,X,Y,Z))/2
zero('full SG2 matrix equals Du*Dv*(kappa I-B)', P-Du*Dv*(kappa*S.eye(4)-op))
QA, QB = S.kronecker_product(pz,Id), S.kronecker_product(Id,pz)
zero('negative Alice angle is local state congruence', P.subs(u,-u)-QA*P*QA)
zero('negative Bob angle is local state congruence', P.subs(v,-v)-QB*P*QB)
su, tv = (1-u*u)/Du, (1-v*v)/Dv
zero('overlap Jacobian', S.Matrix([su,tv]).jacobian([u,v])-
     S.diag(-4*u/Du**2,-4*v/Dv**2))
zero('positive-chart square-root identity Alice', 1-su**2-(2*u/Du)**2)
zero('positive-chart square-root identity Bob', 1-tv**2-(2*v/Dv)**2)

# Independently expand the square-root-free Gram trace formula (SWY Eq.36).
s, t, mu = S.symbols('s t mu')
aa, bb = S.Matrix(S.symbols('a0 a1')), S.Matrix(S.symbols('b0 b1'))
cc = S.Matrix(2,2,S.symbols('c00 c01 c10 c11'))
Gs, Gt = S.Matrix([[1,s],[s,1]]), S.Matrix([[1,t],[t,1]])
scalar = lambda mat: mat[0,0]
w = scalar(aa.T*Gs*aa)+scalar(bb.T*Gt*bb)+S.trace(cc.T*Gs*cc*Gt)
triple = scalar(aa.T*Gs*cc*Gt*bb)
rr = (scalar(aa.T*Gs*cc*Gt*cc.T*Gs*aa)
      +scalar(bb.T*Gt*cc.T*Gs*cc*Gt*bb)
      +scalar(aa.T*Gs*aa)*scalar(bb.T*Gt*bb)
      +cc.det()**2*(1-s*s)*(1-t*t))
delta = S.Poly(S.expand(mu**4-2*w*mu**2-8*triple*mu+w*w-4*rr),s,t)
assert delta.degree(s) == delta.degree(t) == 2
checks.append('determinant trace polynomial has bidegree (2,2)')
zero('determinant leading mixed coefficient', delta.coeff_monomial(s*s*t*t)-16*S.prod(cc))

# The new lemma uses the only possible total-degree-four monomial.
eps, s0, t0, w1, w2 = S.symbols('eps s0 t0 w1 w2')
hij = {(i,j): S.Symbol(f'h{i}{j}') for i in range(3) for j in range(3)}
h = sum(hij[i,j]*s**i*t**j for i in range(3) for j in range(3))
line = S.Poly(S.expand(h.subs({s:s0+eps*w1,t:t0+eps*w2})),eps)
zero('null-direction quartic coefficient', line.coeff_monomial(eps**4)-hij[2,2]*w1*w1*w2*w2)
assert S.Poly(line.as_expr().subs(w1,0),eps).degree() <= 2
assert S.Poly(line.as_expr().subs(w2,0),eps).degree() <= 2
checks.append('axis restrictions have degree at most two')

result = {'scope':'exact symbolic coordinate and coefficient matching only; analytic hypotheses in accompanying proof',
          'sympy_version':S.__version__, 'checks':checks, 'count':len(checks), 'ok':True,
          'H1':'open', 'H2':'open', 'global_exactness':'conditional'}
if '--write-receipt' in sys.argv:
    Path(__file__).with_name('result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
