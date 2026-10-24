"""Part B: reconstruct R6.5, R6.6, R6.7 (remote charts, n0=0) from the C2 parent (RB5)
and check the strictness / determinant claims of c2_short.tex:12397-12405, 12872-12878.
"""
import sympy as sp

u, v, W, X, Y, Z = sp.symbols('u v W X Y Z', real=True)
p, q, zeta, eta = sp.symbols('p q zeta eta', real=True)
n, m, C, A, D, al, K = sp.symbols('n m C A D alpha K_m', real=True)
h, t = sp.symbols('h t', positive=True)
half = sp.Rational(1, 2)

def parent(u=u, v=v, W=W, X=X, Y=Y, Z=Z, n=n, m=m, C=C, A=A, D=D, al=al):
    Du, Dv = 1 + u**2, 1 + v**2
    M = v*W - X + u*v*Y - u*Z
    J = u*W + u*v*X - Y - v*Z
    return (n*Du*Dv*X**2 + m*M**2 + C*Dv*(u*W - Y)**2 + A*Du*Dv*Y**2
            + (D + half)*J**2 + (al + half)*Du*(v*Y - Z)**2
            + half*Dv*(X + u*Z)**2 - half*Du*Dv*Z**2)

P = sp.expand(parent())
state = [W, X, Y, Z]

def gram(expr, roots):
    """symmetric matrix of a quadratic form in roots (roots are symbols)."""
    expr = sp.expand(expr)
    k = len(roots)
    G = sp.zeros(k, k)
    for i in range(k):
        for jj in range(k):
            if i == jj:
                G[i, i] = expr.coeff(roots[i], 2)
            else:
                G[i, jj] = half*expr.coeff(roots[i], 1).coeff(roots[jj], 1)
    # verify
    zv = sp.Matrix(roots)
    assert sp.expand((zv.T*G*zv)[0] - expr) == 0
    return G

Cs = (h + t)/(2*(1 - h))
As = t/(2*h)
Ds = 1/(2*(t - (1 - h)))
resid = {n: 0, al: 0, C: Cs, A: As, D: Ds}

print('==== R6.5: u unbounded (p=1/u), v bounded unit ====')
P5 = sp.expand(sp.cancel(p**2*P.subs(u, 1/p)))
Mrow5 = sp.expand(p*(v*W - X + (1/p)*v*Y - (1/p)*Z))   # = vY - Z + p(vW - X)
G5expr = sp.expand(P5 - (m - K)*Mrow5**2)
assert G5expr.coeff(m) == 0 if G5expr.has(m) else True
print('m eliminated from remainder:', not G5expr.has(m))
G5 = gram(G5expr, state)
R5 = sp.simplify(G5.subs(p, 0).subs(resid))
print('residue matrix at p=0 (n0=alpha0=0):'); sp.pprint(R5)
# explicit claimed decomposition: C*Dv W^2 + A*Dv Y^2 + (D+1/2)(W+vX)^2 + (K+1/2)(vY-Z)^2
claim5 = (Cs*(1 + v**2)*W**2 + As*(1 + v**2)*Y**2 + (Ds + half)*(W + v*X)**2
          + (K + half)*(v*Y - Z)**2)
zv = sp.Matrix(state)
print('residue == C*Dv W^2 + A*Dv Y^2 + (D*+1/2)(W+vX)^2 + (K+1/2)(vY-Z)^2 :',
      sp.simplify(sp.expand((zv.T*R5*zv)[0] - claim5)) == 0)
det5 = sp.factor(R5.det())
print('det R5 =', det5)
minors5 = [sp.factor(R5[:k, :k].det()) for k in range(1, 5)]
print('leading minors:', minors5)

print('\n==== R6.6: v unbounded (q=1/v), u bounded ====')
P6 = sp.expand(sp.cancel(q**2*P.subs(v, 1/q)))
Mrow6 = sp.expand(q*((1/q)*W - X + u*(1/q)*Y - u*Z))  # = W + uY - q(X + uZ)
G6expr = sp.expand(P6 - (m - K)*Mrow6**2)
print('m eliminated from remainder:', not G6expr.has(m))
G6 = gram(G6expr, state)
R6m = sp.simplify(G6.subs(q, 0).subs(resid))
print('residue matrix at q=0:'); sp.pprint(R6m)
# m-row kernel at q=0: W = -u Y ; basis (Y, X, Z)
Bk = sp.Matrix([[-u, 0, 0], [0, 1, 0], [1, 0, 0], [0, 0, 1]])  # columns: Y-dir, X, Z in (W,X,Y,Z)
R6k = sp.simplify(Bk.T*R6m.subs(K, 0)*Bk)
print('residue (K=0) restricted to m-row kernel W=-uY, basis (Y,X,Z):'); sp.pprint(sp.simplify(R6k))
print('kernel det =', sp.factor(R6k.det()))
print('kernel X,Z block det =', sp.factor(R6k[1:, 1:].det()))
det6 = sp.factor(sp.expand(R6m.det()))
print('full residue det (poly in K) =', det6)
print('coefficient of K in det:', sp.factor(sp.expand(R6m.det()).coeff(K, 1)),
      ' K^0 term:', sp.factor(sp.expand(R6m.det()).coeff(K, 0)))

print('\n==== R6.7: both unbounded, roots (W,X,Y,zeta=pZ,eta=qZ) ====')
P7 = sp.expand(sp.cancel(p**2*q**2*P.subs({u: 1/p, v: 1/q})))
Mrow7 = p*W - p*q*X + Y - q*Z
# BA15-type allocation with m retained as (m-K) square
alloc = (n*(1 + p**2)*(1 + q**2)*X**2 + m*(p*W - p*q*X + Y - eta)**2
         + C*(1 + q**2)*(W - p*Y)**2 + A*(1 + p**2)*(1 + q**2)*Y**2
         + D*(q*W + X - p*q*Y - zeta)**2 + al*(1 + p**2)*(Y - eta)**2
         + half*((q*W + X - p*q*Y - zeta)**2 + (1 + p**2)*(Y - eta)**2
                 + (1 + q**2)*(p**2*X**2 + 2*X*zeta - zeta**2)))
print('allocation identity (zeta=pZ, eta=qZ) holds:',
      sp.expand(alloc.subs({zeta: p*Z, eta: q*Z}) - P7) == 0)
G7expr = sp.expand(alloc - (m - K)*(p*W - p*q*X + Y - eta)**2)
print('m eliminated:', not G7expr.has(m))
roots7 = [W, X, Y, zeta, eta]
G7 = gram(G7expr, roots7)
R7 = sp.simplify(G7.subs({p: 0, q: 0}).subs(resid))
print('residue matrix on (W,X,Y,zeta,eta) at p=q=0:'); sp.pprint(R7)
print('(X,zeta) block det =', sp.factor(R7.extract([1, 3], [1, 3]).det()), ' expected D*/2 =', sp.factor(Ds/2))
print('(Y,eta) block det =', sp.factor(R7.extract([2, 4], [2, 4]).det()), ' expected A*(K+1/2) =', sp.factor(As*(K + half)))
print('cross entries between blocks zero:',
      all(R7[i, jj] == 0 for i in (0,) for jj in (1, 2, 3, 4)) and
      all(R7[i, jj] == 0 for i in (1, 3) for jj in (2, 4)))
print('det R7 =', sp.factor(R7.det()))
# literal-endpoint check: state residue of p^2 q^2 P at p=q=0 in (W,X,Y,Z) is singular in Z
R7state = gram(sp.expand(P7 - (m - K)*Mrow7**2), state).subs({p: 0, q: 0}).subs(resid)
print('for comparison, residue in raw state roots (W,X,Y,Z): Z-row =', list(R7state[3, :]))
