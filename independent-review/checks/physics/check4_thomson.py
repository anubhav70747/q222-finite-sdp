"""Check 4: Thomson/flow identity.  Weights a_k and d are independent symbols (pure algebra).
Parametrise t_k = sqrt(r) rho_k / sqrt(c_k) with sum_k rho_k h_k = e4  (then sum t_k g_k = ell, r = 2 d Du Dv).
Then a_k (g_k.psi - t_k s/a_k)^2 = a_k c_k (h_k.psi - r rho_k Z/(c_k a_k))^2 and sum t_k^2/a_k = r sum rho_k^2/(c_k a_k)."""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import numpy as np

a = list(symbols('a1:9', positive=True))
Tg = T_of(a)
Fg = (psi.T*Tg*psi)[0] - r*Z**2                      # F = psi^T T psi - s^2

# general solution of sum rho_k h_k = e4: 4 free parameters f, pivots = (lam, n, A, m)
f = symbols('f1:5')
piv, free = [0, 1, 4, 2], [3, 5, 6, 7]
rho = [None]*8
for k, fk in zip(free, f):
    rho[k] = fk
sol = sp.solve(list(sum((rho[k]*h[k] for k in free), zeros(4, 1)) + sum((symbols(f'p{k}')*h[k] for k in piv), zeros(4, 1)) - e4),
               [symbols(f'p{k}') for k in piv], dict=True)[0]
for k in piv:
    rho[k] = sol[symbols(f'p{k}')]
assert sp.simplify(sum((rho[k]*h[k] for k in range(8)), zeros(4, 1)) - e4) == zeros(4, 1)

flow = sum(a[k]*c[k]*(sp.expand((h[k].T*psi)[0] - r*rho[k]*Z/(c[k]*a[k])))**2 for k in range(8))
cost = r*sum(rho[k]**2/(c[k]*a[k]) for k in range(8))     # sum t_k^2/a_k
print('Thomson identity (symbolic, free flow params f1..f4):',
      sp.cancel(sp.together(Fg - (flow + r*Z**2*(1 - cost)))) == 0)

# --- minimiser = ell^T T^{-1} ell, and F>=0 <=> ell^T T^{-1} ell <= 1, at random rational points (exact)
rng = random.Random(7)
agree, ok_min, n_pt, n_neg = 0, 0, 40, 0
for _ in range(n_pt):
    sub = {u: rand_q(rng, 1, 5), v: rand_q(rng, 1, 5), d: rand_q(rng, 1, 9)}
    sub.update({a[k]: rand_q(rng, 1, 30)*Rational(1, rng.choice([1, 10, 100])) for k in range(8)})
    Tn = Tg.subs(sub); rn = r.subs(sub)
    H = Matrix.hstack(*[h[k].subs(sub) for k in range(8)])          # 4x8
    b = [(c[k]*a[k]).subs(sub) for k in range(8)]
    Tinv = Tn.inv()
    rho_s = [b[k]*(H[:, k].T*Tinv*e4)[0] for k in range(8)]          # stationary point
    assert H*Matrix(rho_s) == e4                                      # feasible
    cost_s = rn*sum(rho_s[k]**2/b[k] for k in range(8))
    ellTinvell = rn*Tinv[3, 3]
    # minimality: for any null vector nu of H, cost(rho*+nu) - cost(rho*) = r sum nu^2/b > 0 (cross term vanishes)
    nu = H.nullspace()[rng.randrange(4)]
    cross = sum(rho_s[k]*nu[k]/b[k] for k in range(8))
    ok_min += bool(cost_s == ellTinvell and cross == 0)
    # sign test:  det(P) = det(T)(1 - ell^T T^-1 ell);  compare with numeric eigenvalues of P = T - ell ell^T
    Pn = np.array((Tn - rn*e4*e4.T).evalf(30).tolist(), dtype=float)
    lam_min = np.linalg.eigvalsh(Pn).min()
    agree += bool(lam_min >= -1e-9) == bool(ellTinvell <= 1)
    n_neg += bool(ellTinvell > 1)
print(f'min cost == ell^T T^-1 ell and first-order cond (cross term=0) : {ok_min}/{n_pt}')
print(f'(P psd) <=> (ell^T T^-1 ell <= 1): {agree}/{n_pt}  (violating points: {n_neg})')
