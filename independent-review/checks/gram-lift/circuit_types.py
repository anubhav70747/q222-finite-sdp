"""Does the determinant's AM-GM circuit structure predict which sources the k = 2 lift certifies?

For each of the 302 test sources (cases_302.json), take the contact of the determinant at d = d*:
the balanced weights w_e = q_e U_c^e1 V_c^e2 on the 8 outer exponents of {-1,0,1}^2.  Decompose w
into minimal circuits (4 opposite pairs, 4 triangles with weights (1,1,2), 4 with weights (1,1,1))
and ask whether the (1,1,1) triangles, which need sixth roots of U = u^2 (Hurwitz), can be avoided.

Result recorded in EXACT_LIFT_CANDIDATE.md: all 9 k = 2 gaps (gaps_k2.json) need (1,1,1) circuits,
but so do most sources that k = 2 certifies, so the scalar circuit structure is not the mechanism.
Usage: python circuit_types.py   (about 5 minutes)
"""
import json, os, statistics
from fractions import Fraction
import numpy as np, sympy as sp
from scipy.optimize import linprog
import gram_lift as gl

HERE = os.path.dirname(os.path.abspath(__file__))
E = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0)]
IDX = {e: i for i, e in enumerate(E)}


def vec(items):
    v = np.zeros(8)
    for e, w in items:
        v[IDX[e]] += w
    return v


PAIRS = [vec([((1, 0), 1), ((-1, 0), 1)]), vec([((0, 1), 1), ((0, -1), 1)]),
         vec([((1, 1), 1), ((-1, -1), 1)]), vec([((1, -1), 1), ((-1, 1), 1)])]
T112 = [vec([((1, 1), 1), ((-1, 1), 1), ((0, -1), 2)]), vec([((1, -1), 1), ((-1, -1), 1), ((0, 1), 2)]),
        vec([((1, 1), 1), ((1, -1), 1), ((-1, 0), 2)]), vec([((-1, 1), 1), ((-1, -1), 1), ((1, 0), 2)])]
T111 = [vec([((1, 1), 1), ((-1, 0), 1), ((0, -1), 1)]), vec([((1, -1), 1), ((-1, 0), 1), ((0, 1), 1)]),
        vec([((-1, 1), 1), ((1, 0), 1), ((0, -1), 1)]), vec([((-1, -1), 1), ((1, 0), 1), ((0, 1), 1)])]


def fit(w, gens, weights=None):
    """Minimal L1 residual of w in the cone of gens; optionally then minimize a weighted circuit mass."""
    A = np.array(gens).T
    m, n = A.shape
    Aeq = np.hstack([A, -np.eye(m), np.eye(m)])
    c1 = np.concatenate([np.zeros(n), np.ones(2 * m)])
    res = linprog(c1, A_eq=Aeq, b_eq=w, bounds=(0, None), method='highs').fun
    if weights is None:
        return res, None
    c2 = np.concatenate([weights, np.zeros(2 * m)])
    r2 = linprog(c2, A_ub=c1.reshape(1, -1), b_ub=[res + 1e-9], A_eq=Aeq, b_eq=w, bounds=(0, None), method='highs')
    return res, r2.fun


def contact_weights(e):
    s = max(e)
    es = [x / s for x in e]
    dsn, z = gl.dstar(gl.oracle(es))
    P = gl.pmat([gl.d] + [sp.Rational(x.numerator, x.denominator) for x in es])
    Q, rem = sp.div(sp.Poly(P.det(method='berkowitz'), gl.u, gl.v),
                    sp.Poly(sp.expand(16 * (1 + gl.u**2)**2 * (1 + gl.v**2)**2), gl.u, gl.v))
    assert rem.is_zero
    co = {(a // 2 - 1, b // 2 - 1): float(sp.Poly(c, gl.d).eval(sp.Float(dsn, 40))) for (a, b), c in Q.terms()}
    q = np.array([max(co.get(x, 0.0), 0.0) for x in E])
    with np.errstate(divide='ignore'):
        lw = np.log(q) + np.array(E, float) @ z
    w = np.exp(lw - lw[np.isfinite(lw)].max())
    w[~np.isfinite(lw)] = 0.0
    return w / w.sum(), z


if __name__ == '__main__':
    cases = json.load(open(os.path.join(HERE, 'cases_302.json')))
    gaps = {(g['fam'], g['idx']) for g in json.load(open(os.path.join(HERE, 'gaps_k2.json')))}
    out = []
    for c in cases:
        e = [Fraction(x).limit_denominator(10**30) for x in c['e']]
        w, z = contact_weights(e)
        res_no111, _ = fit(w, PAIRS + T112)
        res_all, share = fit(w, PAIRS + T112 + T111, np.array([0] * 8 + [3] * 4, float))
        out.append(dict(gap=(c['fam'], c['idx']) in gaps, needs111=res_no111 > res_all + 1e-7, share111=share,
                        name='%s %d' % (c['fam'], c['idx']), zu=float(z[0]), zv=float(z[1]),
                        balance=float(np.abs(np.array(E, float).T @ w).max())))
    g = [o for o in out if o['gap']]
    k = [o for o in out if not o['gap']]
    print('sources:', len(out), ' k=2 gaps:', len(g), ' max balance error: %.1e' % max(o['balance'] for o in out))
    print('k=2 gaps that need (1,1,1) circuits:      %d / %d' % (sum(o['needs111'] for o in g), len(g)))
    print('k=2-certified that need (1,1,1) circuits: %d / %d' % (sum(o['needs111'] for o in k), len(k)))
    for key in ('zu', 'zv'):
        sm = lambda rows: tuple(round(f(o[key] for o in rows), 2) for f in (min, lambda x: statistics.median(list(x)), max))
        print('contact log-%s (min, median, max): gaps %s, certified %s' % (key[1], sm(g), sm(k)))
