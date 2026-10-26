"""exact_primal: pointwise EXACT certificate of mu(theta) > 0  (primal side).  usage: python -I exact_primal.py FAM p1..p5 [--allow-outside]
 1. rationalise theta (limit_denominator 1e6), verify native chart exactly (nd_core.chart_ok); --allow-outside skips the refusal for a source outside the chart (negative controls only)
 2. equilibrated double SDP -> G' (x coords), float affine projection (3 Newton-like rounds) onto A_S x = b_S
 3. exact: x_rat = dyadic(x), exact residual r = b_S - A_S x_rat (Fractions), exact solve (A_S A_S^T) w = r (flint fmpq_mat), x_ex = x_rat + A_S^T w  => A_S x_ex = b_S EXACTLY;
    also all 412 rows checked exactly (dropped rows).  4. exact LDL^T of every G'_b(x_ex): exact PD  =>  mu(theta) > 0 rigorously (strict: G' > 0 on the contact complement)."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import numpy as np, scipy.linalg as sla
import nd_core
from hcore import equilibrated, make_struct, beta_rows, beta_all, cx_vec
import exact_check as ec
import flint
from flint import fmpq, fmpq_mat, arb
flint.ctx.prec = 400

def ball_pd_fr(M):
    """rigorous PD test of an EXACT rational symmetric matrix (list of lists of Fractions): ball Cholesky at 400 bits; True => positive definite (proved)."""
    n = len(M); A = [[arb(fmpq(v.numerator, v.denominator)) for v in row] for row in M]; L = [[None]*n for _ in range(n)]
    for j in range(n):
        s = A[j][j]
        for k in range(j): s -= L[j][k]*L[j][k]
        if not (s > 0): return False
        d = s.sqrt(); L[j][j] = d
        for i in range(j + 1, n):
            s = A[i][j]
            for k in range(j): s -= L[i][k]*L[j][k]
            L[i][j] = s/d
    return True

def exact_b(S, eb):
    b = [Fr(0)]*len(S.keys)
    for k, x in eb.items(): b[S.kidx[k]] = Fr(x)
    return b

def rows_exact_dot(ent_row, x):
    s = Fr(0)
    for c, v in ent_row.items(): s += v*x[c]
    return s

def exact_primal(fam, p, passes=4, verbose=True, th=None, r=None, allow_outside=False):
    th = th or tuple(Fr(float(x)).limit_denominator(10**6) for x in p)
    ok, src, ct = nd_core.chart_ok(fam, th)
    if not ok and not allow_outside: return dict(error='rationalised point outside native chart')
    pf = np.array([float(x) for x in th]); t0 = time.time()
    if r is None: r = equilibrated(fam, pf, passes=passes, check_ok=bool(ok))
    S = make_struct(fam, r['e']); keep = r['keep']; bf = beta_rows(S, beta_all(fam, pf[:, None])[0][0]); cx = cx_vec(S)
    Ax = S.A[keep].tocsr(); bS = bf[keep]; x = (r['y']*cx).copy()
    K = (Ax @ Ax.T).toarray(); Kf = sla.cho_factor(K)
    # --- G-metric (affine-invariant) Newton projection: x <- x + L dH L^T, dH = min-Frobenius solution of A(L dH L^T) = residual  (keeps G' >= (1-|dH|) G'_ipm, so a tiny IPM margin survives)
    from pcsolve import tmap
    gproj_info = []
    for it in range(6):
        res = bS - Ax @ x
        if np.abs(res).max() < 1e-15*max(1.0, np.abs(bS).max()): break
        try:
            Ls = [np.linalg.cholesky(0.5*(G + G.T)) for G in S.unvec(x)]
        except np.linalg.LinAlgError: gproj_info.append('chol-fail'); break
        Ts = [tmap(L) for L in Ls]; Ay = np.hstack([Ax[:, S.offs[bi]:S.offs[bi + 1]].toarray() @ Ts[bi] for bi in range(len(Ls))])*cx[None, :]
        dy, *_ = np.linalg.lstsq(Ay, res, rcond=1e-14); dxH = dy*cx
        dH = [S.unvec(dxH)[bi] for bi in range(len(Ls))]; mdH = max(np.abs(np.linalg.eigvalsh(D)).max() for D in dH)
        x = x + np.concatenate([Ts[bi] @ dxH[S.offs[bi]:S.offs[bi + 1]] for bi in range(len(Ls))]); gproj_info.append((float(np.abs(res).max()), float(mdH)))
    for it in range(2):
        res = bS - Ax @ x; x = x + Ax.T @ sla.cho_solve(Kf, res)
    def lam(xv):
        Gs = S.unvec(xv); return [np.linalg.eigvalsh(G)[0] for G in Gs]
    def lam_sh(xv, tt):
        out = []
        for bi, G in enumerate(S.unvec(xv)):
            m = G.shape[0]; a = np.array([float(z) for z in S.info[bi]['a']]); Sh = np.eye(m) - np.outer(a, a)/(1 + a @ a); out.append(np.linalg.eigvalsh(G - tt*Sh)[0])
        return out
    maxdiag = r['maxdiag']; res_f = np.abs(bS - Ax @ x).max()
    if verbose: print('float: t=%.4e trel=%.4e %s | after G-metric projection %s: float residual %.2e, lam_min(G\') per block %s, lam_min(G\'-0.5 t Sh) %s' % (r['t'], r['trel'], r['status'], gproj_info, res_f, ['%.2e' % v for v in lam(x)], ['%.2e' % v for v in lam_sh(x, 0.5*r['t'])]), flush=True)
    # ---- exact
    eb = ec.exact_beta(fam, th); b = exact_b(S, eb)
    xr = [Fr(float(v)) for v in x]; print('exact: rational x built %.0fs' % (time.time() - t0), flush=True)
    res_ex = [b[i] - rows_exact_dot(S.ent[i], xr) for i in keep]
    rmax = max(abs(float(v)) for v in res_ex); print('exact: residual computed %.0fs' % (time.time() - t0), flush=True)
    n = len(keep)
    Ks = [[Fr(0)]*n for _ in range(n)]; ents = [S.ent[i] for i in keep]
    for i in range(n):
        di = ents[i]
        for j in range(i, n):
            dj = ents[j]; a, bb = (di, dj) if len(di) <= len(dj) else (dj, di); s = Fr(0)
            for c, v in a.items():
                w = bb.get(c)
                if w is not None: s += v*w
            Ks[i][j] = Ks[j][i] = s
    t1 = time.time(); print('exact: Gram K built %.0fs' % (time.time() - t0), flush=True)
    Kq = fmpq_mat(n, n, [fmpq(v.numerator, v.denominator) for row in Ks for v in row]); rq = fmpq_mat(n, 1, [fmpq(v.numerator, v.denominator) for v in res_ex])
    w = Kq.solve(rq); wl = [Fr(int(w[i, 0].p), int(w[i, 0].q)) for i in range(n)]
    print('exact: solve done %.0fs' % (time.time() - t0), flush=True); xe = list(xr)
    for i in range(n):
        wi = wl[i]
        if wi == 0: continue
        for c, v in ents[i].items(): xe[c] += v*wi
    # exact residual over ALL 412 rows
    bad = [i for i in range(len(S.keys)) if rows_exact_dot(S.ent[i], xe) != b[i]]
    if verbose: print('exact projection: |res| before %.2e; time %.0fs; rows with nonzero exact residual after: %d of %d' % (rmax, time.time() - t1, len(bad), len(S.keys)), flush=True)
    Ms = []
    for bi, m in enumerate(S.ms):
        M = [[Fr(0)]*m for _ in range(m)]; off = int(S.offs[bi])
        for q in range(m):
            for pp in range(q + 1):
                v = xe[off + q*(q + 1)//2 + pp]; M[pp][q] = M[q][pp] = v
        Ms.append(M)
    pd = [ball_pd_fr(M) for M in Ms]
    out = dict(th=[str(v) for v in th], t=r['t'], trel=r['trel'], status=r['status'], exact_residual_rows=len(bad), exact_pd_blocks=pd, exact_pd=bool(all(pd)) and not bad,
               lam_float=[float(v) for v in lam(x)], resid_before=rmax, time=time.time() - t0)
    if verbose: print('EXACT PRIMAL CERTIFICATE mu>0: %s  (blocks PD %s)  [%.0fs]' % (out['exact_pd'], pd, out['time']), flush=True)
    return out

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('fam', choices=['C1', 'C2']); ap.add_argument('p', nargs=5, type=float, metavar='p')
    ap.add_argument('--allow-outside', action='store_true', help='do not refuse a rationalised source outside the native chart (negative controls)')
    a = ap.parse_args()
    print(exact_primal(a.fam, a.p, allow_outside=a.allow_outside))
