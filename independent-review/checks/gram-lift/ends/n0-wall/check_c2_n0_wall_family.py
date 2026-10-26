"""check_c2_n0_wall_family.py -- the whole open C2 n=0 wall reduced to an explicit 16-parameter Gram family (sympy + Fractions only, no absolute paths).
Usage: python check_c2_n0_wall_family.py   (about 90 seconds; prints PASS lines, asserts everything)

Setting as in check_c2_n0_wall_k2.py: native chart (c,a,b,t,r), d=1, lam=beta=0, parent F = psi^T P(u,v) psi, contact (u,v,W,X,Y,Z)=(c,ac,b/a,-cr,c(b/a-at),1), wall n=0.

PROVED SYMBOLICALLY ON THE WHOLE WALL (parameters c,a,b,r; t=t_w):
 (S1) the five wall weights m,A,D,alpha,C have the short forms (g,h,j,k,B the chart quantities):
        m = r(1-r)(1+a^2 c^2)/(2 g (b+r)),  A = a^2 r/(2h),  D = r(1+c^2) B/(2 a^2 (1+c^2 r) j),
        alpha = c^2 (1-r) B/(2 (1+c^2 r) k),  C = (1-r) B/(2 a^2 c^2 t (b+r)),   g+k = b+r.
 (S2) a UNIVERSAL SUPPORT of 18 features of the k=2 face lift suffices (found numerically, see NOTES.md):
        block 1 (13 integer-power features)   vW,uW,uvW | X,vX,uvX | Y,vY,uY,uvY | Z,vZ,uZ
        block 2 (3)  sqrt(v) X, sqrt(v) Y, sqrt(v) Z          (reduced rank <= 2)
        block 4 (2)  sqrt(uv) W, sqrt(uv) Z                   (reduced rank <= 1)
      The coefficient equations  Z^T G Z == F  are solved EXACTLY with the weights m,C,A,D,alpha symbolic (so valid for every wall point):
      100 Gram entries, 72 determined (explicit, listed), 28 free.  Identity Z^T G(free) Z == F asserted symbolically.
 (S3) G y_c = 0 on the whole wall: the 18 kernel equations (symbolic in c,a,b,r, DomainMatrix over Q(c,a,b,r)) are CONSISTENT and have rank 12
      => the exact Gram family {Z^T G Z = F, G y_c = 0} over the whole wall is an explicit 16-dimensional affine family (28 - 12).
EXACT AT SAMPLE POINTS (S4): at 22 rational wall points (c,a not necessarily squares; includes the 8 points of check_c2_n0_wall_k2.py) a member of this
      18-feature family is positive definite on every reduced block (exact LDL; kernel exactly <y_c>), i.e. the 18-feature lift (not only the 27-feature one) is strictly admitting.
NOT PROVED: positive definiteness of the family for ALL wall points (the 16-dim LMI is only checked at the sample points)."""
import sys, time, json, os
from fractions import Fraction as Fr
import sympy as sp
from sympy.polys.matrices import DomainMatrix
T0 = time.time()
def say(s): print('[%6.1fs] %s' % (time.time() - T0, s), flush=True)
c, a, b, t, r = sp.symbols('c a b t r', positive=True)

def c2_weights(c, a, b, t, r):
    z = c*c; Ac = 1 + z*r; h = b - a*a*t; g = (1 + z)*b - (1 - r + a*a*z*t); j = 1 + z*r - t
    Dl = ((1 + z)*z*b*b + 2*z*(1 + z)*b*r + (1 + z)*(1 + a*a*z)*b + z*z*(1 - a*a)*r*r + 2*z*(1 + a*a*z)*r + 1 + a*a*z)
    N = (1 + z)*(b + r)*(b + a*a*(1 - t)) - a*a*(1 + a*a*z)*(1 - r)*t
    return dict(n=((b + r)*Ac*(Ac + a*a*z*t) - r*Dl)/(2*r*Dl), m=(1 + a*a*z)*(1 - r)*Ac*(Ac + a*a*z*t)/(2*g*Dl), C=Ac*N/(2*a*a*t*Dl),
                A=a*a*Ac*(b + r)*(Ac + a*a*z*t)/(2*h*Dl), D=(1 + z)*(b + r)*(b - a*a*z*r)*(Ac + a*a*z*t)/(2*a*a*j*Dl),
                al=z*z*(b + r)*N/(2*(1 - z*h)*Dl))
Wsym = c2_weights(c, a, b, t, r)
tw = sp.solve(sp.numer(sp.together(Wsym['n'])), t)[0]
Wwall = {k_: sp.cancel(sp.together(v_.subs(t, tw))) for k_, v_ in Wsym.items()}
assert Wwall['n'] == 0
# ------------------------------------------------------------------------------------------------ (S1)
z = c*c; Ac = 1 + z*r
h = sp.cancel(b - a*a*tw); g = sp.cancel((1 + z)*b - (1 - r + a*a*z*tw)); j = sp.cancel(1 + z*r - tw); k = sp.cancel(1 - z*h); B = b - a*a*z*r
assert sp.cancel(g + k - (b + r)) == 0
nice = dict(m=r*(1 - r)*(1 + a*a*z)/(2*g*(b + r)), A=a*a*r/(2*h), D=r*(1 + z)*B/(2*a*a*Ac*j), al=z*(1 - r)*B/(2*Ac*k), C=(1 - r)*B/(2*a*a*z*tw*(b + r)))
for k_, v_ in nice.items(): assert sp.cancel(sp.together(Wwall[k_] - v_)) == 0, k_
say('PASS (S1): on the wall  m = r(1-r)(1+a^2c^2)/(2g(b+r)), A = a^2 r/(2h), D = r(1+c^2)B/(2a^2(1+c^2 r)j), alpha = c^2(1-r)B/(2(1+c^2 r)k), C = (1-r)B/(2a^2c^2 t(b+r)); g+k=b+r  (exact)')
# ------------------------------------------------------------------------------------------------ (S2)
u, v, W, X, Y, Z = sp.symbols('u v W X Y Z'); Us, Vs = sp.symbols('Us Vs'); PSI = [W, X, Y, Z]
m_, C_, A_, D_, al_ = sp.symbols('m C A D alpha', positive=True)
def parent(u, v, W_, X_, Y_, Z_, n, m, C, A, D, al):
    Du, Dv = 1 + u**2, 1 + v**2
    return (4*n*Du*Dv*X_**2 + 4*m*(v*W_ - X_ + u*v*Y_ - u*Z_)**2 + 4*C*Dv*(u*W_ - Y_)**2 + 4*A*Du*Dv*Y_**2
            + (4*D + 2)*(u*W_ + u*v*X_ - Y_ - v*Z_)**2 + (4*al + 2)*Du*(v*Y_ - Z_)**2 + 2*Dv*(X_ + u*Z_)**2 - 2*Du*Dv*Z_**2)
Fsym = sp.expand(parent(Us**2, Vs**2, W, X, Y, Z, 0, m_, C_, A_, D_, al_))           # n = 0 on the wall
names = 'WXYZ'
def nm(f): i, j, l = f; return {0: '', 1: 's', 2: 'u'}[i] + {0: '', 1: 't', 2: 'v'}[j] + names[l]
BL = [[(0,2,0),(2,0,0),(2,2,0),(0,0,1),(0,2,1),(2,2,1),(0,0,2),(0,2,2),(2,0,2),(2,2,2),(0,0,3),(0,2,3),(2,0,3)],
      [(0,1,1),(0,1,2),(0,1,3)],
      [(1,1,0),(1,1,3)]]
def feat(f): i, j, l = f; return Us**i*Vs**j*PSI[l]
Gs = {}
for bi, blk in enumerate(BL):
    for a1 in range(len(blk)):
        for b1 in range(a1, len(blk)): Gs[(bi, a1, b1)] = sp.Symbol('g_%s_%s' % (nm(blk[a1]), nm(blk[b1])))
ZGZ = sum(gsym*feat(BL[bi][a1])*feat(BL[bi][b1])*(1 if a1 == b1 else 2) for (bi, a1, b1), gsym in Gs.items())
diff = sp.Poly(sp.expand(ZGZ - Fsym), Us, Vs, W, X, Y, Z)
unk = list(Gs.values())
sol = sp.solve([co for _, co in diff.terms()], unk, dict=True)[0]
free = [gsym for gsym in unk if gsym not in sol]
assert len(unk) == 100 and len(sol) == 72 and len(free) == 28
entry = {key: sp.expand(sol.get(gsym, gsym)) for key, gsym in Gs.items()}
ZGZ_sol = sum(entry[key]*feat(BL[key[0]][key[1]])*feat(BL[key[0]][key[2]])*(1 if key[1] == key[2] else 2) for key in Gs)
assert sp.expand(ZGZ_sol - Fsym) == 0
say('PASS (S2): 18-feature Gram family G(free) with Z^T G Z == F identically in (m,C,A,D,alpha,free); 100 entries, 72 determined, 28 free: %s' % [str(f_) for f_ in free])
say('   determined diagonals: ' + ', '.join('%s=%s' % (nm(BL[0][i]) , entry[(0, i, i)]) for i in range(13) if not (entry[(0, i, i)].free_symbols & set(free))))
# ------------------------------------------------------------------------------------------------ (S3)
W_c, X_c, Y_c, Z_c = b/a, -c*r, c*(b/a - a*tw), sp.Integer(1)
psi_c = [W_c, X_c, Y_c, Z_c]
def yvec(bi):
    blk = BL[bi]
    if bi == 0: return [c**(f[0]//2)*(a*c)**(f[1]//2)*psi_c[f[2]] for f in blk]        # integer-power features: u^(i/2) v^(j/2) with i,j in {0,2}
    return [psi_c[f[2]] for f in blk]                                                 # odd blocks: the common factor sqrt(v_c) resp. sqrt(u_c v_c) drops out of G y = 0
wsub = {m_: Wwall['m'], C_: Wwall['C'], A_: Wwall['A'], D_: Wwall['D'], al_: Wwall['al']}
rowsM, rhsV = [], []
for bi, blk in enumerate(BL):
    ys = yvec(bi)
    for f in range(len(blk)):
        e = sp.expand(sum(entry[(bi, min(f, gg), max(f, gg))]*ys[gg] for gg in range(len(blk))))
        coefs = [e.coeff(fr) for fr in free]
        cst = sp.expand(e - sum(co*fr for co, fr in zip(coefs, free)))
        rowsM.append([sp.cancel(sp.together(x.subs(wsub))) for x in coefs]); rhsV.append(sp.cancel(sp.together(cst.subs(wsub))))
K = sp.QQ.frac_field(c, a, b, r)
Mdm = DomainMatrix([[K.from_sympy(x) for x in row] for row in rowsM], (len(rowsM), len(free)), K)
Adm = DomainMatrix([[K.from_sympy(x) for x in row] + [K.from_sympy(-bv)] for row, bv in zip(rowsM, rhsV)], (len(rowsM), len(free) + 1), K)
rkM = Mdm.rank(); rkA = Adm.rank()
assert rkM == rkA == 12, (rkM, rkA)
say('PASS (S3): G y_c = 0 (18 equations, symbolic over Q(c,a,b,r), t=t_w) is consistent with rank 12 => explicit 28-12 = 16-dimensional affine Gram family on the whole wall')
# ------------------------------------------------------------------------------------------------ (S4)
def solve_exact(Aeq, rhs, x0):
    m1, n1 = len(Aeq), len(Aeq[0]); piv = []; basis = []
    for i in range(m1):
        row = list(Aeq[i]) + [rhs[i]]
        for (pc, prow) in piv:
            if row[pc] != 0:
                f = row[pc]/prow[pc]; row = [x - f*y for x, y in zip(row, prow)]
        nz = [q for q in range(n1) if row[q] != 0]
        if nz: piv.append((nz[0], row)); basis.append(i)
        else: assert row[n1] == 0, 'inconsistent kernel equations'
    Ab = [Aeq[i] for i in basis]; bb = [rhs[i] for i in basis]; mb = len(Ab)
    res = [bb[i] - sum(Ab[i][q]*x0[q] for q in range(n1)) for i in range(mb)]
    Mm = [[sum(Ab[i][q]*Ab[jj][q] for q in range(n1)) for jj in range(mb)] + [res[i]] for i in range(mb)]
    for i in range(mb):
        p_ = next(q for q in range(i, mb) if Mm[q][i] != 0); Mm[i], Mm[p_] = Mm[p_], Mm[i]; pv = Mm[i][i]; Mm[i] = [x/pv for x in Mm[i]]
        for q in range(mb):
            if q != i and Mm[q][i] != 0: f = Mm[q][i]; Mm[q] = [x - f*y for x, y in zip(Mm[q], Mm[i])]
    yy = [Mm[i][mb] for i in range(mb)]
    return [x0[q] + sum(Ab[i][q]*yy[i] for i in range(mb)) for q in range(n1)]
def ldl_pd(M):
    n1 = len(M); A1 = [row[:] for row in M]; piv = []
    for k1 in range(n1):
        if A1[k1][k1] <= 0: return False, piv
        piv.append(A1[k1][k1])
        for i in range(k1 + 1, n1):
            f = A1[i][k1]/A1[k1][k1]
            if f:
                for jj in range(k1 + 1, n1): A1[i][jj] -= f*A1[k1][jj]
    return True, piv
Wfun = sp.lambdify((c, a, b, r), [Wwall[k_] for k_ in ('m', 'C', 'A', 'D', 'al')], 'sympy')
twf = sp.lambdify((c, a, b, r), tw, 'sympy')
constf = {}; csf = {}
for key, e in entry.items():
    cs_ = [e.coeff(fr) for fr in free]; cst = sp.expand(e - sum(co*fr for co, fr in zip(cs_, free)))
    constf[key] = sp.lambdify((m_, C_, A_, D_, al_), cst, 'sympy'); csf[key] = [Fr(int(x.p), int(x.q)) for x in cs_]
def certify(pt, p0):
    cv, av, bv, rv = [sp.Rational(x) for x in pt]; tv = twf(cv, av, bv, rv); wv = Wfun(cv, av, bv, rv)
    assert tv > 0 and all(x > 0 for x in wv)
    zv = cv*cv; hv = bv - av**2*tv
    assert hv > 0 and (1 + zv)*bv - (1 - rv + av**2*zv*tv) > 0 and 1 + zv*rv - tv > 0 and 1 - zv*hv > 0 and bv - av**2*zv*rv > 0 and 0 < rv < 1
    psi = [bv/av, -cv*rv, cv*(bv/av - av*tv), sp.Integer(1)]
    ent = {key: [Fr(int(sp.Rational(constf[key](*wv)).p), int(sp.Rational(constf[key](*wv)).q)), csf[key]] for key in entry}
    ys = []
    for bi, blk in enumerate(BL):
        if bi == 0: yy = [cv**(f[0]//2)*(av*cv)**(f[1]//2)*psi[f[2]] for f in blk]
        else: yy = [psi[f[2]] for f in blk]
        ys.append([Fr(int(sp.Rational(q).p), int(sp.Rational(q).q)) for q in yy])
    Aeq = []; rhs = []
    for bi, blk in enumerate(BL):
        for f in range(len(blk)):
            row = [Fr(0)]*28; cst = Fr(0)
            for gg in range(len(blk)):
                cs0, cs_ = ent[(bi, min(f, gg), max(f, gg))]
                cst += cs0*ys[bi][gg]
                for q in range(28): row[q] += cs_[q]*ys[bi][gg]
            Aeq.append(row); rhs.append(-cst)
    p = solve_exact(Aeq, rhs, [Fr(x) for x in p0])
    Gb = []; minpiv = []
    for bi, blk in enumerate(BL):
        n1 = len(blk); G = [[Fr(0)]*n1 for _ in range(n1)]
        for a1 in range(n1):
            for b1 in range(a1, n1):
                cs0, cs_ = ent[(bi, a1, b1)]; val = cs0 + sum(cs_[q]*p[q] for q in range(28)); G[a1][b1] = G[b1][a1] = val
        assert all(sum(G[i][q]*ys[bi][q] for q in range(n1)) == 0 for i in range(n1))         # G y_c = 0 exactly
        yk = max(range(n1), key=lambda q: abs(ys[bi][q])); cols = [q for q in range(n1) if q != yk]
        N = [[Fr(0)]*len(cols) for _ in range(n1)]
        for ci, q in enumerate(cols): N[q][ci] = Fr(1); N[yk][ci] = -ys[bi][q]/ys[bi][yk]
        S = [[sum(N[i][a1]*G[i][jj]*N[jj][b1] for i in range(n1) for jj in range(n1) if N[i][a1] != 0 and N[jj][b1] != 0) for b1 in range(len(cols))] for a1 in range(len(cols))]
        ok, piv = ldl_pd(S); assert ok, ('reduced block not PD', bi, pt)
        minpiv.append(float(min(piv))); Gb.append(G)
    # independent identity check at the point:  Z^T G Z == F  (sympy)
    ZG = sum(sp.Rational(Gb[bi][a1][b1].numerator, Gb[bi][a1][b1].denominator)*feat(BL[bi][a1])*feat(BL[bi][b1])*(1 if a1 == b1 else 2) for bi in range(3) for a1 in range(len(BL[bi])) for b1 in range(a1, len(BL[bi])))
    Fpt = parent(Us**2, Vs**2, W, X, Y, Z, 0, *wv)
    assert sp.expand(ZG - Fpt) == 0
    return minpiv
DATA_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data18.json')
if os.path.exists(DATA_FILE):
    D = json.load(open(DATA_FILE)); assert D['free'] == [str(f_) for f_ in free]
    for pt, p0, mg in D['pts']:
        t1 = time.time(); mp = certify(pt, [Fr(x) for x in p0])
        say('PASS (S4): wall point (c,a,b,r)=%s: 18-feature Gram, Z^T G Z == F, G y_c == 0, 3 reduced blocks PD (min pivots %s) [%.0fs]' % (tuple(pt), ['%.2g' % x for x in mp], time.time() - t1))
say('DONE (partial).  NOT proved: PD of the 16-dim family at every wall point.')
