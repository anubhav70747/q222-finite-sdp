"""check_c2_c_limits.py -- exact (sympy) proofs for the limit problems of the two C2 ends c -> infinity and c -> 0
of the native chart (c,a,b,t,r), d=1, lam=beta=0.  The finite-c perturbation step is NOT proved here.
Run:  python check_c2_c_limits.py   (exact symbolic checks only; prints PASS lines and asserts every identity)
Optional numerical sanity check: add --numeric-sanity. No solver is used.

What is proved here (all identities are exact rational-function / polynomial identities, all positivity statements are
'every numerator and denominator has only positive coefficients in positive coordinates of the open face'):

 c -> infinity  (b = a^2 t + eta/c^2, r = rho/c^2; face = {0<rho<t, rho>t-1, 0<eta<1, a^2 t+eta>1}, chart parameters a,t,eta,rho)
  (I1) exact rescaled parent  F~ = c^-2 F_c(cU,cV,W,X'/c,Y'/c,Z)  is an explicit polynomial in eps = c^-2 (no limit taken)
  (I2) limit weights n,m,A,D,alpha (finite), C = C2 eps^2 + O(eps^3), C2 > 0 on the face
  (I3) limit form  F_0 = F~|_{eps=0,C=0}
  (I4) F_0 = cX (UVX+LX)^2 + cY (UVY+LY)^2 + p2 (VW - a^2 t UZ)^2 + q2 (UW - t VZ)^2   (exact; cX,cY,p2,q2 > 0 explicit rational functions)
  (I5) k=2 Gram of F_0 (features UW,UVX,VZ,VW,UZ,UVY | sqrt(UV)W, sqrt(UV)Z) with one free parameter kappa:
       F_0 = Q_kappa(y) + kappa*UV (W - a t Z)^2 ;  Q_kappa = cX e_X^2 + cY e_Y^2 + [p2 u2^2 + q2 u1^2 - kappa (u1 + t w)(u2 - a t w)]
       and Q_kappa is positive definite on R^6/<y_contact> iff 0 < kappa < kappa_max = 4 a p2 q2/(q2 + a^2 p2)   (explicit LDL)
       => with kappa = kappa_max/2 the k=2 limit Gram has kernel EXACTLY the contact (rank 5 + 1 of 6 + 2 features), on the whole open face.
 c -> 0  (b,a,t,r fixed;  face = {0<r<1, 0<t<1, b+r>1, b>a^2 t}) : the same chain (Z1)-(Z5).
"""
import argparse
import sympy as sp, time
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--numeric-sanity', action='store_true', help='also sample 4000 face points; not part of the proof')
args = parser.parse_args()
T0 = time.time()
def say(s): print('[%6.1fs] %s' % (time.time() - T0, s), flush=True)

# =====================================================================================================  C -> INFINITY
a, t, eta, rho = sp.symbols('a t eta rho', positive=True)
c, eps = sp.symbols('c epsilon', positive=True)
U, V, W, X, Y, Z = sp.symbols('U V W X Y Z')
n, m, C, A, D, al = sp.symbols('n m C A D alpha')


def qf_zero(expr, psi=(W, X, Y, Z)):
    """True iff the form `expr` (quadratic in psi, polynomial in U,V, rational in the parameters) vanishes identically: every psi-Hessian entry cancels."""
    for i in range(4):
        for j in range(i, 4):
            e = sp.diff(expr, psi[i], psi[j])
            if sp.cancel(sp.together(e)) != 0: return False
    return True

def parent(u, v, W_, X_, Y_, Z_, n, m, C, A, D, al):      # d=1, lam=beta=0 parent form of the repository (gram_lift.form)
    Du, Dv = 1 + u**2, 1 + v**2
    return (4*n*Du*Dv*X_**2 + 4*m*(v*W_ - X_ + u*v*Y_ - u*Z_)**2 + 4*C*Dv*(u*W_ - Y_)**2 + 4*A*Du*Dv*Y_**2
            + (4*D + 2)*(u*W_ + u*v*X_ - Y_ - v*Z_)**2 + (4*al + 2)*Du*(v*Y_ - Z_)**2 + 2*Dv*(X_ + u*Z_)**2 - 2*Du*Dv*Z_**2)

# (I1) -----------------------------------------------------------------------------------------------
Ft_formula = (4*n*(U**2 + eps)*(V**2 + eps)*X**2 + 4*A*(U**2 + eps)*(V**2 + eps)*Y**2 + 4*m*(V*W + U*V*Y - U*Z - eps*X)**2
              + (4*D + 2)*(U*W + U*V*X - V*Z - eps*Y)**2 + (4*al + 2)*(U**2 + eps)*(V*Y - Z)**2
              + 4*C/eps*(V**2 + eps)*(U*W - eps*Y)**2 + 2*(V**2 + eps)*(2*U*X*Z + eps*X**2 - Z**2))
Ft_def = parent(c*U, c*V, W, X/c, Y/c, Z, n, m, C, A, D, al)/c**2
assert sp.simplify((Ft_def - Ft_formula).subs(eps, 1/c**2)) == 0
say('PASS (I1): c^-2 F_c(cU,cV,W,X/c,Y/c,Z) = explicit eps-form (eps=c^-2, 4C c^2 = 4C/eps)')

# (I2) -----------------------------------------------------------------------------------------------
def c2_weights(c, a, b, t, r):
    z = c*c; Ac = 1 + z*r; h = b - a*a*t; g = (1 + z)*b - (1 - r + a*a*z*t); j = 1 + z*r - t
    Dl = ((1 + z)*z*b*b + 2*z*(1 + z)*b*r + (1 + z)*(1 + a*a*z)*b + z*z*(1 - a*a)*r*r + 2*z*(1 + a*a*z)*r + 1 + a*a*z)
    N = (1 + z)*(b + r)*(b + a*a*(1 - t)) - a*a*(1 + a*a*z)*(1 - r)*t
    return dict(n=((b + r)*Ac*(Ac + a*a*z*t) - r*Dl)/(2*r*Dl), m=(1 + a*a*z)*(1 - r)*Ac*(Ac + a*a*z*t)/(2*g*Dl), C=Ac*N/(2*a*a*t*Dl),
                A=a*a*Ac*(b + r)*(Ac + a*a*z*t)/(2*h*Dl), D=(1 + z)*(b + r)*(b - a*a*z*r)*(Ac + a*a*z*t)/(2*a*a*j*Dl),
                al=z*z*(b + r)*N/(2*(1 - z*h)*Dl))
def limit_weights(a, t, eta, rho):
    return dict(n=(t - rho)/(2*rho*(t + 1)), m=(1 + rho)/(2*(1 + t)*(a*a*t + eta - 1)), A=a*a*t*(1 + rho)/(2*eta*(1 + t)),
                D=t*(t - rho)/(2*(t + 1)*(1 + rho - t)), al=(a*a*t*(1 + rho) + eta*(t + 1) + rho - t)/(2*(1 - eta)*(t + 1)))
Wt = c2_weights(c, a, a**2*t + eta/c**2, t, rho/c**2)
LW = limit_weights(a, t, eta, rho)
Wte = {k: sp.cancel(sp.together(v_.subs(c, 1/sp.sqrt(eps)))) for k, v_ in Wt.items()}
for k in ('n', 'm', 'A', 'D', 'al'):
    num, den = sp.fraction(Wte[k]); assert sp.simplify(den.subs(eps, 0)) != 0
    assert sp.simplify(Wte[k].subs(eps, 0) - LW[k]) == 0, k
Cs = Wte['C']; assert sp.simplify(Cs.subs(eps, 0)) == 0 and sp.simplify(sp.diff(Cs, eps).subs(eps, 0)) == 0
C2 = sp.factor(sp.simplify(sp.diff(Cs, eps, 2).subs(eps, 0)/2))
assert sp.simplify(C2 - (rho + 1)*(a**2*rho*t + a**2*t + eta*t + eta + rho - t)/(2*a**4*t**2*(t + 1))) == 0
say('PASS (I2): n,m,A,D,alpha -> limit weights (finite, denominators nonzero at eps=0); C = C2 eps^2 + O(eps^3), C2 = %s' % C2)

# (I3) -----------------------------------------------------------------------------------------------
F0 = (4*LW['n']*U**2*V**2*X**2 + 4*LW['A']*U**2*V**2*Y**2 + 4*LW['m']*(V*W + U*V*Y - U*Z)**2 + (4*LW['D'] + 2)*(U*W + U*V*X - V*Z)**2
      + (4*LW['al'] + 2)*U**2*(V*Y - Z)**2 + 4*U*V**2*X*Z - 2*V**2*Z**2)
F0e = Ft_formula.subs({n: LW['n'], m: LW['m'], A: LW['A'], D: LW['D'], al: LW['al']}).subs(C, 0).subs(eps, 0)
assert qf_zero(F0e - F0)
say('PASS (I3): F~ at eps=0 with the limit weights and C-term dropped is F_0')

# (I4) -----------------------------------------------------------------------------------------------
xs = sp.Symbol('xs')
cX = 2*t*(-2*rho + t - 1)/(rho*(t + 1)*(-rho + t - 1))
cY = -2*a**2*t*(rho + 1)*(a**2*t + 2*eta - 1)/(eta*(eta - 1)*(t + 1)*(a**2*t + eta - 1))
LX = -rho*((rho + 1)*U*W + t*(rho - t)*V*Z)/(t*(t - 2*rho - 1))
LY = -eta*((eta - 1)*V*W + a**2*t*(a**2*t + eta)*U*Z)/(a**2*t*(a**2*t + 2*eta - 1))
p2 = 2*(rho + 1)*(a**2*t + eta)/(a**2*t*(t + 1)*(a**2*t + 2*eta - 1))
q2 = 2*(rho + 1)*(t - rho)/(t*(t + 1)*(2*rho + 1 - t))
RHS = cX*(U*V*X + LX)**2 + cY*(U*V*Y + LY)**2 + p2*(V*W - a**2*t*U*Z)**2 + q2*(U*W - t*V*Z)**2
assert qf_zero(F0 - RHS)
say('PASS (I4): F_0 = cX (UVX+LX)^2 + cY (UVY+LY)^2 + p2 (VW-a^2 t UZ)^2 + q2 (UW-tVZ)^2 EXACT')

# positivity on the open face: positive coordinates  rho=r, t=r+x/(1+x), eta=1/(1+y), a^2 t + eta - 1 = s  (all of r,x,y,s>0)
r_, x_, y_, s_ = sp.symbols('r_ x_ y_ s_', positive=True)
def has_one_sign(poly, gens):
    cs = sp.Poly(sp.expand(poly), *gens).coeffs()
    return (all(v > 0 for v in cs) or all(v < 0 for v in cs)), (1 if cs[0] > 0 else -1)
def manifestly_positive_inf(expr, a_sq_only=True):
    if sp.simplify(expr.subs(a, -a) + expr) == 0: expr = sp.simplify(expr/a)      # odd in a: a>0 is a factor, test the even part
    """expr is a rational function of a (even or odd), t, eta, rho.  Substitute the positive coordinates; a = sqrt(a2) is treated by a2-substitution
    after splitting a-odd part (expr = e0(a^2) + a e1(a^2) must be a monomial in a times a function of a^2)."""
    a2 = (s_ + 1 - 1/(1 + y_))/(r_ + x_/(1 + x_))
    sub = {t: r_ + x_/(1 + x_), eta: 1/(1 + y_), rho: r_}
    e = sp.together(expr.subs(a, sp.sqrt(sp.Symbol('A2', positive=True)))).subs(sp.Symbol('A2', positive=True), a2)
    e = sp.cancel(sp.together(sp.simplify(e.subs(sub))))
    num, den = sp.fraction(e)
    ok1, s1 = has_one_sign(num, (r_, x_, y_, s_)); ok2, s2 = has_one_sign(den, (r_, x_, y_, s_))
    return ok1 and ok2 and s1*s2 > 0
# the chart conditions 0<rho<t, rho>t-1, 0<eta<1, a^2 t+eta>1 are exactly r_>0, 0<x_/(1+x_)<1, y_>0, s_>0
pos_list = [('cX', cX), ('cY', cY), ('p2', p2), ('q2', q2), ('C2', C2), ('n', LW['n']), ('m', LW['m']), ('A', LW['A']), ('D', LW['D']), ('alpha', LW['al'])]
for nm, ex in pos_list:
    assert manifestly_positive_inf(ex), nm
say('PASS (I4+): cX, cY, p2, q2, C2 and the limit weights n,m,A,D,alpha are > 0 on the whole open face (all coefficients positive in r,x,y,s)')

# (I5) k=2 Gram -------------------------------------------------------------------------------------
kap = sp.Symbol('kappa', positive=True)
a1_, a2_, a3_, b1_, b2_, b3_ = sp.symbols('a1 a2 a3 b1 b2 b3')       # features UW, UVX, VZ, VW, UZ, UVY
eX = a2_ - rho*((rho + 1)*a1_ + t*(rho - t)*a3_)/(t*(t - 2*rho - 1))        # UVX + LX in feature variables
eY = b3_ - eta*((eta - 1)*b1_ + a**2*t*(a**2*t + eta)*b2_)/(a**2*t*(a**2*t + 2*eta - 1))
u1 = a1_ - t*a3_; u2 = b1_ - a**2*t*b2_; w = a3_ - a*b2_
Qk = cX*eX**2 + cY*eY**2 + p2*u2**2 + q2*u1**2 - kap*(a1_ - a*t*b2_)*(b1_ - a*t*a3_)
yfeat = {a1_: U*W, a2_: U*V*X, a3_: V*Z, b1_: V*W, b2_: U*Z, b3_: U*V*Y}
assert qf_zero(Qk.subs(yfeat) + kap*U*V*(W - a*t*Z)**2 - F0)
say('PASS (I5a): F_0 = Q_kappa(UW,UVX,VZ,VW,UZ,UVY) + kappa*(sqrt(UV)(W - a t Z))^2   (polynomial identity in U,V,W,X,Y,Z; kappa free)')
# contact kernel
yc = {a1_: a*t, a2_: -a*rho, a3_: a, b1_: a**2*t, b2_: 1, b3_: eta}
GE = sp.hessian(Qk, [a1_, a2_, a3_, b1_, b2_, b3_])/2
ycv = sp.Matrix([yc[s] for s in (a1_, a2_, a3_, b1_, b2_, b3_)])
assert (GE*ycv).applyfunc(sp.simplify).is_zero_matrix
assert sp.simplify(Qk.subs(yc)) == 0
say('PASS (I5b): G_E(kappa) (6x6) annihilates the contact vector y_c = (a t,-a rho,a,a^2 t,1,eta); the 1x2 odd Gram kappa*(1,-a t)(1,-a t)^T annihilates (W_c,Z_c)=(a t,1)')
# change of variables to the 5 forms eX,eY,u1,u2,w; kernel of these 5 forms is exactly <y_c>
Hk = p2*u2**2 + q2*u1**2 - kap*(u1 + t*w)*(u2 - a*t*w)
assert sp.cancel(sp.together(sp.expand(p2*u2**2 + q2*u1**2 - kap*(a1_ - a*t*b2_)*(b1_ - a*t*a3_) - Hk))) == 0     # (a1-atb2)=u1+tw, (b1-ata3)=u2-atw
M5 = sp.Matrix([[sp.diff(f, s) for s in (a1_, a2_, a3_, b1_, b2_, b3_)] for f in (eX, eY, u1, u2, w)])
assert M5.rank() == 5 and (M5*ycv).applyfunc(sp.simplify).is_zero_matrix
say('PASS (I5c): (eX,eY,u1,u2,w) are 5 independent linear forms on R^6 vanishing exactly on <y_c>;  Q_kappa = cX eX^2 + cY eY^2 + H_kappa(u1,u2,w)')
# LDL of H_kappa:  H = kappa a t^2 (w + (a u1 - u2)/(2 a t))^2 + Aq u1^2 - kappa/2 u1 u2 + Bq u2^2  (Schur in w), then 2x2 LDL
g1, g2 = t, a*t
Hdet = sp.Matrix([[q2, -kap/2, kap*g2/2], [-kap/2, p2, -kap*g1/2], [kap*g2/2, -kap*g1/2, kap*g1*g2]])
uvw = sp.Matrix([u1, u2, w])
assert sp.simplify(sp.expand((uvw.T*Hdet*uvw)[0] - Hk)) == 0
Aq = q2 - kap*g2/(4*g1); Bq = p2 - kap*g1/(4*g2)
schur = sp.expand(Hk - kap*g1*g2*(w + (g2*u1 - g1*u2)/(2*g1*g2))**2)
assert sp.simplify(schur - (Aq*u1**2 - kap/2*u1*u2 + Bq*u2**2)) == 0
assert sp.simplify(sp.expand(Aq*u1**2 - kap/2*u1*u2 + Bq*u2**2 - Aq*(u1 + (-kap/4)/Aq*u2)**2 - (Bq - kap**2/(16*Aq))*u2**2)) == 0
say('PASS (I5d): exact LDL  H_kappa = kappa a t^2 (w + (a u1 - u2)/(2 a t))^2 + Aq (u1 - kappa u2/(4 Aq))^2 + (Bq - kappa^2/(16 Aq)) u2^2 ,  Aq = q2 - kappa a/4, Bq = p2 - kappa/(4a)')
kmax = 4*a*p2*q2/(q2 + a**2*p2)
kk = kmax/2
Aq_ = sp.simplify(Aq.subs(kap, kk)); Bq_ = sp.simplify(Bq.subs(kap, kk))
assert sp.simplify(Aq_ - q2*(2*q2 + a**2*p2)/(2*(q2 + a**2*p2))) == 0
assert sp.simplify(Bq_ - p2*(q2 + 2*a**2*p2)/(2*(q2 + a**2*p2))) == 0
assert sp.simplify(Aq_*Bq_ - kk**2/16 - p2*q2/2) == 0
assert sp.simplify(Aq.subs(kap, kmax)*Bq.subs(kap, kmax) - kmax**2/16) == 0       # kappa_max is exactly the edge of positive definiteness
say('PASS (I5e): kappa_max = 4 a p2 q2/(q2 + a^2 p2) is the exact edge (Aq Bq = kappa^2/16); at kappa = kappa_max/2: Aq = q2(2q2+a^2 p2)/(2(q2+a^2 p2)), Bq = p2(q2+2a^2 p2)/(2(q2+a^2 p2)), Aq Bq - kappa^2/16 = p2 q2/2 > 0')
for nm, ex in (('kappa/a', kk/a), ('kappa a t^2', kk*a*t**2), ('Aq/q2', Aq_/q2), ('Bq/p2', Bq_/p2)):
    assert manifestly_positive_inf(ex), nm
say('PASS (I5f): kappa = kappa_max/2 > 0 and all pivots (kappa a t^2, Aq, Bq - kappa^2/(16 Aq)) > 0 on the whole open face.  CONSEQUENCE: for every face point the c=infinity limit form F_0 has a k=2 Gram')
say('        (6 + 2 features: blocks (UW,UVX,VZ,VW,UZ,UVY) and (sqrt(UV)W, sqrt(UV)Z)) that is PSD with kernel EXACTLY the contact in each block (strictly feasible in the k=2 limit problem).')
# numerical sanity of the closed form against an independent evaluation (not part of the proof)
if args.numeric_sanity:
    import random
    import numpy as np
    GEf = sp.lambdify((a, t, eta, rho, kap), GE, 'numpy'); kkf = sp.lambdify((a, t, eta, rho), kk, 'numpy')
    random.seed(11); bad = 0
    for _ in range(4000):
        rr = random.uniform(1e-2, 3); w_ = random.uniform(1e-2, 1 - 1e-2); tt = rr + w_; e_ = random.uniform(1e-2, 1 - 1e-2); ss = random.uniform(1e-2, 5)
        # a^2 t + eta - 1 = ss -> a^2 = (ss + 1 - eta)/t, eta = 1 - e_
        aa = ((ss + e_)/tt)**0.5; et = 1 - e_
        ev = np.linalg.eigvalsh(np.array(GEf(aa, tt, et, rr, kkf(aa, tt, et, rr)), float))
        if ev[0] < -1e-9*max(1, ev[-1]) or ev[1] < 1e-12*max(1, ev[-1]): bad += 1
    say('numerical sanity (not a proof): 4000 random face points: G_E(kappa_max/2) PSD with exactly one zero eigenvalue in %d/4000 cases failing' % bad)
    assert bad == 0

# =====================================================================================================  C -> 0
b, t0, r0 = sp.symbols('b t0 r0', positive=True)    # b,t,r of the c=0 chart (a is shared)
tt_, rr_ = sp.symbols('tt rr', positive=True)
c0 = sp.Symbol('c0', positive=True)
def zero_weights_limit():
    nn = b*(1 - rr_)/(2*rr_*(1 + b)); CC = (a*a*(b*(1 - tt_) + rr_ - tt_) + b*(b + rr_))/(2*a*a*tt_*(b + 1)); AA = a*a*(b + rr_)/(2*(b + 1)*(b - a*a*tt_))
    DD = b*(b + rr_)/(2*a*a*(b + 1)*(1 - tt_)); mm = (b*(1 - rr_) - (1 - rr_)**2)/(2*(b**3 + 2*b**2*rr_ - b**2 + b*rr_**2 - b + rr_**2 - 2*rr_ + 1))
    return dict(n=nn, m=mm, C=CC, A=AA, D=DD, al=sp.Integer(0))
LZ = zero_weights_limit()
Wz = c2_weights(c0, a, b, tt_, rr_)
def fin_limit(expr):
    num, den = sp.fraction(sp.cancel(sp.together(expr))); assert sp.simplify(den.subs(c0, 0)) != 0
    return sp.simplify(num.subs(c0, 0)/den.subs(c0, 0))
for k_ in ('n', 'm', 'C', 'A', 'D'):
    assert sp.simplify(fin_limit(Wz[k_]) - LZ[k_]) == 0, k_
al4 = fin_limit(Wz['al']/c0**4)
say('PASS (Z1): c->0 limits of n,m,C,A,D exist (rational, finite, denominators nonzero at c=0) and alpha = c^4 * (finite rational, alpha4 = %s) -> 0' % sp.factor(al4))
# rescaled parent F(cU,cV,W,cX',cY',Z)/c^2 -> F^(0) coefficientwise
Wtz = {k_: v_ for k_, v_ in Wz.items()}
Ftz = sp.expand(parent(c0*U, c0*V, W, c0*X, c0*Y, Z, Wtz['n'], Wtz['m'], Wtz['C'], Wtz['A'], Wtz['D'], Wtz['al'])/c0**2)
Pz = sp.Poly(Ftz, U, V, W, X, Y, Z); F00 = 0
for mon, co in Pz.terms():
    lim = fin_limit(co)
    F00 += sp.simplify(lim)*sp.prod([s_**e_ for s_, e_ in zip((U, V, W, X, Y, Z), mon)])
n_, m_, C_, A_, D_ = [LZ[k_] for k_ in ('n', 'm', 'C', 'A', 'D')]
F00ref = (4*n_*X**2 + 4*m_*(V*W - X - U*Z)**2 + 4*C_*(U*W - Y)**2 + 4*A_*Y**2 + (4*D_ + 2)*(U*W - Y - V*Z)**2 + 2*(X + U*Z)**2 - 2*V**2*Z**2 - 4*V*Y*Z)
F00ref = F00ref.subs({tt_: tt_})
assert qf_zero(F00 - F00ref)
say('PASS (Z2): lim_{c->0} c^-2 F_c(cU,cV,W,cX\',cY\',Z) = F^(0) (all monomials, exact coefficientwise limits)')
# (Z3) Schur decomposition of F^(0):  y=(A1..A3|B1..B3)=(UW,Y,VZ|VW,X,UZ)
muz = sp.Symbol('muz')
A1, A2, A3, B1, B2, B3 = sp.symbols('A1 A2 A3 B1 B2 B3')
Fy0 = (4*C_*(A1 - A2)**2 + 4*A_*A2**2 + (4*D_ + 2)*(A1 - A2 - A3)**2 - 2*A3**2 - 4*A2*A3 + 4*n_*B2**2 + 4*m_*(B1 - B2 - B3)**2 + 2*(B2 + B3)**2)
featz = {A1: U*W, A2: Y, A3: V*Z, B1: V*W, B2: X, B3: U*Z}
assert qf_zero(Fy0.subs(featz) - F00ref)
def comp(Fq, v):
    pl = sp.Poly(sp.expand(Fq), v); c2_ = pl.coeff_monomial(v**2); c1_ = pl.coeff_monomial(v); c0_ = pl.coeff_monomial(1)
    return sp.simplify(c2_), sp.simplify(c1_/(2*c2_)), sp.simplify(c0_ - c1_**2/(4*c2_))
cA2, lA2, R1 = comp(Fy0, A2)           # A2 + lA2  (lA2 linear in A1,A3)
cB2, lB2, R2 = comp(R1, B2)            # B2 + lB2  (lB2 linear in B1,B3)
# the relation multiplier is already absorbed: Fy0 is the polynomial form; R2 must be (after using A1*A3=B1*B3) two squares
P0 = 2*a**2*(b + rr_)*(a**2*tt_ - a**2 - b)/(b*(b + 1)*(2*a**2*tt_ - a**2 - b)); Q0 = -2*(b + rr_)*(rr_ - 1)/(b*(b + 1)*(b + 2*rr_ - 1))
muval = sp.symbols('muval')
R2s = sp.expand(R2 + muval*(A1*A3 - B1*B3))
two = P0*(A1 - b/a**2*A3)**2 + Q0*(B1 - b*B3)**2
sol_mu = sp.solve(sp.Poly(sp.expand(R2s - two), A1, A3, B1, B3).coeffs(), muval, dict=True)
assert len(sol_mu) == 1
assert sp.simplify(sp.expand(R2s.subs(sol_mu[0]) - two)) == 0
say('PASS (Z3): F^(0) = cA2 (Y+lA2)^2 + cB2 (X+lB2)^2 + P0 (UW - (b/a^2) VZ)^2 + Q0 (VW - b UZ)^2   with P0 = %s, Q0 = %s  (exact; multiplier of A1A3-B1B3 solved uniquely: %s)' % (sp.factor(P0), sp.factor(Q0), sp.factor(sol_mu[0][muval])))
# positivity on the c=0 face: r = p/(1+p), t = q/(1+q), b = 1 - r + s, a^2 t = b k/(1+k)   (p,q,s,k > 0)
p_, q_, s2_, k_ = sp.symbols('p_ q_ s2_ k_', positive=True)
def manifestly_positive_zero(expr):
    if sp.simplify(expr.subs(a, -a) + expr) == 0: expr = sp.simplify(expr/a)
    A2s = sp.Symbol('A2s', positive=True)
    bb = 1 - p_/(1 + p_) + s2_; a2 = bb*k_/((1 + k_)*(q_/(1 + q_)))
    e = expr.subs(a, sp.sqrt(A2s)).subs(A2s, a2).subs({b: bb, tt_: q_/(1 + q_), rr_: p_/(1 + p_)})
    e = sp.cancel(sp.together(sp.simplify(e))); num, den = sp.fraction(e)
    ok1, s1 = has_one_sign(num, (p_, q_, s2_, k_)); ok2, s2s = has_one_sign(den, (p_, q_, s2_, k_)); return ok1 and ok2 and s1*s2s > 0
ge1, ge2 = b/a**2, b/a            # (g1,g2) of the c=0 end:  (a1 - (b/a) b3)=u1+g1 w,  (b1-(b/a) a3) = u2 - g2 w,  w = A3 - a B3
kmaxz = 4*ge1*ge2*P0*Q0/(P0*ge1**2 + Q0*ge2**2); kkz = kmaxz/2
for nm, ex in (('cA2', cA2), ('cB2', cB2), ('P0', P0), ('Q0', Q0), ('kappa0', kkz), ('kappa0*g1*g2', kkz*ge1*ge2),
               ('n', n_), ('m', m_), ('C', C_), ('A', A_), ('D', D_)):
    assert manifestly_positive_zero(ex), nm
say('PASS (Z4): cA2, cB2, P0, Q0, kappa0 = kappa_max/2 and the limit weights n,m,C,A,D are > 0 on the whole open c=0 face')
# (Z5) the k=2 Gram at c=0:  F^(0) = Q0form(y) + kappa*UV (W - (b/a) Z)^2 and its LDL
u1z = A1 - b/a**2*A3; u2z = B1 - b*B3; wz = A3 - a*B3; kz = sp.Symbol('kz', positive=True)
Qz = cA2*(A2 + lA2)**2 + cB2*(B2 + lB2)**2 + P0*u1z**2 + Q0*u2z**2 - kz*(A1 - b/a*B3)*(B1 - b/a*A3)
assert qf_zero(Qz.subs(featz) + kz*U*V*(W - b/a*Z)**2 - F00ref)
say('PASS (Z5a): F^(0) = Q_kappa(UW,Y,VZ,VW,X,UZ) + kappa*(sqrt(UV)(W - (b/a) Z))^2  (polynomial identity; contact (W_c,Z_c) = (b/a,1))')
ycz = {A1: b/a, A2: b/a - a*tt_, A3: a, B1: b, B2: -rr_, B3: 1}
assert sp.simplify(Qz.subs(ycz)) == 0
GEz = sp.hessian(Qz, [A1, A2, A3, B1, B2, B3])/2
ycvz = sp.Matrix([ycz[s] for s in (A1, A2, A3, B1, B2, B3)])
assert (GEz*ycvz).applyfunc(sp.simplify).is_zero_matrix
M5z = sp.Matrix([[sp.diff(f, s) for s in (A1, A2, A3, B1, B2, B3)] for f in (A2 + lA2, B2 + lB2, u1z, u2z, wz)])
assert M5z.rank() == 5 and (M5z*ycvz).applyfunc(sp.simplify).is_zero_matrix
Hz = P0*u1z**2 + Q0*u2z**2 - kz*(u1z + ge1*wz)*(u2z - ge2*wz)
assert sp.cancel(sp.together(sp.expand(P0*u1z**2 + Q0*u2z**2 - kz*(A1 - b/a*B3)*(B1 - b/a*A3) - Hz))) == 0
schurz = sp.expand(Hz - kz*ge1*ge2*(wz + (ge2*u1z - ge1*u2z)/(2*ge1*ge2))**2)
Az = P0 - kz*ge2/(4*ge1); Bz = Q0 - kz*ge1/(4*ge2)
assert sp.simplify(schurz - (Az*u1z**2 - kz/2*u1z*u2z + Bz*u2z**2)) == 0
assert sp.simplify(Az.subs(kz, kmaxz)*Bz.subs(kz, kmaxz) - kmaxz**2/16) == 0
Az_ = sp.simplify(Az.subs(kz, kkz)); Bz_ = sp.simplify(Bz.subs(kz, kkz))
assert sp.simplify(Az_ - P0*(2*P0*ge1**2 + Q0*ge2**2)/(2*(P0*ge1**2 + Q0*ge2**2))) == 0
assert sp.simplify(Bz_ - Q0*(P0*ge1**2 + 2*Q0*ge2**2)/(2*(P0*ge1**2 + Q0*ge2**2))) == 0
assert sp.simplify(Az_*Bz_ - kkz**2/16 - P0*Q0/2) == 0
say('PASS (Z5b): same LDL at c=0: kappa_max = 4 g1 g2 P0 Q0/(P0 g1^2 + Q0 g2^2), g1=b/a^2, g2=b/a; at kappa_max/2 the pivots are P0(2P0 g1^2+Q0 g2^2)/(2(P0 g1^2+Q0 g2^2)), ..., all > 0 => the k=2 limit Gram at c=0 has kernel exactly the contact, on the whole open c=0 face')
say('ALL ASSERTIONS PASSED')
