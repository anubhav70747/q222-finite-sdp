# Strict native C1 and the K364 fixed-grid question: parity-resolved contact jets

27 October 2026. **Exact status:** a source-wide nondegenerate contact Gram and a *single* source-only polynomial skew matrix that is PSD on a two-sided open Bob-angle neighbourhood of the native contact are constructed. The unrestricted full parent is matched coefficientwise for every Alice angle and four-state vector. **H1, global PSD for every Bob angle in K364, is not proved and no K364 separator is claimed.** The current result is a new contact-local sufficient mechanism, not a global admission theorem. No numerical sweeps, parameter searches, or GitHub changes are used.

## 1. Exact domain, parent and fixed features

Let `theta=(k,f,g,a,t)` be a **genuine strict real oriented C1 source** in SD1/NC1, with `0<k<=1`, `k,f,g,a,t>0`, `r=1-af^2 in (0,1)`, `theta_B=tg^2 in (0,1)`, `L_c=r>0`, `H_c=1+(1-a)f^2>0`, `Q_c=1-k^2g^2+ak^2f^2g^2>0`, `M=1+k^2r-t>0`, `J=1+k^2-k^2g^2t-a>0`, all other original source denominators/gaps positive, and all six retained weights `m,C,A,D,alpha,beta>0`, with `d>0` and `lambda=n=0`. Normalize `C=2` through the original positive-source scale. We do not assume `a=t=1`, and no source/testing rate condition enters.

The actual stationary contact, as provided by the native dictionary, is

```
u_* = k, v_* = k f g,
psi_* = (W_*,X_*,Y_*,Z_*)
      = (g/f,-k(1-tg^2),kg(1-af^2)/f,1).
```

All four state coordinates are nonzero; `W_*,Y_*,Z_*>0`, `X_*<0`. Put `D_u=1+u^2`, `D_v=1+v^2`, and retain the **whole** C1 polynomial

```
F = 4m(vW-X+uvY-uZ)^2 + 4C D_v(uW-Y)^2 + 4A D_u D_v Y^2
  + (4D+2d)(uW+uvX-Y-vZ)^2
  + (4alpha+2d)D_u(vY-Z)^2 +(4beta+2d)D_v(X+uZ)^2
  - 2d D_u D_v Z^2 .                                        (1)
```

Exactly extract `F=psi^T[P0(v)+uP1(v)+u^2P2(v)]psi`, with `Pj=Pj^T` and source-linear entries. All three matrices are polynomials of degree <=2 in v; the full `D_v F` has Bob degree <=4. For any source-only skew `S(v)` of Bob degree <=4, define

```
H(v) = [[D_v P0(v), D_v P1(v)/2 + S(v)],
        [D_v P1(v)/2 - S(v), D_v P2(v)]].                    (2)
```

Then `H(v)=H(v)^T`, and **for every independent u,v and every unrestricted real state psi**,

```
(psi,u psi)^T H(v) (psi,u psi) = D_v F(u,v;psi).             (3)
```

When H(v) is PSD for *all* v>=0, the previously reviewed univariate matrix-SOS factorization supplies a coefficientwise positive Gram in the **same fixed 40 rows** with `u` exponents 0,1 and Bob exponents 0,1/2,1,3/2,2, all embedded in the fixed `Z_(6,12)` 364-feature library. Equation (3) alone does **not** imply that global positivity.

## 2. Exact 4+4 parity reduction and the final global two-Schur obligation

Group the upper state coordinates into `w=(W,X)` and lower into `y=(Y,Z)`; use indices `0,1,2,3` for W,X,Y,Z. The actual full-parent matrices have `P0(v)=diag(P0_w(v),P0_y(v))`, `P2(v)=diag(P2_w(v),P2_y(v))`, and `P1(v)=[[0,E(v)],[E(v)^T,0]]`. This is a direct polynomial identity from (1), not a state slice.

Choose the **cross-parity source-only skew**

```
S(v) = [[0,T(v)],[-T(v)^T,0]],  T(v) in Mat_2(R[v]).
```

In the fixed 8-coordinate order `(W0,X0,Y0,Z0;W1,X1,Y1,Z1)`, H(v) is exactly the direct sum, after permutation, of

```
H_A(v) = [[A0(v), L(v)+T(v)],
          [(L(v)+T(v))^T, A2(v)]],

H_B(v) = [[B0(v), L(v)^T-T(v)^T],
          [L(v)-T(v), B2(v)]],                                (4)
```

with `A0=D_v P0_w`, `B0=D_v P0_y`, `A2=D_v P2_y`, `B2=D_v P2_w` and `L=D_v E/2`. These two blocks contain every state normal and every state cross term. For v>0 both `A0` and `B0` are **strictly** PD: `det P0_w=8m(2beta+d)D_v v^2>0`, its first diagonal is `4m v^2>0`, and `P0_y` is a positive sum of `(C+A)D_v Y^2`, `alpha(vY-Z)^2`, and the other positive terms in the previously verified u=0 parent.

Therefore the exact global sufficient condition for **this cross-parity 40-feature ansatz** is

```
R_A(v) := A2(v) -(L(v)+T(v))^T A0(v)^(-1)(L(v)+T(v)) >= 0,
R_B(v) := B2(v) -(L(v)-T(v)) B0(v)^(-1)(L(v)-T(v))^T >= 0,
                                for ALL v>0.                 (5)
```

The two inversions in (5) are *only* for the analysis of positivity: they are not SOS feature functions or multipliers in (3). The matrix congruence showing equivalence between (4) and (5) is exact. At v=0 use the original homogeneous endpoint form rather than inverting its rank-deficient `A0` block. The unknown T(v) has four entries, each a polynomial of degree at most four with source-only coefficients. **Simultaneous positivity of both entire 2x2 matrices (5) on (0,infinity) remains unproved.** This is sufficient for H1, not necessary for membership in full K364.

## 3. A genuinely strict rank-six contact Gram, for EVERY source

Put `z=psi_*` and `x_+=(z,u_*z)` in the bilinear coordinate space R^8. The exact Alice sign congruence of (1) gives a second contact zero

```
z_-=(W_*,X_*,-Y_*,-Z_*),  x_-=(z_-,-u_* z_-).
```

The two projective contact points are distinct because `u_*>0` and both upper/lower state subvectors are nonzero. They are the **only** zeros of `F(u,v_*;psi)` on `RP^1 x RP^3`: SG6's strictly convex/coercive Laurent function has its unique minimum at `|u|=u_*` for fixed `v_*`, the three-state principal block has a one-dimensional nullspace, the Alice sign congruence gives the other sign, and the `u=0,infinity` projective faces are strictly positive by the original u^0 and u^2 endpoint squares (v_*>0). SG18-SG19 give a **positive-definite transverse Hessian** at both zeros.

Let `N=span_R{x_+,x_-}` (dimension two) and form the completely explicit **source-rational** orthogonal projector

```
Pi_R = I_8 - Xc*(Xc^T Xc)^(-1)*Xc^T, Xc=[x_+ x_-].        (6)
```

Its only inverse is a 2x2 positive definite *pure-source* matrix. The bilinear contact-annihilator square

```
Q_ref(s,t;psi) = (s psi,t psi)^T Pi_R (s psi,t psi) >=0      (7)
```

vanishes on the rank-one Segre variety precisely at the two physical contact points: a nontrivial linear combination of the two corresponding rank-one 2x4 matrices has rank two, since **both** Alice factors and both state factors are independent. In a compact normalization `s^2+t^2=1`, `||psi||=1`, `Q_ref` is O(dist^2) to these two zeros, while the physical nonnegative fixed-Bob form `D_(v_*) F` has a **positive quadratic lower bound** from the strict native Hessian. On the remaining compact set the form is strictly positive. Hence a sufficiently small **real positive** source-dependent epsilon satisfies globally

```
D_(v_*) F(u,v_*;psi) - epsilon*Q_ref(1,u;psi) >=0           (8)
```

for every homogeneous Alice pair and every state, including negative u and Alice infinity. This is a rigorous compactness/Hessian argument at a **single source and a single Bob angle**, not a compact *native-family* admission cover.

Calderon's 2-by-n biquadratic SOS theorem applies to (8), giving a PSD bilinear Gram `G_rem`. Every root of that SOS vanishes at both physical contacts, hence `G_rem x_+=G_rem x_-=0`. Thus

```
G_c = G_rem + epsilon Pi_R >=0,
ker G_c = N,  rank G_c=6 .                                  (9)
```

Average G_c under the exact Alice-sign involution on the eight bilinear coordinates. It still has kernel exactly N, retains the positive-definite complement, and becomes parity-block diagonal with blocks of order four and **one** zero vector in each. The polynomial coefficient identity of the contact Gram forces its diagonal 4x4 subblocks to equal `D_(v_*) P0(v_*)` and `D_(v_*) P2(v_*)`, and its cross symmetric part to equal `D_(v_*) P1(v_*)/2`. Therefore it is exactly a skew completion of the form (2) at v=v_*.

**Source-algebraic, not testing-dependent choice.** The feasible contact skew lies in a one-dimensional affine family given explicitly below. Its two three-dimensional contact complements have strictly positive definite Gram matrices. All their leading principal minors are polynomial inequalities of degree <=3 in that one free scalar, with source-rational coefficients. They therefore define a nonempty open semialgebraic subset of the real line. Select a point by isolating their finitely many algebraic roots and choosing a point in any feasible component. This is a finite piecewise source-algebraic selection (or a rational selection for each fixed real source), with no numerical SDP or testing-dependent Gram coefficient. The external Calderon theorem establishes nonemptiness; no statement about global Bob v has been used.

## 4. The EXACT affine contact-skew line

Write `w_*=(W_*,X_*)`, `y_*=(Y_*,Z_*)`, `h=u_*`, `D_*=1+v_*^2`. Let

```
U = D_* P0_w(v_*),    V = D_* P2_y(v_*),
L_* = D_* E(v_*)/2,
p = -U*w_*/h - L_* y_*,
q = -L_*^T w_* - h V y_* .                               (10)
```

At a PSD parity contact the two contact-kernel equations give exactly

```
T_* y_* = p,             T_*^T w_* = q,
w_*^T p = y_*^T q .                                      (11)
```

The last equality is guaranteed by the genuine state/u stationarity and the two sign-contact kernels, not assumed for arbitrary event weights. Define `s=w_*^T w_*>0`, `t=y_*^T y_*>0`, `w_perp=(-X_*,W_*)`, `y_perp=(-Z_*,Y_*)`. The **full solution** of (11) is the source-rational affine line

```
T_* = p*y_*^T/t + w_* q^T/s - (w_*^T p)*w_* y_*^T/(s t)
      + tau*w_perp*y_perp^T,             tau in R.         (12)
```

All inversions in (12) are positive source-only norms. The 3x3 principal minors of the two contact-orthogonal block compressions give an exact one-variable criterion for selecting tau with **both** contact blocks PSD and positive on their three-dimensional complements. This is not simply a kernel solution; existence of a strictly feasible tau has been proved by (8)-(9).

## 5. A TWO-SIDED polynomial PSD neighbourhood of the contact

The preceding rank-six Gram is stronger than the prior isolated PSD Calderon anchor. In fact it supports a **single source-algebraic polynomial** `T(v)=T_*+(v-v_*)T_1+(v-v_*)^2 T_2` whose complete 8x8 `H(v)` is PSD for *all Bob v in a two-sided interval* around v_*, and strictly PD when v!=v_* in that interval. The proof is an exact finite Schur argument, not an appeal to continuity of a singular matrix.

For its proof let `x=x_+`, `xbar=x_-`, `N=span{x,xbar}`, and `M=H(v_*)=G_c`. Define the source-only inverse-on-range

```
Msharp = (M+Pi_N)^(-1)-Pi_N,          Pi_N=I-Pi_R.         (13)
```

The matrix `M+Pi_N` is positive definite, so (13) has a strictly positive source-only determinant. This is the ordinary inverse of M on `N^perp`, with zero action on N; no singular pivot or mixed testing expression is inverted.

Define the physical Segre tangent space at x,

```
E = { (dotpsi,h dotpsi + dotu psi_*) : dotpsi in R^4,
                                           dotu in R } . (14)
```

It has dimension five and `E intersect N = R x`, because `xbar` is a different rank-one tensor whose two factors are independent of those of x. Consequently `L=E^perp intersect N^perp` has dimension two.

For any 2x2 source-only matrix U, write `S(U)=[[0,U],[-U^T,0]]` and the 8x8 symmetric coefficient-neutral correction `J(U)=[[0,S(U)],[-S(U),0]]`. Its effect on the contact vector is

```
J(U)x=(h S(U) psi_*,-S(U) psi_*),
(e^T J(U)x)=0 for every e in E.                      (15)
```

The map `U -> J(U)x` has image **all of E^perp** (dimension three): its four source entries act by `U*y_*` and `-U^T*w_*`, constrained only by `w_*^T U y_* = y_*^T U^T w_*`. Its kernel is exactly `R w_perp*y_perp^T`. Moreover, under the parity block split

```
x_A^T J(U)x_A= +2h w_*^T U y_*,
x_B^T J(U)x_B= -2h w_*^T U y_* .                     (16)
```

Native Bob stationarity makes the sum of the two first Bob derivatives zero. Pick the source-rational direction `U0=w_*y_*^T/(s t)`, for which `w_*^T U0 y_*=1`, and choose a scalar multiple in T_1 to make **both** first soft coefficients zero. The remaining two independent effective directions are

```
U1=w_perp*y_*^T,             U2=w_* y_perp^T,             (17)
```

with `w_*^T U1 y_*=w_*^T U2 y_*=0`. Their images `l_i=J(U_i)x` form a basis of `L`.

Let `H_1=dH/dv|_(v_*)`, whose source-only skew derivative is adjustable by U1,U2. With the first soft coefficients already zero, `b=H_1 x` belongs to `N^perp`. Form

```
G_ij = l_i^T Msharp l_j,      b_i = l_i^T Msharp b,
(alpha_1,alpha_2)^T = -G^(-1)*(b_1,b_2)^T.           (18)
```

Here `G` is **strictly positive definite**, because `l_1,l_2` are independent in `N^perp` and `Msharp` is positive there. Every denominator is a positive source-only determinant. Apply this correction to T_1. The resulting `b_opt` is the minimum-`Msharp`-norm vector in its exact two-dimensional affine gauge fiber. Its Euler equations say `Msharp*b_opt` is orthogonal to `L`, hence lies in `E+N`. Since M annihilates N, `b_opt=M e` for some physical tangent `e in E`.

Expand `H(v_*+h_v)=M+h_v H_1+h_v^2 H_2+O(h_v^3)`. The **complete** second-order minimized physical variation over `(dot u,dot psi)` at the strict native contact is

```
kappa_phys = x^T H_2 x - b_opt^T Msharp b_opt .      (19)
```

Indeed the Taylor quadratic along the actual rank-one Segre tangent is
`e^T M e+2 e^T b_opt+x^T H_2 x`; its minimum over E agrees with the full-range Schur minimum because `b_opt in M(E)`. The strict five-variable state/angle Hessian of SG18-SG19 makes that minimum **strictly positive** for a unit Bob displacement. Hence `kappa_phys>0` on the *entire genuine strict native C1 family*, with no special limit, source ratio, or sign assumption beyond native strictness.

The parity blocks are orthogonal. Their two soft second-order Schur coefficients `kappa_A,kappa_B` (computed with `Msharp` on each block) sum exactly to `kappa_phys>0`. A source-only change of the quadratic skew coefficient `T_2 -> T_2+delta U0` shifts these two coefficients by **opposite** `+2h delta,-2h delta`, by (16), without altering their sum or the linear soft jets. Choose

```
delta = (kappa_phys/2-kappa_A)/(2h).                (20)
```

Then both soft second-order Schur coefficients equal `kappa_phys/2>0`. Each original 4x4 contact block is positive definite on its three-dimensional complement. Schur elimination there gives, on the scalar soft line, precisely

```
h_v^2*(kappa_phys/2 + O(h_v)),                      (21)
```

and a complementary 3x3 block remaining PD. Thus **the entire H(v)** is PSD at v=v_* with kernel N and **positive definite** at every noncontact Bob angle in some source-dependent open interval `|v-v_*|<eps_theta` (v>0). All coefficients of T(v) in this proof are source-only rational/algebraic, selected by the finite formulas (12),(18),(20). The test features of (3) remain the same fixed 40 rows.

The conclusion (21) is a local positive-definite matrix polynomial, **not yet a globally sound Gram identity on the whole Bob half-line**. It is stronger than a single PSD anchor and closes the formerly unsupported neighbourhood inference, but it must not be used as H1.

## 6. The exact unresolved global implication

The same source-only polynomial T(v) constructed in §5 need **not** satisfy (5) for remote v; the derivation says nothing about those remaining angles. The previously verified endpoint source-rational choice `S_0,S_1,S_4` gives strict PSD near Bob zero and infinity, but has not been shown to agree with the three contact jets while retaining positivity elsewhere. The current precise sufficient problem is:

> For every genuine strict oriented C1 source theta, construct source-only polynomial T_theta(v) of degree <=4 (or a direct full 364-feature source-only Gram) satisfying BOTH matrix inequalities (5) for every positive real Bob v, with the endpoint homogeneous conditions. A proof must retain the complete full-state parent and cannot divide by a mixed source/testing polynomial as an SOS multiplier.

Proving that would give H1 by the already verified 40-feature half-line matrix-SOS implication and K364 embedding. It has **not** been proved here; failure of this more restrictive cross-parity 40-feature ansatz would **not** refute H1, which allows all 364 fixed features. Neither an exact K364 separator nor a source-wide PSD certificate has been found. The qualitative Q222 existence theorem, the evaluated conditional outer pencil, C2 results, and singular-density/reflection mechanisms are untouched.

## 7. Deterministic exact identity evidence

`verify_contact_geometry.py` independently constructs the unrestricted SG2 C1 parent for seven independent source weights, checks its full Alice coefficient extraction, all off-parity vanishing entries of the 8x8 matrix, the unrestricted scalar equality, the closed affine skew-contact line (12), the entire five-dimensional Segre-tangent annihilation in (15), the opposite soft-pivot transfers (16), and the exact rank-three/two-dimensional derivative gauge action. It runs symbolic polynomial arithmetic only; no source samples, SDP, parameter sweeps, or repeated native/rank baseline calculations.

The nondegenerate rank-six Gram and two-sided positivity theorem are **pen-and-paper consequences** of SG6, SG18-SG19, Calderon's biquadratic SOS theorem, and the stated exact finite-dimensional Schur construction. The checker is not presented as an automatic proof of Calderon, positivity of the native Hessian, or global admission. The original fixed-Bob full-parent, endpoint and interpolation source/checkers are retained unchanged inside the companion ZIP.
