# C1 C14.4: the doubly infinitesimal BN/UN source cell

**Scope.** This is a self-contained coefficient-assembly proof conditional only on the *native source identities* ES1–ES10 / SN3 / the explicit polynomial contact equalities in SN10, SN16 and UN12 (these are reproduced as monic coefficient identities below). It does **not** appeal to a general physical-C1 tensor-SOS closure result. It treats the actual C1 parent with both unrestricted state normals, and the two source cells `S,kappa,nu` positive infinitesimal, `ell=1-nu`, `q=nu/ell`, with `tau=h/q` bounded (BN) or unbounded (UN). Every sign and limit is meant in two separately labelled real closed fields `K_c,K_t`; the commutative tensor ring is `T=K_c tensor_R K_t` and its bounded subring has the real residue homomorphism. The finite exact checker accompanying this note proves the delicate polynomial helper and source-cancellation identities, not the asymptotic inequalities themselves. Sources: supplied exact c1_short.tex lines 13176–13953 and the preceding ES1–ES10 native entry, including the physical signs and gaps.

## 1. Actual parent, scalar, two normals and exact contact division

Normalize the positive CHSH coefficient to `d_B=1`. Use `D_u=1+u^2`, `D_w=1+w^2` and retain **all** original C1 event coefficients `m,C,A,D,alpha,beta` with `D,beta >=0`. With the reciprocal Bob testing coordinate `w`, the original unnormalized parent is

\[
\begin{aligned}
\mathcal F:=F_w/4={}&m(W-wX+uY-uwZ)^2+C D_w(uW-Y)^2+A D_uD_wY^2\\
 &+(D+\tfrac12)(uwW+uX-wY-Z)^2
 +(\alpha+\tfrac12)D_u(Y-wZ)^2\\
 &+(\beta+\tfrac12)D_w(X+uZ)^2-\tfrac12D_uD_wZ^2. \tag{1}
\end{aligned}
\]

Set `x=kappa S`, `delta=w^2-q`, `ell=1-nu`, `zeta=ell/kappa`, `e=kappa h/nu=kappa tau/ell`, and source-scaled weights `c0=SC(1-r), a0=SAr, ds=D/(S q^2), ns=beta/q^2, mu=m/kappa^2`. All `S,kappa,nu,q` are strictly positive infinitesimals; `ell` is a positive unit. Define the exact two independent normals and five scalar rows:

\[
\begin{gathered}
 X_e=X+ewW,\quad N_y=\frac{Y-ruW}{\sqrt{Sr(1-r)}},\qquad
 N=\kappa W-\ell uwZ,\quad N_0=N/\kappa=W/\zeta-uwZ,\\
 G=\frac{uW}{\sqrt S}-\ell\sqrt S\,wZ,\quad
 B=\sqrt S\,(w^2-q)Z,\quad B_u=u(w^2-q)Z,\\
 R=w(u^2-x)Z=\frac{\kappa\sqrt S G-uN}{\ell}
 =\frac{\sqrt S}{\zeta}G-uN_0.
\end{gathered}\tag{2}
\]

No state normalization or source/test scale comparison occurs. With `a=u^2,c=w^2,d=D+1/2`, the zero-normal scalar is **exactly**

\[
 \mathcal F|_{X_e=N_y=0}= A_NW^2+\frac{a}{S}A_P W^2-2KuwWZ+F_ZZ^2, \tag{3}
\]
where
\[
\begin{aligned}
 A_P&=(1+c)[c_0(1-r)+a_0r(1+a)]+S[d(1-r-e)^2c+(\alpha+\tfrac12)(1+a)r^2],\\
 A_N&=m(1+ra+ec)^2+(\beta+\tfrac12)e^2c(1+c),\\
 K&=m(1+ra+ec)+d(1-r-e)+(\alpha+\tfrac12)(1+a)r+(\beta+\tfrac12)e(1+c),\\
 F_Z&=D+\alpha c+\beta a+(m+\alpha+\beta+\tfrac12)ac.
\end{aligned}\tag{4}
\]

This follows by direct substitution `X=-ewW,Y=ruW` in (1), including the negative CHSH packet. Set `J1=K-ell A_P`, `L1=J1-(ell/kappa)A_N`, and
\[
 H=F_Z-\ell^2ScA_P+ac[(\ell/\kappa)^2 A_N-2(\ell/\kappa)J1].
\]
The full actual stationary contact gives `L1(x,q)=H(x,q)=H_a(x,q)=H_c(x,q)=0` and `H(a,0)=D+beta a`. Consequently *finite monic polynomial division* (in `a-x`, then `c-q`) gives unique coefficient polynomials
\[
 \begin{aligned}
 L1&=(a-x)L_a+(c-q)L_c,\\
 H&=S(c-q)^2[d_s+\kappa n_s+cH_{02}]
 +(a-x)(c-q)[-qn_s+cH_{11}]+c(a-x)^2H_{20}. \tag{5}
 \end{aligned}
\]
These equalities use no inverse of `(a-x)`, `(c-q)`, `u`, `w`, or a mixed tensor. The auxiliary quotients are literal differences of polynomials (for example `L_a=(L1(a,c)-L1(x,c))/(a-x)`). The canceled combinations, **not** their potentially singular constituents, are
\[
 C_R=H_{20}+2\ell L_a/\kappa,
 \quad C_G=\kappa(H_{11}-n_s)/\ell,
 \quad C_N=(n_s-H_{11})/\ell-2L_c/\kappa.
\tag{6}
\]
In particular the source identity `e=kappa tau/ell` yields `C_R=-ell mu r(a ell r+2c kappa tau+2ell)` and the SN12 bounded combination for `C_N`. The exact contact relation `uN=kappa sqrt(S)G-ell R` and the universal polynomial conductor
\[
 c(a-x)(c-q)Z^2=\{\kappa G(wB)-N(wB_u)\}/\ell \tag{7}
\]
transform (3) into the **complete scalar identity**
\[
\boxed{F_0=\frac{A_N}{\kappa^2}N^2+A_PG^2+d_sB^2+n_s B_u^2
 +H_{02}(wB)^2+C_RR^2-2\sqrt S L_a GR+C_GG(wB)+C_NN(wB_u).} \tag{8}
\]
One decisive coefficient check is `kappa B^2-q(a-x)(c-q)Z^2=B_u^2-c(a-x)(c-q)Z^2`. It produces the `n_sB_u^2` in (8); deleting that term would destroy the actual `beta` coefficient.

For clarity **both** normal blocks are obtained by coefficient comparison of (1), with
\[
\begin{aligned}
 H_X&=mc+da+(\beta+\tfrac12)(1+c),\\
 H_Y&=rc_0(1+c)+(1-r)a_0(1+a)(1+c)\\
 &\quad+Sr(1-r)[ma+dc+(\alpha+\tfrac12)(1+a)],\\
 H_{XY}&=-uw\sqrt{Sr(1-r)}(m+d).
\end{aligned}\tag{9}
\]
Writing
\[
\begin{aligned}
 A_X&=-m(1+ra+ec)+da(1-r-e)-e(\beta+\tfrac12)(1+c),\\
 T_X&=(\ell/\kappa)cA_X+\beta-D+c(m+\beta+\tfrac12),\\
 Y_P&=(1+c)[-c_0+a_0(1+a)]
 +S[m(1+ra+ec)-dc(1-r-e)+(\alpha+\tfrac12)(1+a)r],\\
 Y_C&=\ell Y_P-ma+d-(\alpha+\tfrac12)(1+a),
\end{aligned}
\]
exact contact vanishing implies
`T_X=c(a-x)X_a+(c-q)((D-beta)/q+cX_c)` and
`Y_C=(a-x)Y_a+(c-q)Y_c` by monic division. The **entire** first-normal terms are
\[
\boxed{\begin{aligned}
L_X={}&w\widehat A_X N+rac{uw\sqrt S}{\ell}\widehat X_aG
 +d_s\sqrt S\,q uB-n_sqB_u+wX_c(wB_u),\\
L_Y={}&\sqrt{r(1-r)}[Y_PG+\sqrt S Y_a R+Y_c(wB)],\\
\widehat A_X&=A_X/\kappa-aX_a/\ell,
\quad \widehat X_a=\kappa X_a=\ell[d(1-r-e)-mr].
\end{aligned}} \tag{10}
\]
Indeed `X=W` elimination first gives `w A_X W+u[mc-d+(beta+1/2)(1+c)]Z` and the exact identity for `T_X` above; the two relations `uN=kappa sqrt(S)G-ell R` and `d_s sqrt(S)q uB-n_sqB_u=(D-beta)B_u/q` give (10). Likewise direct differentiation with respect to `N_y` gives the second row. Thus
\[
 \boxed{\mathcal F=F_0+H_XX_e^2+H_YN_y^2+2H_{XY}X_eN_y+2L_XX_e+2L_YN_y} \tag{11}
\]
is an identity for **all** `W,X,Y,Z`, not only at the contact or at zero normals.

## 2. A precise finite real Gram allocation lemma

All coefficient vectors below are *named state-linear rows*, not presumed independent. A formal Gram on `g=(g_1,...,g_N)` is obtained mechanically: the entry for `g_i^2` is its coefficient; an `a g_i g_j` term contributes `a/2` in both symmetric off-diagonal positions. Replacing rows via the polynomial relations (2), (7) and those below is an exact equality of the complete parent. 

**Pivot lemma.** Suppose primary rows have a bounded real-residue positive definite coefficient block `D`, and helper rows `h_j` have bounded coefficients (including all cross terms). If exact identities let us add *any fixed real* positive diagonal `K_j h_j^2` at a zero-residue perturbation of the primary and already installed helper entries, while retaining separately a finite SOS, then choose `K_j` large enough that the finite real residue matrix has positive Schur complement `diag(K_j)+H-B^T D^{-1}B`. This produces one bounded Gram with strictly positive definite real residue. Its evaluation is tensor SOS by the strict-residue lemma; the retained packets are tensor SOS by their displayed separated roots. The `K_j` may be chosen **after** fixing the labelled source and test field elements; no common cap/feature library is needed for the qualitative criterion.

A helper linear identity `h=\sum_i\varepsilon_i g_i`, with `\rho(\varepsilon_i)=0`, buys a fixed `K h^2` by subtracting `K(\sum_i\varepsilon_i g_i)^2` from reserved primary squares. For a helper *quadratic budget*, use the exact zero identity explicitly. This is finite source/test algebra, never an inverse of a mixed tensor or a Gram chosen from a numerical sample.

## 3. BN: bounded `tau`, including `beta=0` or `D=0`

The **original physical source**, not merely stationarity, provides (BN2): with `b=beta/q`, `z=kappa beta/q^2`,
\[
 \tau,\mu,b,z,d_s\text{ bounded},\quad \tau,\mu\text{ positive units},
 \quad Jd_s+\theta z\sim\tfrac12,\quad\theta=1-h\sim1,\quad
 A_{00}:=\mu\ell^2-b\sim\tfrac12. \tag{12}
\]
Here `J` is a unit; hence **at least one** of `d_s,z` is a positive unit. The source identity
`b=[tau(B_native/S)-ell H0]/[2 theta(P_native/S)]` proves the positive-unit lower bound for `tau` from `beta>=0` and, with `mu ell^2=tau(1+x)(1+sigma r)/(2M(P_native/S))`, proves `A00~1/2` by cancellation. `c0,a0,alpha/S~1/2`, while `m,D,alpha,beta` are infinitesimal. The bounds do not divide by `beta/q^2` and hold through all source/test relative rates.

Define four **actual elementary source/testing monomials**
\[
 A_1=\sqrt\beta\,uZ,\quad A_2=\sqrt\beta\,W/(\zeta\sqrt q),\quad
 A_3=\sqrt\beta\sqrt{x/q}\,wZ,\quad
 A_4=\sqrt\beta\,uwW/(\zeta q\sqrt x).
\]
Their product satisfies `A1 A2 A3 A4=(beta uwWZ/(zeta q))^2`. Consequently
\[
 \Psi=\sum A_i^2-4\beta uwWZ/(\zeta q),\qquad
 \Psi-\tfrac12[(A_1-A_2)^2+(A_3-A_4)^2]\in\Sigma T^2. \tag{13}
\]
For `WZ>=0`, this follows from the two geometric means; for `WZ<0`, the signed product is a positive separated square. This includes literal `beta=0`. Reserve fixed fractions of `E_1=A1-A2`, `E_2=A3-A4` and retain the remaining circuit SOS in (13).

Subtracting exactly `Psi` from (8), define the SN-to-BN coefficient shift (BN8–BN12). Crucially the resulting scalar and the **two unchanged complete normals** become
\[
\begin{split}
 F_0={}&\Psi+A N_0^2+A'_P G^2+d_sB^2+H_{02}(wB)^2
 +C_RR^2-2\sqrt S L_aGR\\
 &+C'_GG(wB)+C'_NN_0(wB_u),\qquad
 L_X={}&\ell w\widehat A_X N_0+rac{uw\sqrt S}{\ell}\widehat X_aG
 +d_s\sqrt S\,q uB-bB_u+wX_c(wB_u),\qquad L_Y={}&\sqrt{r(1-r)}[Y_PG+\sqrt S Y_aR+Y_c(wB)].
\end{split}\tag{14}
\]
The exact canceled coefficients are
\[
 A=\mu\ell^2(1+ra+ec)^2+(\beta+\tfrac12)\tau^2 c(1+c)-b,\quad
 \rho(A)=\rho(A'_P)=\tfrac12,\]
\[
 C'_G=\kappa(H_{11}-2n_s)/\ell,\qquad
 C'_N=-Sa_0\ell^2r+(\beta+m+\tfrac12)\tau^2(c+1/\ell)
      +2\ell\kappa\mu\tau-m\tau^2. \tag{15}
\]
`C_R,L_a,C'_G,C'_N,H02,Y_a,Y_c,X_c,hat(A_X),hat(X_a)` are bounded by the explicit polynomial quotient formulas SN12–SN17 and (12); terms not routed to helpers have **zero residue** due to explicit factors `u,w,sqrt(S),kappa,q`, or the cancellation `Y_P\in ker(rho)`. Note `C'_G,C'_N` need only be bounded: their potentially nonzero real residues are paired with helpers, not claimed zero.

**Order of the two helper purchases:**

1. Put `g=(N0,G,X_e,N_y,E1,E2,sqrt(d_s)B,wB,wB_u)` (omit a literally zero packet). If `d_s` is a unit, `wB=(w/sqrt(d_s))(sqrt(d_s)B)` costs infinitesimally. Otherwise `z` is a unit and the **exact** identity
   \[
     wB=-\frac{w+\sqrt q}{\sqrt z}E_2-\frac{w(w+\sqrt q)}\ell G \tag{16}
   \]
   buys an arbitrary real pivot on `wB` at zero-residue Gram cost.
2. The exact homogeneous identity
   \[
    \delta^2N_0^2+\frac{2\kappa\delta}{\ell}G(wB)
          +2\kappa(wB)^2-(wB_u)^2=(\delta W/\zeta)^2 \tag{17}
   \]
   buys a pivot on `wB_u`: add `K(wB_u)^2`, subtract the first three terms multiplied by `K` from the current bounded Gram, and retain the positive square `K(\delta W/\zeta)^2`. All deducted coefficients have zero residue because `delta,kappa` do. No inverse of `delta` is ever taken.

The potentially dangerous `-bB_u` in the **unrestricted** `X_e` normal is eliminated exactly by
\[
 bB_u=-(w+\sqrt q)\sqrt b\,E_1-b(w+\sqrt q)N_0. \tag{18}
\]
All right-side coefficients have zero residue. The `d_s B` term is `sqrt(d_s)` times its positive packet; literal `D=0` simply removes that term. The `C_R R^2` and `L_a G R` terms contribute zero-residue entries by `R= (sqrt(S)/zeta)G-uN0`. The normal block diagonals satisfy `H_X~H_Y~1/2`, and `H_XY` has zero residue. The primary residue diagonals on `(N0,G,X_e,N_y)` are `(1/2,1/2,1/2,1/2)`; circuit rows and `sqrt(d_s)B` carry their retained positive pivots. Every remaining off-diagonal entry not incident to a purchased helper has zero residue. The two real pivots on `(wB,wB_u)` make the **complete** bounded Gram strictly positive by the pivot lemma. This certifies the exact parent (1), not just its zero-normal graph.

## 4. UN: unbounded `tau`, including `theta` infinitesimal and `D=0`

Here `theta=1-h>0`, `e=kappa tau/ell` is infinitesimal, and `e tau` is **bounded** by the original source sphere. The source quantities `m,D,alpha` are infinitesimal, `c0,a0,alpha/S~1/2`, and
\[
 n:=\beta\theta^2+2Dh+\tfrac12h^2>0,\quad
 z_n=\kappa n/q^2,\quad
 d_s,z_n\text{ bounded},\quad d_s+z_n\sim\tfrac12. \tag{19}
\]
At least one `d_s,z_n` is a unit. The actual `beta>=0` equation rules out literal `beta=0` here: it would imply `tau=ell H0/(B_native/S)`, bounded. Literal `D=0` is allowed and gives a unit `z_n`.

The **critical unbounded radial cancellation**, not a separate estimate of its terms, is
\[
\begin{gathered}
 A_{00}=m\zeta^2-n/q,\quad
 A_{00}=\frac{\theta\ell H_0}{2P_0}
   +\frac{e T_\kappa+h\ell S T_\nu}{2JMP_0},\\
 T_\nu=r-1-\nu-\nu S(1-r)(r+\ell),\quad
 T=\ell JU+JM_*P_0-JM_*\kappa\theta(J-\sigma)-2M_*U,\\
 U=(1+x)(1+\sigma r),\quad M_*=\ell M,\quad
 T=\kappa T_\kappa+\ell\nu S T_\nu,\quad T_\kappa\text{ bounded polynomial}.
\end{gathered}\tag{20}
\]
The exact native relations are `2theta P beta=h B_native-sigma H0`, `2P(JD+x theta beta)=sigma^2 mathcal H`; with `P0=P/S`, `sigma=nu S`, `x=kappa S`, these prove (20) by polynomial subtraction and monic division in `kappa`. All `J,M,P0` are positive units. Hence `A00-theta/2` is infinitesimal even though `m zeta^2` and `n/q` are independently unbounded. Also `n/q=m zeta^2-A00~tau/2` is unbounded, so `q/n` is infinitesimal. When `theta` is infinitesimal, `n~1/2` and `kappa/q` is infinitesimal. The exact `X_e`-coefficient pole likewise cancels:
\[
 t_0=\frac{D+h/2-\beta\theta}{q}
 =\frac Dq+\frac{\ell H_0-\kappa\tau(J-\sigma)}{2P_0}\sim\tfrac12. \tag{21}
\]
Both equations (20)–(21) are independently replayed in the attached finite symbolic checker.

Use the **alternative** unrestricted first normal `N_x=X+huZ` and the same `N_y`. Define `N0,G,B,Bu,R` by (2), and an additional *actual positive circuit*, with `n` replacing `beta` in the four `A_i`. Let `E1=A1-A2,E2=A3-A4,E3=A1-A3`, and `Psi_n= sum A_i^2-4nuwWZ/(zeta q)`. Averaging the two geometric-mean pairings (including the signed-state case) yields
\[
 \Psi_n-\tfrac14(E_1^2+E_2^2+E_3^2)\in\Sigma T^2. \tag{22}
\]
Reserve all three rows at positive fractions. For `WZ<0`, the remainder is a positive quadratic form in the four `A_i`: its matrix is `I-L_path/4`, where the four-vertex path Laplacian has largest eigenvalue `2+sqrt(2)<4`; the remaining signed cross term is a separated positive square. For `WZ>=0`, average the two exact geometric-mean pairings. After exact subtraction, the contact monic division (UN12) and its explicit combined coefficient identities (UN13) yield the **complete scalar**
\[
\begin{aligned}
 F_0={}&\Psi_n+A_{00}N_0^2+\widetilde A_PG^2
 +2u\sqrt S(L_a-C_R/\zeta)N_0G+d_sB^2+dG_h^2\\
 &+C_GG(wB)+C_NN_0(wB_u),\qquad
 \widetilde A_P\sim\tfrac12,\quad
 C_R=-m\zeta^2r(2+ra),\quad
 C_G=2\ell P_1(x)+\ell xa_0r-2z_n/\ell,\quad
 C_N=-Sa_0\ell^2r.
\end{aligned}\tag{23}
\]
The unbounded `C_R` has **already been combined** with the original `A`: `A+a C_R=A00`. No isolated `C_R R^2` is bounded. Instead `C_R/\zeta=-m\zeta r(2+ra)` is infinitesimal and `C_R/\zeta^2=-mr(2+ra)` is bounded, so the displayed coefficients in (23) are bounded, and their non-helper cross terms have zero residue.

The **entire unrestricted** `N_x` coefficient is the printed UN15 row
\[
 L_X=-m\zeta wN_0+d\tau a(wN_0+B_u)
 +uw\sqrt S[d(1-r-e)-mr]G+t_0B_u.
\]
The potentially unbounded `tau` coefficient is *not estimated*; use the exact positive hard-row identity, with `Q=wN0+Bu`:
\[
 \boxed{\begin{aligned}
 G_h&=-\tau u Q+ew\sqrt S\,G,\qquad
 N_x+\tau Q=X+ewW,\\
 dG_h^2+daN_x^2+2d\tau a Q N_x
 &=da(X+ewW)^2-2de\tau uw\sqrt S\,QG+de^2Sw^2G^2.
 \end{aligned}}\tag{24}
\]
Retain the whole nonnegative square `da(X+ewW)^2`, `beta(1+c)N_x^2`, and the last positive `de^2 S w^2G^2` if convenient. Because `e tau` is bounded, the only new signed `QG` coefficient is zero-residue. The remaining `N_x` diagonal is `mc+(1+c)/2~1/2`; the `N_y` diagonal is `~1/2` and its cross has zero residue. The full `N_y` coupling (UN17–UN18) is
\[
L_Y=\sqrt{r(1-r)}[Y_PG+\sqrt S Y_aR+Y_c(wB)],\quad
 Y_P\in\ker\rho,
\]
with bounded `Y_a,Y_c` across both `r` openings. Thus both first normals and the hard-normal block have been explicitly retained.

**Ordered helper allocation.** The finite row list is
\[
 g_{UN}=(G,N_x,N_y,N_0,E_1,E_2,E_3,\sqrt{d_s}B,B_u,wB), \tag{25}
\]
removing any literally zero packet. First, `wB` is free via `wB=(w/sqrt(d_s))sqrt(d_s)B` if `d_s` is a unit; otherwise `z_n` is a unit and
\[
 wB=-\frac{w+\sqrt q}{\sqrt{z_n}}E_2-\frac{w(w+\sqrt q)}\ell G. \tag{26}
\]
If `theta` is a unit, `A00~theta/2` already supplies `N0`. If `theta` is infinitesimal, `A00` has residue zero (possibly signed); obtain a **new** strict `N0` pivot *first* using
\[
 N_0=\sqrt{q/n}(E_2-E_1)-\frac{w-\sqrt q}{\sqrt n}E_3
            +\frac{\sqrt{\kappa/q}}\ell wG. \tag{27}
\]
Every coefficient is infinitesimal since `n~1/2`, `kappa/q` infinitesimal. Having installed `N0`, the exact identity
\[
 B_u=-(w+\sqrt q)\sqrt{q/n}E_1-(w+\sqrt q)N_0 \tag{28}
\]
buys an arbitrary real pivot on `B_u` at zero-residue cost. Since `w` is infinitesimal, `wB_u` is then free. This order is noncircular and uses **no** inverse of `theta`, `D`, or `beta`. All remaining signed terms of (23)–(24) and the entire `N_y` coupling have bounded helper coefficients or zero-residue primary entries. The pivot lemma makes the formal coefficient Gram on (25) strictly positive definite after finitely many real cap choices, including `theta` infinitesimal. The exact SOS packets (22), (24) and positive event excesses restore (1) completely.

## 5. Every independent testing angle, state, source opening and sign

The local BN/UN calculation covers both `u,w` infinitesimal, even when either is literally zero, without division by a testing coordinate or a state coordinate. Both signs of `WZ` are covered by the piecewise exact circuit decompositions. The full parent (1) is quadratic in **all four** unrestricted state amplitudes; neither `X_e,N_y` nor `N_x,N_y` is specialized to zero. Physical sign changes carry negative real test angles to the nonnegative chart. Complex states split into two real quadratic evaluations.

Since `S` infinitesimal and the exact original source gives `c0,a0~1/2`, the positive original events `C,A` are unbounded (`C=c0/[S(1-r)]`, `A=a0/[Sr]`), including the `r` endpoint openings. Thus the established adjacent-`C,A` theorem CE1–CE10 covers every test with Alice not infinitesimal. For infinitesimal Alice, the original unmatched Bob theorem SE10–SE12 covers all projective Bob sine-square residues different from `ell` (here `ell~1`), including Bob zero, infinity and all finite noncritical angles. The sole matched residue has reciprocal Bob `w` infinitesimal and is exactly the present BN/UN chart. Thus all independent projective tests and all states are covered. Source `tau>0` is either bounded (BN, where `beta=0` is allowed) or unbounded (UN, where `D=0` is allowed); no rate subcase is omitted. No formal negative-event source or physical-C1 closure is invoked.

**Theorem (conditional only on the earlier exact native polynomial identities and the already-assigned CE/SE outer tests).** For every actual physical doubly infinitesimal C14.4 source with `kappa,nu` positive infinitesimals, `ell=1-nu`, and any permitted `tau=h/q`, the complete original parent at every independent testing point and every state belongs to `Sigma(K_c tensor_R K_t)^2`. The strict bounded Gram is selected pointwise from (14) or (23), with helper pivots (16)–(18) or (26)–(28). There is **no claim of a fixed feature list** or an evaluated global SDP, and no need for uniform real caps. The broader finite-SDP conclusion depends on the independent global tensor theorem and other local sectors, not on this C1 cell alone.

## Reproduction and limitation

`Q222_C1_BN_UN_exact_checker.py` checks 16 symbolic identities with no random samples or optimizer, including (2)'s radial conductor, (17), BN (16),(18), UN (24),(26)–(28), native `P4`, and the exact source-pole cancellations (20)–(21). The unexpanded source estimates (BN2/UN3), the original contact polynomials SN10/SN16/UN12, and the previously assigned CE/SE outer theorems are **explicit proof dependencies**; this note supplies the missing **complete coefficient/Gram assembly** from them. The separate `Q222_C1_BN_UN_parent_expansion_checker.py` verifies the scalar, both quadratic normals, their mixed hard block and both complete first-normal coefficients **directly from the original seven-event parent** (six more symbolic equalities). Each is a fixed algebraic or ordered-field identity, not a uniformity hypothesis on independent tests.
