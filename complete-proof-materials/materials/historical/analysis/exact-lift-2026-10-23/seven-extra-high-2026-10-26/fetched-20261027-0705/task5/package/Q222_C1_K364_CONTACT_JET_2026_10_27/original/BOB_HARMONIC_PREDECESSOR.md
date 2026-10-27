# Oriented strict C1: an exact full-state Bob-harmonic construction, with a source-wide graded PSD boundary theorem

**27 October 2026.** Read-only analytical result. Source: the original C1 seven-coquotient dictionary SD1/NC1 and full four-state parent SG2/SD2, on the actual strict **oriented** domain `0<k<=1`. This note does **not** prove H1 (`C1 ⊂ K364`) and does not purport to refute it. It derives one concrete globally fixed-grid **sufficient** PSD mechanism, proves its coefficient identity, proves its genuinely simultaneous whole-state PSD near the Bob-zero and Bob-infinity ends for **every** strict native source, and shows how to preserve a source-algebraic PSD contact anchor. The remaining middle-angle matrix positivity is explicitly unsolved.

## Theorem and genuine native hypotheses

Write the original canonical C1 source as

\[
 c=(d,0,0,m,C,A,D,\alpha,\beta),\qquad d,m,C,A,D,\alpha,\beta>0,
\]

with all seven entries given **exactly** by SD1 and the seven original native coquotients in `source_refs/c1_native.py`, normalized by `C=2`. The native parameters \(k,f,g,a,t>0\) obey \(0<k\le1\),

\[
\begin{split}
 L_c&=1-af^2>0,&H_c&=1+(1-a)f^2>0,\\
 Q_c&=1-k^2g^2+ak^2f^2g^2>0,\\
 M&=1+k^2(1-af^2)-t>0,&
 J&=1+k^2-k^2g^2t-a>0,
\end{split}
\]

and the remaining strict native factors and the original \(D,\beta>0\) conditions in SD1/NC1. Their positive nonzero stationary contact is
\[
 (u_c,v_c;\psi_c)=\left(k,kfg;\frac gf,-k(1-tg^2),\frac{kg(1-af^2)}f,1\right).
\]
No coefficient is free beyond these true source formulas. Here \(\psi=(W,X,Y,Z)^T\), \(D_u=1+u^2\), \(D_v=1+v^2\), and
\[
\begin{aligned}
F_c={}&4m(vW-X+uvY-uZ)^2+4CD_v(uW-Y)^2+4AD_uD_vY^2\\
&+(4D+2d)(uW+uvX-Y-vZ)^2+(4\alpha+2d)D_u(vY-Z)^2\\
&+(4\beta+2d)D_v(X+uZ)^2-2dD_uD_vZ^2.
\end{aligned}\tag{1}
\]
All identities and matrices below use the **whole polynomial (1)**, before any state or testing restriction.

### Statement 1: four globally valid positive endpoint forms

Put \(K=4C+4D+2d>0\), \(L=4\beta+2d+4m>0\), and decompose
\[
 F_c(u,v;\psi)=\psi^T\big(P_0(v)+uP_1(v)+u^2P_2(v)\big)\psi.
\tag{2}
\]
Then the *complete* endpoint forms are exactly
\[
\begin{aligned}
\psi^TP_0(v)\psi={}&4m(vW-X)^2+4(C+A)D_vY^2+4D(Y+vZ)^2\\
&+4\alpha(vY-Z)^2+4\beta D_vX^2+2dD_v(X^2+Y^2),\\[1mm]
\psi^TP_2(v)\psi={}&4m(vY-Z)^2+4CD_vW^2+4AD_vY^2+4D(W+vX)^2\\
&+4\alpha(vY-Z)^2+4\beta D_vZ^2+2d[(W+vX)^2+(vY-Z)^2],\\[1mm]
F_c(u,0;\psi)={}&K(uW-Y)^2+L(X+uZ)^2+4AD_uY^2+4\alpha D_uZ^2,\\[1mm]
[v^2]F_c={}&4m(W+uY)^2+4C(uW-Y)^2+4(A+\alpha)D_uY^2\\
&+4D(uX-Z)^2+4\beta(X+uZ)^2+2dD_u(X^2+Y^2).
\end{aligned}\tag{3}
\]
Here \([v^2]\) denotes the leading Bob-angle coefficient, not a value at one state. Direct termwise expansion of (1) proves all four polynomial identities. Every displayed coefficient is positive on the entire genuine strict native domain. Consequently \(P_0(v),P_2(v)\succ0\) for every \(v>0\). At \(v=0\), \(\ker P_0(0)=\mathbb Re_W\), \(\ker P_2(0)=\mathbb Re_X\); these are precisely the compatible-measurement corner degeneracies.

### Statement 2: a genuinely fixed 40-feature *sufficient* route into K364

Let
\[
 \mathcal V(v)=\sum_{j=0}^{4}v^jS_j,\qquad S_j^T=-S_j\in\operatorname{Mat}_{4\times4}(\mathbb R),
\]
with **source-only** entries, and form the symmetric \(8\times8\) polynomial matrix
\[
 \boxed{\mathscr H_c(v)=
 \begin{pmatrix}
 D_vP_0(v)&\tfrac12D_vP_1(v)+\mathcal V(v)\\
 \tfrac12D_vP_1(v)-\mathcal V(v)&D_vP_2(v)
 \end{pmatrix}.}\tag{4}
\]
For the *unrestricted state* and the independently variable Alice angle define \(x=(\psi,u\psi)^T\). Skew symmetry gives the exact coefficientwise identity
\[
 \boxed{x^T\mathscr H_c(v)x=(1+v^2)F_c(u,v;\psi)}\tag{5}
\]
for **every** choice of the skew coefficients. The entries of \(\mathscr H\) have degree at most four in \(v\).

If a source-only choice of these five skew matrices satisfies the genuinely global matrix inequality
\[
 \mathscr H_c(v)\succeq0\qquad\text{for every real }v\ge0,\tag{6}
\]
then it supplies an exact PSD Gram for \((1+v^2)F_c\) in the **fixed 40-feature library**
\[
 \begin{aligned}
 \mathcal Z_{40}={}&\{u^i v^j\psi_\ell: i=0,1,\ j=0,1,2,\ \ell=1,\ldots,4\}\\
 &\cup\{u^i v^{j+1/2}\psi_\ell: i=0,1,\ j=0,1,\ \ell=1,\ldots,4\}.
 \end{aligned}\tag{7}
\]
All forty rows are source independent and belong to \(Z_{(6,12)}\) because their sixth-root exponent indices obey \(0\le i\le6,0\le j\le12\).

**Proof of the implication (6) → fixed PSD Gram.** Apply the univariate real polynomial-matrix spectral factorization theorem to \(\mathscr H_c(t^2)\), which is PSD for every real \(t\) and has degree at most eight. It is a sum of polynomial matrix squares \(R(t)^TR(t)\), with \(\deg R\le4\). Split \(R(t)=R_{\rm even}(t^2)+tR_{\rm odd}(t^2)\), average the two square identities at \(t\) and \(-t\), and obtain
\[
 \mathscr H_c(v)=R_{\rm even}(v)^TR_{\rm even}(v)
       +vR_{\rm odd}(v)^TR_{\rm odd}(v),\quad
 \deg R_{\rm even}\le2,\quad\deg R_{\rm odd}\le1.
\tag{8}
\]
Substitution of \(x=(\psi,u\psi)\) gives squares in exactly the span (7), including singular polynomial matrices. Their constant coefficient vectors depend only on the source. Spectral factorization over the real algebraic closure of the native parameter field, or exact real quantifier elimination, gives real-algebraic source coefficients, not testing-dependent Gram coefficients. Embed their Gram by zero rows into the original 364-feature library. No mixed determinant, source/test denominator or nonaffine multiplier is canceled. This is a **sufficient condition** for H1, not an equivalent characterization of all K364 Grams.

### Statement 3: exact source-rational Bob-zero and Bob-infinity PSD completion

Take the exact endpoint Gram matrices of (3), in the ordered variables \(x=(W_0,X_0,Y_0,Z_0;W_1,X_1,Y_1,Z_1)\):
\[
\begin{aligned}
 x^TG_0x={}&K(W_1-Y_0)^2+L(X_0+Z_1)^2
          +4A(Y_0^2+Y_1^2)+4\alpha(Z_0^2+Z_1^2),\\
 x^TG_\infty x={}&4m(W_0+Y_1)^2+4C(W_1-Y_0)^2
 +4(A+\alpha)(Y_0^2+Y_1^2)\\
 &+4D(X_1-Z_0)^2+4\beta(X_0+Z_1)^2
 +2d(X_0^2+X_1^2+Y_0^2+Y_1^2).
\end{aligned}\tag{9}
\]
Thus \(G_0\succeq0\) with **exact** kernel \(\operatorname{span}\{e_{W_0},e_{X_1}\}\), and \(G_\infty\succ0\) because all eight root components are forced to zero by the positive rows. Let \(B_0=G_0[0{:}4,4{:}8]\), \(B_\infty=G_\infty[0{:}4,4{:}8]\) and set
\[
 S_0=\frac{B_0-B_0^T}{2},\qquad S_4=\frac{B_\infty-B_\infty^T}{2}.
\tag{10}
\]
These are explicit source-rational skew matrices; the endpoint identities (3) ensure the symmetric off-diagonal coefficients are exactly correct. Hence \(\mathscr H_c(0)=G_0\), and \(\lim_{v\to\infty}v^{-4}\mathscr H_c(v)=G_\infty\) regardless of \(S_1,S_2,S_3\). In particular **every** allowed choice is strictly positive for all sufficiently large \(v\) (with a source-dependent cutoff).

The novel boundary-safe choice is
\[
 \boxed{S_1=\tau(e_We_Z^T-e_Ze_W^T+e_Xe_Y^T-e_Ye_X^T),\quad
 \tau=2D+d-2m.}\tag{11}
\]
Choose any source-only skew \(S_2,S_3\) with \((S_2)_{WX}=0\). With \(N=(e_{W_0},e_{X_1})\), the entire \(N\)-to-positive-range first derivative of (4) at \(v=0\) consists only of
\[
 W_0\to-4m(X_0+Z_1),\qquad
 X_1\to(4D+2d)(W_1-Y_0).
\tag{12}
\]
The order-\(v\) \(N\)-to-\(N\) block is identically zero. The order-\(v^2\) direct kernel diagonal is \(\operatorname{diag}(4m,4D+2d)\), with zero mixed term because \((S_2)_{WX}=0\). Eliminating the **strict six-dimensional** range block of \(G_0\) gives the exact full two-soft-mode graded Schur matrix
\[
 \boxed{\Sigma_0=
 \begin{pmatrix}
 \dfrac{4m(2\beta+d)}{2\beta+d+2m}&0\\[2mm]
 0&\dfrac{4C(2D+d)}{2C+2D+d}
 \end{pmatrix}\succ0.}\tag{13}
\]
This is a matrix calculation, not merely two individual diagonal positivity tests. Both off-diagonal terms are proved identically zero. The exact all-source symbolic checker reconstructs the entire \(8\times8\) parent, the two endpoint Grams and both entries of (13), retaining all state normals.

**Rigorous local PSD consequence.** By the Schur formula, after the range block is eliminated, the full kernel Schur is \(v^2(\Sigma_0+O(v))\). The range block remains positive definite as \(v\downarrow0\). Therefore for **every strict real native C1 source** and any \(S_2,S_3\) as above, there is \(\epsilon_c>0\) such that
\[
 \mathscr H_c(v)\succ0\quad(0<v<\epsilon_c),\qquad
 \mathscr H_c(v)\succ0\quad(v>V_c)
\tag{14}
\]
for some finite \(V_c\). This is simultaneous whole-\(u\), whole-four-state PSD for the **same** polynomial Gram, with no testing/source rate comparison. It does **not** assert positivity on the intermediate Bob-angle interval.

### Statement 4: exact contact anchor compatible with both boundary mechanisms

Put \(z=\psi_c\), \(u_*=u_c>0\), \(v_*=v_c>0\), \(D_*=1+v_*^2\). The native state and angle stationarity at the noncommuting contact give
\[
 (P_0+u_*P_1+u_*^2P_2)(v_*)z=0,
\quad z^T(P_1+2u_*P_2)(v_*)z=0,
\quad z^T(P_0+u_*P_1+u_*^2P_2)(v_*)z=0.
\]
Define the **explicit source-rational vector**
\[
 b=-D_*\bigl(P_0(v_*)z/u_*+P_1(v_*)z/2\bigr),
 \qquad z^Tb=-D_*\bigl(F(u_*,v_*;z)/u_*-\partial_uF(u_*,v_*;z)/2\bigr)=0.
\tag{15}
\]
Then
\[
 S_c^{\rm ker}=\frac{bz^T-zb^T}{z^Tz}\in\mathfrak{so}(4)
\tag{16}
\]
is source-rational, and the **entire** \(8\times8\) matrix (4) at \(v_*\) with \(\mathcal V(v_*)=S_c^{\rm ker}\) annihilates \((z,u_*z)\). This is an exact kernel identity; it does *not* claim that this particular rational matrix is PSD.

There is, however, at least one **source-algebraic** \(S_c\in\mathfrak{so}(4)\) making the full \(8\times8\) contact matrix PSD. Indeed, homogenize \(F_c(u,v_*;z)\) in \((1,u)\). Physical global nonnegativity, including negative \(u\) via the exact Alice sign congruence and \(u=\infty\) via (3), makes it a nonnegative **\(2\times4\) biquadratic form**. Calderón's classical theorem (A. P. Calderón, *A note on biquadratic forms*, LAA 7 (1973), 175–177, DOI `10.1016/0024-3795(73)90052-9`) states that every nonnegative \(2\times n\) biquadratic form is SOS of bilinear forms. This yields a PSD \(8\times8\) Gram with the prescribed \(P_0,P_1,P_2\) blocks and some skew cross correction. Its entire contact vector lies in its kernel because the Bell polynomial vanishes there. The feasible skew set is defined by finitely many polynomial inequalities over the field generated by the native parameters; real-closed-field transfer permits an algebraic point or a finitely piecewise semialgebraic choice. **This is a one-Bob-slice theorem, not a PSD selection for all Bob angles.**

For any such chosen PSD contact skew \(S_c\), put
\[
 T=S_c-S_0-v_*S_1-v_*^4S_4,
\quad E_{WX}=T_{WX}(e_We_X^T-e_Xe_W^T),
\]
\[
 \boxed{S_2=(T-E_{WX})/v_*^2,\qquad S_3=E_{WX}/v_*^3.}\tag{17}
\]
These are source-only algebraic coefficients, with \((S_2)_{WX}=0\) and the exact interpolation
\(S_0+v_*S_1+v_*^2S_2+v_*^3S_3+v_*^4S_4=S_c\).
Therefore **one and the same degree-four skew polynomial** has: positive-definite \(\mathscr H\) at all sufficiently small positive \(v\), positive-semidefinite \(\mathscr H(v_*)\), and positive-definite \(\mathscr H\) for sufficiently large \(v\), while retaining coefficient identity (5) for every independent Alice/state coordinate. These are proven PSD regions, not an assertion about the untested interior.

## Exact remaining global positivity implication and H1 status

The still-unproved sufficient matrix statement for this explicit 40-row subcone is
\[
\boxed{\forall\theta\in\Theta_{C1}^{\rm strict,ori},\ \exists\ S_2,S_3\in\mathfrak{so}(4),\ (S_2)_{WX}=0:\quad
 \mathscr H_\theta(v)\succeq0\ \forall v\ge0.}\tag{18}
\]
The polynomial is defined explicitly by (2), (4), (10), (11). For each fixed positive Bob \(v\), \(P_0(v)\succ0\) by (3), so its **exact** positivity condition is the \(4\times4\) Schur inequality
\[
\boxed{D_vP_2(v)
 -\bigl(\tfrac12D_vP_1(v)-\mathcal V(v)\bigr)
 [D_vP_0(v)]^{-1}
 \bigl(\tfrac12D_vP_1(v)+\mathcal V(v)\bigr)\succeq0.}\tag{19}
\]
This inverse is used **only to analyze matrix positivity**, not as a multiplier of the Bell parent, a testing feature, or an affine source constraint. No positivity of (19) has been established for every interior \(v\); anchors and the local graded Schur (13) cannot replace it. This is the first precise unproved implication in this **specific candidate**, which is *stronger* than the original H1. A failure of (18) would not be a counterexample to \(K_{364}\), and neither (18) nor a numerical source example establishes H1. The full H1 remains the original fixed-grid coefficientwise PSD Gram equality involving 364 features.

No source-side reflection transport, singular-density closure, C2 construction, full 128-orbit, or global outer SDP is re-audited here. Their previously verified scopes are retained. No Git write was made.

## Exact verification

`verify_full_parent_bob.py` checks (1)–(5) over seven algebraically independent source weights, including all four unrestricted state coefficients and four endpoint SOS formulas, and verifies the 40 source-independent feature indices.

`verify_boundary_schur.py` checks (9)–(13) by an actual six-dimensional block elimination with five completely free higher skew entries \(S_2\), all six entries of \(S_3\), and all seven source weights independent; it confirms the full two-by-two matrix, not just two diagonal entries.

`verify_anchor_interpolation.py` checks (17), including zero \(WX\) second-order skew, contact matching, endpoint preservation, and full matrix symmetry with no testing inverse. The external Calderón theorem and the strictly positive Schur implications are pen-and-paper ingredients, not silently attributed to checker output.
