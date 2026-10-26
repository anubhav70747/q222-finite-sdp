# Oriented strict native C1: exact fixed Z18 charge transfer and the remaining PSD implication

**26 October 2026.** This is an analytical continuation of the verified C1 16-feature coefficient-map and obstruction package, not an assertion of universal strict-C1 Gram admission. The orientation is **exactly** `0 < k <= 1`. No numerical SDPs, parameter sweeps, guessed Gram positivity, connection-to-closedness argument, or GitHub write is used. The earlier strict `k=4` counterexample is outside the oriented source domain, though it remains a genuine obstruction to the old 16-feature cone on unrestricted A1.

## Theorem 1: genuine native scope and the full parent

Let `k,f,g,a,t` be real with `0<k<=1` and `f,g,a,t>0`, and put

\[
\kappa=k^2,\quad \eta=f^2,\quad \zeta=g^2,\quad
 V_c=\kappa\eta\zeta,\quad e=kf/g,\quad
 L_c=1-a\eta,\quad H_c=1+(1-a)\eta,\quad
 Q_c=1-\kappa\zeta+aV_c,
\]
\[
 R=1+\kappa L_c,\quad T=1+\kappa-\kappa\zeta t,\quad
 M=R-t,\quad J=T-a.
\]

The **strict original-source domain** requires

\[
 L_c,H_c,Q_c,M,J,1-t\zeta,1-\kappa\eta+tV_c>0,
\]

and strictly positive seven native coquotients `c2, cb, c3, c5, c4, calpha, cd` and their defining `P4,Palpha,Pd`, as in the unchanged `original_baseline/c1_native.py` and SD1 (these definitions are part of the exact source, not substitute free weights). With

\[
 d=\frac{2a\eta c_d}{c_3},\ m=\frac{2ae^2c_2}{Mc_3},\ C=2,\ A=\frac{2a\eta^2c_4}{c_3},\
 \alpha=\frac{2ae^{-2}\kappa^3\eta^2c_\alpha}{c_3},
\]
\[
 D=\frac{2ac_5}{Jc_3}-\frac d2>0,\qquad
 \beta=\frac{2ae^2c_b}{tc_3}-\frac d2>0,
\]

this is the genuine strict C1 source `(d,lambda,n,m,C,A,D,alpha,beta)=(d,0,0,m,2,A,D,alpha,beta)`. All denominators are positive on this domain, and the established native entry/stationarity theorem gives the nonzero physical contact

\[
 (u_c,v_c;W_c,X_c,Y_c,Z_c)
 =(k,kfg;g/f,-k(1-tg^2),kgL_c/f,1).
\tag{1}
\]

Use the original unrestricted four-state SG2 parent

\[
\begin{split}
F={}&4m(vW-X+uvY-uZ)^2+4C(1+v^2)(uW-Y)^2
 +4A(1+u^2)(1+v^2)Y^2\\
&+(4D+2d)(uW+uvX-Y-vZ)^2
 +(4\alpha+2d)(1+u^2)(vY-Z)^2\\
&+(4\beta+2d)(1+v^2)(X+uZ)^2
 -2d(1+u^2)(1+v^2)Z^2.
\end{split}
\tag{2}
\]

All coefficient statements below concern (2) **for all independent positive `u,v` and all four unrestricted state coordinates**. The original source identities and full-parent stationarity are in `original_baseline`, and are not silently replaced by just scalar contact or determinant positivity.

## Theorem 2: an exact universal 18-feature coefficient admission and precisely two fractional charges

The fixed, source-independent testing library is

\[
\begin{aligned}
Z_{18}^{\rm int}&=(uW,vW,X,uvX,Y,uvY,uZ,vZ,uvW,vX,uY,vY,Z,uvZ),\\
Z_{18}^{\rm half}&=(\sqrt{uv}W,\sqrt{uv}Z,\sqrt{uv}X,\sqrt{uv}Y).
\end{aligned}\tag{3}
\]

Every feature belongs to the fixed sixth-root grid `Z6` (the half-power exponent is `3/6`). The feature library does **not** depend on any native source. Source-only changes `u=u_c U`, `v=v_c V`, `psi=diag(psi_c)phi`, and division of individual rows by their nonzero contact evaluations are invertible diagonal congruences. No testing denominators or new testing degrees are introduced.

Use the twelve established integer contact differences

\[
\begin{split}
r_0,\ldots,r_6={}&(UW-VZ,\ VW-VZ,\ X-VZ,\ UVX-VZ,\\
&Y-VZ,\ UVY-VZ,\ UZ-VZ),\\
r_7,\ldots,r_{11}={}&(UVW-UVZ,\ VX-UVZ,\ UY-UVZ,\\
&VY-UVZ,\ Z-UVZ).
\end{split}\tag{4}
\]

Here the state letters are the **contact-normalized** `phi`, not the original unnormalized state. Introduce two fractional contact-root squares

\[
 h_W=\sqrt{UV}(W-Z),\qquad h_X=\sqrt{UV}(X-Y).
\tag{5}
\]

The underlying constant exact coefficient map on all contact-difference Gram entries has **46 generated angular/state coefficient keys, 49 entries, rank 40, cokernel 6**. **These 46 generated keys must not be confused with the full 46 parity-compatible bidegree-(2,2) keys**: they consist of the old 42 parity-compatible invariant-support keys plus four extra *parity-incompatible* keys (forced to zero). The four parity-compatible keys still absent from both selected libraries are the constant `W^2,WZ` and the `U^2X^2,U^2XY` coefficients; the genuine native parent has zero coefficients there. The old integer plus WZ-odd submap has rank 36 on its 42 invariant keys, with the previous eight gauges. A four-by-four new coefficient minor on the four new parity-incompatible monomials

\[
 (UV\,WX,UV\,WY,UV\,XZ,UV\,YZ)
\]

has determinant **16** (its columns are the half-root pairs `(0,1),(0,2),(1,1),(2,2)`). Thus the rank is at least `36+4=40`. Every contact-difference product satisfies the four state and two angle contact jets; these six independent rows bound the rank by `46-6=40`. The checker constructs the exact whole matrix and verifies these statements by integer arithmetic. This keeps the dimensional distinction precise: the **entire parity-compatible** polynomial space has 46 keys, but the selected invariant Bell support has 42; the four additional keys of the *generated* extended map are forbidden by the parent symmetry, not those four missing compatible keys.

For the actual C1 parent those four additional parity-incompatible coefficients vanish. Writing the general real symmetric Gram on the three half-root differences

\[
 h_0=\sqrt{UV}(W-Z),\quad h_1=\sqrt{UV}(X-Z),\quad h_2=\sqrt{UV}(Y-Z)
\]

and matching the four forbidden monomials gives the **complete exact odd-block form**

\[
 \Omega=\begin{pmatrix}\gamma&0&0\\0&\delta&-\delta\\0&-\delta&\delta\end{pmatrix},
 \qquad \Omega\succeq0\iff\gamma\ge0,\ \delta\ge0.
\tag{6}
\]

The integer/fractional cross block can be removed by averaging the formal square-root sign automorphism. The integer parity groups and the half WZ/XY sectors separate by averaging the exact simultaneous sign congruence `(U,V,X,Y) -> (-U,-V,-X,-Y)` (and keeping the formal square root of `UV` fixed). Each operation preserves the parent and PSD. At a genuine contact, each block has its corresponding all-ones contact kernel. Thus the contact-reduced description loses no PSD certificate from this 18-feature cone.

**The two new globally valid zero identities** in the contact-normalized polynomial ring are

\[
 \begin{aligned}
 r_1(r_0-r_6)&=UV(W-Z)^2=h_W^2,\\
 (r_2-r_4)(r_3-r_5)&=UV(X-Y)^2=h_X^2.
 \end{aligned}\tag{7}
\]

Consequently every actual strict C1 source has the **exact full-parent coefficient identity**

\[
\boxed{
 \widehat F_\theta
 =r_{0:6}^TH_7(\theta;s,\gamma,\delta)r_{0:6}
 +r_{7:11}^TH_5(\theta;s)r_{7:11}
 +\gamma h_W^2+\delta h_X^2,
}\tag{8}
\]

where `s=(s0,...,s6)` is source-only, and all entries are source-rational except for any algebraically chosen gauges. Take the original rational 36-pivot inverse and seven kernels from the recovered 16-feature proof, retain its eighth gauge `gamma`, and add **one** new integer correction `delta K_X`, where the only nonzero symmetric entries (zero-based indices in `H7`) are

\[
 (K_X)_{23}=-\tfrac12,\quad(K_X)_{25}=+\tfrac12,
 \quad(K_X)_{34}=+\tfrac12,\quad(K_X)_{45}=-\tfrac12.
\tag{9}
\]

Indeed `r^T K_X r=-h_X^2` exactly. Hence the 18-feature **C1 invariant** coefficient map is `42 x 45`, rank 36, with nine independent source-only gauge coordinates. Matching includes the complete first and second state-normal polynomials: (8) is an identity of all coefficients of (2), not just the six contact jets.

The cone of sources admitted by the original fixed `Z18` Gram is globally sound: its defining polynomial equality is linear in the Bell-source coefficients and one finite PSD Gram; positivity on positive testing angles follows from the square representation, and the original sign gauges and homogeneous limits cover all other projective tests. The fixed feature Gram cone is closed by the positive-definite moment-integral trace bound on `(u,v) in [1/2,2]^2` and the state unit sphere. **No assertion that this cone contains all strict C1 sources is made.**

### Why this is a real construction change

In the old 16-feature cone the Fourier/contact inequality imposes a lower bound on the *single* odd charge `gamma`, while a forced H7 principal Schur block bounds that same scalar above. The rigorous unrestricted-`k` obstruction exploits their conflict. In (8) the exact sign-related-contact argument (valid only on the physical Bell slice) becomes a bound on **`gamma+delta`**, while the old forced **upper bound on `gamma` alone** is unchanged by (9). The added source-only positive charge can therefore carry the previously conflicting demand. Explicitly, on the **canonical Bell slice with the four sign-related contact zeros**, put

\[
 S_X=(4\beta+2d)v_c^2X_c^2,\quad S_Y=4Au_c^2Y_c^2,\quad
 S_Z=4\beta u_c^2v_c^2Z_c^2,
\]

and define the three-side deficit

\[
 L_{\triangle}=\max\{0,\sqrt{S_X}-\sqrt{S_Y}-\sqrt{S_Z},
 \sqrt{S_Y}-\sqrt{S_X}-\sqrt{S_Z},
 \sqrt{S_Z}-\sqrt{S_X}-\sqrt{S_Y}\}.
\]

The same four-contact integer-root Fourier argument as in the recovered exact counterexample gives **the necessary inequality**

\[
\boxed{\gamma+\delta\ge L_{\triangle}^2.}\tag{9a}
\]

The old forced Schur upper bound acts on `gamma` alone. Identity (7) transfers a positive cross charge to `delta` while leaving every original full-parent coefficient matched, so the previous `gamma`-only contradiction no longer follows. This removes the specific *logical contradiction* but does **not** establish PSD feasibility of the whole H7/H5 pair, even for the old `k=4` source. The inequality (9a) is used only on the actual four-contact Bell-validity slice, **not** on arbitrary single-contact polynomials. It also does not prove the old 16-feature cone fails on `k<=1`.

## Theorem 3: exact complete Schur admission interface, including singular certificates

Let `H7,H5` be the source-rational affine matrices in (8), constructed without any inverse of testing coordinates. Define the **forced positive** block, with indices local to `H7`,

\[
 P=H_7[\{1,2\}]=\begin{pmatrix}
 4m v_c^2W_c^2&-4mv_cW_cX_c\\
 -4mv_cW_cX_c&2(2\beta+d+2m)X_c^2
 \end{pmatrix},
\]
\[
 \det P=8m(2\beta+d)v_c^2W_c^2X_c^2>0,
\quad p_4=H_7[4,4]=2Y_c^2(2A+2C+2D+d)>0,
\quad p_{11}=H_5[4,4]=4\alpha Z_c^2>0.
\tag{10}
\]

No gauge, including the new `delta`, changes any entry of `P` or either scalar pivot. Put `Tidx=(0,3,4,5,6)` and `Ridx=(0,3,5,6)` (global `H7` indices). Define

\[
 T_7=H_7[Tidx,Tidx]-H_7[Tidx,\{1,2\}]P^{-1}H_7[\{1,2\},Tidx],
\]
\[
\boxed{ S_7=T_7[Ridx,Ridx]
   -T_7[Ridx,\{4\}]\,T_7[\{4\},Ridx]/p_4, }
\tag{11}
\]

where in the second formula `Ridx` means the corresponding positions `(0,1,3,4)` **inside** `Tidx`. For the global `H5` block indexed `(7,8,9,10,11)`, put

\[
\boxed{ S_5=H_5[0{:}4,0{:}4]
       -H_5[0{:}4,4]H_5[4,0{:}4]/p_{11}. }
\tag{12}
\]

Exact block-triangular congruences prove, **without assuming any nonzero eigenvalue of `S7` or `S5`**,

\[
\boxed{
 H_7\succeq0\iff S_7\succeq0,\qquad
 H_5\succeq0\iff S_5\succeq0.
}\tag{13}
\]

The only divisions are by the forced strict positive source blocks in (10). In particular singular final certificates are retained exactly: there is no inverse of an unproved Schur pivot, no pseudoinverse, and no loss of a `gamma`-equality face. The *one remaining finite, globally quantified admission implication* is

\[
\boxed{
 \forall\theta\in\Theta_{\mathrm{C1},\,0<k\le1}^{\rm strict}
 \ \exists\ s\in\mathbb R^7,\ \gamma,\delta\in\mathbb R_{\ge0}:\quad
 \operatorname{diag}(S_7(\theta;s,\gamma,\delta),
 S_5(\theta;s),[\gamma],[\delta])\succeq0.
}\tag{14}
\]

For the original 16-feature cone set `delta=0`. The verified `k>1` counterexample disproves this original condition on unrestricted A1, **but not on the oriented domain in (14)**. The 18-feature implication (14) is also currently **unproved**, as is universal strict-C1 admission to the full sixth-root cone. Neither algebraic matching, full native stationarity, positive forced minors, connectedness nor the known C2 open-end construction proves (14).

### A forced charge bound survives the enlargement

Writing `B=2C+2D+d`, the `(0,0)` entry of `S7` is the exact gauge-independent expression

\[
 (S_7)_{00}=
 \frac{4ABu_c^2W_c^2}{2A+B}
 -\frac{(2\beta+d+2m)\gamma^2}{16m(2\beta+d)v_c^2W_c^2}.
\tag{15}
\]

It is independent of `delta` since the new gauge avoids every pivot/coupling entering this scalar. Thus the old source-rational upper bound on `gamma` remains necessary, but the two-charge extension avoids making it an upper bound on `gamma+delta`. No positivity of the other 4x4 entries follows from this scalar alone.

## Preservation, exact replay, and scope

`c1_z18_extension.py` constructs the fixed 46-by-49 coefficient map, all half-difference columns, the explicit added gauge, source-rational original-frame constructor and exact Schur matrices. `verify_z18_exact.py` proves all source-independent monomial identities, rank bounds and odd-block forcing. `verify_schur_identities.py` checks both block congruences algebraically, independently of any numerical source or SDP. The unchanged `original_baseline` contains the seven actual coquotients, full parent, six-jet stationary verification, 36-pivot coefficient map, and earlier forced positivity proof. The original source files are preserved byte-for-byte.

**Completed:** exact full-state Z18 coefficient admission for all strict C1 sources (indeed all strict A1 sources), fixed nine-gauge parametrization, and complete singularity-safe 4x4+4x4 Schur criterion. **Not completed:** a single choice of those gauges with PSD proved for every oriented strict real native source, or a proven source-wide alternative finite cone. This is an algebraic matrix-feasibility frontier, not an evaluated global affine lift. Conditional singular-source density and the independently verified C2 open-end theorem remain untouched.
