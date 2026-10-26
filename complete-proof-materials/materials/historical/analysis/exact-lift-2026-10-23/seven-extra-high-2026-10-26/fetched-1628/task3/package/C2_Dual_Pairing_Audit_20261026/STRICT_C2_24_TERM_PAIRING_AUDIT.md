# Strict native C2, fixed sixth-root dual pairing: an exact signed 18-channel reduction, and the unresolved PSD step

**26 October 2026.** This is a completed, solver-free algebraic audit of one *globally defined pairing architecture*. It **does not prove** the desired nonstrict inequality (13) and **does not disprove** membership in the full fixed sixth-root library. The earlier independently verified finite-positive-\(c\) open-end theorem and the separately reviewed qualitative spectrahedral-shadow existence theorem are preserved without alteration.

## Domain, independent tests, and the actual full parent

The native parameters are \(\theta=(a,b,c,t,r)\in\Omega\):
\[
 a,b,c,t>0,\quad 0<r<1,\quad
 h=b-a^2t>0,\quad g=(1+c^2)b-(1-r+a^2c^2t)>0,\quad
 j=1+c^2r-t>0,
\]
\[
 1-c^2h>0,\quad b-a^2c^2r>0,\quad N_*>0,\quad N_n>0.
\]
The polynomials \(N_*,N_n,\Delta\), all six native positive weights \(n,m,C,A,D,\alpha\), and the **complete four-state parent** are printed exactly in `c2_fixed_map.py` and the retained `C2_CLEARED_POLYNOMIAL_TARGET.json`. We keep \(d=1,\lambda=\beta=0\). The contact is
\[
 (u_c,v_c;W_c,X_c,Y_c,Z_c)
 =(c,ac;b/a,-cr,c(b/a-at),1).
\]
For this *fixed real source*, contact normalization \(u=cU,v=acV\) and the diagonal state normalization \(\psi=\operatorname{diag}(\psi_c)\phi\) give an identity for **every independently chosen positive real** \(U,V\) and every four-component state \(\phi\). No testing-dependent denominator is introduced. This normalization is **not** claimed to be a pure-testing-field substitution when \(c\) belongs to a separate non-real source field; the statement concerns the explicit real finite Gram-lift problem, not the qualitative separated-tensor proof.

The full target has exactly 24 nonzero monomial coefficient rows. Write \(\widehat F_\theta\) for that normalized polynomial. Its 412-row integer sixth-root contact map and its 24 integer-polynomial, source-denominator-cleared target rows are preserved unchanged alongside this note. The single common denominator \(\mathscr D>0\) and the exact compact dual spectrahedron \(\mathcal Y\) are as stated in the prior verified baseline. **The outstanding exact claim remains**
\[
 \forall \theta\in\Omega,\ \forall y\in\mathcal Y:\quad
 \sum_{k=1}^{24}P_k(\theta)y_k\ge0.\tag{G}
\]
The six contact constraints, exact rank 406, closedness and compactness are assumed from the already checked baseline. This note does not repeat their derivations.

## A smaller exact full-state affine identity

Consider the 14 integer-corner features among the pruned \(Z_6\) block:
\[
 z_i=U^{p_i}V^{q_i}\phi_{\ell_i},\qquad (p_i,q_i)\in\{0,1\}^2,
\]
with \(W\) at \((0,0)\) and \(Z\) at \((1,1)\) omitted by the C2 Newton support. Every feature equals one at the normalized contact. Their 13 differences with the first feature have 91 symmetric Gram coordinates. The exact coefficient map of their products has **78 nonzero rows and rank 72**. Of these, 54 are outside the 24 active target rows. Their exact homogeneous-zero map has rank 54, so the image of the remaining 37-dimensional nullspace on the active rows has rank **18**, exactly the 24 active coefficients modulo the six native contact constraints. An explicit 18-by-18 *rational, source-independent* minor and its 91-by-18 right inverse are stored in `C2_CORNER_AFFINE_18.json`. This proves a sharper **affine** result; it says nothing about PSD.

Among these 14 features, select all two-feature contact squares \((z_i-z_j)^2\) whose **three individual coefficient keys** lie among the 24 supported keys. There are 19 such edges, grouped into three disconnected \(K_4\) graphs and one \(K_2\):
\[
\begin{aligned}
\mathcal C_1&=\{V\phi_W,\phi_X,UV\phi_Y,U\phi_Z\},\\
\mathcal C_2&=\{U\phi_W,UV\phi_X,\phi_Y,V\phi_Z\},\\
\mathcal C_3&=\{UV\phi_W,U\phi_X,V\phi_Y,\phi_Z\},\\
\mathcal C_4&=\{V\phi_X,U\phi_Y\}.
\end{aligned}
\]
The 19 edge-square coefficient columns have exact rank **17**. Their annihilator on the 24 active coefficient rows has dimension seven: the six mandatory contact constraints, and the independent seventh relation
\[
 \mathfrak e(P)=[U^2\phi_Y^2]P-[V^2\phi_X^2]P=0.
 \tag{E1}
\]
Thus individually supported two-root contact squares cannot reproduce the full native parent whenever its \(\mathfrak e\) is nonzero. (The 14-corner **complete** Gram has additional zero-target cancellations, and does reach all 18 dimensions affinely; these two statements must not be conflated.)

The exact parent gives, **before native substitution**,
\[
 \mathfrak e(\widehat F_\theta)=
 \frac{4c^4Ah^2}{a^2}-2a^2c^4r^2(2n+1).
\]
The exact native source fractions imply the simple identity
\[
 2Ah=a^2r(2n+1).
\]
Consequently the complete native coefficient gap has the factored form
\[
 \boxed{\mathfrak e(\widehat F_\theta)=
 2c^4r(2n+1)(h-a^2r)
 =\frac{2c^4(b+r)(1+c^2r)(1+c^2r+a^2c^2t)}{\Delta}(h-a^2r).}
 \tag{E2}
\]
All displayed factors **except** \(h-a^2r\) are strictly positive throughout \(\Omega\). The remaining factor has **both** signs on \(\Omega\), checked at two previously published strict rational native points, and is not a native nonnegative gap. Thus it cannot be assigned one globally fixed nonnegative circuit multiplier.

An exact five-term, full-state *bridge polynomial* is
\[
 \boxed{\mathscr B(U,V;\phi)=
 2V\phi_W\phi_X-V^2\phi_W^2-V^2\phi_X^2
 -U^2\phi_W^2+U^2V^2\phi_W^2.}
 \tag{B1}
\]
It and all its six first derivatives vanish at the normalized contact; \(\mathfrak e(\mathscr B)=1\). By the 18-dimensional rank calculation there is an **exact source-rational, all-angle, full-four-state polynomial identity**
\[
 \boxed{\widehat F_\theta
 =\mathfrak e(\widehat F_\theta)\,\mathscr B
  +\sum_{e\in J_{17}}\gamma_e(\theta)(z_{i_e}-z_{j_e})^2.}
 \tag{B2}
\]
Here \(J_{17}\) is the fixed set of 17 independent edges selected by integer elimination, and **every** \(\gamma_e(\theta)\) is an explicitly specified source-rational linear combination of the 24 normalized target coefficients. Precisely, with \(E_J\in\mathbb Z^{24\times17}\) the selected edge columns and \(I\) the 17 selected row indices, put
\[
 \gamma_J(\theta)= (E_J)_{I,:}^{-1}
   [f_I(\theta)-\mathfrak e(f(\theta))\,b_I],
 \tag{B3}
\]
where \(b\) is the five-term coefficient vector of (B1). `C2_AFFINE_BRIDGE_17_PLUS_1.json` gives the *exact* integer row/edge indices, fixed rational 17-by-17 inverse, and the five nonzero bridge coefficients. This is an **identically valid coefficient identity**, not an SOS: \(\gamma_e\) need not be nonnegative and \(\mathscr B\) is indefinite.

## Precise failure of the natural 18-channel PSD pairing

The edge squares vanish on the four-dimensional space of feature vectors that are constant on each \(\mathcal C_i\). Compress the Gram of \(\mathscr B\) to that space, using the component indicator vectors in the order above. The **exact** compression is
\[
 Q_{\rm comp}(\mathscr B)=
 \begin{pmatrix}
 1&0&-1&0\\0&-1&0&1\\-1&0&1&0\\0&1&0&-1
 \end{pmatrix},\qquad
 \operatorname{spec}Q_{\rm comp}=(2,-2,0,0).
 \tag{B4}
\]
Every 19-edge Gram has **zero** compression on this space. Therefore **for each** native source with \(h-a^2r\ne0\), *every Gram in this exact 19-edge-plus-bridge representation is indefinite*, whatever two graph-edge gauges are chosen. This is a proof that this **particular simple pairing architecture** cannot close (G). It is **not** a separation witness for the full 154-feature cone and not an obstruction to another 14-corner allocation using cancellations outside the active rows. A correct finite-degree proof must include coupled, off-support-canceling PSD-block transfers; the additional fractional rows of \(Z_6\) are available for this purpose. The saved full-domain open-end gauge/tree construction already does so in its own source chart, but it does not establish the all-interior transfer.

## Exact remaining implication

No parameter-dependent collection of PSD Gram blocks with nonnegative native-gap multipliers has been constructed that matches all **412** coefficient equations over \(\Omega\). In particular (B3) has **no proved factor-sign properties** sufficient for PSD, and the bridge (B1) cannot be retained as a positive square. Hence the first *global* unproved condition remains **(G), nonstrictly**. To close it one must exhibit a full 154-feature PSD pairing whose fractional/zero-coefficient transfer terms repair (B4) and every other remaining Schur direction, or prove a different finite globally sound admission cone. Affine compatibility, native physical curvature, strictness/nondegeneracy continuation and representation of each pseudo-moment by a positive testing measure are **not** substituted for (G). No counterexample to the full \(L_6\) cone is claimed.

## Verification and dependencies

`verify_corner_pairing.py` replays all claims (14-feature 72-rank, supported 18-rank, 19-edge rank17, seventh functional (E1), native factor (E2), entire 18-channel coefficient identity (B2) for the complete 24-dimensional stationary-contact target, graph quotient matrix (B4) and the sign change). It uses exact SymPy rational arithmetic, and *no* solver or sampling. The only mathematical inputs are the original full C2 parent/native source formulas and the previously verified 24-key and 154-support JSON data. Its output files are `C2_CORNER_AFFINE_18.json` and `C2_AFFINE_BRIDGE_17_PLUS_1.json`.

The independently checked baseline source map, original full polynomial target and original native source generator are copied verbatim in this folder, so no absent binary intermediates or GitHub access is required. The qualitative finite-SDP existence proof and strict open-end theorem remain separate completed inputs, unchanged.
