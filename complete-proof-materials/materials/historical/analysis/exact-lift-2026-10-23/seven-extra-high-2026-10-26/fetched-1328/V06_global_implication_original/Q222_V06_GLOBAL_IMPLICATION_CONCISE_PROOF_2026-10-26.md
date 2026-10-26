# V06. Conditional qualitative finite-SDP theorem: complete global implication

**Frozen inputs.** Main `8108102c092894c62768b7c8b0008889737bf1e7`; its eight foundational TeX chapters compared for this audit agree byte-for-byte with `332483309586e0fec6e16fc1dd855915b648dca1`. Primary external theorem: C. Scheiderer, *Second-order cone representation for convex sets in the plane*, SIAM J. Appl. Algebra Geom. 5 (2021), 114–139, arXiv:2004.04196v2, **Corollary 3.19** (printed p. 15, PDF page index 14). Sections and formula tags below refer to the current repository. The manuscript's native component proof is **not** being re-proved here.

## The one labelled local input

**(NATIVE)** For every pair of independently labelled real closed fields \(K_c,K_t\supseteq\mathbb R\), every physical canonical source \(c\) of either of the following supports, possessing its actual finite noncommuting stationary zero,
\[
\begin{array}{ll}
\mathsf C_1:&\lambda=n=0,\ d,m,C,A,\alpha>0,\ D,\beta\ge0,\ D+\beta>0;\\
\mathsf C_2:&\lambda=\beta=0,\ d,m,C,A,D,\alpha>0,\ n\ge0,
\end{array}
\]
with all its genuine native source gaps, including all permitted third-omission faces and infinitesimal opening rates, and every **independent** projective pair of tests and arbitrary real state \(\psi\in K_t^4\), the full homogeneous nine-packet parent \(\mathcal F_c\) of (RFA10)/(LC2) belongs to \(\Sigma(K_c\otimes_{\mathbb R}K_t)^2\). Real/imaginary splitting supplies complex states. This is exactly the universal native atom assertion of `07_cover.tex`; neither a fixed feature list nor a common number/rank of squares is part of (NATIVE). The local evidence for (NATIVE) is reviewed separately, and its proof is **not an input** to any lemma below.

**Theorem (V06).** Assume (NATIVE). Then the full eight-coordinate finite-dimensional tensor-product quantum set \(\mathcal Q_{222}\subset\mathbb R^8\) is the projection of a **finite real affine spectrahedron**. All steps after (NATIVE), using the already displayed non-native elementary/closed-face results, are justified below. No explicit pencil, matrix order or semidefinite-extension-degree bound follows.

### Lemma 1. Nine real-qubit atoms, both directions, and first-order transfer

Let \(\Omega=S^1\times S^1\times S^3\) and \(\Phi:\Omega\to\mathbb R^8\) be the polynomial expectation map of (RFA2). Then
\[
\mathcal Q_{222}=\operatorname{conv}\Phi(\Omega)
=\bigl\{\sum_{i=1}^{9}w_i\Phi(\omega_i):w_i\ge0,\ \sum_iw_i=1,\ \omega_i\in\Omega\bigr\}.\tag{1}
\]
Indeed simultaneously dilate both binary effects of each party to projections on the *same* doubled Hilbert space with common isometry; Jordan's lemma decomposes each pair of involutions into blocks of dimension at most two; diagonal blocks of the bipartite state give nonnegative mixture weights; real Bell matrices split every complex pure state \(r+is\) into the real contributions of \(r\) and \(s\). Conversely, each atom is a physical two-qubit behaviour and any finite convex mixture is realized by local direct sums and a matching block-diagonal state. Carathéodory in \(\mathbb R^8\) reduces every mixture to nine atoms, including zero weights.

The right side of (1) is the continuous image of a compact real simplex/sphere product, hence is compact, closed and convex. Its finite polynomial-existential description makes it semialgebraic. **Define** \(Q_K\) for any real closed \(K\supseteq\mathbb R\) by this same nine-atom formula. This is first-order transfer, not an assumption of compactness in the order topology of \(K\). (Source: `02_model.tex`, RFA1–RFA5.)

### Lemma 2. The transferred affine polar has the correct quantifiers

Let \(P_Q=\{f=f_0+\sum_{\nu=1}^8 f_\nu x_\nu\in\mathbb R[x]_{\le1}:f\ge0\text{ on }Q\}\). The defining sentence
\[
\forall q\,[q\in Q\Rightarrow f_0+\sum_\nu f_\nu q_\nu\ge0]
\]
is first-order over the real field because membership of \(Q\) is given by (1). Tarski transfer gives, **for each** real closed \(K\supseteq\mathbb R\),
\[
(P_Q)_K=\{f\in K[x]_{\le1}:f(q)\ge0\ \forall q\in Q_K\}.\tag{2}
\]
Thus no narrower choice of affine inequalities or illicit identification of two coefficient/test factors enters the criterion. (RFA6; TF11.)

### Lemma 3. Every transferred valid affine functional has finitely many certified source atoms

A valid affine \(f\) is nonnegative on all local deterministic behaviours. If \(f\) is nonnegative at every one of the eight PR vertices, it is nonnegative on the entire nonsignalling polytope, whose only vertices are the 16 deterministic and eight PR behaviours. By the finite probability-inequality description and Farkas duality, \(f\) is a nonnegative linear combination of the 16 outcome probabilities. Each has an immediate rank-one separated square for every test.

Otherwise choose a PR vertex at which \(f<0\), apply one of the finitely many local relabelings, and use the affine-basis identity (SG1):
\[
f=d\,\bigl(1-\tfrac12\mathrm{CHSH}\bigr)+4\sum_{e\in\{\lambda,n,m,C,A,D,\alpha,\beta\}}a_e p_e,
\quad d=-f(\mathrm{PR})>0,\quad a_e=f(v_e)/4\ge0.\tag{3}
\]
No assertion that this is an SOS is made here. Its normalized physically valid coefficient body \(\{(d,a_e)\ge0:\sum d+\sum a_e=1\}\) is a compact real convex semialgebraic set of affine dimension at most eight. By finite-dimensional Krein–Milman and Carathéodory, it is a convex combination of at most nine extreme points. The *fixed-nine* conclusion and the midpoint definition of extremality are first-order formulas (the source body's membership is the universally quantified physical-atom test), so this conclusion transfers to all \(K_c\). Multiply by the original positive pure-source total scale to recover \(f\).

The operative extreme-source classification of `04_geometry.tex` is: \(d=0\) gives events; \(d>0\) with noncoercive Laurent support belongs, up to finite symmetries, to one of the two **closed faces** \(\lambda=m=0\) or \(\lambda=\alpha=0\); \(d>0\) with coercive support attains a finite stationary zero and has exactly five or six positive event weights, hence, by support incidence, is a C1/C2 native contact or permitted third omission (SG7–SG10, SG17–SG22). The stationary zero gives global **physicality** from the balanced Laurent identity and the positive three-state block; it does **not** by itself give the required separated full-state SOS. Source/angle gauges and strict contact signs are carried by `05_entry.tex`, EN1–EN16. Every displayed classification and finite decomposition is a first-order real-algebraic implication, so applies to \(K_c\). A nonnegative decomposition of a source with a literal zero coordinate forces that coordinate to vanish in every summand; it is not silently opened.

### Lemma 4. Elementary and closed-face atoms already have exact separated certificates

An outcome event is a pure-testing nonnegative rank-one probability square multiplied by a nonnegative pure-source coefficient. For a closed-face source, SG11 gives the exact univariate quadratic **matrix** Gram criterion, including singular faces. Its applications SG12 (one \(6\times6\) Gram) and SG13 (two \(3\times3\) Grams and one auxiliary) are necessary and sufficient for physicality on their corresponding exclusive-omission faces. Both homogeneous source pencils have explicit strictly feasible points after incorporating the coefficient-sign blocks. This ensures their dual PSD images are closed and every nonnegative testing-evaluation functional has an **attained** dual PSD witness. The existence of these finite-size witnesses is first-order, hence transfers to \(K_t\).

Let \(A_c\succeq0\) be a source pencil over \(K_c\) and \(B_t\succeq0\) its testing dual witness over \(K_t\). Auxiliary coefficients vanish by dual orthogonality. The exact slack equals the separated pairing
\[
\langle A_c,B_t\rangle_\otimes
=\sum_{i,j}\left(\sum_k R_{ki}\otimes S_{kj}\right)^2,
\quad A_c=RR^T,\ B_t=SS^T.\tag{4}
\]
Both factorizations are performed in the respective real closed field, including singular matrices. This proves the closed-face atom certificates, not merely pointwise nonnegativity. (SG11–SG14; TF9.)

### Lemma 5. The strict bounded tensor-Gram mechanism is genuinely algebraic

For any real closed \(K_c,K_t\), put \(T=K_c\otimes_{\mathbb R}K_t\); **no order and no field inverse on \(T\)** are assumed. The bounded pure-field rings \(B_c,B_t\) have standard-part maps to \(\mathbb R\). Their tensor ring \(B=B_c\otimes_{\mathbb R}B_t\subset T\) has residue \(\rho:B\to\mathbb R\), whose kernel is \(\mathfrak m_c\otimes B_t+B_c\otimes\mathfrak m_t\). If \(\rho(a)=\epsilon>0\), write \(a=\epsilon+\sum_j x_j\otimes y_j\) where one factor of each elementary summand is infinitesimal, and split \(\epsilon=\sum_j\epsilon_j\), \(\epsilon_j>0\). When \(x_j\) is infinitesimal and \(M_j>|y_j|\) is real,
\[
\epsilon_j+x_j\otimes y_j=
\frac{(\epsilon_j/M_j+x_j)\otimes(M_j+y_j)
+(\epsilon_j/M_j-x_j)\otimes(M_j-y_j)}2.\tag{5}
\]
All factors on the right are strictly positive in their **own** fields, so their separate square roots give an exact finite SOS in \(T\). If the other factor is infinitesimal, swap roles.

If \(Q\in\operatorname{Sym}_r(B)\) and \(\rho(Q)\succ0\), use a **real** congruence so that \(\rho(Q)=I\), choose a small positive real \(\delta\) with \((r-1)\delta<1\), and write
\[
Q=\sum_i[Q_{ii}-(r-1)\delta]e_ie_i^T+
\sum_{i<j}\left[\frac{\delta+Q_{ij}}2(e_i+e_j)(e_i+e_j)^T+
\frac{\delta-Q_{ij}}2(e_i-e_j)(e_i-e_j)^T\right].\tag{6}
\]
Every scalar coefficient has positive real residue and hence is SOS by (5). Therefore \(g^TQg\in\Sigma T^2\) **for any tensor root vector \(g\in T^r\)**, bounded or not. This is stronger than checking a positive value or a positive-definite real evaluation matrix. (TF1–TF6.)

For a pure-source \(s>0\) and pure-testing \(t>0\), multiplication or division by \(s\otimes t\) preserves SOS, using \(s^{\pm1/2}\otimes t^{\pm1/2}\). No general mixed tensor sum may be divided by or square-rooted. Homogeneous state/angle coordinates must be nonzero before division; zero states and literal projective endpoints are handled by their original identities. (TF8; RFA10–RFA13.)

### Lemma 6. From the exact parent to each independent atom slack

For projective angle pairs \((a_0,a_1),(b_0,b_1)\in K_t^2\setminus\{0\}\) and nonzero \(\psi\in K_t^4\), let \(\Delta_A=a_0^2+a_1^2>0\), \(\Delta_B=b_0^2+b_1^2>0\), \(N=\|\psi\|^2>0\). The complete nine-packet homogeneous polynomial obeys the **identity in** \(K_c\otimes K_t\)
\[
\mathcal F_c=(1\otimes\Delta_A\Delta_BN)\,[f_c,\Phi(a,b,\psi/\sqrt N)].\tag{7}
\]
RFA10–RFA12 check every event and signed CHSH packet, including affine constants. For event/face atoms use Lemma 4; for each native contact atom apply the one labelled assumption (NATIVE). Then multiply every square root by \(1\otimes(\Delta_A\Delta_BN)^{-1/2}\), a **pure-testing** square root. Positive source normalization is likewise restored in \(K_c\). No mixed inversion is used. The zero state is trivial.

To check all entry conditions without strengthening (NATIVE): the auxiliary noncommuting-residue theorem `noncommuting_refined.tex` applies to an **already physical** nonnegative-coefficient source and an independent test where **both** local pairs are noncommuting in residue. It does not cover infinitesimal nonzero incompatibility. At a **literally commuting** pair the CHSH operator satisfies \(H\succeq0\) over the testing field, so CHSH and event terms separately give tensor squares; the literal homogeneous zero/infinite axes are included. For positive half-angles with tested state \(WZ\le0\), the exact pure-testing Schur expression (SG15) makes the CHSH core nonnegative, again an SOS for that *specific* tested state. When \(WZ>0\), \(Z\ne0\) and a pure-testing state normalization gives \(Z=1,W>0\). These elementary alternatives do not purport to certify singular infinitesimal-angle contacts. SG7's determinant/Laurent argument establishes physicality and stationarity, not the separated full-state Gram. Those contacts are exactly why (NATIVE) is required.

### Lemma 7. Two finite decompositions preserve the affine constant and zero weights

Let \(f=\sum_{j=1}^{J}\alpha_j f_j\) be the source decomposition in \(K_c\), with \(J\le9\), \(\alpha_j\ge0\) (plus finitely many event atoms if using the nonsignalling branch). Let \(q=\sum_{i=1}^I w_i\Phi(\omega_i)\) be (1) in \(K_t\), with \(I\le9\), \(w_i\ge0\) and **\(\sum_iw_i=1\)**. If the resulting independent atom slacks are \([f_j,\Phi(\omega_i)]=\sum_k g_{ijk}^2\), then the full source-first tensor evaluation is exactly
\[
[f,q]=\sum_{i,j,k}\bigl[(\sqrt{\alpha_j}\otimes\sqrt{w_i})g_{ijk}\bigr]^2.\tag{8}
\]
The equality of constant terms is \(\sum_{ij}\alpha_jw_i f_{j,0}=\sum_j\alpha_jf_{j,0}\), which uses the testing weight sum **one**. Source conic weights need not sum to one. Neither zero source nor zero testing weights are inverted; they simply contribute zero squares. This proves
\[
\forall K_c,K_t,\ \forall f\in(P_Q)_{K_c},\ \forall q\in Q_{K_t},\quad
[f,q]\in\Sigma(K_c\otimes_{\mathbb R}K_t)^2.\tag{9}
\]

### Lemma 8. Correct qualitative criterion and both pencil inclusions

Choose **two labelled copies of one arbitrary** real closed \(K\supseteq\mathbb R\) in (9). By Lemma 2 the source functional ranges over \((P_Q)_K\), and the test ranges over \(Q_K\). Scheiderer's *Corollary 3.19* states exactly that **for all such** \(K,f,q\), unrestricted finite SOS of the point-first tensor evaluation in \(K\otimes_{\mathbb R}K\) characterizes a finite spectrahedral shadow of any **closed convex semialgebraic** \(Q\). The flip of the two labelled tensor factors is an algebra isomorphism preserving sums of squares, so (9) has precisely the paper's ordering. The distinct *Theorem 3.10* prescribes \(\operatorname{sxdeg}(Q)\le d\) by additionally bounding the tensor rank of each square by \(d\); the **corollary does not** require any common \(d\), square count, polynomial degree, finite NPA level or global feature library.

For clarity, here is a finite pencil with **both inclusions**, contingent only on the finite-dimensional existence conclusion in the *proof* of Corollary 3.19, not an alleged explicit extraction from local charts. That proof supplies an affine real variety \(X\), a real algebraic map \(\phi:X\to\mathbb A^8\), and a finite-dimensional \(U\subset\mathbb R[X]\) such that
\[
Q\subseteq\phi(X(\mathbb R)),\qquad \phi^*(P_Q)\subseteq\Sigma U^2.\tag{10}
\]
Choose a basis \(u_1,\ldots,u_r\) of \(U\), enlarge by \(1\) if convenient, and set
\[
V=\operatorname{span}_{\mathbb R}\bigl(1,\phi_1,\ldots,\phi_8,\{u_iu_j:1\le i\le j\le r\}\bigr)\subset\mathbb R[X].
\]
The proposed fixed spectrahedral lift uses **linear functional variables** \(L\in V^*\):
\[
L(1)=1,\quad L(\phi_\nu)=q_\nu\ (1\le\nu\le8),\quad
M(L)=[L(u_iu_j)]_{i,j=1}^r\succeq0.\tag{11}
\]
The coordinate-ring linear dependencies in \(V\) give finitely many exact linear equations; eliminate them or impose them with opposite scalar PSD blocks. Choosing a basis of the finite real space \(V\) yields one **finite real affine matrix pencil** with finitely many auxiliary variables.

**\(Q\subseteq\mathrm{proj}(11)\).** For \(q\in Q\), (10) gives \(x\in X(\mathbb R)\) with \(\phi(x)=q\). Evaluation \(L=\mathrm{ev}_x|_V\) is feasible, because its moment matrix is \(u(x)u(x)^T\succeq0\).

**\(\mathrm{proj}(11)\subseteq Q\).** If \((q,L)\) is feasible and \(f\in P_Q\), then \(\phi^*f=\sum_j(\sum_i c_{ji}u_i)^2\) by (10). Thus \(f(q)=L(\phi^*f)=\sum_j c_j^TM(L)c_j\ge0\). Since \(Q\) is closed and convex, finite-dimensional affine separation gives \(q\in Q\). No representing measure, moment-completeness assertion, or Slater assumption is being inserted here. In particular (10)'s weaker inclusion \(Q\subseteq\phi(X(\mathbb R))\) suffices for the first direction, while the SOS pullback and separation prove the second.

The theorem follows. The matrices in (11) are **existential** because Corollary 3.19 only establishes existence of \(X,U\). The passage \(Q\) semialgebraic over \(\mathbb Q\) and a fixed, but unspecified, finite pencil size additionally allows algebraic real matrix entries by elementary real-closed-field transfer; it supplies **no evaluated entries or bound on** \(r,N,m\).

## Scope of the verdict

**All global implications are valid conditional on (NATIVE).** The demonstrated non-native prerequisites are (1) the exact real-qubit/nine-atom model and its first-order extension, (2) the canonical extreme-source support reduction, including the zero-event and d=0 branches, and (3) the attained SG12/SG13 closed-face dual certificates. The strict tensor-Gram mechanism, parent normalization, noncommuting-residue entry, pure test alternatives, conic/convex assembly, and the Corollary 3.19 interface match their exact hypotheses. The evidence does not replace an independent audit that (NATIVE) is actually proved at **every** physically permitted C1/C2 source and singular/projective test. Neither an evaluated fixed lift nor a finite-NPA theorem is asserted.
