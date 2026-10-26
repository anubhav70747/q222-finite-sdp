# Q222: coefficient-space density of singular canonical contact sources

**26 October 2026 | focused conditional construction theorem**

**Sources inspected:** frozen source package `07_FIXED_CONE_DENSITY_sources.txt` (SG1–SG22, new strict C1/C2 native domains, and the fixed sixth-root cone); primary original `latex/chapters/05_entry.tex` at public commit `aeea5ef87892c9897a27cb4c288d2f835f7e88df`, especially EN1–EN16. The references to exact PSD admission for *strict* C1/C2 sources below are **hypotheses**, not conclusions. No numerical or SDP search is used.

## 1. Statement, including the essential convex-assembly qualification

Write the canonical source vector as
\[
 c=(d,\lambda,n,m,C,A,D,\alpha,\beta),\qquad
 F_c=dF_H+4\sum_{i=1}^{8}a_i E_i,
\]
where each `E_i` is its fixed nonnegative squared-event polynomial in the complete four-state parent SG2. Define the physical source cone \(\mathscr P\) by positivity on every quantum behaviour and the normalized canonical section
\[
 \mathscr B=\{c\in\mathscr P:\ d,a_i\ge0,\ d+\textstyle\sum_i a_i=1\}.
\]
Let \(\mathscr S_1^\circ\) be the actual strictly positive **six-event** physical native C1 sources, in the full strict domain **A1** of the supplied C1 reduction (\(\lambda=n=0\), \(m,C,A,D,\alpha,\beta>0\)). Let \(\mathscr S_2^\circ\) be the actual strict physical C2 sources in the full **\(\Omega\)** of the supplied C2 reduction (\(\lambda=\beta=0\), \(n,m,C,A,D,\alpha>0\)). Both are considered modulo positive source scaling so that they lie in \(\mathscr B\). Let \(\mathscr S_i\) also allow their **one additional permitted event omission** while retaining all the original strict noncommuting-contact gap hypotheses:
\[
\begin{array}{ll}
\mathscr S_1:\;&m,C,A,\alpha>0,\quad D,\beta\ge0,\ D+\beta>0,\quad \lambda=n=0,\\
\mathscr S_2:\;&m,C,A,D,\alpha>0,\quad n\ge0,\quad \lambda=\beta=0.
\end{array}
\]
These are physical **finite noncommuting contact** sources; the two source-only zero faces with noncoercive support are *not* included in \(\mathscr S_i\).

**Theorem 1 (exact singular-contact entrance).** In the ordinary coefficient topology of \(\mathbb R^9\),
\[
 \boxed{\mathscr S_1\subseteq\overline{\mathscr S_1^\circ},\qquad
 \mathscr S_2\subseteq\overline{\mathscr S_2^\circ}.}\tag{D1}
\]
Indeed every five-event point admits a real-analytic positive one-event opening, with a simultaneous real-analytic adjustment of its Bell coefficient and complete source contact. The other omitted event weights remain **literally zero**. The perturbed sources satisfy the actual strict native inequalities A1 or \(\Omega\), not merely their formal stationarity equations.

Let \(\mathcal K\subseteq\mathscr P\) be **one fixed closed convex globally sound cone in Bell-coefficient space**, defined by a finite source-independent sixth-root Gram map, and assume *conditionally*
\[
\boxed{\mathscr S_1^\circ\cup\mathscr S_2^\circ\subseteq\mathcal K.}\tag{H}
\]
Then (D1) yields \(\mathscr S_1\cup\mathscr S_2\subseteq\mathcal K\), including \(D=0\), \(\beta=0\), and \(n=0\), without any separate singular-source Gram construction.

Let \(\mathscr E\) denote the physical event-ray cone, \(\mathscr F_{12}\) the complete original physical closed face \(\lambda=m=0\) (SG12), and \(\mathscr F_{13}\) the complete original physical closed face \(\lambda=\alpha=0\) (SG13), including all their singular matrices and further omissions. Let \(\mathfrak G\) be a finite set of the required physical relabellings, including all PR orientations and the dihedral source orbits. **Do not assume** that an arbitrarily chosen fixed feature subcone is itself invariant under every relabelling; use its explicitly relabelled coefficient cones.

**Theorem 2 (exact, conditional conic assembly).** Under (H), the *entire* cone of quantum-valid affine functionals satisfies
\[
 \boxed{\mathscr P=\mathscr E+\sum_{g\in\mathfrak G}
  \bigl(g\mathcal K+g\mathscr F_{12}+g\mathscr F_{13}\bigr).}\tag{D2}
\]
The finite Minkowski sum on the right is **closed**. Consequently, conditional on (H) and the existing exact SG12/SG13 pencils and their dual-attainment hypotheses, its dual intersection is a finite SDP representation of the quantum behaviour set. This is a **conditional construction implication**, not a proof of (H), the universal 24-term C2 dual inequality, the C1 eight-gauge PSD inequalities, or any already-unverified explicit global lift.

The closure of the **union** of the strict families, events and face sources is *not* all of \(\mathscr P\): a general strictly positive coefficient vector need not have either forced two-event omission. Finite **convex sums** are indispensable. This distinction is exactly why (D2), rather than union-density, is the correct global theorem.

## 2. Complete curvature and the analytic perturbation

Fix a real point \(c_0\in\mathscr S_1\cup\mathscr S_2\) with **five** positive events, and its positive finite contact \((u_c,v_c,W_c,X_c,Y_c,Z_c=1)\). Choose its original PR orientation and contact signs from EN1–EN2; none is assumed by continuity from an unrelated source chart.

The four Laurent coefficients indexed by the axial exponent vectors \((-1,0),(0,-1),(0,1),(1,0)\) are, from the actual SG6 determinant,
\[
\begin{aligned}
q_{01}&=4(\lambda D t_7t_6+\alpha m t_5t_8)>0,\\
q_{10}&=4(\alpha A t_7t_3+\lambda n t_8t_4)>0,\\
q_{12}&=4(\beta m t_1t_6+CDt_5t_2)>0,\\
q_{21}&=4(\beta n t_3t_2+CA t_1t_4)>0.
\end{aligned}\tag{D3}
\]
The strict signs use only \(m,C,A,\alpha>0\), and \(D+\beta>0\) in C1 or \(D>0\) in C2. In particular, all four remain positive when a third event is **literally absent**. The contact balance of SG7 gives \(\sum_e w_e e=0\), with \(w_e=q_eu_c^{2e_1}v_c^{2e_2}\). Therefore the angular Hessian in logarithmic coordinates is \(\sum_e w_eee^T\succ0\), directly because these four positive axial terms span both directions.

SG3 supplies a positive definite \((W,X,Y)\) state block for \(d_0>0\). Completing these three variables within the source field gives a scalar angular Schur form proportional to the determinant Laurent polynomial (SG19); at the stationary zero its angular Hessian is strictly positive by (D3). Thus the **full five-variable Hessian**
\[
 H_0=\nabla_s^2F_{c_0}(s_0)\succ0,\qquad s=(W,X,Y,u,v)\tag{D4}
\]
is positive definite, including the *state-angle mixed blocks*. There is no need for a uniform Hessian margin along a sequence of different sources.

Write \(E_j(s)=\partial_{a_j}F_c(s)\ge0\) for the actual event polynomial of the **single** additional omission to be opened. The signed core \(F_H=\partial_dF\) is strictly negative at the contact. Indeed, if all positive event amplitudes vanished there, the stationary state-kernel equation would imply \(P_H(u_c,v_c)\psi_c=0\). But SG20 gives \(\det P_H=-16D_u^2D_v^2u_c^2v_c^2\ne0\). Hence some positive event amplitude is nonzero and
\[
0=F_{c_0}(s_0)=d_0 F_H(s_0)+\sum_i a_iE_i(s_0),\qquad F_H(s_0)<0.\tag{D5}
\]

Hold all old event coefficients fixed. Put the new weight \(a_j=\varepsilon\), the other forced omitted weights equal to zero, and solve
\[
\mathscr J(s,d,\varepsilon)=
 (\nabla_sF_{d,a+\varepsilon e_j}(s),\ F_{d,a+\varepsilon e_j}(s))=0.\tag{D6}
\]
The six-by-six Jacobian in \((s,d)\) at \((s_0,d_0,0)\) is exactly
\[
 D_{(s,d)}\mathscr J=
 \begin{pmatrix}H_0&\nabla_sF_H(s_0)\\0&F_H(s_0)\end{pmatrix},
 \qquad\det D_{(s,d)}\mathscr J=\det(H_0)F_H(s_0)\ne0.\tag{D7}
\]
Thus the implicit-function theorem gives an analytic real branch \(\varepsilon\mapsto(s(\varepsilon),d(\varepsilon))\) near zero. In particular,
\[
 d'(0)=-\frac{E_j(s_0)}{F_H(s_0)}\ge0,\qquad
 s'(0)=-H_0^{-1}[\nabla_s E_j(s_0)+d'(0)\nabla_sF_H(s_0)].\tag{D8}
\]
These derivatives are **pure-source calculations**, not Gram multipliers. For each sufficiently small \(\varepsilon>0\), \(d(\varepsilon)>0\), all event weights are nonnegative, the original positive event weights remain positive, and \(u(\varepsilon),v(\varepsilon)>0\) with \(Z=1\). Homogeneity plus (D6) supplies the fourth state derivative. The SG7 balanced-Laurent contact lemma then establishes physicality at **every** angle/state. This step is essential: generic coefficient perturbations of a boundary inequality need not stay valid.

If the opened event already vanishes at the old contact, \(E_j(s_0)=\nabla_sE_j(s_0)=0\); one may even choose \(s(\varepsilon)=s_0\), \(d(\varepsilon)=d_0\), and add its nonnegative square without moving the contact. Formula (D6) handles both possibilities uniformly.

## 3. The original native coordinate entrance, without opening losses

**C1.** In the EN13 chart, set \(k=u_c\), \(S=W_c\), \(r=Y_c/(kS)\), \(\theta=1+X_c/k\), \(\eta=v_c/(kS)\), \(\zeta=Sv_c/k\), \(a=(1-r)/\eta\), \(t=\theta/\zeta\), and \(f=\sqrt\eta,g=\sqrt\zeta\). These invert the supplied native C1 A3 contact exactly. At the five-event boundary EN1–EN2 give
\[
 0<r,\theta<1,\quad
 R_m=k\zeta M>0,\quad R_D=-v_cJ<0,
 \quad W_c+v_cX_c=SE>0,\quad1-v_cY_c=Q>0.\tag{D9}
\]
Thus \(M,J,E,Q,1-t\zeta\) are strict positive source gaps. The identities \(L_c=r\), \(H_c=r+\eta\), \(Q_c=Q\), and \(1-\kappa\eta+tV_c=E\) give the other explicit A1 gaps; \(P_4,P_\alpha,P_d>0\) follow from the printed positive NC2 formulas. The last two A1 inequalities are **exactly** \(D>0\) and \(\beta>0\) after clearing positive native denominators. Opening whichever one was zero therefore places every nearby source (D6) in **the original strict A1**, with the other missing \(\lambda=n=0\) unchanged.

**C2.** Use the exact EN14 inverse
\[
 c=u_c,\ a=v_c/u_c,\ b=aW_c,\ r=-X_c/u_c,\quad
 t=(u_cW_c-Y_c)/v_c.\tag{D10}
\]
Then, identically,
\[
 h=aY_c/u_c>0,\quad g=R_m/u_c>0,\quad j=-R_D/v_c>0,
 \quad 1-c^2h=1-v_cY_c>0,
 \quad b-a^2c^2r=a(W_c+v_cX_c)>0.\tag{D11}
\]
The source geometry gives \(0<r<1\), and all native denominator signs are strict (EN15–EN16). The stationary six-column determinant stays nonzero; the original rational stationary formulas consequently recover the same weights. Its \(N_*\) stays strictly positive, and \(N_n=2r\Delta\,n\). Thus opening \(n=0\) to \(n=\varepsilon>0\) gives **exactly** \(N_n>0\) in the original strict C2 \(\Omega\), while \(\lambda=\beta=0\) stay literal.

This also handles the common five-event C1/C2 overlap \(\lambda=n=\beta=0\), either by opening \(\beta\) on C1 or \(n\) on C2. There is **no** assumption that a limiting state coordinate stays a positive real *unit*; each original finite contact has strictly nonzero coordinates and the entrance is local to that individual coefficient point.

Finally, for \(c_0\in\mathscr B\), normalize the perturbed source by the **positive source-only** scalar
\[
 \overline c(\varepsilon)=\frac{(d(\varepsilon),a+\varepsilon e_j)}{
obreak d(\varepsilon)+\sum_i a_i+\varepsilon}
 \longrightarrow c_0.\tag{D12}
\]
If a source formula uses \(d=1\) (C2), first divide by \(d(\varepsilon)>0\); if it uses \(C=2\) (C1), first scale to that section. Positive homogeneity makes all normalizations equivalent. This proves (D1).

## 4. No hidden singular-contact stratum

SG9 establishes that a **coercive** physical finite-contact source cannot omit an incompatible pair of losing events, and compatible omissions number at most three. The extreme-source lemma in SG10 proves that an extreme normalized coercive source with \(d>0\) has **five or six** positive event weights. EN1–EN16 put every such extreme point, after a physical relabelling, in C1 or C2 with the specified strict native gaps. We have just proved that their five-event boundaries lie in the coefficient closure of the six-event strict loci.

A source with *further omissions* has noncoercive Newton support (it omits an exclusive pair) and falls in the two exact **closed-face** orbits SG12/SG13. Physical projective/commuting zero configurations also give no new extreme class: if an extreme has coercive support, its Laurent polynomial attains its zero at a **finite noncommuting** contact; otherwise it belongs to the noncoercive faces. At \(d=0\) the entire source is a nonnegative combination of event rays. Limits in which a contact coordinate tends to zero or infinity while the coefficients converge are automatically covered by this dichotomy applied to the *limiting coefficient source*. The limiting contact parameters are **not** required to converge.

An arbitrary valid source, including one with seven or eight positive event weights, is not asserted to lie in the closure of one strict family. Rather, the SG extreme decomposition supplies at most **nine** extreme-source atoms in the normalized canonical base; each belongs to the strict-contact closure, a closed face, or an event cone. Its *finite nonnegative sum* lies in the corresponding sum of certificate cones. The claimed universal source decomposition over real closed fields is the first-order finite-atom SG result, not a new compactness assumption over a non-Archimedean field.

## 5. Closedness of the **fixed** Gram cone and of the finite source-cone sum

For the full fixed sixth-root library \(Z_6=(u^{i/6}v^{j/6}\psi_l)\), \(0\le i,j\le6\), the original Gram cone is
\[
 \mathcal K_6=\{c:\ F_c=Z_6^TGZ_6,\ G\succeq0\}.\tag{D13}
\]
Fix once for all the product probability measure on \([1/2,2]^2\times S^3\) with strictly positive density, and form \(M=\int Z_6Z_6^T\,d\mu\). The feature functions are independent: distinct positive fractional powers are linearly independent on the open angle rectangle, and the four state coordinates separate linearly. Consequently \(M\succ0\). If \(c_n\to c\) and \(F_{c_n}=Z_6^TG_nZ_6\), \(G_n\succeq0\), then
\[
 \lambda_{\min}(M)\operatorname{tr}G_n
 \le\operatorname{tr}(MG_n)=\int F_{c_n}\,d\mu,
\tag{D14}
\]
whose right side is a **bounded linear function** of the *original coefficients*. Thus \(\operatorname{tr}G_n\) is bounded; a convergent Gram subsequence proves \(c\in\mathcal K_6\). This proves coefficient-space closedness with no moving-contact-coordinate assumption. The same argument works for any fixed independent-feature sublibrary (after removing linear dependence); a general fixed closed cone \(\mathcal K\) is already assumed closed.

The sum in (D2) also is closed, although arbitrary sums of closed cones need not be. Write any Bell affine form in the deterministic-sign basis as
\[
 f(\chi)=t+\sum_{i=1}^8b_i\chi_i,
 \qquad\chi=(A_0,A_1,B_0,B_1,A_0B_0,A_0B_1,A_1B_0,A_1B_1).
\]
Average over all 16 deterministic assignments of \(A_0,A_1,B_0,B_1\in\{\pm1\}\). The nine sign monomials including 1 are orthonormal. For **every globally sound** source \(f\ge0\) at these vertices,
\[
 t=\mathbb Ef\ge0,\qquad b_i=\mathbb E[\chi_i f],\quad |b_i|\le t.\tag{D15}
\]
Thus \(t=0\) forces \(f=0\), and \(\{f\in\mathscr P:t(f)=1\}\) is bounded. For a convergent sequence \(f_n=\sum_{j=1}^N f_{j,n}\), \(f_{j,n}\) lying in **fixed closed globally sound** cones, \(0\le t(f_{j,n})\le t(f_n)\). Equation (D15) bounds every \(f_{j,n}\); a diagonal subsequence converges term by term to elements of those closed cones. Hence their finite Minkowski sum is closed. Relabelled copies stay closed because relabellings are invertible linear maps of Bell coefficients. The source faces are closed as intersections of \(\mathscr P\) with their exact coordinate hyperplanes (and have the existing exact SG pencils); the event cone is polyhedral and closed.

Nothing here divides by a contact gap, takes a mixed inverse, or assumes a positive **uniform** Gram eigenvalue as a boundary is approached.

## 6. Conditional full-polar equality, exact quantifiers, and what is still open

Every summand on the right side of (D2) is globally sound, so that sum is contained in \(\mathscr P\). Conversely, take any \(f\in\mathscr P\). If it is nonnegative at all eight PR vertices, SG1's nonsignalling vertex argument expresses it as a nonnegative event sum. Otherwise choose a PR orientation with \(d=-f(\mathrm{PR})>0\) and its nonnegative event coordinates, and normalize by \(d+\sum_i a_i\). The compact convex canonical source base is a convex combination of at most nine of its extreme points (SG10). Extreme points at \(d=0\) are event sources; noncoercive points lie on SG12/SG13; coercive extremes have five/six-event finite contacts. The six-event ones belong to \(\mathcal K\) **by hypothesis (H)**; their five-event boundaries belong to \(\mathcal K\) by (D1) and its closedness. Scale and relabel back, proving \(f\) belongs to the right side of (D2). This proves equality.

The native/contact/density assertion is a first-order semialgebraic statement: for every normalized physical five-event source and every \(\varepsilon>0\), there exists a strict normalized native six-event source at squared coefficient distance below \(\varepsilon^2\). It follows over \(\mathbb R\) from (D6) and transfers to any real closed coefficient field. Closedness of a **fixed** finite Gram image is likewise a first-order coefficient statement (quantifying over all positive radii and PSD matrices) and transfers. The SG at-most-nine-atom reduction is already a fixed finite first-order statement. Therefore this proof introduces **no restriction on independent source/testing infinitesimal rates**; with the fixed \(Z_6\) basis all testing/state variables are represented by one original full-parent identity. The globally sound nature of each Gram certificate covers finite angles, signs, literal zero angles, infinity via homogeneous limits, and arbitrary real or complex state coordinates.

If \(\mathcal K\) is the fixed finite Gram cone, its dual is the explicit spectrahedral condition \(\mathcal A^*y\succeq0\). The duals of SG12/SG13 have their separately established finite SDP presentations (SG11–SG14); event duals are linear inequalities. Hence (D2) conditionally yields an explicit **finite hybrid intersection of dual SDP constraints**. Its conclusion is not justified without hypothesis (H), and it does not assert that SG12/SG13 themselves belong to \(\mathcal K\).

**Final boundary of this result.** There is no residual *singular extreme-contact* stratum requiring a new fixed-cone certificate after (H); SG12/SG13 and events are the only other source mechanisms, and arbitrary nonextreme sources enter by finite sums. The two *unproved construction implications* remain precisely the universal strict C1 PSD eight-gauge feasibility and strict C2 PSD dual positivity (the supplied 24-term inequality). This theorem does **not** settle either one and does not establish an unconditional global SDP lift.

## 7. Reproducibility

`verify_native_entrances.py` checks exact forward/inverse C1 and C2 native contact maps, contact-gap identities, and orthogonality of the deterministic coefficient basis, with symbolic parameters and no numerical examples. The analytic implicit-function step and coefficient-space closedness are proved above, not inferred from these identities. The source snapshot `INPUT_07_FIXED_CONE_DENSITY_sources.txt` is reproduced without alterations in the accompanying ZIP, as is the exact checker output.
