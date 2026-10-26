# Reader route to the frozen zero-computer proof that \(\operatorname{sxdeg}(\mathcal Q_{2,2,2})=2\)

**Corrected companion edition — 31 August 2026**  
**Supersedes:** `Q222_SXDEG2_READER_ROUTE_2026-08-30.md`

This document is the shortest route through the proof that I would give to a human referee. It is not a replacement for the displayed determinant and packet calculations in the full parchment. The canonical full proof is the exact 30 August parchment with the six surgical replacements recorded in `Q222_SXDEG2_CORRECTION_PATCH_AND_FREEZE_2026-08-31.md`.

## Theorem

For the complete bipartite quantum behavior set with two binary observables per party,
\[
\boxed{\operatorname{sxdeg}(\mathcal Q_{2,2,2})=2.}
\]
Equivalently, the set has a finite lift by \(2\times2\) positive-semidefinite blocks, hence by Lorentz cones.

The theorem does **not** supply a solver-ready affine Lorentz formulation and does not determine the minimum number \(M_{\rm SOC}\) of Lorentz blocks.

---

# 1. What has to be proved

Let \(K_\star\) be one canonical CHSH star in the normalized Bell polar. Let \(E\) be the compact semialgebraic family consisting of the normalized physical supports, the coordinate supports of the star simplex, and the constant support.

The full parchment proves
\[
P_{K_\star}=\operatorname{cone}(E).
\tag{R1}
\]
Take
\[
S=\operatorname{Ext}(K_\star).
\]
The set \(S\) is semialgebraic because extremality is first-order definable, and compact convexity gives
\[
K_\star=\operatorname{conv}(S).
\tag{R2}
\]
Scheiderer's Theorem 3.10 therefore reduces the upper bound to
\[
\boxed{
\operatorname{sosx}(e^{\otimes}(x))\le2
\quad
\text{for every real closed }R\supseteq\mathbb R,
\ e\in E_R,
\ x\in S_R.}
\tag{R3}
\]
The point and support copies in \(e^{\otimes}(x)\) remain independent throughout.

If the multiplication residue of the tensor is strictly positive, Scheiderer's positive-residue result gives \(\operatorname{sosx}=1\). The only problem is a zero-residue contact.

---

# 2. The universal rank-two mechanism

Let \(Y\) be a smooth one-copy contact variety, \(B=\mathbb R[Y]\), and
\[
J=\ker(B\otimes B\to B)
\]
the diagonal ideal. Put \(\delta f=f\otimes1-1\otimes f\).

## Exact diagonal lemma

Suppose
\[
\tau\in J^2
\]
and the class of \(\tau\) in \(J^2/J^3\) has positive-definite real residue. Then
\[
\boxed{
\tau=\sum_\nu\gamma_\nu\,\delta(f_\nu)^2,
\qquad
\overline{\gamma_\nu}>0.}
\tag{R4}
\]

The proof is exact, not Taylor-truncated. Choose regular coordinates \(y_i\) and a finite list \(z_\nu\) containing \(y_i\) and \(y_i\pm y_j\). Their squared divided differences generate \(J^2\), since
\[
4\delta y_i\delta y_j
=
\delta(y_i+y_j)^2-\delta(y_i-y_j)^2.
\]
Subtract a sufficiently small positive real frame while preserving positivity of the quadratic class. Represent the remaining class modulo \(J^3\), diagonalize its real residue matrix, and rewrite cross terms as divided-difference squares. Only then one has
\[
\tau=\sum_\nu q_\nu\delta(z_\nu)^2+\mu,
\qquad
\mu\in J^3,
\qquad
\overline q_\nu>0.
\]
Since the same squares generate \(J^2\),
\[
J^3=JJ^2=\sum_\nu J\,\delta(z_\nu)^2,
\]
so the whole remainder is absorbed into the coefficients without changing their positive residues. This proves (R4).

Each \(\delta f\) has tensor rank at most two. Thus every smooth contact with positive intrinsic quadratic class is finished.

### Consequence

Mixed Bell normal variables are allowed as **cotangent coordinates** for identifying the intrinsic Hessian. They are never required to be SOS roots. Once positivity is known intrinsically, (R4) manufactures new genuine divided-difference roots of rank at most two.

---

# 3. Bell geometry is reduced to one explicit slack

Fine's theorem gives the local polar as the convex hull of the positivity body \(G\) and the eight CHSH apices \(C_s\). The visible facet \(F_s\) of \(G\) is unique. Indeed, facets of \(G\) are indexed by deterministic and PR vertices of the binary no-signalling polytope. For deterministic \(q\),
\[
\langle q,C_s\rangle\le1,
\]
whereas for a PR vertex \(q_t\),
\[
\langle q_t,C_s\rangle
=
\begin{cases}
2,&t=s,\\
-2,&t=-s,\\
0,&t\ne\pm s.
\end{cases}
\]
Thus only \(q_s\) exposes a facet from which \(C_s\) lies strictly outside. Together with the pairwise midpoint identity, the visible-star lemma gives
\[
\mathcal B_Q
=G\cup\bigcup_{s=1}^{8}K_{\star,s},
\tag{R5}
\]
where \(G\) is polyhedral and all eight stars are affine Bell relabelings of one \(K_\star\).

Use the real half-angle measurement chart
\[
A_0=B_0=Z,
\]
\[
A_1(u)=\frac{(1-u^2)Z+2uX}{1+u^2},
\qquad
B_1(v)=\frac{(1-v^2)Z+2vX}{1+v^2}.
\]
Every physical support of \(K_\star\) pulls back, after a positive pure-support clearing, to the explicit nine-channel matrix
\[
\begin{aligned}
P_\star={}&
4\lambda D_ur_0r_0^T
+4e_1D_uD_vr_1r_1^T
+4e_2r_2r_2^T\\
&+4e_3D_vr_3r_3^T
+4e_4D_uD_vr_4r_4^T
+4e_5r_5r_5^T\\
&+4\alpha D_us_3s_3^T
+4\beta D_vs_2s_2^T\\
&+2d\bigl(r_5r_5^T+D_vs_2s_2^T+D_us_3s_3^T-D_uD_vzz^T\bigr).
\end{aligned}
\tag{R6}
\]
The eight event coefficients and the signed coefficient \(d\) are the only original Bell coefficients filtered in the proof.

---

# 4. Finite contacts

The eight positive event rows form a rank-four matroid. The full parchment prints one nonzero minor for every independent orbit and one determinant for every four-row orbit. The only four-row circuits are
\[
Q_0=(\lambda,e_1,e_2,\alpha),
\]
\[
Q_1=(e_1,e_2,e_4,e_5),\qquad u^2=v^2,
\]
\[
Q_2=(\lambda,e_3,\alpha,\beta),\qquad u^2v^2=1.
\tag{R7}
\]
No five positive rows have a nonzero common state kernel.

For every noncircuit finite contact the determinant remaining after the state-normal directions are removed is a \(B_2\)-Laurent polynomial. In logarithmic variables
\[
\nabla^2\Phi=\sum_\rho M_\rho\rho\rho^T,
\]
so
\[
\det\nabla^2\Phi
=\sum_{\rho<\sigma}M_\rho M_\sigma\det(\rho,\sigma)^2.
\tag{R8}
\]
The full parchment prints the only four possible collinear discriminants and the empty-mask discriminants; all are strictly nonzero. Hence every noncircuit finite zero is a smooth contact with positive complete Hessian and closes by (R4).

## The two hard circuits

For \(Q_1\) and \(Q_2\), the proof never chooses a largest active weight. Let \(w_1,\dots,w_4\) be the four active event weights. An unresolved circuit satisfies
\[
\frac d{w_i}\in\mathfrak m.
\]
Use the literal point-side scale
\[
\rho=A_{\rm raw}+B_{\rm raw}+\sum_i\frac{d^2}{w_i}.
\tag{R9}
\]
Then
\[
\rho/w_i\in B,
\qquad
\frac{d^2}{w_i\rho}\in B,
\tag{R10}
\]
for arbitrary relative orders.

Complete each normal packet separately:
\[
4w_if_i^2+2d\ell_i f_i
=
\left(2\sqrt{w_i}f_i+\frac d{2\sqrt{w_i}}\ell_i\right)^2
-
\frac{d^2}{4w_i}\ell_i^2.
\tag{R11}
\]
Every displayed root contains at most two separated tensors. The residual scalar children factor as
\[
R_{Q_1}=4(\sqrt B\,U-\sqrt A)^2+\eta U,
\tag{R12}
\]
\[
R_{Q_2}=\frac{1+U}{U}
\left[4(\sqrt A\,U-\sqrt B)^2+\eta U\right],
\tag{R13}
\]
with \(\eta\ge0\). Their interior equality zeros have positive tangent curvature. Hence \(Q_1,Q_2\) close without a valuation-rate fan.

---

# 5. The collar theorem is exact cofactor linear algebra

The three historical noncoordinate divisors have exact one-copy gauges
\[
G_F=c^2G_E\circ\phi.
\tag{R14}
\]
Let \(J_F,J_E\) be the value-plus-first-derivative \(5\times6\) contact matrices. Differentiation gives
\[
J_F=
\begin{pmatrix}
c^2&0\\
2c\nabla c&c^2D\phi
\end{pmatrix}
J_E\circ\phi.
\tag{R15}
\]
Therefore every signed maximal cofactor is multiplied by
\[
\det A=c^{10}\det D\phi.
\tag{R16}
\]
For independent point and support copies,
\[
\boxed{
\tau_F(q_p,q_s)
=\lambda_F(q_p)c(q_s)^2
\tau_E(\phi(q_p),\phi(q_s)).}
\tag{R17}
\]

For \(F_1\), the map must be read with its full definition
\[
\boxed{p:=Bs^2-u,}
\]
\[
\phi_1(s,u,a,B)=\left(u,s,\frac{as}{p},\frac1p\right),
\qquad
\det D\phi_1=\frac{s^3}{p^3}.
\]
The exact positive point multipliers on the three open physical gauges are
\[
\lambda_{F_1}=s^3p^7,
\qquad
\lambda_{F_4}=-b^7s^{18},
\qquad
\lambda_{F_3}=r^{26}.
\tag{R18}
\]
Thus every open collar is literally an EVENT contact multiplied by a pure-point factor and a pure-support square. No comparison of point and support infinitesimal rates occurs.

Finite algebraic boundary factors are routed through the strict children printed in the parchment. If a denominator vanishes because the measurement leaves the affine gauge, no quotient is taken; the contact is routed to the one-projective or simultaneous-projective homogeneous analysis.

---

# 6. Projective contacts

There are only three projective levels.

## One projective measurement

After an exact reciprocal chart, write the small measurement coordinate as \(u\in\mathfrak m\), with the other coordinate a unit. The original nine packets force
\[
\tau=\frac{u\mathcal H}{d}\in B^\times,
\qquad
\delta=\frac d{\mathcal H}\in\mathfrak m,
\tag{R19}
\]
where \(\mathcal H\) is the positive harmonic split scale. The four forcing coefficients satisfy
\[
\lambda,e_2,e_5,\alpha\in\frac{d^2}{\mathcal H}B.
\tag{R20}
\]
This is derived from the first radial nonnegative envelope, whose residue derivative has simple Jacobian \(8(1+\bar v^2)>0\). Thus there is an exact point-side critical center.

The normal split is exhaustive: full-normal, one-missing, or both-missing. The first two have explicit positive normal residue matrices; both-missing is impossible. The remaining scalar determinant has a finite nine-entry coefficient table. Its only nonregular stationary support is
\[
\frac\Phi4
=(H_0\tau\theta-Qv)^2+(H_0\tau v\theta-A)^2,
\tag{R21}
\]
with nonzero two-root Jacobian.

After that, exactly three possibilities remain for the delayed split pair:

1. both current, where the missing scalar curvature becomes positive;
2. neither current, where the first scalar class is \(A\tau^2+B\tau^{-2}+C\) and every interior zero has curvature \(8A>0\);
3. exactly one current, the A17 terminal.

No mixed projective coordinate is used as an SOS root.

## A17

Two exact unit-Jacobian divisions derive, rather than assume, the terminal orders. The centered complete tensor satisfies
\[
A_s^4\widetilde\Theta_\beta\in J^2,
\qquad
A_s^2\widetilde\Theta_{e_1}\in J^2,
\tag{R22}
\]
for
\[
J=(C_Xx_0,h_1,h_2,dx_3,drx_4).
\]
The full parchment prints the finite below-cost quotient bases and cancellation triples. The four-dimensional residue block is positive definite, and the **centered** fifth Schur pivots are
\[
\boxed{
q_\beta=\frac4{T^3(1+V^2)}>0,
\qquad
q_{e_1}=\frac4{T^3V^2(1+V^2)}>0.}
\tag{R23}
\]
There are no \(K_3\)-dependent additions in these pivots. The positive infinitesimal critical value \(m=u^3K_3\ge0\) remains in the proof and is restored separately as a pure-side square after centering.

A five-simplex of genuine divided differences handles the fifth direction. The remaining four-dimensional part closes by the diagonal lemma.

---

# 7. Simultaneous projective contacts

The five homogeneous positive EVENT squares have exactly two real base components:
\[
B_0:\quad u_1=v_1=0,\ [y]=[0:1:0],
\]
\[
B_\infty:\quad v_0=0,\ y_3=0.
\tag{R24}
\]
There is no third real projective component.

## \(B_0\)

The exact leading cofactor signs derive the strict chamber. The special factor \(LT^2=1\) is treated before division: one contact-matrix column vanishes, the remaining five columns have rank four, and their only relation has mixed signs. Hence the only nonnegative kernel ray is the exact \(e_1\) coordinate child.

On the strict chamber, after one pure-point congruence the complete \(4\times4\) physical matrix is literally
\[
\begin{aligned}
S_0={}&
\frac{AT(1+T^2)p}{w}(1+U^2)e_0e_0^T
+A\varepsilon_0(1+V^2)\omega\omega^T\\
&+\varepsilon_0p(1+U^2)(1+V^2)e_1e_1^T
+\frac{A\varepsilon_0(1+T^2)wp}{h}U^2e_2e_2^T.
\end{aligned}
\tag{R25}
\]
The scalar coefficients are nonnegative and **separate factorwise**: \(A,T,p,w,h,\varepsilon_0\) are point-side, while the only support dependence is through \(1+U^2\), \(1+V^2\), and \(U^2\). Absorbing these pure-support factors into the displayed support vectors gives literal tensor-square roots of separation rank at most two. The only remaining boundary \(h=0\) is closed by its next positive tiers.

## \(B_\infty\)

Use the reciprocal point chart
\[
s=cr,
\quad u=r^{-1},
\quad a=((J+c^2)r^2)^{-1},
\quad B=\frac{1+Kr^2}{(J+c^2)r^3}.
\tag{R26}
\]
The exact point-boundary factorizations prove that every unresolved row lies in
\[
\boxed{c>0,\qquad J>1,\qquad0\le K<J.}
\tag{R27}
\]
Here \(K=0\) is retained. Exact \(J=1\), \(K=J\), and \(c=0\) are routed before reciprocal division.

A fraction-free Cramer identity in the complete value/first-jet system is
\[
D_0D_1d=2J\zeta N_d\lambda,
\tag{R28}
\]
with
\[
\bar D_0=\bar D_1=\bar N_d=\bar J>0.
\]
Hence the signed packet cannot occur strictly before \(\lambda\). On the genuinely coupled \(\lambda\)-row one gets
\[
d=2z,
\qquad
 e_3=\frac zP,
\qquad
 e_5=\frac{Jz(1+z)}R,
\tag{R29}
\]
with the unit-or-child firewall
\[
\boxed{\bar z,\bar P,\bar R,\bar H>0.}
\tag{R30}
\]
Only after (R30) is established is the conormal Hessian used.

The hard \(2\times2\) conormal determinant factors as
\[
\det Q(j)
=
\frac{16z(1+z)(j-J)^2}{(1-Jz)(Jz+J-1)},
\tag{R31}
\]
and the complete even radial restriction is
\[
F_{\ker}(J,\rho)
=
4Jz(1+z)(J+\gamma^2)\rho^2
+\rho^4R_1(\rho^2).
\tag{R32}
\]
Thus the intrinsic conormal quadratic class is
\[
\frac12\nabla^2F_{\rm con}(0)
=
\operatorname{diag}(A,\alpha_j,\alpha_\rho)
\succ0
\tag{R33}
\]
with **positive residues**, not merely positive valued-field elements. The mixed triangular variables are only cotangent coordinates. The diagonal lemma supplies genuine rank-two divided-difference roots.

---

# 8. Why the descent terminates

Gauge maps and blow-ups are coordinate changes, not descent arrows.

A genuine child either sets at least one of the nine original nonnegative Bell coefficients to zero, or, at fixed active coefficient set, imposes an additional exact coordinate/projective/attachment equation and enters a proper lower-dimensional stratum.

Use
\[
\boxed{\kappa=(N_{\rm act},d_{\rm stratum})}
\tag{R34}
\]
with lexicographic order. Every genuine child strictly decreases \(\kappa\). Since \(0\le N_{\rm act}\le9\) and the one-copy dimensions are finite, no infinite chain exists.

A separate common-zero incidence lemma ensures that the tested support and a point-defining active support can be placed in one common strict contact node or one common proper child. At a zero affine-gauge denominator, the proof routes to the appropriate finite child or homogeneous projective analysis before division. No statement of the form
\[
\text{different specializations}\Rightarrow\text{positive residue}
\]
is used.

Combining the finite, collar, one-projective, A17, \(B_0\), and \(B_\infty\) modules with induction on (R34) proves (R3). Therefore
\[
\operatorname{sxdeg}(K_\star)\le2.
\tag{R35}
\]

---

# 9. Global upper bound

All eight stars are affine Bell relabelings of \(K_\star\), and \(G\) is polyhedral. Scheiderer's finite convex-hull closure gives
\[
\operatorname{sxdeg}(\mathcal B_Q)\le2.
\]
Polarity preserves this upper bound, hence
\[
\boxed{\operatorname{sxdeg}(\mathcal Q_{2,2,2})\le2.}
\tag{R36}
\]

---

# 10. Independent lower bound

In the correlator projection impose
\[
C_{00}=C_{01}=C_{10}=t=\cos\alpha,
\qquad
0\le\alpha\le\pi/4.
\]
Tsirelson's vector representation gives real unit vectors with these scalar products. The spherical triangle inequality yields
\[
d(a_1,b_1)\le3\alpha,
\]
so
\[
C_{11}\ge\cos3\alpha.
\]
Equality is attained by coplanar vectors at angles
\[
0,\ \alpha,\ -\alpha,\ 2\alpha.
\]
Thus the exact boundary contains
\[
\boxed{
C_{11}=4t^3-3t,
\qquad
1/\sqrt2\le t\le1.}
\tag{R37}
\]
This is a curved cubic arc. Hence the quantum set is not polyhedral. Since semidefinite extension degree one is exactly the polyhedral case,
\[
\operatorname{sxdeg}(\mathcal Q_{2,2,2})\ge2.
\tag{R38}
\]
Together with (R36),
\[
\boxed{\operatorname{sxdeg}(\mathcal Q_{2,2,2})=2.}
\]

---

# 11. Hostile-referee verification list

The Bell-specific verification is finite:

1. visible-facet uniqueness, pairwise-midpoint identity, and star decomposition;
2. support-cone separation;
3. exact diagonal lemma and \(J^3\) absorption;
4. event-row minors and the three circuit loci;
5. finite collinear and empty-mask discriminants;
6. \(Q_1,Q_2\) syzygies, scale \(\rho\), packet completion, and Laurent factorizations;
7. the exact gauge/cofactor identity, including \(p=Bs^2-u\) for \(F_1\);
8. terminal \(F_3\) kernel coefficient and genuine one-copy cotangent change;
9. one-projective radial envelope and two-square singular child;
10. A17 unit-Jacobian divisions, cancellation basis, residue block, **K3-free centered pivots**, and five-simplex divided-difference completion;
11. two-component projective base locus;
12. \(B_0\) chamber, factorwise-separated four-square suspension, and \(h=0\) closure;
13. \(B_\infty\) boundary factorizations, Cramer identity, unit firewall, determinant factorization, and radial coefficient;
14. common-zero incidence and lexicographic descent;
15. global Scheiderer assembly and cubic lower bound.

No item asks the reader to run a solver, symbolic algebra system, Newton scan, Pareto search, or production atlas.

---

# 12. Frontier after the theorem

The theorem-level existence statement is separate from
\[
\boxed{\text{printing an explicit finite affine }L_3\text{ lift}}
\]
and from determining the exact minimum \(M_{\rm SOC}\). Neither is proved by the extension-degree theorem itself.
