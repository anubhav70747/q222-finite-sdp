# A consolidated parchment proof that the complete \((2,2,2)\) quantum behavior set has semidefinite extension degree two

**Final self-contained reader version — 30 August 2026**


## Theorem

Let
\[
\mathcal Q_{2,2,2}\subset\mathbb R^8
\]
be the complete bipartite quantum behavior set with two binary observables per party, in marginal-correlator coordinates
\[
(A_0,A_1,B_0,B_1,C_{00},C_{01},C_{10},C_{11}).
\]
Then
\[
\boxed{\operatorname{sxdeg}(\mathcal Q_{2,2,2})=2.}
\]
Equivalently, \(\mathcal Q_{2,2,2}\) admits a finite semidefinite lift using only \(2\times2\) positive-semidefinite blocks, or equivalently Lorentz cones.

This is an existence theorem. It does **not** print a solver-ready finite affine Lorentz system and does not determine the minimum number \(M_{\rm SOC}\) of Lorentz blocks.

### Logical scope

No numerical SDP, computer algebra output, Newton or Pareto enumeration, stored Gram matrix, machine-generated contact atlas, or fixed production-cell table is a premise. The lengthy identities that remain are finite displayed polynomial identities and determinant calculations. They may be checked by hand expansion.

The proof uses three specialized external results only:

1. Fine's characterization of the local \((2,2,2)\) polytope.
2. The standard two-dichotomic-observable qubit reduction, which is also reconstructed below.
3. Scheiderer's semidefinite-extension-degree criterion, positive-residue tensor lemma, convex-hull/polar closure results, and the diagonal-ideal argument underlying his 2025 Lemma 3.12.

---

# 1. Global funnel: from the quantum set to one canonical Bell star

## 1.1 Real two-qubit reduction

For outcomes \(a,b\in\{\pm1\}\),
\[
p(a,b|x,y)=\frac14\bigl(1+aA_x+bB_y+abC_{xy}\bigr).
\]
A binary observable is a Hermitian contraction. For any such contraction \(A\), define the dilation
\[
\widehat A=
\begin{pmatrix}
A&\sqrt{I-A^2}\\
\sqrt{I-A^2}&-A
\end{pmatrix}.
\]
Because \(A\) commutes with \(\sqrt{I-A^2}\), one has \(\widehat A^2=I\), and compression to the first summand recovers \(A\). Thus both observables on each side may be taken projective.

Jordan's two-projector lemma decomposes each local Hilbert space into invariant blocks of dimension at most two. Hence every behavior is a convex combination of behaviors from local blocks of dimensions at most two. On a two-dimensional block, local unitaries place both observables in a common real \(XZ\)-plane. All measured operators are real. Replacing a state \(\rho\) by
\[
\rho_{\mathbb R}=\frac12(\rho+\overline\rho)
\]
preserves all measured expectations, and a real spectral decomposition expresses \(\rho_{\mathbb R}\) as a convex combination of pure real states. Therefore
\[
\mathcal Q_{2,2,2}=\operatorname{conv}\Phi(\mathscr P),
\]
where \(\mathscr P\) is a compact semialgebraic union of real two-qubit projective charts. In particular \(\mathcal Q_{2,2,2}\) is compact, convex, and semialgebraic.

## 1.2 Polar and eight-star decomposition

Let
\[
\mathcal B_Q=\mathcal Q_{2,2,2}^{\circ}.
\]
For each setting pair and outcome pair define a positivity normal \(E_{xy}^{ab}\) by
\[
(E_{xy}^{ab})_{A_x}=-a,\qquad
(E_{xy}^{ab})_{B_y}=-b,\qquad
(E_{xy}^{ab})_{C_{xy}}=-ab.
\]
Then
\[
1-\langle E_{xy}^{ab},q\rangle=4p(a,b|x,y)\ge0.
\]
Put
\[
G=\operatorname{conv}\{E_{xy}^{ab}\}.
\]
For an odd CHSH sign pattern \(s\), put
\[
C_s=\frac12(0,0,0,0,s_{00},s_{01},s_{10},s_{11}).
\]
Fine's theorem gives
\[
\mathcal L^{\circ}=\operatorname{conv}\bigl(G,\{C_s\}\bigr),
\]
where \(\mathcal L\) is the local behavior polytope.

The PR behavior associated with \(s\) exposes exactly the eight positivity normals with \(ab=-s_{xy}\). Their convex hull is a facet \(F_s\) of \(G\), and
\[
\frac18\sum_{E\in\operatorname{vert}(F_s)}E=\frac12C_s.
\]
For distinct odd patterns \(s,t\),
\[
\frac{C_s+C_t}{2}\in G.
\]
We use the following elementary convexity lemma.

**Visible-star lemma.** Suppose each exterior point \(c_i\) sees exactly one facet \(F_i\) of a compact convex set \(G\), and \((c_i+c_j)/2\in G\) for \(i\ne j\). Then
\[
\operatorname{conv}(G,c_1,\ldots,c_N)
=G\cup\bigcup_i\operatorname{conv}(F_i,c_i).
\]
Indeed, in a convex combination of the \(c_i\), either every coefficient is at most \(1/2\), in which case the point is a convex combination of pairwise midpoints, or one coefficient is larger than \(1/2\), in which case the point lies in the corresponding star.

Since
\[
\mathcal L\subseteq\mathcal Q_{2,2,2}\subseteq\mathcal{NS},
\]
polarity gives
\[
G\subseteq\mathcal B_Q\subseteq\mathcal L^{\circ}.
\]
Therefore
\[
\boxed{
\mathcal B_Q=G\cup\bigcup_sK_{\star,s},
\qquad
K_{\star,s}=\mathcal B_Q\cap\operatorname{conv}(F_s,C_s).
}
\]
All eight stars are affine Bell-relabeling images of one canonical star \(K_\star\). It is enough to prove
\[
\operatorname{sxdeg}(K_\star)\le2.
\]

---

# 2. The support cone and Scheiderer's quantitative target

Fix one canonical apex \(C\) and put
\[
x_*=\frac23C.
\]
The average of the visible facet is \(C/2\), so \(x_*\) is strictly inside the star simplex. The CHSH operator
\[
\mathcal C=A_0(B_0-B_1)+A_1(B_0+B_1)
\]
satisfies
\[
\mathcal C^2=4I+[A_0,A_1]\otimes[B_0,B_1],
\]
hence \(\|\mathcal C\|\le2\sqrt2\). It follows that every physical support has a strictly positive value at \(x_*\):
\[
1-\langle x_*,p\rangle\ge1-\frac{2\sqrt2}{3}>0.
\]
Normalize each physical support by this value. Let \(E\) consist of

1. all normalized supports coming from pure real two-qubit projective behaviors;
2. the nine normalized coordinate supports of the star simplex;
3. the constant support \(1\).

Then \(E\) is compact semialgebraic and
\[
K_\star=\{x:e(x)\ge0\ \forall e\in E\}.
\]

## 2.1 Exact support-cone theorem

Let \(P_{K_\star}\) be the cone of affine functions nonnegative on \(K_\star\). Then
\[
\boxed{P_{K_\star}=\operatorname{cone}(E).}
\]

**Proof.** Every normalized generator satisfies \(e(x_*)=1\), so \(\operatorname{cone}(E)\) is closed. Suppose \(f\in P_{K_\star}\setminus\operatorname{cone}(E)\). Strict separation gives a linear functional
\[
L(g)=tg_0+\langle y,g_1\rangle
\]
that is nonnegative on \(E\) and negative on \(f\). The constant support gives \(t\ge0\).

If \(t>0\), then \(y/t\) satisfies every inequality in \(E\), hence lies in \(K_\star\), contradicting \(L(f)<0\). If \(t=0\), then \(x_*+sy\in K_\star\) for all \(s\ge0\), contradicting compactness unless \(y=0\). This also contradicts strict separation. \(\square\)

Conic Carathéodory gives at most nine generators. Since the statement is first-order semialgebraic, it transfers to every real closed extension.

## 2.2 Quantitative criterion

Scheiderer's quantitative criterion says that, for a closed convex semialgebraic set \(K\), semialgebraic \(S\subseteq K\), and semialgebraic \(E\subseteq P_K\) satisfying
\[
K=\operatorname{conv}(S),\qquad P_K=\operatorname{cone}(E),
\]
one has \(\operatorname{sxdeg}(K)\le d\) if and only if, over every real closed extension \(R/\mathbb R\),
\[
\operatorname{sosx}(e^{\otimes}(x))\le d
\]
for all \(e\in E_R\) and \(x\in S_R\).

We take
\[
S=\operatorname{Ext}(K_\star).
\]
This is a semialgebraic subset of \(K_\star\): extremality is the first-order condition that whenever \(x=(y+z)/2\) with \(y,z\in K_\star\), one has \(y=z=x\), and quantifier elimination preserves semialgebraicity. Since \(K_\star\) is compact and finite dimensional, it is the convex hull of its extreme points. Thus Scheiderer's condition (iii) may be tested only at extreme points. Evaluations with strictly positive multiplication residue have \(\operatorname{sosx}=1\) by Scheiderer's positive-residue proposition, so only zero-residue contact evaluations require work.

Our target is therefore
\[
\boxed{
\operatorname{sosx}(e^{\otimes}(x))\le2
\quad
\text{for every real closed }R/\mathbb R,
\ e\in E_R,
\ x\in\operatorname{Ext}(K_\star)_R.
}
\tag{2.1}
\]

---

# 3. The exact diagonal rank-two lemma

Let \(Y\) be a smooth one-copy variety with coordinate ring \(B\), and let
\[
J=\ker(B\otimes_{\mathbb R}B\to B)
\]
be the diagonal ideal. Write
\[
\delta f=f\otimes1-1\otimes f.
\]

## Lemma 3.1: exact diagonal compression

Suppose \(\tau\in J^2\), and its class in \(J^2/J^3\) has strictly positive-definite residue at a real point \(u\). Then, after localizing at functions with positive residue,
\[
\boxed{
\tau=\sum_\nu\gamma_\nu\,\delta(f_\nu)^2,
\qquad
\overline{\gamma_\nu}>0.
}
\tag{3.1}
\]
Consequently every evaluated square root has tensor rank at most two, and
\[
\operatorname{sosx}(\tau)\le2.
\]

**Proof.** Choose regular local coordinates \(y_1,\ldots,y_m\) whose divided differences generate \(J\). Let \(z_\nu\) be the finite list consisting of the \(y_i\) and the sums and differences \(y_i\pm y_j\). Their squared differences generate \(J^2\), because
\[
4\,\delta y_i\,\delta y_j
=
\delta(y_i+y_j)^2-\delta(y_i-y_j)^2.
\tag{3.2}
\]

Choose \(\varepsilon>0\) so small that
\[
\tau'=\tau-\varepsilon\sum_\nu\delta(z_\nu)^2
\]
still has positive-definite class in \(J^2/J^3\). Write
\[
\tau'\equiv\sum_{i,j}p_{ij}\,\delta y_i\,\delta y_j
\pmod{J^3}.
\]
After a constant real linear change of the \(y_i\), the residue matrix \((\overline{p_{ij}})\) is diagonal and nonnegative. Using (3.2), rewrite the remaining off-diagonal terms as squares of the \(\delta(z_\nu)\). Thus
\[
\tau=\sum_\nu q_\nu\,\delta(z_\nu)^2+\mu,
\qquad
\mu\in J^3,
\qquad
\overline{q_\nu}\ge\varepsilon>0.
\tag{3.3}
\]
Since the \(\delta(z_\nu)^2\) generate \(J^2\),
\[
J^3=J\,J^2=\sum_\nu J\,\delta(z_\nu)^2.
\]
Hence
\[
\mu=\sum_\nu h_\nu\,\delta(z_\nu)^2,
\qquad h_\nu\in J.
\]
Substitution into (3.3) gives (3.1), because \(h_\nu\) has zero residue and therefore \(\overline{q_\nu+h_\nu}>0\).

Scheiderer's positive-residue proposition writes each \(\gamma_\nu\) as an SOS of rank-one elementary tensors. Finally
\[
\delta f=f\otimes1-1\otimes f
\]
has tensor rank at most two. \(\square\)

### Corollary 3.2: mixed cotangent coordinates are harmless

To check the hypothesis of Lemma 3.1 it is enough to exhibit **any** basis \(\eta_1,\ldots,\eta_m\) of \(J/J^2\) for which
\[
[\tau]=\eta^{\mathsf T}H\eta,
\qquad \overline H\succ0.
\]
The \(\eta_i\) themselves need not be legal rank-two roots. They are used only to identify the intrinsic quadratic class. Lemma 3.1 then manufactures new genuine divided-difference roots.

This observation is the reason historical weighted Bell normal variables may be used for Hessian calculations without ever being promoted to SOS roots.

### Lemma 3.3: pure-side invariance

Multiplication or division by a nonzero factor lying wholly in the point factor or wholly in the support factor does not increase tensor rank. Thus if
\[
\widetilde\tau=(a_{\rm p}b_{\rm s})^2\tau
\]
and \(a_{\rm p},b_{\rm s}\ne0\), a rank-two SOS for \(\widetilde\tau\) gives one for \(\tau\) by dividing each root by \(a_{\rm p}b_{\rm s}\).

### Lemma 3.4: critical-value splitting

If an exact one-copy support has minimum value \(m\ge0\) on its exact critical graph, subtract the support-side constant \(m\). The centered loss has zero value and zero first jet. The value is restored by the pure-side square \((\sqrt m)^2\). Positive infinitesimal critical values are therefore never silently set to zero.

---

# 4. The canonical Bell slack

Put
\[
D_u=1+u^2,\qquad D_v=1+v^2,
\]
and define the eight positive EVENT row vectors
\[
\begin{array}{c|rrrr}
\lambda&1&v&0&0\\
e_1&0&1&0&0\\
e_2&v&-1&uv&-u\\
e_3&u&0&-1&0\\
e_4&0&0&1&0\\
e_5&u&uv&-1&-v\\
\alpha&0&0&v&-1\\
\beta&0&1&0&u.
\end{array}
\tag{4.1}
\]
Equivalently, with
\[
\begin{aligned}
r_0&=(1,v,0,0)^T,&r_1&=(0,1,0,0)^T,\\
r_2&=(v,-1,uv,-u)^T,&r_3&=(u,0,-1,0)^T,\\
r_4&=(0,0,1,0)^T,&r_5&=(u,uv,-1,-v)^T,\\
s_2&=(0,1,0,u)^T,&s_3&=(0,0,v,-1)^T,\\
z&=(0,0,0,1)^T,
\end{aligned}
\]
the denominator-cleared canonical slack is
\[
\boxed{
\begin{aligned}
P_\star={}&4\lambda D_ur_0r_0^T
+4e_1D_uD_vr_1r_1^T+4e_2r_2r_2^T
+4e_3D_vr_3r_3^T\\
&+4e_4D_uD_vr_4r_4^T+4e_5r_5r_5^T
+4\alpha D_us_3s_3^T+4\beta D_vs_2s_2^T\\
&+2d\bigl(r_5r_5^T+D_vs_2s_2^T+D_us_3s_3^T-D_uD_vzz^T\bigr).
\end{aligned}
}
\tag{4.2}
\]
The coefficients
\[
d,\lambda,e_1,e_2,e_3,e_4,e_5,\alpha,\beta
\]
are nonnegative star coordinates. The first eight packets are literal positive rank-one EVENT squares. The last line is the signed CHSH packet.

For a real realization \((\psi,u,v)\), direct substitution gives
\[
D_uD_v\|\psi\|^2\bigl(1-\langle x,p_{\psi,u,v}\rangle\bigr)
=\psi^TP_\star(x;u,v)\psi.
\tag{4.3}
\]
Thus it suffices to prove the rank-two tensor theorem for the canonical slack, after positive pure-support clearing.

---

# 5. Finite EVENT contacts

The finite classification is short because the eight EVENT rows form a tiny rank-four matroid.

## 5.1 The event matroid

Every set of at most three rows is independent on \(uv\ne0\). Up to the cyclic Bell symmetry, representatives have nonzero minors
\[
\begin{array}{c|c}
(\lambda)&1\\
(\lambda,e_1)&1\\
(\lambda,e_2)&-(1+v^2)\\
(\lambda,e_3)&-uv\\
(\lambda,\alpha)&v\\
(\lambda,e_1,e_2)&uv\\
(\lambda,e_1,e_4)&1\\
(\lambda,e_1,e_5)&-1\\
(\lambda,e_2,e_3)&1+v^2+u^2v^2\\
(\lambda,e_3,\alpha)&-uv^2.
\end{array}
\tag{5.1}
\]
For four rows, the eight bracelet representatives have determinants
\[
\begin{array}{c|c}
(\lambda,e_1,e_2,e_3)&-u\\
(\lambda,e_1,e_2,e_4)&u\\
(\lambda,e_1,e_2,e_5)&-u(1+v^2)\\
(\lambda,e_1,e_2,\alpha)&0\\
(\lambda,e_1,e_2,\beta)&u^2v\\
(\lambda,e_1,e_4,\alpha)&-1\\
(\lambda,e_1,e_4,\beta)&u\\
(\lambda,e_3,\alpha,\beta)&-(uv-1)(uv+1).
\end{array}
\tag{5.2}
\]
The alternating orbit also contains
\[
\det(e_1,e_2,e_4,e_5)=-(u-v)(u+v).
\tag{5.3}
\]
Hence the only four-row circuit loci, up to Bell relabeling, are
\[
\boxed{Q_0=(\lambda,e_1,e_2,\alpha),}
\]
\[
\boxed{Q_1=(e_1,e_2,e_4,e_5),\quad u^2=v^2,}
\]
\[
\boxed{Q_2=(\lambda,e_3,\alpha,\beta),\quad u^2v^2=1.}
\tag{5.4}
\]
No five-row set has a nonzero common kernel. The five dihedral orbits have the following nonzero four-minors:
\[
\begin{array}{c|c|c}
\text{five-set}&\text{four-subset}&\det\\\hline
\lambda,e_1,\beta,e_2,\alpha&\lambda,e_1,\beta,e_2&-u^2v\\
\lambda,e_1,\beta,e_2,e_4&\lambda,e_1,\beta,e_2&-u^2v\\
\lambda,e_1,\beta,\alpha,e_4&\lambda,e_1,\beta,\alpha&-uv\\
\lambda,e_1,\beta,\alpha,e_3&\lambda,e_1,\beta,\alpha&-uv\\
\lambda,e_1,e_2,\alpha,e_3&\lambda,e_1,e_2,e_3&-u.
\end{array}
\tag{5.5}
\]

## 5.2 Finite coefficient filtration

Filter only the nine original nonnegative Bell coefficients. At each step normalize by the sum of the remaining coefficients. Before the first layer containing \(d\), the normalized initial matrix is a positive sum
\[
G=\sum_i\overline w_i r_ir_i^T\succeq0.
\]
At a zero-residue scalar contact, \(\psi^TG\psi=0\), so every strictly earlier positive row annihilates the residue state. Hence their common kernel is one of the cases in (5.1)-(5.5). When the first layer containing \(d\) is reached, retain the complete signed packet and every EVENT row tied with it. No tied positive row is set to zero.

Thus a hard finite contact is an independent mask of size at most three, one of \(Q_0,Q_1,Q_2\), the empty-earlier-mask case, or impossible.

### Hand-check table for the ten independent finite masks

For the ten independent masks in (5.1), one may choose state-normal coordinates by a nonzero pure-support minor. The determinant of that normal change and the determinant of the remaining positive state core are as follows:

\[
\begin{array}{c|c|c|c}
\text{type}&\text{mask}&\text{normal determinant}&\text{positive-core determinant}\\ \hline
S1&\lambda&1+v^2&(1+v^2)^2\\
P1&\lambda,e_1&1&(1+v^2)(1+u^2+v^2)\\
P2&\lambda,e_2&(1+u^2)(1+v^2)^2&(1+v^2)^4\{u^2v^2+(1+v^2)^2\}\\
P3&\lambda,e_3&u^2v^2+v^2+1&(1+v^2)(u^4v^2+u^2v^4+u^2v^2+1)\\
P4&\lambda,\alpha&(1+v^2)^2&(1+v^2)^2\\
T1&\lambda,e_1,e_2&u^2(1+v^2)&u^2v^2+(1+v^2)^2\\
T2&\lambda,e_1,e_4&1&1+u^2+v^2\\
T3&\lambda,e_1,e_5&1+v^2&u^2+(1+v^2)^2\\
T4&\lambda,e_2,e_3&(1+u^2)(1+v^2)(u^2v^2+v^2+1)&\begin{aligned}
&u^6v^4+u^4v^6+4u^4v^4+u^4v^2\\[-1mm]
&+2u^2v^6+5u^2v^4+2u^2v^2+v^6+3v^4+3v^2+1
\end{aligned}\\
T5&\lambda,e_3,\alpha&(1+v^2)(u^2v^2+1)&\{uv(1+v^2)\}^2+(u^2v^2-1)^2
\end{array}
\tag{5.6a}
\]
Every entry is strictly positive on the finite chart \(uv\ne0\). Thus every surviving noncircuit zero has a one-dimensional state kernel, and all state-normal changes invert only pure-support units.

## 5.3 Regular finite contacts

For an independent mask, eliminate the state-normal directions. The remaining determinant is a \(B_2\) Laurent polynomial. In logarithmic coordinates its noncentral exponents lie in
\[
\Phi_{B_2}=\{(\pm1,\pm1),(\pm1,0),(0,\pm1)\},
\]
with nonnegative outer coefficients. At a stationary point
\[
\nabla^2\Phi=\sum_\rho M_\rho\rho\rho^T,
\]
so
\[
\det\nabla^2\Phi
=\sum_{\rho<\sigma}M_\rho M_\sigma\det(\rho,\sigma)^2.
\tag{5.6}
\]
Thus singularity can occur only if all active roots are collinear. Up to symmetry the only opposite-root possibilities have discriminants
\[
4d^2\mathcal R^2,
\quad
16d^2(2\alpha+2\beta+d)^2,
\quad
16d^2(2\beta+d+2e_3)^2,
\quad
64d^2,
\tag{5.7}
\]
where
\[
\mathcal R=
4\alpha\beta+2\alpha d+4\alpha e_3+2\beta d+4\beta e_3+d^2+2de_3>0.
\]
If every noncentral coefficient vanishes, the possible central coefficients are exactly
\[
\begin{gathered}
-2d^3,
-2d^2(2\beta+d),
-2d^2(d+2e_4),
-2d(d+2e_4)(2\beta+d),\\
-2d^2(d+2e_3),
-2d^2(d+2e_2),
-2d(d+2e_2)(d+2e_4),
-2d(d+2e_2)(d+2e_3),\\
-4d^2,
-4d(d+2e_4),
-4d(d+2e_3),
-4d^2,
-4d(d+2e_2),
-8d .
\end{gathered}
\tag{5.7a}
\]
They are all strictly negative for \(d>0\), so rank-zero contact is impossible.

In the empty-mask case the two possible opposite-pair discriminants factor as
\[
4c_{00}c_{22}-\Gamma_+^2=-d^2Q_+^2<0,
\qquad
4c_{02}c_{20}-\Gamma_-^2=-d^2Q_-^2<0,
\tag{5.8}
\]
where, explicitly,
\[
\begin{aligned}
Q_+={}&4\alpha\beta d+8\alpha\beta e_3+8\alpha\beta\lambda
+2\alpha d^2+4\alpha de_3+4\alpha d\lambda+8\alpha e_3\lambda\\
&+2\beta d^2+4\beta de_3+4\beta d\lambda+8\beta e_3\lambda
+d^3+2d^2e_3+2d^2\lambda+4de_3\lambda,
\end{aligned}
\tag{5.8a}
\]
\[
\begin{aligned}
Q_-={}&d^3+2d^2(e_1+e_2+e_4+e_5)\\
&+4d(e_1e_2+e_1e_4+e_1e_5+e_2e_4+e_2e_5+e_4e_5)\\
&+8(e_1e_2e_4+e_1e_2e_5+e_1e_4e_5+e_2e_4e_5).
\end{aligned}
\tag{5.8b}
\]
Every displayed monomial has positive coefficient, so \(Q_\pm>0\) whenever \(d>0\).

Hence every noncircuit finite zero contact has a smooth exact contact graph with positive complete Hessian. If
\[
\tau(q_1,q_2)=e_{q_2}(X(q_1)),
\]
then \(\tau(q,q)=0\) and both first derivatives vanish along the diagonal, so \(\tau\in J^2\). Its intrinsic quadratic class is positive definite. Lemma 3.1 gives
\[
\boxed{\operatorname{sosx}(\mathrm{REG})\le2.}
\tag{5.9}
\]

---

# 6. The three finite circuits

## 6.1 \(Q_0\)

The exact syzygy
\[
r_2=vr_0-(1+v^2)r_1+us_3
\tag{6.1}
\]
shows directly that the signed packet is positive on the circuit kernel. Thus \(Q_0\) closes as a deterministic/coordinate child.

## 6.2 Exact coordinates for \(Q_1\)

Put \(u=t\), \(v=t+q\) and
\[
\begin{aligned}
h_1&=X,\\
h_2&=vA-X+uvY-u,\\
h_3&=Y,\\
h_4&=uA+uvX-Y-v.
\end{aligned}
\]
The Jacobian with respect to \((A,X,Y,q)\) is \(2t\ne0\), and
\[
u(1+v^2)h_1+uh_2-v(1+u^2)h_3-vh_4=v^2-u^2.
\tag{6.2}
\]
Hence
\[
q=
\frac{t(1+v^2)h_1+th_2-v(1+t^2)h_3-vh_4}{t+v},
\tag{6.3}
\]
with a pure-support unit denominator.

## 6.3 Exact coordinates for \(Q_2\)

Put \(u=t\), \(v=(1+q)/t\) and
\[
g_1=A+vX,\quad
g_2=uA-Y,\quad
g_3=vY-1,\quad
g_4=X+u.
\]
The Jacobian is \(-2\), and
\[
uvg_1-vg_2-g_3-uv^2g_4=1-u^2v^2.
\tag{6.4}
\]
Again the circuit defect is an exact local coordinate.

## 6.4 Arbitrary active-weight orders

Let \(w_1,\ldots,w_4\) be the four active positive EVENT weights on an unresolved \(Q_1\) or \(Q_2\) circuit. Positivity of the current matrix and the circuit kernel imply
\[
\frac d{w_i}\in\mathfrak m
\quad\text{for every }i.
\tag{6.5}
\]
Indeed the active positive Gram annihilates the circuit kernel, while the signed packet has nonzero action on that kernel; a current residue of \(d/w_i\) would contradict positive semidefiniteness.

Define
\[
(A_{\rm raw},B_{\rm raw})=
\begin{cases}
(\lambda+\alpha,e_3+\beta),&Q_1,\\
(e_2+e_5,e_1+e_4),&Q_2,
\end{cases}
\]
and the literal positive point-side scale
\[
\boxed{
\rho=A_{\rm raw}+B_{\rm raw}+\sum_{j=1}^4\frac{d^2}{w_j}.}
\tag{6.6}
\]
A complementary positive channel with unbounded \(A_{\rm raw}/w_i\) or \(B_{\rm raw}/w_i\) would close at an earlier positive layer. Hence these ratios are bounded. Also
\[
\frac{d^2/w_j}{w_i}=\frac d{w_i}\frac d{w_j}\in B.
\]
Consequently
\[
\boxed{
\rho/w_i\in B,
\qquad
\frac{d^2}{w_i\rho}\in B.}
\tag{6.7}
\]
These are the only scale bounds needed; no largest weight or discrete valuation class is selected.

Complete each active normal separately:
\[
4w_if_i^2+2d\ell_i f_i
=
\left(2\sqrt{w_i}f_i+\frac{d}{2\sqrt{w_i}}\ell_i\right)^2
-
\frac{d^2}{4w_i}\ell_i^2.
\tag{6.8}
\]
All coefficients here are point-side and \(f_i,\ell_i\) are one-copy support functions, so each completed square root contains at most two separated tensors.

For \(Q_1\) put
\[
U=u^2,
\qquad
A=\frac{\lambda+\alpha}{\rho},
\qquad
B=\frac{e_3+\beta}{\rho},
\]
\[
R_0=\frac{d^2}{4\rho}
\left(\frac1{e_1}+\frac1{e_2}+\frac1{e_4}+\frac1{e_5}\right).
\tag{6.9}
\]
For \(Q_2\) put
\[
U=u^{-2},
\qquad
A=\frac{e_2+e_5}{\rho},
\qquad
B=\frac{e_1+e_4}{\rho},
\]
\[
R_0=\frac{d^2}{4\rho}
\left(\frac1\lambda+\frac1{e_3}+\frac1\alpha+\frac1\beta\right).
\tag{6.10}
\]
Every inverse above is of an active point-side coefficient. Exact zero of one of those coefficients is processed before the division and is a proper coefficient child.

The four normal completions give the exact scalar children
\[
R_{Q_1}(U)=4(1+U)(A+BU)-R_0U,
\tag{6.11}
\]
\[
R_{Q_2}(U)=4\frac{(1+U)^2}{U}(AU+B)-R_0(1+U).
\tag{6.12}
\]
Global nonnegativity for \(U>0\) gives
\[
R_0\le4(\sqrt A+\sqrt B)^2.
\]
Set
\[
\eta=4(\sqrt A+\sqrt B)^2-R_0\ge0.
\]
Then, identically,
\[
\boxed{R_{Q_1}=4(\sqrt B\,U-\sqrt A)^2+\eta U,}
\tag{6.13}
\]
\[
\boxed{
R_{Q_2}=\frac{1+U}{U}
\left[4(\sqrt A\,U-\sqrt B)^2+\eta U\right].}
\tag{6.14}
\]
At an interior equality zero, \(A,B>0\), and in regular tangent coordinates the curvatures are
\[
32\sqrt{AB}>0\quad(Q_1),
\qquad
32\sqrt A(\sqrt A+\sqrt B)>0\quad(Q_2).
\tag{6.15}
\]
Together with the four positive normal directions and the bounded scale relations (6.7), this gives a positive complete five-variable quadratic residue at the exact critical center. Lemma 3.1 closes the centered remainder, while (6.8) already has explicit rank-two roots. Thus
\[
\boxed{\operatorname{sosx}(Q_0),\operatorname{sosx}(Q_1),\operatorname{sosx}(Q_2)\le2.}
\tag{6.16}
\]

---

# 7. Exact Bell gauge transfer: the collar theorem

This is the structural step that removes the old flat-collar obstruction.

Suppose a one-copy support map satisfies
\[
G_F(q)=c(q)^2G_E(\phi(q)),
\qquad c\ne0,
\qquad\det D\phi\ne0.
\tag{7.1}
\]
Form the \(5\times6\) value-plus-first-derivative contact matrix
\[
J_G(q)=
\begin{pmatrix}
G(q)^T\\
\partial_1G(q)^T\\
\vdots\\
\partial_4G(q)^T
\end{pmatrix}.
\]
Differentiation gives
\[
J_F(q)=
\underbrace{\begin{pmatrix}
c^2&0\\
2c\nabla c&c^2D\phi
\end{pmatrix}}_{A(q)}
J_E(\phi(q)).
\tag{7.2}
\]
Hence
\[
\det A(q)=c(q)^{10}\det D\phi(q).
\tag{7.3}
\]
Every signed maximal cofactor of a \(5\times6\) matrix scales by the determinant of a left \(5\times5\) multiplier. Therefore the physical cofactor/contact vectors satisfy
\[
X_F(q)=\lambda_F(q)X_E(\phi(q)),
\tag{7.4}
\]
where the physical orientation is chosen so \(\lambda_F>0\). With independent point and support copies,
\[
\boxed{
\tau_F(q_{\rm p},q_{\rm s})
=
\lambda_F(q_{\rm p})\,c(q_{\rm s})^2\,
\tau_E(\phi(q_{\rm p}),\phi(q_{\rm s})).
}
\tag{7.5}
\]
This is an exact whole two-copy identity.

The three Bell gauges are
\[
\phi_1=(u,s,as/p,1/p),
\qquad
G_{F_1}=p^2G_E\circ\phi_1,
\qquad
\det D\phi_1=\frac{s^3}{p^3},
\]
so
\[
\lambda_{F_1}=s^3p^7>0;
\tag{7.6}
\]
\[
\phi_4=(1/s,w,1/(bs),(As^2+w)/(bs)),
\]
\[
G_{F_4}=b^2s^4G_E\circ\phi_4,
\qquad
\det D\phi_4=-\frac1{b^3s^2},
\]
so on the physical chamber \(b<0\),
\[
\lambda_{F_4}=-b^7s^{18}>0;
\tag{7.7}
\]
and
\[
\phi_3=(W/r,1/r,a,g/r),
\qquad
G_{F_3}=r^6G_E\circ\phi_3,
\qquad
\det D\phi_3=r^{-4},
\]
so
\[
\lambda_{F_3}=r^{26}>0.
\tag{7.8}
\]

Therefore every nonzero collar gauge evaluation is exactly an EVENT evaluation multiplied by a positive pure-point factor and a pure-support square. Arbitrarily different infinitesimal rates of the two copies are irrelevant.

The strict collar boundaries are finite and algebraic. They are explicitly:
\[
\begin{array}{c|c}
\multicolumn{2}{c}{F_1}\\ \hline
a=0&e_4\text{ ray}\\
a=u^2&e_3\text{ ray}\\
u(a+1)=B&e_2\text{ ray}\\
B=0&\text{proper coordinate face}\\
u=0&B_0
\end{array}
\qquad
\begin{array}{c|c}
\multicolumn{2}{c}{F_4}\\ \hline
bw+1=0&e_3\text{ ray}\\
A+b(1+w^2)=0&e_2\text{ ray}\\
w=0&\text{flat }(d,e_2,e_3)\text{ child}\\
A=0&\text{proper coordinate face}\\
b=0&\text{proper coordinate/projective child}
\end{array}
\tag{7.9}
\]
and
\[
\begin{array}{c|c}
\multicolumn{2}{c}{F_3}\\ \hline
a=0&e_4\text{ ray}\\
g=0&e_5\text{ ray}\\
W(a-g)=1&e_2\text{ ray}\\
a=g&\text{terminal }F_3\\
W=0&\text{projective base chart}.
\end{array}
\tag{7.10}
\]
These equations are obtained by setting to zero the numerator or denominator factors in the three displayed birational maps and substituting into the EVENT rows (4.1). Exact zero of a gauge denominator is never treated by division; it is routed before the quotient is taken.

---

# 8. Terminal \(F_3\)

At the only noncoordinate strict \(F_3\) boundary, \(a=g\), write the point copy as
\[
s=Rp,\qquad a=-h+R^2,\qquad g=-h,
\]
and the independent support copy as
\[
t=S(p+V),\qquad W'=W+U,\qquad
a'=-(h+K)+S^2,\qquad g'=-(h+K),
\]
with \(W,h,p>0\). Put
\[
X=S-R,\qquad
Y=U-\frac WhK,\qquad
D=\sigma-\rho,
\]
and center the full radial square at
\[
d=D+\frac\rho W Y.
\tag{8.1}
\]
The degree-four event range form is the literal positive sum
\[
E_4=c_1q_1^2+c_2q_2^2,
\tag{8.2}
\]
where
\[
q_1=Sp(hU-WK),\qquad
c_1=2^{15}W^2h(Wh+1),
\]
\[
q_2=p^2(S^2-R^2),\qquad
c_2=2^{15}W^2h^2(W+h)(Wh+1)^2.
\]
The next range-kernel term has exact center
\[
X_c=\frac{RK}{2h(Wh+1)}-\frac{RV}{p},
\qquad Y_c=0.
\tag{8.3}
\]
After this completion, its effective residue is
\[
E_{6,\mathrm{eff}}
=2^{15}W^2p^2\Delta_{F_3}R^4K^2,
\qquad
\Delta_{F_3}=4W^3h-p^2(W+h).
\tag{8.4}
\]
The oriented contact-normal inequality gives
\[
p^2\le\frac{W^3h}{(Wh+1)(W+h)},
\tag{8.5}
\]
hence
\[
\boxed{
\Delta_{F_3}\ge
\frac{W^3h(4Wh+3)}{Wh+1}>0.
}
\tag{8.6}
\]
Thus the terminal event quadratic class is strictly positive in every unresolved kernel direction. The transfer back to the complete physical slack is also exact. In the orientation used above, the normalized full contact normal has residue factor
\[
c_{\rm norm}=\frac{1}{2^{13}W^2h^2(W+h)(Wh+1)},
\]
and the first full kernel layer is
\[
\boxed{
\mathcal F_{\rm eff}
=
\frac{4p^2\Delta_{F_3}}{h^2(W+h)(Wh+1)}
R^4\rho K^2\ge0.}
\tag{8.6a}
\]
The factor \(\rho\) is pure point-side. On the unresolved strict row it may therefore be removed before taking multiplication residue; the remaining coefficient is strictly positive by (8.6). The deformation of the event normal lies in the already completed range ideal and does not alter this first kernel coefficient. Hence the complete physical contact, not merely the EVENT truncation, has the asserted positive intrinsic kernel class.

To make the tensor typing intrinsic, use the genuine one-copy functions
\[
z_0=\rho W/h,\qquad
z_1=W/h,\qquad
z_2=R^4p^4(1+Wh)/h^2,\qquad
z_3=h.
\tag{8.7}
\]
With
\[
\widetilde X=X+\frac RpV-\frac{RK}{2h(1+Wh)},
\]
the corrected normal variables satisfy
\[
dz_0=(W/h)d,\qquad dz_1=Y/h,
\]
\[
\widetilde X=
\frac{R}{4z_2}dz_2
-
\frac{Rh^2}{4(1+Wh)}dz_1,
\qquad dz_3=K.
\tag{8.8}
\]
This is an invertible triangular cotangent change. Hence the positive terminal quadratic form is the intrinsic quadratic class in genuine one-copy coordinates. Lemma 3.1 gives
\[
\boxed{\operatorname{sosx}(F_3^{\rm term})\le2.}
\tag{8.9}
\]
Loss of a strict unit is a coordinate or projective child.

---

# 9. One-projective contacts

Assume \(0\ne u\in\mathfrak m\) while \(v\) has nonzero residue. Reciprocal observable charts and Bell relabelings reduce every one-projective contact to this form using only pure-support squares. Exact coincidence or opposition of one party's two measurements is a local/polyhedral child, so we work on the open one-projective chart.

Put
\[
c_X=e_1+\beta,
\qquad
c_Y=e_3+e_4,
\]
\[
e_1=(1-k_X)c_X,
\quad \beta=k_Xc_X,
\qquad
e_3=k_Yc_Y,
\quad e_4=(1-k_Y)c_Y,
\]
with \(0\le k_X,k_Y\le1\). Define the harmonic point-side scales
\[
h_L=\frac{e_1\beta}{e_1+\beta}=c_Xk_X(1-k_X),
\qquad
h_R=\frac{e_3e_4}{e_3+e_4}=c_Yk_Y(1-k_Y),
\tag{9.1}
\]
with value zero when the corresponding sum is exactly zero. In the affine state chart \(\psi=(\theta,X,Y,1)^T\), use the exact completed normals
\[
N_X=X+k_Xu,
\qquad
N_Y=Y-k_Yu\theta.
\tag{9.2}
\]
At \(u=0\) the signed packet has positive residue \(2(1+v^2)(X^2+Y^2)\), so a hard zero has \(\bar X=\bar Y=0\). The forcing amplitudes belonging to \(\lambda,e_2,e_5,\alpha\) then have nonzero real residues. Therefore any one of these coefficients strictly earlier than the signed scale closes by positive residue.

### Proposition 9.1: separated projective entry

After all earlier positive packets and exact children are removed, set
\[
v_0=\bar v\ne0,
\qquad
\theta_0=\bar\theta\ne0,
\qquad
\mathcal H=h_L+h_R\theta_0^2.
\tag{9.3}
\]
On a hard branch,
\[
\boxed{
\mathcal H>0,
\qquad
\tau_{\rm test}:=\frac{u\mathcal H}{d}\in B^\times,
\qquad
\delta:=\frac d{\mathcal H}\in\mathfrak m,}
\tag{9.4}
\]
and
\[
\boxed{
\lambda,e_2,e_5,\alpha\in\frac{d^2}{\mathcal H}B.}
\tag{9.5}
\]
Moreover there is an exact point-side normal--\(u\) critical center with the same residue, and after division by \(\mathcal H\) one may write exactly
\[
u=\delta\tau,
\qquad
\lambda/\mathcal H=\delta^2\ell,
\quad e_2/\mathcal H=\delta^2p,
\quad e_5/\mathcal H=\delta^2q,
\quad \alpha/\mathcal H=\delta^2a,
\tag{9.6}
\]
with \(\tau\in B^\times\) and \(\ell,p,q,a\in B_{\ge0}\).

#### Proof
Complete the two positive split normals in the one-copy field. The normal-center correction is later than the first scalar radial layer. The first possibly nonzero scalar face in the independent support radial variable \(u_s\) is therefore a subface of
\[
4(1+v_0^2)\mathcal H u_s^2
-4dv_0\theta_0u_s
+4(\lambda\theta_0^2+e_2v_0^2\theta_0^2+e_5v_0^2+\alpha).
\tag{9.7}
\]
A hard zero cannot have a strictly earlier positive constant or quadratic face. The negative linear term also cannot be the unique earliest term: if it were tied only with the nonnegative constant, the resulting real initial polynomial \(C-4v_0\theta_0t\) becomes negative for a sufficiently large positive real rescaling \(u_s\mapsto tu_s\), contradicting global nonnegativity. Hence the harmonic quadratic and signed linear terms tie, proving \(u\mathcal H/d\in B^\times\). Since \(u\in\mathfrak m\), this gives \(d/\mathcal H\in\mathfrak m\). A forcing coefficient earlier than \(du\) would itself be a positive initial face, so (9.5) follows.

Normalize by the pure-point scale \(d^2/\mathcal H\) and put \(u_s=\delta t\). The multiplication-residue scalar envelope is
\[
\mathscr E_0(t)=
4(1+v_0^2)t^2-4v_0\theta_0t
+4(\bar\ell\theta_0^2+\bar p v_0^2\theta_0^2+\bar qv_0^2+\bar a).
\tag{9.8}
\]
At an unresolved zero it is nonnegative and vanishes, so its derivative has the simple zero
\[
\mathscr E_0'(t)=8(1+v_0^2)t-4v_0\theta_0,
\qquad
\mathscr E_0''(t)=8(1+v_0^2)>0.
\tag{9.9}
\]
For the exact normal-minimized envelope the packet formulas give
\[
\mathscr E'(t)=8(1+v_0^2)t-4v_0\theta_0+\delta R(t),
\qquad R(t)\in B.
\tag{9.10}
\]
The derivative equation therefore has unit Jacobian at the residue zero. The real-closed-field implicit-function theorem gives an exact point-side critical center with the same residue. This proves (9.4)-(9.6). \(\square\)

## 9.1 Actual normal packets

Complete the actual positive channel pairs before introducing any mixed normal coordinate. A model identity is
\[
\begin{aligned}
e_1(1+u^2)X^2+\beta(X+uB)^2
={}&\bigl(e_1(1+u^2)+\beta\bigr)
\left(X+\frac{\beta uB}{e_1(1+u^2)+\beta}\right)^2\\
&+
\frac{e_1(1+u^2)\beta}{e_1(1+u^2)+\beta}u^2B^2.
\end{aligned}
\tag{9.11}
\]
The denominator lies wholly in one tensor factor.

On a full-normal row every active normal coefficient \(c_i\) satisfies \(d/c_i\in\mathfrak m\). The exact coefficient matrix in the actual one-copy amplitudes is therefore a positive diagonal matrix plus a zero-residue perturbation and is strictly diagonally dominant. Pointwise,
\[
\begin{aligned}
f^TMf={}&
\sum_i\left(m_{ii}-\sum_{j\ne i}|m_{ij}|\right)f_i^2\\
&+\sum_{i<j}\frac{|m_{ij}|+m_{ij}}2(f_i+f_j)^2
+\sum_{i<j}\frac{|m_{ij}|-m_{ij}}2(f_i-f_j)^2.
\end{aligned}
\tag{9.12}
\]
Every root is a real sum or difference of actual one-copy amplitudes, multiplied only by pure-side square roots. Hence each root has tensor rank at most two.

On a one-missing row, the remaining nondominant sector is a genuine point-side \(2\times2\) PSD packet:
\[
af^2+2bfg+cg^2
=
a(f+(b/a)g)^2+\frac{ac-b^2}{a}g^2.
\tag{9.13}
\]
Again the first root contains only two separated functions. The both-missing row is incompatible with the scalar stationarity equation.

## 9.2 Scalar singularity

After the exact normal completion, the scalar child can be written
\[
\begin{aligned}
M(\tau,v,\theta)={}&
(\ell+pv^2+h_Y(1+v^2)\tau^2)\theta^2-\tau v\theta\\
&+a+qv^2+h_X(1+v^2)\tau^2.
\end{aligned}
\tag{9.14}
\]
The determinant coefficient table has nine entries:
\[
\begin{array}{c|l}
(2,2)&4h_Xh_Y\\
(2,1)&8h_Xh_Y\\
(2,0)&4h_Xh_Y\\
(1,2)&4(h_Xp+h_Yq)\\
(1,1)&4ah_Y+4\ell h_X+4h_Xp+4h_Yq-1\\
(1,0)&4(ah_Y+\ell h_X)\\
(0,2)&4pq\\
(0,1)&4(ap+\ell q)\\
(0,0)&4a\ell.
\end{array}
\tag{9.15}
\]
The same \(B_2\) Cauchy-Binet argument as in Section 5 shows that, up to left-right symmetry, the only nonzero singular stationary pattern is
\[
\ell=p=h_X=0,\qquad a,q,h_Y>0.
\tag{9.16}
\]
Write \(a=A^2\), \(q=Q^2\), \(h_Y=H_0^2\). The exact residue scalar becomes
\[
\boxed{
\frac\Phi4
=(H_0\tau\theta-Qv)^2+(H_0\tau v\theta-A)^2.
}
\tag{9.17}
\]
The Jacobian of the two square arguments with respect to \((\theta,v)\) is
\[
2H_0Q\tau v\ne0.
\tag{9.18}
\]
Hence the two squares are exact local normal coordinates.

On their exact zero curve,
\[
Q=\frac1{2H_0(1+v^2)},
\qquad
A=Qv^2,
\]
\[
\theta=\frac{Qv}{H_0\tau}.
\tag{9.19}
\]
Every delayed scalar packet restricts to a Laurent polynomial of the form
\[
R(\tau)=A_1\tau^2+B_1\tau^{-2}+C_1,
\qquad A_1,B_1\ge0.
\tag{9.20}
\]
If \(R\ge0\) and has an interior zero, then \(A_1,B_1>0\), \(A_1\tau_0^4=B_1\), and
\[
R''(\tau_0)=8A_1>0.
\tag{9.21}
\]
Thus the both-current and neither-current rows have positive scalar curvature after the full exact normal center is lifted. Lemma 3.1 applies to the complete centered contact.

The only remaining exactly-one-current branch is A17. Exact loss of any normal scale is a proper child.

Therefore every non-A17 one-projective contact has \(\operatorname{sosx}\le2\).

---

# 10. A17

A17 is the unique exactly-one-current continuation of the singular one-projective scalar. It is also the only terminal at which the exact critical value can survive as a positive infinitesimal. We therefore derive its scale pattern from the original critical equations before using the terminal square-ideal calculation.

Put \(u=r^2\), \(S=1+V^2\). On the strict chart the point-side quantities
\[
T,V,C_X,Z,1-Z
\]
are units with positive residues and \(0<\bar Z<1\). Normalize the current right split by
\[
e_3=Z^{-1},
\qquad
e_4=(1-Z)^{-1},
\tag{10.1}
\]
and let \(\chi=1\) in the \(\beta\)-current orientation and \(\chi=0\) in the \(e_1\)-current orientation.

### Lemma 10.1: terminal entry from two unit-Jacobian divisions

After all earlier positive packets and exact children are removed, there is an exact full critical center with the same residue as the tested contact. At that center the first exact orders are
\[
d\in uB,
\qquad \overline{d/u}=2TS,
\tag{10.2}
\]
\[
\lambda,e_2\in u^2\mathfrak m,
\qquad
e_5,\alpha,m\in u^2B,
\tag{10.3}
\]
and the delayed member \(c\) of \(\{e_1,\beta\}\) satisfies
\[
c=u^2N_2,
\qquad N_2\in B.
\tag{10.4}
\]

#### Proof of the first orders
The two current positive range squares have unit residue. Their exact normal equations give
\[
X_p+\chi u\in uB,
\qquad
Y_p-(1-Z)uTV\in uB.
\]
The first nonzero \(u\)-stationarity layer of the complete nine-channel slack is
\[
8T^2V^2S\,u-4TV^2d,
\]
which forces (10.2). After completing the two unit range directions, the first kernel Schur layer is
\[
\begin{pmatrix}
4(Su^2+\lambda+V^2e_2)&-2Vdu\\
-2Vdu&4(V^2e_5+\alpha)
\end{pmatrix}.
\tag{10.5}
\]
The exact affine state equations, with the fourth equation equal to the retained critical value \(m\), give (10.3). Finally the exact completion of the current/delayed split pair has residual form
\[
R_{\rm split}=u^2c\,U_{\rm split},
\qquad U_{\rm split}>0,
\quad U_{\rm split}\in B_p^\times.
\tag{10.6}
\]
If \(c/u^2\) were unbounded, this positive layer would close earlier, proving (10.4).

Now write
\[
d=uD_1,
\quad
\lambda=u^2\ell_2,
\quad
e_2=u^2m_2,
\quad
e_5=u^2Q_2,
\quad
\alpha=u^2R_2,
\quad
m=u^2\eta_2,
\]
and \(\psi=(TV,ux_1,uy_1,1)^T\). Define the exact normalized five-equation map
\[
\begin{aligned}
\mathcal T_0&=u^{-2}(P_\star\psi)_0,
&\mathcal T_1&=u^{-1}(P_\star\psi)_1,\\
\mathcal T_2&=u^{-1}(P_\star\psi)_2,
&\mathcal T_3&=u^{-2}((P_\star\psi)_3-m),\\
\mathcal T_4&=u^{-2}\psi^T(\partial_vP_\star)\psi.
\end{aligned}
\tag{10.7}
\]
At \(u=\ell_2=m_2=0\) its residue system has the unique solution
\[
(x_1,y_1,D_1,Q_2,R_2)
=
(-\chi,(1-Z)TV,2TS,T^2,T^2V^4+\eta_2/4),
\tag{10.8}
\]
and
\[
\boxed{
\det D_z\mathcal T
=-\frac{2^{10}C_XV^2S^2}{Z(Z-1)}\ne0.}
\tag{10.9}
\]
Thus unit-Jacobian division controls all first corrections. The remaining exact \(u\)-stationarity equation has linear part
\[
8TV^2[Su-T(\ell_2+V^2m_2)]
\quad(\beta),
\]
or
\[
8T[Su-TV^2(\ell_2+V^2m_2)]
\quad(e_1).
\]
If \(u\) were later than \(\ell_2+m_2\), division by their sum would leave a strictly negative real combination of nonnegative residues, impossible. Hence
\[
\ell_2,m_2\in uB,
\quad
\boxed{L_3=\lambda/u^3\in B,
\qquad M_3=e_2/u^3\in B.}
\tag{10.10}
\]
Compare the two nonnegative quantities \(u\) and \(\eta_2\). If \(\eta_2\) is earlier, the normalized critical value is the literal positive scalar \(1\), so the branch closes. On the hard branch
\[
\boxed{K_3=\eta_2/u\in B_{\ge0},
\qquad m=u^3K_3.}
\tag{10.11}
\]
A second application of the same unit-Jacobian system, after dividing the corrections by \(u\), gives bounded coefficients \(D_2,Q_3,\widehat R_3,X_3,\widehat Y_3\) and the exact defect congruences
\[
D_2-2T\left[L_3+V^2M_3+T\{V^2-S^2Z(1-Z)\}\right]\in uB,
\tag{10.12}
\]
\[
T(L_3+V^2M_3)-S\in uB\quad(\beta),
\qquad
TV^2(L_3+V^2M_3)-S\in uB\quad(e_1).
\tag{10.13}
\]
These are precisely the extra factors needed in the finite terminal cancellation ledger. \(\square\)

Consequently, in the \(\beta\)-current orientation, the exact terminal point data may be written
\[
\begin{aligned}
d&=2TSu+D_2u^2,&\lambda&=L_3u^3,&e_2&=M_3u^3,\\
e_1&=N_2u^2,&\beta&=C_X,\\
e_5&=T^2u^2+Q_3u^3,&
\alpha&=T^2V^4u^2+(\widehat R_3+K_3/4)u^3,
\end{aligned}
\tag{10.14}
\]
with the mirrored orientation obtained by interchanging the current split member. The critical value is (10.11).

### Exact relative graph

Keep the support copy independent. Put
\[
x_4=A_s-TV,
\qquad
x_3=v_s-V,
\]
\[
h_2=u_s-u+\frac{u}{A_s}x_4,
\qquad
h_1=\frac{Y_s-Y_p-(1-Z)TVh_2}{Z(1-Z)}.
\tag{10.15}
\]
In the \(\beta\)-orientation set
\[
x_0=X_s-X_p+h_2-\frac{u}{A_s}x_4,
\]
and in the \(e_1\)-orientation set
\[
x_0=X_s-X_p.
\tag{10.16}
\]
The inverse change is triangular and uses only the support-side unit \(A_s^{-1}\) and point-side unit \([Z(1-Z)]^{-1}\). Hence these five relative binomials generate the full local contact ideal. With \(u=r^2\) define
\[
\boxed{J=(C_Xx_0,h_1,h_2,dx_3,drx_4).}
\tag{10.17}
\]

Center the loss:
\[
\widetilde\Theta=\Theta-m.
\tag{10.18}
\]
Use the independent relative ideal \(J\) from (10.17). After clearing the pure-support squares \(A_s^4\) or \(A_s^2\), the only possible below-cost first-jet terms occur in a finite list. Their channel coefficients reduce to the cancellation triples
\[
(2Z+1,-2(Z-1),-3),
\quad
(Z+2,-(Z-1),-3),
\quad
(1,0,-1)
\tag{10.20}
\]
for the \(e_1\)-orientation, and
\[
\begin{gathered}
(4Z+1,-4(Z-1),-5),\\
(3Z+2,-3(Z-1),-5),\\
(2Z+3,-2(Z-1),-5),\\
(Z+4,-(Z-1),-5),\\
(1,0,-1),
\end{gathered}
\tag{10.21}
\]
for the \(\beta\)-orientation, together with
\[
Z-(Z-1)-1=0.
\tag{10.22}
\]
Every triple sums to zero. Hence
\[
\boxed{
A_s^4\widetilde\Theta_\beta\in J^2,
\qquad
A_s^2\widetilde\Theta_{e_1}\in J^2.
}
\tag{10.23}
\]
The common four-direction residue block is
\[
A_4=
\begin{pmatrix}
\dfrac{4S}{C_X}&0&0&0\\
0&4Z(1-Z)S&0&0\\
0&0&4T^2V^2S&\dfrac{2TV(V^2-1)}S\\
0&0&\dfrac{2TV(V^2-1)}S&\dfrac1S
\end{pmatrix}.
\tag{10.24}
\]
Its lower \(2\times2\) determinant is
\[
\frac{16T^2V^4}{S^2}>0,
\]
so \(A_4\succ0\). The centered fifth Schur pivots are
\[
\boxed{
q_\beta=\frac4{T^3S}>0,
\qquad
q_{e_1}=\frac4{T^3V^2S}>0.
}
\tag{10.25}
\]

To make the fifth-direction tensor rank manifest, use genuine one-copy functions
\[
F_1=X+\chi u,
\qquad F_2=Y,
\qquad F_3=Au,
\qquad F_4=Auv,
\qquad D_i=\delta F_i.
\tag{10.26}
\]
Here \(A\) is the affine state coordinate on the one-copy A17 chart (so its support-copy value is \(A_s\)); it is unrelated to the scalar symbols \(A,B\) used in the finite circuit normalizations.
Let \(g=drx_4\). After a triangular pure-side change of the first four directions, the exact centered form is
\[
F=D^TAD+2rg\,b^TD+\gamma g^2,
\qquad
\overline A\succ0,
\quad
\overline\gamma>0.
\tag{10.27}
\]
For
\[
k=dr/u,
\qquad
f_\lambda=u+\sum_i\lambda_iF_i,
\]
and
\[
t(\lambda)=
-(d/u)
(A_s\lambda_1,A_s\lambda_2,1+A_s\lambda_3,A_s\lambda_4),
\]
one obtains the exact identity
\[
\boxed{
g+r\,t(\lambda)^TD=-A_sk\,\delta f_\lambda.}
\tag{10.28}
\]
Each right-hand side is a genuine divided difference multiplied by a pure-side factor.

Choose five real slopes \(\lambda^{(j)}\) forming a small affine simplex around the residue barycentric slope. The affine-system determinant is
\[
5\varepsilon^4(dA_s/u)^4,
\tag{10.29}
\]
a product of pure-side units. The exact weights solve
\[
\sum_j\theta_j=\gamma,
\qquad
\sum_j\theta_jt_j=b,
\tag{10.30}
\]
and satisfy
\[
\overline{\theta_j}=\overline\gamma/5>0.
\tag{10.31}
\]
Subtracting the five squares leaves
\[
D^T\left(A-r^2\sum_j\theta_jt_jt_j^T\right)D,
\tag{10.32}
\]
whose residue is still \(\overline A\succ0\). Lemma 3.1 closes this four-dimensional residual in genuine divided differences. Finally restore the critical value by the pure-side square
\[
mA_s^4=(\sqrt m\,A_s^2)^2
\quad\text{or}\quad
mA_s^2=(\sqrt m\,A_s)^2.
\tag{10.33}
\]
Therefore
\[
\boxed{\operatorname{sosx}(A17)\le2.}
\tag{10.34}
\]
Loss of \(T,V,C_X\) as units or loss of \(0<\overline Z<1\) is an exact coordinate/projective child.

---

# 11. The simultaneous projective base locus

Use homogeneous measurement coordinates \([u_0:u_1]\), \([v_0:v_1]\) and the projective EVENT-kernel vector \([y_1:y_2:y_3]\). After removal of common positive factors, the five positive homogeneous EVENT coefficients are
\[
\begin{aligned}
G_1&=(u_0^2+u_1^2)v_0^2(v_0^2+v_1^2)y_3^2,\\
G_4&=(u_0^2+u_1^2)v_0^2(v_0^2+v_1^2)y_1^2,\\
G_5&=u_0^2v_0^2(v_0y_1+v_1y_2)^2,\\
G_3&=(v_0^2+v_1^2)(u_0v_0y_1+u_1v_1y_3)^2,\\
G_2&=(u_1v_1v_0y_1-u_1v_0^2y_2-u_0(v_0^2+v_1^2)y_3)^2.
\end{aligned}
\tag{11.0}
\]
Their simultaneous real zero set is immediate. If \(v_0=0\), then \(v_1\ne0\). The two equations \(G_3=G_2=0\) force \(y_3=0\): if \(u_1\ne0\), \(G_3=0\) does so, while if \(u_1=0\), projectivity gives \(u_0\ne0\) and \(G_2=0\) does so. This gives
\[
\boxed{B_\infty:\ v_0=0,\ y_3=0.}
\tag{11.1}
\]
If \(v_0\ne0\), \(G_1=G_4=0\) give \(y_1=y_3=0\). Projectivity forces \(y_2\ne0\). Now \(G_2=0\) gives \(u_1=0\); projectivity of \([u_0:u_1]\) then gives \(u_0\ne0\), and \(G_5=0\) gives \(v_1=0\). Hence
\[
\boxed{B_0:\ u_1=v_1=0,\ [y]=[0:1:0].}
\tag{11.2}
\]
The alternatives are disjoint and exhaustive. There is no third real simultaneous-projective component.

---

# 12. Endpoint \(B_0\)

Use the separated point chart
\[
s_p=T,
\qquad
u_p=r^2,
\qquad
a_p=Ar^4,
\qquad
B_p=Lr^2,
\]
and independent support variables \((U,V)\) in the full physical suspension. Put
\[
x=T^2,
\qquad
w=1-Lx,
\qquad
p=1-A-Lx,
\qquad
h=1-L(1+x).
\tag{12.0}
\]
Then \(w=A+p\). We first derive, rather than assume, the strict point chamber. Let \(X=(X_d,X_{e_1},\ldots,X_{e_5})\) be the oriented EVENT cofactor vector. Direct substitution of the \(B_0\) chart gives
\[
\left.r^{-12}X\right|_{r=0}
=C_0(0,0,0,A,p,0)^T,
\qquad
C_0=2^{13}T^9(1+T^2)^3(LT^2-1)^3(LT^2+L-1).
\tag{12.0b}
\]
The exact factor \(LT^2=1\), equivalently \(w=0\), is processed before any division by \(w\). There \(B=u/T^2\), the primitive square \(G_1\) and its complete first jet vanish, so the corresponding contact-matrix column \(\mathcal C_1\) is zero. The remaining columns satisfy the exact relation
\[
-2\mathcal C_0+\mathcal C_2+\mathcal C_3=0,
\tag{12.0c}
\]
while the \(4\times4\) minor on rows \((0,1,2,3)\) and columns \((2,3,4,5)\) is
\[
-2^{11}A^2T^5u^7(1+T^2)^2(1+Au^2)
(AT^2u^2-1)(2Au^4+Au^2+u^2+1)\ne0.
\tag{12.0d}
\]
At the projective endpoint \(u=r^2\in\mathfrak m\), every displayed factor is nonzero. Hence the five nonzero columns have rank four and their unique relation has signs \((- ,+,+,0,0)\), so no nonzero multiple is componentwise nonnegative. The only nonnegative kernel direction is the zero column \(\mathcal C_1\), namely the exact \(e_1\) coordinate child.

Away from \(w=0\), put
\[
m=L(x+2)-1,
\qquad
q_0=x+2-Lx^2.
\]
After removal of one common positive factor, the exact leading cofactor signs are
\[
\begin{array}{c|c|c}
\text{channel}&\text{order in }r&\text{remaining factor}\\ \hline
 d&16&2Aw^2ph\,q_0\\
 e_1&16&ALx^2w^2ph\\
 e_2&16&Aw^2pm\\
 e_3&12&Aw^3h\\
 e_4&12&w^3ph\\
 e_5&20&A(1+x)w^4ph.
\end{array}
\tag{12.0e}
\]
The \(e_3,e_4\) entries show that \(A\) and \(p\) have the same sign. They cannot both be negative: then \(w=A+p<0\), \(L=(1-w)/x>0\), \(h=w-L<0\), and \(m=L(x+2)-1>0\). Making the \(e_3,e_4\) entries nonnegative would reverse the common orientation, while the \(e_2\) factor is positive and would become negative. Thus
\[
A>0,
\qquad
p>0,
\qquad
w=A+p>0.
\]
If \(h<0\), then \(L>1/(1+x)\), hence \(m>1/(1+x)>0\), and the same \(e_3,e_4\) versus \(e_2\) sign contradiction occurs. Therefore \(h\ge0\). The exact boundary \(h=0\) is treated below and has no unresolved nonzero support. On the remaining strict row \(h>0\), the \(e_3,e_4\) entries fix the orientation. Then \(e_1\) forces \(L\ge0\); equality would give \(m=-1\) and make \(e_2\) negative, so \(L>0\). Finally \(e_2\) gives \(m\ge0\). In particular the unresolved open chamber has
\[
0<A<1,
\qquad
T>0,
\qquad
L,p,w,h>0,
\qquad
m\ge0,
\tag{12.0f}
\]
where \(A<1\) follows from \(p=1-A-Lx>0\). No residue-unit assumption on \(h\) is made: it may be positive infinitesimal until its exact boundary is selected.

The exact transverse comparison has one nonchild tie,
\[
\chi_p=-r^4\varepsilon_0,
\qquad
\varepsilon_0\in B_{p,\ge0}.
\tag{12.0a}
\]
If the transverse term is earlier, the literal positive \(\lambda\)-packet closes the branch; if it is later, it belongs to the next actual coefficient layer. Thus (12.0a) is the only unresolved transverse case.

After a pure-point radial congruence and one common positive point clearing, the exact full \(4\times4\) strict-transform matrix has a literal four-square decomposition. With
\[
\omega=(U,-1,0,-UV)^T,
\]
one has
\[
\boxed{
\begin{aligned}
S_0={}&
\frac{AT(1+T^2)p}{w}(1+U^2)e_0e_0^T
+A\varepsilon_0(1+V^2)\omega\omega^T\\
&+\varepsilon_0p(1+U^2)(1+V^2)e_1e_1^T
+\frac{A\varepsilon_0(1+T^2)wp}{h}U^2e_2e_2^T.
\end{aligned}
}
\tag{12.1}
\]
Every coefficient is nonnegative and lies wholly in the point factor. Hence (12.1) is already a tensor SOS, with rank-one or rank-two roots after the trivial support-function split.

The only dangerous open-chamber boundary is \(h=0\). Its exact next tiers are not discarded. Put
\[
\zeta=\frac{1-A(1+T^2)}{1+T^2}\ge0.
\]
After common positive clearing, the first \(h=0\) cofactor tier is proportional to
\[
A(0,0,\zeta,A,\zeta,0).
\tag{12.2}
\]
The next positive coefficients include
\[
[e_d]_{r^{20}}\propto A^2(3T^2+2)\zeta,
\qquad
[e_1]_{r^{20}}\propto A^2T^4\zeta,
\qquad
[e_5]_{r^{24}}\propto A^2\zeta.
\tag{12.3}
\]
If \(A=0\) we are on the exact \(e_4\) child. If \(\zeta=0\), then \(p=0\), the exact \(e_3\) child. If \(A,\zeta>0\), the leading positive EVENT squares force all nonzero support coordinates to vanish, while the later \(e_5\) square remains strictly positive, a contradiction. Hence no unresolved nonzero support remains on \(h=0\).

Therefore
\[
\boxed{\operatorname{sosx}(B_0)\le2.}
\tag{12.4}
\]

---

# 13. Endpoint \(B_\infty\)

This is the endpoint where positive infinitesimal pivots must be distinguished from positive residues.

## 13.1 Reciprocal chart, transverse endpoints, and exact coefficient firewall

Use the reciprocal point chart
\[
s_p=cr,
\qquad
u_p=r^{-1},
\qquad
a_p=((J+c^2)r^2)^{-1},
\qquad
B_p=\frac{1+Kr^2}{(J+c^2)r^3},
\tag{13.0a}
\]
with \(0\ne r\in\mathfrak m_{>0}\). The first reciprocal cofactor tier gives
\[
J\bigl(2(J-1)(J-K),0,K(J-1),J-K,0,J(J-1)(J-K)\bigr),
\tag{13.0b0}
\]
up to one common positive point factor and the fixed physical orientation. Here and below let
\[
Y=(Y_d,Y_{e_1},Y_{e_2},Y_{e_3},Y_{e_4},Y_{e_5})
\]
denote the exact reciprocal cofactor vector after removal of that common positive point factor. Hence the first-tier inequalities are \(J\ge1\) and \(0\le K\le J\). They are not yet exact child statements. The exact boundary routing is as follows.

On \(J=1\), put
\[
b_1=K(1+c^2r^2)-1.
\tag{13.0b1}
\]
For \(K>0\), direct factorization of the exact reciprocal cofactors gives
\[
Y_{e_1}=P_1b_1,
\qquad
Y_{e_3}=-P_3b_1,
\qquad P_1,P_3>0,
\tag{13.0b2}
\]
where
\[
P_1=Kc^6r^6(1+Kr^2)(1+(1+c^2)r^2)(1+r^2-Kc^2r^4),
\]
\[
P_3=(1+(1+c^2)r^2)(1+r^2-Kc^2r^4)^2
(1+r^2+c^2r^2(1-K)).
\]
Thus \(b_1\ne0\) gives incompatible cofactor signs. If \(b_1=0\), then
\[
K=(1+c^2r^2)^{-1}
\]
and all reciprocal cofactors except \(Y_{e_2}\) vanish exactly, while
\[
Y_{e_2}=-\frac{c^2r^2(1+r^2)(1+(1+c^2)r^2)^2}
{(1+c^2r^2)^3}<0.
\tag{13.0b3}
\]
After the common orientation reversal this is exactly the \(e_2\) coordinate ray. The remaining subcase \(J=1,K=0\) is exactly the \(e_3\) ray,
\[
Y=(0,0,0,(1+r^2)^2(1+(1+c^2)r^2)^2,0,0).
\tag{13.0b4}
\]

On \(K=J\), put
\[
a_1=1-J+Jc^2r^2.
\tag{13.0b5}
\]
Exact factorization gives
\[
Y_{e_1}=Q_1a_1,
\qquad
Y_{e_3}=-Q_3,
\tag{13.0b6}
\]
with
\[
Q_1=Jc^6r^6(1+Jr^2)(1+(J+c^2)r^2)(1+Jr^2-Jc^2r^4)>0,
\]
\[
Q_3=Jc^2r^2(1+(J+c^2)r^2)(1+Jr^2-Jc^2r^4)^2
\{J+Jr^2+c^2r^2(1-J)\}>0.
\]
If \(a_1>0\), the two displayed coefficients have incompatible signs. If \(a_1=0\), direct substitution in the remaining four cofactors gives zero identically and, after the common orientation reversal, this is exactly the \(e_3\) coordinate ray. If \(a_1<0\), set
\[
q=J(1-c^2r^2)>1.
\]
The exact \(e_2\) cofactor factors as
\[
\boxed{
Y_{e_2}=q(q-1)U_2,
\qquad
U_2=\frac{1+(J+c^2)r^2}{1-c^2r^2}\,\Pi(r),}
\tag{13.0b7}
\]
where
\[
\begin{aligned}
\Pi(r)={}&J+(J^2+J)r^2+(J^2-Jc^4+c^4)r^4\\
&+(Jc^4-2J^2c^4)r^6
+(J^2c^6-J^2c^4-Jc^6)r^8.
\end{aligned}
\tag{13.0b8}
\]
The denominator \(1-c^2r^2\) is a pure-point unit and
\[
\overline{U_2}=\bar J>0.
\tag{13.0b9}
\]
Hence \(Y_{e_2}>0\) while \(Y_{e_3}<0\), again impossible. Every assertion in this subcase is now a finite displayed factorization.

The exact boundary \(c=0\) is processed before \(c^{-1}\) is introduced and is the already covered one-projective chart. By contrast, \(K=0\) with \(J>1\) is **not** deleted: it remains part of the reciprocal chamber. Consequently every unresolved reciprocal point lies in
\[
\boxed{c>0,\qquad J>1,\qquad0\le K<J.}
\tag{13.0b}
\]
Positive infinitesimal \(c,J-1,J-K\) remain explicit pure-point scales and are handled by finite coefficient layers. Every inverse in (13.0a) is therefore pure point.

For the transverse full-contact coordinate use one raw point variable \(\varepsilon\) and put
\[
\chi=\frac{J}{J+c^2}\varepsilon,
\qquad
\zeta=\frac{\varepsilon}{c},
\qquad
\eta=1-J\zeta.
\tag{13.0c}
\]
At \(\varepsilon=0\), all six non-\(\lambda\) full cofactors vanish identically. In the chosen nonnegative orientation
\[
X^{\rm full}|_{\varepsilon=0}
=\Lambda_\lambda(r,c,J,K)(1,0,0,0,0,0,0)^T,
\qquad \Lambda_\lambda\ge0.
\tag{13.0d}
\]
If \(\Lambda_\lambda>0\), this is the exact \(\lambda\) coordinate ray; if \(\Lambda_\lambda=0\), the coefficient filtration continues in the already retained tangent EVENT tensor. The first radial coefficient is
\[
\left.r^9X^{\rm full}\right|_{r=0,\,\zeta=0}
=\frac{2^{16}J^3c^{10}(J-1)(J-K)}{(J+c^2)^8}
(1,0,0,0,0,0,0)^T.
\tag{13.0e}
\]
At the opposite endpoint the first radial tier is only an \(e_5\) ray,
\[
\left.r^9X^{\rm full}\right|_{r=0,\,\zeta=1/J}
=\frac{2^{16}J^3c^{10}(J+1)(J-K+1)}{(J+c^2)^8}
(0,0,0,0,0,0,1)^T,
\tag{13.0f}
\]
so it must not be declared an exact child. Retain that \(e_5\) square and restrict the next \(r^2\) tier to its exact kernel. With support coordinates \(U,V\), a kernel basis is
\[
C_5=
\begin{pmatrix}
1&0&0\\
0&1&0\\
U&UV&-V\\
0&0&1
\end{pmatrix}.
\tag{13.0g}
\]
After one common positive point normalization, the next coefficient vector is
\[
(\lambda,d,e_1,e_2,e_3,e_4,e_5)
=(J^2q,2Jq,0,JK,q,0,0),
\qquad q=J-K+1.
\tag{13.0h}
\]
For \(M_2=C_5^TP_\star C_5\), the leading principal minors are
\[
\Delta_1=4J(1+U^2)\{J(J-K+1)+(J+1)U^2V^2+KV^2\},
\tag{13.0i}
\]
\[
\Delta_2=16J(J+1)(J+U^2V^2)(1+U^2)(1+V^2)q
(J+KV^2+U^2V^2),
\tag{13.0j}
\]
\[
\Delta_3=64J^2V^2(J+1)^2(1+U^2)^2(1+V^2)^2q^2
(J+KV^2+U^2V^2).
\tag{13.0k}
\]
Hence \(M_2\succ0\) for \(V\ne0\). At \(V=0\) its unique physical kernel is \((0,-U,0,1)^T\). The first later \(e_1\) coefficient is genuinely positive: all coefficients through order \(r^5\) in \(r^9X^{\rm full}_{e_1}\) vanish, while
\[
\boxed{
[r^6]\!\left(r^9X^{\rm full}_{e_1}\right)_{\zeta=1/J}
=\frac{2^{16}Jc^{14}(J+1)(J-K+1)}{(J+c^2)^7}>0.}
\tag{13.0l}
\]
On the retained kernel its support square is
\[
4U^2(1+U^2).
\tag{13.0m}
\]
This is strictly positive for \(U\ne0\); if \(U\ne0\) but \(\bar U=0\), the pure-support \(U\)-strict transform has residue \(4\); if \(U=0\) exactly, one is on the exact one-measurement coordinate child. Thus \(\zeta=1/J\) is closed without promoting its leading ray to an exact ray. Consequently every unresolved coupled row satisfies
\[
0<\zeta<1/J.
\tag{13.0n}
\]

Now filter only the original nonnegative Bell coefficients. On every smooth reciprocal contact row the complete \(6\times7\) value/first-jet system satisfies the exact Cramer identity
\[
\boxed{D_0D_1\,d=2J\zeta\,N_d\,\lambda,}
\tag{13.0o}
\]
where \(D_0,D_1,N_d\) are the explicit pure-point polynomials listed in Appendix A and
\[
\boxed{\bar D_0=\bar D_1=\bar N_d=\bar J>0.}
\tag{13.0p}
\]
Therefore, whenever \(\lambda\ne0\),
\[
\frac d\lambda=2\zeta\,U_d,
\qquad
U_d:=\frac{JN_d}{D_0D_1}\in B_{\rm point}^{\times},
\qquad \bar U_d=1.
\tag{13.0q}
\]
Since \(0\le\zeta\le1/J\), the signed coefficient cannot occur strictly before \(\lambda\). If \(\lambda=0\), (13.0o) and the unit property force \(d=0\) exactly. Hence every coefficient layer strictly before \(\lambda\) consists only of literal positive EVENT packets. No support-side determinant or mixed ratio is used in this ordering.

On a genuinely coupled \(\lambda\)-row, normalize \(\lambda=1\) and put
\[
z=\varepsilon/c,
\qquad
P=J(1+z)-1,
\qquad
R=1-Jz,
\qquad
H=2Jz+J-1=P+Jz.
\tag{13.1}
\]
The complete contact equations give
\[
\boxed{
d=2z,
\qquad
e_3=\frac zP,
\qquad
e_5=\frac{Jz(1+z)}R,
\qquad
0<z<1/J.
}
\tag{13.2}
\]
The coefficient filtration gives the unit-or-child firewall:
\[
\boxed{
\overline z,\overline P,\overline R,\overline H>0
}
\tag{13.3}
\]
on the strict coupled row. Indeed:

- if \(\overline z=0\), the positive \(\lambda\) packet is earlier and the zero stratum is the already covered EVENT tangent;
- if \(\overline P=0\), the positive \(e_3\) packet is earlier;
- if \(\overline R=0\), the positive \(e_5\) packet is earlier.

Thus no positive infinitesimal is falsely treated as a residue unit.

## 13.2 Complete conormal Hessian

After the earlier literal positive packets are removed, the coupled \(2\times2\) conormal block is
\[
Q(j)=
\begin{pmatrix}Q_{11}&Q_{12}\\Q_{12}&Q_{22}\end{pmatrix},
\]
where
\[
Q_{11}=2(d+2e_3+2e_5+2),
\]
\[
Q_{12}=2(d-2e_3j+2e_3+2e_5),
\]
\[
Q_{22}=-2(2dj-d-2e_3j^2+4e_3j-2e_3-2e_5).
\tag{13.4}
\]
Substituting (13.2) gives the positive pure-point pivot
\[
\boxed{
A=Q_{11}
=
\frac{4(2Jz+J-1)}{(1-Jz)(Jz+J-1)}.
}
\tag{13.5}
\]
The determinant factors exactly as
\[
\boxed{
\det Q(j)=
\frac{16z(1+z)(j-J)^2}{(1-Jz)(Jz+J-1)}.
}
\tag{13.6}
\]
Therefore
\[
(\xi,\eta)Q(j)(\xi,\eta)^T
=
A\left(\xi+\frac{Q_{12}(j)}A\eta\right)^2
+
\alpha_j\bigl((j-J)\eta\bigr)^2,
\tag{13.7}
\]
where
\[
\boxed{
\alpha_j=\frac{4z(1+z)}{2Jz+J-1}>0.
}
\tag{13.8}
\]
Only the pure-point unit \(A\) is inverted.

For the independent reciprocal support copy write its two measurement scales as
\[
U_s=\rho^{-1},
\qquad
V_s=\gamma\rho,
\tag{13.8a}
\]
where \(\rho\) is the support radial parameter and \(\gamma\) is the remaining finite support parameter. On the exact kernel at \(j=J\), substituting (13.8a) into the **complete untruncated** cleared physical slack gives the even radial restriction
\[
F_{\ker}(J,\rho)
=
4Jz(1+z)(J+\gamma^2)\rho^2
+
\rho^4R_1(\rho^2).
\tag{13.9}
\]
Hence
\[
\boxed{
\alpha_\rho=4Jz(1+z)(J+\gamma^2)>0.
}
\tag{13.10}
\]
The complete cleared conormal scalar is even in \(\rho\) for every nearby \(j\), so
\[
\partial_\rho\partial_jF_{\rm con}(0)=0,
\qquad
\partial_\rho\partial_nF_{\rm con}(0)=0.
\]
The first-face completion also gives \(\partial_n\partial_jF_{\rm con}(0)=0\). The completed symbol below is used only to read the intrinsic Hessian. At the real contact let
\[
c_0=\overline{Q_{12}(J)/A}\in\mathbb R.
\]
Then the genuine one-copy local coordinate \(n_0=\xi+c_0\eta\) has the same first differential at the residue contact as the moving completion \(\xi+Q_{12}(j)\eta/A\). The two Hessians are related by an invertible real triangular congruence. Therefore it is legitimate to compute in the convenient exact triangular variables
\[
n=\xi+\frac{Q_{12}(j)}A\eta,
\qquad
q=j-J,
\qquad
s=\rho,
\]
the complete intrinsic quadratic class is
\[
\boxed{
\frac12\nabla^2F_{\rm con}(0)
=\operatorname{diag}(A,\alpha_j,\alpha_\rho),
}
\tag{13.11}
\]
and every diagonal entry has strictly positive residue by (13.3).

The variables \((n,q,s)\), together with the earlier actual EVENT normal classes, form a basis of the local contact cotangent space. They are used only to certify the intrinsic quadratic form. Lemma 3.1 then replaces them by genuine divided-difference roots of tensor rank at most two.

All exact rank drops or coordinate boundaries are coefficient/coordinate children already covered. Therefore
\[
\boxed{\operatorname{sosx}(B_\infty)\le2.}
\tag{13.12}
\]

---

# 14. Clean coverage and finite descent

There are only three measurement types:

1. both measurements finite;
2. exactly one measurement projective;
3. both measurements projective.

The preceding sections classify them as follows.

- Finite: REG, \(Q_0,Q_1,Q_2\), or a coefficient/coordinate child.
- One-projective: full-normal, one-missing, both-current smooth, neither-current Laurent, A17, or a proper child. Both-missing is impossible.
- Double-projective: exactly \(B_0\) or \(B_\infty\), plus their proper children.
- The noncoordinate divisors \(F_1,F_4,F_3\) are not new tensor problems: on their open loci they are exact EVENT gauges by Section 7. Their strict boundaries are coordinate children, terminal \(F_3\), or projective nodes already listed.

Nested coordinate restrictions commute. For the deeper determinant children, with \(H=D+2E_2\),
\[
\begin{aligned}
G|_{E_1=E_3=0}&=DH(D+2E_4),\\
G|_{E_1=E_4=0}&=H(D^2+2DE_3-4E_3E_5),\\
G|_{E_3=E_4=0}&=H(D^2-4E_1E_5),\\
G|_{E_1=E_3=E_4=0}&=D^2H.
\end{aligned}
\tag{14.1}
\]
Every pair restriction agrees with \(D^2H\) on the triple overlap.

## 14.1 Common-zero incidence lemma

The point and support copies are independent, but a zero-residue tensor must still be placed in one common one-copy contact node. This is a consequence of the displayed Bell geometry, not an extra atlas assumption.

**Lemma 14.1.** Let \(R/\mathbb R\) be real closed, let \(x\in(K_\star)_R\) be extreme, and let \(e,g\in E_R\) satisfy
\[
\overline{e(x)}=0,
\qquad
g(x)=0.
\]
Assume \(g\) is an active generator from the normal cone of \(x\). Then, after passing only through the exact coordinate restrictions and one-copy gauges already printed above, the residue support copies \(\bar e,\bar g\) lie either on one common strict contact node or on one common proper child. Their infinitesimal rates need not agree.

**Proof.** If \(g\) is a coordinate generator, its coordinate vanishes at \(x\), and every other nonnegative active support restricts to that same coordinate face. So suppose \(g\) is an evaluation support.

On a finite strict node, the logarithmic circuit Hessian is positive definite and the Bell slack has a simple kernel. Hence the measurement contact and projective kernel ray are unique; the normalized active evaluation support is unique. If one of the strict factors in the finite determinant table vanishes, the corresponding exact coefficient or amplitude is zero, so both active supports restrict before division to the named coefficient/coordinate child. The three circuit branches \(Q_0,Q_1,Q_2\) have unit-Jacobian support coordinates, so the same conclusion holds there.

On \(F_1,F_4,F_3\), Section 7 transports both copies by the same one-copy birational EVENT gauge. Thus common-node incidence is reduced exactly to the preceding EVENT/finite statement. A zero gauge denominator is, by definition, one of the printed strict boundary children; no quotient is taken there.

On a one-projective strict node, Proposition 9.1 gives the unique exact radial critical center and the actual normal packets have positive current coefficients. The full-normal, one-missing, both-current and Laurent children are therefore unique strict contacts. Exactly-one-current is A17, whose two unit-Jacobian divisions give its common terminal graph. Loss of a required unit is an exact child before division.

It remains only to check the simultaneous-projective base. The required routing is read from three signed-cofactor formulas. At \(B_0\), after the substitution of Section 12 and removal of the common positive power,
\[
\left.r^{-12}X_{B_0}\right|_{r=0}
=C_0(0,0,0,A,1-A-LT^2,0)^T,
\tag{14.2a}
\]
where
\[
C_0=2^{13}T^9(1+T^2)^3(LT^2-1)^3(LT^2+L-1).
\tag{14.2a'}
\]
On the strict \(B_0\) chamber this has fixed nonzero sign. The exact factor \(LT^2=1\) is processed before division and gives the \(e_1\) coordinate child; the factor \(LT^2+L-1=0\) is the \(h=0\) boundary closed in Section 12. Hence the unresolved exceptional support lies in \(\operatorname{cone}\{e_3,e_4\}\). If the remaining incidence factor vanishes, the next cofactor layer is
\[
\left.r^{-16}X_{B_0}\right|_{r=0}
=2^{13}AT^9(0,0,\zeta,A,\zeta,0)^T,
\qquad
\zeta=\frac{1-A(1+T^2)}{1+T^2},
\tag{14.2b}
\]
so it lies in \(\operatorname{cone}\{e_2,e_3,e_4\}\).

At \(B_\infty\), with the reciprocal endpoint substitution,
\[
\left.r^7X_{B_\infty}\right|_{r=0}
=
\frac{2^{13}J^2c^9}{(J+c^2)^7}
\begin{pmatrix}
2(J-1)(J-K)\\[1mm]0\\[1mm]K(J-1)\\[1mm]J-K\\[1mm]0\\[1mm]J(J-1)(J-K)
\end{pmatrix}.
\tag{14.2c}
\]
Thus on the strict reciprocal chamber \(J>1,0\le K<J\), the \(e_1\)- and \(e_4\)-coordinates vanish and the endpoint lies in the already named \(F_1\cap F_4\) child. No coordinate ray is inferred from the leading vector alone. Exact \(J=1\) is routed by (13.0b1)--(13.0b4) to \(e_2\) or \(e_3\); exact \(K=J\) is routed by (13.0b5)--(13.0b7) either to \(e_3\) or to a sign-incompatible row. The exact boundary \(c=0\) was routed to the one-projective chart before reciprocal division, while \(K=0,J>1\) remains in the strict reciprocal analysis.

If two residue supports specialize to different base components, (14.2a) and (14.2c) show that their first nonzero positive coordinate packets belong to the displayed coordinate faces. Since both nonnegative supports vanish at the same residue point, every coordinate carrying positive residue in either packet vanishes there. The union of those exact coordinate equations defines a proper common coordinate child. Coefficients with zero residue are not discarded: the finite coefficient filtration is repeated inside that child until either a positive coordinate packet appears or an exact zero coefficient lowers \(N_{
m act}\). If the two supports specialize to the same base component, the strict endpoint theorem or one of its exact zero-coefficient exits applies.

These cases exhaust finite, one-projective and double-projective measurement geometry. No inference of the form “different specializations imply positive residue” has been used. \(\square\)

## 14.2 The well-founded invariant

Blow-ups and the exact gauge maps of Section 7 are coordinate changes, not descent arrows. Every genuine nonterminal child does one of two things:

1. it sets at least one of the nine original nonnegative Bell coefficients to zero; or
2. without restoring any deleted coefficient, it imposes an additional exact coordinate/projective/attachment equation, hence passes to a proper lower-dimensional stratum.

Define
\[
\boxed{
\kappa=(N_{\rm act},d_{\rm stratum})
}
\tag{14.3}
\]
with lexicographic order, where \(N_{\rm act}\in\{0,\ldots,9\}\) is the number of original Bell coefficients not yet deleted and \(d_{\rm stratum}\) is the real dimension of the current exact one-copy stratum after all already imposed equations.

A coefficient child strictly lowers \(N_{\rm act}\). At fixed \(N_{\rm act}\), a proper coordinate/projective/attachment child strictly lowers \(d_{\rm stratum}\). No deleted coefficient is ever reintroduced, and no imposed exact stratum equation is removed along a child arrow. Hence every nonterminal child strictly lowers \(\kappa\) in the finite set
\[
\{0,\ldots,9\}\times\{0,\ldots,d_{\max}\}.
\]
Therefore no infinite child chain or cycle exists.

This is the entire termination argument. No valuation fan, Newton depth, or projective-depth convention is needed.

---

# 15. Universal rank-two star theorem

Fix a real closed extension \(R/\mathbb R\), a generator \(e\in E_R\), and a point \(x\in\operatorname{Ext}(K_\star)_R\).

If the multiplication residue of \(e^{\otimes}(x)\) is strictly positive, Scheiderer's positive-residue proposition gives \(\operatorname{sosx}=1\).

Suppose the residue is zero. Because \(x\) is extreme, it lies on the boundary and its normal cone contains a nonzero affine support \(f\) with \(f(x)=0\). The supporting-hyperplane assertion transfers to \(R\): for the fixed compact semialgebraic family \(K_\star\), existence of a nonzero affine functional nonnegative on \(K_\star\) and vanishing at a boundary point is a first-order semialgebraic statement. Write \(f\) as a finite nonnegative conic combination of generators from \(E_R\). Since every generator is nonnegative at \(x\), at least one generator \(g\) occurring with positive coefficient satisfies \(g(x)=0\). If \(g\) is a coordinate support, we are already on its exact coordinate child. Otherwise \(g\) is an evaluation support and parametrizes the point copy through the contact equations.

Apply Lemma 14.1 to the residue supports \(\bar e,\bar g\). It gives either one common strict contact node or one common proper child. Lift the two copies independently into that common one-copy chart; only the residue node is shared, while all point/support infinitesimal rates remain independent. The tested support \(e\) may have a positive infinitesimal exact critical value even though its residue is zero. If so, split that value off by Lemma 3.4. The centered loss then enters one of the finite/projective types in Section 14.

At REG and every strict smooth node, the centered loss lies in \(J^2\) and has positive intrinsic quadratic class, so Lemma 3.1 gives rank-two roots. The exceptional nodes are closed explicitly:

- \(Q_0,Q_1,Q_2\): Section 6;
- \(F_1,F_4,F_3\): exact EVENT transfer, Section 7;
- terminal \(F_3\): Section 8;
- one-projective singularities: Section 9;
- A17: Section 10;
- \(B_0\): Section 12;
- \(B_\infty\): Section 13;
- coordinate/polyhedral/rotated-SOC children: literal rank-one/rank-two squares.

If a proper child is reached, induction on \(\kappa\) applies. Section 14 proves that the induction terminates.

Therefore
\[
\boxed{
\operatorname{sosx}(e^{\otimes}(x))\le2
\quad
\text{for all }R,e,x\text{ occurring in (2.1).}
}
\tag{15.1}
\]
Scheiderer's quantitative criterion and the support-cone theorem now yield
\[
\boxed{\operatorname{sxdeg}(K_\star)\le2.}
\tag{15.2}
\]

---

# 16. Global Bell assembly

All eight stars are affine Bell relabelings of \(K_\star\), so each has semidefinite extension degree at most two. The central body \(G\) is polyhedral. Scheiderer's finite convex-hull bound gives
\[
\operatorname{sxdeg}(\mathcal B_Q)\le2.
\tag{16.1}
\]
His polar bound then gives
\[
\boxed{\operatorname{sxdeg}(\mathcal Q_{2,2,2})\le2.}
\tag{16.2}
\]

---

# 17. Independent lower bound

Restrict the correlator projection to
\[
C_{00}=C_{01}=C_{10}=t=\cos\alpha,
\qquad
0\le\alpha\le\frac\pi4.
\tag{17.1}
\]
Tsirelson's vector representation supplies real unit vectors with
\[
C_{xy}=a_x\cdot b_y.
\]
Hence
\[
d(a_0,b_0)=d(a_0,b_1)=d(a_1,b_0)=\alpha.
\]
The spherical triangle inequality gives
\[
d(a_1,b_1)\le3\alpha,
\]
so
\[
C_{11}\ge\cos3\alpha.
\tag{17.2}
\]
Equality is attained by coplanar vectors at angles
\[
0,\quad\alpha,\quad-\alpha,\quad2\alpha.
\]
Therefore the exact boundary contains
\[
\boxed{
C_{11}=\cos(3\arccos t)=4t^3-3t,
\qquad
\frac1{\sqrt2}\le t\le1.
}
\tag{17.3}
\]
This is a genuinely curved cubic arc. A linear projection of a polytope is a polytope, so \(\mathcal Q_{2,2,2}\) is not polyhedral. Since semidefinite extension degree one is exactly the polyhedral case,
\[
\operatorname{sxdeg}(\mathcal Q_{2,2,2})\ge2.
\tag{17.4}
\]
Combining (16.2) and (17.4),
\[
\boxed{
\operatorname{sxdeg}(\mathcal Q_{2,2,2})=2.
}
\qquad\square
\]

---

# 18. Dependency firewall

The proof above does **not** use any of the following historical moves:

- “different specializations imply positive residue”;
- positivity in a non-Archimedean field as a substitute for positive residue;
- a weighted or mixed Bell Gram generator as an SOS root;
- the generic filtered-collar/R2.4 implication;
- identification of the independent point and support tensor copies;
- a moving rank-six reciprocal Schur root;
- setting the A17 positive infinitesimal critical value to zero;
- a machine sparse replay as the proof of A17 \(J^2\)-membership;
- a Newton/Pareto fan as the coverage theorem;
- a fixed production atlas as a hypothesis of Scheiderer's pointwise criterion.

Historical computer algebra used during discovery or independent checking is treated as nonlogical. Every identity used as a premise above is finite and displayed in this document, or is one of the explicitly cited external theorems.

---

# 19. Reader's verification checklist

A hostile human audit can proceed in the following order.

1. Check the Bell star decomposition using Fine, the visible facets, and midpoint identity.
2. Check the support-cone separation argument.
3. Check Lemma 3.1 line by line. In particular, note that the \(J^3\) remainder appears **only after** the residual quadratic class has been represented modulo \(J^3\).
4. Verify the EVENT row minors (5.1)-(5.5).
5. Verify the four discriminants (5.7) and the empty-mask factorizations (5.8).
6. Verify the \(Q_1,Q_2\) syzygies, actual scale \(\rho\), packet completion, and Laurent squares (6.2)-(6.15).
7. Verify the one-line matrix identity (7.2), its determinant (7.3), and the three gauge factors (7.6)-(7.8).
8. Verify terminal \(F_3\) positivity and full-slack transfer from (8.2)-(8.8).
9. Verify the one-projective entry theorem (9.3)-(9.10), then the actual normal packets and scalar alternatives.
10. Verify the A17 terminal-entry divisions (10.2)-(10.13), the finite cancellation lists (10.20)-(10.23), residue matrix (10.24), fifth pivots (10.25), and five-simplex divided-difference completion (10.28)-(10.32).
11. Verify the two-component projective base locus (11.1)-(11.2).
12. Verify the \(B_0\) chamber from the leading cofactor formulas (12.0b)-(12.0f), including the exact \(LT^2=1\to e_1\) exit; then expand the four squares in (12.1), check the transverse scale (12.0a), and verify the \(h=0\) exit (12.2)-(12.3).
13. Verify the reciprocal chart, the exact point-boundary routing (13.0b0)-(13.0b9), and the transverse endpoint closure (13.0c)-(13.0n); then check the exact Cramer firewall (13.0o)-(13.0q) using Appendix A, the determinant factorization (13.6), and the even radial Hessian (13.9)-(13.11).
14. Verify the common-zero incidence lemma, especially the endpoint cofactor routings (14.2a)-(14.2c), and then check that every genuine child lowers \((N_{\rm act},d_{\rm stratum})\).
15. Apply Scheiderer's criterion and closure theorems.
16. Check the cubic lower-bound slice.

This checklist contains no request to run a program.

---

# 20. Self-containment and provenance note

The proof above is intentionally written so that the *logical* argument can be checked from finite displayed algebra. The finite mask table, the explicit empty-mask polynomials \(Q_\pm\), the projective entry envelope, and the two A17 unit-Jacobian divisions have been reproduced here rather than cited as machine output. Historical project computations may be used as independent discovery or cross-checking records, but none is a premise of the argument above. A reader may verify every Bell-specific identity by ordinary determinant expansion, substitution, completion of squares, and differentiation of the displayed nine-channel slack.



# Appendix A. Exact reciprocal Cramer factors

This appendix makes the signed/\(\lambda\) ordering in (13.0o) directly checkable. Put
\[
x=r^2,
\qquad z=\zeta.
\]
The first factor is
\[
D_0=J+x(Jc^2z+Jc^2-Kc^2)-Kc^4x^2.
\tag{A.1}
\]
Write
\[
D_1=J+xD_{1,1}+x^2D_{1,2}+x^3D_{1,3}+x^4D_{1,4},
\tag{A.2}
\]
with
\[
D_{1,1}=J^2+Jc^2+J-Kc^2+c^2,
\tag{A.3}
\]
\[
D_{1,2}=J^2c^2z+J^2c^2+J^2-JKc^2z-2JKc^2
+Jc^2z+2Jc^2-Kc^4+c^4,
\tag{A.4}
\]
\[
D_{1,3}=c^2\bigl(J^2z+J^2-JKc^2z-2JKc^2-JK
+Jc^2z+Jc^2+K^2c^2-Kc^2\bigr),
\tag{A.5}
\]
\[
D_{1,4}=Kc^4(-J+Kc^2-c^2).
\tag{A.6}
\]
Likewise
\[
N_d=J+xN_1+x^2N_2+x^3N_3+x^4N_4+x^5N_5+x^6N_6,
\tag{A.7}
\]
where
\[
N_1=J^2z+2J^2+Jc^2z+2Jc^2+J-Kc^2+c^2,
\tag{A.8}
\]
\[
\begin{aligned}
N_2={}&J^2c^2z^2+4J^2c^2z+3J^2c^2+J^2z+2J^2
-JKc^2z-3JKc^2\\
&+Jc^4z+Jc^4+2Jc^2z+3Jc^2-2Kc^4+c^4,
\end{aligned}
\tag{A.9}
\]
\[
\begin{aligned}
N_3=c^2\bigl(&J^2c^2z^2+2J^2c^2z+J^2c^2+J^2z^2+4J^2z+3J^2
-3JKc^2z\\
&-5JKc^2-2JK+2Jc^2z+2Jc^2+K^2c^2-Kc^4-2Kc^2\bigr),
\end{aligned}
\tag{A.10}
\]
\[
N_4=c^4\bigl(J^2z^2+2J^2z+J^2-2JKc^2z-2JKc^2
-2JKz-4JK+2K^2c^2-2Kc^2\bigr),
\tag{A.11}
\]
\[
N_5=-Kc^6(2Jz+2J-Kc^2-K),
\qquad
N_6=K^2c^8.
\tag{A.12}
\]
Substitution in the complete seven-channel value/first-jet equations and ordinary fraction-free elimination give
\[
D_0D_1d-2JzN_d\lambda=0.
\tag{A.13}
\]
No division is required to check (A.13): after inserting the canonical slack (4.1), both sides are polynomials in the reciprocal point variables and equality follows by direct expansion. Since \(x=r^2\in\mathfrak m\), all nonconstant terms in (A.1), (A.2), and (A.7) have zero residue, and hence
\[
\bar D_0=\bar D_1=\bar N_d=\bar J>0.
\tag{A.14}
\]
This proves the pure-point unit claim used in Section 13 without any computational certificate.


# References

1. A. Fine, *Hidden Variables, Joint Probability, and the Bell Inequalities*, Phys. Rev. Lett. **48**, 291-295 (1982).
2. L. Masanes, *Extremal quantum correlations for N parties with two dichotomic observables per site*, arXiv:quant-ph/0512100.
3. C. Scheiderer, *Second-order cone representation for convex sets in the plane*, Math. Programming / arXiv:2004.04196. In particular Proposition 1.6, Corollary 1.9, Proposition 3.5, and Theorem 3.10.
4. C. Scheiderer, *Smooth hyperbolicity cones are second-order cone representable*, arXiv:2509.17121 (2025). The proof of Lemma 3.12 supplies the exact diagonal-ideal absorption used in Lemma 3.1 above.
