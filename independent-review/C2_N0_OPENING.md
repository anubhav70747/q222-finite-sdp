# Exact C2 wall-to-strict opening in a fixed half-power space

27 October 2026. This note supplies a local PSD completion for the fractional carrier introduced in the [2039 package](../complete-proof-materials/materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-2039/task3/package/PROOF_AND_OPEN_GATE.md). Every C2 n=0 wall point with a positive definite Gram on all four reduced blocks of the larger 27-feature wall support has an explicitly described neighborhood of strict n>0 admission in **28 fixed half-power features**. Multiplication by 1+v² embeds the result in **49 fixed K364 features**.

The theorem does not assume that the selected carrier can be subtracted globally. It proves a PSD residual near the stated wall points by a fixed coefficient right inverse and a quantitative positive-definiteness bound. Whole-wall PSD selection, arbitrary strict C2 admission H2, and C1 admission H1 remain open. The original derivation and replay were performed together; the subsequent independent checks are linked below.

Subsequent independent [analytical review](../complete-proof-materials/materials/historical/analysis/c2_n0_opening_analytic_verification_2026_10_27.txt) and [bounded exact replay](../complete-proof-materials/materials/historical/analysis/c2_n0_opening_exact_verification_2026_10_27.txt) support the stated conditional neighborhood theorem. The replay verifies the already archived positive step 2^-24 without repeating the original source-search loop.

## Source, contact and quantifiers

Use the complete SG2 C2 parent F_theta, d=1, lambda=beta=0, with its actual native parameters theta=(a,b,c,t,r), and all six native weights n,m,C,A,D,alpha. All tests u,v and all four state coordinates are independent. At a wall point theta0, n=0; all five retained weights and the other native gaps are strictly positive. Require c,a,b,r,h=b-a²t>0 and the genuine stationary contact

\[
(u_c,v_c;\psi_c)=(c,ac;b/a,-cr,ch/a,1).
\]

Every component of this contact is nonzero. Normalize using pure-source scalars:

\[
u=cU,\quad v=acV,\quad
\psi=\operatorname{diag}(b/a,-cr,ch/a,1)\Phi,
\qquad\Phi=(W,X,Y,Z)^T.
\]

The moving contact is now U=V=1, Phi=(1,1,1,1). These rescalings preserve PSD and the fixed half-power dictionary. They introduce no testing denominators.

Define the four wall feature blocks

| Block | Fixed normalized features |
|---|---|
| 0 | VW, UW, UVW, X, VX, UVX, Y, VY, UY, UVY, Z, VZ, UZ |
| 1 | U sqrt(V) W, sqrt(V) X, sqrt(V) Y, U sqrt(V) Y, sqrt(V) Z |
| 2 | sqrt(U) V W, sqrt(U) V X, sqrt(U) Y, sqrt(U) V Y, sqrt(U) Z |
| 3 | sqrt(UV) W, sqrt(UV) X, sqrt(UV) Y, sqrt(UV) Z |

Every feature equals one at the normalized contact. In each block replace the features by differences from its first feature. Denote these difference columns by rho_b. Their lengths are (12,4,4,3).

**Wall hypothesis.** At theta0 there are four matrices Q_b^0 strictly positive definite such that

\[
\widehat F_{\theta_0}=\sum_{b=0}^3\rho_b^TQ_b^0\rho_b.
\tag{1}
\]

This is positivity on the entire contact complement in every one of the four larger wall blocks. A certificate positive only on the smaller 18-feature support does not automatically satisfy this hypothesis. The archived larger-wall checker supplies eight instances with the required positivity; the standalone checker below reconstructs one of them from its rational proposal.

## A fixed linear map exhausts the residual constraints

Let M map the 104 independent symmetric entries of the four Q_b to the coefficients of sum rho_b^T Q_b rho_b. Off-diagonal entries contribute twice. Products have 72 distinct keys (i,j,l,k), representing U^i V^j Phi_l Phi_k with l<=k. Let J be the six-row map that evaluates the derivatives in U,V,W,X,Y,Z at the normalized contact.

The evaluated rational data prove

\[
M\in\mathbb Q^{72\times104},\qquad
\operatorname{rank}M=66,\qquad
JM=0,\qquad\operatorname{rank}J=6.
\tag{2}
\]

In particular **im M = ker J**. A specified 66-by-66 minor of M has determinant 2^54, and a specified six-by-six minor of J has determinant 2. These two minors and JM=0 establish the ranks without an unspecified generic specialization. All source parameters have been removed from these maps.

The selected M minor uses only the first 78 columns, belonging to the integer block. Its inverse defines a constant rational E in Q^(104 by 72), supported in those first 78 rows, with

\[
MEM=M,\qquad ME r=r\quad(r\in\ker J).
\tag{3}
\]

Thus **only the integer block needs a coefficient correction**. No exceptional native denominator or parameter-dependent pivot is required. The [fixed data](checks/gram-lift/ends/n0-opening/fixed_map.json) give every feature, coefficient key, upper-entry position, pivot row/column and inverse entry. The standalone checker reconstructs all of them from feature products.

## The contact-compatible carrier removes the missing diagonal

Choose the source-rational q=a/(1+a), and write

\[
\begin{split}
\widehat C_\theta=4nc^4r^2\big[&(-UX+qVX+(1-q)VW)^2\\
&+2qUV(W-X)^2\big].
\end{split}\tag{4}
\]

This is the two-square carrier before multiplication by 1+v². It is PSD for U,V>0 and arbitrary Phi, vanishes at the moving contact, and has all six contact derivatives zero. At n=0 it vanishes identically. It is a positive-quadrant half-power SOS, not a separate Bell source.

Let

\[
r_\theta=\operatorname{coeff}(\widehat F_\theta-\widehat C_\theta).
\]

The following statements are coefficientwise identities, checked with the six weights and a,b,c,h,r,q formal:

1. The carrier removes the forced U²X² coefficient exactly.
2. Its UVX² term cancels internally.
3. Every remaining coefficient belongs to the same 72-key support of M, including the required negative mixed WX couplings.

For actual native stationary sources, the six full-parent derivatives vanish. Since the carrier also has zero contact jets,

\[Jr_\theta=0.\tag{5}\]

Equation (2) therefore solves **all** residual coefficient conditions, including state normals and mixed angular/state terms. Positive raw diagonal budgets are not substituted for this calculation.

## Explicit PSD neighborhood theorem

Write x0 for the 104 upper entries of the wall matrices in (1), and r0=Mx0. Set

\[
x(\theta)=x_0+E(r_\theta-r_0),
\qquad Q_b(\theta)=\operatorname{unpack}_b x(\theta).
\tag{6}
\]

For every actual native source near theta0, (3) and (5) give the exact full identity

\[
\widehat F_\theta=\widehat C_\theta+
\sum_{b=0}^3\rho_b^TQ_b(\theta)\rho_b.
\tag{7}
\]

This identity is valid before any positivity argument. All coefficients depend only on the source and the fixed wall certificate. The last three matrices remain exactly Q_b^0 in normalized coordinates.

For an explicit sufficient positivity bound put

\[
\mu_b=\|(Q_b^0)^{-1}\|_\infty^{-1}>0,
\qquad C_b=\max_{j\text{ in block }b}\sum_k|E_{jk}|.
\]

The infinity norm here is the maximum absolute row sum. Since the inverse is symmetric, its largest eigenvalue is at most this norm, hence Q_b^0 >= mu_b I. The stored right inverse has

\[(C_0,C_1,C_2,C_3)=(141,0,0,0).\]

If epsilon=||r_theta-r0||_infinity, every correction entry in block b has magnitude at most C_b epsilon. Consequently its spectral norm is at most d_b C_b epsilon, where d_b=(12,4,4,3). In particular the explicit condition

\[
\|r_\theta-r_0\|_\infty<
\varepsilon_0:=\frac{\mu_0}{2\cdot12\cdot141}
\tag{8}
\]

implies Q_0(theta) >= mu_0 I/2; the other three blocks retain their wall positivity. The native coefficient formulas and all normalization factors are continuous at theta0, so (8) holds throughout a relative neighborhood of that wall point. The stationary-source event-opening argument supplies nearby physical strict sources; the exact rational instance below also exhibits one directly.

For every source in this neighborhood with n>=0, both terms on the right of (7) are PSD. The wall features together with the single new integer feature UX contain every root in (7), so the **same 28-feature half-power space** works throughout the neighborhood. Source-dependent contact shifts and rescalings change only Gram coefficients.

Multiplying (7) by 1+v² uses the original roots and their v multiples. Their union contains 49 distinct fixed features, all with Alice exponent at most 1 and Bob exponent at most 2, in half-integer steps. They embed directly in Z_(6,12), giving K364 admission. The original SG2 sign congruences and homogeneous continuity extend full-parent positivity to the other testing quadrants and projective endpoints; the carrier alone need not have those symmetries.

This proves the neighborhood theorem. It also supplies one fixed 28-feature cone for the union of these neighborhoods over all wall points satisfying (1). It supplies neither a common metric radius at singular wall limits nor PSD wall certificates where (1) is unproved.

## Exact delivered instance

The reconstructed baseline has

\[
\theta_0=(a,b,c,t,r)=
(1,22/13,1,2581/2210,2/5).
\]

The checker repairs the archived rational proposal exactly using E, verifies all four wall reduced blocks by rational LDL, and computes a positive rational epsilon0. It proves epsilon0 > 9/1000000 exactly, so coefficient distance less than 9/1000000 is an explicit sufficient neighborhood condition at this seed. It then certifies

\[
\theta_1=(1,22/13,1,2581/2210+2^{-24},2/5),
\qquad n(\theta_1)=221/19042140160>0.
\]

All native gaps and five other weights are strictly positive. All coefficient equations and every reduced LDL pivot are checked exactly. The rational Gram vectors, weights, eigenvalue lower bounds and coefficient radius are in [exact_results.json](checks/gram-lift/ends/n0-opening/exact_results.json).

Run from any directory with SymPy installed:

```bash
python independent-review/checks/gram-lift/ends/n0-opening/verify_n0_opening.py
```

The checker uses no solver, sample sweep, network or unavailable runtime object. Its [seed proposal](checks/gram-lift/ends/n0-opening/seed_wall_proposal.json) is only a starting vector; exact coefficient matching and positivity are rechecked.

## Remaining scope

The local coupled completion is now proved wherever the larger wall Gram is positive definite on its contact complement. It is not enough to know only that a wall source has some singular Gram. Universal PSD selection on the whole n=0 wall is still open, as is admission away from these neighborhoods and the previously verified strict open c-end regions. The first-u-jet/frozen-top obstruction and the global half-power counterexample retain their scopes. H1/H2, and therefore unconditional exactness of the evaluated global SDP, remain open.

The physical point is that an independent n-event square cannot be added to a tight source: it fails to vanish at the shifted contact. The contact-zero carrier supplies the compulsory new diagonal without spending that zero. After it is extracted, stationarity is exactly the remaining affine compatibility condition, and the strict wall Gram absorbs a sufficiently small coefficient correction.
