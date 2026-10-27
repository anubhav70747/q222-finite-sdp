# Sun–Wang–Yu: exact matching to the C1/C2 contact geometry

27 October 2026. Base: `1b293832ffeff46ebbee834b89554481d70f544e`.

**The matching is complete for every real canonical quantum-valid source with d>0 and a finite noncommuting contact.** This includes both native families and their permitted zero-weight faces. It gives one short alternative proof of the positive contact Hessian already established in SG18–SG19. It does not prove universal fixed-library admission or change the open H1/H2 status.

The external input is Liang-Liang Sun, Jie Wang and Sixia Yu, [*Exact Identification of Quantum Correlations in the Minimal Bell Scenario*, arXiv:2610.34364v1](https://arxiv.org/html/2610.34364v1), submitted 25 October 2026. We use their Proposition 4 and Lemma 6. The Hessian consequence below is our deduction; it is not a theorem quoted from that paper. This is a targeted comparison and hypothesis audit, not an independent certification of their entire paper.

## 1. Match the affine functional, including its constant

Use the source and event ordering of [SG1–SG2](../latex/chapters/04_geometry.tex):

\[
f=d\left[1-\frac{E_{00}-E_{01}+E_{10}+E_{11}}2\right]
  +4\sum_i a_i p_i,
\qquad (a_i)=(\lambda,n,m,C,A,D,\alpha,\beta).
\]

Write

\[
f=\kappa-\mathcal B,\qquad
\mathcal B=\sum_x a_x^{B}A_x+\sum_y b_y^{B}B_y+\sum_{xy}c_{xy}E_{xy}.
\]

The exact coefficient map is

\[
\begin{aligned}
\kappa&=d+\lambda+n+m+C+A+D+\alpha+\beta,\\
(a_0^B,a_1^B)&=(A+\alpha-\lambda-n,\ C+D-m-\beta),\\
(b_0^B,b_1^B)&=(n+\beta-C-A,\ m+\alpha-\lambda-D),\\
(c_{00},c_{01},c_{10},c_{11})
&=\left(\frac d2+n+A,-\frac d2-\lambda-\alpha,
         \frac d2+C+\beta,\frac d2+m+D\right).
\end{aligned} \tag{1}
\]

For d>0 and nonnegative event weights, the four correlator coefficients have signs (+,−,+,+). In particular

\[
c_{00}c_{01}c_{10}c_{11}
=-\left(\frac d2+n+A\right)\left(\frac d2+\lambda+\alpha\right)
  \left(\frac d2+C+\beta\right)\left(\frac d2+m+D\right)<0. \tag{2}
\]

Thus frustration follows from the canonical orthant alone. No assumption about strictly positive contact probabilities is needed. If f is quantum-valid and has a quantum zero, then κ is exactly the quantum maximum of B, not merely an upper bound.

The [native entrance domains EN1](../latex/chapters/05_entry.tex) match as follows. In each row d>0 and the contact is finite and noncommuting.

| Family | Original coefficient hypotheses | Consequence for (1) |
|---|---|---|
| C1 | λ=n=0; m,C,A,α>0; D,β≥0; D+β>0 | (+,−,+,+), including D=0 or β=0 |
| C2 | λ=β=0; m,C,A,D,α>0; n≥0 | (+,−,+,+), including n=0 |

This uses the established native entrance and stationary-validity results, not arbitrary rational parameters outside their domains. The oriented C1 restriction 0<k≤1 is unnecessary for this contact theorem. More strongly, both families have λ=0: the all-positive deterministic assignment has f=4λ=0. Physicality and the noncommuting contact therefore give L_B=Q_B=κ. The strict-quantum-advantage theorem is inapplicable to these native sources; the frustrated-case lemmas are the required inputs.

## 2. Match the measurements and cleared state matrix

Let Z and X be the real Pauli matrices. The actual gauge of SG2 is

\[
A_0=B_0=Z,\quad
A_1=\frac{(1-u^2)Z+2uX}{1+u^2},\quad
B_1=\frac{(1-v^2)Z+2vX}{1+v^2}.
\]

For u,v>0 its overlaps are

\[
s=\frac{1-u^2}{1+u^2},\qquad t=\frac{1-v^2}{1+v^2},\qquad
J=\frac{\partial(s,t)}{\partial(u,v)}
=\operatorname{diag}\left(-\frac{4u}{(1+u^2)^2},-\frac{4v}{(1+v^2)^2}\right). \tag{3}
\]

They lie in (−1,1), and J is invertible. The off-diagonal measurement entries equal the positive square roots used in the external canonical gauge. Put Δ=(1+u²)(1+v²), M=κI−B(s,t), and let P be the state matrix of SG2. Then, identically in all source coefficients and all testing/state variables,

\[
P(u,v)=\Delta M(s(u),t(v)),\qquad
F(u,v,\psi)=\psi^TP\psi=\Delta\|\psi\|^2 f(q). \tag{4}
\]

There is no source-dependent multiplier or state substitution hidden in (4). Negative finite angle gauges reduce to the positive ones by local diagonal state congruences: conjugating with diag(1,−1) on Alice flips u, and on Bob flips v. These invertible coordinate changes preserve contact Hessian definiteness.

## 3. A short polynomial lemma

**Lemma.** Let h(s,t) have separate degree at most two and be nonnegative near an interior zero z. Suppose its s²t² coefficient is negative and it has no horizontal or vertical zero line through z. Then its Hessian at z is positive definite.

**Proof.** At the minimum, the gradient vanishes and the Hessian is positive semidefinite. For a nonzero null direction w, the restriction h(z+εw) has zero constant, linear and quadratic terms. Nonnegativity for both signs of ε forces its cubic coefficient to vanish. If both components of w are nonzero, its quartic coefficient is [s²t²]h times w₁²w₂²<0, a contradiction. If one component is zero, the restriction has degree at most two and vanishes identically, giving a forbidden axis-parallel zero line. ∎

For δ(s,t)=det M(s,t), global validity gives M≥0 and δ≥0 on the overlap square. Sun–Wang–Yu Proposition 4 gives separate degree at most two and [s²t²]δ=16c₀₀c₀₁c₁₀c₁₁. Their Lemma 6 excludes axis-parallel zero lines when all four coefficients are nonzero. Equations (1)–(2) therefore satisfy every hypothesis of the lemma, including at tied local and quantum bounds. Hence

\[
D^2_{s,t}\delta(s_c,t_c)\succ0. \tag{5}
\]

Global validity is an input here. A stationary point of an arbitrary nonphysical source is not sufficient. In the native application, the earlier SG7 stationary-validity argument supplies this input without using SG18–SG19.

## 4. Match the full five-variable Hessian exactly

The native SG3 already proves that the Z=0 principal block H₃ of P is positive definite for d>0 and uv≠0. Thus at a contact P has a simple kernel with Z≠0. Rescale the state to Z=1, write z=(W,X,Y), and partition P into (H₃,b;bᵀ,c). The exact completion is

\[
F=(z+H_3^{-1}b)^TH_3(z+H_3^{-1}b)+R(u,v),\qquad
R=\frac{\det P}{\det H_3}
  =\frac{\Delta^4}{\det H_3}\,\delta(s(u),t(v)). \tag{6}
\]

At the contact δ=0 and Dδ=0. Differentiating (6) therefore gives

\[
D^2_{u,v}R\big|_c
=\frac{\Delta_c^4}{\det H_{3,c}}\,
 J_c^T D^2_{s,t}\delta\big|_c J_c\succ0. \tag{7}
\]

The local change z↦z+H₃⁻¹b has invertible Jacobian. At its zero the full Hessian is congruent to diag(2H₃,D²R). Both blocks are positive definite. This proves

\[
D^2_{W,X,Y,u,v}F\big|_{Z=1,\,c}\succ0. \tag{8}
\]

The power Δ⁴ in (6)–(7) is required: δ is a four-dimensional determinant. A factor Δ would be incorrect. The homogeneous six-state/angle formulation has a radial state null direction; (8) concerns the five-dimensional projective chart, not an unrestricted homogeneous Hessian.

## 5. What can be shortened, and what must remain

Equations (1)–(8) provide one alternative to the SG18–SG19 contact-curvature derivation for both native families, including their permitted third omissions. They avoid a separate eigenvalue perturbation argument and do not need the external state-simplicity lemma.

They also recover finite-contact coercivity without the original separate support exclusions. In the SG5 notation,

\[
\delta(s(u),t(v))=\frac{16u^2v^2}{D_u^2D_v^2}L(u^2,v^2).
\]

The multiplier is positive at contact, and the changes between positive (u,v), overlaps (s,t), and log squared angles are locally invertible. Since L and its first derivatives vanish there, (5) gives a positive definite logarithmic Hessian of L. By SG5–SG6 that Hessian is ∑ₑ wₑ eeᵀ, where wₑ=qₑ U_c^{e₁}V_c^{e₂}>0 on the active exponents, and stationarity gives ∑ₑ wₑ e=0. Positive definiteness makes the active exponents span the plane. Their positive balance then places zero in the interior of their convex hull: a supporting functional through zero would vanish on every active exponent, contradicting that span. This is exactly the finite-contact coercivity conclusion preceding SG18. It also makes the Laurent function strictly convex in log squared angles and gives uniqueness of its positive contact.

Thus this is a short alternative proof of finite-contact coercivity, uniqueness and full curvature. The ensuing source classification and deletion arguments still need those conclusions and their own remaining steps. The existing manuscript proof is preserved. This note does not remove it, substitute a new manuscript PDF, or replace the singular C1/C2 tensor-SOS cases. Over the reals, (8) is an algebraic implication with global PSD quantified over the testing variables. Clearing positive denominators expresses that implication in the language of ordered fields, so its real-closed-field transfer is legitimate. This remains a single-field Hessian statement; it is not a separated source/testing tensor-SOS certificate. Logarithms above abbreviate U∂_U and V∂_V; no logarithmic function is required for that transferred Hessian statement.

Curvature is source-wise. Its least eigenvalue and the size of a valid neighborhood need not have a positive lower bound as the source approaches a wall. No uniform degree, full-testing Gram, source-affine certificate, or K364 membership follows from (8).

## 6. Related work and the remaining frontier

Sun–Wang–Yu Theorem 1 gives a one-real-two-qubit-plus-local membership characterization. Its nonlinear feasibility formulation does not supply a fixed affine SDP. For affine f, validity on that decomposition is still equivalent to validity on all real two-qubit behaviors: local behaviors are already in their convex hull. The polar positivity obligation therefore remains.

The distinct construction target is one fixed finite affine semidefinite representation. This matching establishes neither a priority claim nor a completed representation. The previously reported failure of all finite NPA levels would be compatible with a different finite lift; an alternative lift need not use an NPA word-moment space.

The [frontier](../FRONTIER.md) is unchanged in its universal obligations: H1 and H2 remain open; the evaluated L18 outer SDP has an unconditional quantum inclusion and a conditional reverse inclusion. L19 is a reviewed conditional local opening. L20/L21 retain their stated scopes and outstanding independent review. The next substantive step is sufficient universal admission in one fixed library, or a proved adequate alternative. Contact nondegeneracy alone does not supply that step.

## Verification record

The [standalone exact checker](checks/sun-wang-yu-match/verify_match.py) checks the entire source-to-Bell coefficient map, the deterministic saturation, the full SG2 state-matrix identity, both angle sign congruences, the overlap Jacobian, the determinant polynomial coefficient from the Gram trace formula, and the null-direction quartic identity. Its [receipt](checks/sun-wang-yu-match/result.json) records exact algebra only. The analytic proof above supplies positivity; a script PASS is not universal admission. See the [bounded independent review](SUN_WANG_YU_CONTACT_REVIEW.md) for the hypothesis and scaling audit. Scheduled research tasks remain paused.
