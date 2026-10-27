# Q222: comparison with exact geometry and semidefinite representation results

27 October 2026. This comparison follows the [preservation audit](../reconciliation/2026-10-27/preservation-audit/README.md). Existing manuscripts, proofs, matrices and result packages are preserved. Two bounded GPT-6-Luna/medium research passes supplied source comparisons; the lead checked the central statements against primary papers and corrected the distinctions below. This is a targeted comparison of nine relevant works, not an exhaustive priority search or a new audit of every cited proof.

Our two results must remain distinct: the reviewed qualitative finite-SDP existence argument with V01–V06 dependencies, and the evaluated L18 outer SDP whose reverse inclusion still requires H1/H2 or adequate alternative coverage. No comparison below changes either status.

## Which sets and representations are being compared?

The full behavior has eight coordinates: four marginal expectations and four correlators. Its four-correlator projection is a different convex set. A parametrization or polynomial membership test is also different from a fixed affine semidefinite lift. A finite level of NPA is one particular semidefinite construction, not the class of all finite lifts.

| Primary work and exact locator | Established object/result | Relation to this project |
|---|---|---|
| [Sun–Wang–Yu, arXiv:2610.34364v1](https://arxiv.org/html/2610.34364v1), 25 October 2026; Theorem 1, Theorem 4, Proposition 4, Lemma 6, Appendix E | Full behaviors: one pure real two-qubit behavior mixed with a local behavior. Frustrated-function geometry includes tied local/quantum bounds. | Direct overlap in full-set geometry. Our [completed match](SUN_WANG_YU_CONTACT_MATCH.md) gives a common contact proof. Their nonlinear membership criterion supplies neither K364 admission nor an evaluated fixed affine SDP. |
| [Barizien–Bancal, arXiv:2406.09350v1](https://arxiv.org/html/2406.09350v1), 13 June 2024; *Nature Physics* 21, 577–582 (2025); Theorem 1, Appendix A | Full behaviors: exact nonlocal extremal-point and realization conditions using the lifted arcsine/full-alternation characterization. | The continuous extreme-point description already characterizes the generating geometry of the same full quantum set. It does not provide the required fixed affine Gram lift. |
| [Masanes, quant-ph/0512100v1](https://arxiv.org/pdf/quant-ph/0512100), 13 December 2005; Theorem 1 | For N parties with two binary observables per party, extreme behaviors admit pure N-qubit/projective realizations. | Supports finite-dimensional extreme-point reduction. Finite-dimensional realization of generators alone is not finite SDP representability of their convex hull. |
| [Le Phuc Thinh–Varvitsiotis–Cai, arXiv:1809.10886v1](https://arxiv.org/html/1809.10886v1), 28 September 2018; *PRA* 99, 052108 (2019); Sections III–IV, Theorem 1 | Correlators: elliptope projection for arbitrary setting counts; an angle-inequality characterization when min(n,m)≤2, including TLM for (2,2). | The correlator projection already has a finite PSD-completion formulation. The angle inequalities are linear after arccos of the given correlators, not linear inequalities in the original correlators. Marginal consistency is absent. |
| [Navascués–Pironio–Acín, arXiv:0803.4290v1](https://arxiv.org/html/0803.4290v1), 31 March 2008; Sections 4.1–4.2 and 6.1 | Moment SDPs complete in the limit for the commuting-measurement model; rank-loop conditions can certify finite stopping for particular inputs. | Neither asymptotic convergence nor pointwise finite stopping implies one globally exact finite level. Our reported no-finite-NPA result and a different finite lift are compatible. |
| [Bamps–Pironio, arXiv:1504.06960v1](https://arxiv.org/abs/1504.06960v1), 27 April 2015; *PRA* 91, 052111; Sections II–III | Tight operator SOS decompositions for the singly tilted CHSH family and applications to self-testing. | Exact certificates for a restricted family are established precedents. They do not cover all affine supports or the doubly tilted limiting family used in our obstruction. |
| [Scheiderer, arXiv:2004.04196v2](https://arxiv.org/html/2004.04196v2), 28 January 2021; Theorem 3.10 and Corollary 3.19 | For closed convex semialgebraic K, universal tensor-SOS evaluation over every real closed extension characterizes finite semidefinite representability. Theorem 3.10 additionally controls tensor rank/block size. | This is the qualitative criterion used in the retained proof. It does not evaluate K364 or the L18 reverse inclusion. Its planar second-order-cone theorem cannot simply be applied to the eight-dimensional behavior set. |
| [Dritschel, arXiv:1811.06005v4](https://arxiv.org/html/1811.06005v4), revised 9 August 2024; Section 4, Theorem 4.1 | A nonnegative matrix-valued two-torus polynomial of degree (d1,d2) has at most 2d2 matrix squares, with factor degrees (max(dA,m),2d2−1), where m is finite for that input. | Fixed square count is not fixed feature dimension. The first factor degree is not controlled uniformly by input degree. The [comparison corollary](ANALYTICAL_COMPARISON_BRIDGES.md#2-a-small-square-count-does-not-give-a-uniform-trigonometric-library) proves divergence on our existing physical limiting family. |
| [Scheiderer, *Spectrahedral shadows*, arXiv:1612.07048v3](https://arxiv.org/abs/1612.07048v3), revised 4 December 2017; introduction and counterexample results | Some compact convex semialgebraic sets have no finite semidefinite representation. | Compactness, semialgebraicity and finite parametrization alone are insufficient. Those examples are not a nonrepresentability theorem for Q222. |

For clarity, the familiar correlator SDP asks for real a,b such that

\[
\begin{pmatrix}
1&a&E_{00}&E_{01}\\
a&1&E_{10}&E_{11}\\
E_{00}&E_{10}&1&b\\
E_{01}&E_{11}&b&1
\end{pmatrix}\succeq0.
\]

It certifies only the four correlators. Specifying four additional nonzero marginals is not handled by this matrix. This distinction prevents the correlator result from being misidentified as a solution of our full-behavior construction.

## What the comparison lets us simplify

1. **Contact geometry:** use the exact frustrated-source matching already published. It supplies finite-contact coercivity, uniqueness and full curvature in one argument; retain the classification, singular tensor identities and their hypotheses.
2. **Constructive membership:** eliminate the local behavior and its mixture weight by the finite interval test in [Bridge 1](ANALYTICAL_COMPARISON_BRIDGES.md). This follows from Sun–Wang–Yu and local-polytope facets, includes both endpoint weights, and avoids restricting the local mixture to two vertices. It remains a nonlinear search over the real-qubit family.
3. **Choose the right fixed-library route:** [Bridge 2](ANALYTICAL_COMPARISON_BRIDGES.md#2-a-small-square-count-does-not-give-a-uniform-trigonometric-library) shows why a uniform ordinary trigonometric degree bound cannot replace the fractional library along the existing doubly tilted family. This does not exclude a nonsmooth fractional construction or an arbitrary finite lift.

The [compact strict-family theorem](COMPACT_NATIVE_CONTACT_COVER.md) remains useful on each stated compact subset. Compactness does not turn its source-dependent grid size into one grid for all sources approaching walls. Conversely, a genuinely uniform finite real Gram identity would directly define a finite affine source cone; it need not first be converted into the tensor criterion. These are alternative proof routes, and only the applicable route's obligations should be imposed.

## Positioning and open work

The project's own completed comparison is also stronger than a finite-lift existence statement: the [previously delivered common-tilt planar package](../extracted-lifts/common-tilt-planar/README.md), now preserved in Git, gives an evaluated (5,4)-block membership lift for W=(marginal sum, CHSH), on the same region as the imported infinite-NPA-depth directions. Its 32 physical images give exact support on 32 planes and strict strengthening of every fixed NPA relaxation; equality with the full quantum set does not follow. The [original coverage comparison](../extracted-lifts/common-tilt-planar/COVERAGE_COMPARISON.txt) distinguishes this from the single-marginal projection, balanced C1 source cone, C2 equal-/unequal-angle cones, and signed asymmetric fixed-input SDPs. The latter still do not evaluate one common affine lift of W3=(A0,B0,CHSH). This preservation pass reuses the stated prior proof/review evidence and does not claim a fresh audit of all these results.

The comparison does not support claiming the first exact characterization, the first finite parametrization of full behaviors, or the first finite SDP for the correlator projection. The representation claim must specify the **full behavior set and an arbitrary finite affine lift**, distinguish qualitative existence from an evaluated exact pencil, and compare against the checked versions above. This bounded pass does not establish priority over all literature.

The immediate constructive obligations are unchanged: universal oriented C1 and strict C2 admission, or another adequate fixed sound family; then apply the retained closure/orbit argument and prove the reverse inclusion for the actual matrices. L20/L21 independent reviews remain outstanding. None of the cited theorems or the two bridges supplies missing H1/H2 coverage.
