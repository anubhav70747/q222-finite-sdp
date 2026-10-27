# Exact-lift frontier — 27 October 2026

The complete rational **outer SDP** is evaluated. Its unconditional exactness remains open. The reviewed qualitative finite-SDP existence theorem is a separate established input with the V01–V06 dependencies retained. L19 supplies independently reviewed [C2 n=0 wall-to-strict opening](independent-review/C2_N0_OPENING.md). L20 now adds an [exact global C1 source slice and a five-native-parameter admission tube](independent-review/C1_BALANCED_SLICE.md). L21 extends the [conditional source opening to all three zero-weight walls](independent-review/ZERO_WEIGHT_OPENING.md). Universal H1/H2 remain open.

This update incorporates upstream through [`e37c26d`](https://github.com/anubhav70747/q222-finite-sdp/commit/e37c26d15ac2bb11d94c90b98f7a2dd152b6ee19), preserving the independent L19 review, guarded real-factor checker and new C1 scalar-Schur reformulation. [Fourteen earlier exact checker runs](reconciliation/2026-10-27/verification/latest_checks.json) and [the two new exact certificate checks](reconciliation/2026-10-27/completion-advance/verification.json) passed on their stated contracts. The two new runs are same-assistant checks, not independent review. L19 has independent [analytical](complete-proof-materials/materials/historical/analysis/c2_n0_opening_analytic_verification_2026_10_27.txt) and [bounded exact](complete-proof-materials/materials/historical/analysis/c2_n0_opening_exact_verification_2026_10_27.txt) reviews of its conditional local theorem.

## Sun–Wang–Yu matching completed

The [new exact comparison and proof](independent-review/SUN_WANG_YU_CONTACT_MATCH.md), based on main `1b29383`, matches their canonical Bell operator to the full SG2 parent. Every canonical d>0 source has correlator signs (+,−,+,+). Both native C1/C2 families have tied local and quantum bounds, so the relevant external inputs are Proposition 4 and Lemma 6, not the strict-advantage theorem. The resulting short argument proves finite-contact coercivity, uniqueness and the complete five-variable contact Hessian, including the allowed zero-weight faces. It is an alternative proof of existing geometry, with [bounded independent review](independent-review/SUN_WANG_YU_CONTACT_REVIEW.md) and [13 exact matching checks](independent-review/checks/sun-wang-yu-match/result.json).

Their one-real-two-qubit-plus-local membership description does not supply fixed affine SDP admission. H1/H2 and global reverse inclusion remain open. This update preserves the original manuscript, singular proofs and the outstanding L20/L21 reviews. Scheduled research tasks remain paused.

## Preservation and broader comparison

The [Git preservation audit](reconciliation/2026-10-27/preservation-audit/README.md) confirms no files were deleted and no existing manuscript, PDF, extracted lift, archived proof or checker was changed between the audited baselines and `ac8861a`. Original checklist entries and component coverage remain recoverable; the audit identifies the index/status changes separately.

The [nine-work literature comparison](independent-review/LITERATURE_COMPARISON.md) distinguishes full-behavior geometry, correlator SDPs, NPA convergence and arbitrary finite lifts. Two [analytical consequences](independent-review/ANALYTICAL_COMPARISON_BRIDGES.md) eliminate the local mixture variables by a facet interval test and derive trigonometric-degree divergence along our existing doubly tilted family from the retained smooth-feature obstruction. These simplify membership search and rule out an unsuitable uniform-degree route; they do not supply H1/H2. The [nine exact checks](independent-review/checks/comparison-bridges/result.json) cover their stated algebra and examples. Original proof and manuscript files are unchanged.

## What would complete the current candidate

The separately completed [common-tilt planar lift and 32 physical projections](extracted-lifts/common-tilt-planar/README.md) are now preserved in Git with their original proof, actual matrices, checkers and delivery archives. They give exact support on the stated 32 planes and, using the imported nontermination theorem, strictly strengthen every fixed NPA level. They do not establish full-set membership or H1/H2. The common affine construction for W3 also remains open. This publication imports existing local results; no original manuscript or certificate was rewritten.

Let

\[
K_{364}=\{s\ge0:(1+v^2)F_s=z_{6,12}^TPz_{6,12},\ P\succeq0\},
\qquad z_{6,12}=(u^{i/6}v^{j/6}\psi_l)_{i=0..6,j=0..12,l=1..4}.
\]

Equality is coefficientwise for the complete four-state SG2 parent and independent testing angles. A sufficient route consists of two remaining universal statements:

| Input | Exact required statement | Status |
|---|---|---|
| H1 | Every genuine physical strict C1 native source with 0<k<=1 belongs to K364 | Open |
| H2 | Every genuine physical strict C2 source in its complete native domain belongs to K364 | Open |

An unweighted Z18 or Z6 certificate is sufficient but stronger than direct weighted K364 admission. Another adequate fixed finite family of sound cones could also complete a lift. H1/H2 are sufficient contracts for this route, not asserted necessary conditions for exactness of every possible orbit construction.

Once H1/H2 hold, use Kbar=K364+sigma(K364). The verified native reflection covers the opposite C1 orientation; coefficient-space density supplies permitted singular contacts; event and SG12/SG13 cones cover the other classified extremes. The full 128-element physical orbit absorbs sigma. Attained source duality then supplies the reverse inclusion for the already delivered global matrices. No separate uniform limiting-Gram bound is required on this route after universal strict admission in the same fixed cone/sum.

The source-dual target is nonnegative pairing against **every** PSD truncated moment. On a genuine contact source, the trace-normalized minimum is at most zero by contact evaluation; membership makes that minimum exactly zero. Individual pairings need not vanish. General PSD pseudo-moments cannot be replaced by physical evaluations or assumed to have representing measures.

## Current constructive results and their limits

| Result | Completed statement | Remaining statement |
|---|---|---|
| C1/C2 contact matching, L22 | Common proof of finite-contact coercivity, uniqueness and full Hessian using the exact SG1/SG2 map and the frustrated-case lemmas | Uniform all-testing fixed-library admission; no new H1/H2 coverage |
| C1 balanced slice and source tube, L20 | Exact order-five source/behavior pencils on the stated slice; a rational global Z18 interval certificate and explicit five-parameter all-testing tube | The rest of the native C1 domain; independent review of the new result |
| Zero-weight openings, L21 | Fixed carriers and coefficient maps for C2 n, C1 beta and C1 D opening from full wall-complement PD Grams | Whole-wall PD selection and coverage outside the local source neighborhoods |
| C1 Z18 reduction | Nine source-only gauges, two exact 4x4 Schur residuals | Simultaneous PSD selection on the actual oriented native domain |
| C1 Bob endpoints | One full-state polynomial matrix is PD near v=0 and v=infinity for every strict source | A globally PSD polynomial within a fixed adequate library |
| C1 scalar-Schur gate | Four scalar sign conditions equivalent to the sufficient 40-feature mechanism; exact endpoint and conditional contact factors | One common source-only quartic making all four polynomials nonnegative for every positive Bob angle; no universal choice is proved |
| C1 contact jet | Rank-six contact Gram and two-sided Bob-neighborhood positivity for every strict oriented source, including k=1 | Join contact and endpoint mechanisms at all intervening Bob angles; the separate local recipes cannot simply be patched without controlling degree |
| C2 open c ends, L12 | Exact finite-c admission near both strict open ends, uniform on compact interior face sets; half-power features suffice there | Arbitrary finite c and coupled approaches to face walls |
| C2 n=0 wall | Universal coefficient identities; finite exact PSD witnesses; smaller generic 16-parameter family | Whole-wall PSD selection; exceptional specializations in the smaller family remain controlled only within their stated scope |
| C2 n=0 opening, L19 | A PD Gram on all four larger wall contact complements gives an explicit source neighborhood of finite-n admission in28 fixed half-power features (49 weighted K364 features) | Wall points lacking that PD premise, and the rest of the strict native domain |
| Global L18 SDP | All rational coefficients, 128 actions, orthants, face/event blocks and affine equations; quantum-to-feasible inclusion | H1/H2 or sufficient alternative coverage for unconditional reverse inclusion |

The C1 contact-jet result is a **testing-angle neighborhood for each source**. L19, the L20 coefficient tube and L21 are **source neighborhoods with full all-testing identities**. Their distinct hypotheses and quantifiers remain essential.

## The local C2 opening now completed

An independent n-event square fails to vanish at the moving physical contact. The fractional carrier replaces it by contact-zero roots. After pure-source contact normalization, subtracting that carrier leaves 72 coefficient keys. The fixed wall contact-difference map has 104 symmetric variables and rank66; its image is exactly the kernel of the six state/angular contact jets.

A constant rational right inverse therefore supplies every residual coefficient. It changes only the 12x12 integer contact block; all three fractional contact blocks can remain fixed. A quantitative norm bound preserves PSD near a wall certificate PD on all four contact complements. The proof includes the negative WX couplings created by the carrier rather than discarding them.

Delivered: [proof](independent-review/C2_N0_OPENING.md), [standalone SymPy checker](independent-review/checks/gram-lift/ends/n0-opening/verify_n0_opening.py), [fixed matrix data](independent-review/checks/gram-lift/ends/n0-opening/fixed_map.json), [exact wall and strict certificates](independent-review/checks/gram-lift/ends/n0-opening/exact_results.json). The verified strict example has n=221/19042140160>0. Whole-wall or general H2 admission is not claimed.

## New exact C1 slice and quantified source tube

On d=4m, C=m, D=A-m, beta=alpha-m, lambda=n=0, canonical orthant constraints and A*alpha>=2m^2 are necessary and sufficient for validity. The resulting source LMI has blocks(2,1,1,1), total order5, and an evaluated behavior dual. Necessity has an explicit maximally entangled witness with orthogonal measurement pairs; sufficiency follows from the full boundary curve and nonnegative event-pair additions.

Eight fixed rational matrices, each greater than I/4, certify the whole native interval k=f=g=1, a=t=1-r, 1/3<=r<=1/2 in Z18. The fixed contact right inverse then gives an all-testing certificate whenever the actual normalized coefficient vector is within eta^2/246960 of that curve, for r in[1/3+eta,1/2-eta]. This is a genuine open region in all five native parameters. [The proof, evaluated data and standalone checker](independent-review/C1_BALANCED_SLICE.md) state every restriction. Same-assistant exact checking is complete; independent review is not claimed.

The [general zero-weight theorem](independent-review/ZERO_WEIGHT_OPENING.md) supplies the beta and D carriers missing from L19. The residual maps are both72x104 of rank66, with corrections confined to the integer block. Wall relative-PD is still a hypothesis.

## Routes already excluded or insufficient

* The unrestricted C1 Z16 target is refuted; oriented Z18 and K364 are not excluded by that counterexample.
* The selected one-direction C2 scalar completion fails at a strict seed for every real tau. All lower polynomial gauges also fail when its old complete top allocation is frozen. Those are restricted-ansatz failures.
* The mixed determinant identity represents det(B)F and does not give fixed source-affine admission for F after cancellation.
* The fractional carrier alone is a positive-quadrant half-power SOS, not a Bell-source cone element. L19 supplies a residual only under its local PD-wall hypothesis.
* The unweighted global half-power cone fails at an exact strict source, consistently with its proved local coverage.
* The angular order13 pencil has a quartic source map. It does not replace a source-affine full-state lift.

## Remaining completion checklist

1. Prove H1 and H2, or state and prove adequate alternative fixed finite orbit coverage. Preserve actual native inequalities, independent tests and all state normals. Universal strict-source coverage includes every strict point arbitrarily near a wall.
2. Verify the new PSD mechanism and all equality/saturation cases independently. L19's independent review is complete for its stated all-four-wall-complement PD premise. L20/L21 have new exact same-assistant certificates and proofs; independent review remains outstanding. The remaining universal constructions still require their own proofs and review.
3. Apply the established reflection/density/closed-sum theorem and conclude equality with the full Bell-validity cone; discharge reverse inclusion for L18.
4. Integrate the final theorem, domains, formulas and both inclusions into the article and companion. A shorter complete replacement proof remains an additional objective; local results do not replace the imported long native argument.
5. Rebuild changed manuscript PDFs, reproduce any changed matrices, synchronize the final source/data/Overleaf release and then close L08–L10. No manuscript PDF was changed or rebuilt in this update.

The matrices have 897 PSD blocks of total order51088, 423168 auxiliary scalars and3456 affine equalities. Their size and materialization do not establish admission or computational practicality. The remaining universal positivity statements are substantive research obligations and remain explicitly open.

The latest task3 delivery adds a pointwise verifier, not an admission theorem. Use the [reviewed real-factor variant and its scope](complete-proof-materials/STRICT_SOURCE_REDUCTIONS_UPDATE.md), which rejects the original program's imaginary-factor/negative-weight false positive. The original package remains preserved.

The [reviewed scalar-Schur reformulation](complete-proof-materials/materials/historical/analysis/c1_global_schur_gate_verification_2026_10_27.txt) sharpens the remaining sufficient C1 target. The endpoint skew is forced; the contact determinant factors also need both individual kernels and both zero soft derivatives. This is an exact reformulation of one sufficient mechanism, not new H1 coverage or a K364 obstruction.
