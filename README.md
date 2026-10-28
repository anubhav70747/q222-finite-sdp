# Finite semidefinite representation of the simplest quantum correlation set

**Infinite NPA depth, finite lift:** the [two-page note](planar-npa-lift/note/main.pdf), [PRL Letter](planar-npa-lift/article/main.pdf) and [complete proofs](planar-npa-lift/proofs/main.pdf) give an exact finite SDP for the full marginal-sum/CHSH quantum region, although every fixed NPA level misses some of its supporting lines. The main article proves why equal overlaps suffice in every nonclassical support direction; the resulting quartic matrix gives explicit blocks of orders 5 and 4. The [package](planar-npa-lift/) includes the author-approved abstract, reproducible sources, a 16-page complete proof companion and exact data for all 32 relabelled projections. The new physical proof shows that two real pure equal-overlap strategies suffice for every planar point, with simultaneous SDP witness ranks at most two. The Letter is three pages including references; the short note remains two pages.

**Dependency repair, 28 October:** the [R5.5 and R6.5–R6.7 correction](complete-proof-materials/DEPENDENCY_REPAIR_2026_10_28.md) restores proof bodies omitted from the compressed companion. Two bounded analytical reviews verify the full-parent extraction, its actual R6 source-domain transfer and the stated R6/R7 projective endpoint joins. Their 50 exact checks support those named identities. The historical global reviews remain on file; this repair is not a new exhaustive audit. H1/H2 and global candidate exactness remain open.

**Complete material publication, 26 October:** [proofs, data and exhaustive checklist](complete-proof-materials/). This includes the original 215-component source archive, all 45 C1 and 49 C2 case records with direct proof links, the complete historical analytical reviews, and newly fetched proof/checker packages. The [26 October reconciliation](complete-proof-materials/RECONCILIATION.md#final-qualitative-conjunction-26-october-13-utc) records the qualitative-existence review using those body reviews; its use with the compressed companion includes the dependency repair above. The [conditional global outer SDP is now fully evaluated](complete-proof-materials/CONDITIONAL_GLOBAL_PENCIL_UPDATE.md). Its exactness remains conditional on direct oriented C1 and strict C2 admission into K364; the verified density and orbit arguments then cover the singular source limits. These admission/exactness obligations remain open as L08-L10.

**Current frontier:** [the completion map](FRONTIER.md) records the remaining universal C1/C2 admission proofs and the evaluated global matrices. The [C2 wall-to-strict opening theorem](independent-review/C2_N0_OPENING.md) now gives an explicit PSD neighborhood in 28 half-power features at every n=0 wall point with a Gram positive on all four contact complements. It does not establish whole-wall or general strict C2 admission.

**New global C1 result:** the [balanced source slice](independent-review/C1_BALANCED_SLICE.md) has an exact rational SDP of total order5, a fixed global Z18 certificate, and an explicit five-native-parameter admission tube. The [zero-weight opening theorem](independent-review/ZERO_WEIGHT_OPENING.md) now handles the C1 beta=0 and D=0 walls as well as C2 n=0, conditional on a Gram PD on all four larger wall complements. These results enlarge the proved admission regions; universal H1/H2 remain open.

**New exact construction:** [L12 and its proof package](complete-proof-materials/CHECKLIST.md) give a fixed 16-feature PSD certificate near both strict open C2 ends. The complete finite-c identity and graded Schur positivity have been independently checked; the result is uniform on compact subsets of the stated face interior.

The central difficulty lies where measurements become compatible. There,
the state matrix loses rank while the Bell coefficients and measurement
angles can approach their limits at independent rates. The proof follows
the surviving state directions and completes the resulting quadratic forms.

- [Article](pdf/Q222_Paper.pdf): the mechanism and global argument, seven pages.
- [Complete proofs](pdf/Q222_Complete_Proofs.pdf): the 347-page technical companion.
- [Overleaf sources](Q222_Overleaf.zip).
- [Explicit representations](pdf/Explicit_Lifts.pdf) and [matrix data](extracted-lifts/).

The companion retains the complete singular-case argument. The explicit
construction has advanced from that argument toward one exact affine pencil
for the full eight-coordinate quantum set:

- **Assembly.** Given an exact canonical pencil, eight adjoint copies give
  the full lift. The sixth-root candidate for that pencil has an explicit
  rational strictly feasible point.
- **Reduction.** Every extreme Bell inequality is an event, lies on one of two
  closed faces with exact pencils, or is a strict C1/C2 contact source. The
  construction therefore reduces to admitting those two contact families in
  one fixed finite collection of cones ([route](independent-review/EXACT_LIFT_ROUTE.md)).
- **Coverage.** Every compact strict native subfamily has a finite cone, whose
  size depends on the subfamily
  ([compact cover](independent-review/COMPACT_NATIVE_CONTACT_COVER.md)). The
  sixth-root lift admits the hardest sources found, with exact certificates,
  and its limit problems at the `c → 0` and `c → ∞` ends are solved in closed
  form on whole faces ([sextic lift](independent-review/EXACT_LIFT_CANDIDATE.md)).
  Universal strict admission in the same fixed closed cone remains open.
  The verified coefficient-space density theorem then handles the singular
  sources; no separate uniform limiting-Gram estimate is required on that route.
- **Explicit blocks.** Exact cones cover singular source families on explicit
  testing regions, as described below.

The construction below resolves a singular source family on an explicit testing region.

On the review branch, an [exact order-seven block](extracted-lifts/c2-equal-angle-cone/)
covers two independently varying parameters of the C2 contact family.
Every feasible source is one tight contact source plus two nonnegative
event terms. This gives the source and behavior matrices directly and
strictly improves `1+AB` while retaining every quantum behaviour.

An [order-six block](extracted-lifts/c2-unequal-angle-curve/) now covers a
compact curve with unequal measurement angles. Its exact moment
representation needs at most three contact points. The
[general contact equations](independent-review/C2_CONTACT_GEOMETRY.md)
identify an additional fourth-degree relation, the next ingredient for a
fixed lift of their full three-parameter convex hull.

## One cancellation across the singular source limits

The [uniform contact proof](extracted-lifts/canonical-completion/uniform-contact/Uniform_Contact_Certificate.pdf)
has a three-page main argument and eleven pages of analytical appendices.
Two square completions contain coefficients that become singular in different
source limits. The weight `phi = p^2 / (1-p+p^2)` combines them before either
coefficient is bounded. Keeping the entire first-normal coupling preserves
the cancellation. Both state directions are then controlled by one positive
matrix.

This proves a fixed certificate for the broader stationary R1 source family
on `|u| <= 9/256`, `|v| <= 1/16`, including vanishing event weights and
unbounded source ratios. Its [evaluated rational system](extracted-lifts/canonical-completion/uniform-contact/)
has total PSD order 283. The proof gives every square, its fixed function
span, the affine equations and a finite-grid proof of closedness.

The cap contains every native contact in this source box. A joint bound
on the source coefficients makes that enlargement possible without adding
matrix rows. A ten-row matrix covers `|v| >= 1/16` at every Alice angle.
The [joined system](extracted-lifts/canonical-completion/uniform-contact/joined_certificate.json)
has order 293. The adjoining R1 testing strip `|u| > 9/256`, `|v| < 1/16`
is the next region to join.

An
[explicit rational example](extracted-lifts/canonical-completion/uniform-contact/Annulus_Obstruction.md)
shows why the earlier separated arc matrix fails for this source family:
it retains positive squares where the full Bell slack vanishes.
The enlarged contact certificate now includes this example.

## From canonical certificates to the full quantum set

The [three-page assembly proof](extracted-lifts/canonical-completion/full-polar/Q222_Polar_Assembly.pdf)
reduces the full quantum polar to eight output relabellings of one canonical
cone. An exact canonical pencil gives the eight-coordinate lift through
eight adjoint copies. The input to assemble is a single evaluated
certificate family admitting every canonical source at every testing angle.
The constructions above and the [review route](independent-review/EXACT_LIFT_ROUTE.md)
cover growing parts of the source space.

The earlier [contact completion](extracted-lifts/canonical-completion/contact-completion/)
and [regional matrices](extracted-lifts/canonical-completion/fixed-c2-gates/Fixed_C2_Completions.pdf)
already give an order-635 construction valid at all measurements for their
smaller native source family. The long companion carries the complete
singular-case argument, and these explicit constructions make growing parts
of it concrete.

The [canonical reduction](extracted-lifts/canonical-completion/Q222_Canonical_Completion.pdf)
also turns validity into a common completion of two four-by-four matrices
at every Bob angle. Its simpler fixed-selector candidates are now understood:
the [order-33 half-line cone](extracted-lifts/canonical-completion/halfline/)
and the [order-49 ramified cone](extracted-lifts/canonical-completion/ramified-completion/)
are strict sub-cones, with [exact counterexamples](extracted-lifts/canonical-completion/ramified-obstruction/).
The [fixed polynomial multiplier obstruction](extracted-lifts/canonical-completion/polynomial-obstruction/)
explains why increasing a separable polynomial degree does not repair the
global construction. Their proofs and matrices remain available.

## Other explicit results

| Construction | Representation |
|---|---|
| Tilted-CHSH projection | Exact order 4, without auxiliaries; both inclusions |
| Balanced C1 source slice | [Exact order 5](independent-review/C1_BALANCED_SLICE.md), with the stated linear relations and canonical source orthant |
| Canonical angular validity | [Exact order 13](extracted-lifts/global-weighted/angular-order13/) in Laurent coefficients; their Bell-source map is quartic |
| Determinant-weighted state polynomial | Fixed certificate of order 133, with 104 functions |
| Specified local source cones | Exact pencils of orders 62, 669 and 731 |
| Full-set size bound | Total PSD order at least 9 |

The [scientific data guide](extracted-lifts/README.md) gives each domain,
matrix convention and reproduction command. It also retains the larger
outer candidates and the reverse inclusions proved for them so far.

## Build

Compile each entry in `latex/` twice with pdfLaTeX:

```text
main.tex
Supplement.tex
Explicit_Lifts.tex
Uniform_C1_Certificate.tex
```

The bibliography is included. The independent completion notes have their
own editable sources beside their matrices. Earlier versions remain in Git
history; the current repository contains the scientific sources, formulas,
matrix arrays and exact reproduction code.
