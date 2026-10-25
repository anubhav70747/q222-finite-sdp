# Finite semidefinite representation of the simplest quantum correlation set

The central difficulty lies where measurements become compatible. There,
the state matrix loses rank while the Bell coefficients and measurement
angles can approach their limits at independent rates. The proof follows
the surviving state directions and completes the resulting quadratic forms.

- [Article](pdf/Q222_Paper.pdf): the mechanism and global argument, seven pages.
- [Complete proofs](pdf/Q222_Complete_Proofs.pdf): the 343-page technical companion.
- [Overleaf sources](Q222_Overleaf.zip).
- [Explicit representations](pdf/Explicit_Lifts.pdf) and [matrix data](extracted-lifts/).

The companion retains the complete singular-case argument. **An exact affine
pencil for the full eight-coordinate quantum set has not yet been extracted.**
The construction below resolves a singular source family on an explicit testing region.

On the review branch, an [exact order-seven block](extracted-lifts/c2-equal-angle-cone/)
covers two independently varying parameters of the C2 contact family.
Every feasible source is one tight contact source plus two nonnegative
event terms. This gives the source and behavior matrices directly and
strictly improves `1+AB` while retaining every quantum behaviour.

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
has order 293 and leaves the R1 testing strip
`|u| > 9/256`, `|v| < 1/16` unresolved.

An
[explicit rational example](extracted-lifts/canonical-completion/uniform-contact/Annulus_Obstruction.md)
shows why the earlier separated arc matrix fails for this source family:
it retains positive squares where the full Bell slack vanishes.
The enlarged contact certificate now includes this example.

## From canonical certificates to the full quantum set

The [three-page assembly proof](extracted-lifts/canonical-completion/full-polar/Q222_Polar_Assembly.pdf)
reduces the full quantum polar to eight output relabellings of one canonical
cone. An exact canonical pencil gives the eight-coordinate lift through
eight adjoint copies. The unresolved requirement is a single evaluated
certificate family admitting every canonical source at every testing angle.

The earlier [contact completion](extracted-lifts/canonical-completion/contact-completion/)
and [regional matrices](extracted-lifts/canonical-completion/fixed-c2-gates/Fixed_C2_Completions.pdf)
already give an order-635 construction valid at all measurements for their
smaller native source family. The long companion is retained because a
short, complete global replacement has not yet been established.

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
| Canonical angular validity | Exact order 49 in the Laurent coefficients; their Bell-source map is quartic |
| Determinant-weighted state polynomial | Fixed certificate of order 133, with 104 functions |
| Specified local source cones | Exact pencils of orders 62, 669 and 731 |
| Full-set size bound | Total PSD order at least 9 |

The [scientific data guide](extracted-lifts/README.md) gives each domain,
matrix convention and reproduction command. It also retains the larger
outer candidates and their proved or unresolved reverse inclusions.

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
