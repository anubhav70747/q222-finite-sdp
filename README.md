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
The following construction isolates the remaining step in a small matrix problem.

## The completion approach

The [canonical reduction](extracted-lifts/canonical-completion/Q222_Canonical_Completion.pdf)
eliminates Alice's angle. A Bell slack is valid exactly when two four-by-four
matrices admit a common completion J at every Bob angle v. This includes
zero source weights, all state signs and projective endpoints.

Choosing J as a polynomial on the full real line fails. A rational source
and a two-coordinate negative quadratic form give a short exact obstruction.
The [multiplier theorem](extracted-lifts/canonical-completion/polynomial-obstruction/Q222_Fixed_Polynomial_Multiplier_Obstruction.pdf)
shows that no fixed nonnegative separable polynomial multiplier repairs
this failure for every valid source. Its limiting polynomial is Motzkin plus
a constant.

On a half-line, square-root factors change the problem. The
[half-line construction](extracted-lifts/canonical-completion/halfline/Q222_Halfline_Completion.pdf)
represents quadratic completions by four PSD blocks of orders 8, 8, 4, 4.
With the source-sign constraints its total order is 33, with 44 auxiliaries.
Its matrices are affine in all nine canonical source coefficients. The
projected cone is closed, and an explicit rational witness certifies the
source that defeats every global polynomial completion.

Every source admitted by this pencil is valid. The remaining implication
is whether every valid canonical source admits such a half-line completion.
The exact arrays, rational witness and short proofs are together in
[canonical-completion/](extracted-lifts/canonical-completion/).

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
