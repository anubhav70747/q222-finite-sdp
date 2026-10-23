# Fixed polynomial multipliers fail

The [proof](Q222_Fixed_Polynomial_Multiplier_Obstruction.pdf) supplies a rational
family of strictly valid Bell slacks. Every fixed product p(u)q(v) of nonzero
polynomials nonnegative on the real line fails to turn all sufficiently small
members of this family into polynomial sums of squares. The threshold depends
on p and q and is not evaluated.

This excludes every fixed polynomial Bob multiplier in the canonical
two-matrix completion. The critical polynomial has an explicit square-root
certificate on a half-line; algebraic testing functions remain available.

The JSON and pullback file give the complete source family, four-state
substitution, Bell matrix and exact coefficients. Both checkers use SymPy:

```text
python check_fixed_q_obstruction.py
python check_multiplier_obstruction_independent.py
```

The original recovered archive has SHA256
`5f77931c5a41e9c48108a6cb39d3ddbc92de937a51c410f559dbcf8c890e8bb4`.
`ORIGINAL_SHA256SUMS.txt` records its original members; current public hashes
are in the parent `PROVENANCE.json`. The public proof has been reorganized
around the separable-multiplier theorem and its completion corollary.
Compile `source/fixed_polynomial_multipliers.tex` twice with pdfLaTeX.
