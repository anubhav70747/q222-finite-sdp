# Fixed functions at the singular contact

At a zero of a positive quadratic form, every square in a positive decomposition
must vanish. The native contact couples two state directions, so completing
them together is essential. Their combined positive coefficients control the
remaining terms even when one event weight is zero. Uniform matrix bounds then
give a finite family of testing functions with fixed square-root arguments.

- [Weak-contact proof](Contact_Completion.pdf), with its [editable source](Contact_Completion.tex).
- [Strong-contact continuation](Strong_Contact_Completion.pdf), with its [editable source](Strong_Contact_Completion.tex).
- [Reduced rational coefficient system](contact_affine_parity_reduced.json).
- [Unreduced coefficient system](contact_affine_independent.json), retained for comparison.

The two proofs cover the native source family defined in the
[regional completion note](../fixed-c2-gates/Fixed_C2_Completions.pdf), with
`0 < H < 1`, `0 < z,y <= 2^-64 H`, `0 < r < 1`,
`0 < t < 1+zr`, and `n >= 0`. Source parameters and testing variables remain
independent. Nonnegative additions to the eight event coefficients are allowed.

The contact construction covers `|u| <= 2^-12`, `|v| <= 2^-25`. The previously
derived order-20 and order-10 matrices cover the remaining testing regions,
including the projective endpoints. Their intersection with this coefficient
system has total PSD order **635**, with **4024 auxiliary variables** and
**619 coefficient equations**.

This gives a valid-source cone containing the displayed native family and its
nonnegative event additions. Universal admission of all canonical sources
remains open; these matrices do not yet give the exact full quantum set.

## Reading the matrices

The source order is `(d, lambda, n, m, C, A, D, alpha, beta)`. Impose their
nonnegativity and positive semidefiniteness of every Gram block. A record
`[row, block, i, j, value]` contributes `value * G_block[i,j]` to that row.
Rows and matrix indices are one-based; off-diagonal product columns already
include the factor two. Source records form the right-hand side.

There are 63 Gram blocks of total order 596. The nine source-sign blocks give
contact order 605. The coefficient system has 4021 free symmetric Gram entries,
20610 Gram records, and 120 source records. Appending the unchanged
[`alice_order20.coo.json`](../fixed-c2-gates/alice_order20.coo.json) and
[`bob_order10.coo.json`](../fixed-c2-gates/bob_order10.coo.json), with the
Alice matrix's three free auxiliaries, gives the joined system above.
The equation count is not a claim about row rank.

Averaging the eight formal radical-sign changes removes cross-character
entries without changing the projected cone. The unreduced system retains
44 Gram blocks, 4697 Gram variables, and 948 coefficient equations.
`original_block_and_indices` supplies the exact embedding of each reduced
block in the original matrices. Source-dependent formulas in the proofs
choose Gram witnesses; they are not additional nonlinear SDP constraints.

## Reproduce

The affine arrays use only Python's standard library:

```text
python reconstruct.py
python reduce_parity.py
```

The symbolic identities and coefficient estimates use SymPy 1.14.0:

```text
python weak/check_weak.py
python weak/check_spans.py
python weak/check_source_bounds.py
python strong/scripts/strong/extract_lII.py
python strong/scripts/strong/verify_native_contact.py
python strong/check_strong_algebra.py
python strong/scripts/strong/bound_lII.py
python strong/scripts/strong/monomial_budgets.py
python strong/scripts/strong/export_source_witness.py
```

The scripts produce their intermediate pickle files locally from the included
symbolic expressions. The distributable matrices and estimates are plain JSON.
No numerical optimization or source sampling is used.

Compile each TeX note twice with pdfLaTeX. The complete source formulas,
grouped Gram entries, monic divisions, and rational coefficient bounds provide
the mathematical derivation independently of the check outputs.
