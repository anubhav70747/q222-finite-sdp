# Polynomial half-line selectors fail; an algebraic enlargement

[The proof](Q222_Halfline_Order33_Counterexample.pdf) gives an all-positive
rational valid source outside the quadratic half-line completion cone.
A rational dual separates it with every auxiliary coefficient cancelling
exactly. The growth bound on the completion matrix makes this an obstruction
to polynomial selectors of every degree in the native Bob variable.

Allowing the selector to be quartic in sqrt(v) repairs the example. The
resulting rational affine pencil has total order 49, with block orders
12, 8, 12, 8 and nine scalar source signs, 148 auxiliaries, and no remaining
coefficient equations. Its 40 fixed scalar functions are defined on a
compact algebraic cover, including every projective endpoint.

The proved relation is C33 strictly contained in Cram, contained in the
canonical valid-source cone. Universal equality for Cram remains open.
This order-49 pencil is distinct from the Laurent-coefficient cone elsewhere
in this repository.

[Endpoint and contact constraints](Q222_Selector_Contacts.pdf) give a
two-page reduction of the remaining selection problem. On the principal
singular face, ramification adds only four coefficients; an interior
contact leaves one scalar freedom in the selector value. The note also
derives the necessary derivative conditions and covers literal zero faces.

Run `python reproduce.py` from this directory. The standard-library checks
verify the exact separation, both feasible repair certificates, every
affine coefficient, the source-validity identity and the fixed-function
data. No numerical solver is used. The `data/` directory contains the
literal matrices and auxiliaries; `proof/` contains editable LaTeX.

The recovered archive has SHA256
`2a02c671858f83e2c11f6a8d05ea6a3ddb70004f98658a6a297362d6b0596079`.
The public copy corrects one stale endpoint-document pointer and streamlines
the exposition. Its matrix data and mathematical statements are unchanged.
`SHA256SUMS.json` records the current public files.

## Closedness and exact dualization

Appendix B gives a positive moment weight annihilating all 148 auxiliary
matrices and a rational strictly feasible point with margin I49/100.
The projected source cone is closed. Its eight-orientation outer
approximation has total order 392 and 640 free auxiliary variables.
The assembly theorem makes equality with Q222 equivalent to universal
coverage by this canonical cone.

Additional independent checks, from this directory:

```text
python independent/check_order33_counterexample.py
python independent/independent_ramified_check.py
python independent/check_ramified_regularity.py
```

The first two use standard-library fractions; the regularity checker also
uses SymPy for the two exact rank computations. The independent feature
check expands all 157 affine basis identities for both signs directly
against the homogeneous physical state polynomial.
