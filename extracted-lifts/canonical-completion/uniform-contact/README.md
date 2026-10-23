# A uniform certificate near compatible measurements

Two positive coefficients can vanish along different limits of the source family.
Completing the squares separately introduces singular denominators. The weight
`phi = p^2 / (1 - p + p^2)` combines the decompositions first. The factor
`1-phi` cancels `1-p`; the bound `p^2 <= 32 Sigma` controls `phi/Sigma`.
Their common denominator is at least `3/4`. Keeping the entire
first-normal coupling preserves the cancellation that makes this possible.

The resulting certificate covers the stationary source family with
`0 < B, jB, T <= 2^-10`, `0 < p < min(j,1)`, `0 < H < 1` and `n >= 0`.
There is no bound on `j` or `jT`. The testing region is
`|u| <= 2^-12`, `|v| <= 2^-25`, with every real state and both state signs.
The source parameters and testing variables are independent.

The proof retains both state directions, the literal `n = 0` case,
nonnegative event additions and finite limits of conic combinations of the
displayed sources. Fixed testing functions turn the square identities into
an affine semidefinite system. Source-dependent formulas choose its Gram
matrices; the constraints themselves are linear in the nine source
coefficients and the free Gram entries.

The complementary testing regions and the other canonical source families
remain necessary for an exact lift of the full quantum correlation set.
The earlier native-source construction and its regional matrices remain
available in [contact-completion](../contact-completion/).
An [exact annulus zero](Annulus_Obstruction.md) shows why its earlier
separated arc construction cannot be extended unchanged to the larger family.

Read [the proof](Uniform_Contact_Certificate.pdf) for the cancellation and
its coefficient estimates. The evaluated construction has 49 Gram blocks,
total PSD order 283, 2269 Gram variables and 554 coefficient equations.
These counts describe the fixed-cap certificate cone, not a full-set lift.

## Matrices and reproduction

- `extraction/fixed_family.json` lists every monomial vector and rational weight.
- `extraction/affine_system.json` contains the evaluated affine constraints.
- `extraction/root_catalog.json` gives every square used to derive those vectors.
- `source_dictionary.py` gives the stationary source coefficients.

The source order is `(d,lambda,n,m,C,A,D,alpha,beta)`. Indices start at zero.
Gram entries are actual symmetric entries; the factor two for an
off-diagonal entry is already included in the equation records. Source
records are on the right-hand side. The polynomial target is
`(1+v^2) S_c`, after the state rotation in the proof.

With Python 3 and SymPy, run `python reproduce.py` from this directory.
The checks use exact rational and polynomial arithmetic. They verify the
original Bell slack, both state charts, the scalar remainders, all function
supports, the rational equation arrays and the finite-grid closedness argument.
Compile `Uniform_Contact_Certificate.tex` twice with pdfLaTeX for the proof.
All inputs needed by this folder are included here.
