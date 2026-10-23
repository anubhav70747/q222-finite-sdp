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
`|u| <= 9/256`, `|v| <= 1/16`, with every real state and both state signs.
The source parameters and testing variables are independent.
Every native contact lies strictly inside this cap. Its enlargement uses
the joint estimate `(j+1)/(j+1-p) * (1-phi+2phi/j) <= 4`, which retains
the dependence between two source coefficients while bounding their sum.

The proof retains both state directions, the literal `n = 0` case,
nonnegative event additions and finite limits of conic combinations of the
displayed sources. Fixed testing functions turn the square identities into
an affine semidefinite system. Source-dependent formulas choose its Gram
matrices; the constraints themselves are linear in the nine source
coefficients and the free Gram entries.

A second, ten-row matrix covers `|v| >= 1/16` for every Alice angle.
The two constraints share their source coordinates and nine sign blocks,
giving an explicit system of total order 293. Its remaining R1 testing
region is `|u| > 9/256`, `|v| < 1/16`. This strip and the other canonical
source families remain necessary for a lift of the full quantum set.
The earlier native-source construction and its regional matrices remain
available in [contact-completion](../contact-completion/).
An [exact annulus zero](Annulus_Obstruction.md) shows why the earlier
separated arc construction cannot be extended unchanged to the larger family.
The enlarged contact certificate now covers this zero.

Read [the proof](Uniform_Contact_Certificate.pdf) for the cancellation and
its coefficient estimates. The evaluated construction has 49 Gram blocks,
total PSD order 283, 2269 Gram variables and 554 coefficient equations.
These counts describe the cap certificate. The ten-row Bob matrix adds no
free variables or coefficient equations.

## Matrices and reproduction

- `extraction/fixed_family.json` lists every monomial vector and rational weight.
- `extraction/affine_system.json` contains the evaluated affine constraints.
- `extraction/root_catalog.json` gives every square used to derive those vectors.
- `source_dictionary.py` gives the stationary source coefficients.
- `bob-complement/bob_complement_order10.coo.json` gives the fifteen rational
  upper-triangular records of the second matrix.
- `joined_certificate.json` specifies the common source coordinates,
  constraints and domain of the order-293 system.

The source order is `(d,lambda,n,m,C,A,D,alpha,beta)`. Indices start at zero.
Gram entries are actual symmetric entries; the factor two for an
off-diagonal entry is already included in the equation records. Source
records are on the right-hand side. The polynomial target is
`(1+v^2) S_c`, after the state rotation in the proof.

With Python 3 and SymPy, run `python reproduce.py` from this directory.
The checks use exact rational and polynomial arithmetic. They verify the
original Bell slack, both state charts, the scalar remainders, all function
supports, the rational equation arrays, the finite-grid closedness argument
and the ten-row Bob identity and source bounds.
Compile `Uniform_Contact_Certificate.tex` twice with pdfLaTeX for the proof.
All inputs needed by this folder are included here.
