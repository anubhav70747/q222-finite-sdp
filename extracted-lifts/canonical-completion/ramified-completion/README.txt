Q222: exact half-line counterexample and one fixed algebraic enlargement
======================================================================

PROVED

The quadratic half-line completion cone C33 is a proper subcone of the full
nonnegative nine-coordinate canonical valid-slack cone V9.

An all-positive rational counterexample, in coordinate order
(d,lambda,n,m,C,A,D,alpha,beta), is

  cbar = (1, 2^-34, 2^-14, 2^-24, 1/576, 161, 1/2, 76, 1/250).

The integer separating functional has coefficients

  (-23803618, 906393494300170, 2744379904, 1378590360912,
   5935589516, 20425, 2600028, 91248, 278104928).

It is nonnegative on every C33 certificate and has exact value

  -5417798011632274183 / 9663676416000

at cbar. The six rational dual moments and all four exact LDL factorizations
are supplied. The witness also excludes the previously established contact
source cstar. The principal-minor growth bound excludes every polynomial
half-line selector, since every such selector necessarily has degree <=2.

A single fixed enlargement repairs both sources. Use v=16*s^2, s>=0, and
J(s) polynomial of degree <=4. The four Gram blocks have orders 12,8,12,8.
Together with nine scalar source signs this gives an evaluated affine pencil
of order 49, with 148 auxiliary variables and no coefficient equations.

  C33 is strictly contained in C49_ram, which is contained in V9.

NOT CLAIMED

The equality C49_ram=V9 remains unproved. This is not an affine lift of the
full eight-coordinate Q222. This order-49 completion pencil is distinct from
the previously solved order-49 Laurent-coefficient cone. The existing complete
companion is preserved and is not claimed to be replaced by this note.

READ

  Q222_Halfline_Order33_Counterexample.pdf
  proof/halfline_counterexample.tex

The six-page proof gives the exact separating matrices, source validity,
strict enlargement, affine parameterization and compact testing cover.

VERIFY

  python3 reproduce.py

Python 3 standard library only; no solver, CAS, network, or tolerance is used.
The script verifies both sources, the rational dual and repair certificates,
all affine basis matrices and all exported auxiliary values. It also checks
the retained small determinant/AM-GM source identity. Runtime is short.

MATRIX AND FUNCTION CONVENTIONS

data/ramified_order49_pencil.json contains nine source coefficient matrices
and 148 auxiliary coefficient matrices. Sparse entries are upper-triangular
[row,column,rational] triples with zero-based indices. Reflect off-diagonal
entries. The constant matrix is zero. The block offsets are:

  source signs 0..8; A_G at 9; B_G at 21; A_H at 29; B_H at 41.

There are 411 nonzero coefficient entries across all 157 basis matrices.
data/ramified_order49_source_auxiliaries.json and
data/ramified_order49_cbar_auxiliaries.json give actual feasible values for
all 148 auxiliaries. The PSD certificates are the referenced rational Grams.

data/linear_coefficient_equations.json explicitly gives the 80 source-linear
equations in 228 symmetric Gram entries before elimination. The evaluated
pencil parametrizes their full affine solution space.

data/compact_fixed_functions.json lists all forty bounded algebraic testing
functions, their source-independent domain, algebraic relations, multiplier,
uniform bounds and endpoint values. data/fixed_functions.json gives their
equivalent finite-angle formulas. All source signs and state signs remain.

The compact cover uses a0^2+a1^2=1, a0>=0, rho^2+tau^2=1, rho,tau>=0,
epsilon=+/-1, and kappa^2=tau^8+256*rho^8 with kappa>0. It satisfies
kappa^2>=1/8. Every feature is bounded by a number less than 3 on unit states.
The represented target is a1^2 times the physical homogeneous slack; the
fixed multiplier is positive on a dense testing set.

DETERMINISTIC BUILDERS

  python3 data/build_exact_separator.py
  python3 data/build_ramified_pencil.py
  python3 data/build_function_and_equation_data.py

These recreate the explicit rational certificate, affine pencil and function/
equation specifications from literal data and formulas. They perform no
optimization. The repair certificate is supplied directly as exact rational
Gram matrices and facial congruences; verify_ramified.py checks every entry.

To rebuild the PDF, run pdflatex twice on proof/halfline_counterexample.tex
from the proof directory. Only standard LaTeX packages are used.

PROVENANCE

Public input commit: b2bfe64cd81a951b7bbe33a84b9dedd5d4ba162c
Repository: https://github.com/anubhav70747/q222-finite-sdp
The three files in reused/ are unchanged copies of the established small
source-validity inputs. No earlier lift matrices were regenerated.

Numerical calculations served only to discover candidate rational matrices.
Every result claimed here is certified by the supplied exact rational data,
explicit algebraic identities and independent standard-library verifiers.
