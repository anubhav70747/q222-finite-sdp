Q222: exact canonical completion and polynomial-selector obstruction

RESULT
------
The full nine-source canonical slack has an exact reduction to two 4 x 4
matrix inequalities with one Bob-dependent 2 x 2 selector. Canonical validity
is exactly: for every Bob test v there exists a selector J satisfying both
inequalities. All four state coordinates, zero source weights, and projective
endpoints are retained.

For C1, every globally polynomial selector reduces necessarily to three
parameters. The already certified rational extreme source admits none.
Its exact contact forces a unique candidate with a principal minor at v=1
strictly less than -103/80. This also proves that the reduced polynomial
Phi_c(u^2,v;a,b) is not a polynomial sum of squares.

STATUS
------
No global affine lift is obtained. The order-8 data are POINTWISE pencils:
the universal Bob quantifier has not been eliminated. The obstruction does
not exclude a fixed algebraic testing cover, an additional testing-only
multiplier, or another completion of the original slack. It must not be
read as a multiplier-removal theorem.

CONTENTS
--------
Q222_Canonical_Completion.pdf
    Complete short proof and precise remaining implication.
source/canonical_completion.tex
    Editable LaTeX source.
POINTWISE_SELECTOR_ORDER8.coo.json
    Two 4 x 4 blocks; nine nonnegative source parameters and four free
    pointwise selector parameters; 108 full symmetric integer COO entries.
C1_POLYNOMIAL_SELECTOR3_ORDER8.coo.json
    The derived three-parameter polynomial restriction; 80 integer entries.
    The supplied valid source disproves its adequacy.
POLYNOMIAL_SELECTOR_COUNTEREXAMPLE.json
    Exact contact linear equations, unique selector, negative minor and
    factored determinant. All values are rational strings or integers.
ESTABLISHED_VALID_SOURCE.json
    Accepted source vector and contact, copied verbatim from public commit
    5fc0a0487bb94d7bdfa9929edb5f4f6f3ff50240, file
    extracted-lifts/outer-candidates/strict_outer_physical_source.json.
    Its physical validity is an established input from the independently
    verified order-13744 certificate. That matrix archive is not rebuilt.
check_selector_certificate.py
    Independent verifier of every new matrix coefficient, contact equations,
    uniqueness, and exact negative minor. Python standard library only.
check_alice_quadratic_reduction.py
    Symbolic verification of the complete nine-source polynomial identities.
build_selector_data.py
    Rebuilds only the new order-8 matrices and obstruction certificate.
verification.txt
    Recorded output of the two delivered verifiers.
SHA256SUMS.txt
    Checksums for all other package files.

REPRODUCTION
------------
From this directory, the main independent check requires only Python 3:

    python3 check_selector_certificate.py

With SymPy 1.14.0 installed, verify the full-state identities:

    python3 check_alice_quadratic_reduction.py

Optionally rebuild only the NEW data:

    python3 build_selector_data.py
    python3 check_selector_certificate.py

Compile the proof (TeX Live with standard AMS packages):

    pdflatex -jobname=Q222_Canonical_Completion source/canonical_completion.tex
    pdflatex -jobname=Q222_Canonical_Completion source/canonical_completion.tex

COO CONVENTION
--------------
An entry [k,j,r,s,a] means that a*parameter[k]*v^j is added to matrix
entry (r,s). Indices are zero-based. Both off-diagonal copies are stored;
do not symmetrize or double them again. Each JSON file declares its exact
source/selector order and the scope of its quantifiers.

The source and testing fields remain independent. The real-closed-field
selector theorem is a pointwise assertion in one field; it does not assert
that its selector belongs to the source/testing tensor product.
