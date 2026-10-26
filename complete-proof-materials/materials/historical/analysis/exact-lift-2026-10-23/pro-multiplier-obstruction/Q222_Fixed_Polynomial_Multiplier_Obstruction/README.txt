Q222: FIXED POLYNOMIAL MULTIPLIER OBSTRUCTION

Proved result
-------------
For every fixed nonzero polynomial q(v)>=0 on the real line, the accepted
two-LMI canonical completion fails for all sufficiently small members of
the explicit rational source family in EXACT_SOURCE_AND_FUNCTIONS.json.
All nine source coefficients are strictly positive for 0<T<=1/100.
The sources are valid for all states, both signs of each measurement sine,
commuting tests and projective endpoints. The d=0 endpoint is retained.

More generally, no fixed product p(u)q(v) of nonzero globally nonnegative
polynomials makes every canonical valid slack polynomial SOS. The argument
is analytical and uses exact polynomial identities. It does not use an SDP
solution or an identification with NPA certificates.

The half-line quadratic selector remains unresolved. Its proposed C1 form
is necessary, with eight coefficients; reflection retains sign(v) unless
r0=0 is independently justified. The small rational Gram in the JSON file
certifies ONLY the critical polynomial M(x,y)+2 on y>=0. It is not a Gram
certificate for the full canonical slack or an affine lift of Q222.

Contents
--------
Q222_Fixed_Polynomial_Multiplier_Obstruction.pdf  Complete proof and scope.
source/fixed_polynomial_multiplier_obstruction.tex  Editable LaTeX.
EXACT_SOURCE_AND_FUNCTIONS.json  Nine source functions, 4x4 Bell matrix,
  four-component state, half-line selector functions and critical Gram.
CRITICAL_PULLBACK.txt  Every coefficient of the exact polynomial pullback.
check_fixed_q_obstruction.py  Exact identity and data verification.
verification.txt  Recorded verification output.
SHA256SUMS.txt  Checksums of the above files.

Reproduction
------------
Requires Python 3 and SymPy (checked with SymPy 1.14.0):

  python3 check_fixed_q_obstruction.py

The check verifies the canonical coefficient map, determinant identity,
cubic bound identity, all pullback coefficients, the machine-readable
matrix/function data, and the half-line reflection and Gram identity.
It uses exact rational arithmetic throughout.

To rebuild the PDF, run the following command twice from this directory:

  pdflatex -halt-on-error -interaction=nonstopmode -jobname=Q222_Fixed_Polynomial_Multiplier_Obstruction source/fixed_polynomial_multiplier_obstruction.tex

Scope and provenance
--------------------
Input repository commit: 918523fc022bb7d1b0bf5304fcd485342528456f.
The accepted two-LMI theorem, previous obstruction and order-49 Laurent
lift are used as inputs and are not regenerated. The 343-page companion
is unchanged. The critical Motzkin idea is credited to A. Chaturvedi,
arXiv:2607.14569v1. This package proves the new multiplier obstruction with
an explicit rational cubic upper bound and polynomial full-state test.

The threshold for a general prescribed q is proved to exist, but is not
numerically evaluated. The rational sequence T=1/n^2 eventually excludes
each prescribed q. No globally adequate algebraic family or affine pencil
is claimed. The exact unresolved half-line implication is displayed in
the final section of the proof.
