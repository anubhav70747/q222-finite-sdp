FIXED TESTING FUNCTIONS FOR THE FULL-R1 TESTING CAP

Scientific scope
The rational affine system defines a closed, evaluated certificate cone
K283 in source order (d,lambda,n,m,C,A,D,alpha,beta). It is sound for all
states and signed testing angles |u|<=2^-12, |v|<=2^-25. It admits the
native source box 0<B,T,jB<=2^-10, 0<p<min(j,1), 0<H<1, n>=0, with
the actual rational dictionary in the accompanying source appendix;
there is no bound on j or jT. It also admits the eight event rays and
their nonnegative sums with admitted sources. This is not a certificate
for all testing angles and is not an exact global Q222 lift.

Files
fixed_family.json: explicit 49 monomial vectors, rational polynomial
weights, chart assignment, source coefficient conventions and fixed grid.
root_catalog.json: every certificate row/support (306 records), constant
G7/G3 blocks and exact four scalar loss polynomials. Source symbols in a
row are coefficients, not testing variables or function arguments.
affine_system.json: 554 rational coefficient equations, 2707 Gram records
and 132 source records; 2269 actual upper-triangular Gram variables.
Offdiagonal symmetric entries contribute twice. All indices start at0.
The 49 Gram block orders sum to274; nine source-sign scalar blocks give
total PSD order283. The two chart orders are143 and131.

Reproduction (Python3 with SymPy)
  python build_support.py
  python compile_affine.py
  python check_factored_assembly.py
  python check_fixed_grid.py
Each script reads/writes only this directory. No optimizer, numerical
sampling, received pickle or search over degree is used. The grid checker
uses only the Python standard library. The factored assembly avoids an
unnecessary single expansion of the entire source/state expression.

Exact logical dependencies
The small primitive circuit identities, both W19 sign charts, permuted
W19, and every testing support are checked by build_support.py. The eight
event rays and correct rotated target columns are checked independently
by compile_affine.py. The source-loss vectors and simultaneous entire
scalar residual are checked by check_factored_assembly.py. Its full
native-parent interface is the independently verified identity recorded
in ../normal_core_checks.json and ../full_r1_enlarged_gram.py.
The latter's current private imports are explicitly:
  ../../contact-cap-independent/weak-recovery/received_generator.py
  ../../contact-cap-independent/reconstruct.py (covariance prefix only).
For portable staging the covariance check can instead use the nine
explicit target columns in compile_affine.py, as all four sign actions
preserve these polynomials. The generator dependency supplies the actual
native rational source dictionary; none of the four scripts above imports
it. The exact source/Schur and four simultaneous budget inequalities are
in ../normal_core_checks.json, ../normal_and_ru_budgets.txt and
../check_normal_and_ru_budgets.py/.json.

Human-readable proof
../proof/Fixed_Testing_Functions.tex defines the support prescription,
affine equations, sign coverage and positive-grid closedness argument.
It is included by ../proof/Uniform_Contact_Certificate.tex, whose main
text and other appendices provide the source formulas, parent/core and
uniform inequalities. This file deliberately gives actual order283;
no unavailable announced order299 data is assumed or reconstructed.

Closedness
Every block has distinct monomials with individual degrees at most
(4,3,1,1,1,1). The explicit320-point grid per chart is Cartesian with
sizes(5,4,2,2,2,2), and every block weight is strictly positive there.
Interpolation gives a positive definite rational evaluation form for
each block. Nonnegative evaluated contributions bound every Gram trace
for a convergent source sequence; a bounded subsequence preserves all
encoded coefficient equations. No mixed-field inverse is used.
