FIXED FINITE C1 FAMILY AND EVALUATED LOCAL PENCIL
23 October 2026

WHAT IS PROVED

The complete stationary C1 reciprocal parent has a fixed finite algebraic
SOS certificate throughout the native source domain

  d=1, 0<x,sigma<=2^-30, 0<z<=1,
  0<b,h<1, r=b*z<1, E>0, M>0, J>0, D>=0, beta>=0,

with no upper bound on H=h/sigma and with either sign of lambda.

The independent testing domain is

  abs(u)<=2^-24, abs(w)<=2^-12, all four real state coordinates.

Complex states follow by real and imaginary splitting. The permitted D=0 and
beta=0 source faces, lambda=0, testing axes, and zero state coordinates are
included. Literal x=0, sigma=0, z=0, other source regimes, and C2 are not
claimed covered.

The positive-lambda construction uses exactly the recovered twelve-root
identity and all retained elementary, angular, and circuit roots. The Gram
has only seven nonconstant testing monomials. The numerical source envelope
128 and delta=2^-12 eliminate point-selected TF1 constants. A finite geometric
sum reduces all the principal gate radicals to two fixed algebraic variables.

The complementary six-root identity uses fourteen fixed radicals with
constant 1, and its uniform nonpositive-lambda source theorem holds on the
larger x,sigma<=2^-10 box.

ACTUAL MATRIX

  root/C1_UNIFORM_LOCAL_ORDER669.coo.json

This is an evaluated rational affine pencil representing exactly

  K = (K_POS intersection K_NEG) + K_SIX.

Its visible coordinates are (d,m,C,A,D,alpha,beta). POS and NEG refer to the
two state-sign charts, not to the sign of lambda. The first summand comes
from the positive-lambda proof; K_SIX comes from the six-root proof.

  Visible source coordinates: 7
  Auxiliary variables:       2694
  PSD blocks:                120
  Total matrix order:        669
  Largest block order:       12
  Nonzero rational records:  7716

The exact block formula is

  L669(c,z,y_plus,y_six) = diag(L508(z,y_plus), L161(c-z,y_six)).

The source split z contributes seven of the auxiliaries. The 508-order
intersection has 2021 auxiliaries; the 161-order component has 666.

Each COO record [block,i,j,variable,rational] specifies an actual upper
triangle matrix coefficient. All indices are zero-based. Mirror off-diagonal
entries exactly once. Visible variables are 0..6; remaining variables follow
the file's auxiliary_order. The constant matrix is zero.

No source substitution or symbolic operation is required to read this file
as the rational matrices A_j and B_k of L(c,y)=sum(c_j A_j)+sum(y_k B_k).

The positive source construction and exact elimination prove both directions
for THIS FIXED CERTIFICATE CONE. The cone is not asserted equal to all locally
nonnegative sources. This is not a primal lift of full eight-coordinate Q222,
and its feasibility does not establish global Bell validity. The uncovered
source and testing charts still require a uniform fixed-family adequacy proof.

PROOF AND FUNCTIONS

  proof/Fixed_Family_Extraction.pdf
  proof/Fixed_Family_Extraction.tex and its six input fragments

The 17-page proof is self-contained: source formulas, full parent identity,
all twelve principal roots, both state-normal conductors, all retained
squares, explicit source estimates, fixed radicals, six-root coefficients,
both local cone inclusions, and the remaining global obligation.

The exact feature functions and their product maps are in:

  uniform_tf1/C1_POS_FULL_FIXED_FEATURES.json
  uniform_tf1/C1_NEG_FULL_FIXED_FEATURES.json
  uniform_tf1/C1_POS_FULL_COEFFICIENT_MAP.coo.json
  uniform_tf1/C1_NEG_FULL_COEFFICIENT_MAP.coo.json
  candidate_boundary/six_root_general_formulas.json
  candidate_boundary/six_root_structured_product_map.json

The original unscaled twelve-root matrix, eight principal coefficient
matrices, and complete retained roots are in c1_fixed_family/*.json.
The native five-parameter rational source is specified in
root/NATIVE_SOURCE_MODEL.json and in Section 2 of the proof.

The compact positive-state algebraic cover uses
  u=delta^2*A_t^8, w=delta*B_t^4, W=omega^4, Z=zeta^4,
  xi^2+A_t^4=1, eta^2+B_t^4=1,
  omega^8+zeta^8+X^2+Y^2=1.
The negative chart uses Z=-zeta^4. Its principal source map includes the
explicit three column sign changes printed in the proof and in
uniform_tf1/README.txt. Retained negative-root coefficients already include
that sign substitution. The positive pure-testing normalization is
kappa_t^2*(1+u^2)*(1+w^2)=1.

VERIFY THE DELIVERED MATRICES

From this bundle's root directory run:

  python3 root/verify_final_pencils.py

Only the Python standard library is required. This verifier imports no
generator and no archived checker. It reads the delivered pencils without
modifying them and writes root/FINAL_EXACT_REPLAY.json.

It independently rebuilds all feature products and the full canonical parent,
checks the final coefficient maps, recomputes exact ranks 281,221,180, verifies
all 2687 component free Gram coordinates, and checks the final source split
through all 682 component coefficient equations. The seven source-splitting
coordinates bring the final auxiliary count to 2694.

A separate exact rational six-root fixture can be replayed with:

  python3 candidate_boundary/verify_six_root_pencil.py

The fixture is in the general guard domain and has lambda>0. It is not offered
as an example from the lambda<=0 subdomain. It contains 331 rational
upper-triangle source-Gram entries.

REGENERATE FROM THE PRINTED IDENTITIES

Install the optional symbolic dependency in a chosen Python environment:

  python3 -m pip install -r requirements.txt

Then run, in this order:

  python3 c1_fixed_family/build_core.py
  python3 c1_fixed_family/build_retained.py
  python3 c1_fixed_family/verify_native_parent.py
  python3 c1_fixed_family/verify_reassembly.py
  python3 c1_fixed_family/verify_retained.py
  python3 uniform_tf1/compile_principal.py
  python3 uniform_tf1/eliminate_principal.py
  python3 uniform_tf1/compile_full_charts.py
  python3 uniform_tf1/combine_signs.py
  python3 candidate_boundary/compile_six_root_chart.py
  python3 candidate_boundary/compile_six_root_pencil.py --structured
  python3 root/compile_combined.py
  python3 root/verify_final_pencils.py

The symbolic programs were checked with SymPy 1.14.0. The fixed-family
coefficient compilers and final replay use exact rational arithmetic.
The six-root generator also contains an optional merged 157-order alternative;
the PRIMARY delivered and combined cone is the structured 161-order one.
No numerical SDP, floating-point positivity check, or CHSH test is used.

To rebuild the PDF, run pdflatex twice in proof/:

  pdflatex -interaction=nonstopmode -halt-on-error Fixed_Family_Extraction.tex
  pdflatex -interaction=nonstopmode -halt-on-error Fixed_Family_Extraction.tex

PROVENANCE

Repository: https://github.com/anubhav70747/q222-finite-sdp
Commit:     0c1ce2e15ed067a6845464e2468559cea5274f2c

The exact current seven-page source is preserved as
reference/Uniform_C1_Certificate.tex. Input hashes are in
root/INPUT_PROVENANCE.json. These source identities, rather than a guessed
larger Gram frame, were the basis of the extraction.

The publication file hashes and recovered-source hashes are recorded in
PROVENANCE.json. README.md describes the corrected intermediate row bound.
The mathematical checks are the separate exact verifiers above.

