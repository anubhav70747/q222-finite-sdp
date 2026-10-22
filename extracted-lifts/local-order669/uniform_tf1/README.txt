FIXED C1 FAMILY AND EVALUATED RATIONAL SOURCE PENCILS
23 October 2026

Result and scope
---------------
This directory removes the variable-M TF1 radicals from the entire lambda>0
bounded-positive-z C1 parent on a fixed testing rectangle. It derives one finite
algebraic family from the twelve principal roots and all retained event/circuit/
angular roots. Both independent state normals are retained. The uniform source
box proved in ../root/uniform_source_bounds.tex has

  0<x,sigma<=2^-30, 0<z<=1, native source inequalities, D,beta>=0, lambda>0.

No bound on H=h/sigma or positive lower bound on z,b,1-b,h,1-h is needed.
The test box is |u|<=2^-24, |w|<=2^-12, with all normalized real or complex
states. The exact native formulas are in ../c1_fixed_family and in the public
Uniform_C1_Certificate source. Literal D=0 and beta=0 are allowed. Literal zero
openings x,sigma,z and the other C1/C2 source families are not asserted covered.

This is NOT a full-Q222 lift. The actual affine pencils represent membership in
the displayed fixed local Gram cones exactly. Their source domains cover a
nontrivial continuous native source family; they do not prove global coverage.

Principal finite family
-----------------------
The fixed numerical bounds are delta=2^-12 and M=128. The principal 12x12 Gram
has exactly eight coefficient monomials in (sqrt(u),w), including the constant.
A fixed matrix reserve splits each of the seven nonconstant monomials. Finite
geometric sums reduce all test radicals to two gate variables. On the algebraic
cover

  u=delta^2*a^8, w=delta*btest^4,
  W=omega^4, |Z|=zeta^4,
  xi^2+a^4=1, eta^2+btest^4=1,
  omega^8+zeta^8+X^2+Y^2=1,

the principal functions, event roots, angular roots and circuit roots are
polynomials. All powers delta^(j/4) used for normalization are rational.

Actual lambda>0 full-parent matrices
-----------------------------------
C1_POS_FULL_SOURCE_PENCIL.coo.json:
  52 PSD blocks, total order268, maximum order12;
  7 canonical source coordinates,1087 auxiliaries;
  168 distinct functions,281 independent coefficient equations.

C1_NEG_FULL_SOURCE_PENCIL.coo.json:
  49 PSD blocks, total order240, maximum order11;
  7 canonical source coordinates,934 auxiliaries;
  143 distinct functions,221 independent coefficient equations.

C1_BOTH_SIGNS_LOCAL_SOURCE_PENCIL.coo.json:
  101 PSD blocks,total order508,maximum order12;
  source order(d,m,C,A,D,alpha,beta),2021 auxiliaries;
  5326 nonzero rational affine entries.

Every COO tuple is [block,row,column,variable_index,rational_coefficient]. Only
row<=column is stored; mirror off-diagonal entries with the same value. The first
seven variables are the source coordinates; all remaining variables are free
Gram auxiliaries. There is no constant term. Thus these are fully evaluated
homogeneous affine pencils. They require no source formula, SDP solve or
symbolic compilation to read as matrices.

C1_*_FULL_FIXED_FEATURES.json gives every testing monomial, each block, each
retained source vector and scalar weight. C1_*_FULL_COEFFICIENT_MAP.coo.json
prints the coefficient equations and the linear canonical-parent target map.
The block-diagonal pencil is an exact representation of the intersection of
the two fixed certificate cones. Feasible Grams imply local slack positivity;
all uniformly admitted native sources have explicit feasible source Grams.

Verification
------------
The coefficient maps and affine elimination use only Python's standard-library
Fraction arithmetic. Every output equation is independently rebuilt from the
final pencil. The positive and negative charts pass281+221 exact target-equation
checks. The source reassembly, radical signs, both state normals and all retained
root identities are checked independently in ../c1_fixed_family. No numerical
CHSH test is used for adequacy. No whole-Q222 reverse inclusion is claimed.

Reproduction (preserve the sibling directory layout)
----------------------------------------------------
  python3 compile_principal.py
  python3 eliminate_principal.py
  python3 compile_full_charts.py
  python3 combine_signs.py

Inputs are the JSON exports in ../c1_fixed_family. The programs do not import or
execute archive verification scripts. Outputs are deterministic exact rational
matrices, feature lists, coefficient maps, and check receipts. The first two
commands also expose a separate principal-only cone(order192,280 output
coefficient coordinates,968 auxiliaries), useful for checking the uniform
matrix split independently of the retained roots.

Mathematical proof
------------------
UNIFORM_TF1_PROOF.tex gives the fixed matrix reserve, explicit source guards,
finite monomial radical elimination, algebraic test cover, pencil construction,
and a precise explanation of why the pointwise TF1 argument alone does not
supply a finite global family. Complete parent and retained root formulas are in
../c1_fixed_family/c1_construction_appendix.tex. Uniform source constants are in
../root/uniform_source_bounds.tex.

Negative-state source congruence
--------------------------------
The matrix source_matrix_negative_state_sign in
../c1_fixed_family/C1_ROOT_MATRIX.json is denoted A_minus here. It sets the
twelfth source-root column (the common circuit centre G) to zero, but its
testing-span entries still use the actual signed state coordinate Z.
Consequently the substitution Z=-zeta^4 must also be included when constructing
a source witness in the negative chart. Let A_minus_11 be its first eleven
columns, and use the first eleven monomials of

  chi=(X,Y,zeta^4,omega^4,a^8*zeta^4,btest^4*zeta^4,
       a^8*omega^4,a^8*btest^4*omega^4,a^16*omega^4,
       a^12*omega^4,a^4*omega^4,a^4*btest^2*omega^2*zeta^2).

Then, with delta=2^-12, the complete normalized negative-chart source map is

  R_minus=A_minus_11 * diag(1,1,-1,1,-delta^2,-delta,delta^2,
                            delta^3,delta^4,delta^3,delta).

For every principal source PSD block H_m the associated fixed-feature Gram is
R_minus^T * H_m * R_minus. The positive chart uses all twelve columns and

  R_plus=A_plus * diag(1,1,1,1,delta^2,delta,delta^2,delta^3,
                      delta^4,delta^3,delta,delta^(3/2)).

In contrast, the negative-chart retained roots exported in
../c1_fixed_family/C1_RETAINED_ROOTS.json already include the signed-Z
substitution. Do not apply a second sign change to those coefficients. Their
quarter-angle monomials are normalized by u^(1/4)=2^-6*a^2 and
w^(1/4)=2^-3*btest. These conventions are already implemented in the compiler;
this paragraph clarifies how to reconstruct the source witnesses and does not
change any matrix entry.

Standalone final replay
-----------------------
After the sibling candidate_boundary and root exports are present, run

  python3 ../root/verify_final_pencils.py

This standard-library-only verifier reads the exported fixed features,
coefficient maps and final matrices. It independently rebuilds the feature
products, the complete seven-coordinate canonical parent, every affine
coefficient equation, exact ranks and free-coordinate parametrizations, the
two-sign intersection, and the final source-sum split with the structured
six-root cone. It writes ../root/FINAL_EXACT_REPLAY.json. The final combined
pencil has order669,120 blocks of maximum order12,7 source coordinates and
2694 auxiliaries. All682 final component target equations are checked directly
against the final combined pencil. No pencil or input export is modified.
