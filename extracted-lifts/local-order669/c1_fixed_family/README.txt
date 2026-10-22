Q222: exact twelve-root C1 extraction
23 October 2026

Authority: anubhav70747/q222-finite-sdp, commit
0c1ce2e15ed067a6845464e2468559cea5274f2c,
latex/Uniform_C1_Certificate.tex. Inputs were read without modification.

This directory supplies the exact algebraic extraction from the bounded-z,
arbitrary-H C1 certificate. It does not claim a full-set Q222 affine lift.

DELIVERED IDENTITIES

* C1_ROOT_MATRIX.json: source-only 12 by 12 matrices expressing the twelve
  principal roots in the user's proposed common testing span. Both state-sign
  charts are printed. There are 24 nonzero entries on the positive chart.
* C1_PRINCIPAL_GRAM.json: eight exact symmetric source coefficient matrices
  in test monomials 1,w,w^2,s,s^2,s^2*w,s^4,s^4*w, s=sqrt(u). Upper-triangle
  storage has 104 entries. Off-diagonal numbers are matrix entries, not doubled
  polynomial coefficients. The original full four-state parent is preserved.
* C1_RETAINED_ROOTS.json: every elementary/angular/circuit coefficient matrix
  on both state-sign charts. Positive: 33 fixed retained testing monomials,
  30 roots. Negative: 27 monomials, 31 roots. Each root has separation rank
  at most three. Circuit roots are printed before multiplication by qC.
* source_definitions.json: algebraic/rational source definitions used by the
  matrices. Native p,l,kappa,C,A,alpha,D,beta are defined by P1 in the authority.
* c1_construction_appendix.tex: complete human-readable derivation, principal
  Gram, all nineteen elementary roots, eight/eleven angular roots, three/one
  circuit roots, exact boundary and sign routing, and fixed circuit multiplier.

INDEPENDENT VERIFICATION

All programs below use exact SymPy algebra; version used was 1.14.0. A normal
Python environment with SymPy installed is sufficient. No archived Python
packages or upstream checker code are imported.

  python3 build_core.py
  python3 build_retained.py
  python3 verify_native_parent.py
  python3 verify_reassembly.py
  python3 verify_retained.py

The native-parent check verifies all 17 nonzero pre-source polynomial
coefficients after the generic rational native P1 substitution. This is an
identity in every state coordinate and both independent testing angles, not
a finite collection of evaluations.

The reassembly check verifies all 12 source-matrix rows, both full state-normal
conductors, and the entire twelve-root Gram plus retained expression against
SP6. The retained verifier reads the actual exported matrices and squares
every root on both sign charts, checking each complete elementary polynomial,
angular certificate, and circuit remainder. The receipts are respectively
NATIVE_PARENT_VERIFIED.json, REASSEMBLY_VERIFIED.json, and
RETAINED_ROOTS_VERIFIED.json.

The later fixed-testing compiler combines these data with explicit uniform
source bounds; that compiler and its evaluated rational pencil are separate
deliverables. The source bounds and compiler are mathematical obligations,
not assumptions inferred from these algebraic identity checks. None of the
verification above is a proof-assistant formalization.

BOUNDARIES AND SCOPE

No state coordinate is inverted in the final formulas. The testing axes and
W=0/Z=0 are polynomial specializations, while D=0, beta=0, and c=0 are legal
source specializations without cancellation. The lambda<=0 source branch
has its own six-root construction supplied separately. Literal x=0, sigma=0,
z=0, and all source/test charts outside the explicit finite box are not added
by taking limits; global fixed-family coverage remains a separate task.
