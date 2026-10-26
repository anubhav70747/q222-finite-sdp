Q222 REGENERATED COMPILER COMPONENT -- 22 OCTOBER 2026

This folder contains two distinct mathematical components.

1. EXACT TILTED-CHSH PROJECTION EXTRACTION

The projection t=a0, s=E00+E01+E10-E11 is exactly
  T = {(t,s): |t|<=1, t^2+s^2/4<=2}.

tilted_family_exact_extraction.tex gives a complete analytical proof of both
inclusions, the operator upper-bound identity and explicit attaining models,
fixed source-independent testing functions, separately labelled tensor-field
certificates, and all tilt endpoints. The evaluated lift has N=5,m=1; its
eliminated version has N=4,m=0. Their complete exact coefficient files contain
14 and 12 upper-triangular nonzero records, respectively. All coefficients
are encoded as rational pairs in Q(sqrt(2)), with sqrt(2)>0.

Run:
  python3 build_and_check_tilted_extraction.py

This standard-library-only program freshly regenerates all TILTED_*.json
files and CHECK_TILTED_EXTRACTION.json. It checks 38 named exact algebraic
conditions. TILTED_MANIFEST.json binds the executable helper, checker and
generated data with SHA-256.

SCOPE: These pencils describe T, or its inverse image in the eight original
coordinates. They do not describe the full quantum set. The exact rational
point (0,0,0,0,1,1,1/4,-1/4) satisfies every probability constraint and every
standard relabeled single-marginal tilted-CHSH inequality but is nonquantum.
The proof and TILTED_SCOPE_COUNTEREXAMPLE.json contain the exact obstruction.

2. CONDITIONAL EFFECTIVE EXTRACTION AND UNSOLVED QE INPUT

effective_algebraic_extraction.tex proves the rational physical membership
formula using at most nine real two-qubit atoms, including all finite complex
and mixed quantum models and binary measurements. Its 81 quantified variables
satisfy 36 equations and 9 weak inequalities, of maximum degree 5.

The note proves conditional real-algebraic coefficient descent and an exact
search procedure if a finite lift of the searched type exists. It does not
freshly verify the historical full-Q222 finite-lift existence premise.
Arbitrary finite semidefinite-lift existence does not guarantee an order-2
block lift. Neither existence premise is asserted by this component.

Run:
  python3 generate_q222_qe_problem.py

The program generates Q222_SOC_M3_UNSOLVED_INPUT.wl and its exact JSON
counterpart. The fixed three-block candidate has 162 unknown coefficients
and 9 auxiliary variables. CHECK_Q222_QE_ENCODING.json records 20 fresh
encoding checks. QE_ENCODING_MANIFEST.json binds the exact input files.

No quantifier elimination, feasibility decision, or coefficient sampling has
been run. Loading the Wolfram Language file only constructs the formulas.
The candidate remains UNSOLVED. Global extraction criteria E1--E3 remain OPEN.

REPRODUCIBILITY AND REVIEW SCOPE

Both programs use Python >=3.9 and the Python standard library only. Exact
arithmetic is implemented with integers and fractions.Fraction; there are no
floating-point checks or external solver dependencies. The checker receipts
describe finite algebraic and encoding checks, not proof-assistant formal
verification. The full analytical inequalities and field/quantum arguments
are in the two TeX documents. Source citations include direct primary-source
URLs. Bamps--Pironio's operator SOS is Eq. (27), not Eq. (25).

This is newly regenerated work. It is not a byte-identical restoration of an
unavailable historical 549-page R21 archive, and no old PASS is carried over.
