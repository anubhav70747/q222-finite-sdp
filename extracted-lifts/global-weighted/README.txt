Q222: a global certificate for the multiplied slack
==================================================

Source commit: fb58c232e374dadf3cbb5f23e7e7413d429de884
Repository: https://github.com/anubhav70747/q222-finite-sdp

WHAT IS PROVED

For the complete canonical Bell parent F, with d and all eight event
weights nonnegative, let H be its W,X,Y principal state block. Then

  F is quantum-valid
    if and only if
  (det H) F belongs to one explicitly evaluated fixed Gram cone.

The testing cover is compact and algebraic, independent of every source
coefficient. It includes all projective tests and zero event weights.
The state completion is a polynomial adjugate identity for all four
state coordinates. The fixed angular certificate uses the twelve
balanced circuits already present in the source geometry.

The indexed construction has 189 functions. Exact linear compression
gives 104 independent functions and 27 PSD blocks:

  22 blocks of order 4; four blocks of order 6; one block of order 21.

Total matrix order: 133. There are 535 symmetric Gram coordinates,
261 product coefficient equations of rank 189, and 346 free directions.
The eliminated rational pencil has 6079 upper-triangle coefficient
records. Its evaluated quartic source map has 7472 monomial records.

WHAT IS STILL OPEN

The 189 visible pencil coordinates are polynomial coefficients, not the
eight coordinates of Q222. Their map from the nine canonical Bell/source
coefficients has homogeneous degree FOUR. Thus the result is an exact
quartic pullback of a fixed spectrahedral shadow, not an affine lift of
Q222. Replacing those quartic monomials by independent variables would
lose the proved equivalence.

The first obstruction to the direct cancellation route is explicit:
even on a compact, uniformly strictly positive physical source/test
rectangle, exact inverse/whitening factors require infinitely many
linearly independent testing functions. This rules out that route,
not alternative unmultiplied SOS constructions or finite lifts.

MAIN DELIVERABLES

proof/global_weighted_certificate.pdf and .tex
  Seven-page complete proof of the multiplied-slack theorem, its two
  inclusions, the evaluated pencil, and the limited whitening obstruction.

weighted/WEIGHTED_FIXED_FUNCTIONS.json
  Actual 104 polynomial functions on the cover; reconstruction of all
  189 indexed functions; exact rational maps for each PSD block.

weighted/WEIGHTED_GRAM_PRODUCT_MAP.coo.json
  Evaluated integer 261-by-535 product matrix, with monomial rows and
  all off-diagonal factors specified.

weighted/WEIGHTED_CONE_ORDER133.coo.json
  Actual rational homogeneous pencil in 189 visible coefficient
  variables and 346 free variables. Upper-triangle entries are mirrored.

weighted/CANONICAL_SOURCE_QUARTIC_MAP.json
  The essential nonlinear source input to that pencil, fully expanded.

weighted/multiplier_schema.json
  Eight original rows, twelve integer circuit rays, twelve evaluated
  rational 21-by-21 PSD matrices, and positive dyad certificates.

weighted/DETERMINANT_SOURCE_COEFFICIENTS.json
  All nine expanded quartic determinant coefficients. Their complete
  identity was checked independently from the original affine parent.

REQUESTED BOUNDED-R BRIDGE

bridges/bounded_R/ proves the unchanged complete C1 construction for
every real R>=1, with B=1+R and

  epsilon_R = 2^(-40) B^(-8),   delta_R = 2^(-20) B^(-2).

The source domain is 0<x,sigma<=epsilon_R, 0<z<=R and all native gaps,
with D,beta>=0 and unrestricted H=h/sigma. Independent testing variables
satisfy sqrt(|u|),|w|<=delta_R. BE3--BE9, the e-bound, SP6 and both state
normals are unchanged. This gives the stated bounded-positive-z transfer
bridge. The testing family can depend on R; no R-independent family is
inferred. Literal zero source denominators are not assigned values.

The included rational negative minor is only an obstruction to keeping
the unchanged main Gram and restoring radial reserves. A separately
retained angular square can repair its scale; it is not a negative
physical slack or an impossibility result for fixed testing families.

COMPLEMENTARY ACTUAL-SOURCE ANNULUS

bridges/compact_Bob/ gives a six-function identity eliminating the
chosen real Bob residue from the leading completion. A complete
four-state reserve, source-only clipping, and a fixed algebraic grid
then give 65 functions per chart for the actual C1 parent on an
explicit complementary annulus. The source admits every compatible
radial opening and literal D=0,beta=0. The finite bounds on L and the
omitted double corner are stated explicitly. The exact negative defect
below the three-event wall explains why that cap cannot be extended
over the missing corner. These are supplements to the main global
multiplier result, not an asserted full-set completion.

REUSE OF THE PREVIOUS EXTRACTION

reused/C1_UNIFORM_LOCAL_ORDER669.coo.json is copied byte-for-byte from
the pinned public commit. Its SHA256 is in PROVENANCE.json. It was not
recompiled. The correction to the w^2 row bound, 5 in place of 4, is
accepted; the old bound128 and I/6 reserve are unchanged. The new R
admission is not identified with an unchanged old coefficient embedding.

REPRODUCTION

Python 3.10+ and SymPy 1.14.0 suffice. Install requirements.txt, then run:

  python reproduce.py

This checks the new circuit identities, all 189 original functions,
the compressed basis, every one of the 535 products using an independent
Groebner computation, every pencil coefficient, and the entire quartic
source identity in all nine source variables and all ten state
coefficients. It also checks both supplements. All arithmetic is exact.
The saved independent receipt includes hashes of the evaluated inputs.

To regenerate the NEW matrices before checking them:

  python reproduce.py --rebuild

The order-669 file is only hash-checked, in either mode.

For the main proof PDF, run pdflatex twice on
proof/global_weighted_certificate.tex. The supplement TEX files are
standalone. No external archive, solver, or network lookup is needed to
run the mathematical checks in this package.


PUBLICATION NOTE

The seven-page proof has been edited for exposition; its hypotheses,
identities and matrix data are unchanged. Both supplementary TEX notes now
have matching PDFs. The bounded-R domain was checked against the original
complete square formulas: sqrt(|u|)<=delta_R is essential, and the altered
gate scale does not establish admission to the old evaluated order-669
pencil. RECOVERED_SHA256SUMS.json describes the original received ZIP;
SHA256SUMS.json describes this publication. The original archive SHA256 and
the current file hashes are recorded in PROVENANCE.json.
