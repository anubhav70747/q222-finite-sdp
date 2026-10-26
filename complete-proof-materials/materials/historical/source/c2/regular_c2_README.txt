REGULAR-C2 EXACT LOCAL EXTRACTION -- REGENERATED 2026-10-22

Scope
-----
This is a complete finite coefficient extraction for the declared local
normalized regular-C2 contact graph and independent testing box. It is not a
global lift of the full eight-coordinate quantum set. Global E1-E3 remain open.

The complete seven-source parent, all exact domain conditions, and the full
analytic proof are printed in regular_c2_local_proof.tex and its seven-page PDF.
The source order is (d,n,m,C,A,D,alpha), with d=Z=1 for the local certificate.
The source anchor is (1,7/8,9/8,1/18,11/8,9/8,1/18).
The contact anchor, in (u,v,W,X,Y) order, is (1,1,1,-1/10,1/10).
The exact radius is 510101005/350806159660744.

Reproduce
---------
Python 3.10 or later, standard library only; no numerical package is needed.
Run without -O, because assertions are part of the verification.

  python regular_c2_build.py
  python regular_c2_check.py

The first command deterministically recreates the three coefficient JSON files
and a builder verification receipt. The second uses a separately implemented
Fraction dictionary engine, imports no builder code, independently reconstructs
the full parent, and writes regular_c2_independent_verification.json.

If the sibling diagonal/ and inputs/ directories are present, the optional
cross-review corroboration runs as follows:

  python regular_c2_cross_review.py

It imports only arithmetic from the independent local checker, not the diagonal
generator. Its analytical review findings are distinguished explicitly from the
exact coefficient checks in regular_c2_diagonal_cross_review.json.

Recompile the standalone proof with two passes:

  pdflatex -interaction=nonstopmode -halt-on-error regular_c2_local_proof.tex
  pdflatex -interaction=nonstopmode -halt-on-error regular_c2_local_proof.tex

Complete exact data
-------------------
regular_c2_frames.json
  Six contact-equation polynomials, the contact Jacobian and nonzero minor,
  fifteen upper-triangular Taylor entries before and after recentering, H0,
  rational LDL data, all 25 frame polynomials and fixed linear forms, all positive
  shifted factors/remainders, and all 785 block labels/radicands.

regular_c2_feature_list.json
  All 1570 fixed features, including repetitions, each written as the square
  root of a rational testing polynomial times a rational polynomial multiplier.
  All 785 source Gram matrices are given exactly as U*v*v^T. The factorization
  retains exact polynomial source weights and vectors; it is not numerical.

regular_c2_product_map.json
  Complete 237-monomial indexing, 27747 nonzero symmetric product records, and
  all 295 nonzero records for the full seven-source target. Upper off-diagonal
  entries carry multiplier 2 in coefficient matching. The compiler has 2355 scalar
  upper-triangular Gram variables and 237 coefficient equations.

Record conventions
------------------
All rational coefficients are reduced strings. All exponent lists are ordered.
Source polynomials use 11 centered variables: the six source coordinates followed
by the five contact displacements. Testing polynomials use the five testing
displacements from the fixed anchor. The H_delta polynomials use the same first
eleven variables followed by five actual differences z-theta. H_centered and
gamma polynomials use the five independently centered testing variables instead.

The 785 blocks consist of 735 mixed positive products and 50 pure remainders.
The 237-monomial number is an ambient support count, not a claim of independence
of the listed features or of full rank of their product map. Maximum product
degree is 6.

Verification limits
-------------------
Fresh exact reconstruction checks all mathematical identities as complete
coefficient dictionaries. The final Taylor remainder has 1193 nonzero monomials.
All 50 strict remainder bounds and all 735 mixed-factor sign bounds are checked.
The original contact graph is retained, and its source Jacobian is nonsingular
throughout the declared contact box by an exact row-norm bound below 1.

Receipts distinguish encoding checks from mathematical identity groups. Counts
are overlapping checks, not counts of independent theorems. The independent
implementation and separate-agent review come from the same model-assisted
system. They are neither external peer review nor proof-assistant verification.
No historical PASS status is used in place of the current exact replay.
