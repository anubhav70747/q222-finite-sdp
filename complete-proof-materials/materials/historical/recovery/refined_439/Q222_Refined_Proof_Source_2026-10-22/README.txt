Q222 — A REFINED COMPLETE PROOF
22 October 2026

Main document: Q222_Refined_Complete_Proof.pdf
Main source:   Q222_Refined_Complete_Proof.tex

BUILD

From this directory, with a standard TeX Live installation:

  pdflatex -interaction=nonstopmode -halt-on-error Q222_Refined_Complete_Proof.tex
  pdflatex -interaction=nonstopmode -halt-on-error Q222_Refined_Complete_Proof.tex
  pdflatex -interaction=nonstopmode -halt-on-error Q222_Refined_Complete_Proof.tex

The project uses standard LaTeX packages listed in its preamble. No original
archive, external image, computer algebra output, network call, or private
dependency is needed for compilation. The bibliography is embedded.

READING

Chapter 1 states the qualitative finite-lift theorem and gives its short
global argument. Chapters 2–7 prove its reductions and assembly. Appendices
A–G contain the complete exact certificates and their local hypotheses.
Equation tags such as RC or SC are local to the named certificate section;
hyperlinked section references identify their intended occurrence. The PDF
bookmarks include all certificate sections.

The certificate appendix is an essential part of the proof. The opening
alone is not a self-contained proof of the canonical certificate theorem.
No numerical SDP matrix-size bound or order-two block bound is asserted.

VERIFICATION

Read Q222_Verification_Record_2026-10-22.txt for the result, refinements,
scope and limitations. The checks/ directory contains 340 named targeted
exact checks and three supplemental Boolean support tables. These checks
corroborate the algebra; they are not a formal verification of the theorem.

The audit/ directory records original component identities and hashes,
source scoping, final reference rewrites, analytical review scopes, and
validation. Review notes describe intermediate states: the two small core
corrections and final dependency repairs they request are incorporated in
the delivered proof. They are retained as an audit trail, not open issues.

The complete recovered 798-page inputs were preserved. This refinement
removes repeated and superseded routes and adds new complete arguments;
it does not claim that the proof's length is globally optimal.
