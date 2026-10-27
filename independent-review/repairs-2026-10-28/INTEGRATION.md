# Integration of the R5.5 and R6.5–R6.7 dependency repair

Baseline: `d01e5e089c5b4ad73667e922eec5b62b5d82485d`.

The compressed C2 companion retained the component labels but omitted the proof bodies. Restoring a reference alone would not resolve the issue: the original R5 source theorem uses different source hypotheses, and the ordinary R6.5 state matrix ceases to be strictly positive at the `(infinity,0)` corner when `n0=0`.

The integrated [proof](../../latex/certificates/c2_remote_dependencies.tex) supplies two subsections at the original destinations:

- **C.17.5:** the whole-parent centered-event identity R5.5, the complete WA matrix identity, a real Schur-cap lemma, the actual unbounded-a R6 source bounds, the full finite-remote kernel form and its positive-definiteness proof. The remainder is explicitly in the bounded separated tensor ring, so the existing Lemma 3.2 applies.
- **C.17.6:** the three exact homogeneous parents and R6.5–R6.7, the full centered corrections, and a common endpoint construction retaining all positive event excesses. The exceptional `(infinity,0)` case has its own complete polynomial identity in four roots. Real caps cover the stated R6/R7 regimes even when the actual A or D coefficient is unbounded.

The source-field and testing-field factors stay separate. No mixed inverse, mixed square root, source/testing scale comparison, or theorem transfer across different source domains is used. Positive determinants are accompanied by proofs of full-kernel or full-matrix positive definiteness and by the bounded tensor-Gram step. Literal reciprocal zeros are handled by the polynomial charts themselves.

The [R5.5 review](r55/REPORT.md) and [endpoint review](r65-r67/REPORT.md) retain their original bounded scopes. Their exact checkers verify 29 and 21 algebraic identities respectively; the analytical positivity and domain arguments are in the written proof. The integration changes only their placement, cross-references, and minor prose, with the original reviewer proposals retained alongside the checkers.

## Delivery verification

- The current `Supplement.tex` compiled with local pdfLaTeX through stable cross-references: 347 PDF pages. Final log has no undefined references, multiply defined labels, or overfull boxes. Existing PDF-bookmark warnings from mathematical headings remain.
- The nine immediate dependency and Gram-lemma labels each resolve uniquely. The restored text is on numbered proof pages 339–343 (PDF pages 342–346).
- All restored pages and the final reference page were rendered and inspected for clipping, missing glyphs, and equation layout.
- `Q222_Overleaf.zip` contains the new fragment and agrees byte-for-byte with the corresponding current source files. Its existing README is preserved. The old archive also lagged behind earlier tracked C1/C2 source edits; refreshing it introduces no new C1 source revision in this commit.
- A fresh extraction of the archive compiled independently to 347 pages, with text identical on every page to the published PDF and stable references.
- The built-in compiler encountered a local MiKTeX access error and an incompatible fallback engine. The successful builds above used the installed pdfLaTeX executable directly.

These checks establish delivery consistency, not mathematical truth. The written repair closes the named missing dependencies under the stated R6/R7 hypotheses. It is not a new audit of every contact conductor, all 49 C2 leaves, or the global existence theorem. Earlier reviews and original component archives are preserved as dated evidence.

Universal H1/H2 admission into the fixed explicit certificate cone is still unproved. Consequently reverse inclusion and exactness of that global SDP candidate remain open; this repair does not change them or the separate planar result.
