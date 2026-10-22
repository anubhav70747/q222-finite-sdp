# Finite semidefinite representation of the simplest quantum correlation set

Manuscript and analytical proofs, 23 October 2026.

- [Article](pdf/Q222_Paper.pdf): 9 pages.
- [Complete proofs](pdf/Q222_Complete_Proofs.pdf): 348 pages.
- [Uniform C1 certificate](pdf/Uniform_C1_Certificate.pdf): a self-contained reading copy of the new local lemma, including the tensor sum-of-squares construction.
- [Overleaf source package](Q222_Overleaf.zip).

The new C1 argument combines two nonnegative source weights whose sum is a
positive unit. It controls both state-normal terms and the angular term
through the same allocation, including the permitted zero weights. It
replaces eight sections of the preceding proof. The remaining singular
source regimes and complementary measurement charts are retained in full.

## Build

The editable sources are in `latex/`. Compile `main.tex`, `Supplement.tex`
or `Uniform_C1_Certificate.tex` with pdfLaTeX, twice to resolve references.
The bibliography is included in the source. No numerical solver is required
for the written derivations.
