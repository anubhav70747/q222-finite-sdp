Q222 ANALYTICAL EXTRACTION REFINEMENT - REGENERATED 2026-10-22

Start with the report PDF, then CHECKLIST.txt. The full original sources,
coefficient data and runnable exact checkers are included in this archive.

WHAT WAS REGENERATED
* Complete-parent regular-C2 local fixed features and rational product maps.
* Exact tilted-CHSH projection, all tilt parameters, original-coordinate
  projection-cylinder coefficients, and a nonquantum scope witness.
* Separated-coefficient and denominator-free smooth-diagonal local theorems.
* The displayed B0 coefficient module, with the frozen factor-typing repair.
* Restricted event-contact cofactor data.
* The literal Schur-remainder obstruction and its limited scope.
* A rational full-quantum membership formula and a conditional exact
  coefficient-search specification. Its QE candidate is UNSOLVED.

WHAT IS STILL OPEN
E1: a complete evaluated global feature list and actual global dimensions.
E2: an evaluated affine pencil exactly projecting to the full quantum set.
E3: both full-set inclusions for that same fixed pencil.
The required restoration, atlas and coverage data are itemized in
GLOBAL_EXTRACTION_OBLIGATIONS.json. No new historical-existence certification,
global matrix-size bound, or full-set order-two representation is asserted.

DIRECTORIES
inputs/       Hash-verified retained source inputs and selected prior context.
recovered/    The unchanged earlier continuation ZIP (including its sources).
c2/           Local regular-C2 extraction and literal Schur-remainder audit.
diagonal/     Local theorems, B0 and restricted contact data.
compiler/     Exact tilted coefficients and UNSOLVED global QE input.
publication/  Combined report source and documented source transformations.
verification/ Fresh isolated replay, cross-review and document QA records.

VERIFY
  python3 verify.py
checks every manifest-bound file and the retained-source hashes.
  python3 -m pip install -r requirements.txt
  python3 verify.py --replay quick
also runs six exact construction/verification suites and three separate review
suites in an isolated temporary copy, and requires
every manifest-bound regenerated output to remain byte-identical. The local
C2 and tilted/QE suites use the standard library; diagonal and Schur suites
require SymPy 1.14.0. No network access is used by the checkers themselves.

REBUILD THE REPORT
  python3 build_release_content.py
  cd publication
  pdflatex -interaction=nonstopmode -halt-on-error main.tex
  pdflatex -interaction=nonstopmode -halt-on-error main.tex
The source uses standard LaTeX mathematical and typography packages. The
delivered PDF was rendered and inspected. A rebuilt PDF may differ in binary
metadata even when its mathematical content is identical.

This is a newly regenerated refinement, not an exact byte restoration of the
unavailable prior refinement or a reproduction of the historical 549 pages.
All previous releases are preserved as source artifacts; none is overwritten.
