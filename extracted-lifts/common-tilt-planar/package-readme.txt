Read proof/main.pdf for the six-page extended proof and explicit lift.
RESULTS_AND_BRANCHES.txt maps the completed and remaining implications.

This extends the supplied symmetric support-epigraph construction to the
full planar behavior projection (sum of marginals, CHSH), by allowing all
real homogeneous support coefficients and dualizing with attainment.
It also extends the asymmetric fixed-input construction to all tilt signs.

Reproduction (Python3 with SymPy):
  python verify_symmetric.py
  python verify_membership_and_asymmetric.py

Compile proof/main.tex twice with pdfLaTeX. No project-specific libraries,
numerical solver, external files, or downloaded local dependencies are needed.

The original archive is preserved separately. The proof retains the prior
NPA nontermination theorem as an imported result; exactness of the wider
eight-coordinate evaluated global pencil remains subject to the universal
C1/C2 admission targets identified in the branch map.
