Q222: exact obstruction to the ramified quartic completion
========================================================

Result: C_ram is a proper subcone of the full nonnegative canonical
valid-slack cone. The all-positive rational source

  (d,lambda,n,m,C,A,D,alpha,beta)
  = (1,1/576,161,2^-24,2^-34,2^-14,1/2,1/250,76)

is valid by an exact canonical symmetry of the already certified cbar,
and is separated from C_ram by the exported integer functional and four
strictly positive-definite integer Hankel matrices of orders 12,8,12,8.

This is a refutation of universal quartic coverage, not a refutation of
all finite lifts. The exact physical transpose moves the accepted Bob
quarter-power functions to Alice. The two certificate cones are
incomparable; global adequacy of their union or sum is not asserted.

Read the PDF or proof/ramified_obstruction.tex for the proof. Every one
of the ten integer 4x4 moment matrices is printed in the PDF appendix.
The complete machine-readable certificate is data/certificate.json.

Verify with Python 3 (standard library only):

  python3 verify.py

The verifier uses fraction-free positive leading minors for the dual
Hankel matrices; derives all nine functional coefficients directly from
the accepted P,N,E,K definitions; checks all 148 auxiliary contractions
against the unchanged public COO data; proves the physical transport
as a polynomial identity in homogeneous angle coordinates; and checks
the exact value and all six first derivatives at the transported contact.
No SDP solver, floating point, SymPy, or matrix regeneration is used.

Rebuild the PDF (pdflatex):

  cd proof
  pdflatex -interaction=nonstopmode -halt-on-error ramified_obstruction.tex
  pdflatex -interaction=nonstopmode -halt-on-error ramified_obstruction.tex

Established inputs reused unchanged:

  Commit 16b37f6deb198fa48134145b4d12b7cc7d37be36
  https://github.com/anubhav70747/q222-finite-sdp

  reused/ramified_order49_pencil.json is the accepted evaluated pencil.
  SHA256 5428ca18e35f2d2468c2560eba3109cf26b6c98edc1dcca9c7a108aef80d7421

  reused/ESTABLISHED_VALID_SOURCE.json is the accepted rational contact
  source; it is used only for the transported contact check. Its original
  validity proof and the admission of cbar in C_ram are accepted inputs.

All prior successful artifacts and the complete companion are preserved.
