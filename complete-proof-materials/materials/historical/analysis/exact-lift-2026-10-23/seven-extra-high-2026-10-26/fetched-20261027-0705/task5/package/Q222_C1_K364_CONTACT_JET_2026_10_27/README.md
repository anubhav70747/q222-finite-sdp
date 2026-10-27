# Q222 oriented C1 K364: exact source-only contact-jet construction

This is a constructive, **contact-local** full-state PSD result and explicit first-order/second-order Schur mechanism, not universal K364 admission. The proof retains the genuine native C1 coquotients and full parent. No numerical SDP search, parameter sweep, global-matrix recompilation, or GitHub write was made.

Files:
- `PROOF.md`: complete mathematical derivation, full 4+4 decomposition, rank-6 Calderon contact anchor, source-only contact skew line, full Schur/contact-jet positivity theorem, exact unresolved global matrix inequality.
- `COORDINATOR_THEOREM.md`: short theorem and precise frontier.
- `source_only_contact_jet.py`: deterministic exact formula from native source and a source-algebraic admissible contact-gauge parameter tau to the quadratic Bob polynomial and its strict Schur margin.
- `verify_contact_geometry.py`: exact polynomial checks for full-parent parity, affine contact skew, tangent annihilation, opposite soft charge and two-dimensional source skew freedom.
- `verify_two_schur.py`: exact matrix congruences proving the full 4x4 parity blocks reduce to **two 2x2** Schur matrix inequalities.
- `verify_contact_jet_exact.py`: one exact rational regression using a known genuine strict native point. It verifies the local jet algorithm; the source-wide theorem is mathematical, not inferred from it.
- `verification.log`: combined exact replay output.
- `original/`: unchanged native source dictionary and prior verified full-parent Bob-harmonic proof/checker.
- `MANIFEST.json`: SHA-256 and byte counts of every included file.

The source-algebraic contact parameter tau is obtained by solving six **strict 3x3 leading-principal-minor inequalities in one scalar** after substituting the complete exact native weights. Existence throughout the genuine domain is proved analytically from Calderon + the contact annihilator and Hessian. Singular full finite-lift matrices elsewhere are not excluded; this theorem uses strict contact-complement Grams only as a local construction device.

External import: A. P. Calderon, *A note on biquadratic forms*, Linear Algebra Appl. 7 (1973), 175–177, DOI 10.1016/0024-3795(73)90052-9.
