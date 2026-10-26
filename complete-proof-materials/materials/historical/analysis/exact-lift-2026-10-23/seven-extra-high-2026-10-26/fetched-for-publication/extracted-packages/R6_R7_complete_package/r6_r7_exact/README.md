# Exact R6/R7 native C2 coefficient audit

- `R6_R7_COEFFICIENT_PROOF.md`: complete self-contained source bounds; scalar, both normal and hard formulas; bounded conductor with quantified divisibilities; residues and projective joins; scope and noncircular imported identities.
- `verify_arrays.py`: 12 exact whole-parent coefficient identities for independent state normals.
- `verify_fast_divisibility.py`: direct polynomial proof of the `F0/b` and `AN/(bs)` coquotients and their constant residues.
- `check_jets.py`, `check_b_unit_jets.py`: exact native source expansion coefficients for the two R6 scales and R7.
- `verify_source_jets_table.py`, `verify_unit_b_RU.py`, `check_RB17_coeffs.py`: low-weight full scalar coefficient matches from the source, for R6/R7 and both b types.
- `check_conductor.py`: RI14--RI15 polynomial identities and contact-gradient determinant.
- `check_remote.py`: original full parent determinant on the finite remote `M=0` kernel (RB24).
- `verification.log`: saved outputs from completed exact checks. `check_RB17_coeffs.py`, `check_conductor.py` and `check_remote.py` appear after the first interrupted batch and were separately completed exactly.

These tests use symbolic algebra only. No SDP solver, numerical sweep, Monte Carlo run or approximate positivity inference was used. The proof relies on the independently checked CA/SE local conductor, the strict bounded-Gram tensor lemma WA10, and listed exact BE/projective selected-frame polynomial identities, not on any global Q222 theorem.

For reproduction, run each script with Python 3 and SymPy. Some source-jet checks are intentionally separate to avoid constructing a giant universal expansion.
