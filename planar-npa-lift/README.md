# Infinite NPA depth at minimal semidefinite complexity

The full quantum region of the CHSH value and the marginal sum
\(M=\langle A_0\rangle+\langle B_0\rangle\) has an explicit affine lift with
PSD blocks of orders 5 and 4. Every fixed NPA level projects strictly
outside that same region. Its semidefinite extension degree is 2.

Equal measurement overlaps determine every supporting line. At fixed
product of the two overlap parameters, the determinant is affine in their
sum; the arithmetic–geometric mean inequality selects equal overlaps in
the nonclassical directions. The remaining directions have local bounds.
One univariate quartic matrix then gives a fixed finite lift. Actual states
supply its membership variables directly, including at the endpoints.

Every point of this region is a mixture of at most **two real pure
equal-overlap two-qubit strategies**. Both membership matrices can
simultaneously have rank at most two; at an extreme point, at most one.
These are ranks of available witnesses. The explicit block orders remain
5 and 4; the optimal extension degree 2 follows separately from the planar
representation theorem.

- [PRL Letter](article/main.pdf): three pages including references.
- [Two-page note](note/main.pdf): the central insight and proof.
- [Complete proofs](proofs/main.pdf): 16 pages, including every membership
  matrix, both inclusions, the two-strategy realization, all signed support
  branches and the entire earlier NPA argument.
- [Complete Overleaf package](NPA_Finite_Lift_Overleaf.zip).
- [Accepted abstract](ABSTRACT.txt), reproduced verbatim in the Letter.
- [Exact data for all 32 projections](checks/ALL32_EXACT_PROJECTION_CONSTRAINTS.json).

## Open in Overleaf

Upload the ZIP and compile **main.tex** for the Letter. Select
**supplement.tex** for the complete proofs or **note.tex** for the short
note. Use pdfLaTeX. Bibliographies are included in the sources; no external
figures, bibliography files or private source paths are needed. The three
documents in article/, proofs/ and note/ also compile independently.

## What the proof covers

| Result | Exact scope | Complete proof |
|---|---|---|
| Common-tilt quantum region | Every real support direction, all boundary points; one 5/4-block affine membership lift | F1–F2 |
| Physical realizations | At most two real pure equal-overlap two-qubit strategies per planar point; simultaneous witness ranks at most two | F2 |
| Minimal block order | Extension degree 2 by the planar theorem and nonpolyhedrality | F3 |
| Independently signed tilts | Exact fixed-input support SDP for every real pair of tilts, including axes, critical and commuting branches | F4 |
| Relabelled projections | All 32 marginal–CHSH planes; exact support after intersection with any NPA level | F5 |
| Infinite NPA depth | Divergent relative support error near the local threshold at every fixed level; every fixed finite word list | S1–S5 |

The signed two-tilt support construction depends on the chosen tilts. The
full eight-coordinate quantum-set lift and universal C1/C2 admission remain
a separate project; see the repository [frontier](https://github.com/anubhav70747/q222-finite-sdp/blob/main/FRONTIER.md).

## Reproduce

Run `python build.py` to build all three PDFs using pdfLaTeX and RevTeX 4.2.
Intermediates are written to build/. To check the exact algebra, install
SymPy and run:

    python checks/verify_symmetric.py
    python checks/verify_membership_and_asymmetric.py
    python checks/verify_state_witness.py
    python checks/build_orbit_constraints.py

These check the 17 original polynomial identities; the adjoint matrices,
strict anchor and signed formulas; the physical rank-one witnesses,
normalization, endpoints and equal-angle Motzkin factorization; and the
32 projections (64 PSD blocks, total order 288, 512 auxiliaries and 96
affine equations). The universal positivity and planar chord arguments
are proved analytically in the companion.

[Verification record](VERIFICATION.json) · [Source provenance](PROVENANCE.md)
