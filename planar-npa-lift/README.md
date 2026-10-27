# Infinite NPA depth at minimal semidefinite complexity

The full quantum region of the CHSH value and the marginal sum
\(M=\langle A_0\rangle+\langle B_0\rangle\) has an explicit affine lift with
PSD blocks of orders 5 and 4. Every fixed NPA level projects strictly
outside that same region. Its semidefinite extension degree is 2.

This is a sequel to the NPA Letter: the earlier result supplies infinite
depth, and the new proof supplies the finite representation. Its central
step is elementary. At fixed product of the two overlap parameters, the
determinant is affine in their sum; the arithmetic–geometric mean inequality
reduces the nonclassical support directions to equal overlaps. The remaining directions have local
bounds. One univariate quartic matrix, affine in all support coefficients,
then gives the complete finite lift.

- [Two-page note](note/main.pdf) — the central insight and proof.
- [Two-page PRL sequel draft](article/main.pdf) and [abstract](ABSTRACT.txt).
- [Complete proofs](proofs/main.pdf) — 13 pages, including all matrix entries,
  both membership inclusions and the entire NPA nontermination argument.
- [Portable source package](NPA_Finite_Lift_Overleaf.zip).
- [Exact data for all 32 projections](checks/ALL32_EXACT_PROJECTION_CONSTRAINTS.json).

The abstract is author-approved and reproduced verbatim. The revised Letter
includes the analytical proof of the finite lift in its two-page main text.
The source authorship follows the existing PRL manuscript.

## What the proof covers

| Result | Exact scope | Complete proof |
|---|---|---|
| Common-tilt quantum region | Every real support direction, all boundary points; one 5/4-block affine membership lift | F1–F2 |
| Minimal block order | Extension degree 2 by the planar theorem and nonpolyhedrality; the explicit displayed lift has orders 5 and 4 | F3 |
| Independently signed tilts | Exact fixed-input support SDP for every real pair of tilts, including axes, critical and commuting branches | F4 |
| Relabelled projections | All 32 marginal–CHSH planes; exact support after intersection with any NPA level | F5 |
| Infinite NPA depth | Divergent relative support error near the local threshold at every fixed level; extension to every fixed finite word list | S1–S5 |

The signed two-tilt construction is a support calculation with coefficients
depending on the chosen tilts. The full eight-coordinate quantum-set lift
and its universal C1/C2 admission questions remain a separate project;
see the repository [frontier](https://github.com/anubhav70747/q222-finite-sdp/blob/main/FRONTIER.md).

## Reproduce

Each of note/main.tex, article/main.tex and proofs/main.tex is self-contained.
Upload the ZIP to Overleaf and select the desired main document, or run:

    python build.py

This uses pdfLaTeX (with RevTeX 4.2 for the article) and writes intermediates
to build/. To check the exact algebra, install SymPy and run:

    python checks/verify_symmetric.py
    python checks/verify_membership_and_asymmetric.py
    python checks/build_orbit_constraints.py

The first program checks 17 polynomial identities. The second checks the
adjoint membership matrices, strict anchor and signed asymmetric formulas.
The third enumerates the physical relabellings and regenerates 64 PSD
blocks of total order 288, 512 auxiliaries and 96 affine equations.
These checks support the analytical proofs; they do not replace their
positivity, compactness or quantifier arguments.

[Verification record](VERIFICATION.json) · [Source provenance](PROVENANCE.md)
