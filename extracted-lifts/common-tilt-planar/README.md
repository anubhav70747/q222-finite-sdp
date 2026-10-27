# Exact common-tilt projection and 32 physical images

Preserved and published on 27 October 2026 from the previously delivered local `minimal-complexity-extension-2026-10-27` package. These results were outside Git; they were not lost or removed from the repository. This publication copies the original proof, data, checkers, reviews and both delivery archives without changing their bytes.

## Completed result and scope

The six-page [proof](proof/main.pdf) ([TeX](proof/main.tex)) gives an exact membership lift for

\[
W=\{(\langle A_0\rangle+\langle B_0\rangle,\mathrm{CHSH}):p\in\mathcal Q_{222}\}.
\]

Its evaluated matrices have PSD blocks of orders 5 and 4, sixteen auxiliary scalars and three affine equations. The homogeneous support construction covers every real common tilt and CHSH coefficient. Nonpolyhedrality together with Scheiderer's planar theorem gives semidefinite extension degree two; explicit order-two blocks are not evaluated here. The prior infinite-NPA-depth theorem is an imported input applying to support directions of this same region.

The [32-projection data](ALL32_EXACT_PROJECTION_CONSTRAINTS.json) transport this lift through all distinct marginal–CHSH planes. They add 64 PSD blocks, total order 288, 512 auxiliaries and 96 affine equations to a chosen finite NPA relaxation. The [coverage proof and comparison](COVERAGE_COMPARISON.txt) establish exact support on those planes and strict improvement over each fixed NPA level using the imported nontermination theorem. A combined solver instance and numerical benchmark are not supplied.

The signed asymmetric construction is exact for each fixed pair of tilts. Its nonlinear dependence on those tilts does not give one evaluated affine lift of the three-coordinate image W3. Neither these exact projections nor their intersection establish full eight-coordinate membership. Global H1/H2 admission remains open.

## Original materials and prior review

* [Results and remaining branches](RESULTS_AND_BRANCHES.txt).
* [Symmetric checker](verify_symmetric.py), [membership and asymmetric checker](verify_membership_and_asymmetric.py), and [orbit-data builder](build_orbit_constraints.py).
* [Prior symmetric proof review](historical/symmetric/PROOF.md), [its exact checker](historical/symmetric/check_exact.py), [recorded check results](historical/symmetric/CHECK-RESULTS.json), [branch map](historical/RESULT_BRANCH_MAP.md), and [signed extension record](historical/signed_extension.txt).
* Original [finite-lift package](archives/Finite_Lift_Extension_and_Branch_Map.zip) and [32-projection package](archives/Expanded_Coverage_32_Projections.zip).

The prior delivery state and reviews are historical evidence, with their original paths and scopes retained. The symmetric report does not independently establish the imported NPA theorem or every asymmetric extension. This import verifies preservation; it does not claim a new complete mathematical audit. Existing full-Q222 manuscripts, proof packages and global matrices are untouched.

## Reproduction and preservation

The original package instructions are in [package-readme.txt](package-readme.txt). The proof compiles with pdfLaTeX; the two exact mathematical checkers require Python and SymPy. The orbit builder uses the standard library. The existing delivery reported proof review and rendered-PDF inspection; those unchanged operations were not repeated just to archive the same bytes.

The [import receipt](IMPORT_VERIFICATION.json) records all 17 copied file hashes. Both ZIPs pass CRC checks, all 13 archive members match the corresponding current local files, and the proof TeX/PDF match the earlier delivery hashes. Run `python verify_import.py` to verify the published copies and archives. The local two-page manuscript candidate remains a separate delivery and is not silently substituted for this proof.
