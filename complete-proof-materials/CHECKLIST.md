# Complete material checklist — 26 October 2026

This checklist links to the proofs, coefficients, matrix data, exact checkers and recorded results themselves. A checked box below means that the identified material and its stated result are on file; it is not a fresh certification of the entire theorem.

The 22 October full analytical review supports qualitative finite-SDP existence. It covers 215 indexed components, 45 C1 cells and 49 C2 leaves. The current-version interface checks and final source/global conjunction are now reconciled in V01-V06 using those complete body reviews. Constructing one evaluated global pencil remains a separate unfinished task.

- [All 215 components](COMPONENTS.md)
- [Every C1 cell](C1_CASES.md) and [every C2 leaf](C2_CASES.md)
- [All published material](INVENTORY.md), [machine-readable checklist](CHECKLIST.json), and [file hashes](MANIFEST.json)
- [What the earlier work already completed](RECONCILIATION.md)
- [Recovery limits](RECOVERY_LIMITS.md) and [the next verification pass](VERIFICATION_PLAN.md)
- [Library audit and actionable remaining steps](REMAINING_STEPS.md)

## Foundations, full local proof and delivery

- [x] **P01 — Real two-qubit reduction and compact nine-atom physical model.** prior complete analytical review.
  Evidence: [02_model.tex](../latex/chapters/02_model.tex); [lead_completion_review.txt](../complete-proof-materials/materials/historical/analysis/lead_completion_review.txt).
  Completion condition/scope: Preserve the real-state reduction, compactness and convex mixture argument.

- [x] **P02 — Universal tensor-SOS criterion and algebraic coefficients.** prior review; targeted criterion recheck available.
  Evidence: [03_tensor.tex](../latex/chapters/03_tensor.tex); [global_tensor_criterion_audit_2026_10_26.txt](../complete-proof-materials/materials/historical/analysis/global_tensor_criterion_audit_2026_10_26.txt); [global_criterion_adversarial_check.txt](../complete-proof-materials/materials/historical/analysis/c2_completion/global_criterion_adversarial_check.txt).
  Completion condition/scope: Retain the field quantifiers and both inclusions; no fixed square count is required for qualitative existence.

- [x] **P03 — Canonical source reduction, all omission supports, stationary contact validity.** prior review; targeted geometric recheck available.
  Evidence: [04_geometry.tex](../latex/chapters/04_geometry.tex); [canonical_coverage_short_audit_2026_10_26.txt](../complete-proof-materials/materials/historical/analysis/canonical_coverage_short_audit_2026_10_26.txt); [CANONICAL_SOURCE_REDUCTION.md](../independent-review/CANONICAL_SOURCE_REDUCTION.md).
  Completion condition/scope: Include both closed-face orbits, event rays, finite contacts and permitted third omissions.

- [x] **P04 — Native source/contact entrance.** prior complete analytical review.
  Evidence: [05_entry.tex](../latex/chapters/05_entry.tex); [en_native_entry_audit_2026_10_26.txt](../complete-proof-materials/materials/historical/analysis/en_native_entry_audit_2026_10_26.txt).
  Completion condition/scope: Match actual source parameters to all contact equations, independently of testing variables.

- [x] **P05 — Noncommuting residue theorem A6 and strict-Gram mechanism.** prior complete analytical review.
  Evidence: [noncommuting_refined.tex](../latex/chapters/noncommuting_refined.tex); [lead_completion_review.txt](../complete-proof-materials/materials/historical/analysis/lead_completion_review.txt).
  Completion condition/scope: Use the full event matching/induction and state/projective joins, not only positive residues.

- [x] **P06 — All 45 C1 source/test/state cells.** prior complete analytical review.
  Evidence: [All C1 cases](C1_CASES.md); [audit_report.txt](../complete-proof-materials/materials/historical/analysis/c1_completion/audit_report.txt); [refined_delta_audit.txt](../complete-proof-materials/materials/historical/analysis/c1_completion/refined_delta_audit.txt).
  Completion condition/scope: Preserve every listed hypothesis, complete state form, opening and complementary test.

- [x] **P07 — All 49 C2 source/test leaves.** prior complete analytical review.
  Evidence: [All C2 cases](C2_CASES.md); [final_audit.txt](../complete-proof-materials/materials/historical/analysis/c2_completion/final_audit.txt); [refined_delta_audit.txt](../complete-proof-materials/materials/historical/analysis/c2_completion/refined_delta_audit.txt).
  Completion condition/scope: Preserve strong-angle branches, weak-source R1–R10, both normals, literal n=0 and projective joins.

- [x] **P08 — Global source cover and finite assembly.** prior complete analytical review.
  Evidence: [07_cover.tex](../latex/chapters/07_cover.tex); [final_audit.txt](../complete-proof-materials/materials/historical/analysis/final_audit.txt); [coverage_union.json](../complete-proof-materials/materials/historical/analysis/coverage_union.json).
  Completion condition/scope: Keep local implications under their native hypotheses and the universal criterion. Coverage counts do not prove those implications.

- [x] **P09 — Refined BA, CT, DS/CJ and BN/UN changes.** prior delta reviews complete.
  Evidence: [refined_delta_audit.txt](../complete-proof-materials/materials/historical/analysis/c1_completion/refined_delta_audit.txt); [refined_delta_audit.txt](../complete-proof-materials/materials/historical/analysis/c2_completion/refined_delta_audit.txt).
  Completion condition/scope: Use generic BA identities for R8/R10 with their own source bounds; preserve all unit testing ratios in CT.

- [x] **P10 — Original 215 source bodies and 215 presentation copies.** materially delivered; source identity checked.
  Evidence: [Component-by-component proof index](COMPONENTS.md); [COMPONENT_INDEX.json](../complete-proof-materials/materials/historical/recovery/full_r21/records/COMPONENT_INDEX.json); [PROOF_DEPENDENCIES.json](../complete-proof-materials/materials/historical/recovery/full_r21/records/PROOF_DEPENDENCIES.json).
  Completion condition/scope: Check original source and transformed presentation hashes separately; preserve all original path bindings in the archive.

- [x] **P11 — Original full archive and 798-page recovered volume.** materially delivered.
  Evidence: [Q222_R21_Rederived_Proof_Package_2026-10-22.zip](../complete-proof-materials/materials/historical/recovery/full_r21/archives/Q222_R21_Rederived_Proof_Package_2026-10-22.zip); [Q222_Analytical_Existence_Proof_RECOVERED_2026-10-22.pdf](../complete-proof-materials/materials/historical/recovery/full_r21/Q222_Analytical_Existence_Proof_RECOVERED_2026-10-22.pdf); [Q222_Complete_Package.zip](../complete-proof-materials/materials/historical/delivery/Q222_Complete_Package.zip).
  Completion condition/scope: Preserve byte identity and distinguish these historical versions from current manuscript sources.

- [x] **P12 — Current paper, companion, sources and matrix data.** already published; indexed here.
  Evidence: [Q222_Paper.pdf](../pdf/Q222_Paper.pdf); [Q222_Complete_Proofs.pdf](../pdf/Q222_Complete_Proofs.pdf); [Q222_Overleaf.zip](../Q222_Overleaf.zip); [All evaluated representations](../extracted-lifts).
  Completion condition/scope: Keep sources, matrices, exact coefficient systems, code and results together.


## Explicit lift construction

- [x] **L01 — Exact closed-face pencils SG12/SG13 and symmetry copies.** proved in supplied source; prior review.
  Evidence: [04_geometry.tex](../latex/chapters/04_geometry.tex); [ASSEMBLY_REFINEMENT.md](../independent-review/ASSEMBLY_REFINEMENT.md).
  Completion condition/scope: Use the exact globally sound face cones directly; no sixth-root inclusion for them is needed.

- [x] **L02 — Canonical-to-full eight-copy polar assembly.** proved conditional construction.
  Evidence: [full_polar_assembly.tex](../extracted-lifts/canonical-completion/full-polar/full_polar_assembly.tex).
  Completion condition/scope: Insert a single exact canonical pencil; both inclusions are already established for that input.

- [x] **L03 — Compact strict native source families.** scoped complete construction.
  Evidence: [COMPACT_NATIVE_CONTACT_COVER.md](../independent-review/COMPACT_NATIVE_CONTACT_COVER.md).
  Completion condition/scope: The finite size depends on the compact subfamily; this does not supply uniform singular-limit admission.

- [x] **L04 — Explicit regional R1 and C2 certificates.** scoped constructions and matrices delivered.
  Evidence: [R1 data and source](../extracted-lifts/canonical-completion/uniform-contact); [C2 regional completions](../extracted-lifts/canonical-completion/fixed-c2-gates); [r1_full_proof_audit_2026_10_26.txt](../complete-proof-materials/materials/historical/analysis/r1_full_proof_audit_2026_10_26.txt).
  Completion condition/scope: Retain the stated source and testing domains when joining regional cones.

- [x] **L05 — Equal-angle and unequal-angle C2 exact subfamilies.** scoped complete representations.
  Evidence: [c2-equal-angle-cone](../extracted-lifts/c2-equal-angle-cone); [c2-unequal-angle-curve](../extracted-lifts/c2-unequal-angle-curve).
  Completion condition/scope: Retain their exact parameter domains and source/behavior inclusions.

- [x] **L06 — Sixth-root library specification, soundness and closedness.** candidate properties established in supplied analysis.
  Evidence: [EXACT_LIFT_CANDIDATE.md](../independent-review/EXACT_LIFT_CANDIDATE.md); [ASSEMBLY_REFINEMENT.md](../independent-review/ASSEMBLY_REFINEMENT.md).
  Completion condition/scope: These properties make it a valid candidate cone; they do not prove every native source belongs to it.

- [x] **L07 — Closed-form C2 c=0 and c=infinity limit Grams.** exact endpoint result; local symbolic replay completed.
  Evidence: [C2_END_GRAM_LEMMA.md](../independent-review/C2_END_GRAM_LEMMA.md).
  Completion condition/scope: Keep the distinct state rescaling at the two ends. Endpoint positivity alone does not give finite-c positivity.

- [ ] **L08 — Finite-c complete-state certificate in one fixed library.** OPEN for arbitrary finite native c; exact neighborhoods of both strict open C2 ends verified in L12.
  Evidence: [EXACT_LIFT_CANDIDATE.md](../independent-review/EXACT_LIFT_CANDIDATE.md); [Finite-c recovered record](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1207/C2_finite_c_recovery/C2_Finite_c_Recovery_Record.md); [Missing runtime dependencies](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1207/C2_finite_c_recovery/MISSING_DEPENDENCIES.txt); [Original open-end theorem and proof](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1407/FINITE_C_THEOREM_AND_PROOF.md); [Independent analytic PSD verification](materials/historical/analysis/finite_c_open_end_analytic_verification_2026_10_26.txt); [Exact source and coefficient verification](materials/historical/analysis/finite_c_open_end_code_verification_2026_10_26.txt).
  Completion condition/scope: Give a Gram or equivalent exact certificate for the complete parent at all finite native c, with actual coefficient matching and positivity.

- [ ] **L09 — Uniform admission through every singular source limit.** OPEN construction.
  Evidence: [ASSEMBLY_REFINEMENT.md](../independent-review/ASSEMBLY_REFINEMENT.md); [COMPACT_NATIVE_CONTACT_COVER.md](../independent-review/COMPACT_NATIVE_CONTACT_COVER.md).
  Completion condition/scope: Prove every remaining C1/C2 source lies in one fixed finite sum of globally sound cones; include zero weights, independent rates and testing endpoints.

- [ ] **L10 — Final evaluated global affine pencil and both inclusions.** OPEN construction.
  Evidence: [full_polar_assembly.tex](../extracted-lifts/canonical-completion/full-polar/full_polar_assembly.tex).
  Completion condition/scope: After L08/L09 or a proved alternative, publish every matrix, auxiliary variable, coefficient equation and an independently replayable two-inclusion proof.

- [x] **L11 — NPA and simpler feature-library limitations.** proved/scoped obstruction evidence delivered.
  Evidence: [SMOOTH_FEATURE_OBSTRUCTION.md](../independent-review/SMOOTH_FEATURE_OBSTRUCTION.md); [NATIVE_MULTIPLIER_OBSTRUCTION.md](../independent-review/NATIVE_MULTIPLIER_OBSTRUCTION.md); [C2_HALF_POWER_OBSTRUCTION.md](../independent-review/C2_HALF_POWER_OBSTRUCTION.md).
  Completion condition/scope: Distinguish failure of one fixed NPA/feature cone from impossibility of an arbitrary finite SDP lift.


- [x] **L12 — Exact C2 finite-c admission near both strict open c-ends.** VERIFIED on the stated strict real face domain, uniformly on compact interior subsets.
  Evidence: [Original open-end theorem and proof](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1407/FINITE_C_THEOREM_AND_PROOF.md); [Independent analytic PSD verification](materials/historical/analysis/finite_c_open_end_analytic_verification_2026_10_26.txt); [Exact source and coefficient verification](materials/historical/analysis/finite_c_open_end_code_verification_2026_10_26.txt); [Verified delivery ZIP (assembled from original files)](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1407/C2_open_end_verified_delivery.zip).
  Completion condition/scope: For every fixed a,t,eta,rho>0 with rho<t, 0<t-rho<1, eta<1 and a^2*t+eta>1, the complete native parent has an exact PSD Gram on 16 fixed features for sufficiently small positive epsilon=c^-2. The gauge, full 2x2 tail Schur complement, constant rational all-orders correction and inversion are verified. Coupled approaches to source walls and arbitrary finite c are outside this result.

- [x] **L13 — Strict C1 fixed-span coefficient and eight-gauge reduction.** VERIFIED affine reduction and exact PSD equivalence; universal feasibility OPEN.
  Evidence: [Original C1 eight-gauge proof](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1522/task5/package/Q222_C1_exact/PROOF.md); [Original complete C1 package](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1522/task5/original-package.zip); [Independent C1 analytical review](materials/historical/analysis/c1_strict_reduction_analytic_verification_2026_10_26.txt); [Exact data replay and limitations](materials/historical/analysis/strict_interior_reductions_code_verification_2026_10_26.txt).
  Completion condition/scope: On the genuine strict real C1 domain A1, the fixed 16-feature invariant map has shape 42x44, rank36 and eight free gauges. Exact admission to this subcone is equivalent to H7(theta;s)>=0, H5(theta;s)>=0 and s7>=0. The 42 rows are the support generated by these features, not all 46 parity-compatible monomials. Forced positive diagonals and partial minors do not solve full PSD feasibility. Zero-weight walls are excluded.

- [x] **L14 — Strict C2 fixed dual inequality and exact coefficient data.** VERIFIED affine/dual reduction; nonnegative dual inequality OPEN.
  Evidence: [Original C2 dual reduction](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1522/task3/package/interior/STRICT_C2_INTERIOR_THEOREM_AND_OPEN_GATE.md); [Original complete C2 package](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1522/task3/original-package.zip); [Independent C2 analytical review](materials/historical/analysis/c2_strict_reduction_analytic_verification_2026_10_26.txt); [Exact data replay and limitations](materials/historical/analysis/strict_interior_reductions_code_verification_2026_10_26.txt).
  Completion condition/scope: On strict real six-event Omega with n>0, the contact-reduced 154-feature map has 412 rows,3021 columns and rank406. Six exact jets give affine reachability. Closedness and the compact normalized dual give PSD admission iff the displayed cleared 24-term evaluation is nonnegative for every dual point. Data and rank were exactly replayed. Physical curvature does not prove the inequality; source walls and uniform singular limits remain outside the result.

## Reconciliation and verification

- [x] **V01 — Reconcile the older complete review with later narrow C1/C2 reviews.** RESOLVED current-version native source-coverage and acyclic dependency conjunction.
  Evidence: [Source-coverage conjunction](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1328/V01_THEOREM_AND_PROOF.md); [Exact dependency crosswalk](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1328/FINITE_DEPENDENCY_CROSSWALK.csv); [Independent source and hypothesis verification](materials/historical/analysis/v01_conjunction_verification_2026_10_26.txt).
  Completion condition/scope: The current source partitions, whole-parent testing joins and V02-V04 interfaces supply exactly V06 NATIVE for arbitrary independently labelled real closed fields, actual physical C1/C2 sources, all independent projective tests and states, and permitted zero faces. This retains the earlier complete body reviews; it is not a new derivation of unchanged leaves.

- [x] **V02 — Thin-R10 formal-base and conductor interfaces in the latest source.** RESOLVED for the frozen thin-R10 version/interface task.
  Evidence: [Five-interface analytical reconciliation](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1207/thin_R10_reconciliation/THIN_R10_INTERFACE_RECONCILIATION.md); [Independent native-source verification](materials/historical/analysis/thin_r10_interface_verification_2026_10_26.txt); [r10_source_and_formal_base_audit_2026_10_26.txt](../complete-proof-materials/materials/historical/analysis/r10_source_and_formal_base_audit_2026_10_26.txt); [ca_se_conductor_structural_review_2026_10_26.txt](../complete-proof-materials/materials/historical/analysis/ca_se_conductor_structural_review_2026_10_26.txt); [r10_all_tests_join_audit_2026_10_26.txt](../complete-proof-materials/materials/historical/analysis/r10_all_tests_join_audit_2026_10_26.txt); [final_audit.txt](../complete-proof-materials/materials/historical/analysis/c2_completion/final_audit.txt).
  Completion condition/scope: The quarter-parent stress, both normals, native CA/SE frames, normalized two-stage SE division, negative-beta reserve and all-test imports compose at 332483309586e0fec6e16fc1dd855915b648dca1. This closes this interface check, not all leaves or an evaluated global lift.

- [x] **V03 — IC/LII current-source and full-parent proof reconciliation.** RESOLVED for the specified native IC and LII sectors.
  Evidence: [IC/LII version-and-proof audit](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1207/IC_LII_3324833_VERSION_AND_PROOF_AUDIT.md); [Independent native-source verification](materials/historical/analysis/ic_lii_interface_verification_2026_10_26.txt); [Recovery limits](RECOVERY_LIMITS.md).
  Completion condition/scope: The independent comparison and targeted source verification found no correction needed. Earlier inaccessible chat history remains unrecovered; it supplies no substantiated defect. Other C2 sectors and global construction remain separate.

- [x] **V04 — Latest BN/UN, R2, R3/R4, R6/R7 and UL readbacks.** RESOLVED for the specified current-version normalization, domain and import interfaces.
  Evidence: [Full frozen-version crosswalk](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1248/Q222_V04_FROZEN_INTERFACE_CROSSWALK_2026-10-26.md); [Independent native-source verification](materials/historical/analysis/v04_interface_verification_2026_10_26.txt); [Original evidence and checker package](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1248/V04_crosswalk_original_evidence.zip).
  Completion condition/scope: Corrected BN/UN factors, restricted reflected k=1 subdomains, R3/R4/R6/R7 replacement routing and UL whole-parent restoration agree with current native sources. Existing source bounds and previously reviewed imports remain explicit. Final V01/V06 conjunction is recorded below; L08-L10 construction remains separate.

- [x] **V05 — Identity, link, archive and delivery coverage checks.** performed for this publication.
  Evidence: [Integrity verifier](verify_materials.py); [Verification receipt](VERIFICATION.json).
  Completion condition/scope: Verify bytes, all component hashes, source/dependency availability and all generated checklist links; do not call this a new mathematical proof.

- [x] **V06 — Final independent mathematical verification of the integrated package.** RESOLVED integrated qualitative-existence implication on the reviewed source.
  Evidence: [Eight-lemma global proof](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1328/V06_global_implication_original/Q222_V06_GLOBAL_IMPLICATION_CONCISE_PROOF_2026-10-26.md); [Independent global verification and notation clarification](materials/historical/analysis/v06_global_conjunction_verification_2026_10_26.txt); [Original V06 proof package](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1328/V06_global_implication_original.zip).
  Completion condition/scope: V01 supplies the exact NATIVE input. Independent verification checks the transferred physical model and valid source body, separated tensor Gram, affine constant and zero weights, primary Corollary 3.19, and both existential lift inclusions. The review shorthand for the normalized source body must retain physical validity as the current manuscript does. L08-L10 remain open; no evaluated pencil or size bound follows.
