# Library audit and remaining steps

26 October 2026. Scientific baseline: `1ad5125dbfc9c56881c283724b04994e3d50ec00`.

The recovered library supports the qualitative finite-SDP existence argument, with the current source interfaces reconciled in V01–V06. That argument already implies the existence of some universal finite algebraic certificate space. The remaining construction problem is to identify an explicit space or finite collection of cones, prove that it admits every required source, and evaluate its affine matrices. A short complete replacement for the long contact proof remains a separate objective.

This audit read the substantive historical reviews and current verification notes, checked the scope of the actual construction bodies, compared the original “Recover manuscript” conversation with its recovered mathematical records, and inspected the two active derivations once. It is a reconciliation of the available library, not a fresh rederivation of all 215 components or a claim to have recovered every historical online file. The baseline and all original research files remain unchanged.

## Completed work retained

| Item | What is established | Evidence and precise limit |
|---|---|---|
| Physical model | Real two-qubit reduction, compactness and the nine-atom convex model | [Model](../latex/chapters/02_model.tex), [lead review](materials/historical/analysis/lead_completion_review.txt). Arbitrary-dimensional quantum realizations are reduced before the real-algebraic argument. |
| Canonical sources | Events, two closed-face types and native contact sources, including allowed omission supports | [Geometry](../latex/chapters/04_geometry.tex), [source reduction](../independent-review/CANONICAL_SOURCE_REDUCTION.md). |
| Full C1 local argument | The 45 source/test/state cells, full normals, opening faces and testing complements | [Complete C1 review](materials/historical/analysis/c1_completion/audit_report.txt), [case links](C1_CASES.md). These are substantive body reviews, with their stated shared dependencies. |
| Full C2 local argument | The 49 leaves, including R8 joint-D, thin R10, both normals and literal n=0 | [Complete C2 review](materials/historical/analysis/c2_completion/final_audit.txt), [case links](C2_CASES.md). |
| Current-version conjunction | V01–V06 reconcile the current native hypotheses and the existential global implication | [V01 verification](materials/historical/analysis/v01_conjunction_verification_2026_10_26.txt), [V06 verification](materials/historical/analysis/v06_global_conjunction_verification_2026_10_26.txt). Physical validity remains part of the source body. |
| Abstract finite library | The tensor-SOS criterion implies a finite algebraic certificate space and an existential finite pencil | [Eight-lemma global argument](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1328/V06_global_implication_original/Q222_V06_GLOBAL_IMPLICATION_CONCISE_PROOF_2026-10-26.md). Its basis, entries and size are not evaluated. |
| Closed-face cones | SG12/SG13 and their finite symmetry copies have exact pencils | [Geometry](../latex/chapters/04_geometry.tex), [assembly refinement](../independent-review/ASSEMBLY_REFINEMENT.md). No new proof of their membership in Z6 is required for a hybrid lift. |
| Full-set assembly | An exact canonical cone gives the full set through eight orientations; both inclusions are proved | [Assembly proof](../extracted-lifts/canonical-completion/full-polar/full_polar_assembly.tex). Exact canonical input and the stated duality conditions are required. |
| Compact strict source families | Every specified compact strict native family has one finite globally sound certificate cone | [Compact cover](../independent-review/COMPACT_NATIVE_CONTACT_COVER.md). Its grid size depends on the compact family. |
| Candidate fixed library | Soundness, closedness, a rational strictly feasible point and constant contact-map rank in the established strict native contact charts are established | [Candidate](../independent-review/EXACT_LIFT_CANDIDATE.md), [rank and assembly](../independent-review/ASSEMBLY_REFINEMENT.md). None alone proves universal PSD admission. |
| C2 open-end neighborhoods | Exact all-orders PSD Grams on 16 fixed features near both strict open ends; uniform on compact interior face sets | [Original theorem](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1407/FINITE_C_THEOREM_AND_PROOF.md), [analytic verification](materials/historical/analysis/finite_c_open_end_analytic_verification_2026_10_26.txt), [exact checks](materials/historical/analysis/finite_c_open_end_code_verification_2026_10_26.txt). This is stronger than an endpoint identity, but excludes arbitrary finite c and coupled source-wall limits. |
| Restricted explicit lifts | Exact equal-angle cone, unequal-angle curve and other restricted source constructions | [Equal-angle cone](../extracted-lifts/c2-equal-angle-cone/), [unequal-angle curve](../extracted-lifts/c2-unequal-angle-curve/), [order-635 construction](../extracted-lifts/canonical-completion/contact-completion/). Their stated source domains remain essential. |
| Library delivery | Original source/presentation bodies, dependency bindings, packages and verification evidence are available | [Inventory](INVENTORY.md), [source map](SOURCE_MAP.json), [verification receipt](VERIFICATION.json). Delivery checks establish identity and availability; the mathematical reviews supply the scientific conclusions. |

The original recovery conversation's unequal-angle result agrees with the recovered [geometry verification](materials/historical/analysis/gpt6_unequal_angle_geometry_verification_2026_10_26.txt) and [order-six curve review](materials/historical/analysis/unequal_curve_order6_review_2026_10_26.txt). Its three-parameter hull remains open. No overlooked full-hull or full-Q222 construction was found there.

## Remaining checklist

These are completion conditions, not a requirement to execute every possible proof route. R01–R05 can be discharged by one sufficiently strong admission theorem. They refine L08–L10 in the [material checklist](CHECKLIST.md).

- [ ] **R01 — Specify the explicit finite construction.** Fix the testing functions, source-to-coefficient maps and allowed additional globally sound cones. Their finite dimensions must be independent of the source. The sixth-root library is a candidate; a proved finite hybrid alternative is sufficient. **Maps to L08/L09.**

- [ ] **R02 — Prove complete strict C1 admission.** For every physically valid strict C1 source, supply an exact PSD Gram or equivalent certificate in the chosen construction. Preserve the complete four-state parent, independent testing angles and all coefficient equations. A rational affine preimage or positive diagonal masses alone do not close this item. **Maps to L09, with the strict-family part made explicit.**

- [ ] **R03 — Prove complete strict C2 admission.** Cover arbitrary finite native c and the source parameters outside the verified open-end neighborhoods. Retain the known obstruction to universal use of the smaller half-power library. **Maps to L08.**

- [ ] **R04 — Include all singular source limits in the same finite construction.** Include permitted zero weights, additional omissions, compatible-contact limits, vanishing state coordinates and simultaneous openings at independent relative rates. Full strict-source admission plus a proved density statement and closedness may settle these limits together. Otherwise supply genuine neighborhood certificates, not only limiting Grams. **Maps to L09.**

- [ ] **R05 — Establish complete testing and canonical-source coverage.** Every cone used globally must certify all independent measurements and states, including projective endpoints. Join regional testing certificates for the same source before global use. Then apply the existing canonical classification and include event and closed-face orbits. **Maps to L09.**

- [ ] **R06 — Instantiate the proved duality and assembly.** Establish the required closedness, strict feasibility or explicit facial reduction for the actual completed blocks. Insert the exact canonical cone into the existing eight-orientation or hybrid assembly. **Maps to L10.**

- [ ] **R07 — Evaluate the affine pencil.** Deliver every matrix in A0 + sum(q_i A_i) + sum(y_j B_j), each exact coefficient, all auxiliary equations, block orders and variable counts. Supply deterministic generation instructions. If algebraic constants occur, identify them exactly. **Maps to L10.**

- [ ] **R08 — Verify both inclusions for the delivered matrices.** Instantiate quantum-to-feasible and feasible-to-quantum arguments for the actual evaluated system. Check the decisive new identities, PSD arguments and universal coverage independently. Reuse unchanged reviewed dependencies. **Maps to L10.**

- [ ] **R09 — Obtain a genuinely shorter complete proof.** Replace the long native-contact argument by proved general mechanisms, preserving every singular case and dependency. [PEN_AND_PAPER](../independent-review/PEN_AND_PAPER.md) still imports I5 from the long companion; the seven-page exposition does not itself replace that proof.

- [ ] **R10 — Integrate the final manuscript.** Incorporate the completed admission theorem, exact domains, constructive formulas and final lift into a connected article and pedagogical companion. Preserve the distinction between finite NPA impossibility and the more general finite lift. L12 is already published as evidence; it also needs integration into the eventual final constructive exposition.

- [ ] **R11 — Rebuild and reproduce the changed package.** Compile and inspect the revised PDFs, verify references and source inclusion, reproduce the new exact matrix construction in a clean environment, and synchronize matrix data, article, companion and Overleaf ZIP. Current delivery checks remain completed; new edits require new checks.

- [ ] **R12 — Publish verified completion.** Push the revised sources, exact data and evidence to main; verify remote identity and links; close L08–L10 and R01–R11 against their actual deliverables. Pause the recurring proof checks only when the exact lift and short complete proof are verified.

## Immediate dependency and assignment map

| Work | Present assignment or next dependency |
|---|---|
| R02 strict C1 | Existing Extra High C1 assignment; preserve its current run. |
| R03 strict C2 | Existing Extra High C2 assignment; preserve its current run. |
| R01/R04/R05 | Consolidate completed outputs before choosing direct global admission, finite boundary-neighborhood coverage or a proved alternative. |
| R06–R08 | Require an exact canonical construction from R01–R05. |
| R09 | Can develop alongside the new mechanism, but closure requires a complete replacement for the imported long argument. |
| R10–R12 | Integrate and release only the proved final construction; keep existing original evidence available. |

At **20:53 IST** on 26 October both assignments supplied complete normal answers and original packages. Independent analytical reviews and exact data replay validate the C1 eight-gauge coefficient reduction (L13) and the C2 compact-dual reduction (L14). Neither proves the remaining universal PSD inequality. The next work therefore starts from the explicit reduced problems; R02 and R03 remain open. The earlier interrupted C2 state has ended.

## Route-dependent work and work already superseded

- **Boundary-neighborhood route:** construct finitely many all-testing certificates on true neighborhoods of the singular source–contact boundary, prove a fixed finite span in every neighborhood, and show the remainder lies in finitely many compact strict native charts. The existing compact-interior theorem then handles that remainder. This is a proposed route, not a theorem established by this audit.
- **Direct dual route:** prove nonnegative evaluation of every normalized dual PSD functional on every native parent. A zero margin can still give valid admission. A global strict-Gram or flat-extension theorem is stronger than necessary. This route remains a proposed algebraic target.
- **Continuation route:** if used, prove the actual positivity/nondegeneracy continuation step. Constant rank, openness, connectedness and closedness of the feasible cone do not by themselves establish it.
- **R1 regional route:** the [order-283/293 construction](../extracted-lifts/canonical-completion/uniform-contact/) leaves the testing strip |u| > 9/256, |v| < 1/16. Complete it only if that regional construction is retained as a global component. Another admission theorem may bypass it.
- **Three-parameter unequal-angle hull:** remains a possible constructive subroute. Its exact compact curve is already completed; the larger hull need not be separately solved if another construction covers it.
- **Known cones at extreme rays:** an actual extreme native ray cannot be split into nonproportional globally valid event/face/contact summands. It must be admitted directly by a retained cone.
- **Closed faces:** retain SG12/SG13 directly. Their inclusion in the sixth-root cone is optional.
- **NPA and size optimization:** one fixed NPA level can be retained for its proved coverage. Finding the smallest lift or extending numerical sweeps is not required for exactness.
- **Earlier missing intermediates:** the new self-contained L12 files supersede the lost-runtime issue for that theorem. Historical missing pickles are provenance limits, not a current L12 dependency.
- **Historical releases and additional replays:** see [recovery limits](RECOVERY_LIMITS.md). The separate 549-page release, damaged nested archive and interrupted external multipart extraction remain unrecovered; no missing dependency of the indexed current proof was identified from them. Additional historical checker replays are separate reproducibility work, not a reason to reopen V01–V06.

## Audit finding

No substantive contradiction was found between the completed body reviews, the V01–V06 conjunction and the scoped construction results. No overlooked explicit universal fixed-library theorem was found in the inspected library material. The remaining scientific work is explicit universal admission and extraction (L08–L10), followed by a shorter complete proof and integrated delivery. Earlier narrow “not checked in this task” statements do not erase the complete reviews. Current status must follow the actual theorem domains, not package names, file counts or interim chat claims.
