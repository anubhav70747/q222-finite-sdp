# Q222 comparison and additive review package

Baseline: `ca4fb5e4acaadb234aa3f6a523b630e25ca99502`, 26 October 2026, commit “Verify the EL and NU endpoint certificates with complete state dependence”.

The current repository states the finite-SDP-lift existence theorem. This package accepts that source theorem, reconciles the R20 checkpoint with the later manuscript, and adds two deductions with exact computations. It does not certify every argument in the long companion afresh, and does not extract the global eight-coordinate affine pencil.

## What can be reviewed and integrated

1. **Angular coefficient cone:** `extracted-lifts/global-weighted/angular-order13/` gives a complete order-13 affine coefficient pencil, its integer arrays, proof including zero-support closure, and executable exact checks. On the Bell-source orthant it uses order 11; on either C1/C2 contact support, order nine and at most six squares. It represents the determinant coefficient cone. Substitution of the quartic Bell-source map is still nonlinear.
2. **Exposed Bell sources:** `independent-review/CONTACT_RAY_EXPOSURE.md` proves that every finite physical C1/C2 contact ray is exposed, including the allowed third omissions. A concrete mixture of one contact experiment and two deterministic neighbours exposes it in the full Bell-validity cone. Its source reconstruction, two exact examples, and limitations are stated.
3. **Checklist reconciliation:** 77 named rows distinguish historical R20 acceptance, the current manuscript's statements, and the narrower verification performed here.
4. **Inventories:** every repository file and every located project/title-search archive record has an explicit disposition. Downloaded or read content is not labelled mathematically verified.

## Corrections to earlier notes

The October 24 order-13 note did not exhibit the rotated weighted-triangle identity needed to preserve diagonal Gram subblocks. Merely reflecting its first weighted triangle generates only two of the four weighted triangles; directly rotating that identity creates forbidden within-subblock cross entries. The new proof supplies
\[
U^2V^2+U^2+2V-4UV
=(U-\sqrt V)^2+(UV-\sqrt V)^2+2U\sqrt V(\sqrt V-1)^2.
\]
All twelve circuits and their fixed diagonal blocks now pass exact coefficient checks.

The conditional atlas note also needs a precise distinction between variables. In the angular pencil \(U=u^2,V=v^2,t=\sqrt V=|v|\). The ramified order-49 completion's published `fixed_functions.json` instead has \(16\tau^4=|v|\). Its quarter-powers, normalization and testing multiplier define a different full-state feature library. Neither the angular order reduction nor the 36-half-power obstruction, by itself, establishes admission or exclusion of the sixteen-chart ramified atlas.

The two original October 24 notes and the audit are preserved under `archive/`; these corrections govern their use in the new package. The order-940 source-atlas / order-7520 behaviour proposal remains conditional. Its contact coverage and a full evaluated array export were not established or replayed here.

## R20 versus the current proof

R20 left C14.1a/b, C14.4 and weak-C2 R1–R4 residuals. The present C1 chapter has the all-\(h/\sigma\), all-Alice-angle completion and the five normalized small-\(L\) source cells. The C2 chapter now contains complete R1, R2 and joint R3/R4 sections. These are later claimed closures, so R20's OPEN labels must not be copied as the present theorem status.

The current review remains more qualified than a full fresh proof seal. It explicitly isolates the thin-R10 TL6-to-eight-root bounded-matrix implication, remaining source domains and their joins as verification targets. The final C2 paragraph still says “the only remaining task” despite its preceding claimed coverage. This is an editorial/interface issue to reconcile after checking the cited hypotheses. This package records the issue and does not remove it by assertion.

Existence and explicit extraction have different interfaces. Scheiderer's criterion is fed by the universal separated tensor-SOS premise; it does not require a separate uniform local square-count theorem. Extracting one evaluated affine pencil still requires fixed finite source/test features and both global inclusions. E1–E3 remain unresolved in the current repository.

## What was actually inspected

- Repository inventory: **577 files**, 70,738,279 bytes in the pinned tree. **546** nonbinary contents fetched: 315 text/source files and 231 successfully parsed JSON files. These are structural reads, not 546 proof certifications.
- **Three** large JSON files could not be fetched; their paths and Git blob identities are recorded. **28** PDF/ZIP outputs were inventoried separately. The Overleaf ZIP was retrieved and its twenty member records enumerated before the workspace connection was lost; its complete byte comparison against the source tree was not finished.
- Archive inventory: **689 distinct provider records**, the union of all 590 project-folder records and 296 Q222-title search results. **391** records were read to the end as text. **35** had no readable text extraction. Two large JSON records were read through line 81,000 only. **261** records remain metadata-only, principally binary/other artifacts. Distinct records with equal names or sizes were not assumed byte-identical.
- The handoff ZIP contents could not be exhaustively materialized or unpacked: the transfer helper encountered DNS failure, and the execution environment then became unavailable. This is why the historical archive comparison remains incomplete.

The two partial records are `C1D_LOCAL_AFFINE_PENCIL.coo.json` (247,070 lines) and `full_slack_level_11.json` (157,285 lines). The inventory records exact read boundaries. A full archive-byte and dependency audit remains outstanding. Reading historical “locked” status files is not a substitute for that audit.

## Executed verification

- The existing `independent-review/checks/el-nu-endpoints/run_exact_checks.py` completed with exit zero and “ALL EXACT CHECKS PASSED”. Its eight modules check source identities, sphere relations, full scalar graphs, both state normals, hard blocks, budgets and complementary identities on their printed scope.
- The new angular checker was evaluated as the shipped pure JavaScript function: all twelve circuit supports/certificates, all block-span and fixed-diagonal conditions, and all 22 visible/auxiliary coefficient basis identities passed in exact integer arithmetic.
- The new contact checker was evaluated as the shipped pure JavaScript function: both six-jet systems, selected minors 69120 and 184320, positive losing-event probabilities, three deterministic neighbours, and CHSH values \(103/37,46/17,107/39,8/3\) passed in exact rational arithmetic.
- No new numerical SDP claim is made. The remainder of the repository scripts were inspected as source/data, not all executed. Prior audit receipts are historical evidence, not re-labelled as this session's runs.
- No final LaTeX build/render was performed. pdfLaTeX was absent before the environment became unavailable. This package adds Markdown, data and code and does not change the existing manuscript or PDFs.

## Package use

This is an additive overlay against the pinned revision. It contains no changes to existing tracked scientific files. The historical notes are labelled as such; the new corrected note is the integration target.

The prepared Git object can be reviewed or cherry-picked into a local branch. No repository branch is updated by this preparation. Check the base revision, inspect the diff, and run the two standard-library JavaScript checkers in Node before integrating in another environment. The complete historical/full-parent audit is an explicitly separate unfinished obligation, not a requirement silently declared satisfied here.

See `CHECKLIST_CROSSWALK.md`, `PROVENANCE.json`, `REPLAY_REPORT.json`, `BLOCKERS.md`, and the two inventories.
