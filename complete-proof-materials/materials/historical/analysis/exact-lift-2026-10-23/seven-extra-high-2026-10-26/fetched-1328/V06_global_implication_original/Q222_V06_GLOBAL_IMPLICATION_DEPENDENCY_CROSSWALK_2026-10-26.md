# V06 dependency, quantifier and version crosswalk

**Audit outcome:** valid conditional implication. **Only labelled unverified premise for V06:** `(NATIVE)` as precisely stated in the accompanying proof. Its comprehensive source/test coverage is reserved to the separate final native review. The already printed geometric, closed-face, and noncommuting lemmas were inspected at the global interface; their hypotheses were not replaced by coverage counts.

## Frozen source evidence

- Repository main `8108102c092894c62768b7c8b0008889737bf1e7`: `latex/chapters/01_theorem.tex`, `02_model.tex`, `03_tensor.tex`, `04_geometry.tex`, `05_entry.tex`, `06_mechanisms.tex`, `07_cover.tex`, `noncommuting_refined.tex` **all match text at** `332483309586e0fec6e16fc1dd855915b648dca1`. The user-provided `04_QUALITATIVE_GLOBAL_IMPLICATION_SOURCES_1252.txt` reproduces the first, second, third, sixth, seventh and noncommuting chapters and previous scoped audits.
- External primary criterion: C. Scheiderer, *Second-order cone representation for convex sets in the plane*, SIAM J. Appl. Algebra Geom. 5 (2021), 114–139; `https://arxiv.org/pdf/2004.04196`; Cor. 3.19, PDF p. 15; Thm. 3.10, PDF p. 12; Thm. 2.10 and the finite-algebraic-space construction. The title refers to the main planar application, but **Cor. 3.19 itself is dimension independent**.
- Earlier global audits in the uploaded bundle: `analysis/global_tensor_criterion_audit_2026_10_26.txt`, `analysis/c2_completion/global_criterion_adversarial_check.txt`, `analysis/lead_completion_review.txt`, `analysis/final_audit.txt`. Their verdicts are evidence, not independent assumptions. V02 thin-R10, V03 IC/LII, V04 changed native routing are version/interface reviews and do not replace `(NATIVE)`.

## Inputs and outputs, without hidden promotion

| ID | Genuine input | Output used globally | Exact reference | Limit |
|---|---|---|---|---|
| MODEL | Two binary POVMs, arbitrary finite local dimensions | Polynomial real two-qubit atom map; at most nine convex atoms, compact/convex/semialgebraic real Q | RFA1–RFA5; `02_model.tex` | No Hilbert-space claim over non-Archimedean K; extend the polynomial formula |
| POLAR | Nine-atom first-order formula | \((P_Q)_K\) equals all K-affine functions nonnegative on \(Q_K\) | RFA6; TF11 | Transfer applies to *all* real-closed K over R |
| SOURCE | f valid; NS positivity and PR vertices; canonical SG1; stationary and support results | Nonnegative events or ≤9 extreme canonical atoms: events, closed faces, native C1/C2 including third omissions | SG1, SG5–SG10, SG17–SG22, EN1–EN16 | Extreme decomposition transferred via finite first-order midpoint predicate; no K-compactness |
| EVENT | Source event coefficient \(a_e\ge0\) | Pure-source times pure-testing rank-one square | SG1–SG2, RFA12 | Literal zero valid |
| CLOSED | Source in two noncoercive exclusive-omission faces | SG12/SG13 exact LMI; attained PSD test dual; separated trace SOS | SG11–SG14 | Strict feasible full pencil needed for attained dual and closed image |
| NATIVE | **Labelled assumption:** actual SG10 contact, physical and all permitted boundaries | Full homogeneous parent SOS at **all** independent projective tests and arbitrary states over \(K_c\otimes K_t\) | LC2 / `07_cover.tex` canonical theorem | **Only independent native-coverage check not performed in V06** |
| GRAM | One bounded tensor coefficient Gram with PD real residue, after legal retained packets and source/test normalization | Exact finite tensor SOS on arbitrary roots | TF1–TF8; CC1–CC6; UA1–UA3 | A numerical/pointwise PSD test does not suffice; no mixed inverse |
| NONCOMMUTING | Physical source, all original event weights nonnegative; **both** testing measurement pairs noncommuting **in residue** | Full-state SOS for these independent test atoms | NC1–NC26; `noncommuting_refined.tex` | Not tests with nonzero *infinitesimal* incompatibility; cannot replace NATIVE |
| COMMUTING | One measurement pair literally commuting | Testing-field PSD CHSH core plus source-weighted event squares | RFA10, literal commuting paragraph | Not a near-commuting test |
| SCALAR-SIGN | Affine positive angles, tested WZ≤0, original canonical signs | Pure-testing nonnegative CHSH core plus events | SG15 | Source physicality is still required for other state sectors; not a general scalar determinant certificate |
| PARENT | Full homogeneous RFA10/LC2 SOS, nonzero test angle pair and state | Independent atom slack by inverse **pure-testing** square root | RFA9–RFA13, TF8 | Never invert a mixed tensor or a zero homogeneous coordinate |
| ASSEMBLE | Finite nonnegative source atom weights; ≤9 testing weights summing to 1 | Universal f-first tensor SOS; zeros harmless | M2 / TF10 | Source weights need not sum to one; test weights must |
| CRITERION | Closed convex semialgebraic Q; SOS for every K, f, q | Some finite spectrahedral shadow; no fixed rank requirement | Scheiderer Cor. 3.19 | Thm. 3.10 concerns *specified* bound d |
| PENCIL | Finite X,U from Cor. 3.19 proof | Fixed affine LMI in finite moments; both inclusions by eval/separation | V06 proof Lemma 8, Thm. 2.10 | Size and entries not evaluated; no finite NPA level |

## Three fragile implications checked explicitly

1. **No confusion about extremes and boundaries.** An arbitrary physical canonical source is decomposed in the compact **real** normalized cone, with \(d=0\) giving event atoms and closed faces/noncoercive cases retained. The fixed-nine/extremality assertion is first-order, so decompositions transfer to any \(K_c\). Strict positivity in the SG10 support means positivity **in the coefficient field**, not necessarily positive real residue. Literal \(D=0,\beta=0,n=0\) are included in `(NATIVE)`. An infinitesimal positive coefficient must not be reset to zero.
2. **No confusion about positive matrices.** Closed-face dual factorizations take place separately in \(K_c\) and \(K_t\). Local strict-Gram SOS requires a complete coefficient identity and bounded coefficients with **positive-definite real residue**. Merely \(P_c(u,v)\succeq0\) over one common field, or a determinant circuit SOS, would not certify \(K_c\otimes K_t\).
3. **No confusion about schematic extraction.** Cor. 3.19 supplies a finite algebraic space \(U\) only **existentially**. The fixed moment pencil built from \(\{1,\phi_i,u_iu_j\}\) has both inclusions, but no evaluated matrices follow from the native root lists. A fixed grid Z6, a common tensor-rank bound, or a finite NPA level is an additional theorem, not this criterion.

## Smallest readable manuscript architecture

**Main proof, approximately 5–7 pages:** (i) display eight Bell coordinates, parent RFA10 and real nine-atom formula; (ii) state the transferred polar and SG10 source reduction with exact event/noncoercive/native cases; (iii) prove TF1–TF2 in less than one page; (iv) state exactly `(NATIVE)` and point to complete scoped proofs; (v) prove parent restoration + zero-weight double assembly; (vi) quote Cor. 3.19 precisely and give the two-inclusion finite-moment pencil argument. Do not present an evaluated lift.

**Appendix A (foundations):** one-page Naimark/Jordan/real split; SG1 non-signalling/Farkas; SG5–SG10/SG17 source extrema and EN signs; SG11–SG14 closed-face LMIs and their attained duals.

**Appendix B (mechanisms):** TF1/TF2, SG15 sign/commuting homogeneous endpoints, full NC noncommuting-residue theorem, CC1–CC6 monic conductor and real allocation lemmas. Include hypotheses about **bounded coquotients** and the source/test factor split.

**Appendices C/D (independently reviewed local proof bodies):** native C1 and C2 full-parent sector proofs, exact local coefficient arrays, both normals, literal zero faces, all-testing joins, formal thin-R10 negative-beta paired-reserve and its physical source restoration where used. No global assertion is allowed as an import to a leaf that is used to prove that assertion.

**Separate artifact, not an appendix claimed complete:** evaluated fixed affine lift or global feature/pencil extraction. This is distinct from Cor. 3.19's qualitative existence conclusion.

## Final scope

No new unresolved *global implication* was found. The proof's conditional premise `(NATIVE)` remains exactly the source-coverage task being reviewed separately. No reliance is placed on a 215-component count, 45/49 charts, PASS labels or numerical samples as a substitute for `(NATIVE)`.
