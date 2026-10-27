# Complete proof materials and checklist

This collection publishes the mathematical material behind the Q222 work: original source archives, individual proof bodies, refined sources, coefficient and matrix data, symbolic checkers, recorded outputs and analytical reviews.

Start with the [complete checklist](CHECKLIST.md). It links the foundations, every local case, the explicit constructions and the remaining verification tasks to their actual evidence.

- [215 proof components](COMPONENTS.md), with original and presentation sources separately linked.
- [45 C1 cases](C1_CASES.md) and [49 C2 cases](C2_CASES.md), retaining their source, testing and state hypotheses.
- [Material inventory](INVENTORY.md), including the original full R21 archive and the newly fetched complete proof/checker packages.
- [Existing explicit representations and matrix data](../extracted-lifts/), [current paper sources](../latex/) and [current PDFs](../pdf/).
- [Status reconciliation](RECONCILIATION.md), [recovery limits](RECOVERY_LIMITS.md) and [next verification pass](VERIFICATION_PLAN.md).
- [Machine-readable checklist](CHECKLIST.json), [source-to-publication map](SOURCE_MAP.json), [file hashes](MANIFEST.json) and [delivery verification](VERIFICATION.json).

The full 22 October analytical review supports the qualitative existence theorem. Later scoped reviews raised narrower questions about particular versions and interfaces. Those questions are recorded without treating previously reviewed cases as absent. The current-version qualitative conjunction is reconciled in V01–V06 using those body reviews. L18 now supplies the complete evaluated conditional global outer SDP. Its unconditional exactness still requires universal oriented C1 and strict C2 admission into the fixed cone, followed by the proved density and assembly implication. See the [current frontier](../FRONTIER.md).

Historical files retain their original wording and dates, including superseded constructions, unsuccessful attempts and earlier status claims. Their presence is evidence of the work performed, not endorsement of every claim in every version. The reconciliation and checklist identify how to use them.

Run `python complete-proof-materials/verify_materials.py` from the repository to check the delivered bytes, original component/dependency hashes and checklist links. This verifier performs no SDP solves and makes no mathematical correctness inference from file integrity.
