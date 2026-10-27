# Preservation audit after the Sun–Wang–Yu update

27 October 2026. Audited head: [`ac8861a`](https://github.com/anubhav70747/q222-finite-sdp/commit/ac8861aab22778ee06dbbf79c1784898e1e550f8).

**No tracked file was deleted, and no manuscript, PDF, archived proof, extracted lift or existing exact-checker file was changed in the audited ranges.** This conclusion compares Git blobs, independently of Windows checkout line endings. The [script](verify_preservation.py) and [complete receipt](result.json) make the two comparisons reproducible.

| Baseline | Existing files unchanged | New files | Deleted files |
|---|---:|---:|---:|
| Immediately before the matching update, `1b29383` | 2,484 of 2,491 | 4 | 0 |
| Original complete-material publication, `3324833` | 1,930 of 1,943 | 552 | 0 |

The latest update changed only `FRONTIER.md`, the two checklist formats, the manifest, the delivery receipt, the verification plan and the independent-review index. **Every pre-existing checklist item is identical to `1b29383`; L22 was added.** All 215 component records, 45 C1 cell records and 49 C2 leaf records are identical at both baselines and the audited head.

Across the longer range, the 13 modified files are indices, status/reconciliation records, and `EXACT_LIFT_CANDIDATE.md`. That candidate note has added passages only: none of its earlier lines was deleted. The original manuscript, mathematical packages and their artifacts remain byte-identical. The changed earlier checklist entries are L08–L10 and V01–V04/V06: their later evidence and scope are retained, and no item identifier was removed. Both baselines are ancestors of the audited head and remain recoverable in history.

The manifest diff looked larger because 14 pre-existing text entries gained an explicit `canonical_lf` flag. Their recorded SHA-256 hashes and byte lengths did not change. Those values already matched the Git blobs; the flag makes the verifier handle Windows CRLF checkouts consistently. The five hashes changed by the matching update belong to the edited frontier/checklist/plan/index documents, not to prior proof data.

This audit establishes preservation through the stated commit. It does not certify every mathematical claim afresh or cover untracked material never present in these Git snapshots. The separate reviewed qualitative V01–V06 conjunction is retained, and the explicit construction still has open H1/H2 obligations. No rollback or restoration was necessary.

Reproduce with `python reconciliation/2026-10-27/preservation-audit/verify_preservation.py`. The default head is frozen to the audited publication; `--head <commit>` permits a later explicitly scoped comparison. Scheduled browser monitoring was not resumed.
