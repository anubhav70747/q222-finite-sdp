# Complete exact conditional K364 global SDP package

**Generated and verified 26 October 2026.** The included sources contain exact rational matrices and are self-contained with SymPy installed. The *outer SDP is unconditional*; exactness for Q222 awaits separate universal strict oriented C1 and strict C2 direct admissions to K364.

Read `PROOF_AND_MATRIX_SPECIFICATION.md` for theorem, every matrix indexing and two-inclusion proof. `COORDINATOR_SHORT_THEOREM.md` is the short handoff.

Materialized data:
- `K364_FULL128_FACTORED_SPARSE_PENCIL.json` (and gz): complete 364-feature basis, 3250 coefficient keys, all 66430 Gram columns, sparse T, SG1 basis, 128 action matrices, every target, face source PSD bases, copy offsets.
- `K364_GLOBAL_ALL_3456_AFFINE_EQUALITIES.json` (and gz): each scalar equality as literal rational coefficients and absolute variable indices.
- `K364_ALL_EVENT_LMI_ENTRIES.json`: every physical-event diagonal LMI.
- `K364_STRICT_RATIONAL_SOURCE_AND_GRAM.json`: exact 196/364 strict rational witness.

Scripts, in dependency order:
- `generate_exact_pencil.py` (stdlib + SymPy)
- `assemble_global_pencil.py` (stdlib)
- `verify_exact_pencil.py` (stdlib + SymPy, invokes the two local modules)

Run `python generate_exact_pencil.py && python assemble_global_pencil.py --check-structural && python verify_exact_pencil.py`.
Optional `python assemble_global_pencil.py --stream-all-lmi` unfolds every copied LMI entry to compressed JSONL. Not needed for reproducibility since the shared matrices and all 128 explicit copy offsets already evaluate every entry uniquely.

No numerical optimization, GitHub writes or sampling. The geometric/cone sources are the uploaded `batch1739_task7_verified_sources.txt` and the checked frozen main source.
