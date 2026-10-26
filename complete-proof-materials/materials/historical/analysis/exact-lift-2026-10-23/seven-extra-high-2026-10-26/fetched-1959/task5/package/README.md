# C1 K364 direct-admission frontier, 27 October 2026

This package is a **complete exact symbolic proof of a source-wide local-in-Bob full-state PSD theorem**, not a proof of the global K364 admission H1.

- `PROOF.md` gives the original strict C1 native hypotheses, complete Bell parent, four exact positive endpoint SOS identities, a fixed 40-row globally sound *conditional* Gram mechanism, the exact source-rational Bob-zero skew, strictly positive **joint 2x2 Schur**, and a source-algebraic exact contact interpolation preserving both projective endpoints.
- `verify_full_parent_bob.py` verifies the four-state coefficient polynomial and all four endpoint identities, the generic 8x8 skew coefficient identity, and feature grid inclusion.
- `verify_boundary_schur.py` verifies the full 6+2 graded Schur calculation with arbitrary remaining higher skew entries, allowing singular/contact limits. **No numerical sources or SDP.**
- `verify_anchor_interpolation.py` verifies interpolation and source-only denominators for arbitrary skew contact data.
- `verification.log` is the deterministic replay output (all pass).
- `source_refs/c1_native.py` is the unchanged original seven-coquotient native C1 source and full parent normalization.
- `source_refs/PREVIOUS_Z18_SCOPE.md` is the untouched preceding proof for context; it does not prove H1.

To replay: `python run_exact_checks.py` (requires SymPy). No GitHub writes.

**Exact unresolved inequality:** global `H_theta(v) >= 0` for every real `v>0` with one source-only degree-four skew polynomial; this is sufficient for H1 but not necessary for membership in the full 364-feature cone. Both the universally oriented strict C1 K364 admission H1 and strict C2 H2 remain separate open statements. No counterexample to K364 is asserted.
