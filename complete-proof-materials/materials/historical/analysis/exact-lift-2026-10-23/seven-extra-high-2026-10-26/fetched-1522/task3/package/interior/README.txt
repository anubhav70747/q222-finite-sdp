Q222 native strict C2 finite-sixth-root interior audit, 2026-10-26

PROVED: exact native source/parent, six full-state/angle contact jets,
154-feature support-face coefficient map, rank 406, rational affine Gram
preimage, complete polynomial target clearing, closed Gram-image cone,
exact normalized PSD dual criterion, strict physical angular curvature.

NOT PROVED: global strict C2 PSD admission in Z6, full nondegenerate
continuation, arbitrary finite native c, any singular source wall, or full
Q222 fixed rational lift. The missing inequality is (13) in the proof note.

Dependencies: Python 3; SymPy for symbolic files; NumPy for the optional
reproduction of the exact 406-mod-65521 rank certificate. No solvers,
numerical sweeps, GitHub connection, or absent binary intermediates.

Minimal exact verification from THIS directory:
    python c2_fixed_map.py --check
    python verify_native_jets.py
    python verify_rank_mod65521.py
    python verify_angular_curvature.py
    python generate_cleared_targets.py

Exact witness verification (supplied rational witness, no search):
    python verify_exact_certificates.py --primal MY_GRAM.json
    python verify_exact_certificates.py --dual MY_DUAL.json

To generate an affine-only rational preimage from a native rational source,
read generate_affine_preimage.py. Such a preimage has NO PSD guarantee and
rational 406-by-406 elimination can be slow.

C2_154_CONTACT_MAP.json: all 412 row keys, all 3021 sparse Gram coefficient
columns and six native contact-annihilator rows.
C2_RANK406_MODULAR_PIVOTS.json: original row/column indices of a 406-minor
that is provably nonzero modulo the prime 65521.
C2_CLEARED_POLYNOMIAL_TARGET.json: exactly 24 native integer-polynomial
nonzero target rows after multiplying by one positive source denominator.

Original independently verified open-end theorem files are copied, with
unchanged bytes, to preserved-open-end/. The previous theorem is strictly
local at the C2 c-end faces, not all finite native c.

Public current research repository cited by user (no GitHub writes):
https://github.com/anubhav70747/q222-finite-sdp/tree/1ad5125dbfc9c56881c283724b04994e3d50ec00
See independent-review/ASSEMBLY_REFINEMENT.md and
independent-review/checks/gram-lift/contact_rank.py there.
