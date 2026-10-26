STRICT C2 DUAL MOMENT PAIRING AUDIT, 2026-10-26.

1) STRICT_C2_24_TERM_PAIRING_AUDIT.md: full analytical signed-affine result,
   exactly specified attempted pairing, native signs and decisive failure
   of the chosen pairing architecture. Global nonstrict delta>=0 UNPROVED.
2) verify_corner_pairing.py: standalone exact deterministic checker and
   generator. Requires only SymPy and bundled files.
3) C2_CORNER_AFFINE_18.json: source-independent fixed 18x18 rational
   inverse and 91x18 contact Gram data, AFFINE ONLY.
4) C2_AFFINE_BRIDGE_17_PLUS_1.json: full 17-edge plus one bridge
   coefficient recipe, AFFINE ONLY; bridge has indefinite 4x4 compression.
5) C2_154_CONTACT_MAP.json, C2_CLEARED_POLYNOMIAL_TARGET.json and
   c2_fixed_map.py: verbatim copies of prior independently checked baseline.

Verify:
 python verify_corner_pairing.py

There are no SDP solves, numerical sweeps, testing-dependent inverses,
missing .pkl data, or claims of C2 interior L6 membership.
The fixed 154-feature sixth-root dual moment inequality remains open.
