# Strict oriented C1: exact fixed-Z18 coefficient and PSD frontier

Original analytical theorem: `PROOF.md`.

New source-independent mechanism: the exact zero identity
`(r2-r4)(r3-r5)=u*v*(X-Y)^2`, which adds a second source-only half-power gauge.
Exact map: 46 generated coefficient keys, 49 formal coordinates, rank 40;
after forbidden parity terms are removed, 42 native invariant keys, 45 coordinates,
rank 36 and 9 gauges. The full native parent is matched exactly.

Exact positive pivots reduce the remaining Gram admission to two 4x4
matrices with gamma,delta>=0. This is a conditional PSD criterion, not a
universal PSD proof. In particular the k>1 obstruction is not asserted on
oriented 0<k<=1.

Run from this directory:

    python verify_z18_exact.py
    python verify_schur_identities.py
    python original_baseline/verify_c1_native_stationarity.py
    python original_baseline/verify_c1_map_exact.py

The checkers perform deterministic SymPy integer/rational polynomial operations.
No numerical sweep, approximate PSD matrix, optimizer or repository mutation.
