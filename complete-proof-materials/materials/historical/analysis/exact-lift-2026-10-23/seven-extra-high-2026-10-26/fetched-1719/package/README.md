# C1 reflection transport and fixed cone, reproducible artifacts

Frozen input: `C1_REFLECTION_SOURCE_PACKET.txt` (provided with task). No GitHub writes.

- `C1_REFLECTION_TRANSPORT_AND_FIXED_CONE_THEOREM.md`: original complete conditional proof, full source map, 16 reflected features, fixed 364-row cone, closedness, endpoint validity and exact dual-adjoint assembly.
- `verify_reflection_transport.py`: deterministic SymPy exact assertions, no numerical sweeps. Run `python verify_reflection_transport.py`.
- `transport_sparse_maps.json`: generated fixed source/Bell maps, 364 feature indices, 32 signed reflected-feature entries, 196-to-364 integer embeddings, and native map.
- `EXACT_CHECKS.log`: saved exact-checker stdout.
- `C1_REFLECTION_SOURCE_PACKET.txt`: frozen user-supplied source excerpt.

**Conditional inputs NOT checked here:** universal oriented `0<k<=1` strict C1 `Z16` PSD, and universal strict C2 `Z6` PSD. The published source-contact density theorem, SG geometry/entry and the independent strict full-`Z6` feasibility point are cited as established imports. No claim of unconditional explicit SDP exactness.
