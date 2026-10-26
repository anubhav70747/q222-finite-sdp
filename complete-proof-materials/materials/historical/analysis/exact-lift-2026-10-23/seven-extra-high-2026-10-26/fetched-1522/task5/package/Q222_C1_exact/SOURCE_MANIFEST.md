# Immutable source and exact dependency manifest

The analytical results in `PROOF.md` use only the following primary identities, plus symbolic algebra derived directly from them:

1. Original physical native C1 entry and complete four-state source (SD1--SD11): [c1_short.tex at 1ad5125dbfc9c56881c283724b04994e3d50ec00](https://github.com/anubhav70747/q222-finite-sdp/blob/1ad5125dbfc9c56881c283724b04994e3d50ec00/latex/certificates/c1_short.tex). The `k,f,g,a,t` coquotients and state contact are transcribed without numerical substitutions in `c1_native.py`.
2. Original complete canonical source parent and global stationarity-to-validity theorem: [04_geometry.tex](https://github.com/anubhav70747/q222-finite-sdp/blob/1ad5125dbfc9c56881c283724b04994e3d50ec00/latex/chapters/04_geometry.tex), and corresponding source parent in `latex/certificates/c1_short.tex`.
3. Current fixed sixth-root candidate and its *conjectural* universal C1/C2 contact admission: [EXACT_LIFT_CANDIDATE.md](https://github.com/anubhav70747/q222-finite-sdp/blob/1ad5125dbfc9c56881c283724b04994e3d50ec00/independent-review/EXACT_LIFT_CANDIDATE.md).
4. Contact-rank/affine assembly and caution about connectedness without closed nondegeneracy: [ASSEMBLY_REFINEMENT.md](https://github.com/anubhav70747/q222-finite-sdp/blob/1ad5125dbfc9c56881c283724b04994e3d50ec00/independent-review/ASSEMBLY_REFINEMENT.md).
5. Independently verified **strict C2 open-end only** construction at [commit 1ad5125](https://github.com/anubhav70747/q222-finite-sdp/commit/1ad5125dbfc9c56881c283724b04994e3d50ec00). The present C1 proof reuses its *idea* of source-only gauge and fixed image-inverse correction, not its C2 matrices, source signs or PSD conclusion.
6. User supplied consolidated evidence `/mnt/data/05_C1_FIXED_ADMISSION_SOURCES_1418.txt`, also preserved verbatim in this package under `sources/`.

The proof is valid independently of any numerical SDP evidence and does not rely on the user's new C2 theorem to assert C1 positivity. `c1_small_map.py`, `verify_c1_map_exact.py` and the other checker sources regenerate all claimed integer and rational identities, with Python3 and SymPy.

No GitHub mutation was made.
