# Native C2 R2 source and all-testing audit

Pinned input: `anubhav70747/q222-finite-sdp`, commit `8960b3fead3ce6b185bb1387d148fb87d05d3d23`.

`R2_C2_COEFFICIENT_AND_ASSEMBLY_PROOF.md` is the complete scoped argument. Its R2 theorem uses the **exact algebraic** RFN/SH/SC/QI/SJ/R2O lemmas with precisely the hypotheses listed in section 1. No full C1 or full C2 physicality theorem is used on a formal negative-beta source. The independent UL proof handles positive-unbounded relative current. The last section gives an exact positive native R8 source, proving why the currently checked leaf list alone does **not** close global C2.

Checks use exact SymPy rational polynomial operations, with no solver and no numerical sweep:

- `verify_r2.py`: 8 independent original-source TU/NC4 polynomial checks, including all six native/formal coefficient differences.
- `verify_rl_fast.py`: compact RL RSR5 zero-p numerator and all four native RFN21–22 source balances by **polynomial numerators**, not heavy rational expansion.
- `verify_full_state.py`: generic full four-state RFN7 readback, both hard diagonals and cross normal, and the general harmonic shift identity.
- `verify_circuit.py`: SH1, SH2, RL17, TU17, and the SH4 numerator factorization.
- `verify_outer_caps.py`: all R2O21–R2O24 finite-cap leading/minor identities.
- `verify_r8_scope.py`: one formal positive rational-Puiseux R8 source, not a numerical search.
- `verification.log`: successful exact output from each script. (Scripts were run individually; a combined shell invocation exceeded its execution time limit after the first scripts, so no unverified aggregate-success claim is made.)

Each script can be invoked separately with `python filename.py`. Dependencies: SymPy 1.14, Python 3.11+.

**Status:** R2 local coefficient and full-state proof established relative to specifically named algebraic conductors; native C2 end-to-end depends on independently checked R8 and separate one-angle/strong/two-weak leaves, and full Q222 existence additionally depends on C1 and global assembly. No fixed-degree source-independent lift is claimed.
