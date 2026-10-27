# C1 K364: source-wide nondegenerate contact jet, not H1

**Exact scope, 27 October 2026.** For every genuine strictly positive native C1 source with oriented `0<k<=1` and the unmodified full four-state parent, there exists a **source-algebraic 2x2 polynomial**

`T_theta(v)=T_c+(v-v_c)T_1+(v-v_c)^2 T_2`

such that the full coefficientwise identity

`(1+v^2)F_theta(u,v,psi)=(psi,u psi)^T H_theta(v)(psi,u psi)`

holds at every pair of testing angles and every state; its `8x8` matrix `H_theta(v)` splits into **two complete 4x4 parity blocks**, has **rank 6** at `v=v_c`, and is **positive definite at every v in a two-sided punctured neighbourhood** of `v_c`. It is not known to be PSD for all `v>0`, so it does not certify H1 or K364 membership.

**Mechanism:** Calderon's 2xn SOS theorem plus the *unique* native C1 contact zeros and strict contact Hessian produce a rank-six source-only contact Gram (not merely a PSD anchor). At contact its cross-parity skew belongs to an exact 1-parameter affine family `T_c=T_0+tau w_perp y_perp^T`. A source-only `2x2` positive definite inverse solves the linear Bob-jet minimization. The strict SG18-SG19 Hessian gives a positive scalar total second-order Schur curvature. One more source-only skew entry splits that margin equally between the two soft modes, proving two-sided PSD by an ordinary full matrix Schur calculation. See PROOF.md §§3–5.

**Remaining exact sufficiency condition:** choose one source-only *global* polynomial `T(v)` of degree at most four making **both** complete 2x2 Schur residuals `R_A(v),R_B(v)` of PROOF.md Eq. (5) PSD for every `v>0` (with homogeneous Bob endpoints), or find a different entire coefficientwise 364-feature PSD Gram. This is stronger than H1; its failure does not refute K364. No universal gauge selection, no global Gram, and no exact separating functional were found in this assignment.

**Exact evidence:** `verify_contact_geometry.py`, `verify_two_schur.py` and `verify_contact_jet_exact.py` are symbolic/rational, solver-free; the last verifies one already-known native seed as a targeted regression, not a source survey. The universal positivity proof is mathematical and uses the earlier verified SG6, SG18–SG19, original native source classification, and Calderon as an explicit imported theorem. The previous endpoint and conditional 40-feature proofs are preserved unchanged under original/.
