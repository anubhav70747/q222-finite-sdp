# Is the ramified cone plus its mirror an exact lift?

The repo reduces an exact lift of Q222 to an exact semidefinite description of
the 9-dimensional canonical cone `C_can` of valid sources
`c = (d, λ, n, m, C, A, D, α, β)`. The order-49 ramified cone `C_ram` is known to
be strictly smaller, because its Alice↔Bob mirror source `c†` is valid but
excluded. The natural repair is the order-98 cone `C_ram + πC_ram`, whose
eight-copy assembly would be an order-784 candidate for Q222. Here
`π` swaps λ↔C, n↔A and α↔β, and fixes d, m and D.

**Finding.** `C_ram + πC_ram` is strictly smaller than `C_can`. An exact
rational certificate proves this for one source (next section).
High-precision solves show the gap is wide on multiscale sources.

## Exact certificate

The rational source, in the order `(d, λ, n, m, C, A, D, α, β)`, is

    c_sep = (481/10⁶, 1, 123/10⁵, 11/(125·10⁶), 1, 123/10⁵, 17/(25·10⁶), 89/10⁹, 89/10⁹),

a rationalized version of "sym 49" below with `d ≈ 0.88·d*`. It is π-symmetric.
[`checks/lift-experiment/verify_cert.py`](checks/lift-experiment/verify_cert.py)
checks the following in exact rational arithmetic, from
[`cert_data.json`](checks/lift-experiment/cert_data.json) and the repo's pencil.
All 17 checks pass.

**(A) c_sep is valid.** The determinant is divided exactly into the
`u²v²·L(u²,v²)` form. Five AM–GM circuits with rational coefficients cover the
outer coefficients. Each circuit minimum `Π(c_e/β_e)^{β_e}` has a rational lower
bound `m_k`, proved by comparing integer powers. Then
`h + Σ m_k = 1.845·10⁻⁷ > 0`, so `L > 0` for all `U, V > 0`.

Exact LDLᵀ shows `P_c` is positive definite at one point in each sign quadrant.
`det P_c` never vanishes inside a quadrant, so `P_c ≻ 0` on every open quadrant,
and `P_c ⪰ 0` everywhere by closure. This part does not rely on the paper's
"P ⪰ 0 ⇔ det ≥ 0" lemma.

**(B) c_sep is not in C_ram + πC_ram.** Rational 49×49 matrices `Y` and `Y′`
follow the pencil's block structure. Their four blocks are positive definite by
exact LDLᵀ, and their nine scalar entries are nonnegative. Both annihilate all
148 auxiliary matrices exactly. Together they define one functional `ℓ` with
`ℓ_k = Tr(A_k Y)` and `ℓ_{π(k)} = Tr(A_k Y′)`. So `ℓ ≥ 0` on `C_ram` and on
`πC_ram`, hence on their sum. But `ℓ(c_sep) = −5.548·10⁻⁵ < 0`.

To run it: `python verify_cert.py <repo>/extracted-lifts/canonical-completion/ramified-completion/data/ramified_order49_pencil.json`. It takes about 5 s.

## Method

- **Validity.** For a vector `e` of event weights, `d*(e)` is the largest CHSH
  weight with `(d, e) ∈ C_can`. Because the 4×4 state matrix has one violating
  direction, validity is `h_c + φ(q_c) ≥ 0`, where
  `φ(q) = min over z ∈ ℝ² of Σ q_e·exp(e·z)` is a convex minimization (see
  `PEN_AND_PAPER.md`, steps 3 and 5). The test agreed with a brute-force
  minimum-eigenvalue scan of `P_c(u,v)` at 28 of 28 points placed at
  `d*(1 ± 0.3%)`.
- **Cones.** `d_ram(e)` and `d_sum(e)` are the largest `d` with `(d, e)` in `C_ram`
  and in `C_ram + πC_ram`. They are computed by semidefinite programs on the
  repo's pencil `ramified_order49_pencil.json`.
- **Reconstruction check.** The pencil was reconstructed and checked: it is
  exactly positive semidefinite at both of the repo's feasible points.
  `d_ram(c̄) = d*(c̄)`, `d_ram(c†) = 0.8454·d*`, and `d_sum(c†) = d*`, as
  expected.
- **Precision.** 302 event vectors were screened in float64: random,
  π-symmetric, the C1 face and the C2 face. Float solvers are unreliable here, so
  every candidate gap was re-solved with a 192–320-bit log-barrier method, with
  primal iterates checked in ball arithmetic.

## Results

High-precision ratios are given as brackets of width about 1.5·10⁻⁸. Event order
is `(λ, n, m, C, A, D, α, β)`.

| source | event weights | d* | d_sum / d* | d_ram / d* |
|---|---|---|---|---|
| sym 36 | (0.26924, 5.5224, 6.09e-7, 0.26924, 5.5224, 3.49e-6, 2.34e-7, 2.34e-7) | 0.031173 | **0.58827** | 0.55552 |
| C2 2 | (0, 4.60e-4, 0.7947, 46.03, 1.30e-4, 33.54, 14.63, 0) | 0.66085 | 0.73432 | 0.73432 |
| C1 7 | (0, 0, 0.08741, 1.518, 4.42e-4, 0.4237, 5.114, 5.00e-7) | 0.13132 | 0.75946 | 0.75946 |
| sym 49 | (18.341, 0.022574, 1.62e-6, 18.341, 0.022574, 1.25e-5, 1.63e-6, 1.63e-6) | 0.010036 | **0.76031** † | 0.72648 † |
| random 33 | (see data) | 0.021824 | 0.96963 | 0.96963 |
| random 91 | (see data) | 0.084907 | 0.98397 | 0.98397 |
| C2 38 | (see data) | 0.047534 | 0.99030 | 0.91228 |

† Dual matrix verified positive semidefinite in 256-bit ball arithmetic, with a
final dual residual of about 10⁻⁷⁵.

In each gap case, a brute-force eigenvalue check confirms that `P_c ⪰ 0` holds
at the midpoint of `(d_sum, d*)`, and fails just above `d*`. So these sources are
valid but lie outside `C_ram + πC_ram`. Sixteen other high-precision solves
returned `d_sum = d*` to within 1.4·10⁻⁸.

## What it says about the exact lift

- **The order-784 candidate is not exact.** Its eight-copy assembly is strictly
  larger than Q222.
- **The missing functions are about scale, not symmetry.** In four of the seven
  gap cases the mirror adds nothing (`d_sum = d_ram`). Every gap source has
  event weights spread over 5–7 orders of magnitude. Selectors built from fixed
  powers of Bob's angle alone cannot follow such multiscale contacts.
- **A testable next family.** The separated AM–GM certificates of the existence
  proof use cube roots of `u²` and `v²`, and the square roots of those cube
  roots, in both angles (`PEN_AND_PAPER.md`, step 5). The ramified family uses
  only `√v`. A candidate testing space should add `u^{1/3}`- and
  `v^{1/3}`-type functions in both angles. It should be tested on the seven
  sources above before anything else. This is a hypothesis, not a result.

## Caveats

- Apart from the exact certificate, the brackets in the table come from
  high-precision primal–dual solves, with duality gaps of at most 10⁻⁷.
- The C1 and C2 face samples replace exact zeros by 10⁻¹⁸.
- The float screen is not reliable across all 302 vectors, so lower ratios may
  exist among the unverified ones.
- The exact certificate takes the repo's pencil file as the definition of
  `C_ram`.

The high-precision results are recorded in
[`checks/lift-experiment/hp_results.jsonl`](checks/lift-experiment/hp_results.jsonl).
