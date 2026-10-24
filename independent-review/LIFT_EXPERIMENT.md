# Is the ramified cone plus its mirror an exact lift?

The repo reduces an exact lift of Q222 to an exact semidefinite description of
the 9-dimensional canonical cone `C_can` of valid sources
`c = (d, λ, n, m, C, A, D, α, β)`. The order-49 ramified cone `C_ram` is known to
be strictly smaller, because its Alice↔Bob mirror source `c†` is valid but
excluded. The natural repair is the order-98 cone `C_ram + πC_ram`, whose
eight-copy assembly would be an order-784 candidate for Q222. Here
`π` swaps λ↔C, n↔A and α↔β, and fixes d, m and D.

**Finding.** `C_ram + πC_ram` is strictly smaller than `C_can`, and by a wide
margin on multiscale sources. This is established numerically, with
high-precision solves and one dual certificate in ball arithmetic. It is not yet
an exact rational certificate.

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

- Only sym 49 has a dual certificate. The other brackets come from
  high-precision primal–dual solves, with duality gaps of at most 10⁻⁷.
- The C1 and C2 face samples replace exact zeros by 10⁻¹⁸.
- The float screen is not reliable across all 302 vectors, so lower ratios may
  exist among the unverified ones.
- Turning the sym 49 or sym 36 separation into an exact rational certificate
  needs two pieces. The first is a rational dual matrix that is exactly positive
  semidefinite. The second is an exact AM–GM validity certificate for a rational
  `d` between `d_sum` and `d*`. Both are routine follow-ups.

The high-precision results are recorded in
[`checks/lift-experiment/hp_results.jsonl`](checks/lift-experiment/hp_results.jsonl).
