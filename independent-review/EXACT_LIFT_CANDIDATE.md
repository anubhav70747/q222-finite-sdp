# A candidate exact lift: the sextic Puiseux–Gram lift

This note proposes an explicit semidefinite lift of the canonical cone `C_can`,
and hence of Q222, and records the numerical evidence for it. The lift is a
single fixed SDP: one Gram matrix in a fixed basis of fractional monomials in
the two measurement angles. It has no charts and no infinitesimals.

Its exactness is a conjecture. It is supported by 302 test sources, including
every source that defeats the earlier candidates, but it is not proved.

## The lift

Let `c = (d, λ, n, m, C, A, D, α, β) ≥ 0` be a source and `F_c(u,v,ψ)` its parent
form (see `PEN_AND_PAPER.md`, step 3). For an integer `k ≥ 1`, set

    Z_k = { u^{i/k} · v^{j/k} · ψ_l  :  0 ≤ i, j ≤ k,  l = 1..4 },

and define

    L_k = { c : there is G ⪰ 0 with F_c(u,v,ψ) = Z_kᵀ G Z_k for all u, v > 0 }.

- **It is an SDP.** Distinct fractional monomials are linearly independent on
  `(0,∞)²`, so the identity is a finite set of linear equations in `(c, G)`. So
  `L_k` is the projection of a spectrahedral cone.
- **It is sound.** A Gram certificate gives `F_c ≥ 0` for `u, v > 0`. The exact
  symmetries `(u,Y,Z) → (−u,−Y,−Z)` and `(v,X,Z) → (−v,−X,−Z)` extend this to
  every quadrant, and continuity covers the axes and endpoints. Hence
  `L_k ⊆ C_can`.
- **It is nested.** `L_k ⊆ L_{k′}` whenever `k` divides `k′`.
- **Its size.** A sign symmetry splits the Gram matrix into blocks.

| k | basis | Gram blocks | free auxiliaries per copy |
|---|---|---|---|
| 1 | {1, u}⊗{1, v}⊗ψ | 4, 4, 4, 4 | n/a |
| 2 | {1, √u, u}⊗{1, √v, v}⊗ψ | 16, 8, 8, 4 | 128 |
| 6 | sixth-root grid | 64, 48, 48, 36 | 4608 |

The full sixth-root system has 5098 symmetric Gram variables and 490
coefficient equations. On either open C1 or C2 support face, forced zero
rows reduce it to blocks of orders 52, 36, 36, 30: 3175 Gram variables,
412 equations and 2763 free auxiliaries. The smaller count applies only
after the source face has been fixed.

**Conjecture.** `L₆ = C_can`.

## What is proved

**Proposition 1 (soundness).** `L_k ⊆ C_can` for every k. See the bullets
above.

**Proposition 2 (closedness).** Each `L_k` is a closed convex cone.

*Proof.* Let μ be the uniform probability measure on `[½, 2]² × S³`. Distinct
fractional monomials `u^{i/k}v^{j/k}` are real-analytic functions with distinct
exponents. They are therefore linearly independent on any open subset of
`(0,∞)²`, so `M = ∫ Z Zᵀ dμ` is positive definite. If `F_c = Zᵀ G Z`, then
`Tr G ≤ Tr(G M)/λ_min(M) = ∫ F_c dμ / λ_min(M)`, and the right-hand side is
linear in `c`. Take a convergent sequence `c_n → c` with certificates `G_n`. The
`G_n` are then bounded. A subsequence converges to some `G ⪰ 0`, and the linear
identity passes to the limit. ∎

**Proposition 3 (reduction to extreme sources).** `L₆ = C_can` if and only if
`L₆` contains every extreme point of the compact base
`B = {c ∈ C_can : d + λ + n + m + C + A + D + α + β = 1}`.

*Proof.* `L₆ ⊆ C_can` is closed and convex. By Krein–Milman, `B` is the closed
convex hull of its extreme points. ∎

**Proposition 4 (events).** Every single-event source lies in `L₁ ⊆ L₆`.

*Proof.* An event term is `c_k·(h_k·ψ)²`, where `c_k` is one of
`1, D_u, D_v, D_uD_v`, and `h_k·ψ` has degree at most 1 in each of `u` and `v`.
Expanding `D_u = 1² + u²` and `D_v = 1² + v²` writes it as a sum of squares of
elements of the span of `Z₁`. ∎

By the classification of extreme sources (paper, Ch. 4), the conjecture is
therefore **equivalent** to one statement:

> **Remaining lemma.** The extreme sources of the two closed faces, and every
> C1 and C2 contact source, have a Gram certificate in `Z₆`.

This lemma is over the reals and involves no infinitesimals. The contact sources
form explicit rational families, parametrized by their contact point
(`c2_short.tex:263–292`, `c1_short.tex:696–748`). A proof could take the form
of a Gram matrix whose entries are rational functions of those parameters,
verified symbolically.

The determinant's circuit types alone do not predict which sources `L₂`
certifies (see "Why sixth roots, physically"). A construction using the
circuits would also have to recover the state matrix and its coefficients.
The numerical comparison does not rule out such a construction.

## Exact certificates at extreme sources

Four genuine boundary sources have exact rational certificates in `L₆`:
two C1 sources (C1A, C1B) and two C2 sources (C2A, C2B). Each is a contact
source from the paper's rational families, with `d = d*` and every non-forced
weight strictly positive. The contact angles are chosen to be sixth powers of
rationals, which makes the contact vector rational.

The checker ([`checks/gram-lift/boundary/`](checks/gram-lift/boundary/), 56/56
exact checks) proves three things for each source:

1. **Valid, on the boundary.** `F_c ≥ 0`, by a weighted AM–GM bound on the
   determinant plus exact positive-definite points in every quadrant, and
   `F_c` vanishes with zero gradient at the contact.
2. **Exact Gram identity.** `F_c` equals a Gram form in `Z₆` exactly, as
   Puiseux polynomials.
3. **Positive definite after reduction.** First remove the 42 monomials forced
   out by the source support face. In each remaining parity block, remove the
   contact vector `z_c`. The resulting four Gram matrices are positive definite.

At all four sources the certificate has the same structure:

- **Kernel.** In the 154-function support-face basis there is exactly one
  contact vector in each parity block. Embedded in the full 196-function
  basis, the Gram matrix also has the 42 forced null directions, so its rank
  is 150 and its kernel has dimension 46.
- **Rank.** After contact reduction, the coefficient map has 412 rows and
  exact rank 406. The six missing dimensions are the four state derivatives
  and two angle derivatives; homogeneity makes the value condition
  redundant. The exact rank calculation and its extension to all nonzero
  contacts are in [`ASSEMBLY_REFINEMENT.md`](ASSEMBLY_REFINEMENT.md).

## A continuation theorem

For a support pattern σ (C1: λ = n = 0; C2: λ = β = 0), let `Θ_σ` be the open
orthant of its six remaining weights, and for `w ∈ Θ_σ` let `c(w) = (d*(w), w)`.
`Θ_σ` is connected, and the contact `x_c(w)` depends analytically on `w`
(paper, Ch. 5). Fix the 154-function support-face basis described above.
Call a certificate for `c(w)` **nondegenerate** if its four Gram matrices
are positive definite after their contact vectors have been removed.
The reduced coefficient map has rank `412 − 6 = 406` throughout the native
charts with nonzero contact coordinates, by monomial scaling and the exact
rank calculation. This definition is relative to the support face. It
cannot use the full 196-function basis: for example `λ=0` forces the row
of the constant feature `W` to vanish before any contact is imposed.

**Lemma (openness).** If `c(w₀)` has a nondegenerate certificate, so does
`c(w)` for every `w` near `w₀`.

*Proof.* Every reduced Gram form vanishes with zero gradient at `x_c(w)`, so the
image of the reduced map lies in the space `C_w` of forms that do. The six
contact conditions are independent at `w₀` (rank `N − 6`), hence near `w₀`, so
`dim C_w = N − 6`. By lower semicontinuity the rank stays at least `N − 6`, so
it equals `N − 6` and the image is all of `C_w`. Now `F_{c(w)} ∈ C_w`, because
`c(w)` has a contact at `x_c(w)`. Choose continuous bases of the contact
complements and write the reduced coefficient map as `A(w)`. Its
constant-rank pseudoinverse corrects the existing positive definite
certificate `G₀` by
`G(w)=G₀+A(w)†(F_{c(w)}−A(w)G₀)`. The correction tends to zero at `w₀`,
so each reduced Gram remains positive definite nearby. ∎

**Theorem (conditional).** Assume:

> **(ND)** For σ = C1 and σ = C2, every `w ∈ Θ_σ` that is a limit of points
> with nondegenerate certificates itself has a nondegenerate certificate.

Then every C1 and C2 contact source lies in `L₆`.

*Proof.* The set `S_σ` of `w ∈ Θ_σ` with a nondegenerate certificate is open by
the lemma, closed in `Θ_σ` by (ND), and nonempty by the seeds C1A and C2A.
Since `Θ_σ` is connected, `S_σ = Θ_σ`. The C1 and C2 sources with `D = 0`,
`β = 0` or `n = 0` are limits of these sources, so they lie in the closed cone
`L₆` too. ∎

With Propositions 1–4, **`L₆ = C_can` follows from (ND) together with the
extreme sources of the two closed faces.** The closed faces are noncoercive,
with no finite contact, so the continuation argument does not cover them.
Numerical tests support membership on 18 of 20 sampled sources; the other
two results are inconclusive at the available precision. Membership of the
closed faces in `L₆` is unnecessary for an exact hybrid lift: their existing
exact pencils can be added directly, as explained in
[`ASSEMBLY_REFINEMENT.md`](ASSEMBLY_REFINEMENT.md).

**The rank half of (ND) holds everywhere [proved].** The reduced coordinates
`p_i = z_i/z_{c,i} − z_{i0}/z_{c,i0}` are invariant under the torus rescaling
`u → u_c·u′`, `v → v_c·v′`, `ψ_l → ψ_{c,l}·ψ′_l`: the basis picks up a diagonal
factor that cancels in each ratio. C1/C2 contacts have all four amplitudes
nonzero (paper, strict sign lemma). So the reduced coefficient map at any
contact equals the map at the standard point `(1, 1; 1, 1, 1, 1)`, up to a
diagonal rescaling of coefficients. Its rank is therefore constant on each
family. In the 154-function support-face basis it is `412 − 6 = 406`, proved
exactly in [`checks/gram-lift/contact_rank.py`](checks/gram-lift/contact_rank.py)
(see [`ASSEMBLY_REFINEMENT.md`](ASSEMBLY_REFINEMENT.md)).

What (ND) still needs is strict positivity of the reduced Gram matrix.

**The forced zeros, physically.** `f(v_i) = 4c_i`, so a zero event weight means
the inequality is tight on a local deterministic strategy. Such strategies sit
at commuting angles, the corners of the Newton square. Newton pruning removes
exactly these, and the seeds show nothing else is forced besides the contact
vector.

**A fixed-cone formulation.** After the rescaling, every contact parent `F̃`
lies in one fixed space `C_std` of dimension `N − 6`. Let `K_std` be the fixed
cone of reduced Gram forms at the standard point.

- Exactness on the C1/C2 families is equivalent to `F̃ ∈ K_std` for every
  contact source.
- By continuation, it is enough that no `F̃` meets `∂K_std`.
- A boundary point needs a nonzero positive pseudo-moment `Λ` with `Λ(F̃) = 0`.
  Point evaluations inside the domain are excluded, because the contact is the
  unique zero there. Corner evaluations are removed by pruning.

So the remaining statement is a **sextic flat-extension property**: every
positive pseudo-moment on the pruned sixth-root basis that annihilates a contact
parent is the point evaluation at its contact.

**Survey of (ND).** 351 C1/C2 contact sources were tested, in bulk, near chart
walls, and along 8-decade gap scans.

- **Rank.** It is `N − 6` at every source, with a singular-value gap above
  10¹³.
- **Strict positivity.** It is certified with margin, at least 10 times the
  solver noise, on 49 of 81 C2 and 24 of 79 C1 sources. No failure was
  confirmed.
- **Walls.** At the walls `g`, `j`, `1 − zh` and `Qc`, the margin scales
  linearly with the gap. The product of the margin and the diverging weight is
  constant to about 3 digits, so this is a normalization effect, not a loss of
  positivity. At the `h` and `Lc` walls the margin stays flat down to gaps of
  10⁻⁸.
- **Unresolved.** Double precision could not decide sources with contacts
  very near commuting measurements (`v_c ≲ 10⁻⁴`), nor the walls `r → 0`,
  `r → 1` and `c → ∞`.

**A route to a complete proof (not carried out).** As a contact approaches
commuting measurements, the rescaled parent tends to a limit with one live
angle. For a single variable, a positive semidefinite matrix polynomial on a
half-line is a square-root Gram form, so these limit regimes should be provable
exactly. A complete proof would combine three pieces:

1. a rigorous interval-arithmetic covering of the compact middle region;
2. exact univariate certificates for the limit regimes;
3. uniform margins that join the two.

This is substantial work, and it is where the paper's 31,000 lines also live.

**The near-commuting end is a singular perturbation.** Along the C2 family with
Bob's contact angle `v_c = a·c → 0`, the native formulas give the following
asymptotics:

- the weights `C` and `D` grow like `a⁻²`, and `A` vanishes like `a²`;
- `n`, `m` and `α` tend to finite limits;
- the contact amplitudes `W_c` and `Y_c` grow like `a⁻¹`.

After rescaling to the contact, the leading order `a⁻⁴` of the parent is a
single square, `(C₀ + D₀)(cb)²(u′W′ − Y′)²`, where `C₀ = lim a²C` and
`D₀ = lim a²D`. It vanishes at the contact. So the limit is not a
nondegenerate one-angle problem: it collapses to rank one, and every direction
transverse to it has to be certified at lower orders in `a`.

Handling the ends of the contact families therefore needs a multi-scale
analysis. The first plan was this: near each end, find the exponents with which
the Gram entries must scale, then check that the rescaled limit problem is
strictly feasible. That is the tropical version of the companion's chart
analysis, and there are many ends: `a → 0`, `r → 0`, `r → 1`, `c → ∞` and
their combinations. The `a → 0` end imposes the following exact constraints on
such a construction.

**The `a → 0` end, analysed.** On C2, two blocks of every PSD Gram certificate
in `Z_k`, at every positive integer k, are fixed by single monomials of `F_c`:

- on the pair `(uZ, uvY)`, `B₁ = [[K, −K], [−K, K + 4A]]`, where
  `K = 4m + 4α + 2d`;
- on the pair `(uvX, vZ)`, `B₂ = [[4D + 2d + 4n, −4D], [−4D, 4D]]`.

Indeed, setting `W=X=Y=0` gives `F/Z²=4α+Ku²+4Dv²`. Its Newton polygon
forces every `Z` feature with `i+j>k` to have zero Gram row: a square containing
such a feature would have an uncancellable highest-weight monomial outside
this polygon. The coefficients of `u²vYZ` and `uv²XZ` then each have a unique
surviving feature pair. This proves the displayed blocks are fixed. The
argument uses `β=0`; with a nonzero `u²v²Z²` coefficient, the same pruning
does not follow.

On C2 their relative determinants are `ρ₁ = 4A/(K + 4A)` and
`ρ₂ = (d + 2n)/(2D + d + 2n)`. At this end `A ~ a²` and `D ~ a⁻²`, so both are
of order `a²` along the stated strict native paths. Each relative determinant
is invariant under diagonal congruence, and its block, normalized to unit
diagonal, has least eigenvalue `1−√(1−ρ_i)`. Relating this to the reduced Gram
margin also requires controlling the projection off the moving contact kernel.
A small principal minor alone does not provide that control: the projector
`I−wwᵀ/(wᵀw)`, with `w=(1,1,ε)`, has nonzero eigenvalues `1,1`, although
its first two-coordinate block has determinant `ε²/(2+ε²)`.

A numerical survey at three interior points P1–P3, with `a = 10⁻¹ … 10⁻⁸`,
reports agreement with `t(a) = 1 − √(1 − ρ₁) ≈ c·a²`, with
`c = 0.198, 0.124, 0.063`. At `a = 0` both rescaled blocks have rank exactly
one. The null vector of `B₁` is `c²b·e(u′Z′) + e(u′v′Y′)`, because
`v_cY_c/Z_c → c²b`.

This is a loss of scale, not of positivity:

- In the sheared basis `(uZ − uvY, uvY)`, `B₁` is `diag(K, 4A)` (when β = 0),
  which is positive for every `a > 0`. `B₂` has the same structure, a large
  rank-one part `4D(uvX − vZ)²` plus an O(1) remainder.
- The contact itself does not degenerate. Its Hessian keeps four positive O(1)
  eigenvalues, one eigenvalue of order `a⁻⁴`, and the exact zero.
- What degenerates is the source weights: the weight `A` of event (−+|00)
  vanishes and the weight `D` of event (−+|11) diverges. In the `d=1`
  normalization there is no finite source at `a=0`. A projective normalization
  reaches a boundary with additional zero coefficients.

`B₁` is one minor of a larger forced structure, with a physical reading. At
`u → ∞`, that is `θ = π`, Alice's two observables are opposite and so jointly
measurable. There the leading coefficient of the C2 parent is

`[u²]F = 4C(1+v²)W² + 4n(1+v²)X² + K(vY − Z)² + 4A(1+v²)Y² + (4D+2d)(W + vX)²`.

At `A = 0` it vanishes on the whole curve `W = X = 0, Z = vY`: the Bell
functional is then tight on a continuum of classical strategies, not only at the
contact. Only pairs of basis elements on the `θ = π` face contribute to `[u²]F`.
So on the span of the moment vectors of that curve, every Gram certificate
equals `4A(1+v²)` along the curve, and by Cauchy–Schwarz it is `O(A)`
throughout the span: evaluate at k+1 fixed distinct positive values of v to
obtain a basis, then bound its mixed Gram entries by their diagonal values.
For k = 6 the span is 7-dimensional (`uv^{j/6}Y`, `j = 0, …, 6`, with
`uZ` tied to `uvY`). Along the stated paths these seven directions have Gram
weight `O(a²)` in a fixed unscaled basis, for every certificate, and
`B₁` is the 2×2 minor at `j = 6`.

**Classical states are nearly tight.** On product states `ψ = W(1, v, u, uv)`,
which are `(1, u) ⊗ (1, v)` up to scale, the parent is exactly

`F = 2W²(1+u²)(1+v²)(2Au² + 2λ(1+v²) + v²(d + 2n + 2β(1+u²)))`.

On C2 this is `2W²D_uD_v(2Au² + (d+2n)v²)`, which is `O(a²)` near `v ~ a`. At
k = 4–6 the optimal dual of the margin problem is (to about 10⁻³) a point
evaluation on this product cone.

**Facial reduction, three independent attempts.** Each attempt was then checked
by a separate reproduction at another point. All three end at "unclear"; no
intrinsic obstruction was found.

- **Complete the square** (with `s = a²`, `P = uW′ − Y′` and `Q = v′(uX − Z)` in
  contact-scaled variables). Exactly,
  `F = λ₀s⁻²(P + sκQ)² + s⁻¹(4C₀v′² + c₁)P² + G₀ + O(s)`.
  - `G₀` has closed-form coefficients in `(c, b, t, r)`. On a grid it is ≥ 0,
    and its only interior zero is the contact.
  - **Forced kernel.** Seven faces of its Newton polytope have real zero
    families, which force 4 contact vectors, `B₁`, `B₂` and 5 Z-edge elements
    `e(i, 6−i, Z)`. Interior-point duals suggested 6–8 more. The complete
    kernel is given below.
  - **Transversality is proved only for `B₁` and `B₂`.** That is, the `a²`
    corrections are shown to be positive on those two directions only.
- **Exact shears.** Three exact shears remove every block fixed by single
  monomials that is `O(a²)`; the smallest remaining minor is 0.62–0.89 at
  P1–P3. The margin still decays, about like `a^{1.9}` over `a = 10⁻²…10⁻³`.
  - **Mechanism.** In the limit form with the stiff square removed, the scaled
    smallest eigenvalue falls like `v′⁻²` across a wedge of large `v′`. Only the
    `a`-dependent terms cut this off, at `v′ ~ 1/a`.
- **Boundary layer.** There is no single intermediate scale. Along `v = a^h`
  the dominant balance is `u ~ a^{−(1−h)}`, a continuum of scales. A frozen
  `a^{j/6}` ansatz shows no plateau. One positive lead: shearing every
  `(6, j, X)`–`(0, j, Z)` pair raises the margin at `a = 10⁻²` by a factor of
  about 2 (reproduced at P1–P3).

**A proposed complete kernel for the limit.** The endpoint calculation proposes
that all forced directions come from the product states
`ψ = w(1, v′, uv′)` with `Y = uW`. The reported restriction of the transformed
limit form is `2w²(1+u²)(p u² + q v′²)`, where `p, q > 0`, with Newton polygon
`hull{(0,2), (2,0), (4,0), (2,2)}`.

The argument takes any Gram certificate `G` and works in three steps:

1. Push `G` forward to the image monomials of the family. The result is
   positive semidefinite and represents that restriction.
2. So any image monomial outside half the polygon gets zero weight.
3. Hence the sum of the basis elements that map to that monomial lies in
   `ker G`.

The reported count is 57 forced directions plus the 4 contact vectors, 61 of
the 126. Verifying this count requires the transformed limit form, the explicit
monomial map and its exact rank. With these 61 directions imposed, numerical
solves report margins 0.16, 0.18 and 0.09 at P1–P3, with
solver status "Solved" and residual 10⁻¹³. The earlier figure of about 10⁻⁴,
reached by removing only 17–19 directions, was an artefact. The directions it
missed had diagonals from 10⁻⁹ to 10⁻¹.

**Proposed grading for transversality.** The candidate scaling of a forced
direction is `s^{k/6} = a^{k/3}`, for k = 1, …, 6. The proposed value of k
comes from the lower convex envelope of the order in `s` over a hexagon, with
class counts 12, 11, 10, 9, 8 and 7. A full finite-`a` construction is still needed.

- **The 7 classes of order `s` are exactly positive.** `B₁` carries `4A₁s`.
  `B₂` and `v′X` carry `s(4n₀+2)`. The five Z-edge elements lie on an edge
  whose polynomial is a positive binomial.
- **The other 50 classes enter at fractional orders,** `a^{1/3}` to `a^{5/3}`.
  This is consistent with the `a^{2/3}` rate seen earlier.
- **No sign obstruction on the product family.** There the form is strictly
  positive at finite `a`.

**Not settled.** Two things are missing:

- the six-level graded basis `a^{k/3}`, checked on the full finite-`a`
  problem;
- the Schur correction for the directions outside the product family.

Double-precision solves at `a ≤ 3·10⁻³` broke down, so this was not tested
numerically.

The resolved numerical cases support (ND); the unresolved cases and the full
native domain remain open. Closedness of `L₆` passes established membership
to finite limits, but supplies no admission proof for this endpoint. The checker
below verifies the corner identities, the projective-face identity and three
rational path asymptotics. The proposed 61-dimensional kernel and graded Schur
completion require further derivation:
[`checks/gram-lift/ends/check_a0_corners.py`](checks/gram-lift/ends/check_a0_corners.py).

So the conjecture now rests on a single nondegeneracy property of the Gram
certificates, (ND), plus the closed faces; a hybrid lift needs only (ND), since
the closed faces can use their exact pencils. (ND) asks that no accidental
degeneracy appears as one moves through each connected family of contact
sources. The four exact seeds and the 302-source sweep are consistent with it:
there, margins grow in proportion to the distance from the boundary.

If it holds, the repo's full-polar assembly (eight relabelled copies) turns it
into an explicit lift of Q222. The PSD order is 8·196 = 1568, plus 72 scalar
sign blocks. Each canonical pencil has 4608 free Gram coefficients; the
eight-copy dual assembly imposes the corresponding adjoint constraints.
A rational strictly feasible point is now explicit:
`c*=(1/21,1,1,1,1,1,1,1,1)`. Its full positive definite Gram matrix is given
in [`checks/gram-lift/strict_feasibility.json`](checks/gram-lift/strict_feasibility.json).
Strict feasibility settles this assembly hypothesis; exactness still needs
the contact admission theorem.

## Admission across the contact families

The open step is universal contact admission, `μσ(c) ≥ 0` for every native
C1/C2 contact source (see [`EXACT_LIFT_ROUTE.md`](EXACT_LIFT_ROUTE.md)). The
ends of the charts, an adversarial search and an exact pointwise test now
give the following picture.

**Exact pointwise certificates.** [`checks/gram-lift/admission/`](checks/gram-lift/admission/)
proves `μσ > 0` at a given rational source in four steps:

1. Solve the equilibrated SDP in double precision.
2. Project the solution exactly, in rationals, onto all 412 coefficient
   equations.
3. Check that no coefficient equation has a nonzero residual.
4. Prove each of the four reduced blocks positive definite by 400-bit ball
   Cholesky.

It certifies the hardest sources found:

- the three worst C1 sources, whose double-precision margins are only
  5·10⁻¹¹ to 10⁻⁹;
- the C2 corner `b ≈ 1 − r`, down to a gap of 0.00125;
- the minimizers of 14 wall scans.

A control source with `n < 0` is rejected, and the dual side then gives an
exact separator.

**No counterexample.** The search ran 69 SDPs and about 1.2·10⁵ surrogate
samples. No source returned `μ < 0`. Every apparent negative fell into one of
two groups:

- a double-precision artefact, with source weights spanning 10¹² to 10³³;
- a dual vector that failed exact repair.

Two spots remain undecided because their weights span more than 10¹¹: C1
sources with `g ≲ 0.012`, and one C2 source near the `g` wall.

**The ends.** Each analysis was reproduced independently where marked.
[`checks/gram-lift/ends/check_c2_ends.py`](checks/gram-lift/ends/check_c2_ends.py)
checks the exact C2 facts.

| End | Exact mechanism | Status |
|---|---|---|
| C2, `r → 0` | `n = n₀/r + O(1)` and the other weights stay finite. Every coefficient containing X is exactly `O(r)`, so one uniform `√r` scale on the X features suffices | closes numerically; margin flat down to `r = 10⁻⁶` (reproduced) |
| C2, `r → 1` | `n(r=1) = −c²(h + a²)/(2(c²(a²+b+1)+1)) < 0`, so the chart ends at an `n = 0` wall first | the `n = 0` wall closes numerically (reproduced); `r → 1` is reachable only through the corner `c² = O(1 − r)`, which stays open |
| C2, `c → ∞` and `c → 0` | Along the torus orbit through the contact, `F_c(O(L)) = 4C(u_cv_cW_c)²(L − 1/L)² + [(4n+4m+2)X_c² + (4A+4C+4D+2)Y_c²](1 − 1/L²)²` for every C2 source. The limit problem lives on 16 charge-0 features, and the other 138 are forced | closes in the interior numerically; margin flat to `c = 10⁴` and `10⁻⁴` (reproduced); the corners stay open |
| C2, `a → 0` | forced kernel from product states, graded `a^{k/3}` (above) | finite-`a` step not verified |
| C1 | On product states `F = 2W²D_uD_v(2Au² + v²(d + 2βD_u))`, which is `O(k²)` just above the contact: an intrinsic near-tightness | most noise-limited sources are numerical; weights spanning ≥ 10¹³ exceed double precision; the `M` and `J` walls are open |

**A covering of the middle region is out of reach with these certificates.**

- **Per-box affine certificates** are rigorous but reach only about 1.3%
  half-width, which means 10⁹ to 10¹⁶ boxes.
- **One linear map from the source data to Gram matrices** is valid on boxes
  of about 10%. The rigorous positivity check on such a box, however,
  certifies only about 0.3–0.4%, which means 10¹¹ to 10¹⁸ boxes.

A proof over the whole continuum therefore needs certificates that are
explicit functions of the source parameters. The paper's charts are of this
kind, and so are the explicit cones of `EXACT_LIFT_ROUTE.md`.

**The frontier.** [`COMPACT_NATIVE_CONTACT_COVER.md`](COMPACT_NATIVE_CONTACT_COVER.md)
covers every compact strict native subfamily by one finite cone. Joining
finitely many closed sound cones is then routine (assembly lemma). What
remains is a finite set of fixed cones, each covering a neighbourhood of one
stratum of the singular boundary of the native domains. That boundary is
large. A first inventory of C1 (scripts not yet in this folder) found:

- **Six genuine facets, proved.** They are `Lc = 0`, `M = 0`, `J = 0`,
  `Qc = 0`, `D = 0` and `β = 0`. Each one contains an exact rational point
  where every other defining inequality is strict.
- **Exact weight behaviour at the walls.** At the walls `A ~ 1/Lc`,
  `m ~ 1/M`, `D ~ 1/J` and `α ~ 1/Qc`, with the contact fixed. The coordinate
  ends send the contact to an edge or to the corner `(0,0)`.
- **About 225 candidate corner strata, numerically.** By number of meeting
  walls they split 12, 47, 77, 61, 23 and 5, at a resolution of 10⁻⁴.

C2 has a comparable structure. Along `r → 0`, two corners keep a flat margin:
with `h → 0` it stays about 0.2, and with `c → 0` about 0.05–0.08. With
`1 − zh → 0` the margin falls linearly in the gap, while `α` diverges. Only
the four C2 ends above have known mechanisms, and none is yet a proof.

## Evidence

The ratio `d_k(e)/d*(e)` compares the largest CHSH weight the lift admits with
the true validity limit, for event vectors `e`. The 302 test vectors are those
of `LIFT_EXPERIMENT.md`: 150 random, 50 π-symmetric, 50 on the C1 face, 50 on
the C2 face, and c̄, c†.

| k | result |
|---|---|
| 1 | Falls short badly: 0.39, 0.12 and 0.34 at the three hardest sources. This is the polynomial obstruction. |
| 2 | Reaches d* on 293 of 302 sources, including c_sep, the source exactly certified to lie outside `C_ram + πC_ram`. It has 9 confirmed gaps (high precision), the worst at 0.652 (C2 face); c̄ and c† give 0.971. |
| 3 | Closes 4 of the 9 gaps. |
| 6 | No source below 0.999·d* among all 302. On the 9 former gaps, the margin at d = (1−δ)·d* is strictly positive for δ = 10⁻², 10⁻³ and 10⁻⁴ and scales in proportion to δ, as expected if d₆ = d*. |

So `L₂` and the repo's ramified cone are incomparable: each certifies a source
the other misses. `L₆` contains `L₂`, and none of the tested sources separates
it from `C_can`.

## Why sixth roots, physically

- **Parity.** Γ = Y⊗Y commutes with every correlator and anticommutes with every
  marginal (`checks/physics/check7_parity.py`). Correlations alone split into
  two decoupled 2×2 sectors, which is Tsirelson's elliptope. The marginals are
  exactly the coupling between the sectors.
- **Soft modes.** The CHSH operator's spectrum is 1 ± √(1 ± sin a sin b). At
  compatible measurements each sector has one soft mode, near ∓½ sin a sin b,
  and the other two modes stay gapped at eigenvalue ≥ 1. All the hard behaviour
  is a two-level problem at the scale sin a·sin b ∝ uv.
- **Scale covariance.** Rescaling the angles multiplies each fractional monomial
  by a constant. A Gram certificate in monomials therefore covers every scale
  regime at once. This is why one SDP can replace the companion's chart-by-chart
  analysis. The half-powers in both angles together supply √u·√v ∝ √(uv), the
  soft-mode amplitude. Ramification in one angle at a time misses it, which is
  why the ramified cone plus its mirror fails.
- **Hurwitz, as motivation only.** The determinant's AM–GM circuits with
  weights (1,1,1) are x⁶ + y⁶ + z⁶ − 3x²y²z² in sixth roots. Hurwitz (1891)
  showed this is a sum of squares, and circuits with weights (1,1,2) need only
  square roots. That suggested the ⅙ grid, and the same substitution turns the
  Motzkin polynomial behind the repo's polynomial obstruction into a Hurwitz form.
  **It is not the mechanism.** All 9 sources that defeat `L₂` need (1,1,1)
  circuits at their contact, but so do 252 of the 293 sources that `L₂`
  certifies (`checks/gram-lift/circuit_types.py`). Neither the contact location
  nor the spread of scales separates the two groups. Whatever distinguishes `L₂`
  from `L₆` lives in the 4×4 matrix structure, not in the scalar determinant.

## What exactness would give

**An explicit lift.** If the conjecture holds, `L₆` is an explicit,
finite-size SDP representation of Q222 with rational data.

**A short proof.** The statement "every valid source has a positive semidefinite
Gram in `Z₆`" is first-order over ℝ, so Tarski transfer carries it to every real
closed field. There, `F_c = (Lᵀ Z₆)ᵀ(Lᵀ Z₆)` with a Cholesky factor `L` over
the source field and pure test-side monomials `Z₆`, so the hard lemma of
`PEN_AND_PAPER.md` follows at once. In fact Scheiderer's criterion is no longer
needed, because the lift is explicit.

**One lemma left.** What remains to be proved is a single statement over the
reals, with no infinitesimals. For the extreme sources (events, the two closed
faces, and the C1/C2 contact sources), `F_c` has a Gram certificate in `Z₆`.

## Caveats and next steps

- **The evidence is numerical.** Clarabel was combined with a δ-margin test. The
  confirmed k = 2 gaps use a 192–320-bit solver. A high-precision solve at
  k = 6 would take about 10 hours per case, so none was run.
- **Fourteen sources have indirect numerical support.** Their k = 6 margins
  are at numerical noise level, while their high-precision k = 2 ratios are
  at least 0.99999999 and `L₂ ⊆ L₆`. This supplies numerical evidence, not
  exact certificates for those sources.
- **Remaining proof.** Closedness, four exact contact certificates, constant
  rank on the native contact charts and a global strictly feasible point
  are established. Persistence of positive definite reduced Grams, (ND),
  remains unproved. The closed faces can instead use their existing pencils.
- **Sampling limits.** 302 samples, including the adversarial ones from earlier
  experiments, are not a proof. A targeted search for sources that minimize
  `d₆/d*` is the next test.
- **Exact certificates.** Exact rational Gram certificates at the 9 former gap
  sources would turn "reaches 0.999" into "contains these sources" exactly.
  That's cheap, as for `checks/lift-experiment/`.

Data and a standalone checker are in
[`checks/gram-lift/`](checks/gram-lift/).
