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
family, and equals `N − 6` by the exact seeds.

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

So the conjecture now rests on a single nondegeneracy property of the Gram
certificates, (ND), plus the closed faces. (ND) asks that no accidental
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
