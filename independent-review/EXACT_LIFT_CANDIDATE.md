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
| 6 | sixth-root grid | 64, 48, 48, 36 | 2763 |

**Conjecture.** `L₆ = C_can`.

If it holds, the repo's full-polar assembly (eight relabelled copies) turns it
into an explicit lift of Q222. The PSD order is 8·196 = 1568, plus 72 scalar
sign blocks and about 22,000 auxiliary variables. The assembly needs a strictly
feasible point, which is still to be checked.

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
- **Hurwitz.** The determinant's AM–GM circuits with weights (1,1,1) are
  x⁶ + y⁶ + z⁶ − 3x²y²z² in sixth roots. Hurwitz (1891) showed this is a sum of
  squares. Circuits with weights (1,1,2) need only square roots. The common grid
  is therefore ⅙, and k = 2 is expected to fail on sources whose certificate
  needs a (1,1,1) circuit, which is what the sweep shows. The same substitution
  turns the Motzkin polynomial behind the repo's polynomial obstruction into a
  Hurwitz form.

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
- **Fourteen sources were settled indirectly.** They have noise-level k = 6
  margins and are settled by inclusion: their k = 2 ratios are at least
  0.99999999 by high precision, and `L₂ ⊆ L₆`.
- **Untested properties.** Closedness of `L₆`, strict feasibility for the
  assembly, and the boundary sources themselves are not yet verified.
- **Sampling limits.** 302 samples, including the adversarial ones from earlier
  experiments, are not a proof. A targeted search for sources that minimize
  `d₆/d*` is the next test.
- **Exact certificates.** Exact rational Gram certificates at the 9 former gap
  sources would turn "reaches 0.999" into "contains these sources" exactly.
  That's cheap, as for `checks/lift-experiment/`.

Data and a standalone checker are in
[`checks/gram-lift/`](checks/gram-lift/).
