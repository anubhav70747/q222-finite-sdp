# The finite lift of Q222, from first principles

**Theorem.** The quantum correlation set of two parties with two binary
measurements each, Q222 ⊂ ℝ⁸, is the linear projection of a spectrahedron.

This note gives the proof in the order a physicist would tell it. It starts
from the kinematics of two measurements, then the structure of Bell
inequalities, then positivity at fixed measurements, and finally uniformity
across all measurements, which is where the difficulty lives. Every step is
marked:

- **[here]**: proved in this note, by hand;
- **[paper]**: proved in the paper or its companion, and audited in this review;
- **[checked]**: verified exactly by a script in [`checks/`](checks/);
- **[conjecture]**: supported by evidence, not proved.

Notation: `q = (a₀, a₁, b₀, b₁, E₀₀, E₀₁, E₁₀, E₁₁)`; `u = tan(a/2)`,
`v = tan(b/2)`; `D_u = 1+u²`, `D_v = 1+v²`.

---

## I. Kinematics of two measurements

**1. Real qubits [paper, Ch. 2].** Two ±1 observables generate an algebra whose
centre contains `½{A₀,A₁}` and `[A₀,A₁]²`. Jordan's lemma diagonalizes this
centre. Every block has dimension at most two, and on it

    A₀ = Z,   A₁ = cos a·Z + sin a·X,   ½{A₀,A₁} = cos a,   [A₀,A₁]² = −4 sin²a.

So a measurement pair is characterized by one angle, and its **incompatibility**
is `sin a`. The observables are real, so complex states mix real ones.
Carathéodory's theorem in ℝ⁸ then gives `Q222 = conv Φ(S¹×S¹×S³)`, with at most
nine atoms.

**2. Parity and the two-mode picture [checked: `physics/check7_parity.py`,
`physics/check8_two_mode.py`].** Let `Γ = Y⊗Y`. It commutes with every
correlator `A_x⊗B_y` and anticommutes with every marginal `A_x`, `B_y`. In the
Bell basis it splits ℝ⁴ into two planes: `{Φ⁺, Ψ⁻}` (Γ = +1) and `{Φ⁻, Ψ⁺}`
(Γ = −1). Write the state as two complex amplitudes,
`z₊ = ⟨Φ⁺|ψ⟩ + i⟨Ψ⁻|ψ⟩` and `z₋ = ⟨Φ⁻|ψ⟩ + i⟨Ψ⁺|ψ⟩`. Then every Bell
operator `f₀ + Σ f_xy A_xB_y + Σ f_x A_x + Σ g_y B_y` has the form

    f₀(|z₊|² + |z₋|²) + Re(W₊ z̄₊²) + Re(W₋ z̄₋²) + 2 Re(Ω z₊ z̄₋) + 2 Re(N z₊ z₋),

with

    W₊ = Σ f_xy e^{i(β_y − α_x)},   W₋ = Σ f_xy e^{i(α_x + β_y)},
    Ω  = Σ f_x e^{iα_x},           N  = Σ g_y e^{−iβ_y}.

Read physically, this is a two-mode quadratic Hamiltonian:

| term | origin | role |
|---|---|---|
| W₊ | correlators, angle differences | single-mode squeezing |
| W₋ | correlators, angle sums | single-mode squeezing |
| Ω | Alice's marginals | beam-splitter coupling |
| N | Bob's marginals | two-mode squeezing |

**Corollary (Tsirelson–Landau–Masanes) [here].** Each mode alone has energies
`f₀ ± |W|`. A correlation-only Bell expression is therefore valid exactly when
`f₀ ≥ |W₊|` and `f₀ ≥ |W₋|` for all angles. For CHSH,

    |W| ≤ 2|cos((α₀−α₁)/2)| + 2|sin((α₀−α₁)/2)| ≤ 2√2.

**All of the difficulty of Q222 comes from the marginals**, which are the only
couplings between the two modes.

## II. Bell inequalities

**3. CHSH plus penalties [here; paper, Ch. 4].** Let `H = 1 − S/2` with
`S = E₀₀ − E₀₁ + E₁₀ + E₁₁`. The PR box and its eight deterministic neighbours
are affinely independent. So every slack that is negative at the PR box is
uniquely

    f = d·H + 4 Σᵢ cᵢ pᵢ,   d = −f(PR) > 0,   cᵢ = f(vᵢ)/4 ≥ 0,

where `pᵢ` is the probability of the event lost by neighbour `vᵢ`. Each event
`(a,b|x,y)` is the projector `(1+aA_x)(1+bB_y)/4`, a rank-one product state, and
it contributes to every term of the two-mode form. A slack that is nonnegative
at every PR box is a sum of event probabilities, which are already squares.
Up to symmetry, the extreme sources are:

- events;
- two closed faces with small explicit lifts;
- contact sources of type **C1** (λ = n = 0) or **C2** (λ = β = 0) [paper].

**4. One violating direction [here; checked: `physics/check6_tsirelson.py`].**
Landau's identity gives `Ŝ² = 4 + [A₀,A₁]⊗[B₀,B₁] = 4 + 4 sin a sin b·Γ`. The
spectrum of `I − Ŝ/2` is

    1 ± √(1 + sin a sin b)  (one sector),    1 ± √(1 − sin a sin b)  (the other).

- **Exactly one eigenvalue is negative.** It is `1 − √(1 + |sin a sin b|)`,
  about `−½|sin a|·|sin b|`.
- **Penalties cannot create more.** They are positive semidefinite, so the
  Bell operator keeps at most one negative direction.
- **Hence positivity is one scalar inequality**, `det ≥ 0`, because the
  three-dimensional transverse block is positive definite [paper].

At compatible measurements each sector has one **soft mode**, at about
`∓½ sin a sin b`, while the two **hard modes** stay at eigenvalue ≥ 1. The
only hard regime is the classical limit. There, the violation lives in one soft
mode, the marginals couple it to the other, and validity is the stability of
this two-level system. The violation depth is a product of two
incompatibilities, and each can vanish at its own rate. That is the physical
origin of the two independent scales in step 7.

## III. Positivity at fixed measurements

**5. Effective resistance [here; checked: `physics/check4_thomson.py`].** Write
the parent form, `F = D_uD_v·⟨ψ|f̂|ψ⟩`, as `F = ψᵀTψ − (ℓ·ψ)²`. Here
`T = Σₖ aₖ gₖgₖᵀ` collects the positive squares, and `ℓ = √(2dD_uD_v)·|11⟩` is
the violating direction. For any flow `Σₖ tₖgₖ = ℓ`:

    F = Σₖ aₖ (gₖ·ψ − (tₖ/aₖ)·s)² + s²·(1 − Σₖ tₖ²/aₖ),   s = ℓ·ψ.

The minimal cost `Σ tₖ²/aₖ` is `ℓᵀT⁻¹ℓ` (Thomson's principle). So **a Bell
inequality holds at given measurements exactly when the violating direction has
effective resistance ≤ 1 through the network of losing events**, whose
conductances are the Bell coefficients.

**6. The determinant, and AM–GM as force balance [here; checked:
`determinant/`, `circuits/`].**

    det P = 16 D_u² D_v² u² v² · L(u², v²),
    L(U,V) = h + Σ_{e ∈ {−1,0,1}²∖0} q_e U^{e₁} V^{e₂},

where all eight outer `q_e ≥ 0` and only the centre `h` can be negative.

- **Balance.** At a stationary contact `(U_c, V_c)`, the weights
  `w_e = q_e U_c^{e₁} V_c^{e₂}` satisfy `h = −Σ w_e` and `Σ w_e e = 0`. This is
  a force balance in the plane of log half-angles: each losing event pushes the
  optimal measurement one way, and self-testing is the unique equilibrium.
- **Circuits.** The balanced weights split into the 12 minimal circuits of the
  3×3 square: opposite pairs, and triangles with weights (1,1,1) or (1,1,2).
- **Squares.** Each circuit becomes squares through `p + 1/p − 2`, the cubic
  AM–GM identity and the quartic one.
- **Hurwitz (1891).** x⁶ + y⁶ + z⁶ − 3x²y²z² is a plain sum of squares. So
  every circuit is a sum of squares of fractional monomials in `u^{1/6}` and
  `v^{1/6}`. The (1,1,1) circuits need sixth roots of `U = u²`, which is
  `u^{1/3}`. The (1,1,2) circuits need fourth roots, which is `u^{1/2}`.

## IV. Uniformity: the heart of the proof

**7. What has to hold [paper; Scheiderer, Cor. 3.19].** A finite lift exists
exactly when the following holds. For every valid slack f with coefficients in a
real closed field `K_c`, and every behaviour q in a real closed field `K_t`,
`[f,q]` must be a sum of squares in `K_c ⊗ K_t`. Infinitesimals in the two
fields describe the inequality and the experiment degenerating at unrelated
rates. Roots and inverses are therefore allowed only inside one field.

Fractional monomials are the natural currency here. Rescaling the angles
multiplies `u^{i/6}v^{j/6}` by a constant, so a certificate written in monomials
holds at every scale at once. The circuit certificates of step 6 are of this
kind.

**8. A uniform identity [here; checked: `physics/check12_decomp_sos.py`,
`physics/check3_cauchy_binet.py`].**

    det(T) · F = (Pψ)ᵀ adj(T) (Pψ) + det(P) · s²,   adj(T) = Σ_{|S|=3} a_S w_S w_Sᵀ.

The first term is a sum of squares, with pure source weights and pure test-side
vectors. The second is a sum of squares by step 6. So `det(T)·F` is certified
for every contact source with no case analysis. The only obstruction left is the
multiplier `det(T) = Σ_S a_S det[g_S]²`, which mixes the two fields.

**9. The hard lemma 𝒞 [paper; audited].** For every C1 or C2 contact source and
every test, `F` itself is a sum of squares in `K_c ⊗ K_t`.

- **How the companion proves it.** It partitions by scale regime, about 31,000
  lines. On each regime it rescales by pure quantities until the residual Gram
  matrix has a strictly positive real part, keeps unbounded penalties as whole
  squares, and joins regimes with division-free conductor identities.
- **Audit.** Every sketched step was re-derived in this review and holds; the
  defects are editorial.
- **The split is necessary [checked:
  `physics/check_c2_to_c1_counterexample.py`].** A C2 source's formal C1 image
  can fail positivity, so C1 and C2 cannot be merged.
- **Physically,** the regimes are the scales of the two-level soft-mode problem
  of step 4.

**10. A chart-free route [conjecture; checked numerically: `gram-lift/`].** The
sextic Puiseux–Gram lift asks for `F_c = Zᵀ G Z` with `G ⪰ 0` in the fixed
basis `Z = {u^{i/6} v^{j/6} ψ_l}`. Scale covariance (step 7), Hurwitz's identity
(step 6) and the soft-mode amplitude `√(uv)` (step 4) all point to this basis.
On 302 test sources it never falls below 0.999 of the validity limit, while
half-roots fail on 9. See [`EXACT_LIFT_CANDIDATE.md`](EXACT_LIFT_CANDIDATE.md).

If it is exact, "every valid source has such a G" is a statement over ℝ. Tarski
transfer and a Cholesky factor over `K_c` then give 𝒞 at once, and the lift is
explicit.

## V. Assembly [here]

A valid slack is a nonnegative combination of at most nine extreme sources, and
a behaviour is a mixture of at most nine qubit atoms. If each pair is certified,
`[Σ aⱼfⱼ, Σ wᵢqᵢ] = Σ ((√aⱼ ⊗ √wᵢ)·hᵢⱼₖ)²`, where `Σ wᵢ = 1` preserves the
constant term. Scheiderer's criterion then gives a finite lift, and Tarski
transfer gives one with real algebraic entries. ∎

## VI. The explicit lift

| object | status |
|---|---|
| Tilted-CHSH projection, `[[√2+t, s/2],[s/2, √2−t]] ⊕ (1−t) ⊕ (1+t) ⪰ 0` | exact [paper; checked] |
| Canonical reduction: valid ⇔ ∀v ∃J with two 4×4 matrices PSD | exact [paper; checked] |
| Full-polar assembly: an exact canonical pencil of order n gives Q222 at order 8n | exact [paper; checked] |
| Fixed-function cones: polynomial, half-line (33), ramified (49) | strictly too small [paper] |
| Ramified cone plus its mirror (98; 784 assembled) | strictly too small, by an exact rational certificate [checked: `lift-experiment/`] |
| Sextic Puiseux–Gram lift `L₆` (Gram 196 per copy) | candidate exact lift [conjecture] |

## The whole proof in one paragraph

Two measurements per party reduce to one angle each, and a Bell operator becomes
a two-mode quadratic Hamiltonian. Correlators squeeze each mode; marginals
couple the two. The CHSH part has a single violating direction whose depth is
the product of the two parties' incompatibilities. Event penalties stabilize
it, with Bell coefficients acting as conductances. At the optimal measurements,
the penalties balance like forces, and AM–GM, in Hurwitz's sum-of-squares form,
certifies positivity. What remains is to keep these certificates valid as the
inequality and the experiment approach the classical limit at unrelated rates.
The paper does this regime by regime. A single fixed basis of sixth-root
monomials, which are blind to scale, appears to do it in one stroke.
