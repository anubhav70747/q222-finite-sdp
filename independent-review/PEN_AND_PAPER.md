# The finite lift of Q222: a proof you can follow

**Theorem** (the paper's main result). The set Q222 ⊂ ℝ⁸ of quantum
correlations of two parties, each with two binary measurements, is the linear
projection of a spectrahedron. The matrices can be chosen with real algebraic
entries.

This note proves the theorem from five ingredients, [I1]–[I5].

- **[I1]** is a published theorem of Scheiderer, used as stated.
- **[I2]–[I4]** come from the paper, and their short proofs are given below.
- **[I5]** is the hard lemma. It is proved in the paper's companion by a
  case analysis of about 31,000 lines, which this review audited. Its structure
  is outlined below, and Part B explains what it does. A short proof of [I5]
  is the goal of the chart-free candidate route in Part C.

Everything else is proved here in a few lines.

**The idea, in three sentences.** A finite lift exists exactly when every Bell
inequality's positivity can be certified by sums of squares that stay valid
however the inequality and the experiment degenerate, independently of each
other. Both sides reduce to finitely many atoms: qubit experiments, and a few
extreme Bell inequalities, each of which is CHSH plus penalties on the eight
outcomes the PR box never produces. What remains is one lemma: such
certificates exist for the extreme inequalities that touch the quantum set at a
point. Its physics is a single violating direction stabilized by those
penalties, and it carries all of the paper's long appendix.

Part A is the proof. Part B explains why it works and what the hard lemma is
about; none of Part B is used in Part A. Part C records related results that are
not used.

Tags: **[I1]–[I5]** imported, stated below; **[checked: …]** verified exactly by
a script in [`checks/`](checks/).

## Notation

| symbol | meaning |
|---|---|
| `q = (a₀, a₁, b₀, b₁, E₀₀, E₀₁, E₁₀, E₁₁)` | a behaviour: marginals and correlators |
| `θ`, `φ` | Alice's and Bob's measurement angles; `u = tan(θ/2)`, `v = tan(φ/2)` |
| `D_u = 1+u²`, `D_v = 1+v²` | denominators of the half-angle parametrization |
| `σ_z`, `σ_x`, `ε = [[0,−1],[1,0]]` | real Pauli matrices; `ε` is the real form of σ_y |
| `ψ = (W, X, Y, Z)` | a real two-qubit state, amplitudes in the basis `|00⟩, |01⟩, |10⟩, |11⟩` |
| `Φ^± = (|00⟩ ± |11⟩)/√2`, `Ψ^± = (|01⟩ ± |10⟩)/√2` | the Bell basis |
| `H = 1 − S/2`, `S = E₀₀ − E₀₁ + E₁₀ + E₁₁` | the CHSH slack |
| PR box | `E₀₀ = E₁₀ = E₁₁ = 1`, `E₀₁ = −1`, zero marginals; there `S = 4` and `H = −1` |
| `c = (d, λ, n, m, C, A, D, α, β) ≥ 0` | a canonical source: CHSH weight `d` and eight event weights |
| `K_c`, `K_t` | real closed fields containing ℝ, carrying the source and the test |
| `𝕋 = K_c ⊗_ℝ K_t` | the ring where certificates live; `Σ𝕋²` is its sums of squares |
| `[f, q] = f₀⊗1 + Σ_ν f_ν⊗q_ν` | the separated evaluation of an affine function `f` at `q` |

The eight event weights name the outcomes `(ab|xy)` that the PR box never
produces. These are the outcomes with `ab = −E_xy` at the PR box. Each of the
eight deterministic boxes adjacent to the PR box produces exactly one of them.

| weight | λ | n | m | C | A | D | α | β |
|---|---|---|---|---|---|---|---|---|
| event `(ab\|xy)` | (++\|01) | (+−\|00) | (+−\|11) | (−+\|10) | (−+\|00) | (−+\|11) | (−−\|01) | (+−\|10) |

## Imported results

- **[I1] Scheiderer's criterion** (SIAM J. Appl. Algebra Geom. 5 (2021),
  Cor. 3.19). A closed convex semialgebraic set `K ⊂ ℝⁿ` is a spectrahedral
  shadow if and only if, for every real closed field `R ⊇ ℝ`, every affine `f`
  nonnegative on `K_R` and every `q ∈ K_R`, `[f, q]` is a sum of squares in
  `R ⊗_ℝ R`.
- **[I2] The model** (paper, Ch. 2). Q222 is the convex hull of the behaviours of
  real two-qubit atoms, with `A₀ = B₀ = σ_z`, `A₁ = cos θ·σ_z + sin θ·σ_x` and
  `B₁ = cos φ·σ_z + sin φ·σ_x`. Every behaviour is a mixture of at most nine
  atoms, by Jordan's lemma, the real/imaginary split and Carathéodory. This is a
  first-order formula over ℚ, so it holds verbatim over every real closed field.

  *Proof.*
  1. **Dilate.** Dilate each party's two binary effects to projections on one
     common space. For an effect E, use
     `[[E, √(E(I−E))], [√(E(I−E)), I−E]]` on `ℋ⊕ℋ`, with the isometry
     `h ↦ (h, 0)`.
  2. **Split into blocks (Jordan).** Two ±1 observables A₀, A₁ split the space
     into invariant blocks of dimension at most two. The unitary `A₀A₁` is
     conjugated to its inverse by A₀, so its eigenvectors pair up. On a
     two-dimensional block, a real change of basis gives `A₀ = σ_z` and
     `A₁ = cos θ·σ_z + sin θ·σ_x`. One-dimensional blocks are commuting ±1
     values.
  3. **Reduce to atoms.** Expectations depend only on the diagonal blocks of
     the state, so every behaviour is a mixture of block atoms. Because the
     observables are real, a complex vector `r + is` contributes the real
     states `r` and `s`: the cross terms cancel.
  4. **Bound the number of atoms.** Carathéodory in ℝ⁸ gives at most nine
     atoms.
  5. **Converse.** Every atom is a quantum behaviour. ∎
- **[I3] Extreme sources** (paper, Ch. 4). A slack that is nonnegative at every
  PR box is a nonnegative combination of event probabilities. Otherwise, after
  relabelling outcomes, which preserves Q222, it is negative at the PR box above.
  It then has the unique form `f = d·H + 4Σ c_i p_i`, and is a nonnegative
  combination of at most nine sources of these kinds:
  - events (single outcome probabilities);
  - extreme sources of the two noncoercive closed faces, `λ = m = 0` and
    `λ = α = 0`;
  - contact sources of type C1: `λ = n = 0`; `d, m, C, A, α > 0`; `D, β ≥ 0`;
    `D + β > 0`;
  - contact sources of type C2: `λ = β = 0`; `d, m, C, A, D, α > 0`; `n ≥ 0`.

  A *contact* is a finite stationary zero with `uv ≠ 0`:
  `P(u_c,v_c)ψ_c = 0` and `ψ_cᵀ ∂_u P ψ_c = ψ_cᵀ ∂_v P ψ_c = 0`, where `P` is
  the state matrix of Insight 1. The classification is first-order, so it holds
  over every real closed field.

  *Proof (paper, Ch. 4).*
  1. **Slacks nonnegative at every PR box.** The no-signalling polytope has 16
     deterministic and 8 PR vertices. If `f ≥ 0` at all of them, then `f ≥ 0`
     on the polytope, and Farkas's lemma writes it as a nonnegative combination
     of outcome probabilities.
  2. **At most one negative PR box.** Any two PR boxes have a local midpoint,
     where `f ≥ 0`. So `f` is negative at no more than one of them, and a
     relabelling makes it the one above.
  3. **The canonical form.** That PR box and its eight deterministic neighbours
     `v_i` are affinely independent and span the no-signalling space. So `f` is
     determined by its nine values there, which gives `d = −f(PR) > 0` and
     `c_i = f(v_i)/4 ≥ 0`.
  4. **Extreme rays.** The canonical sources form a closed cone in ℝ⁹ with a
     compact base, so every source is a combination of at most nine extreme
     rays. An extreme ray with `d = 0` is an event.
  5. **Coercive supports.** For `d > 0`, call the support *coercive* when the
     outer terms of the determinant (Insight 5) grow in every direction of
     log-angle space. An extreme source with coercive support must attain zero;
     otherwise opposite small perturbations stay valid. That zero is a contact.
     It imposes six homogeneous linear equations on the nine coefficients, so
     extremality leaves at most six positive event weights. The omitted events
     must be pairwise compatible on the eight-cycle of events, which leaves at
     least five. Up to the dihedral symmetries of the cycle, the possible
     patterns are exactly C1 and C2.
  6. **Noncoercive supports.** These form the two closed faces. ∎
- **[I4] Closed faces** (paper, Ch. 4, SG11–SG14). On the face `λ = m = 0` the
  parent is a Gram form built from one 3×3 matrix polynomial in `v` (SG12). On
  the face `λ = α = 0` it uses two 3×3 blocks with one free parameter (SG13).
  Both are strictly feasible, and conic duality makes every extreme source of
  these faces a separated sum of squares (SG14).

  *Proof idea.* On each face the parent depends on one angle through a matrix
  polynomial of low degree. A univariate positive semidefinite matrix
  polynomial is a Gram form, so validity is an explicit LMI with a strictly
  feasible point. Duality then writes the parent as `Σ X_ij·Y_ij`, with `X ⪰ 0`
  over `K_c` and `Y ⪰ 0` over `K_t`, and such a sum is a separated sum of
  squares: factor X over K_c and Y over K_t. ∎
- **[I5] The hard lemma 𝒞** (paper, Proof.tex Prop. 4 and the C1/C2 appendices;
  audited in this review). Let `K_c` and `K_t` be any real closed fields
  containing ℝ. For a C1 or C2 source over `K_c` with a contact, the parent
  form `F_c` (Insight 1) is a sum of squares in `𝕋`, at every projective
  testing pair and every real state over `K_t`.

---

## Part A. The proof

**Step 1. What suffices.** By [I1] it is enough to show `[f, q] ∈ Σ𝕋²` for all
real closed `K_c, K_t ⊇ ℝ`, every valid `f` over `K_c` and every `q ∈ Q222`
over `K_t`. Taking `K_c = K_t = R` (two labelled copies) gives [I1]'s
condition.

Certificates must be identities in the ring `𝕋`, which is not a field. Roots
and inverses are therefore available only for positive elements of one factor;
for example, `1⊗1 + x⊗y` with `x`, `y` infinitesimal has no inverse in `𝕋`.

**Step 2. Reduce to atoms.** By [I3] over `K_c`, after a relabelling that maps
sums of squares to sums of squares, write `f = Σ_j r_j f_j` with `r_j ≥ 0` and
extreme `f_j`. By [I2] over `K_t`, write `q = Σ_i w_i q_i` with `w_i ≥ 0`,
`Σ w_i = 1` and qubit atoms `q_i`. Bilinearity gives

    [f, q] = Σ_{i,j} (r_j ⊗ w_i)·[f_j, q_i].

The constant term matches because `Σ w_i = 1`. If `[f_j, q_i] = Σ_k σ_{ijk}²`
with `σ_{ijk} ∈ 𝕋`, then `[f, q] = Σ_{i,j,k} ((√r_j ⊗ √w_i)·σ_{ijk})²`, since
`√r_j ∈ K_c` and `√w_i ∈ K_t`.

**Step 3. Atoms reduce to parent forms.** For a qubit atom with half-angles
`u, v` and state `ψ`, `F_c(u,v,ψ) = D_uD_v‖ψ‖² · [f_c, q(u,v,ψ)]`. The
prefactor `D_uD_v‖ψ‖²` is a positive element of `K_t`, so `[f_c, q] ∈ Σ𝕋²`
exactly when `F_c ∈ Σ𝕋²`. Angles with `u` or `v ∈ {0, ∞}` are covered by the
homogeneous form of the parent.

**Step 4. Each extreme source.**
- **Events.** An event term is the square of a linear form times a positive
  element of `K_t`.
- **Closed faces.** These are [I4].
- **Contact sources.** These are [I5], which covers every test.

**Step 5. Conclusion.** By Steps 2–4, `[f, q] ∈ Σ𝕋²` in every case. By Step 1
and [I1], Q222 is a spectrahedral shadow. For fixed sizes, "this LMI projects
onto Q222" is a first-order sentence over ℚ in the matrix entries. A real
solution therefore implies a real algebraic one (Tarski). ∎

---

## Part B. Why it works

Nothing below is used in Part A. It explains the objects in [I3]–[I5], and why
[I5] is the hard part.

**Insight 1. The Bell operator [checked: `physics/check6_tsirelson.py`,
`physics/check12_decomp_sos.py`].** For a canonical source, the parent form
(paper, Proof.tex eq. insight-parent) is

    F_c = 4λD_u(W+vX)² + 4nD_uD_vX² + 4m(vW−X+uvY−uZ)² + 4CD_v(uW−Y)²
        + 4AD_uD_vY² + (4D+2d)(uW+uvX−Y−vZ)² + (4α+2d)D_u(vY−Z)²
        + (4β+2d)D_v(X+uZ)² − 2d·D_uD_v·Z².

It equals `D_uD_v·⟨ψ| d(I − Ŝ/2) + 4 Σ c_i Π_i |ψ⟩`, where `Ŝ` is the CHSH
operator for `S`, and each `Π_i` projects onto the product state of event `i`,
so it has rank one. Write `F_c = ψᵀPψ = ψᵀM₊ψ − s²`:

- `M₊ = Σ_k κ_k g_k g_kᵀ` collects the eight positive squares, with
  conductances `κ = (4λ, 4n, 4m, 4C, 4A, 4D+2d, 4α+2d, 4β+2d)`;
- `ℓ = √(2d·D_uD_v)·(0,0,0,1)` is the violating direction, and `s = ℓ·ψ`.

**Insight 2. Two coupled modes [checked: `physics/check7_parity.py`,
`physics/check8_two_mode.py`].**
- **Parity.** Let `Γ = ε⊗ε`. Then Γ commutes with every correlator `A_x⊗B_y`
  and anticommutes with every marginal `A_x`, `B_y`.
- **Two modes.** Let θ_x and φ_y be the setting angles, so `θ₀ = φ₀ = 0`,
  `θ₁ = θ` and `φ₁ = φ`. Set `z₊ = ⟨Φ⁺|ψ⟩ + i⟨Ψ⁻|ψ⟩` and
  `z₋ = ⟨Φ⁻|ψ⟩ + i⟨Ψ⁺|ψ⟩`; for real ψ these simply repackage the four
  amplitudes. Then every Bell operator
  `f₀ + Σ f_xy A_xB_y + Σ f_x A_x + Σ g_y B_y` reads

      f₀(|z₊|² + |z₋|²) + Re(W₊ z̄₊²) + Re(W₋ z̄₋²) + 2Re(Ω z₊ z̄₋) + 2Re(N z₊ z₋),

  with `W₊ = Σ f_xy e^{i(φ_y−θ_x)}`, `W₋ = Σ f_xy e^{i(θ_x+φ_y)}`,
  `Ω = Σ f_x e^{iθ_x}` and `N = Σ g_y e^{−iφ_y}`.
- **Physical reading.** Correlators squeeze each mode, through angle
  differences in one and angle sums in the other. Alice's marginals act as a
  beam splitter, and Bob's as two-mode squeezing.
- **Corollary (Tsirelson–Landau–Masanes).** Each mode alone has energies
  `f₀ ± |W|`. A correlation-only expression is therefore valid exactly when
  `f₀ ≥ |W₊|` and `f₀ ≥ |W₋|` at all angles. For CHSH,
  `|W| ≤ 2|cos(θ/2)| + 2|sin(θ/2)| ≤ 2√2`.
- **The point.** The marginals are the only coupling between the modes, and
  all of Q222's difficulty comes from them.

**Insight 3. One violating direction [checked: `physics/check7_parity.py`].**
Landau's identity gives `Ŝ² = 4I + [A₀,A₁]⊗[B₀,B₁] = 4I + 4 sin θ sin φ·Γ`.
So `I − Ŝ/2` has eigenvalues `1 ± √(1 + sin θ sin φ)` in one parity sector and
`1 ± √(1 − sin θ sin φ)` in the other.

- **The violation.** At most one eigenvalue is negative, about
  `−½|sin θ|·|sin φ|`. The penalties are positive semidefinite, so `P` keeps at
  most one negative direction.
- **Why det ≥ 0 suffices.** Positivity of `P` is equivalent to `det P ≥ 0`
  only because, for `d > 0` and `uv ≠ 0`, the `(W,X,Y)` block of `P` is
  positive definite (paper, SG3). At `sin θ sin φ = 0`, the measurements of one
  party commute and the CHSH slack is directly nonnegative.
- **The classical limit.** At compatible measurements each sector has one soft
  mode, near `∓½ sin θ sin φ`, while the two hard modes stay at eigenvalue ≥ 1.
  Near this limit, validity is the stability of two soft modes coupled by the
  marginals.

**Insight 4. Effective resistance [checked: `physics/check4_thomson.py`].** For
any currents `τ` with `Σ_k τ_k g_k = ℓ`,

    F = Σ_k κ_k (g_k·ψ − (τ_k/κ_k)·s)² + s²·(1 − Σ_k τ_k²/κ_k).

The least dissipation `Σ τ_k²/κ_k` is `ℓᵀM₊⁻¹ℓ` (Thomson's principle, for
`M₊ ≻ 0`). So, at fixed measurements, a Bell inequality holds exactly when the
violating direction has effective resistance at most 1 through the network of
penalized events.

**Insight 5. The determinant is an AM–GM object [checked: `determinant/`,
`circuits/`].**

    det P = 16 D_u² D_v² u² v² · L(u², v²),
    L(U, V) = η₀ + Σ_{e ∈ {−1,0,1}²∖0} η_e U^{e₁} V^{e₂}.

Each outer `η_e` is a product of event weights and positive linear forms, for
example `η₍₁,₁₎ = 4Cβ(2A+2α+d+2m)(2D+d+2λ+2n)`. Only the centre `η₀` can be
negative (paper, SG5–SG6).

- **Contact gives a balanced minimum.** At a contact, the `(W,X,Y)` block is
  positive definite, so `adj P ∝ ψ_cψ_cᵀ`. Jacobi's formula then gives
  `L = ∂_U L = ∂_V L = 0`. With `μ_e = η_e U_c^{e₁} V_c^{e₂}` this reads
  `η₀ = −Σ μ_e` and `Σ μ_e·e = 0`.
- **Convexity makes it global.** `Σ η_e e^{e·x}` is convex in
  `x = (log U, log V)`, so the stationary point is the global minimum. Hence
  `L ≥ 0` everywhere: a single contact certifies validity at every angle. The
  balance `Σ μ_e e = 0` is a force balance in log-angle space.
- **Circuits.** A balanced nonnegative vector peels, greedily, into the 12
  minimal circuits of the 3×3 square: opposite pairs, and triangles with
  weights (1,1,1) or (1,1,2).
- **Squares.** Each circuit is a sum of squares through

      p + 1/p − 2 = (√p − 1/√p)²,
      a³ + b³ + c³ − 3abc = ½(a+b+c)·[(a−b)² + (b−c)² + (c−a)²],
      a⁴ + b⁴ + c⁴ + d⁴ − 4abcd = (a²−b²)² + (c²−d²)² + 2(ab−cd)²,

  applied to the monomials `U_c⁻¹⊗U` and `V_c⁻¹⊗V`. These are positive pure
  tensors, so their roots are taken separately in each field.
- **Conclusion.** `L`, and hence `det P`, is a separated sum of squares at
  every contact source.

**Insight 6. Everything except one multiplier [checked:
`physics/check12_decomp_sos.py`, `physics/check3_cauchy_binet.py`].** From
`P = M₊ − ℓℓᵀ` and the matrix determinant lemma,
`det P = det M₊·(1 − ℓᵀM₊⁻¹ℓ)`. Expanding `(Pψ)ᵀM₊⁻¹(Pψ)` then gives

    det(M₊)·F = (Pψ)ᵀ adj(M₊) (Pψ) + det(P)·s².

By Cauchy–Binet, `adj(M₊) = Σ_{|Λ|=3} κ_Λ w_Λ w_Λᵀ`. Here `Λ` runs over 3-sets
of events, `κ_Λ` is the product of their conductances (a positive element of
`K_c`), and `w_Λ` is the vector of 3×3 minors of their event vectors, an
element of `K_t⁴`. So the first term is a sum of squares in `𝕋`, and so is the
second, by Insight 5.

**`det(M₊)·F` is therefore certified with no case analysis.** [I5] is exactly
the removal of the multiplier `det M₊ = Σ_{|Λ|=4} κ_Λ det[g_Λ]²`, which mixes
the two fields and has no inverse in `𝕋`.

**Insight 7. Why [I5] is long.** Two independent degenerations meet at the
multiplier.

- **The source and the test** approach the boundary at unrelated rates. This
  is why the fields are separate.
- **Inside the test,** Alice's and Bob's incompatibilities, `sin θ` and
  `sin φ`, vanish at unrelated rates. This is the depth `½ sin θ sin φ` of
  Insight 3.

The companion proves [I5] regime by regime, about 31,000 lines. On each regime
it rescales by pure quantities until the remaining Gram matrix has a strictly
positive real part. It then applies its residue lemma (paper, Ch. 3, TF1–TF2):
a bounded Gram matrix in `𝕋` whose real residue is positive definite is a sum
of squares. This review audited the regimes and re-derived the steps that are
only sketched; the defects found are editorial.

One natural shortcut fails. The paper's map OP4 (`c2_short.tex:8796–8825`)
sends a C2 source to a formal C1 source at the same contact, but the resulting
parent can be nonpositive [checked: `physics/check_c2_to_c1_counterexample.py`].

---

## Part C. Toward an explicit lift (not used in Part A)

**A chart-free candidate route to [I5] (conjecture).** The sextic Puiseux–Gram lift asks
for `F_c = Zᵀ G Z` with `G ⪰ 0`, in the fixed basis `Z = {u^{i/6} v^{j/6} ψ_l}`,
which is unchanged in shape by rescaling. Proved so far:

- it is sound and closed;
- it is exact if and only if it contains the extreme sources;
- four genuine C1/C2 boundary sources have exact rational certificates;
- a continuation theorem reduces the C1/C2 families to one nondegeneracy
  property, (ND);
- it has an explicit rational strictly feasible point and coefficient rank
  406 on both contact families;
- the hardest sources found have exact certificates, and the limit problems
  at the C2 ends `c → 0` and `c → ∞` are solved in closed form on whole faces.

Numerically, every tested source reaches at least 0.999 of the validity limit.
If (ND) holds, it replaces [I5] and gives an explicit lift. See
[`EXACT_LIFT_CANDIDATE.md`](EXACT_LIFT_CANDIDATE.md).

**The explicit lift.** The tilted-CHSH projection has an exact order-4 lift.
The canonical reduction and the eight-copy assembly are exact [paper].

| candidate cone for the canonical sources | status |
|---|---|
| fixed-function cones (polynomial, half-line 33, ramified 49) | strictly too small [paper] |
| ramified cone plus its mirror (98) | strictly too small, by an exact rational certificate [checked: `lift-experiment/`] |
| sextic lift `L₆` | candidate exact lift: sound, closed, strictly feasible; exact certificates at the hardest sources |

## The proof in one paragraph

Scheiderer turns "finite lift" into "certificates that survive independent
degenerations of the inequality and the experiment". Qubit atoms and finitely
many extreme Bell inequalities reduce this to parent forms of CHSH plus event
penalties. Events and the two closed faces are immediate. For an inequality
touching the quantum set at a point, the CHSH part has a single violating
direction, stabilized by the penalties like current through a resistor network.
The determinant balances like forces, and AM–GM certifies it everywhere at
once. What remains is to remove one mixed multiplier; that is the companion's
hard lemma. Part C records a candidate chart-free way to do it in one stroke, with a
scale-blind basis of sixth-root monomials.
