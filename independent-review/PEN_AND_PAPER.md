# A pen-and-paper route through the finite-lift proof

The existence proof splits into a short core that can be checked by hand and one
isolated lemma that carries the long appendices. This note gives the core in
eight steps, each with the computation a reader can do on paper, and states the
remaining lemma precisely. Every displayed identity was also checked by script
(`checks/determinant`, `checks/circuits`, `checks/physics`).

## 0. What has to be shown

Q222 ⊂ ℝ⁸ is the set of behaviours `q = (a₀,a₁,b₀,b₁,E₀₀,E₀₁,E₁₀,E₁₁)` from
finite-dimensional states and binary measurements. Scheiderer's criterion
(SIAM J. Appl. Algebra Geom. 5, 2021, Cor. 3.19) reduces "Q222 is the projection
of a spectrahedron" to one algebraic statement. For every real closed field
`K ⊇ ℝ`, every affine `f ≥ 0` on `Q_K` and every `q ∈ Q_K`:

    [f, q] = f₀⊗1 + Σ fᵥ⊗qᵥ  is a sum of squares in  K ⊗_ℝ K.

Infinitesimals in the left factor describe a Bell inequality approaching a
degenerate limit, and those in the right factor an experiment doing the same,
at an unrelated rate. A certificate may take roots and inverses only of
elements of one factor. It may never divide by a quantity that mixes the
inequality with the experiment.

## 1. Real qubits suffice

By Jordan's lemma, two binary observables per party decompose into blocks of
size at most two, in which `A₀ = B₀ = Z`, `A₁ = cos a·Z + sin a·X` and
`B₁ = cos b·Z + sin b·X`. The observables are real, so a complex state is a
mixture of its real and imaginary parts. By Carathéodory, Q222 is the convex
hull of at most nine points `Φ(a, b, ψ)` with `ψ ∈ S³`. Write `u = tan(a/2)`,
`v = tan(b/2)`, `D_u = 1+u²` and `D_v = 1+v²`. (Chapter 2 of the companion.)

## 2. Every Bell inequality is CHSH plus event penalties

Let `S = E₀₀ − E₀₁ + E₁₀ + E₁₁` and `H = 1 − S/2`, so `H = −1` at the PR box.
The PR box and its eight deterministic neighbours are affinely independent, so
every slack that is negative at the PR box has a unique form

    f = d·H + 4 Σᵢ cᵢ·pᵢ,   d = −f(PR) > 0,   cᵢ = f(vᵢ)/4 ≥ 0,

where `pᵢ` is the probability of the event lost by neighbour `vᵢ`:

| weight | λ | n | m | C | A | D | α | β |
|---|---|---|---|---|---|---|---|---|
| event `(ab\|xy)` | (++\|01) | (+−\|00) | (+−\|11) | (−+\|10) | (−+\|00) | (−+\|11) | (−−\|01) | (+−\|10) |

A slack nonnegative at every PR box is a nonnegative combination of
probabilities, which are already squares. (Chapter 4.)

## 3. The Bell operator has one violating direction

For a state `ψ = (W,X,Y,Z)` in the basis `|00⟩,|01⟩,|10⟩,|11⟩`, after clearing
denominators, `F = D_u D_v ‖ψ‖² f(q)` is a quadratic form `ψᵀPψ` with

    P = Σₖ aₖ gₖgₖᵀ − ℓℓᵀ,    ℓ = √(2d D_u D_v)·(0,0,0,1).

Each event contributes `aₖgₖgₖᵀ`, a product projector. The weights are
`a = (4λ, 4n, 4m, 4C, 4A, 4D+2d, 4α+2d, 4β+2d)`, and the vectors are

    g_λ = √D_u (1,v,0,0)        g_n = √(D_uD_v) (0,1,0,0)
    g_m = (v,−1,uv,−u)          g_C = √D_v (u,0,−1,0)
    g_A = √(D_uD_v) (0,0,1,0)   g_D = (u,uv,−1,−v)
    g_α = √D_u (0,0,v,−1)       g_β = √D_v (0,1,0,u).

The terms proportional to `d` add up to `d·D_uD_v·(I − B/2)`, with the CHSH
operator `B = A₀B₀ − A₀B₁ + A₁B₀ + A₁B₁`. Landau's identity gives
`B² = 4I + 4 sin a sin b·(Y⊗Y)` with the real `Y = [[0,−1],[1,0]]`. Hence the
eigenvalues of `I − B/2` are

    1 ∓ √(1 ± sin a sin b).

Exactly one is negative when `sin a sin b ≠ 0`. Its depth `√(1 + |sin a sin b|) − 1`
is about half the product of Alice's and Bob's incompatibilities. Those two
factors vanish independently as each party's measurements become compatible.
This is the physical origin of the two independent infinitesimals and of every
singular case in the proof.

Restricted to `Z = 0`, the positive part is positive definite whenever `d > 0`
and `uv ≠ 0`. So `P ⪰ 0` holds exactly when `det P ≥ 0`. (Proof.tex §2.)

## 4. Validity is an effective-resistance bound

Let `T = Σ aₖgₖgₖᵀ` and `s = ℓ·ψ`. For any "flow" `t ∈ ℝ⁸` with
`Σ tₖgₖ = ℓ`, expanding the squares gives

    F = Σₖ aₖ (gₖ·ψ − (tₖ/aₖ)·s)²  +  s²·(1 − Σₖ tₖ²/aₖ).

By Thomson's principle, the smallest cost `Σ tₖ²/aₖ` is the effective resistance
`ℓᵀT⁻¹ℓ` of the violating direction, in a network whose edges are the event
vectors with conductances `aₖ`. So, for fixed measurements:

- a Bell slack is valid exactly when `ℓᵀT⁻¹ℓ ≤ 1`;
- any flow of cost at most 1 is an explicit sum-of-squares certificate;
- the matrix determinant lemma gives `det P = det T·(1 − ℓᵀT⁻¹ℓ)`.

## 5. The determinant is an AM–GM object

    det P = 16 D_u² D_v² u² v² · L(u², v²),
    L(U,V) = h + Σ_{e ∈ {−1,0,1}² ∖ 0} q_e U^{e₁} V^{e₂}.

The eight outer coefficients are polynomials in the weights with nonnegative
coefficients. The corners factor, for example
`q₍₁,₁₎ = 4Cβ(2A+2α+d+2m)(2D+d+2λ+2n)`. Only the centre `h` can be negative.
Both `det T` and `ℓᵀ adj(T) ℓ` contribute to all nine monomials, and their
negative parts cancel everywhere except at the centre.

At a stationary contact `(U_c, V_c)`, the weights `wₑ = qₑ U_c^{e₁} V_c^{e₂}`
satisfy `h = −Σ wₑ` and `Σ wₑ·e = 0`. A balanced nonnegative vector splits into
minimal circuits. On the 3×3 square there are twelve:

- four opposite pairs;
- four triangles made of a corner and the two axis points opposite it, with
  weights (1,1,1);
- four triangles made of two adjacent corners and the axis point opposite,
  with weights (1,1,2).

Three identities turn each circuit into squares:

    p + p⁻¹ − 2 = (√p − 1/√p)²
    a³ + b³ + c³ − 3abc = ½(a+b+c)[(a−b)² + (b−c)² + (c−a)²]
    a⁴ + b⁴ + c⁴ + d⁴ − 4abcd = (a²−b²)² + (c²−d²)² + 2(ab−cd)²

Write `X = U_c⁻¹ ⊗ U` and `Y = V_c⁻¹ ⊗ V`. These are pure tensors, so their
cube and fourth roots are taken separately in each field. Hence `L`, and with
it `det P`, is a separated sum of squares at every contact source. Physically,
stationarity is a balance of forces in the plane of log half-angles. (Chapter 4.)

## 6. One identity certifies the whole form, up to one multiplier

With `y = Tψ − s·ℓ`, the matrix determinant lemma gives

    det(T) · F = yᵀ adj(T) y + det(P) · s².

To prove it, use `yᵀT⁻¹y = ψᵀTψ − 2s² + s²·ℓᵀT⁻¹ℓ`. Multiply by `det T`, add
`det(P)s² = det(T)·s²·(1 − ℓᵀT⁻¹ℓ)`, and the result is `det(T)·(ψᵀTψ − s²)`.

By Cauchy–Binet for the adjugate:

    adj(T) = Σ_{|S| = 3} a_S · w_S w_Sᵀ.

Here `a_S` is a product of three Bell coefficients, which is pure and
nonnegative. `w_S` is the generalized cross product of three event vectors,
which is pure on the test side. So `yᵀ adj(T) y` is a separated sum of squares
for every source and test. By Step 5, so is `det(P)·s²` at every contact
source.

**Hence `det(T)·F` is a separated sum of squares for every contact source, with
no case analysis.**

## 7. The one hard lemma: removing det(T)

Also by Cauchy–Binet, `det T = Σ_{|S|=4} a_S·det[g_S]²`, a quantity that mixes
the inequality with the experiment. `K⊗K` is not a field, so this multiplier
cannot be divided out. The criterion of Step 0 forbids exactly this division.
In the fraction field, every valid slack is a sum of squares by Artin's
theorem; the whole content of the criterion is the absence of mixed
denominators.

**Lemma 𝒞 (canonical assertion).** For every C1 or C2 contact source over `K_c`
and every test over `K_t`, `F` itself is a separated sum of squares.

In flow language (Step 4), the optimal flow `t*ₖ = aₖ·gₖᵀT⁻¹ℓ` carries the
forbidden denominator `det T`. A flow supported on four events `S` is pure,
because `det[g_S]` lives on the test side alone. But at a contact the cost
must be exactly 1, so the flow must be the optimal one. The optimal flow is
unique, and at a generic contact it uses every event with positive weight. So
no fixed choice of four events works throughout.

The companion proves 𝒞 by a partition into charts by the relative scale of
source and test. On each chart, conductor identities and strict-residue Gram
matrices (Lemma TF2) produce the certificate. The review found this argument
complete in substance; see `README.md` for the editorial defects and the steps
it re-derived.

**A route to a short proof.** Find, uniformly over sources and tests, a flow
`t ∈ (K_c ⊗ K_t)⁸` with `Σ tₖgₖ = ℓ` and `1 − Σ tₖ²/aₖ` a separated sum of
squares. This is open.

One constraint on any such route is known. The C1/C2 split cannot be removed.
The source `c = 1/4, v_c/u_c = 2, b = 3/2, t = 1/13, r = 1/32` in the native C2
coordinates of `c2_short.tex` has formal C1 image with `π ≈ 0.27`, and that
image has a negative Hessian direction at the contact. So C2 does not reduce to
a C1 theorem with negative β, except in the thin-n slices the companion already
uses.

## 8. Assembly

Every valid slack is a nonnegative combination of at most nine extreme
sources, and every behaviour is a mixture of at most nine qubit atoms. If
`[fⱼ, qᵢ] = Σₖ hᵢⱼₖ²`, then

    [Σ aⱼfⱼ, Σ wᵢqᵢ] = Σ (√aⱼ ⊗ √wᵢ · hᵢⱼₖ)²,

because `Σ wᵢ = 1` preserves the constant term. Scheiderer's criterion then
gives a finite lift. Tarski transfer gives one with real algebraic entries.

## What a reader checks, and how

| Step | By hand | By script |
|---|---|---|
| 1 Real qubits | Jordan blocks, real/imaginary split | none needed |
| 2 CHSH plus penalties | affine independence of nine vertices | none needed |
| 3 One violating direction | Landau's identity, a 4×4 spectrum | `checks/physics` |
| 4 Effective resistance | expand one square | `checks/physics` |
| 5 Determinant and AM–GM | three identities, list of 12 circuits | `checks/determinant`, `checks/circuits` |
| 6 Uniform identity | determinant lemma, Cauchy–Binet | `checks/physics` |
| 7 Lemma 𝒞 | the chart argument, as a reference | companion appendices; `checks/r6-remote-charts`, `checks/thin-layer` |
| 8 Assembly | one line | none needed |
