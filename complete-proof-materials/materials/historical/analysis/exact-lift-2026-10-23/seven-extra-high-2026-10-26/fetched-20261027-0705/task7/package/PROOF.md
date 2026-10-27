# K364 full moment-dual: exact quantifier refutation and finite corrected gate

**Scope.** The already materialized main `8d23998` global pencil and the genuine native C1/C2 source maps are inputs. No fixed-cone matrix is rebuilt, no SDP is solved, and no numerical search is used. The original question asks for `⟨Tc(θ),y⟩=0` for **all** positive truncated moments. This is mathematically impossible. The necessary and sufficient condition for the proposed PSD admission is `⟨Tc(θ),y⟩≥0` for all such moments. This note establishes the falsehood of equality, gives a compact and rationally nondegenerate universal dual test section, and isolates the sharp surviving global assertion. It does NOT prove H1/H2 or disprove K364 admission.

## 1. Actual feature and moment coordinates

Use `z_{ij,l}=u^{i/6}v^{j/6} ψ_l` with `i=0,…,6`, `j=0,…,12`, `l∈{W,X,Y,Z}`, giving 364 fixed rows. Let `y_{IJ,ab}` be a real linear functional on all 3250 monomials `u^{I/6}v^{J/6}ψ_aψ_b`, `I=0,…,12`, `J=0,…,24`, `a≤b` in the coefficient convention of the supplied compact T map. Set

`M(y)_{(i,j,l),(i',j',l')} = y_{i+i',j+j',min(l,l'),max(l,l')}`.

The symmetry/off-diagonal factors are incorporated in the coefficient map `T`, not in this moment entry. Thus for any `P⪰0`, `⟨A(P),y⟩=tr(P M(y))`. In particular the closed Gram-image cone `C=A(S⁺₃₆₄)` has dual `C*={y:M(y)⪰0}`.

## 2. Exact counterexample to the proposed universal equality

Evaluate at the pure test `u=v=1` and pure state `ψ=e_W`, then divide by 91:

`y_{IJ,ab}=1/91` for each `0≤I≤12,0≤J≤24, (a,b)=(W,W)`, and zero for every other key.

There are precisely `7·13=91` W-feature rows, all equal one at this point; hence `M(y)=zzᵀ/91⪰0`, `rank M=1`, `tr M=91/91=1`. By the original SG2 parent,

`⟨Tc,y⟩ = ( (1+1²) F_c(1,1;e_W) )/91 = (16 λ+8m+16C+8D+4d)/91`.

This is strictly positive for every genuine strict oriented C1 and strict C2 source because `m,C,d>0` (and D>0 for strict C2). For the genuine strict C2 seed `(d,n,m,C,A,D,α)=(1,1/8,3/8,1/6,5/8,3/8,1/6)` with `λ=β=0`, the exact value is `38/273 > 0`. For the supplied physical C1 k=4 test source, it equals `2316936438928/6385857947469 > 0`. That latter source is an additional algebraic check, not used to assert H1's oriented condition.

**Therefore:** The requested universal equality is false. Equality is only possible for dual supporting functionals annihilating the particular parent, not for all `M(y)⪰0`.

## 3. Trace-normalized section is compact in every one of the actual 3250 coordinates

Define `Y={y:M(y)⪰0, tr M(y)=1}`. Every key `(I,J,a,b)` is an entry of M: choose `i=min(I,6)`, `i'=I−i` and `j=min(J,12)`, `j'=J−j`, and feature state labels `a,b`. The resulting indices fall into the allowed feature ranges. Hence `y→M(y)` is injective. A unit-trace PSD matrix has diagonal entries in `[0,1]` and off-diagonal entries in `[-1,1]` by its 2×2 minors. Every y-coordinate therefore lies in `[-1,1]`. Y is a closed bounded semialgebraic section and is compact. It is nonempty by §2. This argument uses **all zero-target coefficient keys**: they are part of y and appear in M just as the nonzero source columns do.

## 4. A rational strictly positive moment anchor

For each `r=1,…,7`, `s=1,…,13`, `l∈{W,X,Y,Z}`, take the rank-one evaluation matrix `z(r^6,s^6;e_l) z(r^6,s^6;e_l)ᵀ`. Their sum `M₀` is **positive definite**. Indeed, the 364 evaluation vectors are the columns of `(V₇⊗V₁₃⊗I₄)` in a fixed order, where the entries are `V₇[r,i]=r^i` and `V₁₃[s,j]=s^j`; both Vandermonde determinants are products of nonzero positive integer differences. The resulting moments are exactly

`y⁰_{IJ,ab}=δ_ab (Σ_{r=1}^7 r^I)(Σ_{s=1}^{13} s^J)`.

Thus y⁰ is integral and `M(y⁰)≻0`. Normalizing by its positive integer trace yields a rational element `y⁰_norm∈Y` with full-rank M. This is an explicit truncated-moment construction, **not** an unproved representation theorem for arbitrary moments.

## 5. Corrected separation and rational-witness theorem

For any fixed physical strict source c and `f=Tc`, closedness of the already independently verified coefficient Gram image gives

`c∈K364  ⇔  min_{y∈Y} ⟨f,y⟩≥0`.

No condition that `⟨f,y⟩` vanish for every y is involved. A negative minimizer supplies a separator. If c has rational entries and the minimum is negative, one can choose a **rational** negative separator with a strictly positive moment matrix: mix any negative y with `ε y⁰_norm` for sufficiently small positive rational ε, obtaining `M≻0` with the pairing still negative; approximate its 3250 coordinates by rationals, retaining `M≻0` and negative pairing; normalize trace to one. No moment representing-measure theorem is invoked. Consequently K364 failure at an arbitrary real strict native parameter automatically has some sufficiently nearby rational strict native parameter and a rational negative dual witness: the strict source domain is open and rational functions of its parameters are continuous. At the oriented C1 face k=1, approach by k<1 while preserving strict native signs.

For a rational witness, **all** the moment-PSD conditions are checkable by exact rational symmetric-matrix elimination, and the source pairing uses only the 120 nonzero fixed T entries. Conversely, the absence of a rational negative witness has not been proved. Neither the native stationary equations, which enforce the six contact first jets in source-adapted variables, nor scalar determinant nonnegativity proves positive pairing against every PSD truncated pseudo-moment. The demonstrated rank-one evaluations have representing measures, but arbitrary members of Y need not; replacing them by evaluations is logically invalid.

### Remaining exact assertions

- **H1:** `∀θ∈Ω_C1,0<k≤1, ∀y∈Y`, `Σ_a c_a(θ) Σ_key T[key,a] y_key ≥0`.
- **H2:** `∀θ∈Ω_C2, ∀y∈Y`, the same inequality.

Both retain all four state coordinates, all 3250 coefficient keys (including target zeros), source-only stationarity and the actual native inequalities. These are the two genuine candidate-specific gates. The failure of one preselected C2 syzygy correction does not supply a negative y in Y. This note is a corrected moment-dual research result, not a resolution of H1 or H2.

**Reproducibility:** `verify_exact_dual.py` reads the frozen user-supplied native compact T columns, checks the exact rank-one moment, its trace, the source pairings, all 3250-key reachability ranges and the 120 nonzero source coefficients. `EXACT_EVIDENCE.json` is the deterministic output. No SDP, sweep, or Git write.
