# Independent review of the finite-lift proof

This folder records an outside review of the existence proof in `latex/` and of
the explicit-lift material in `extracted-lifts/`. It adds no claims of its own to
the paper. It breaks the argument into stages, states the idea that organizes it,
lists what was checked and what was found, and follows the explicit lift as far
as it currently goes.

The proof sources under `latex/` were reviewed as they stand at commit `ee85fd1`;
they are unchanged from the snapshot the review started on. Three folders under
`extracted-lifts/canonical-completion/` (`uniform-contact`, `contact-completion`,
`fixed-c2-gates`) appeared during the review and were not audited.

[`q222-lift-anatomy.html`](q222-lift-anatomy.html) is the same review as a single
page, with two diagrams. GitHub shows HTML as source, so download it and open it
in a browser.

[`PEN_AND_PAPER.md`](PEN_AND_PAPER.md) is the whole proof from first
principles, in physics order. Each step is marked as proved there, proved in the
paper, checked by script, or conjectural. Along the way:

- Bell operators become two-mode quadratic forms.
- Landau's identity gives a single violating direction.
- Validity at fixed measurements is an effective resistance.
- AM–GM appears as a force balance.
- The single hard lemma is isolated.

## Verdict

- **Existence proof: sound wherever it was tested.** Every core lemma held under
  adversarial review. The key determinant identity was re-derived symbolically for
  all weights. The audit of the ~31,000 appendix lines found editorial defects but
  no mathematical gap.
- **Sketched steps.** The two load-bearing steps that are only sketched in the text
  were tested separately, and all of them hold. These are the thin-layer reuse of
  C1 certificates on a negative-β source, R10 coverage of the Z = 1 contact chart,
  and the HI radical-division pivots (`c2_short.tex:6849–6909`).
- **Explicit lift: no exact pencil for all eight coordinates yet.** The repo
  reduces the problem exactly to one 9-dimensional canonical cone. Every fixed
  selector family tried so far is provably too small, so the full-set pencils of
  order 264 and 392 are strictly larger than Q222.

## The organizing insight

A finite lift needs positivity that ignores how fast things degenerate.
Scheiderer's criterion (2021, Cor. 3.19) asks for every valid Bell slack, at
every quantum point, to be a sum of squares in `K_c ⊗ K_t`. Infinitesimals in the
two fields encode the inequality and the experiment collapsing at independent
rates, so a certificate may only take roots and inverses inside one field.

Q222 needs only two certificates that survive this.

1. **Strict real residue.** A bounded Gram matrix whose standard-part matrix is
   positive definite is a tensor sum of squares (`03_tensor.tex`, TF1–TF2).
2. **Monomial AM–GM.** Monomials `U_c⁻¹ ⊗ U` are pure tensors, so their roots
   split across the two fields.

Q222's geometry makes the second one sufficient for the determinant.

- **One negative square.** Every nonlocal slack has a single negative square, the
  CHSH core `−2d·D_u·D_v·Z²`, and three positive transverse directions, so the
  4×4 state matrix is positive semidefinite exactly when its determinant is
  nonnegative.
- **One interior point.** The determinant is a Laurent polynomial in `(u², v²)`
  on the 3×3 exponent square. Its eight boundary coefficients are products of
  event weights and never negative, so only the centre can be negative.
- **Contacts are AM–GM.** At a stationary contact the weights balance and split
  into 12 circuits, each an exact AM–GM sum of squares.

The appendices then recover the full quadratic form from the determinant. They do
this chart by chart, without the forbidden mixed inverse, until engine 1 applies.

## Stages of the proof

| Stage | What it establishes | Where | Status |
|---|---|---|---|
| Model | `Q222 = conv Φ(S¹×S¹×S³)`, ≤ 9 real two-qubit atoms | `02_model.tex` | proved |
| Criterion | finite lift ⇔ `[f,q]` is a sum of squares in `K⊗K` | `Proof.tex` §1, arXiv:2004.04196 | statement and hypotheses match the source |
| Engine | strict-residue Gram ⇒ tensor SOS | `03_tensor.tex` | holds |
| Sources | extreme slacks are events, two closed faces, or C1/C2 contacts | `04_geometry.tex` | holds; all 256 omission patterns brute-forced |
| Determinant | `det P = 16D_u²D_v²u²v²·L(u²,v²)`, 8 outer coefficients ≥ 0 | SG5–SG6 | asserted in text; re-derived symbolically |
| AM–GM | balanced weights split into 12 circuits | `Proof.tex` §2 | holds; circuit list complete |
| Full form | canonical assertion 𝒞 at every scale and chart | App. A, B, C; `noncommuting_refined.tex` | audited; see findings |
| Assembly | mixtures (M2), Scheiderer, Tarski ⇒ algebraic entries | `Proof.tex` §4, `07_cover.tex` | holds |

## Findings in the companion

Editorial defects. In each case the content is present elsewhere or was rebuilt.

- **Stale "remaining task".** `c2_short.tex:14268–14278` ends with "the only
  remaining task is the contact chart u,v,w_test>0 infinitesimal, Z=1". That chart
  is already proved for every R10 source. R10.1 covers n/e² not infinitesimal, the
  thin layer C.17.3 covers n/e² positive infinitesimal, and n = 0 maps to C1
  through OP4.
- **R6.5–R6.7 never printed.** These remote-chart matrices are cited as printed
  (`c2_short.tex:12397–12423`, `12872–12878`) but appear nowhere, including the
  PDF. They were rebuilt from the parent, and all three claims hold exactly (see
  `checks/r6-remote-charts`). The stated reason for R6.5 is incomplete, and R6.7
  needs the redundant roots `pZ`, `qZ`.
- **Bare labels.** Nine `\label` lines with no content close `c2_short.tex`. The
  whole-R10 theorem promised at line 13888 is never stated as one result.
- **HI leaf contradiction.** `c2_short.tex:7074` says the HI leaf is "narrowed but
  not closed". The proof printed earlier in C.9.2 (`6712–6915`) contradicts this.
- **Undefined tags.** NSR3–5, WS3, UZ4–6 and CL7 are never defined. Their content
  is printed under other headings, or not needed.
- **Stale C1 phrases.** These include "subject to the exact checks … recorded
  alongside" and "revision-6", plus sentences missing their antecedent. The
  mathematics is printed and was recomputed.
- **Inconsistent copies.** The two copies of `Explicit_Lifts.tex` differ. The
  README gives the C1 contact pencil four blocks, but its file has five (order 59).

Steps printed only as sketches:

- **Thin-layer import** (C.17.3). It maps a C2 source to a formal C1 point with
  β < 0 and reuses the CA/CS/SE constructions. The C1 text uses β ≥ 0 in exactly
  two places (the SE τ bound and the SE4 ν reservation), and the text replaces
  both explicitly. TL3, TL5, TL6 and the zero value and jet at contact were
  checked exactly. **Holds.**
- **HI radical division** (`c2_short.tex:6849–6909`). Positive pivots are asserted
  without a printed Gram matrix. A 10×10 Gram on the roots N_c, N_m, E₁–E₄, F_p,
  F_v, R₁₃, R₂₃ has a positive-definite residue for every r₀ ≥ 0, including R₀ = 0
  and n ≈ 0. **Holds.** The averaging weights are left unspecified in the text;
  any positive choice works.
- **IW16 arrays** (`c1_short.tex:3676–3687`). The residues rely on second-order
  cancellations. They were recomputed and are correct.
- **Big identities.** Large identities (EN8, EN16, NC14) are said to follow by
  "direct multiplication", and the repo ships no scripts for them.

## The explicit lift

The existence proof yields matrices only through compactness, so the finite
testing space it uses is never named. The repo's exact route has two steps.

1. **Canonical reduction.** Eliminating Alice's angle makes a source valid exactly
   when, for every Bob angle v, some 2×2 `J` makes two 4×4 matrices
   `G_J = [[P+K, J],[Jᵀ, N+K]]` and `H_J = [[P+E, K−E−J],[K−E−Jᵀ, N+E]]`
   positive semidefinite (`canonical_completion.tex:98–163`).
2. **Full-polar assembly.** A strictly feasible exact pencil of order n for that
   9-dimensional cone gives Q222 as
   `{q : ∃ Y_g ⪰ 0, ℬ*Y_g = 0, 𝒜*Y_g = M_gᵀ(1,q), g ∈ {±1}³}`, of order 8n
   (`full_polar_assembly.tex:137–225`).

The missing piece is the exact canonical pencil. Each fixed family fails at an
explicit rational source:

| Family for J(v) | Order | Witness `(d,λ,n,m,C,A,D,α,β)` |
|---|---|---|
| polynomial on ℝ | 8 | C1-face source `(359,0,0,75/16,60,416,41/2,831/2,151/2)`; the forced J gives a minor of −2645/4 |
| half-line quadratic | 33 | `c̄ = (1, 2⁻³⁴, 2⁻¹⁴, 2⁻²⁴, 1/576, 161, 1/2, 76, 1/250)` |
| ramified `√v`, `v^{3/2}` | 49 | `c† = πc̄ = (1, 1/576, 161, 2⁻²⁴, 2⁻³⁴, 2⁻¹⁴, 1/2, 1/250, 76)`, the Alice↔Bob mirror |

**Available now.**

- **Tilted-CHSH projection.** With `t = a₀` and `s = E₀₀+E₀₁+E₁₀−E₁₁`, the
  projection is exactly `[[√2+t, s/2],[s/2, √2−t]] ⊕ (1−t) ⊕ (1+t) ⪰ 0`.
- **Size bound.** Any lift has total PSD order at least 9.
- **Local pencils.** Local source-cone pencils of orders 62, 669 and 731.
- **New since the snapshot (not audited here).** A fixed certificate of order 283
  for the stationary R1 family on `|u| ≤ 9/256`, `|v| ≤ 1/16`, joined to order 293
  with a Bob complement. An order-635 construction valid at all measurements for
  the native C2 source envelope. The repo itself states that global canonical
  coverage remains open.

**The symmetric repair also fails.** The order-98 cone `C_ram + πC_ram`
restores the Alice↔Bob symmetry and would assemble to order 784. It is strictly
smaller than `C_can`. An exact rational certificate proves it: a valid rational
source is separated by a rational functional built from exactly positive
semidefinite dual matrices. High-precision solves show the cone reaches only
0.588 of the valid CHSH weight on a π-symmetric multiscale source. See
[`LIFT_EXPERIMENT.md`](LIFT_EXPERIMENT.md).

In four of the seven gap cases the mirror adds nothing, so the missing testing
functions concern scale rather than symmetry.

**A candidate exact lift.** [`EXACT_LIFT_CANDIDATE.md`](EXACT_LIFT_CANDIDATE.md)
proposes the sextic Puiseux–Gram lift: one Gram matrix in the fixed basis
u^{i/6}v^{j/6}ψ, with no charts. With half-powers (k = 2) it is not exact; 9 of
302 test sources fall short. With sixth roots (k = 6) no tested source falls
below 0.999 of the validity limit, including every source that defeats the
earlier candidates.

Part of exactness is now proved:

- the lift is sound and closed;
- exactness reduces to the extreme sources;
- four genuine C1/C2 boundary sources have exact rational certificates;
- a continuation theorem reduces the C1/C2 families to a single
  nondegeneracy property, (ND).

What remains open is (ND) and the two closed faces.

## Checks in this folder

Every script uses exact rational or symbolic arithmetic unless noted. Run Python
with SymPy 1.14 (as pinned in the repo), from the script's own folder.

| Folder | What it checks | Run | Result |
|---|---|---|---|
| `determinant/` | `det P` factorization; 8 outer coefficients have only nonnegative coefficients in the weights; only the centre is indefinite | `python det_symbolic.py`, `ruby det_check.rb` | holds for all weights |
| `circuits/` | the 12 minimal circuits on `{−1,0,1}²∖0` and the two AM–GM identities; SG9 (coercive ⇔ omitted events pairwise compatible) over all 256 supports | `python circuits.py`, `python sg9.py` | 12 circuits; 0 mismatches |
| `canonical-completion/` | `F_c = D_u·Φ_c`; the split identity for every J; Alice↔Bob transport; Motzkin half-line identity; the −2645/4 minor; order-49 kernel count (assembles to 392 and 640) | `python independent_checks.py ../../../extracted-lifts/canonical-completion` | all true |
| `r6-remote-charts/` | rebuilds the unprinted R6.5–R6.7 and checks the source residues (RB4) and each strictness or determinant claim | `python src_check.py`, `python remote_charts.py`, `python r66_minors.py` | all claims hold |
| `thin-layer/` | TL3 differences, R10.2, TL2, TL5, TL6, and zero value and first jet of both parents at the OP4 contact | `python tl3.py`, `python tl6.py` | all residuals 0 |
| `tilted-chsh/` | the order-4 projection pencil against its set on a grid, and on random two-qubit behaviours (floating point) | `ruby chsh_check.rb` | 0 mismatches, 0 violations |
| `lift-experiment/` | exact rational proof that `C_ram + πC_ram ≠ C_can`: an AM–GM validity certificate for a rational source, and an exactly positive semidefinite rational dual that excludes it; plus the high-precision results | `python verify_cert.py ../../../extracted-lifts/canonical-completion/ramified-completion/data/ramified_order49_pencil.json` | 17/17 PASS |
| `gram-lift/` | the Puiseux–Gram lift `L_k`: standalone checker for d_k/d* and the δ-margin test, the symmetry check, the k = 2 and k = 6 sweep results over 302 sources, and the test showing the scalar circuit structure does not predict `L₂` | `python gram_lift.py --k 6 --case <name>`, `python symcheck.py`, `python circuit_types.py` | k = 2: 9 gaps; k = 6: none below 0.999 |
| `gram-lift/boundary/` | exact rational `L₆` certificates at four genuine C1/C2 boundary (contact) sources: validity, exact zero at contact, exact Gram identity, positive definite after removing the contact vector | `python verify_boundary.py` | 4/4 sources, 56/56 checks |
| `physics/` | identities behind `PEN_AND_PAPER.md`. Landau's identity and the spectrum of I − Ŝ/2; the parity Γ = Y⊗Y that commutes with correlators and anticommutes with marginals (`check7_parity.py`); the two-mode normal form of every Bell operator in the Bell basis (`check8_two_mode.py`); events as product-state projectors; the Thomson flow identity and the effective-resistance test; Cauchy–Binet for adj(T) and det(T); `det(T)·F = (Pψ)ᵀadj(T)(Pψ) + det(P)s²`; which determinant coefficients the CHSH term reaches; a C2 source whose formal C1 image is not positive | `python check6_tsirelson.py`, `check4_thomson.py`, `check3_cauchy_binet.py`, `check12_decomp_sos.py`, `check5_monomials.py`, `check_c2_to_c1_counterexample.py` | all true; counterexample confirmed |

The repo's own checkers were also run on Python 3.9 with SymPy 1.14, and all
passed. Three needed a small compatibility shim:

- `angular-order49/check_angular_pencil.py` uses `newline=` in `write_text`
  (Python 3.10+).
- `contact-selector/check_c1_contact_selector.py` uses `sys.set_int_max_str_digits`
  (Python 3.11+).
- `tilted-projection/build_and_check_tilted_extraction.py` imports its sibling
  module, which fails under `python -I`.

## Method and limits

The chapters, the appendices, the extraction chain and the shipped data were each
read in full or by targeted search. Each load-bearing claim was then checked
against the source text, with every computable identity recomputed exactly
(scripts in `checks/`). Each apparent gap was re-examined before being recorded.

Not covered:

- a line-by-line re-derivation of every appendix leaf;
- the noncommuting-residue theorem (Thm 6.2) in detail;
- the three canonical-completion folders added during the review.
