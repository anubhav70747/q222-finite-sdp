# Independent review of the finite-lift proof

The [Sun–Wang–Yu matching](SUN_WANG_YU_CONTACT_MATCH.md) gives a short alternative proof of finite-contact coercivity, uniqueness and full C1/C2 curvature. It includes the full affine coefficient map, measurement/state scaling, a [bounded independent review](SUN_WANG_YU_CONTACT_REVIEW.md), and [exact symbolic checks](checks/sun-wang-yu-match/result.json). The native bounds are tied; their frustrated-case lemmas apply. This completes the matching, not universal fixed-library admission.

This folder records an outside review of the existence proof in `latex/` and of
the explicit-lift material in `extracted-lifts/`, followed by new constructions
on `main`. The review breaks the argument into stages, states the idea that organizes it,
lists what was checked and what was found, and follows the explicit lift as far
as it currently goes.

The [existence-proof check](EXACT_LIFT_ROUTE.md#the-original-existence-argument-has-a-different-requirement)
separates the original universal tensor-SOS argument from extraction of one
fixed finite pencil. For the original argument it isolates the thin-R10
bounded-matrix implication as the key verification step. The
[explicit bounded conductor](SE_BOUNDED_CONDUCTOR.md) and the
[two-power conductor](TWO_POWER_CONDUCTOR.md) now supply that eight-root
matrix under their stated hypotheses. Its thin-R10 interface reconciliation
is complete in [V02](../complete-proof-materials/CHECKLIST.md), including the
ordinary–small C2 branch and the CA/SE base. The [EL/NU endpoint proof](C1_EL_NU_ENDPOINTS.md) now
supplies the complete scalar, both state normals, bounded coefficient
construction, and complementary testing charts for those two physical C1
source domains. Its eight symbolic modules have been replayed locally.
The completed 22 October C1/C2 reviews and original proof bodies have now been
[recovered and published](../complete-proof-materials/). They already cover the
remaining source domains and global assembly. The later reports below are
scoped rechecks, not evidence that those earlier proofs were absent. The
[final conjunction](../complete-proof-materials/RECONCILIATION.md#final-qualitative-conjunction-26-october-13-utc)
reconciles the current source coverage with the independently checked global
implication. The [remaining construction](../complete-proof-materials/VERIFICATION_PLAN.md)
concerns universal strict-source PSD admission and exactness of the already evaluated conditional global pencil. The [current frontier](../FRONTIER.md) and [C2 n=0 opening theorem](C2_N0_OPENING.md) give the latest scopes.

New scoped all-testing constructions are [the balanced C1 source cone and native tube](C1_BALANCED_SLICE.md), with exact order-five source/behavior pencils, and [contact-compatible opening at all three zero-weight walls](ZERO_WEIGHT_OPENING.md). Their complete rational data and standalone checkers are included. Universal H1/H2 and independent review of these new results remain open.

For the current construction route, read
[`EXACT_LIFT_ROUTE.md`](EXACT_LIFT_ROUTE.md): a fixed NPA baseline, the
seven-root reduction forced by two classical zeros, the existing closed-face
pencils, and the precise dual inequality to be proved for the contact sources.
The new [equal-angle C2 cone](../extracted-lifts/c2-equal-angle-cone/) covers
two independent contact parameters with total PSD order seven. It contains
the earlier order-three curve, whose extreme sources lie outside `1+AB`.
Both source and behavior matrices have complete elementary proofs and exact
reproduction code.

The stronger [13-word NPA exclusion](C2_W2_EXCLUSION.md) now shows that
the entire finite equal-angle maximally entangled C2 family, including
`n=0`, is missed even after all local words of length two are retained. A separate
[unequal-angle certificate](C2_W2_UNEQUAL_CERTIFICATE.md) gives one extreme
ray admitted by that larger NPA space. The [contact-scalar construction](CONTACT_SCALAR_GRAM.md)
also has a positive-definite three-by-three Gram at every finite native
C1/C2 contact. The [compact-family construction](COMPACT_NATIVE_CONTACT_COVER.md)
completes the full state form on a fixed finite testing grid for any
compact strict native subfamily, with a size that depends on the subfamily.
Universal strict admission in one fixed cone remains open; the verified density theorem then supplies the singular sources. The `r → 0`
end and the zero-weight walls reduce to explicit limit problems, certified
exactly at sample points, and the `c → 0` and `c → ∞` limit problems are
solved in closed form on whole faces ([sextic lift](EXACT_LIFT_CANDIDATE.md)).
The [polynomial obstruction](C2_POLYNOMIAL_OBSTRUCTION.md) proves that
this completion requires more than ordinary polynomial squares for
generic C2 sources; it also reduces a polynomial coupling choice to
twelve free source coefficients.
For C1, a [triangle inequality](C1_POLYNOMIAL_OBSTRUCTION.md) excludes
ordinary polynomial SOS on a nonempty open set of strict native sources.
Its rational counterexample has all retained event coefficients positive.
The same source also shows why a local square completion must preserve
the other deterministic zeros: a natural coupling allocation leaves
a negative remainder there.
The [common-multiplier obstruction](NATIVE_MULTIPLIER_OBSTRUCTION.md)
goes further: one fixed separable polynomial multiplier cannot cover both
native families. On all physical relabelings, even a nonseparable multiplier
positive at the commuting origin fails. Thus a fixed rational testing span
regular there cannot complete the native families.
The [half-power obstruction](C2_HALF_POWER_OBSTRUCTION.md) also excludes
the fixed 36-feature candidate on an open part of strict native C2.
A rational contact and 90 integer moments give an exact positive
compressed dual with negative full-parent pairing. Together these
obstructions motivate the sixth-root features, whose admission is the target.
The [fractional-feature analysis](SIXTH_ROOT_FEATURES.md) separates the scalar
Motzkin limit, which already admits half-powers, from the full state-matrix
problem. A separate exact biquadratic example shows why retaining thirds
and halves need not replace all sixth-root functions.
More generally, the [smooth-feature theorem](SMOOTH_FEATURE_OBSTRUCTION.md)
extends the cubic-scale NPA obstruction to every fixed finite angular-state
Gram library with a polynomial Peano expansion through order three at the
commuting origin. Its coefficients may vary freely with the source. The
common testing frame is essential; the theorem does not exclude arbitrary
finite SDP lifts.
For the strict native C2 counterexample, a further [necessary condition](SIXTH_ROOT_NECESSARY_CROSS_TERMS.md)
forces complementary fractional cross terms in any sixth-root Gram. Its
exact rational bounds constrain the sixth-root construction and leave
admission to be decided.

The [unequal-angle equations](C2_CONTACT_GEOMETRY.md) expose an additional
fourth-degree contact relation and two exact obstructions to simpler
characterizations. An [order-six curve cone](../extracted-lifts/c2-unequal-angle-curve/)
has explicit source and behavior pencils and an elementary three-atom proof.
The full unequal-angle hull is the next extension of this construction.

The [fixed-fiber construction](C2_FIXED_FIBERS.md) now gives an order-four
source pencil for each fixed parameter pair in the maximally entangled C2
family. The [full-family obstruction](C2_MAXIMAL_NPA_OBSTRUCTION.md) proves
that every finite native contact in that family requires certificates beyond
`1+AB`, including unequal angles and `n=0`. The
[C1 obstruction](C1_NPA_OBSTRUCTION.md) and
[partially entangled C2 obstruction](C2_NONMAXIMAL_NPA_OBSTRUCTION.md)
identify further sources for the additional cones to cover. Together they
specify what the fixed finite set of blocks for these whole families must
contain.

The proof sources under `latex/` were reviewed as they stand at commit `ee85fd1`.
On the review branch, two certificate passages were later expanded to display
their conductor steps: the SE local proof in `c1_short.tex` and the
ordinary–small six-row certificate in `c2_short.tex`. Three folders under
`extracted-lifts/canonical-completion/` (`uniform-contact`, `contact-completion`,
`fixed-c2-gates`) appeared during the review and were not audited.

[`q222-lift-anatomy.html`](q222-lift-anatomy.html) is the same review as a single
page, with two diagrams. GitHub shows HTML as source, so download it and open it
in a browser.

[`PEN_AND_PAPER.md`](PEN_AND_PAPER.md) is the proof you can follow, in three
parts:

- **Part A, the proof.** Five short steps from five imported results, each
  stated exactly: Scheiderer's criterion, the qubit model, the extreme-source
  classification, the closed faces, and the hard lemma.
- **Part B, why it works.** Not used in Part A:
  - Bell operators as two coupled modes;
  - Landau's single violating direction;
  - validity as an effective resistance;
  - AM–GM as a force balance;
  - one identity that certifies everything except a single mixed multiplier.
- **Part C, the explicit lift.** The sextic lift and the progress toward an
  explicit pencil.

## Verdict

- **Existence proof: sound wherever it was tested.** Every core lemma held under
  adversarial review. The thin-R10 bounded-matrix implication is being
  verified separately (see the existence-proof check above). The key determinant identity was re-derived symbolically for
  all weights. The audit of the ~31,000 appendix lines found editorial defects but
  no mathematical gap.
- **Sketched steps.** The two load-bearing steps that are only sketched in the text
  were tested separately, and all of them hold. These are the thin-layer reuse of
  C1 certificates on a negative-β source, R10 coverage of the Z = 1 contact chart,
  and the HI radical-division pivots (`c2_short.tex:6849–6909`).
- **Explicit lift: reduced to one canonical cone and advancing.** The repo
  reduces the problem exactly to one 9-dimensional canonical cone. The review
  branch reduces that cone to the strict C1/C2 contact families, covers every
  compact part of them, solves the `c`-end limit problems in closed form, and
  reduces `r → 0` and the walls to limit problems. The exact pencil for all
  eight coordinates is the construction target. Exact counterexamples show
  which fixed selector families are too small: the full-set pencils of order 264 and 392 are strictly larger than
  Q222, which is why the sixth-root features are used.

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

The exact canonical pencil is the piece under construction. Exact rational
sources show which simpler fixed families are too small:

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
  the native C2 source envelope. Global canonical coverage is the target of
  the review route.

**The symmetric repair is also too small.** The order-98 cone `C_ram + πC_ram`
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
u^{i/6}v^{j/6}ψ, with no charts. With half-powers (k = 2) it admits 293 of
the 302 test sources. With sixth roots (k = 6) no tested source falls
below 0.999 of the validity limit, including every source that defeats the
earlier candidates.

Part of exactness is now proved:

- the lift is sound and closed;
- exactness reduces to the extreme sources;
- four genuine C1/C2 boundary sources have exact rational certificates;
- a continuation theorem reduces the C1/C2 families to a single
  nondegeneracy property, (ND).

The two closed faces use their exact pencils, so the decisive step is
admission of the strict contact families, through (ND) or the dual criterion
of the route.

The [assembly refinement](ASSEMBLY_REFINEMENT.md) makes this target more
precise. Continuation takes place in the 154-function support-face basis,
after 42 forced null directions have been removed. Its coefficient rank
is exactly 406 throughout the native contact charts. The full candidate
also has a rational strictly feasible point, with all event weights one
and CHSH weight `1/21`. Both statements now have exact reproduction code.
For a hybrid lift, the two closed faces can use their existing exact
pencils, so their membership in the sixth-root cone is unnecessary.
The next step to prove is global admission of the contact sources; the
nondegeneracy condition (ND) is a sufficient form of it, concerning positivity
of the reduced Grams.

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
| `gram-lift/strict_feasibility.py` | exact positive definite full196 Gram at `c*=(1/21,1,...,1)`; source signs and coefficient identity | `python strict_feasibility.py` | rational witness, no solver |
| `gram-lift/contact_rank.py` | six independent contact conditions and matching modular rank; both native support faces | `python contact_rank.py` | rank406 for both families |
| `gram-lift/ends/` | the `a → 0` end of C2: the two Gram blocks fixed by single monomials, their determinants, and their exact order `a²` at three interior points (so no diagonal rescaling gives a uniform margin there) | `python check_a0_corners.py` | exact identities pass; limits match the SDP survey |
| `gram-lift/ends/check_c2_ends.py` | the other C2 ends: `n(r=1) < 0`, so the chart stops at an `n = 0` wall; the orders as `r → 0`; the torus-orbit identity behind the `c → 0` and `c → ∞` ends, for every C2 source | `python check_c2_ends.py` | exact identities pass |
| `gram-lift/ends/check_c2_r0_face_k2.py` | the `r → 0` limit form of C2: no integer-power (`k = 1`) Gram exists at any face point, and exact `k = 2` Grams positive on each contact-kernel complement at 9 rational face points | `python check_c2_r0_face_k2.py` | all assertions pass (a few seconds) |
| `gram-lift/ends/check_c2_n0_wall_k2.py` | the C2 wall `n = 0` as the explicit graph `t = t_w(c,a,b,r)`, with exact rational `k = 2` Grams (one contact kernel per parity block) at 8 wall points | `python check_c2_n0_wall_k2.py` | all assertions pass (about 20 s) |
| `gram-lift/ends/check_c1_zero_walls_k2.py` | the C1 walls `β = 0` and `D = 0` as explicit graphs, with exact rational `k = 2` support-face Grams (one contact kernel per parity block) at 17 wall points | `python check_c1_zero_walls_k2.py` | all assertions pass (about 30 s) |
| `gram-lift/ends/n0-wall/` | a generic rational 16-parameter family of 18-feature C2 `n = 0` Grams satisfying the identity and contact kernels; exact PSD members positive on the three contact-kernel complements at 22 rational wall points; exceptional specializations and universal PSD remain open | `python check_c2_n0_wall_family.py` | all assertions pass (about 80 s) |
| `gram-lift/ends/check_c2_c_limits.py` | the `c → ∞` and `c → 0` limit forms of C2 as explicit sums of squares, and a one-parameter `k = 2` Gram whose kernel is exactly the contact on the whole open limit face (explicit LDL, all pivots positive); the finite-`c` perturbation step is not included | `python check_c2_c_limits.py` | all assertions pass (about 90 s) |
| `gram-lift/admission/` | exact pointwise proof of strict sixth-root admission `μσ > 0` at rational C1/C2 sources: double-precision SDP, exact rational projection onto all 412 coefficient equations, 400-bit ball Cholesky of the four blocks; with a negative control | `python run_points.py` (needs python-flint) | hardest known sources certified; control rejected |
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
