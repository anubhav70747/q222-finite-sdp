# Use NPA, then certify the missing boundary

The construction should retain the inequalities that a small NPA relaxation
already proves and add certificates for the boundary it misses. The new
mathematics belongs in that second step. The final lift joins the corresponding
finite PSD blocks on the same eight correlation coordinates.

The [equal-angle C2 cone](../extracted-lifts/c2-equal-angle-cone/) now covers
two independent contact parameters with an exact order-seven block.
Its source matrices are `[[m,A−m],[A−m,C]]` and `[[D,A−D],[A−D,α]]`,
with `A−m,A−D,m+D−A≥0`, `d=4A−2m−2D`, `n=m+D−A` and `λ=β=0`.
Every feasible source decomposes into one tight source and the nonnegative
excesses of `C` and `α`. Four AM–GM terms prove validity at every test.

This cone includes the [quadratic curve](../extracted-lifts/c2-quadratic-family/),
whose sources lie outside `1+AB`. The order-seven block therefore strictly
sharpens `1+AB`; there is no need to append the contained order-three curve
block separately. Both constructions avoid the large state-polynomial Gram.

The [unequal-angle contact equations](C2_CONTACT_GEOMETRY.md) now give a
symmetric three-parameter family and a fourth-degree identity that fixes a
direction missed by its quadratic relations. An
[exact order-six cone](../extracted-lifts/c2-unequal-angle-curve/) covers one
compact unequal-angle curve. Its proof reconstructs every feasible moment
vector from at most three contact points. The convex hull with all three
parameters varying remains unresolved.

## A concrete starting relaxation

Fix the nine words

`W = {I,A0,A1,B0,B1,A0B0,A0B1,A1B0,A1B1}`.

Their moment matrix gives the NPA `1+AB` relaxation. Here
`Ax²=By²=I` and observables of different parties commute. Event probabilities
have certificates within this list: the event projector
`P=(I+aAx)(I+bBy)/4` satisfies `P=P†P` and belongs to `span W`.

For the conventional CHSH expression `S=A0(B0+B1)+A1(B0−B1)`, direct expansion gives

`2√2 I−S = [(A0−(B0+B1)/√2)²+(A1−(B0−B1)/√2)²]/√2`.

The tilted-CHSH slacks also have exact SOS certificates in this same
nine-word space for `0≤α<2`; see [Bamps and Pironio, Sec. II, Eqs. (6), (27)–(28)](https://arxiv.org/html/1504.06960v1).
These are concrete exact families. They do not classify which of the full
C1/C2 source families this relaxation certifies. The relevant uniformity is
one word list for the whole construction, as opposed to convergence through
successively larger lists; see [Navascués, Pironio and Acín](https://arxiv.org/abs/0803.4290).

| Sources or tests | Available certificate | Remaining task |
|---|---|---|
| Positive event rays | Event-projector squares in `1+AB` | None |
| CHSH and one-sided tilted CHSH | Exact SOS in `1+AB` | Retain globally; these familiar slacks cannot contribute to the native contact sources below |
| The two canonical closed faces | SG12: one order-6 Gram; SG13: two order-3 Grams | Retain their finite relabeling orbits |
| Equal-angle C2 contacts, two independent parameters | Exact order-7 source and behavior pencils, including the NPA-hard quadratic curve | Retain finite physical relabeling orbits; the curve block is already included |
| A compact unequal-angle C2 curve | Exact order-6 source and behavior pencils, including all interval endpoints | Extend from this curve to independent contact parameters |
| Open C1/C2 contact sources | Four exact sixth-root certificates; constant coefficient rank | Prove admission of every source left by the retained cones |
| Zero-weight limits of admitted contacts | Closedness of the fixed sixth-root cone | Establish admission before taking the limit |
| The R1 contact cap and Bob complement | Exact regional certificates | Cover the remaining tests before using this as a global cone |

The closed-face pencils can be added directly. They need not be represented
again in the sixth-root basis. Their derivations and dual attainment are in
[`04_geometry.tex`, SG12–SG14](../latex/chapters/04_geometry.tex).

## Two classical zeros reduce the NPA test to seven roots

Every native C1 source has `λ=n=0`, and every native C2 source has
`λ=β=0`. The repository's canonical orientation is
`H=I−(A0B0−A0B1+A1B0+A1B1)/2`, a relabeling of CHSH.
In the canonical expansion `f=dH+4Σ ai pi`, each event is one
at its own deterministic CHSH-facet neighbor and zero at the other seven.
The omitted coefficients therefore give two exact zeros of the slack.
In the order `(A0,A1,B0,B1)`, they are

`vλ=(1,1,1,1), vn=(1,−1,−1,−1), vβ=(−1,1,−1,1)`.

This immediately limits what the familiar NPA certificates can contribute.
A tilted-CHSH slack with `0≤α<2` is strictly positive at every deterministic
strategy, since

`√(8+2α²)−(2+α) = (2−α)²/[√(8+2α²)+2+α] > 0`.

Every valid summand of a slack with a classical zero must vanish there.
Consequently, neither the ordinary Tsirelson slack nor a subcritical
tilted-CHSH slack can enter a valid decomposition of these native sources
with positive weight. This applies to their relabelings too. Any further
native coverage by `1+AB` requires different SOS certificates.

The same two zeros make that search smaller. Each square root in a nine-word
SOS must vanish at both deterministic assignments. Solving these two linear
conditions gives the following seven-root bases. Each root has its own
nonpivot word with coefficient one, so the seven roots are independent:

```text
C1: A0−I, A1B0−I, A1B1−I,
    B0−A1, B1−A1, A0B0−A1, A0B1−A1.

C2: A1−I, B1−I, A0B0−I, A1B1−I,
    B0−A0, A0B1−A0, A1B0−A0.
```

Write the appropriate column as `Rσ`. A native source belongs to the
`1+AB` certificate cone exactly when

`f_c = Rσ† Q Rσ`, with `Q⪰0`,

as an identity in the universal observable algebra. To prove the converse
direction, a PSD nine-word Gram annihilates both deterministic evaluation
vectors because its quadratic form is zero on each. Its range lies in
their common orthogonal complement, whose basis is exactly the displayed
list. Compressing the Gram to this basis gives the seven-by-seven matrix.
All coefficient identities remain imposed. This reduces the membership
problem; it does not assert that every native source passes it.

The coefficient equations are explicit. In the order of the nine words
above, the target coefficients are

```text
I:     d+λ+n+m+C+A+D+α+β
A0:    λ+n−A−α          A1:  m−C−D+β
B0:   −n+C+A−β         B1:  λ−m+D−α
A0B0: −d/2−n−A         A0B1: d/2+λ+α
A1B0: −d/2−C−β         A1B1:−d/2−m−D.
```

Expand `Rσ†QRσ`, commute Alice's words past Bob's, and cancel adjacent
equal observables. Match the nine coefficients above and set every other
word coefficient to zero. This specifies the finite linear equations
without evaluating the operators at a particular strategy.

## Why an extreme inequality must be covered directly

Let `C` be a pointed cone and `c` generate an extreme ray. If
`c=c1+...+cr` with every `ci∈C`, then each nonzero `ci` is a nonnegative
multiple of `c`. This is the definition of an extreme ray, applied to each
summand and the sum of the others. Consequently, a sum of sound certificate
cones covers this ray only if at least one summand cone covers the ray.

Every extreme canonical ray with `d>0` is also extreme in the full cone
`P` of quantum-valid affine slacks. Indeed, the eight canonical event
coefficients of any valid slack are nonnegative: each is one quarter of
its value at the corresponding deterministic facet vertex. Since
`d=−f(PR)`, the canonical cone is exactly `P ∩ {d≥0}`.

Suppose `f=g+h`, with `g,h∈P` and `d(f)>0`. For sufficiently small
`ε>0`, both `f±=(1±ε)g+(1∓ε)h` remain valid and have positive `d`.
They lie in the canonical cone and have average `f`. Canonical extremality
forces them, and hence `g,h`, to be proportional to `f`. This applies to
all positive-core C1/C2 extreme rays, including zero-event faces, and
allows the other cones to use arbitrary Bell coefficients. It uses no
strict inequality at a particular quantum contact.

The sum of finitely many closed, globally sound certificate cones is closed,
even outside one nonnegative canonical orthant. Write a slack as
`f=t+b·q`. Averaging over the 16 deterministic assignments gives
`|b_i|≤t`: its constant is the average of its nonnegative values, and each
other coefficient is their signed average. If a sequence of sums converges,
each summand has a bounded constant and hence bounded coefficients.
Subsequence limits give the required decomposition. Thus closure does not
remove the extreme-ray requirement. The
[assembly proof](ASSEMBLY_REFINEMENT.md) also proves closedness of the fixed
NPA source cone directly from `tr G=t`.

There is a related test before any proposed decomposition: if a valid slack
vanishes at a strategy, every valid summand vanishes there too. Their values
are nonnegative and sum to zero. A generic positive NPA certificate cannot
therefore be subtracted from a tight inequality while preserving validity.

## Only the strict contact families need new admission proofs

The finite boundary cases follow from the open families once the certificate
cones are fixed and closed. This removes the need to reconstruct a new lift
at each third omission.

Here is the source reduction behind that statement. The cleared determinant
has the form

`det P = 16 D_u² D_v² u²v² L(u²,v²)`,

`L(U,V)=h+Σ_e q_e U^{e₁}V^{e₂}`, with `q_e≥0` and
`e∈{−1,0,1}²\{0}`. The coefficients and the positive three-state block are
given in [SG3–SG6](../latex/chapters/04_geometry.tex).
If zero lies inside the active exponent polygon, the logarithmic Laurent
sum is coercive. An extreme source then has a finite zero: a strictly positive
minimum permits small opposite source perturbations while preserving all
coefficients' signs and global positivity.

At a zero, the four state-kernel equations and two angular derivatives are
six homogeneous linear source equations. More than six positive event
weights would leave a nonradial perturbation preserving them. Such a
perturbation remains globally valid: stationarity gives
`Σ_e q_e U_c^{e₁}V_c^{e₂} e=0`, and weighted AM–GM makes `L` nonnegative
at every independent test. Hence an extreme contact has at most six
positive event weights.

The lower bound comes from the eight event cycle
`(λ,D,C,A,α,m,β,n)`. The product formulas in SG6 show that a noncoercive
support omits a pair of cyclic distance three or four. Physical relabeling
puts it on one of the two exact closed faces. Otherwise omitted events are
pairwise compatible, extend to a deterministic assignment, and number at
most three. The surviving extreme contacts therefore have five or six
positive event weights: precisely C1 or C2 and their third-omission faces.

To open a third omission, fix the other event weights and replace the missing
weight by `ε>0`. Normalize the contact state by `Z=1`, put
`z=(W,X,Y,u,v)`, and solve `(∇_zF,F)=0` for `(z,d)`. At the original
contact the Jacobian is

```
[ Hess_z F    ∂_d∇_z F ]
[     0         F_H    ].
```

The contact is nondegenerate: `Hess_z F≻0` by the positive state block and
the angular Hessian `Σ_e q_e U_c^{e₁}V_c^{e₂} eeᵀ≻0`.
The core matrix is nonsingular at noncommuting angles. If every positive event
row vanished on the contact state, the full kernel equation would force the
core to annihilate that state. Hence some event contributes, and `F_H<0`
follows from `F=0` and `d>0`. These facts are proved in
SG17–SG21. The Jacobian is invertible, so the
implicit-function theorem gives an analytic branch of strict contacts with
`d>0` and positive angles. Stationarity proves their global validity as
above. This opens `D=0` or `β=0` in C1, or `n=0` in C2, while retaining
the family's two prescribed omissions and an actual contact throughout.

Consequently, a fixed finite union of closed, globally sound certificate
cones covering **all strict six-event C1 and C2 contacts** covers these
finite boundary contacts as well. Take the limit in one cone along a
subsequence. Together with the event and closed-face cones, the extreme-source
classification then supplies every canonical source.

Escaping contacts also have a short characterization. For convergent normalized
sources with `d` bounded away from zero and coercive limiting support, there
are `δ,a>0` such that

`Σ_e q_e exp(e·t) ≥ δ exp(a‖t‖)`.

At a contact this sum equals `−h`, which stays bounded. Thus both angles
stay bounded above and away from zero; the positive state block also bounds
the state with `Z=1`. An escaping sequence must approach an event combination
or one of the noncoercive closed faces.

This reduction needs pointwise admission throughout the two strict families
in **one fixed finite collection of cones**. It requires no uniform positive
Gram margin as a source approaches a boundary. A separate finite cone chosen
for each source would not suffice.

## One fixed dual problem for each contact family

At an open native contact, first remove the 42 features excluded by the
source support, then the contact vector in each of the four parity blocks.
Rescale the nonzero contact coordinates to one. Monomial scaling is invertible
and gives one fixed coefficient map for each family:

`Aσ : S^51 ⊕ S^35 ⊕ S^35 ⊕ S^29 → Vσ`, with `dim Vσ=406`.

There are 3,021 symmetric matrix entries. The previously established rank
406 makes this map surjective. Write `bσ(c)` for the coefficients of the
complete parent form after this normalization. Normalization is used to
compare certificates; the proposed lift still uses its original fixed
testing functions.

Define

`Dσ = { y∈Vσ* : Aσ*y ⪰ 0, tr(Aσ*y)=1 }`,

`μσ(c) = min_{y∈Dσ} <y,bσ(c)>`.

Then `μσ(c)≥0` is equivalent to sixth-root admission, and `μσ(c)>0`
is equivalent to a positive definite Gram on the contact complements.
A negative value supplies an exact dual separator.

Here is the complete finite-dimensional argument. Integration over a fixed
positive testing rectangle and the state sphere gives a positive definite
moment matrix for each independent reduced feature list. Hence bounded
coefficient vectors have bounded PSD Gram preimages. The image of the PSD
cone is therefore closed. Its dual is `{y:Aσ*y⪰0}`. Surjectivity of `Aσ`
makes its adjoint injective, so the trace-one section `Dσ` is compact.
It is nonempty: the same integration functional has a positive definite
adjoint. Separation proves the nonnegative criterion. A surjective linear
map sends the interior of a convex cone onto the interior of its image;
this proves the strict criterion.

Equivalently,

`μσ(c) = max { ε : Aσ(G)=bσ(c), G⪰εI }`.

Any symmetric solution of the linear equations becomes strictly feasible
after choosing a sufficiently negative `ε`. The compact dual section bounds
the objective. The moment bound above gives attainment of the primal maximum.
The dual formula also makes `μσ` a continuous, concave function of the
normalized coefficient vector. It need not be concave in the native source
parameters, since the normalization depends on their contact.

The matrices in this test are the **compressed** adjoints on the contact
complements. Positivity of an uncompressed moment matrix is a stronger
restriction on the dual variable and can miss an obstruction. Physical
evaluation functionals alone likewise do not exhaust this dual cone.

## The shortest remaining target

Two explicit constructions now isolate the state-matrix step. At every native
C2 contact, a [uniform identity](C2_MULTIPLIED_CERTIFICATE.md) writes
`det(H3) F` as at most 183 squares from 676 fixed third-root features.
The positive multiplier depends on both source and testing variables, so
it cannot be cancelled to obtain the affine lift. Independently, the
[contact scalar](CONTACT_SCALAR_GRAM.md) has at most three squares in the
fixed span `{1,u,v,uv}`, without a multiplier. Its three-by-three Gram is
positive definite at every finite native C1/C2 contact, including the
permitted zero-event faces, so its source-only inverse is available.
The three couplings to the remaining state coordinates still require
a uniform completion.

The C1 family also needs additional certificates. Two globally valid
[rational C1 extreme sources](C1_NPA_OBSTRUCTION.md), including one with
`beta=0`, violate a necessary linear inequality for the seven-root NPA
Gram. Its integer dual matrix is positive definite; both global validity
proofs reduce to AM--GM. A further rational contact has an invertible
six-word coefficient matrix. Rationality of the contact action extends
this exclusion to an open dense subset of the strict five-dimensional
C1 contact parameter domain. Thus `1+AB` alone cannot supply the C1 part
of the construction. The selected exceptional algebraic sets remain unclassified.

For maximally entangled C2 sources, a change of contact coordinates gives
an [exact order-four source pencil at every fixed pair of parameters](C2_FIXED_FIBERS.md).
Its nontrivial block is only two-by-two. Every source in that cone is a
sum of at most two contact rays, and the same formulas give a unique
inverse and branch condition for individual generators. Combining the
continuum of these cones into one fixed affine lift remains unresolved.

The need for additional certificates holds on the entire finite maximally
entangled C2 native family: [six forbidden-word equations](C2_MAXIMAL_NPA_OBSTRUCTION.md)
have a strictly positive determinant throughout its contact domain. This
excludes `1+AB` also at unequal angles and on `n=0`, without an exceptional
interior locus. Singular contact limits require separate analysis.

An explicit partially entangled extreme source also lies outside `1+AB`.
Its global validity follows from one weighted AM--GM determinant identity;
six forbidden-word coefficients give an invertible matrix of determinant
`-405`, excluding every nine-word Gram. The same obstruction holds on a
relative Zariski-open dense subset of the strict native C2 interior.
The [complete rational proof](C2_NONMAXIMAL_NPA_OBSTRUCTION.md) includes
the contact, all coefficients and a full-polar extremality argument.
Thus the additional cones must cover generic changes in Schmidt coefficients,
not only the maximally entangled surface. Exceptional parameter loci are
not classified by this generic argument.

The direct source-cone construction now covers a whole two-parameter
surface of extreme rays. The next extension must retain further independent
contact parameters, including unequal measurement angles and changes in
Schmidt coefficients. The order-seven cone does not yet cover those
general C1/C2 sources.

For the sixth-root route, prove `μσ(c)≥0` for every still-required C1/C2 source.
Strict positivity would give a convenient continuation proof, but it is not
required for exactness. Sources with zero margin may have valid singular
certificates. Alternatively, retain another finite, globally sound certificate
cone that covers any rays the sixth-root cone misses.

Once every remaining extreme ray and its limits is covered, the
[assembly lemma](ASSEMBLY_REFINEMENT.md) gives the full lift by intersecting
the dual constraints. The eventual short proof can then follow four steps:
reduce to the two contact families, prove their admission, include their
limits and the closed faces, and join the finite PSD blocks. Universal
contact admission is still the unresolved step; the dual reformulation
does not establish its sign.
