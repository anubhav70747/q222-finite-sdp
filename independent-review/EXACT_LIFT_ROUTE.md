# Use NPA, then certify the missing boundary

The construction should retain the inequalities that a small NPA relaxation
already proves and add certificates for the boundary it misses. The new
mathematics belongs in that second step. The final lift joins the corresponding
finite PSD blocks on the same eight correlation coordinates.

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

For finitely many closed cones in the nonnegative canonical source orthant,
the sum is closed. Indeed, if their sums converge, each summand is bounded
coordinatewise by the total. Subsequence limits remain in the respective
cones and give a decomposition of the limit. Thus closure does not remove
the extreme-ray requirement in this setting. This argument uses canonical
nonnegative coordinates; it is not an assertion that arbitrary projected
cone sums are closed.

There is a related test before any proposed decomposition: if a valid slack
vanishes at a strategy, every valid summand vanishes there too. Their values
are nonnegative and sum to zero. A generic positive NPA certificate cannot
therefore be subtracted from a tight inequality while preserving validity.

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

For the sixth-root route, prove `μσ(c)≥0` for every required C1/C2 source.
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
