# Exact certificates for the remaining sources

A shorter proof can retain the bounds already obtained exactly and concentrate
on the Bell inequalities they leave unresolved. The additional construction
then has a precise job: supply the missing supporting inequalities. This
also allows the existing closed-face pencils to be used directly, without
forcing them into the sixth-root Gram basis.

## How the pieces combine

Let `P` be the cone of affine functions nonnegative on Q222. Let `N_k`
be a cone of sound certificates from one fixed finite NPA level, and let
`H_1,...,H_r` be further globally sound certificate cones. Suppose their
duals have finite SDP representations, with attainment established by
strict feasibility in the appropriate face or an explicit facial reduction.

**Assembly lemma.** If

`P = closure(N_k + H_1 + ... + H_r)`,

then

`Q222 = {q : (1,q) ∈ N_k* ∩ H_1* ∩ ... ∩ H_r*}`.

**Proof.** A vector is nonnegative on a sum of cones precisely when it is
nonnegative on every summand. Taking a closure leaves the dual unchanged.
The right side therefore consists of the points satisfying every quantum-valid
affine inequality, and equals Q222 by separation. Join the finite SDP blocks
using the same eight coordinates `q`; their auxiliary variables remain
independent. This gives one finite lift. ∎

It suffices to prove the cone equality on all extreme rays of a compact
normalized base, including their limiting faces. Every additional source
must be certified at **all** independent measurement angles and states.
The order-293 cap/Bob construction still leaves a testing strip, so it cannot
yet be inserted as an unrestricted globally sound hard-source cone.

NPA's good numerical performance motivates this division of work. Exact
coverage requires one fixed word list and an analytical certificate on the
claimed family. A different level for each source supplies no uniform finite
lift. The distinction is already present in the original
[NPA convergence and finite-certification theorem](https://arxiv.org/abs/0803.4290).
Exact SOS bounds are known for useful families, including
[tilted CHSH](https://arxiv.org/abs/1504.06960).

For the canonical cone, two parts of this assembly are available. SG12 in
`latex/chapters/04_geometry.tex` represents the face `λ=m=0` by a six-row
PSD Gram; SG13 represents `λ=α=0` by two three-row PSD blocks. Each also
has seven remaining nonnegative source coordinates. Their strict points and
the exact duality argument SG14 are printed there. Add their finite physical
relabeling orbits as separate cones. If the C1/C2 contact rays and their
relabelings are admitted by the sixth-root construction, these face cones and
the event rays complete the canonical cone. Proving that the closed faces
themselves belong to `L_6` is unnecessary for this hybrid theorem.

## A rational strictly feasible point

The sixth-root candidate has a strictly feasible point with particularly
simple source coefficients:

`c* = (d,λ,n,m,C,A,D,α,β) = (1/21,1,1,1,1,1,1,1,1)`.

Here is an explicit construction. In the 16-function corner basis
`z_1={1,u,v,uv}⊗(W,X,Y,Z)`, expand the eight unit-weight event terms by
`D_u=1²+u²` and `D_v=1²+v²`. Their 18 square roots give a rational Gram
`Q_0`. They span the corner basis: the `n,A` terms give all `X,Y`
monomials; the other six events successively give `W,uW,uvW,Z,uZ,uvZ,vW,vZ`.
Thus `Q_0` is positive definite.

The following two matrices make its eigenvalue bound easy to reproduce:

```text
B1 = [ 1  0  1  0 ]       B2 = [ 3  2 -1 -1 ]
     [ 0  1  0  0 ]            [ 2  3 -1 -2 ]
     [ 1  0  3  1 ]            [-1 -1  1  1 ]
     [ 0  0  1  1 ]            [-1 -2  1  3 ]
```

After grouping parity classes, `Q_0` is four times the direct sum of
`B1,B2` and sign-conjugate copies of `B2,B1`. Direct inversion gives
`tr(B1^-1)=6`, `tr(B2^-1)=9/2`, hence `tr(Q_0^-1)=21/4` and
`Q_0 ⪰ (4/21)I`. The signed CHSH Gram `Q_d` satisfies `Q_d ⪰ -2I`:
its only negative terms are `-2Z²,-2u²Z²,-2v²Z²,-2u²v²Z²`.
Consequently

`Q_c := Q_0 + Q_d/21 - I/21 ⪰ I/21`.

To fill the sixth-root basis, set `r=u^(1/6)` and let `H` be the symmetric
seven-by-seven matrix with diagonal `(1,1/2,1/2,1/2,1/2,1/2,1)` and
entries `H[i-1,i+1]=H[i+1,i-1]=-1/4` for `i=1,...,5` (zero-based indices).
All other entries vanish. Each interior coefficient cancels, giving

`(1,r,...,r^6) H (1,r,...,r^6)^T = 1+r^12`.

The even and odd blocks of `H` are irreducibly diagonally dominant, with
positive diagonal and a strict row in each block, so `H ≻ 0`. Let `J`
embed the corner monomials into the full 196-function basis. Then

`G* = J Q_c J^T + (H⊗H⊗I_4)/21 ≻ 0`

and `Z_6^T G* Z_6 = F_{c*}` exactly. All nine source coordinates are
strictly positive. This proves the Slater condition for the full candidate,
independently of its conjectured exactness.

[`strict_feasibility.py`](checks/gram-lift/strict_feasibility.py) verifies the
identity in rational arithmetic and generates the complete
[196-row Gram witness](checks/gram-lift/strict_feasibility.json), with 692
nonzero upper-triangle entries. It uses no SDP solver.

## Continuation on the correct source face

The four exact contact certificates use 154 functions. The missing 42
functions are forced out by the source support before the contact kernel is
removed. This reduction matters: when `λ=0`, the constant `W²` coefficient
is zero, so positivity forces the full Gram row of `W` to vanish. Removing
only a contact vector from the full 196-function matrix cannot give a
positive definite matrix.

The appropriate bases follow directly from the Newton polygons of the
four diagonal state coefficients. Write a feature as `u^(i/6)v^(j/6)ψ_l`,
with `0≤i,j≤6`. In both open C1 and C2 families keep `i+j≥6` for `W`.
For C1 also keep `j≥i` for `X`; the `Y,Z` grids are complete. For C2
the `X,Y` grids are complete and keep `i+j≤6` for `Z`. Each basis has
parity blocks of orders `(52,36,36,30)`.

After removing one contact vector in each block, the coefficient map has
**rank 406 at every contact in these native charts**. To see why the rank
does not change with the contact, rescale `u,v,W,X,Y,Z` by their nonzero
contact values. This is an invertible diagonal change of both the feature
basis and the coefficient monomials, and sends the contact to all ones.
The two angle derivatives and four state derivatives give six independent
linear conditions on the 412 coefficients. Every product of two functions
vanishing at contact satisfies these conditions. Thus the rank is at most
406. Exact integer elimination modulo the prime 65521 gives rank 406 for
both source faces, proving the matching lower bound over the rationals.
[`contact_rank.py`](checks/gram-lift/contact_rank.py) constructs both maps,
checks the six annihilators exactly and records the
[nonzero modular pivots](checks/gram-lift/contact_rank.json).

The native contact formulas have nonzero state coordinates in the open
families. Zero openings belong to the limiting faces and must be treated
by limits or their own certificates.

Constant rank proves local continuation from each positive definite reduced
Gram. It does not prove that this positivity persists through the entire
connected source family: a sequence of positive definite Grams may converge
to a singular Gram while the coefficient map retains its rank. The remaining
nondegeneracy assertion must exclude that loss of positivity, or supply a
different certificate at such a point. Closedness of `L_6` alone does not
do so.

This isolates the proof target. The known face pencils and the strictly
feasible assembly can be retained. The new analytical work is the admission
of every remaining contact source, possibly after using one finite NPA
level to remove families it already certifies exactly. The short guide's
imported hard lemma and the numerical sweep do not yet provide that
admission theorem.
