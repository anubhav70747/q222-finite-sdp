# R6.5--R6.7 endpoint dependency repair

Repository inspected at `d01e5e089c5b4ad73667e922eec5b62b5d82485d`.
This work changes only this repair directory. No manuscript, archived source,
status document, or existing proof script was edited.

## Verdict

The missing endpoint argument is recoverable from the actual four-state parent.
The original R6.5--R6.7 identities are correct with their stated coefficient
hypotheses. Their use at `n0=0` is valid except that R6.5 is nonstrict when the
other testing angle also has zero residue. That corner needs the independent
`(infinity,0)` four-root allocation. The proposed TeX supplies its complete
polynomial, not just a determinant or a state restriction.

An additional simplification proves all these endpoint joins for both R6 and
R7: cap each of the actual `m,C,A,D` coefficients at a positive real lower
bound, retain all four excess event squares and the whole actual `n,alpha`
events, and apply the existing strict bounded-Gram lemma to the complete capped
homogeneous polynomial. This does not need bounded actual `A` or `D`, or the
centered extraction. Its four-/five-root residues are exactly the RU28--RU30
algebra with `C` also capped; all are strict except the same `(infinity,0)`
corner, which the explicit allocation closes. These claims concern the
projective endpoint joins, not the local contact conductor or global theorem.

## Actual dependencies read

- Current master parent: `latex/certificates/c2_short.tex:11-25`.
- Original R6.5--R6.7 and its narrower bounded-a source theorem:
  `complete-proof-materials/materials/historical/recovery/full_r21/components/source/pd5fa361b4623.tex`.
- Exact centered event identity and homogeneous core:
  `.../source/p8da9ea6ed9e5.tex`, R5.2--R5.5 and R5.16.
- Original exceptional corner root frame:
  `.../source/p69995d829bb2.tex:86-99`.
- Actual R6 unit-b source bounds and its consumers:
  `.../source/pb60eba754641.tex`, RB1--RB5 and endpoint paragraphs.
- Actual R6 infinitesimal-b source and consumers:
  `.../source/p75b395956985.tex`, source scope and endpoint paragraph.
- Thin-layer source and consumers:
  `.../source/pc82cd677f4b9.tex`, RT1--RT3 and endpoint paragraph.
- Actual R7 source bounds and own endpoint formulas:
  `.../source/p405d96527133.tex`, RU1--RU3 and RU28--RU30.
- The fetched coefficient proof, Sections 1 and 4, as a secondary comparison:
  `complete-proof-materials/materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-for-publication/R6_R7_COEFFICIENT_PROOF.md`.
- Existing tensor-Gram proof: `latex/chapters/03_tensor.tex:77-115`, label
  `p25af2c58da7e:lem:gram`. It allows unbounded or dependent root entries,
  provided the *matrix entries* lie in the bounded separated ring and the
  entrywise real residue is positive definite.

Current consumers at `c2_short.tex:12412-12438` and `12887-12896` refer to the
original R6 formulas, but their target near line 14301 is only a label.
The R7 formulas are displayed at approximately `13400-13425`; their exceptional
corner still imports the BE frame. The proposed subsection defines the two
existing `pd5fa361b4623` aliases and supplies the common algebra explicitly.

## Source-domain transfer

| Source intersection | Coefficients relevant to endpoint transfer |
| --- | --- |
| Original bounded-a R6 | `n0>0`, `m` unbounded, `mg` bounded; `C,A,D` positive bounded units and `alpha0=0`. |
| Unbounded-a R6, unit b | `n=O(s^2)` so `n0=0`; literal `n=0` allowed. RB4 gives positive real lower bounds on `m,C,A,D`. |
| Unbounded-a R6, infinitesimal b | Actual `n>0` but still `n0=0`. The same RB4 bounds hold without a lower bound on b. |
| R6 thin layer | RT2--RT3 give `n0=0` and the same positive lower bounds. |
| Unbounded-a R7 | Actual `m,A` unbounded, `C0=t0/2>0`; when `t0>1`, `D0=1/[2(t0-1)]>0`; when `t0=1`, actual D is unbounded. All four positive real caps exist. |

The R6 formulas are not imported together with the narrower bounded-a theorem.
Their polynomial derivation is independent of that source theorem. In
particular, it is legitimate to set the bounded coefficient residue `n0=0`
without taking a limit of a theorem that excludes literal `n=0`.

The R7 boundaries `h0=0` and `j_gap,0=0` are handled by caps of the actual
unbounded coefficients. Actual h and the native gap remain strictly positive;
the repair does not claim a source theorem on their literally zero faces.
Likewise actual `D=0` is not within these R6/R7 source hypotheses.

## Endpoint checks and strictness

1. **u unbounded, v bounded with nonzero residue.** Exact `p^2 P` has the full
   four-state support `(X,W,Y,Z)`. The centered correction has a factor p and
   vanishes in residue. The R6.5 kernel determinant at `n0=0` is
   `A0*C0*(D0+1/2)*v0^2*(1+v0^2)^2>0`. The full residue adds the positive
   `(K+1/2)(vY-Z)^2` row. The capped proof gives the same strict four-row rank.
2. **v unbounded, u bounded (including zero).** In exact `q^2 P`, the centered
   correction does not vanish identically: its residue is
   `2*kappa*W*(W+uY)`. It does vanish on the limiting m-row kernel. R6.6 is
   strict there, with determinant
   `D0*(n0+1/2)*(1+u0^2)^2*[C0*(1+u0^2)+A0+1/2]>0`.
   A real cap on the transverse squared row makes the full residue strict.
   Alternatively the capped whole-state residue is manifestly a sum of five
   positive square rows with common kernel zero, for every real u0.
3. **Both angles unbounded.** Exact `p^2 q^2 P` is a polynomial on the five
   roots `(X,W,Y,pZ,qZ)`; its canceled core is printed. Centered corrections
   have factor p and vanish in residue. R6.7 has block determinants
   `D0*(n0+1/2)` and `A0*(K+1/2)`, positive at `n0=0`. The capped proof is
   equally strict. The root entries need not be independent.
4. **u unbounded, v infinitesimal, including the literal endpoint.** The
   ordinary n0=0 residue has an X kernel. After retaining exact positive event
   packets, the full selected polynomial has the exact frame
   `(H=pX+Z,W,Y,T=vX)` with polynomial coefficients in p,v. Its residue is
   `H^2/2+(W+T)^2/2+Cc*W^2+Ac*Y^2`, determinant `Ac*Cc/4>0`.
   No inverse of p or v is used to establish support.

All polynomial identities hold for either sign of the finite testing angle.
Literal reciprocal zeros are values of the polynomial charts. For nonzero p
or q, homogeneity is restored by a pure-testing square only. A real state is
unrestricted; for a complex state the two real component identities add.

## Tensor-SOS step and limits of this audit

Real determinant positivity alone would not close the tensor-SOS claim.
Here the exact root identities supply the *complete* finite matrix, including
every state cross term. In the centered proof its entries are bounded
separated tensors because `(m-K)g` is bounded. In the capped endpoint proof all
entries are even bounded pure-testing elements: the only surviving event
coefficients are positive real constants. The existing strict bounded-Gram
lemma is therefore applicable. Restoring the retained event packets is exact
and uses only nonnegative pure-source and pure-testing factors.

`check_homogeneous.py` verifies 21 targeted exact polynomial identities and
determinants with SymPy 1.14.0. Its output is `exact-checks.json`, always beside
the script. It checks all three complete homogeneous parents, centered
extraction in each chart, R6.5--R6.7 kernel/full residues and determinants,
the exceptional complete allocation and determinant, and all three common
capped endpoint residues. There are no samples or parameter sweeps.

This audit does not establish source nonemptiness, all native source formula
derivations, the finite-contact conductor, the finite-remote join, a fixed
affine lift, or the global finite-SDP theorem. The R5 repair supplies the
separate real-cap paragraph referenced by the proposed TeX as
`lem:c2-r55-real-cap`; the existing bounded tensor-Gram lemma is reused.
The proposed fragment has not been inserted into or compiled as the manuscript.

## Files delivered

- `r65_r67_dependency.tex`: proposed integrated subsection with the original
  component aliases, R6.5--R6.7, all exact homogenizations, and the complete
  common capped endpoint proof including `(infinity,0)`.
- `check_homogeneous.py`: bounded reproducible exact verifier.
- `exact-checks.json`: 21 exact checks and explicit scope limitations.
