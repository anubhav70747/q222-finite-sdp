# Zero-computer theorem lock for \(\mathcal Q_{2,2,2}\) — corrected freeze

**Freeze date:** 31 August 2026  
**Lock status:** `LOCK_ZERO_COMPUTER_THEOREM_LEVEL_CORRECTED`  
**Canonical claim:**
\[
\boxed{\operatorname{sxdeg}(\mathcal Q_{2,2,2})=2.}
\]

## 1. Canonical frozen object

The theorem proof is frozen as a **composite exact object**:

1. Base parchment, exact bytes:
   - `Q222_SXDEG2_CONSOLIDATED_PARCHMENT_FINAL_2026-08-30.md`
   - SHA-256: `e268a6ae1c11dd02b5819064fa7e8187173074050d92f3204e703dd46f574731`
   - 2,415 lines.
2. Surgical correction overlay:
   - `Q222_SXDEG2_CORRECTION_PATCH_AND_FREEZE_2026-08-31.md`
   - SHA-256 recorded in `SHA256SUMS_Q222_SXDEG2_FROZEN_2026-08-31.txt`.
3. Corrected human reader route:
   - `Q222_SXDEG2_READER_ROUTE_FROZEN_2026-08-31.md`
   - SHA-256 recorded in the same checksum file.

The base bytes plus the six exact replacements in item 2 are the authoritative full parchment. The reader route is explanatory and subordinate to that composite proof.

This composite convention is intentional: it prevents an unnecessary reserialization of the 2,415 already audited lines and makes every post-audit change explicit.

## 2. Repairs incorporated

Exactly six presentation/self-containment repairs are frozen.

### R1. Bell visible-facet uniqueness

The proof now explicitly records that facets of \(G\) are indexed by deterministic and PR vertices of the binary no-signalling polytope, with
\[
\langle q,C_s\rangle\le1
\]
for deterministic \(q\), and
\[
\langle q_t,C_s\rangle\in\{2,-2,0\}
\]
for PR vertices according as \(t=s,-s,\) or otherwise. Hence \(C_s\) sees exactly \(F_s\).

### R2. Canonical real measurement chart

Before the nine-channel slack, the proof restores
\[
A_0=B_0=Z,
\qquad
A_1(u)=\frac{(1-u^2)Z+2uX}{1+u^2},
\qquad
B_1(v)=\frac{(1-v^2)Z+2vX}{1+v^2}.
\]
Thus the direct-substitution identity for \(P_\star\) is self-contained.

### R3. \(F_1\) gauge typing

The collar map is frozen as
\[
p:=Bs^2-u,
\qquad
\phi_1(s,u,a,B)=\left(u,s,\frac{as}{p},\frac1p\right),
\]
with
\[
\det D\phi_1=\frac{s^3}{p^3},
\qquad
\lambda_{F_1}=s^3p^7>0.
\]
The previous omission of the definition of \(p\) is removed.

### R4. Zero gauge denominators

The proof now distinguishes two boundary types before division:

- finite algebraic boundary factors route through the explicit strict children;
- affine-chart denominator zero routes to the homogeneous one-projective or simultaneous-projective analysis.

No reciprocal/birational quotient is taken on either boundary.

### R5. \(B_0\) factor typing

After the literal four-square suspension, the proof no longer says that every scalar coefficient is wholly point-side. The correct statement is that each coefficient separates into a nonnegative pure-point factor and a nonnegative pure-support square/elementary factor. Absorbing the latter into the support vector preserves separation rank at most two.

### R6. A17 reader-route pivots

The stale reader-route equation with \(K_3\)-dependent pivot additions is deleted. The only frozen centered pivots are
\[
\boxed{
q_\beta=\frac4{T^3(1+V^2)},
\qquad
q_{e_1}=\frac4{T^3V^2(1+V^2)}.
}
\]
The critical value \(m=u^3K_3\ge0\) remains exact and is restored separately as a pure-side square. Thus this correction does not reinstate the historical error of setting a positive infinitesimal critical value to zero.

## 3. Re-audit verdict on affected equations

The affected equations have the following status.

| Object | Frozen status |
|---|---|
| Unique visible facet of \(C_s\) | PASS, proof restored |
| Pairwise-midpoint visible-star lemma | PASS, unchanged |
| Half-angle chart \(A_1(u),B_1(v)\) | PASS, restored |
| Nine-channel slack \(P_\star\) | PASS, unchanged |
| \(F_1\): \(p=Bs^2-u\) | PASS |
| \(\det D\phi_1=s^3/p^3\) | PASS |
| \(\lambda_{F_1}=s^3p^7\) | PASS |
| Zero-denominator routing | PASS, wording disambiguated |
| \(B_0\) four-square identity | PASS, unchanged |
| \(B_0\) tensor-factor interpretation | PASS, corrected |
| A17 \(A_4\succ0\) | PASS, unchanged |
| A17 \(q_\beta=4/[T^3(1+V^2)]\) | PASS |
| A17 \(q_{e_1}=4/[T^3V^2(1+V^2)]\) | PASS |
| A17 critical-value restoration | PASS, unchanged |

No corrected item alters the local-node classification or the induction invariant.

## 4. Frozen theorem chain

The canonical upper-bound chain remains
\[
P_{K_\star}=\operatorname{cone}(E),
\qquad
K_\star=\operatorname{conv}\operatorname{Ext}(K_\star),
\]
followed by Scheiderer's quantitative criterion and the universal estimate
\[
\operatorname{sosx}(e^\otimes(x))\le2.
\]
The local zero-residue nodes are exhausted by
\[
\mathrm{REG},\ Q_0,Q_1,Q_2,\ F_1,F_4,F_3,\ F_3^{\rm term},
\]
\[
\text{one-projective},\ A17,\ B_0,\ B_\infty,
\]
or a strict child. Smooth nodes use exact diagonal compression; exceptional nodes use the printed exact squares/gauges or route to a child. Finite descent is controlled by
\[
\kappa=(N_{\rm act},d_{\rm stratum}).
\]
Therefore
\[
\operatorname{sxdeg}(K_\star)\le2,
\quad
\operatorname{sxdeg}(\mathcal B_Q)\le2,
\quad
\operatorname{sxdeg}(\mathcal Q_{2,2,2})\le2.
\]
The attained cubic correlator boundary
\[
C_{11}=4t^3-3t,
\qquad 1/\sqrt2\le t\le1,
\]
is nonpolyhedral, so
\[
\operatorname{sxdeg}(\mathcal Q_{2,2,2})\ge2.
\]
Hence
\[
\boxed{\operatorname{sxdeg}(\mathcal Q_{2,2,2})=2.}
\]

## 5. Forbidden-dependency screen

The corrected frozen proof has no logical dependency on:

- numerical SDP;
- computer algebra output;
- stored machine Gram matrices;
- Newton/Pareto enumeration;
- a machine-generated contact atlas;
- a fixed production-cell table;
- the invalid generic filtered-collar/R2.4 implication;
- mixed/weighted Bell coordinates as SOS roots;
- identification of independent point/support tensor copies;
- the inference “different specializations imply positive residue”;
- positivity in a non-Archimedean field in place of positive residue;
- deletion of the positive infinitesimal A17 critical value.

Historical computations may remain as discovery or independent checking records only.

## 6. Supersession rules

For theorem-level \(\operatorname{sxdeg}=2\) work after 31 August 2026:

- **Use** the 30 Aug base parchment only together with the 31 Aug correction overlay.
- **Use** `Q222_SXDEG2_READER_ROUTE_FROZEN_2026-08-31.md` instead of the 30 Aug reader route.
- **Do not cite or reuse** the old reader-route R23 \(K_3\)-dependent pivots.
- Older fail-closed audits predating the 30 Aug analytical completion are historical frontier records, not the current theorem status.

## 7. Scope deliberately left open

The freeze is an existence theorem for finite \(2\times2\) PSD / \(L_3\) lifts. It does not provide a literal solver-ready affine lift and does not determine the exact minimum \(M_{\rm SOC}\).

The next independent constructive gates are therefore:

1. materialize a finite uniform rank-two Bell atlas and compile a literal affine \(L_3\) lift;
2. only after that lift is frozen, settle the minimum Lorentz-block count.

---

# FINAL LOCK

\[
\boxed{
\texttt{Q222\_SXDEG2 = THEOREM\_LOCKED\_ZERO\_COMPUTER\_CORRECTED}
}
\]

\[
\boxed{\operatorname{sxdeg}(\mathcal Q_{2,2,2})=2.}
\]
