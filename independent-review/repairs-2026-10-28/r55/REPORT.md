# R5.5 dependency repair

Baseline: `d01e5e089c5b4ad73667e922eec5b62b5d82485d`.
Scope: the centered `m` extraction and the finite remote R6 consumer,
not the global Q222 finite-SDP theorem. Only this repair directory was
intentionally edited. The proposal has not itself been integrated or compiled
by this reviewer.

## Finding and precise repair

The original argument was omitted at the cited destination. In the baseline
`latex/certificates/c2_short.tex:14299-14300`, both R5 component labels have
no following mathematical body. The retained paragraph at 8102-8117 already
contains the centered row and correction matrix, but does not print the
four-state event identity R5.5 or its finite remote source-domain argument.

The recovered complete source is
`complete-proof-materials/materials/historical/recovery/full_r21/components/source/p8da9ea6ed9e5.tex`.
I read this whole file and the immediate complete WA table, source formula,
tensor-Gram proof, and three R6 consumer files. R5.2-R5.5 occur at original
lines 41-77. They give a correct exact full-parent extraction. The equation
is universal event algebra and therefore transfers beyond the R5 source
domain. The R5 source theorem does not transfer: that theorem has bounded
`a` and infinitesimal `r`, whereas the present R6 consumers have unbounded
`a` and `r` of residue one.

`r55_dependency.tex` supplies a self-contained subsection for the restored
component destination. It includes the contact-centered row, full event
extraction, exact matrix extraction, a real cap proof, the actual unbounded-a
R6 source estimates, the complete four-state finite remote kernel form, its
strictness proof and legal separated tensor-Gram step. No other R5 stratum
is reconstructed or admitted by this proposal.

## Exact verified scope

1. R5.2 is a four-state polynomial identity, and its row vanishes at the
   original stationary contact. It follows from `g=1-a(1-r-hB)`, not a
   comparison between source and testing scales.
2. R5.5 subtracts the **whole** square `(m-K)(M-gvW)^2` from the complete
   event parent. Every original `n,C,A,D,alpha` packet and the entire signed
   core remain in the remainder. No restriction to the contact graph,
   Hessian, state kernel, or scalar soft form is made at this step.
3. All 15 independent entries of the unnormalized complete WA matrix obey
   `partial_m G=L_m L_m^T+g H_m`. This verifies R5.3-R5.4 with both normal
   rows and all hard entries, holding all other actual source coefficients
   fixed.
4. A subtle distinction is essential: the formal WA-table substitution
   `G|_{m=K}` is generally not the Gram of the event substitution
   `P|_{m=K}`. The WA table has already used the original contact balances.
   Both formulas correctly express the actual parent after the same whole
   excess square is removed. Neither claims a modified physical source.
5. For the actual unbounded-a R6 source, use `s=1/a`, `b=a^2B` and
   `r=1-s(1+hsb-g)`. Original physicality forces `b` bounded; after exact
   rational cancellation the joint residues are
   `mg -> 1/2`, `n/s^2 -> [1-b_0(h_0+t_0)]/2`,
   `A -> t_0/(2h_0)`, `C -> (h_0+t_0)/(2l_0)`,
   `D -> 1/[2(t_0-l_0)]`, and
   `alpha/(s^3b^2) -> (h_0+t_0)/2`. Thus `n_0=alpha_0=0` and all
   other uncapped event weights are bounded. This covers unit and
   infinitesimal `b` and literal `n=0`. No ratio such as `s/g`, `g/b`, or
   a source/testing ratio is assumed.
6. For bounded tests with `v_0 != 0`, the centered correction has residue
   `vWM_0`, which vanishes on `M_0=0`. The complete kernel determinant
   RT10/RB24 equals exactly
   `D_u^2 D_v^2(h+t)/(8hl(t-l)) * [V^2+t(h+t)U^2V+hUV^2]`.
   The `Z=0` restriction is positive definite from the original A event
   and the nonnegative signed-core squares. Since `V>0`, the determinant
   is positive. These two facts prove kernel positive definiteness. A
   sufficiently large **real** cap then makes the full four-state residue
   positive definite; unbounded source `m` permits every such real cap.
7. The remainder coefficients lie in `B_c tensor_R B_t`, because they are
   finite sums of bounded pure-source times bounded pure-testing factors.
   State roots can be unbounded. The restored square uses the pure-source
   root `sqrt(m-K)`. No mixed inverse or mixed radical is introduced.
8. The focused checker also verifies the homogeneous centered-correction
   terms at the single and double reciprocal endpoints. This is only
   correction algebra, not an independent proof of those endpoint joins;
   their complete proof is assigned to the separate endpoint repair.

## Immediate dependencies and destinations

| Object/use | Baseline/current destination | Original evidence | Proposed readable destination |
| --- | --- | --- | --- |
| Complete parent | `c2_short.tex:8005-8014`, and RB5 at 12152-12160 | `pb60eba754641.tex:59-67` | Existing `eq:c2-ba-parent` |
| Complete WA table | `c2_short.tex:7584-7700` | Entire `p9ec3d6bbd717.tex`; table at original 74-179 | Existing WA3-WA8, plus new `eq:c2-r55-matrix` |
| Centered row/extraction | `c2_short.tex:8102-8117`; empty R5 destination 14299-14300 | `p8da9ea6ed9e5.tex:41-77` | New `eq:c2-r55-centered-row`, `eq:c2-r55-event`, `eq:c2-r55-matrix`, `eq:c2-r55-correction` |
| Unit-b finite remote use | `c2_short.tex:12398-12410` | `pb60eba754641.tex:309-320` | New `eq:c2-r55-event` and `eq:c2-r55-r6-kernel` |
| Imported-identity claim | `c2_short.tex:12434-12436` | `pb60eba754641.tex:357-358` | Point directly to `eq:c2-r55-event`, not an empty component label |
| Infinitesimal-b finite remote use | `c2_short.tex:12884-12886` | `p75b395956985.tex:439-441` | Same whole-parent identity and actual R6 finite remote proof |
| Thin-layer source consumer | Recovered source, not a separate active body destination | `pc82cd677f4b9.tex:173-185` | Same event identity transfers; its source assumptions are a subdomain of the actual R6 bounds verified here |
| Bounded tensor-Gram lemma | `latex/chapters/03_tensor.tex:77-112`, included by `Supplement.tex` | Entire `p25af2c58da7e.tex` | Existing `p25af2c58da7e:lem:gram`; no duplicate lemma required |

The earlier reviewer note that the Gram label was missing was based only on
searching `c2_short.tex`; it is corrected by the complete-document location
above. The current complete lemma is present and adequate.

## Why the bounded tensor-Gram lemma applies

The complete proof was read in both the original and current tensor chapter.
For any bounded tensor of positive real residue, its finite separated
expansion is a positive real reserve plus terms with an infinitesimal in
at least one bounded pure factor. Split the reserve into finitely many real
parts. For infinitesimal pure `x` and bounded pure `y`, choosing real
`L>|y|` makes each pure factor in

`epsilon_j + x tensor y = [(epsilon_j/L+x) tensor (L+y) + (epsilon_j/L-x) tensor (L-y)]/2`

positive in its own field. Their separate field square roots give tensor
squares. Swapping factors handles the other case. A real congruence reduces
a strict matrix residue to the identity; the fixed real diagonal reserve
decomposition TF6 then applies the scalar argument entry by entry. The
argument proves a finite SOS, and never assumes that an arbitrary
positive-residue tensor is a unit or takes a mixed inverse. This matches
the finite-remote remainder's stated coefficient-ring hypotheses.

## Checks and limits

Run `py -3.14 -X utf8 independent-review/repairs-2026-10-28/r55/check_r55.py`.
The focused SymPy checker completes with 29 exact identities recorded in
`exact-checks.json`. It contains no numerical sweep, solver, baseline replay
or identification of the two labelled fields. It checks the named algebra;
the positivity and source-domain deductions are the written proof above.

No further mathematical gap was found in the **isolated extraction and
this finite remote R6 use**, once the actual source estimates, complete
kernel restriction, real cap and tensor-Gram step are supplied. The
extraction by itself proves neither positivity nor tensor SOS. This review
does not certify a contact conductor, all endpoint charts, exhaustive
stratum coverage, a finite-dimensional SDP representation, a universal
quantifier over every source, or the global theorem. Those implications
remain outside this bounded task.
