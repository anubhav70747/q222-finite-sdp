# Branch map: minimal-complexity extension attachment vs Q222 finite-SDP proof

Read-only source map, 27 October 2026. The attached `original/main.tex` and `original/supplement.tex` were read as claims, not treated as proof by their labels. The established Q222 boundary below follows the current repository's FRONTIER, CHECKLIST, RECONCILIATION and VERIFICATION_PLAN, together with their cited source proofs. No new derivation or verifier was run here.

## Results and exact domains

| Result | Mathematical object and domain | Evidence/dependencies | Status and exact boundary |
|---|---|---|---|
| Symmetric doubly tilted support epigraph; affine lift with PSD orders 5 and 4 | `E_sym={(alpha,lambda): 0<=alpha<=1, lambda>=omega_Q(alpha)}`. One affine SDP in visible variables `(alpha,lambda)` and Gram variables; all measurement angles are universally quantified then eliminated. The local bound implies `lambda>=2+2alpha`, the domain used in the determinant/characteristic argument. | Attachment SM: two-qubit Jordan reduction; 4x4 slack has nonnegative first three shifted characteristic coefficients on the local-bound range; determinant reduction; global diagonal identity `D(x,y)-D(sqrt(xy),sqrt(xy))=64(4-alpha^2)xy(sqrt(x)-sqrt(y))^2`; tangent-half-angle maps remaining interval to real line; degree-four 3x3 matrix polynomial; univariate matrix SOS with parity split into 5x5 and 4x4 Gram blocks and 16 affine coefficient equations. | The stated lift is for the support epigraph and is structurally single-lift across alpha. Its theorem relies on the cited Jordan/real two-qubit model reduction and univariate matrix SOS theorem. It is not a lift of all behaviors or the complete Q222 set. Keep this claim independent of the older finite-SDP existence argument for Q222. |
| `sxdeg(E_sym)=2` | Same planar convex semialgebraic epigraph. | Upper bound from planar semidefinite extension theorem (`sxdeg<=2`); lower bound from nonpolyhedral boundary using `omega_Q(1-T)=4-2T+(4/3)T^3+O(T^4)`. | Conditional on the cited planar theorem and cited exact boundary expansion. Not a statement about extension degree of Q or W3. The explicit (5,4) lift is not itself a minimal-block construction; the planar theorem supplies the upper bound 2. |
| Symmetric finite-NPA nontermination comparison | Fixed finite standard NPA level L, limit `alpha -> 1` of NPA-minus-quantum gap divided by `(1-alpha)^3`. | The attachment imports the prior result `+infinity` and Motzkin obstruction from ChaturvediNPA2026; true quantum violation coefficient `4/3` imported/cited. | Not proved within this attachment's displayed SDP derivation. It is a cited input. Do not describe the new lift as a proof of NPA nontermination. |
| General asymmetric doubly tilted support optimization | `h_{alpha,beta}=CHSH+alpha A0+beta B0`, nonnegative tilts. Local/supercritical region `alpha+beta>=2` has local value and a `1+AB` SOS. In subcritical region, Jordan measurement square is reduced to weighted ray `(4-alpha^2)(1-c)=(4-beta^2)(1-d)` after party exchange; clipped edge is separately nonnegative. | SM critical/supercritical SOS; determinant identity and first-three-coefficient condition above local value; minimize `A x+B p/x` at fixed `xy=p`; clipped-edge factor signs; positive diagonal congruence; degree-six 4x4 matrix polynomial on `[0,1]`; interval matrix factorization. | The optimization reduction is asserted global on the stated nonnegative-tilt domain. The finite SDP orders 16 and 12 are for each fixed pair `(alpha,beta)`; coefficients of the affine SDP depend on the fixed ratio `rho=(4-alpha^2)/(4-beta^2)`. This is not one affine lift with `(alpha,beta)` as free visible variables. Degenerate edges/tilts require the stated limiting/critical branches; do not infer uniformity across parameter pairs from fixed-input factorization. |
| Algebraic boundary for general tilts | Fixed tilt pair; boundary value lambda. | Quartic `B^2 D(z,rho z)` has a double root at an interior boundary; eliminate root and derivative; symmetric specialization factors to the cited quartic. | Derivation is presented in SM but source equation/root-domain and elimination details are compressed. Treat as a separate boundary-result proof obligation from the finite SDP; no behavior projection follows. |
| Full behavior projection `W3` | Three-dimensional image `(A0 marginal, B0 marginal, CHSH)` of Q. | Attachment explicitly observes nonlinear weighted ratio obstructs promotion of fixed-input SDP to one affine lift of W3. | Open in attachment. The support-function/evaluation results for each fixed Bell functional do not supply a common lift of W3. |
| Qualitative finite-SDP existence for complete Q222 | Complete eight-dimensional quantum behavior set Q / normalized source body, via universal real-closed-field tensor-SOS criterion and source/contact cover. | Existing established route: physical two-qubit reduction and compact nine-atom model; universal tensor-SOS criterion; canonical source reduction and stationary contacts; noncommuting residue theorem; 215 source components, all 45 C1 cells, all 49 C2 leaves; V01-V06 conjunction; full-polar assembly given exact canonical input. | Reviewed qualitative existence theorem remains established, independently of current evaluated candidate. Retain scope and hypotheses from proofs; counts and reconciliation are source navigation, not substitutes for leaf implications. This is existence, not currently an evaluated single global affine pencil. |
| Current explicit global outer SDP candidate | Full Q222 Bell-validity cone via evaluated rational matrices: 897 PSD blocks, total order 51088, 423168 auxiliary scalars, 3456 affine equalities; full 128 physical copies, orthants, event/face blocks and source adjoints. | L18 conditional pencil; quantum-to-feasible inclusion checked. Reverse inclusion follows from source duality only after universal canonical strict admissions. | Evaluated outer SDP is not yet proved exact. Its endpoint requirement is K364 direct admission for every genuine physical strict C1 source (`0<k<=1`, H1) and every genuine physical strict C2 source on the complete native domain (H2), or adequate alternative fixed globally sound cone coverage. |
| C1 and C2 branches for exactness | Native strict source families. C1 has two orientations; C2 includes finite `c`, n-wall and coupled approaches. | C1: fixed-span Z16/Z18 coefficient reductions and exact Schur reductions; unrestricted full-A1 Z16 admission is refuted, oriented Z18 universal admission open. Balanced C1 slice/tube (L20) is only a five-parameter neighborhood; C1 endpoints/contact mechanisms are separate and unjoined globally. C2: exact admission near both strict open c-ends (L12), local n=0-to-positive opening under all-four wall-complement PD (L19), but arbitrary finite c/coupled wall approaches remain open. | Neither local certificates, endpoint bounds, contact jets, partial minors, map rank, nor samples prove H1/H2. Current exactness is explicitly conditional. New attachment does not fill either source-admission branch. |
| Closed faces, event rays, reflection and singular density | Boundary strata of Q and orbit assembly in behavior coefficient space. | Exact face cones SG12/SG13 and event cones are retained. L17 reflection transports strict C1 opposite orientation into reflected K364; use full physical 128-element orbit (equivalently eight sector representatives plus canonical stabilizer). L15 proves that strict admission into one fixed closed globally sound cone implies permitted singular-boundary coverage by coefficient-space closure and finite hybrid assembly. | These are downstream conditional implications: once H1/H2 (or alternative strict common-cone coverage) hold, no separate uniform limiting-Gram bound is needed. Density does not establish strict admission. The full-polar/eight-copy construction only assembles an exact canonical input; it does not produce the missing input. |

## Missing implications, stated as proof obligations

1. To derive the symmetric `(5,4)` lift: establish/cite the two-qubit Jordan model and exact characterization of the support epigraph; justify determinant iff PSD on `lambda>=2+2alpha`; use the global diagonal identity to reduce the complete square; apply univariate matrix SOS with parity decomposition; match all 16 coefficients. The attachment gives these ingredients in SM, with detailed algebraic checks concentrated there. For the extension-degree equality also invoke planar upper-bound theorem and the nonpolyhedral boundary expansion.
2. To derive a general-tilt fixed-input SDP: handle the local/supercritical branch; prove determinant reduction in the subcritical physical domain; prove weighted-ray minimization including clipped edges; check congruence is positive/invertible in the interior; apply `[0,1]` matrix SOS; separately state limits at degenerate parameter values. None of this makes `(alpha,beta)` free affine coordinates.
3. To turn fixed-tilt values into a lift of W3, a common finite affine representation in the three behavior coordinates is still needed. The nonlinear `rho` dependence is the stated obstruction; no such implication is supplied.
4. To make the current Q222 outer pencil exact: prove H1 and H2 over their complete native domains (or another sufficient fixed finite sound-cone cover), then invoke L15 closure/density, event and face cones, orbit assembly, and L18 source-dual reverse inclusion. V01-V06 qualitative existence does not prove H1/H2 for this particular finite library.
5. For H1, prove direct `K364` admission for all physical strict C1 sources, including arbitrary near-wall points and both orientations (reflection plus common-cone sum handles the latter only after base orientation admission). Existing balanced tube, endpoints, contact jets, Z18 reductions and scalar-Schur gates do not cover/join the full domain.
6. For H2, prove direct `K364` admission for all strict C2 native sources, including arbitrary finite c and coupled approaches to zero-weight walls. Open-end theorem and n=0 local opening are scoped neighborhoods only.

## Dependency graph (plain text)

ATTACHMENT SYMMETRIC:
Jordan/two-qubit reduction + support epigraph definition
  -> local-bound range and 4x4 slack determinant criterion
  -> global two-angle diagonal reduction
  -> one-variable quartic matrix positivity
  -> univariate matrix-SOS theorem + parity split + 16 coefficient identities
  -> one affine (5,4)-block lift of E_sym
  -> [planar sxdeg theorem] upper bound 2
  -> [cubic boundary expansion] nonpolyhedral -> lower bound 2

ATTACHMENT GENERAL TILTS:
Jordan block + local/supercritical SOS
  -> subcritical determinant formula
  -> weighted ray + clipped edge proof
  -> fixed `(alpha,beta)` degree-six matrix polynomial
  -> interval matrix-SOS factorization
  -> fixed-input (16,12)-block SDP
  -X-> common affine lift of W3 (nonlinear rho; attachment explicitly leaves open)

EXISTING Q222 GLOBAL CANDIDATE:
V01-V06 qualitative theorem (215 components / 45 C1 / 49 C2)
  -> qualitative finite-SDP existence for Q (separate result)

H1: all strict physical C1 -> K364 admission -- OPEN
H2: all strict physical C2 -> K364 admission -- OPEN
  -> Kbar = K364 + reflected K364 (L17; full 128-action orbit)
  -> L15 coefficient-space singular density / closed hybrid sum
       + event rays + SG12/SG13 closed-face cones
  -> exact canonical source-cone input
  -> L18 source-dual reverse inclusion
  -> equality for evaluated global Q222 pencil (conditional; not yet established)

L18 rational evaluated outer pencil -> Q-to-feasible inclusion (established)
H1/H2 are the missing reverse-inclusion inputs; attachment's E_sym and fixed-pair tilted SDPs have no arrow to H1/H2 or to W3.
