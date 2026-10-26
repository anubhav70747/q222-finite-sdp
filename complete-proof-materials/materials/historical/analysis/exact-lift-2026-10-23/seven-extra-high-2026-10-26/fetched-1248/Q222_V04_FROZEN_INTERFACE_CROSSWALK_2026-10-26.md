# V04 frozen interface check: C1 BN/UN, C2 R2/UL, R3/R4, R6/R7

**Frozen reference:** `anubhav70747/q222-finite-sdp`, main `332483309586e0fec6e16fc1dd855915b648dca1`; `latex/certificates/c1_short.tex`, `latex/certificates/c2_short.tex`. Date 26 October 2026. Source evidence: user-supplied `05_V04_VERSION_INTERFACES_SOURCES_1218.txt`, its corrected/original rendered proofs and the retained independent 22 October C2 audit; direct read-only check of the frozen native TeX. **No repository mutation, numerical sweeps, new polynomial checker, or fresh global proof.**

## Verdict

**No operative version/interface mismatch found in the four specified changes.** There is one *historical* BN/UN typo, already repaired in frozen main. Each replacement has an exact and narrower importing domain, retains the original entire parent and its two unrestricted state normals, and imports only an algebraic identity/full-parent selected-frame theorem within its stated hypotheses. This is a *composability crosswalk*, not an independent redo of the 45 C1 / 49 C2 local chains. No statement about the unresolved full qualitative conjunction or an evaluated uniform finite SDP pencil follows.

Notation: source coefficients in a real closed field `K_c`, independent test/state in `K_t`; a separated certificate lives in `K_c \otimes_R K_t`. Roots may be unbounded. Every Gram matrix coefficient **must** lie in the bounded separated coefficient ring and have a positive-definite real residue after the displayed packets have been retained. Relative source/test rates are unrestricted. Where a row's coefficient is literally zero, omit it rather than invert it.

## I. C1 BN/UN: normalization already correct in frozen main

**Source scope:** Physical C14.4, `S,\kappa,\nu` positive infinitesimals, `\ell=1-\nu` a unit, `q=\nu/\ell`, `\zeta=\ell/\kappa`, native gaps and `D,\beta\ge0`. BN has bounded `\tau=h/q`; UN unbounded `\tau`, with `e=\kappa\tau/\ell` infinitesimal. Testing Alice `u` and reciprocal Bob `w` infinitesimal (literal zeros allowed), with unrestricted state coordinates. The complementary CE/SE leaves are external with their independent angle conditions.

The *frozen* `c1_short.tex:13529` already states

\[
 N=\kappa W-\ell uwZ,\qquad
 \boxed{N_0=N/\ell=W/\zeta-uwZ},\qquad
 R=w(a-x)Z=\frac{\sqrt S}{\zeta}G-uN_0.
\]

The old proposal's `N_0=N/\kappa` is incompatible with both displayed equalities; it is **not** current. Specifically, `((A_N/\kappa^2)N^2)=\zeta^2 A_N N_0^2`, which is the BN8/BN9 coefficient `A=\zeta^2A_N-b`, and `N(wB_u)=\ell N_0(wB_u)` in the complete scalar SN14. The circuit subtraction changes the latter coefficient to

\[
 C'_N=-H_{11}-2\ell L_c/\kappa,
\]

exactly BN11–BN13; the original `n_sB_u^2` is spent by the circuit, not omitted. The complete BN scalar is

\[
\begin{aligned}
F_0={}&\Psi+A N_0^2+A'_P G^2+d_sB^2+H_{02}(wB)^2+C_RR^2\\
&-2\sqrt S L_aGR+C'_G G(wB)+C'_N N_0(wB_u).
\end{aligned}
\]

The complete **first normals**, not only graph restrictions, remain

\[
 L_X=\ell w\widehat A_XN_0+\frac{uw\sqrt S}{\ell}\widehat X_aG+d_s\sqrt S\,q u B-bB_u+wX_c(wB_u),\quad
 L_Y=\sqrt{r(1-r)}(Y_PG+\sqrt S\,Y_aR+Y_c(wB)).
\]

The original hard entries are `H_X=mc+da+(\beta+1/2)(1+c)`, `H_Y` and `H_XY` of SN18. The BN pivot spending order is **first** `wB` by BN14, **then** `wB_u` by the normalized SN8/BN15, **then** the remaining `bB_u` coupling by BN17. On the zero face `\beta=0`, the second circuit chart is not invoked. On `D=0`, the `d_s` row is omitted and the complementary sphere coordinate supplies the helper.

UN uses the same corrected `N_0` but **different** `N_x=X+huZ` and actual circuit `\Psi_n`. Before any bounded estimate the radial terms are combined using

\[
 A+aC_R=A_{00},\quad
 R=(\sqrt S/\zeta)G-uN_0,
\]

so `UN14` has bounded coefficients, including `C_R/\zeta` and `C_R/\zeta^2`. The corrected original-source identities give

\[
 A_{00}=\frac{\theta\ell H_0}{2P_0}
 +\frac{eT_\kappa+h\ell ST_\nu}{2JMP_0},\quad A_{00}-\theta/2\in\mathfrak m_c;
\qquad t_0=\frac{D+h/2-\beta\theta}{q}\sim\frac12.
\]

The one potentially unbounded X-normal term is removed by the *entire-state* identity

\[
 dG_h^2+daN_x^2+2d\tau a(wN_0+B_u)N_x
 =da(X+ewW)^2-2de\tau uw\sqrt S\,(wN_0+B_u)G+de^2Sw^2G^2.
\]

Retain the positive `da(X+ewW)^2`, `\beta(1+c)N_x^2`, circuit remainder, and nonnegative exchanges; `H_X` then has unit residue `1/2`, `H_Y` residue `1/2`, mixed hard residue zero. UN orders `wB`, then `N_0` from UN21 if `\theta` infinitesimal (otherwise native `A_{00}`), then `B_u` by UN20, finally `wB_u`; no inverse of `\theta` or a circular pivot. Literal `D=0` is included, while `\beta=0` is *excluded* from the unbounded-`\tau` physical source by the actual beta equation (but is included in BN).

**Join:** CE1–CE10 only for noninfinitesimal Alice; SE10–SE12 only for unmatched homogeneous Bob residue `\operatorname{res}(s_B^2)\ne\bar\ell=1`; the residual reciprocal-Bob local corner is BN/UN. These are algebraic/selected-parent inputs, not an invocation of all-C1 physical validity. **V04 action:** no edit required to current definition or downstream formulas.

## II. C2 R2 and UL: the reflected `\bar k=1` simplification is *not* global

**Source scope:** Native physical R2: `t_s` unbounded, `B_s,a_sB_s,t_sB_s` infinitesimal, all gaps, `n\ge0`. The correct local split is bounded `a_s` (RL); unbounded `a_s` with infinitesimal paired parameter `p` (RL); unbounded `a_s` with `p` a unit and `\lambda_n=n/(t_sB_s)` bounded (TU); or unbounded `a_s`, unit `p`, and positive unbounded `\lambda_n` (UL). No source is transferred across this split by continuity.

For `e=1-r` and `g_s=1-a_se+ha_sB_s>0`, one has `a_se<1+ha_sB_s\sim1`. Thus exactly

\[
 k=\frac{1+a_sr}{1+a_s}=1-\frac{a_se}{1+a_s}.
\]

**Only when `a_s` is unbounded** does this imply `\bar k=1`, since `a_se=O(1)` while `1+a_s` is unbounded. In UL/TU coordinates `z=1/a_s`, `w=a_se`, the same identity is

\[
\boxed{k=1-\frac{zw}{1+z}},\qquad \bar k=1.
\]

The original R2O full-parent reflection is an exact homogeneous identity, not a source-sector equivalence:

\[
 G^TG=D_uD_vI,\quad F_{c'}(v,u;G\psi)=D_uD_vF_c(u,v;\psi),\quad
(\lambda',n',m',C',A',D',\alpha',\beta')=(\alpha,m,n,C,D,A,\lambda,\beta).
\]

`R2O10` covers reflected Bob noninfinitesimal, its 2-by-2 bounded current Gram having residue `I_2`. For reflected Bob infinitesimal and reflected Alice noninfinitesimal, the *complete four-state* R2O11 Gram has residue

\[
 \operatorname{diag}(1,2,2,\bar k^2-\bar c_u^2),\qquad c_u=(1+u^2)^{-1/2}.
\]

With `\bar k=1`, the final pivot is positive because reflected Alice is noninfinitesimal (`\bar c_u^2<1`, including projective infinity). No finite cap or R2O13–24 is needed **in these three unbounded-`a_s` subdomains**. For bounded `a_s` with `\bar k<1`, retain R2O13–R2O24 and their four real-cap opening cases, including the two-current normals and kernel Schur margins. The bound `C\ge c_*>0` in the separate middle `bounded-a` proposition is not being used to replace R2O's bounded-`a_s`, unbounded-`t_s` case.

For TU, the original separate `n` packet `4nD_uD_vX^2/S` is retained as a contact-centered whole square **plus its exact signed remainder** (`TU17`); it is omitted only when `n=0`. For UL, `\lambda=n/(qz)>0` is unbounded so **`n>0`**; the same current is completed with `J_3,V_3`, full normals, and the toric pivots in Part IV. Formal same-contact C1 source coefficients may have `\beta^0<0`; only the generic SC/QI polynomial identities are imported, never their physical C1 SOS theorem.

**Composability:** R2O10/11 use the actual reflected full parent and the bound `Sg_*\sim k/2>0`, not an unproved full-R1 or full-C2 theorem. This is exactly the source-specific simplification, not permission to delete the bounded-`a_s` real-cap theorem.

## III. C2 R3/R4 and R6/R7: genuine bodies vs historical aliases

The frozen source partition (`c2_short.tex:184–261`) first splits the directed weak-II chart by the original native `t_s` scale. For `t_s` a positive bounded unit:

* `a_s` **bounded**: invoke `prop:c2-bounded-a` (`c2_short.tex:7950–...`), **only** with its printed premises `B,x=aB,y=tB` infinitesimal, `g,j,1-txr>0`, `n\ge0`, and `C\ge c_*>0`. The directed **exactly-one-weak** orientation supplies `C` strong, hence such a positive real `c_*` exists. The **two-weak** complementary-event sources are routed by their separate theorem and are **not** smuggled through the proposition. On the unit-`e` case, `ae<1+hx` forces bounded `a`; actual source identity gives `nr` a positive real residue and excludes `n=0` *only on that unit-`e` subdomain*. Its parent proof uses the complete five-root WA matrix, the exact `D` syzygy `P=(D-d_0)L_D^2+z^T(G|_{D=d_0}+(D-d_0)jH_D)z`, and the centered `m` packet. Because `Dj` and `(m-K)g` are bounded, **both first normals and the hard block** are legal. For infinitesimal `e`, use the different bounded-`a` opening bounds, not the unit-`e` residues.
* `a_s` **unbounded**: native `ae<1+hx` forces `e=1-r` infinitesimal. The distinct remaining source leaves are R3 (`h` unit, `1-h` infinitesimal) and R4 (`h,1-h` units, `j` infinitesimal), admitted by EN/R34 at `c2_short.tex:11619–12092`; R6 (`h,1-h,j` units, `g` infinitesimal), with physical `n\ge0` forcing bounded `b=a^2B`, split into unit-`b` RB (`12093–12444`) and infinitesimal-`b` RI (`12445–12908`); and R7 (`h,g` infinitesimal) by RU (`12909–13435`) with both `b` alternatives and the `t_0=1` opening. R8/R10 remain separately routed and are not reassigned to these leaves.

At `c2_short.tex:14295–14304`, identifiers `p78dfeedee062`, `p1c9902e9a6f2`, `pd5fa361b4623`, `p9a223c96c696` are **bare compatibility labels**, not retained proof bodies. The actual statement and its domain must be taken from the bounded-`a` proposition and the four EN/RB/RI/RU bodies. This is an alias/version interface, **not a lost proof**.

Unchanged local content (not rederived for V04): the EN/R3/R4 parent uses exact RFN7 with unrestricted `N_0,M_0`, the four native polynomial stationary balances `S_N=S_0=S_1=S_2=0`, the correlated `n` completion and RFN18–20 helper budgets. Weak `m` is handled through `\mu=m(1-p)` in the local frame and `\widetilde m=m/(1+m)` remotely; `n=0` omits the current row. In R6/R7, the actual native source cancels `U_F\in(b), U_N\in(sb)`, giving bounded scalar, **both linear normals and all three hard coefficients**, before the CC/SE monic conductor. Exact filtered corrections have zero residue, and the toric plus coupled blocks are strictly positive. The conductor is used only for its *verified* raw-support/parity/divisibility premises, not via an all-C1 physical conclusion.

**Remote imports with their real caps:** EN21–EN25 assume one of `C,D` unbounded and the other plus `A` have strictly positive **real** residues; no positive unit `m` is needed. R6 RB/RI have `m` unbounded and `C,A,D` positive real lower bounds. R7 RU has `m,A` unbounded, `C` a positive unit and `D` positive unit or unbounded (if `t_0=1`). This validates RB22/RB24, RU25–RU30, the original BE selected frame, and literal projective endpoints as whole-state identities, preserving `n,\alpha` as original positive event packets. The local raised-degree conductor itself is **not** evaluated at infinity. No circular import and no domain missing at the historical aliases.

## IV. UL 24-root finite *local* coefficient algorithm

**Exact domain:** C2 R2 with both native `a_s,t_s` unbounded and relative current `\lambda=n/(qz)>0` unbounded, `S=1/t_s`, `z=1/a_s`, `q=a_s(t_sB_s)`, `w=a_s(1-r)`, `R=1+Sqh`, `M=R-w>0`, and independent positive infinitesimal testing `u,v`. Its source equation gives

\[
 \tau=q/(zw)\in\mathfrak m_c,\quad \lambda\tau\sim\tfrac12,\quad
 L=\lambda q\sim zw/2\in\mathfrak m_c,\quad n=\lambda qz\sim z^2w/2>0,\quad
 M_1=nS/z^2\sim Sw/2\in\mathfrak m_c,\quad H_0\asymp z.
\]

The complete original four-state parent has the exact normalized identity (`SC3`, `TU17`, `UL5`):

\[
 F/S=G_{\rm SC3}+\Pi_m+\Pi_v+\frac{4n}{S}D_uD_vX^2,\quad
 \frac{4n}{S}D_uD_vX^2
 =4D_uD_v(J_3-d_N\widehat N)^2+\text{the exact signed scalar/first-normal remainder},
\]

with `\widehat N=\sqrt{H_0}N`, `N_y=\sqrt\sigma(D_uY_b-\rho_0H)`, `J_3=\sqrt\lambda(U_*-B_*)`, `d_N=\sqrt{nz/H_0}` and `V_3=\sqrt{\lambda\Gamma}(\sqrt{U_*}-\sqrt{B_*})`. The algebraic identity

\[
 J_3^2=4V_3^2+\lambda(\sqrt{U_*}-\sqrt{B_*})^4
\]

permits **retaining** `\Pi_{n,1}=2D_uD_v(J_3-d_N\widehat N)^2`, `\Pi_{n,2}=\lambda D_uD_v(\sqrt{U_*}-\sqrt{B_*})^4`, and the exact QI8 positive `\Pi_E` as whole separated SOS packets. Together with `\Pi_m,\Pi_v`, these reconstruct the original `n` term, including its first normal; no negative packet is declared independently positive.

The complete hard block **after these packets**, ordered `(\widehat N,N_y,J_3,V_3)`, is (full second-derivative extraction from SC3, including the `D_bJ_b^2` and core contributions)

\[
\mathcal H=\begin{pmatrix}
1+\frac{zu^2v^2(4D_bS+2)}{H_0}+2D_uD_vd_N^2 &-\frac{2\sqrt{Sz}\,uv(1+2D_bS)}{\sqrt{H_0\sigma}}&-2D_uD_vd_N&0\\
*&\frac{4D_v(c_b+a_b/D_u)+4D_bS^2+2S(1+v^2/D_u)+4\alpha_bzS^2v^2/D_u}{\sigma}&0&0\\
*&*&D_uD_v&0\\
*&*&*&4D_uD_v
\end{pmatrix},
\qquad \boxed{\overline{\mathcal H}=\operatorname{diag}(1,4,1,4)\succ0}.
\]

The **two complete first-normal changes** are the actual-minus-formal TU9/UL10 and UL12 polynomial expressions, with the original paired current kept. They are bounded after exact source-opening cancellations

\[
 1-\rho_0=(1-h)/R,\quad
 \rho_0\widetilde a=-\frac{(1-e)(1+Sq)P_A}{RP},\quad
 (1-\rho_0)\widetilde c=-\frac{(1-e)P_C}{RP},
\]

and the monic `q,U=u^2,V=v^2` differences. Their only apparently unbounded coefficient `w/(c/g)` is replaced by

\[
 \frac{wg}{c}\,vB_*C_*=wzC_*^2-\frac{wz}{\sqrt n}vC_*J_3,
 \qquad wz/\sqrt n=O(\sqrt w).
\]

The complete `UL10` first normal, `UL12` harmonic normal, and the bare-`\Gamma` families have their source contact ratios `UL16–UL18` of zero residue, so the SC17–SC20 **finite contact/toric division** is applicable. It is strictly polynomial/coefficientwise: clear source/testing *unit* denominators, divide monically in `q,U,V`, then in `u-u_c,v-v_c` and `\sqrt u-\sqrt{u_c},\sqrt v-\sqrt{v_c}`; use the printed toric identities to replace differences by existing roots. No mixed inverse, hidden source-opening inverse, or appeal to the *physical* formal-C1 theorem appears.

The exact root list is **16** SC16 toric roots (four `A_i-\Gamma`, three pairings `Q`, six variance rows `D`, three fourth-root pairings `P`), **four** QI roots `(E_*,F_*,R_E,R_F)`, and **four** hard/current roots `(\widehat N,N_y,J_3,V_3)`, totalling 24. The leading strict form is `2T_4+(2-q)E_*^2+2F_*^2+4R_E^2+4R_F^2+\widehat N^2+4N_y^2+J_3^2+4V_3^2`, with strict toric pivots `\ge1,2/3,2/3,4/3` on the SC16 subblocks. All remaining coefficients are computed in a deterministic finite matrix-entry algorithm and have **bounded zero-residue** differences, conditional on the explicitly stated QI/SC polynomial identities, which the retained full audit had inspected. Hence a per-source/test separated tensor SOS is obtained. **Neither the root library nor its 24-by-24 parameter-dependent Gram is an evaluated uniform affine source pencil.** It cannot be used to claim a global finite SDP lift.

Complementary UL tests use the verified exact reflected parent and **only the unbounded-`a_s`** `\bar k=1` version of R2O10/11 from Part II. Literal axes and infinity are original homogeneous-parent charts, not limits of the toric local roots.

## V. No-circularity and status ledger

* Allowed common algebraic inputs: original native SG/EN source/contact polynomial identities; finite SC/QI/CC/SE **polynomial** conductors with their raw-support and opening hypotheses verified in each leaf; RFN18–20 **polynomial SOS** identities; bounded-positive-residue separated tensor lemma; original signed-core/BE/EN/RB/RI/RU/WA and homogeneous current identities under the displayed source/test caps. No whole C1/C2 source theorem is imported to prove itself.
* Kept literal faces: BN `D=0` or `\beta=0`; UN `D=0` (but not `\beta=0`); R2 RL/TU `n=0` and UL `n>0` by hypothesis; R34/R67 `n=0` where permitted. Other source openings are ordered-field positive, not silently set to a positive real lower bound.
* Preserved historical 22 October verdict: the 45 C1/49 C2 operative chains and earlier readbacks are **not erased** by a narrower later assignment or by compatibility labels. V04 checks only the four changed interfaces; the full qualitative universal-conjunction verdict and evaluation of a global affine SDP remain separate tasks.
* **Result:** no genuine fresh mismatch requiring an analytical repair. The sole historical `N/\kappa` error is corrected already in frozen main. No new checker is needed, since no changed identity is disputed.

## Immutable source and artifact pointers

Frozen TeX: `https://github.com/anubhav70747/q222-finite-sdp/blob/332483309586e0fec6e16fc1dd855915b648dca1/latex/certificates/c1_short.tex` and analogous `c2_short.tex`. Frozen archived proofs and exact checkers reside under `complete-proof-materials/materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/` at this commit; see `V04_FROZEN_SOURCE_MANIFEST.md` in the ZIP. The original proof texts in the user's attachment are extracted **verbatim**, without silent algebraic edits, in `V04_original_evidence/`.
