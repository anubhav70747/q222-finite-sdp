# Q222: an exact, fully indexed conditional K364 global SDP pencil

**Frozen source:** main `3d193cc25ae2fcdfd3817c0bf5ebd54dfada7b55` (26 October 2026). **Status:** finite rational global **outer** SDP is materialized unconditionally. Equality with \(\mathcal Q_{222}\) is **conditional** on (H1) direct membership in K364 for every strict physical C1 source whose actual oriented native contact has \(0<k\le1\), and (H2) direct membership in K364 for every physical strict C2 six-event source in its exact native domain. Neither assertion is supplied or proved by this package. No universal Z6 membership assumption is used. The V01–V06 qualitative-existence review is a separate theorem.

## 1. Coefficient spaces and exactly specified maps

Write
\[
\bar q=(1,A_0,A_1,B_0,B_1,E_{00},E_{01},E_{10},E_{11})^T,\qquad
c=(d,\lambda,n,m,C,A,D,\alpha,\beta)^T.
\]
The physical affine functional is \(f_c(q)=c^TE^T\bar q\), where the **actual numeric rational** \(9\times9\) matrix \(E\), in the above orders, is `canonical_E_source_to_Bell_9x9` in the attached JSON. Its first column is
\((1,0,0,0,0,-1/2,1/2,-1/2,-1/2)^T\); its remaining columns are the eight vectors \(4p_{ab|xy}=1+aA_x+bB_y+abE_{xy}\) in SG1 event order. It is invertible over \(\mathbb Q\).

Let \(F_c(u,v;W,X,Y,Z)\) be the **complete** SG2 parent with all nine source coefficients and the signed CHSH core. In the fixed feature basis
\[
Z_{364}=(u^{i/6}v^{j/6}\psi_\ell)_{0\le i\le6,\,0\le j\le12,\,\ell=0,1,2,3},\qquad \psi=(W,X,Y,Z),
\]
use exactly the 3,250 keys
\[
\mathsf S=\{(I,J,a,b):0\le I\le12,0\le J\le24,0\le a\le b\le3\},\quad
m_\kappa=u^{I/6}v^{J/6}\psi_a\psi_b.
\]
All indices are **zero-based**, sorted first by \(I\), then \(J\), then the lexicographic ten unordered state pairs. Every list and index is materialized in `K364_FULL128_FACTORED_SPARSE_PENCIL.json`.

The matrix \(T\in\mathbb Q^{3250\times9}\) is stored as nine sparse *columns*, `T_source_columns`, with **120 nonzero entries**, and independently verified to be
\[
(Tc)_\kappa=[m_\kappa]\big((1+v^2)F_c\big).
\]
No source native parameters occur in \(T\). In particular the factor \(1+v^2\) is retained; it cannot be dropped after a reflection.

The exact map \(\mathsf A:\mathrm{Sym}_{364}\to\mathbb R^{3250}\) is stored as **all 66,430** upper-triangular independent Gram-variable columns. For a feature row \(a=(i,j,\ell)\) and row \(b=(i',j',\ell')\), let
\[
k(a,b)=(i+i',j+j',\min(\ell,\ell'),\max(\ell,\ell')).
\]
Its upper Gram coordinate \(P_{ab}\), \(a\le b\), contributes **1** if \(a=b\) and **2** if \(a<b\) to the key \(k(a,b)\). Accordingly the exact rational adjoint has the especially simple materialized upper entries
\[
\boxed{\big(\mathsf A^*y\big)_{ab}=y_{k(a,b)},\quad0\le a\le b<364.}
\]
`Gram_variable_upper_entries_map_A` and `Gram_dual_moment_upper_entries_map_A_star` list **every** key, row pair and multiplicity. The entries of \(\mathsf A^*y\) at \((b,a)\) are its symmetric copies. No placeholder map is left unspecified.

The canonical source cone is exactly
\[
K_{364}=\big\{c\ge0:\exists P\succeq0,\ Tc=\mathsf A(P)\big\}.
\]
Its unconditional global soundness follows by dividing the positive certificate \((1+v^2)F_c=Z_{364}^TPZ_{364}\) by the **pure testing** positive \(1+v^2\). The original signed-angle congruences cover both signs; polynomial continuity and SG2 homogeneous normalization cover literal axes and projective infinity, and real/imaginary state splitting covers complex states. Thus no source-dependent mixed multiplier or state restriction is used.

A strictly feasible source point is **explicitly constructed and serialized** as `K364_STRICT_RATIONAL_SOURCE_AND_GRAM.json`: \(c_*=(1/21,1,1,1,1,1,1,1,1)\), with \(P_{364,*}\succ0\) and exact \(Tc_*=\mathsf A(P_{364,*})\). The 364-Gram has 1,334 nonzero upper entries. It comes from the published 196-row rational event-root Gram, positive \(7\times7\) univariate matrix \(H\), and its two shifted embeddings \(Z_6,vZ_6\), whose row ranges jointly span all 364 features. The checker verifies the complete parent identity and exact positive principal minors; no solver is called. The moment-integral trace proof independently closes \(K_{364}\) in source coefficient space.

## 2. The **full 128** physical relabeling group, not the SG8 stabilizer

A group element is an ordered seven-tuple \((s,p_A,p_B,a_0,a_1,b_0,b_1)\), where \(s,p_A,p_B\in\{0,1\}\) and all four final outcome signs are \(\pm1\). First exchange parties when \(s=1\), then exchange settings when \(p_A,p_B=1\), and then apply the four outcome signs. Its exact \(9\times9\) signed permutation matrix \(Q_g\) maps homogeneous behaviours \(\bar q\mapsto Q_g\bar q\). For example when \(s=0\),
\[
A'_x=a_x A_{x\oplus p_A},\ B'_y=b_y B_{y\oplus p_B},\quad
E'_{xy}=a_xb_yE_{x\oplus p_A,y\oplus p_B};
\]
when \(s=1\), interchange the two original parties on the right. The JSON stores the explicit signed row permutation for each of the 128 actions, its complete canonical-target matrix, and its canonical-source matrix
\[
S_g=E^{-1}Q_g^TE,\qquad ES_g=Q_g^TE.
\]
The required actual target for orbit copy \(g\) is
\[
\boxed{b_g(q)=E^TQ_g\bar q=S_g^T E^T\bar q.}
\]
`groups` stores all 128 matrices, `group_targets` the three explicit target matrices for each group. There are exactly **eight distinct PR correlations**, each with a 16-element stabilizer, listed along with its group IDs in `PR_sectors`. The identity action has index **15**, and the verified C1 source reflection \(\sigma\) has index **24**. The latter has source permutation
\((d,\lambda,n,m,C,A,D,\alpha,\beta)\mapsto(d,n,\lambda,C,m,\alpha,\beta,A,D)\).

**Important:** for an action moving the reference PR vertex, the canonical source-coordinate matrix \(S_g\) need not preserve \(\mathbb R_+^9\). The global constraints act on **the entire source cone's Bell-coefficient image**: evaluate the unaltered canonical source-dual cone at the transformed target \(b_g(q)\). There is no illegitimate orthant projection after relabeling.

## 3. The two complete closed-face source pencils and their target injections

The original SG12 closed face omits \(\lambda,m\), so its seven-variable source injection \(I_{12}:\mathbb R^7\to\mathbb R^9\) retains indices
\[
\boxed{I_{12}=(0,2,4,5,6,7,8).}
\]
The homogeneous **source** pencil is the original explicit \(6\times6\) matrix
\[
M_{12}(x,s)=\begin{pmatrix}K_0(x)&K_1(x)/2+S(s)\\K_1(x)/2-S(s)&K_2(x)\end{pmatrix}\succeq0,
\]
where \(S\) is \(3\times3\) skew with three free coordinates and \(K_0,K_1,K_2\) are precisely the diagonal/off-diagonal matrices printed in SG12. All seven exact \(6\times6\) rational coefficient matrices \(M_{12,j}\), all three exact skew bases \(N_{12,t}\), and the \(9\times7\) injection are in `face_templates.SG12`. Every upper entry is stored as `[row,col,"numerator/denominator"]`.

The SG13 closed face omits \(\lambda,\alpha\), with injection
\[
\boxed{I_{13}=(0,2,3,4,5,6,8).}
\]
Its original source pencil has two \(3\times3\) PSD blocks \(B_1(x,\xi),B_2(x,\xi)\), one common free scalar \(\xi\), and seven nonnegative retained weights. All fourteen source basis matrices and both \(\xi\) bases are in `face_templates.SG13`, with the same exact rational upper-entry convention. They reproduce SG13 including the opposite signs of \(\xi\).

The dual constraints for orbit \(g\) are **fully expanded** by `assemble_global_pencil.py`:
\[
\begin{aligned}
I_{12}^Tb_g(q)&=\big(\operatorname{tr}(M_{12,j}Y_{12,g})\big)_{j=1}^7+r_{12,g},\\
\operatorname{tr}(N_{12,t}Y_{12,g})&=0\quad(t=1,2,3),\\
Y_{12,g}&\succeq0,\quad r_{12,g}\ge0,\\
I_{13}^Tb_g(q)&=\big(\operatorname{tr}(B_{1,j}Y_{13,g}^{(1)})+\operatorname{tr}(B_{2,j}Y_{13,g}^{(2)})\big)_{j=1}^7+r_{13,g},\\
\operatorname{tr}(B_{1,\xi}Y_{13,g}^{(1)})+\operatorname{tr}(B_{2,\xi}Y_{13,g}^{(2)})&=0,\\
Y_{13,g}^{(1)},Y_{13,g}^{(2)}&\succeq0,\quad r_{13,g}\ge0.
\end{aligned}
\]
In independent symmetric upper-triangular variables, every matrix trace uses factor **two** for off-diagonal entries. The corresponding **17** affine face equalities per group are serialized individually, with their absolute scalar-variable indices. The omitted coordinates have **no** incorrect additional target signs or constraints.

Source Slater witnesses, verified exactly: SG12 uses \(d=n=D=\alpha=\beta=1,C=A=2,S=0\), with leading minors \((14,112,448,2576,28000,62500)>0\); SG13 uses \(d=m=D=1,n=C=A=\beta=2,\xi=-5\), and its two blocks have minors \((4,72,234)\) and \((22,72,234)\). Each seven-variable source orthant is strictly positive as well. These strict source points prove attainment of the displayed face duals by finite-dimensional separation; no unproved abstract PSD projection-closedness is used.

## 4. The actual 128-copy global LMI and copy/offset tables

For every \(g=0,\dots,127\), take independent auxiliary variables
\[
y_g\in\mathbb R^{3250},\quad r_g\in\mathbb R_+^9,\quad
Y_{12,g}\in\mathbb S_+^6,\ r_{12,g}\in\mathbb R_+^7,\quad
Y_{13,g}^{(1)},Y_{13,g}^{(2)}\in\mathbb S_+^3,\ r_{13,g}\in\mathbb R_+^7.
\]
Write \(M_g(y_g)=\mathsf A^*y_g\), with the **66,430 actual upper entries** from §1. Impose
\[
\boxed{b_g(q)=T^Ty_g+r_g,\quad M_g(y_g)\succeq0,\quad\operatorname{diag}(r_g)\succeq0,}
\]
and all the SG12/SG13 equations/PSD inequalities in §3.

The last required LMI is **all 16 event inequalities**, not an eight-event subset:
\[
\boxed{\operatorname{diag}_{x,y\in\{0,1\},a,b=\pm1}
\big(1+a A_x+b B_y+ab E_{xy}\big)\succeq0.}
\]
Their exact 16 affine rows are stored in `physical_event_inequalities`, and their absolute matrix positions in `K364_ALL_EVENT_LMI_ENTRIES.json`.

The **single block-diagonal global LMI** is indexed as follows. Each of 128 copies contributes the seven blocks

| Per-copy block | PSD order | Independent auxiliary scalars |
|---|---:|---:|
| \(\mathsf A^*y_g\) | 364 | 3,250 |
| \(\operatorname{diag}(r_g)\) | 9 | 9 |
| \(Y_{12,g}\) | 6 | 21 |
| \(\operatorname{diag}(r_{12,g})\) | 7 | 7 |
| \(Y_{13,g}^{(1)}\) | 3 | 6 |
| \(Y_{13,g}^{(2)}\) | 3 | 6 |
| \(\operatorname{diag}(r_{13,g})\) | 7 | 7 |
| **Per-copy total** | **399** | **3,306** |

Plus one order-16 diagonal event block. The 128 copies have **27 explicit affine equalities per copy**: 9 K adjoints, 7 SG12 target adjoints, 3 SG12 skew constraints, 7 SG13 target adjoints, 1 SG13 \(\xi\) constraint. Totals:

\[
\boxed{\begin{gathered}
\text{PSD blocks}=897,\quad\text{PSD total order}=51{,}088,\\
\text{auxiliary scalar variables}=423{,}168,\quad
\text{Bell variables}=8,\quad\text{total}=423{,}176,\\
\text{affine equations}=3{,}456,\quad\text{event inequalities}=16.
\end{gathered}}
\]

`group_copy_offsets` contains all **128 genuine block offsets**, variable offsets and equality offsets. `K364_GLOBAL_ALL_3456_AFFINE_EQUALITIES.json` contains **every complete affine row** in the convention `constant + Σ entries[idx]*variable[idx] = 0`, with the physical eight \(q\) variables occupying indices 0–7. `Gram_dual_moment_upper_entries_map_A_star` and the face upper-basis matrices **are the actual sparse LMI matrices**, not undefined operator names. `assemble_global_pencil.py` deterministically substitutes offsets to generate any complete copy or, with `--stream-all-lmi`, to expand **all 8,510,224 upper-triangular PSD entries** as a compressed JSONL file. The checked factored data already specify every entry exactly, so the enormous repeated full expansion is optional. The only encodings are integers or exact rational strings `numerator/denominator`.

## 5. Both inclusions, with different logical statuses

**Unconditional outer inclusion:** Every \(c\in K_{364}\) is a globally valid Bell inequality by its full-state Gram and pure-testing multiplier. All SG12/SG13 face sources are globally valid by their published exact source pencils; every event inequality is valid. Physical relabelings preserve validity. Consequently every \(q\in\mathcal Q_{222}\) satisfies every target dual membership condition and the complete global LMI. All dual representations are **attained**: the strict \(K_{364}\) witness, the strict face source witnesses and closed Gram-image trace bounds justify exactly the displayed equations, orthants, PSD matrices and face auxiliary cancellations. Thus
\[
\boxed{\mathcal Q_{222}\subseteq\mathcal Q_{364}^{\rm outer}\quad\text{unconditionally}.}
\]

**Conditional reverse inclusion:** Suppose (H1) and (H2). Every oriented strict C1 source lies in \(K_{364}\); every opposite strict C1 source belongs to \(\sigma K_{364}\), by the verified native transport and its positive source normalization; every strict C2 source lies in \(K_{364}\). The finite sum \(\bar K=K_{364}+\sigma K_{364}\) is closed by the verified deterministic coefficient bound for **sound** Bell sources, not by a general theorem on sums of closed cones. The published L15 source-coefficient density theorem places the allowed third-omission contacts in \(\bar K\). The original SG9/SG10 extreme-source reduction puts every remaining valid extreme ray in \(\bar K\), one of the two SG12/SG13 closed-face cones, or an event ray, up to physical symmetry. The exposed-ray theorem ensures that this is actual individual extreme-source admission, not a spurious convex decomposition between different contact rays. Since \(\sigma\) is group element 24,
\[
\sum_{g\in\Gamma}g\bar K
=\sum_{g\in\Gamma}gK_{364},\qquad |\Gamma|=128.
\]
The finite sound cone sum is closed by the same deterministic coefficient bound. It equals the entire Bell-validity cone. Taking duals gives
\[
\boxed{(\mathrm{H1})+(\mathrm{H2})\Longrightarrow
\mathcal Q_{364}^{\rm outer}=\mathcal Q_{222}.}
\]
There is **no claim** that H1/H2 hold, that \(\sigma K_{364}=K_{364}\), or that all nonoriented strict sources lie individually in \(K_{364}\). There is also no claim that the older original-196-row \(Z_6\) candidate is exact.

## 6. Exact verification and reproduction

The package contains:

- `generate_exact_pencil.py`: produces the complete source matrix \(T\), the **entire** fixed map \(\mathsf A\) and its adjoint, all 128 actions, every target/injection, both SG12/SG13 primal source base matrices and dual face constraints, 16 event rows, and all copy offsets, as one rational JSON.
- `assemble_global_pencil.py`: produces all 3,456 literal affine rows and the 16 event-block entries; `--stream-all-lmi` expands the factored PSD matrices entrywise, without changing the mathematical model.
- `verify_exact_pencil.py`: reconstructs the nine SG2 target columns independently with symbolic algebra, checks all 66,430 Gram entries/adjoints, the group closure/PR orbit/sigma, all face matrices and Slater witnesses, the full 3,456 assembled affine rows, q=0 feasibility, and the strictly positive rational source Gram.
- `K364_FULL128_FACTORED_SPARSE_PENCIL.json` and `.json.gz`: all rational matrix bases, exactly indexed.
- `K364_GLOBAL_ALL_3456_AFFINE_EQUALITIES.json` and `.json.gz`: every equality, no symbolic placeholders.
- `K364_STRICT_RATIONAL_SOURCE_AND_GRAM.json`: 196/364 upper-triangle rational witnesses, exact Qc and H, and their rational positivity minors.
- `K364_ALL_EVENT_LMI_ENTRIES.json`: actual 16 diagonal LMI entries.
- `EXACT_VERIFICATION_REPORT.json`, `GENERATION_LOG.txt`, `ASSEMBLY_LOG.txt`, `EXACT_CHECKS.log`: reproducibility records. `requirements.txt`: SymPy, standard Python library otherwise.

Run in the extracted package directory:

```bash
python generate_exact_pencil.py
python assemble_global_pencil.py --check-structural
python verify_exact_pencil.py
```

All checks use integer/rational identities and finite Sylvester determinants. No solver, numerical scan, approximate eigenvalue test, tool-driven search, or repository write occurs.
