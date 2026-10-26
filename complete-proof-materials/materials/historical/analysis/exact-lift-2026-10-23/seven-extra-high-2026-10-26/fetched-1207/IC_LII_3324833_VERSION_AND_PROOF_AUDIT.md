# Q222: IC/LII frozen-version and full-parent interface audit

**Frozen version:** `332483309586e0fec6e16fc1dd855915b648dca1`. **Scope:** the physical native C2 weak-IH (IC) and strong-II (LII) source cells only. The current source is `latex/certificates/c2_short.tex`, component `pc00b84929b12` (IC), endpoint component `p18dbedfff317` (WIE), remote `pf95b1c5fc32f`, `p0fb32cd2f1bc`, and `pe16499ccc4c6` (LII), remote `p590112f102bd`. The recovered proof `Q222_C2_IC_LII_Full_Parent_Proof_2026-10-26.md` and its eight checker files are the exact original Git blobs, verified by Git's blob SHA-1, not newly transcribed substitute proofs. The earlier complete C2 review (`final_audit.txt`, 22 October) is preserved as a **previously reviewed import**, not silently upgraded into a global theorem.

## Verdict and quantified statement

**Verdict (a):** no mathematical discrepancy or newly uncovered admission gap was found in the compared IC/LII formulations. The current source and original recovered proof compose, subject to the explicit native-contact and strict bounded-tensor-Gram lemmas already reviewed in their own proofs. Precisely, for arbitrary separately labelled real-closed extensions \(K_c,K_t/\mathbb R\), a physical native C2 source in either of the following domains and every independent real or complex projective test have a complete-parent tensor SOS in \(K_c\otimes_{\mathbb R}K_t\):

* **IC weak-IH:** \(z=u_c^2\) source-infinitesimal, \(y=v_c^2\) source-unbounded, \(B=H+L\) source-unbounded, original native strict gaps, \(C\) weak, \(n\ge0\), source-shape \(\rho=B/y\) with residue in \([0,1]\), including \(\rho_0=0,1\), \(H\to0,1\), and literal \(n=0\).
* **LII strong-II:** \(z,y\) source-infinitesimal, \(B=H+L\) a **positive bounded unit**, original native strict gaps, \(n\ge0\), \(m,C,A,D,\alpha>0\), with source openings \(r,t=L/y,j=1+zr-t,1-H\) arbitrarily small, including literal \(n=0\).

Neither cell is the complete C2 cone. There is no conclusion here about full qualitative finite-existence assembly or any explicit uniform rational/sixth-root lift.

## Equation-level version crosswalk

| Recovered original proof | Frozen current native TeX | Audited equation-level fact |
|---|---|---|
| (1), (2) | `eq:c2-master-parent`, IC3 | The original full source parent \(F=4P\), versus the reciprocal-Bob homogeneous \(\widehat F=w^2F=4w^2P(u,1/w)\). IC derivatives are defined as \(4L_i=\partial_i\widehat F\). LII works throughout at \(F=4P\). No unnoticed quarter-parent factors. |
| (3), (5)–(6), (24) | WI1–WI3, GE1–GE2, DI7, LII1 | Put \(R=zr\), \(L=yt\), \(S=1+R+L\). WI and GE are the **same original native six-weight dictionary** in these variables; they have different, noninterchangeable source-scale conclusions. |
| (7)–(18) | IC4–IC12 | Graph \(X=-ruZ,Y=HwZ\); five contact roots \(R_0,R_1,R_2,T,E\); complete graph identity IC7; both independent normal rows IC10–IC12 and their complete hard block. |
| (19)–(20) | WIE1–WIE2 | \(\rho_0=0\) requires \(s=r/\rho\) bounded and replacement \(T_\rho=\rho T\); bounded first-normal coquotient \(f_X/\rho\) and strict \((R_2,T_\rho)\) block. |
| (21)–(23) | WIE3–WIE6 | \(\rho_0=1\) requires the entire \(4D(J+eZ)^2\) packet, with \(e=1+zr-\rho+H/y\), be **removed before** bounded-Gram admission. WIE4 and WIE5 are exact identities for the graph and full state, respectively. |
| (25)–(28) | LII2–LII6 | Four-state triangular graph, C/D harmonic completion, and **exactly four** retained nonnegative packets. Their entire remainder is LII6. |
| (29)–(31) | LII7–LII8 | The sole formal \(R_hW\) coefficient is \(uk\mathcal R(U,V)\); \(KR_h=uvWQ_n-sU(V-y)W^2\) is an exact identity on the dependent root; substitution removes that coefficient without dividing by a mixed element. |
| (32)–(35) | LII9–LII16 | Two full state-normal rows, bounded monic contact quotients, scalar support and allocation on six roots. |
| (36)–(37) | LII17–LII20 | Strict residue block and the exact free \(T^2\) pivot from LII20. |
| IV | OA3–OA9, IR1–IR8, RR1–RR13 | Distinct **original homogeneous-parent** outer/testing frames; not limits of raised-degree contact roots. |

## Three structural verifications

### IC: source bounds, full normals, and factor discipline

The WI1 source formulas imply \(BC,Bm/z^2,AH/B,\alpha(1-H)/B\) are positive bounded units and \(nr\) bounded for the full weak-IH source. In the interior \(0<\rho_0<1\), \(D\) is bounded, \(N=4n+2\ge2\), and the source-normalized roots are
\[
E_X=\sqrt N\,(X+ruZ),\qquad
E_Y=\sqrt{B/[H(1-H)]}\,(Y-HwZ).
\]
The graph identity IC7 is determined *exactly* by its four-monomial residual support: the difference of the original graph and candidate equals \(AuwWZ+(B_1u^2+B_2w^2+B_3)Z^2\). Every candidate root vanishes at the native source contact. The \(W\), \(u\), \(w\) derivatives and then contact value force \(A=B_1=B_2=B_3=0\), with no source/test identification in the certificate. IC8 and IC8a are true polynomial zero identities. Their second identity supplies the strict radial \(+\rho T^2\) entry.

For the **whole parent**, set \(4L_X=\partial_X\widehat F|_{X=-ruZ,Y=HwZ}\) and similarly \(4L_Y=\partial_Y\widehat F|_{\mathrm{graph}}\); the complete linear terms are \(4L_X E_X/\sqrt N\) and \(4L_Y E_Y/\sqrt{B/[H(1-H)]}\), not the half derivatives of \(P\). IC10–IC12 write these exact linear forms on the same five roots. In particular, monic source/test-polynomial divisions show
\[
P_X=(U-z)f_X+(V-1/y)g_X,\quad
P_Y=(U-z)f_Y+(V-1/y)g_Y,
\]
with zero contact remainders. All explicit multipliers in IC11–IC12 have zero residue: \(m\sqrt B/z=(Bm/z^2)z/\sqrt B\to0\), \(f_Y/B,g_Y/B\) are bounded at both \(H\)-openings, and the remaining terms contain independent testing \(u,w\) or a small source factor. The complete quadratic hard block has residue \(\operatorname{diag}(1,2/(1+\rho_0))\), not a scalar determinant substitute. The scalar \((R_2,T)\) block in IC9 has strictly positive determinant by \(0\le r_0\le\rho_0/(1+\rho_0)\). This establishes one bounded strict **seven-root coefficient matrix**.

For \(\rho_0=0\), the source regularization uses \(\tau=1/B\), \(s=r/\rho\), \(y=(\tau\rho)^{-1}\). Normalized WI1 rational expressions give \(k/(B\rho)\to2s_0\), \(nr/\rho\to(1-s_0)/2\), \(D/\rho^2\to1/2\) simultaneously, without a forbidden ratio such as \(z/\rho\). Scaling \(T_\rho=\rho T\) makes \(f_X/\rho=-2(nr/\rho)(1+V)-s(2D+1)\) bounded. WIE2 adds a unit \(T_\rho^2\) pivot and yields \(\begin{psmallmatrix}2&-s_0\\-s_0&1+2s_0\end{psmallmatrix}\succ0\).

For \(\rho_0=1\), \(D\) can be unbounded. The graph identity WIE4 is the exact subtraction of \(4D Q^2\), \(Q=R_2-rT-\sqrt{H/B}E=J+eZ\) on the graph. A concise direct check is to collect all \(D\)-terms of IC7: the \(k\)-contributions combine by IC8 into \(8DHr(U-z)(V-1/y)Z^2\), exactly the cross term \(8Dr\sqrt{H/B}\,TE\) of \(4DQ^2\). This is the otherwise missing WIE4 cancellation. In the **full** state, WIE5 expresses \(R_2\) through \(F_D=\sqrt D(J+eZ)\), \(T,E,E_X,E_Y\) with bounded coefficients. The remaining normal arrays are the D-removed IC11–IC12 with constant correction \(De\), not bare \(D\), and \(De\to1/4\) by WIE3. The complete seven-row residue is WIE6: \(\operatorname{diag}(1,1-r_0,1+r_0-r_0^2,1,4,1,1)\succ0\), \(0\le r_0\le1/2\). All original source event boundaries, including \(n=0\), remain admissible.

### LII: bounded packet subtraction and actual six-root Gram

Write \(c=CL\), \(d=Dyj\), \(t=L/y\), \(j=1+zr-t\), \(s=(1+zr)/B\), \(q=D/(C+D)\), \(g=yCD/(C+D)=cd/(cj+dt)\). From DI7 and LII1, \(B,s,c,d,g,m/z^2,A/y\) are positive bounded units, \(nr,\alpha(1-H)\) bounded, \(t+j=1+zr\) a unit. Thus \(cj+dt\) is a unit even when \(t\) or \(j\) is infinitesimal, and \(\ell=y(C+D)=c/t+d/j\) is bounded below by a positive real (possibly unbounded). Choose \(0<\ell_*\le\ell\) and \(0<\alpha_*\le\alpha\) real.

LII4 follows by expanding
\[
Q/k=(t/B)K-N-qR_h,\quad
J/k=(1-q)R_h-N-(j/B)K,\quad \eta=k(N+qR_h).
\]
For instance the \(KN\) coefficient is \(8(d-c)/B\), the \(R_h^2\) coefficient is \(4g\), and the \(K^2\) coefficient is \(4(ct+dj)/B^2\). The \(4CVQ^2\) identity must be included separately; otherwise its \(V\)-dependent normal row is missing. The full remainder LII6 subtracts exactly
\[
4(\ell-\ell_*)N^2+4CV\eta^2+4nD_uD_v\xi^2
 +4(\alpha-\alpha_*)D_u(v\eta-\zeta)^2,
\]
all whole nonnegative squares, zero in value and first jets at the native contact. Expanding the \(n\)- and \(\alpha\)-packets leaves precisely \(nr\), \(n z^2r^2\), and \(\alpha(1-H)\), \(\alpha(1-H)^2\); there is no naked unbounded \(n\) or \(\alpha\) coefficient. The remaining \(m,A\) packets have \(m=z^2\widehat m\), \(A=y\widehat A\) with bounded \(\widehat m,\widehat A\).

The only formal \(R_hW\) term has coefficient \(uk\mathcal R(U,V)\), and LII7 replaces \(KR_h\) by an exact polynomial without \(R_hW\). The printed coefficient \(\mathcal R\) is bounded term by term using \(c,d,g,q,s\), \(nr\), and \(\alpha(1-H)\). Coefficient extraction then gives LII8, containing the entire quadratic hard block \(Q_0(N,\xi,\zeta,R_h)\), both linear state-normal polynomials, and scalar \(\phi(U,V)\). There is no state restriction hidden here.

The derivative zeros and parity yield LII9–10, \(P_\xi\in(U,z)\), \(\deg_U\phi\le2\), and \(\phi(0,V)=V(d_s+e_sV)\), where \(d_s/z^2\), \(e_s/z^2\) are bounded by DI7. The bounded coefficient-divisibility mechanism is:

1. Divide \(P_\xi\) by the monic \(U-z\); its residue at \(U=z\) contains an exact \(z\) factor because \(P_\xi\in(U,z)\) (including the formal \(-8m=-8z^2\widehat m\)). Then divide the bounded quotient by \(V-y\). The other state-normal rows require only the two plain monic divisions.
2. Subtract \(V(d_s+e_sV)(U-z)^2/z^2\) from \(\phi\). The remainder is \(U[a(V)U+b(V)]\). Its value and \(U\)-derivative at \((z,y)\) give \(a(y)=b(y)=0\), so \(a=(V-y)A_s\), \(b=(V-y)B_s\). The \(V\)-jet yields \(B_s(y)+zA_s(y)=0\), making \(C_s=(B_s+zA_s)/(V-y)\) a **bounded monic quotient**. This is LII13.
3. The five explicit polynomial identities LII15 convert the resulting \(U-z,V-y\) terms into a bounded quadratic on \((N,\xi,\zeta,R_h,R,T)\). No mixed inverse, no source/test rate comparison, no discarded normal.

The hard residue LII17 has \(\operatorname{diag}(4\ell_*,2,4\alpha_*,4g_0)\) on the four original roots. The only extra nontrivial finite block before the last \(T\)-row is \(\begin{psmallmatrix}2&2/[h(1+h)]\\2/[h(1+h)]&2/[h^2(1+h)]\end{psmallmatrix}\), whose Schur pivot is \(2/[h(1+h)^2]>0\). This factor of **two** comes from the \(F=4P\) normalization and the full \(\xi W\) coefficient; it is not the half-scaled \(P\) coefficient. LII20 is an exact zero identity; its only nonzero residue entry is \(+\Lambda T^2\), so sufficiently large fixed real \(\Lambda\) absorbs all bounded couplings incident to \(T\). This yields the strict six-row Gram for the **entire** original parent after restoration.

### Complementary tests and exact imports

* **IC:** The \(u,w\) double-infinitesimal chart is IC/WIE. The open \(u\)-unit, \(w\)-infinitesimal strip is the separately reviewed **OA3–OA9** exact original-parent Schur argument, including \(\rho_0=0,1\). Finite Bob and both reciprocal extremes are **IR1–IR8**, with original \(A,\alpha\) rows retained and the \((P,E_X)\) exceptional block determinant exactly 8. Source-entry hypotheses for these joins follow from WI2–WI3, not from an assumed C2 global theorem.
* **LII:** The original homogeneous parent RR3 supplies all testing \(v\) with noninfinitesimal projective sine coordinate. For \(v\) infinitesimal and \(u\) noninfinitesimal, RR5 is the **exact harmonic identity** \(CD_vE^2+D(E+vK)^2=N^2+\gamma V^2+Cv^2E^2\), where \(\gamma=g\) is a unit. The exact toric syzygy RR6 permits a finite real choice of the cross coefficient; RR10–RR12 give strict finite blocks if Alice has finite nonzero residue, while RR13 is the direct three-square endpoint when Alice is at infinity. All retained \(n,m,A,\alpha\) event terms are restored.

The remote chains above are **already-reviewed imports**, with their precise source bounds and relevant complete-parent identities checked for applicability; they are not new all-C2 assertions. Both sectors include sign relabelings and real/imaginary splitting. Literal axes are covered by polynomial identities, infinity by original homogeneous parents, never by specializing raised-degree local roots.

## Integrity and reproducibility limits

The frozen Git tree contains `PROOF.md`, `README.txt`, `REPLAY_OUTPUT.txt`, `run_checks.py`, and **eight** original `Q222_*.py` checker modules in the original IC_LII package. Their exact file bytes are restored from the supplied combined-source attachment and their **Git blob SHA-1** hashes agree with the pinned repository blobs. The reproduced source excerpts reflect current `c2_short.tex` at the pinned commit. The historical replay log records eight passes; no original checker was rerun here because comparison disclosed **no disputed equation**, as requested. A checker validates specified identities but cannot prove boundedness, positivity or exhaustive routing: those arguments are provided above.

**Status:** composable IC/LII full-parent local and all-test certificates **within these exact native domains and reviewed imports**. No additional defect has been inferred from another conversation's filename or absence. The qualitative global Scheiderer criterion, all other actual local leaves and the explicit fixed-size global lift remain separate obligations.
