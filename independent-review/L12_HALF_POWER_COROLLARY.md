# L12 already admits the strict open C2 ends in the half-power library

This is a corollary of the newly published L12 theorem, not a new proof of its coefficient or PSD assertions.

Let
\[
Z_k=\{u^{i/k}v^{j/k}\psi_\ell:0\le i,j\le k,\ \ell=W,X,Y,Z\},\qquad
L_k=\{c:F_c=Z_k^TGZ_k,\ G\succeq0\}.
\]
The L12 proof uses fourteen even rows
\[
UW,\ UVX,\ VZ,\ VW,\ UZ,\ UVY,\ UVW,\ Z,\ UX,\ VX,\ UY,\ VY,\ X,\ Y
\]
and two odd rows \(\sqrt{UV}W,\sqrt{UV}Z\).
Every row has angular exponents in
\(\{0,\tfrac12,1\}^2\). Thus these sixteen functions are already a sublist of \(Z_2\), which itself embeds into \(Z_6\).

The L12 substitutions multiply each angular variable and each state coordinate by a nonzero source-only factor. A monomial row therefore becomes a nonzero source-only multiple of the same original angular/state monomial. Diagonal congruence transports its PSD Gram back into the original \(Z_2\). Its positive source-only scalar multiplier is absorbed into that Gram. No new angular denominator or fractional exponent is introduced.

Consequently L12 proves:
for fixed positive real \(a,t,\eta,\rho\) with
\[
0<\rho<t,\quad 0<t-\rho<1,\quad 0<\eta<1,\quad a^2t+\eta>1,
\]
and native parameters
\(\varepsilon=c^{-2}\), \(b=a^2t+\eta\varepsilon\), \(r=\rho\varepsilon\),
the exact full native C2 source belongs to \(L_2\) for every sufficiently small positive \(\varepsilon\). The statement is uniform on compact subsets of this strict face domain.

The other end also lies in \(L_2\). Under the exact parent inversion, multiplication of a reciprocal row by \(uv\) maps
\[
u^{i/2}v^{j/2}\ \longmapsto\
u^{(2-j)/2}v^{(2-i)/2},
\]
with the stated state permutation. This preserves the entire \(Z_2\) grid.

The accompanying checker verifies all sixteen exponent pairs and their inverted-grid membership exactly. It does not re-run L12's source matching, all-orders right inverse or graded Schur proof. Those remain the explicit input theorem and its published verification notes.

This does not contradict the known strict-native C2 source outside \(L_2\): it admits only the sufficiently near strict open-end neighborhoods. It gives no parameter-independent numerical cutoff, no uniform passage to source walls or coupled noncompact degenerations, and no arbitrary-finite-\(c\) or global-Q222 admission. Its useful refinement is that L12's local progress requires half-powers only; genuine sixth-root features are not used in these displayed certificates.

Sources: [current L12 checklist](../complete-proof-materials/CHECKLIST.md), [the original full theorem and proof](../complete-proof-materials/materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1407/FINITE_C_THEOREM_AND_PROOF.md), and [the global half-power exclusion](C2_HALF_POWER_OBSTRUCTION.md).
