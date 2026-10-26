# K364 exact conditional pencil: short theorem

At main `3d193cc25ae2fcdfd3817c0bf5ebd54dfada7b55`, define the genuine fixed globally sound closed cone
\(K_{364}=\{c\ge0:\operatorname{coeff}[(1+v^2)F_c]=\operatorname{coeff}[Z_{6,12}^{T}PZ_{6,12}],\ P\succeq0\}\).
The ZIP materializes its rational `3250 x 9` source map with **120** nonzero entries, its entire `3250 x 66430` symmetric Gram map as **66,430 explicitly indexed single-key columns**, all 128 physical group action and target matrices, all SG12/SG13 face source/dual matrix bases, the full 3456 affine equality rows, PSD copy offsets, and the 16 event inequalities.

For each group member \(g\), set \(b_g(q)=E_{SG1}^TQ_g(1,q)\). The exact global LMI imposes \(b_g=T^Ty_g+r_g\), \(A^*y_g\succeq0\), \(r_g\ge0\); SG12 requires \(I_{12}^Tb_g=M_{12,x}^*Y_{12,g}+r_{12,g}\), \(M_{12,s}^*Y_{12,g}=0\), \(Y_{12,g}\succeq0,r_{12,g}\ge0\); SG13 likewise uses two \(3\times3\) matrices and \(\xi\)-cancellation. Add all 16 event inequalities. Exact dimensions: **897 PSD blocks, total order 51,088; 423,176 scalar variables including 8 Bell coordinates; 3,456 affine equalities**. Every block is explicitly defined by sparse rational data and copy offsets.

**Unconditional:** \(\mathcal Q_{222}\subseteq\mathcal Q_{364}^{outer}\) by soundness and attained dual representation; the explicit positive-definite source point \((1/21,1,\ldots,1)\) is serialized and checked.

**Conditional:** Assuming (H1) every genuinely physical oriented strict C1 source \(0<k\le1\) is in **K364** and (H2) every genuine strict six-event C2 source is in **K364**, then \(K364+\sigma K364\) covers strict contact rays, L15 admits their singular coefficient limits, and SG9/SG10 + SG12/SG13/events close the source cone. The full 128 group includes \(\sigma\), so its orbit sum needs **only K364 copies**, giving \(\mathcal Q_{364}^{outer}=\mathcal Q_{222}\). Neither H1 nor H2 is proved or inferred from the 196-row Z6 candidate. No Git writes.
