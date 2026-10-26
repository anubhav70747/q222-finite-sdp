# An exact order-13 pencil for the nine-term angular Laurent cone

**Research note, 24 October 2026.** This result concerns the *angular determinant-coefficient cone*, not an affine representation of the full eight-coordinate Bell set. It sharpens the order-49 angular SOCP in [the Q222 repository](https://github.com/anubhav70747/q222-finite-sdp/tree/main/extracted-lifts/global-weighted/angular-order49). The quartic coefficient map from canonical Bell sources is not linearized here.

## Statement and evaluated pencil

Let \(E=\{-1,0,1\}^2\setminus\{(0,0)\}\), and write

\[
 L(U,V)=h+\sum_{(i,j)\in E}q_{ij}U^iV^j,\qquad q_{ij}\ge0,
 \qquad U,V>0.
\]

For the matrix description, index the nine coefficients of \(G=UVL\) by \(g_{ab}\), \(0\le a,b\le2\): set \(g_{11}=h\) and \(g_{ab}=q_{a-1,b-1}\) otherwise. Let

\[
R=\operatorname{diag}(g_{00},g_{01},g_{02}),\qquad
S=\operatorname{diag}(g_{20},g_{21},g_{22}).
\]

Choose the following thirteen free real variables: the nine entries of \(J\in\mathbb R^{3\times3}\), the three entries of a symmetric \(C=\begin{psmallmatrix}c_{00}&c_{01}\\c_{01}&c_{11}\end{psmallmatrix}\), and \(\rho\in\mathbb R\). Define the symmetric \(3\times3\) matrix

\[
B=\begin{pmatrix}
g_{10}-2J_{00}&-J_{01}-J_{10}-c_{00}/2&\rho\\
*&h-2(J_{02}+J_{11}+J_{20})-2\rho-2c_{01}&-J_{12}-J_{21}-c_{11}/2\\
*&*&g_{12}-2J_{22}
\end{pmatrix}.
\]

Then \(L(U,V)\ge0\) for all \(U,V>0\), **with its eight outer coefficients nonnegative**, if and only if the following **single affine pencil** is feasible:

\[
\boxed{\operatorname{diag}\!\left(
\begin{pmatrix}R&J\\J^T&S\end{pmatrix},
B,C,[g_{10}],[g_{12}]
\right)\succeq0.}
\tag{A13}
\]

The PSD block orders are \((6,3,2,1,1)\), total order **13**. Its coefficient entries are rational and affine in the nine visible Laurent coefficients and thirteen auxiliaries. For canonical Bell sources \(c\ge0\), all eight outer coefficients are already nonnegative, so the last two scalar blocks can be omitted: the relevant angular feasibility system then has total order **11**, conditional on that source orthant.

## Exact coefficient identity

Set \(t=\sqrt V\), \(v=(1,t,t^2)^T\), \(w=(1,t)^T\). With \(M=\begin{psmallmatrix}R&J\\J^T&S\end{psmallmatrix}\), the polynomial identity encoded by the displayed five entries of \(B\) is

\[
\begin{aligned}
UVL(U,V)&=(v,Uv)^T M(v,Uv)
  +U v^TBv+Ut\,w^TCw.\tag{I}
\end{aligned}
\]

Here the notation \((v,Uv)\) means vertical concatenation. Indeed the constant and quadratic coefficients of \(U\) are \(\sum_{j=0}^2g_{0j}t^{2j}\) and \(\sum_{j=0}^2g_{2j}t^{2j}\). The coefficient of \(U\) is \(g_{10}+ht^2+g_{12}t^4\); its coefficients at \(t^0,t^1,t^2,t^3,t^4\) give precisely the five linear equations solved in the definition of \(B\).

If (A13) holds, each term on the right of (I) is nonnegative for \(U,t>0\), hence \(L\ge0\). The two scalar blocks enforce the remaining coefficient signs \(g_{10},g_{12}\ge0\); the other six signs are implied by the diagonal blocks \(R,S\) of \(M\).

## Necessity using only twelve circuits

The repository establishes the elementary balance reduction: a nonnegative Laurent polynomial with one potentially negative constant and nonnegative outer coefficients decomposes, after normalization at a positive minimum and with a limiting argument for missing coefficients, into positive monomials, a nonnegative constant, and the twelve primitive balanced circuits of \(E\). The circuits comprise four opposite pairs, four triangles of type \((1,1,1)\), and four of type \((1,1,2)\).

For an opposite pair \(e,-e\), the complete cleared circuit is one square:

\[
UV\big(U^{e_1}V^{e_2}+U^{-e_1}V^{-e_2}-2\big)
=\big(U^{(1+e_1)/2}V^{(1+e_2)/2}-U^{(1-e_1)/2}V^{(1-e_2)/2}\big)^2.
\]

For the first triangle type, the representative is

\[
\begin{split}
U^2V^2+U+V-3UV
={}&(UV-\sqrt V)^2+U(\sqrt V-1)^2\\
&+2U\sqrt V(\sqrt V-1)^2.
\end{split}\tag{T3}
\]

The representative of the second type is

\[
U^2V^2+V^2+2U-4UV
=V^2(U-1)^2+2U(V-1)^2.\tag{T4}
\]

All other triangle orientations follow on replacing \(U,V\) by \(U^{\pm1},V^{\pm1}\), then multiplying by the positive monomial square that changes the normal-form multiplier to \(UV\). Crucially, every ordinary square in (T3)–(T4) is linear in the fixed feature list \((1,t,t^2,U,Ut,Ut^2)\); every square carrying a factor \(U\) is linear in \((1,t,t^2)\), and every square carrying \(Ut\) is linear in \((1,t)\). Moreover, the squares in the six-dimensional block each connect exactly one monomial with \(U\)-degree zero and one with \(U\)-degree one. Thus the two diagonal \(3\times3\) subblocks can be fixed to the positive outer coefficients \(R,S\) without losing any circuit certificate.

Multiplying the normalized identities by positive source-only weights and replacing \(U/U_*\), \(V/V_*\) by their monomials preserves this fixed testing span. Positive residual monomials are single squares in the same span. Hence every strictly supported valid Laurent polynomial has the Gram representation (I) with the fixed diagonal blocks.

For a boundary polynomial, replace each outer coefficient by \(q_e+\varepsilon\), retain \(h\), and take \(\varepsilon\downarrow0\). The resulting Gram matrices can be chosen with fixed diagonal blocks as above. Integrate (I) over any rectangle \(1\le U,t\le2\). The three moment matrices obtained from the linearly independent fixed feature vectors are positive definite. Since the integrals of the polynomials are uniformly bounded, the traces of all three PSD Gram blocks are bounded. A subsequence converges and preserves the affine coefficient equations. The limit gives (A13) even with zero outer coefficients or an infimum at infinity. Thus the cone projection is closed and the equivalence is exact.

## Two short rational certificates

For the equal-weight circuit \(UV+U^{-1}+V^{-1}-3\), the nonzero coefficients of \(G=UVL\) are \(g_{22}=g_{10}=g_{01}=1\), \(h=-3\). Choose \(J_{12}=-1\) and all other \(J\)-entries zero, \(\rho=0\), \(C=2\begin{psmallmatrix}1&-1\\-1&1\end{psmallmatrix}\). The resulting \(B=\begin{psmallmatrix}1&-1&0\\-1&1&0\\0&0&0\end{psmallmatrix}\), and all the blocks in (A13) are PSD. Identity (I) reduces to (T3).

For the weighted circuit \(UV+U^{-1}V+2V^{-1}-4\), the nonzero coefficients are \(g_{22}=g_{02}=1\), \(g_{10}=2\), \(h=-4\). Choose \(J_{22}=-1\), other \(J\)-entries zero, \(\rho=-2\), \(C=0\). Then \(B=2\begin{psmallmatrix}1&0&-1\\0&0&0\\-1&0&1\end{psmallmatrix}\succeq0\) and (I) becomes (T4).

**Executed checks:** independent exact-rational polynomial expansion verified all twelve circuit identities and their membership in the common three block spans. The polynomial products of these blocks have exactly fifteen distinct coefficient monomials; their five independent cross-coefficient equations lead to (A13). For the published C1 and C2 contact seeds, exact-rational circuit decomposition produced order-11 Gram matrices of rank six and traces \(230544\) and \(4775/384\), respectively. These checks support the identities; the necessity proof above supplies global coverage.

## Order nine and six squares on the genuine Bell contact families

The source geometry gives a further reduction which does not hold for an arbitrary nine-term Laurent polynomial. For family C1, \(\lambda=n=0\), the determinant formulas (SG6) force \(g_{00}=g_{20}=0\). For family C2, \(\lambda=\beta=0\), they force \(g_{00}=g_{22}=0\). In (A13), every zero diagonal entry of the positive-semidefinite block \(M\) annihilates its entire corresponding row and column. Delete these two zero rows and columns. The nontrivial angular Gram block sizes are now **\((4,3,2)\)**, total **order nine**, with no scalar sign blocks needed because the Bell source orthant already enforces all outer signs. In the three-event-omission subfamilies, a further outer coefficient is zero and the size drops to at most eight.

At a finite contact \((U_c,V_c)\), all three summands of (I) are nonnegative and have sum zero. Each PSD Gram therefore annihilates its own nonzero evaluation vector, \((v,U_cv)\) (restricted to four surviving coordinates), \(v\), and \(w\), respectively, at \(v=(1,\sqrt{V_c},V_c)\), \(w=(1,\sqrt{V_c})\). Their ranks are at most \(3,2,1\), respectively. **Every C1/C2 extreme source therefore admits a determinant certificate using at most six weighted squares in a fixed eleven-feature family, and with a Gram matrix of order nine.** A third omitted event improves the bound to five squares. This is a structural bound, not a numerical observation: it follows from the coefficient support, (A13), and saturation at contact. The published C1/C2 seed decompositions attain rank six.

## Relation to the full Q222 problem

The canonical source coefficients \(c=(d,\lambda,n,m,C,A,D,\alpha,\beta)\ge0\) enter \(h(c),g_{ij}(c)\) through degree-four polynomials defined by the state determinant (Chapter 4 of the repository). Substitution into (A13) creates a polynomial, **not affine**, matrix inequality in \(c\). Therefore this note does not extract the missing affine eight-coordinate quantum-set lift. It provides a smaller exact angular coefficient cone and a fixed eleven-function tensor-SOS certificate for the determinant term. The remaining global problem is to linearize/convexify this quartic source pullback, or to construct an affine source pencil by a different state-normal completion.