# Canonical extreme-source geometry for Q222: exact proof and audit

**Scope.** SG1--SG10 and EN1--EN16 in the supplied `04_geometry.tex` and `05_entry.tex`, together with the explicit native stationary identities SD1, SD9--SD11 and the corresponding native C2 formulas. This proves the classification and its algebraic entry interface, independently of any SDP-existence theorem or tensor-SOS completion. Source parameters and physical tests are not identified. The finite-order-field Carathéodory transfer used downstream is assumed as already audited.

**Reproduction.** Run `python Q222_SG_EN_geometry_exact_checker.py` and `python Q222_C2_native_stationarity_check.py`. Both require only SymPy. They verify the full generic determinant SG5--SG6, SG15/SG20, SG9 support logic, EN6/EN12, both full-parent physical reflections, EN7 Cramer's identity, EN8, EN16, both seeds, C1 SD9--SD10, C2 native contact equations, and SG12--SG13 limiting determinants/identities. The checker is an exact polynomial/rational equality test, not parameter sampling.

## Theorem: exact canonical geometric classification

Write
\[
H=1-\frac{E_{00}-E_{01}+E_{10}+E_{11}}2,
\qquad f=dH+4\sum_{i\in\{\lambda,n,m,C,A,D,\alpha,\beta\}}c_ip_i.
\]
For a quantum-valid affine functional, either it is a nonnegative linear combination of event probabilities (after allowing all eight physical PR orientations), or one can choose a PR orientation with \(d>0\) and \(c_i\ge0\). Inside that canonical cone:

1. \(d=0\) gives only event rays as extreme rays.
2. Every valid \(d>0\) source with noncoercive Laurent support lies, up to the physical dihedral relabelings of the eight event cycle, on **one of two closed faces** \(\lambda=m=0\) or \(\lambda=\alpha=0\). The *entire valid portions* of these faces, including singular boundaries and further zero coefficients, have the exact finite pencils SG12 or SG13 below.
3. Every extreme valid \(d>0\) source with coercive Laurent support has a finite, noncommuting stationary zero and has exactly five or six positive event weights. Up to those physical relabelings its support is
\[
\begin{array}{ll}
\mathrm{C1}:&\lambda=n=0,\ m,C,A,\alpha>0,\ D,\beta\ge0,\ D+\beta>0,\\
\mathrm{C2}:&\lambda=\beta=0,\ m,C,A,D,\alpha>0,\ n\ge0.
\end{array}
\tag{T}
\]
The **only** C1 third omissions are \(D=0\) and \(\beta=0\); the **only** C2 third omission is \(n=0\). The latter is a C1 third-omission orbit.
4. At every contact in (T), the physical signs EN2 hold *strictly*, including all allowed third-omission boundaries. EN13 or EN14--EN15 therefore gives positive native coordinate gaps; the native rational stationary coefficients equal the original source coefficients by the nonzero determinants EN8 or EN16. Conversely, choosing any native parameters with all these gaps, the **original source coefficient signs** (including \(d>0\) and the permitted \(D,\beta,n\ge0\)), and the exact native stationary coefficient formulas gives a globally quantum-valid contact source. Its stationary linear system has one-dimensional kernel, so that ray is extreme in the canonical valid cone. This converse does **not** claim that arbitrary positive parameters satisfying only a subset of the gap inequalities yield nonnegative source coefficients.

The theorem does **not** prove separated tensor SOS or an SDP representation. It isolates exactly the physical sources to which the independent local theorems must apply.

## Lemma 1. The exact Laurent determinant and stationary-to-global physicality

Put \(D_u=1+u^2,D_v=1+v^2\). The full cleared expectation is
\[
\begin{split}
F={}&4\lambda D_u(W+vX)^2+4nD_uD_vX^2+4m(vW-X+uvY-uZ)^2\\
 &+4CD_v(uW-Y)^2+4AD_uD_vY^2
 +(4D+2d)(uW+uvX-Y-vZ)^2\\
 &+(4\alpha+2d)D_u(vY-Z)^2
 +(4\beta+2d)D_v(X+uZ)^2-2dD_uD_vZ^2.
\end{split}\tag{1}
\]
Let \(P(u,v)\) be its symmetric 4-state coefficient matrix. On \(Z=0\), the \(d\)-part contains
\[
2d[(uW+uvX-Y)^2+D_uv^2Y^2+D_vX^2],\tag{2}
\]
which is positive definite in \((W,X,Y)\) whenever \(d>0,uv\ne0\): vanishing gives \(X=Y=W=0\). Positive event terms only improve it. Therefore \(P\) has three strictly positive eigenvalues; its remaining eigenvalue is nonnegative exactly when its determinant is nonnegative. A zero has a one-dimensional kernel with \(Z\ne0\).

With \(U=u^2,V=v^2\), direct Cauchy--Binet expansion yields
\[
\det P=16D_u^2D_v^2u^2v^2
\left(h+\sum_{(i,j)\ne(1,1)}q_{ij}U^{i-1}V^{j-1}\right).
\tag{3}
\]
Define \(t_1=d+2(\lambda+n+D)\), \(t_2=d+2(\alpha+m+A)\), \(t_3=d+2(\lambda+C+D)\), \(t_4=d+2(\alpha+\beta+m)\), \(t_5=d+2(\lambda+\beta+n)\), \(t_6=d+2(\alpha+C+A)\), \(t_7=d+2(\beta+n+m)\), \(t_8=d+2(C+A+D)\). Every \(t_j>0\), and the eight outer coefficients are *exactly*
\[
\begin{array}{ll}
q_{00}=4\lambda\alpha t_7t_8,&q_{01}=4(\lambda Dt_7t_6+\alpha mt_5t_8),\\
q_{02}=4mDt_5t_6,&q_{10}=4(\alpha At_7t_3+\lambda nt_8t_4),\\
q_{12}=4(\beta mt_1t_6+CDt_5t_2),&q_{20}=4nAt_3t_4,\\
q_{21}=4(\beta nt_3t_2+CAt_1t_4),&q_{22}=4\beta Ct_1t_2.
\end{array}\tag{4}
\]
The central \(h=\det P(1,1)/256-\sum q_{ij}\) need not be positive. Formula (3)--(4) is checked as one generic symbolic identity by the accompanying script, without specializing any weight.

Suppose \(d>0,c_i\ge0\), and at a finite \(u_c,v_c>0\) we have the full stationarity equations
\[
P(u_c,v_c)\psi_c=0,
\quad \psi_c^TP_u(u_c,v_c)\psi_c
=\psi_c^TP_v(u_c,v_c)\psi_c=0.
\tag{5}
\]
At this point \(\operatorname{adj}P=\gamma\psi_c\psi_c^T\), \(\gamma>0\) up to the positive normalization of \(\psi_c\). Jacobi's determinant derivative therefore shows \(L=L_U=L_V=0\) at \(U_c=u_c^2,V_c=v_c^2\). Write \(E=\{-1,0,1\}^2\setminus\{0\}\) and \(w_e=q_eU_c^{e_1}V_c^{e_2}\ge0\). Then
\[
h=-\sum_e w_e,\qquad\sum_e w_e e=0,
\quad
L(U,V)=\sum_e w_e\left[(U/U_c)^{e_1}(V/V_c)^{e_2}-1\right].
\tag{6}
\]
The nonnegative balanced weight vector \(w\) decomposes into at most eight minimal positive balanced packets: four opposite pairs of weights \((1,1)\), four triangles with weights \((1,1,1)\), and four triangles with weights \((1,1,2)\), all up to lattice rotations/reflections. To see existence of this decomposition, choose a minimal subset whose convex hull contains zero, subtract the minimum multiple of its uniquely determined positive integer circuit weights, and repeat; each subtraction eliminates an active exponent. No limits occur.

For positive Laurent monomials \(p_i\), each packet has product \(\prod_i p_i^{b_i}=1\). The elementary identities
\[
 p+p^{-1}-2=(\sqrt p-p^{-1/2})^2,
\]
\[
 a^3+b^3+c^3-3abc=\tfrac12(a+b+c)[(a-b)^2+(b-c)^2+(c-a)^2],
\]
\[
 a^4+b^4+c^4+d^4-4abcd=(a^2-b^2)^2+(c^2-d^2)^2+2(ab-cd)^2
\tag{7}
\]
prove each circuit nonnegative (in the cubic take \(a=p_1^{1/3}\), etc.; in the quartic take \(a=p_1^{1/4},b=p_2^{1/4},c=d=p_3^{1/4}\)). Thus \(L(U,V)\ge0\) for every positive \(U,V\), over any real closed field. Equations (2)--(3) give \(P(u,v)\succeq0\) for every nonzero signed \(u,v\). At \(u=0,\infty\) or \(v=0,\infty\), one measurement pair becomes commuting and its CHSH core has operator norm at most 2, so \(H\succeq0\); all event terms remain positive. The projective homogenization of (1) includes these limits without dividing by zero. Two-reflection decomposition and real/imaginary state splitting extend validity to arbitrary quantum strategies.

This proves **stationarity (5) plus original coefficient nonnegativity implies global physicality**. Conversely, a physical source with a finite zero is PSD at all states and angles; the zero gives \(P\psi=0\) and the two angle derivatives vanish because the expectation has a minimum. Hence (5) is also necessary for a finite physical contact. *The determinant's separated tensor SOS alone does not imply an SOS for the full state matrix; that different question is deliberately not used here.*

## Lemma 2. Physical finite contact forces coercive Laurent support

Let \(E_+=\{e:q_e>0\}\). At a stationary zero, (6) gives \(\sum w_e e=0\). If zero is not in the **interior** of \(\operatorname{conv}E_+\), a nonzero separating functional has nonnegative scalar product with every active \(e\); balance forces all active exponents to lie on one line through the origin. The only possibilities are the four opposite pairs in \(E\), or no exponent at all.

For the horizontal pair \(q_{01},q_{21}>0\), the former arises from \(\lambda D\) or \(\alpha m\), and the latter from \(\beta n\) or \(CA\). Each of the four combinations forces one of \(\lambda n,CD,\beta m,\alpha A\) positive, activating an off-line coefficient \(q_{10}\) or \(q_{12}\). The vertical pair \(q_{10},q_{12}>0\) similarly forces one of \(\alpha m,CA,\beta n,\lambda D\) and activates \(q_{01}\) or \(q_{21}\). Thus neither axial pair is possible.

For the diagonal pairs, SG6 forces respectively
\[
(q_{00},q_{22}>0):\quad\lambda,\alpha,\beta,C>0,
\qquad
(q_{02},q_{20}>0):\quad n,m,A,D>0,
\]
with all off-line weights zero. These two positive-event matchings are related by one physical cycle rotation. In the second case, substituting \(u=v=1\), \(\psi=(1,-\varepsilon,\varepsilon,1)\) into the **whole** parent (1) gives
\[
F=16\varepsilon[(d+n+m+A+D)\varepsilon-d]<0
\]
when \(\varepsilon=d/[2(d+n+m+A+D)]\). This contradicts physicality. The first matching is therefore also excluded by relabeling.

If all eight outer coefficients vanish, the twelve positive-product pairs in (4) cannot simultaneously have both endpoints positive. Their eight maximal permitted sets of positive event weights are the physical rotations of \(\{\lambda,m,A\}\), with no more than three events; this follows directly by taking maximal independent sets of those twelve pair constraints (the checker verifies all 256 possible support subsets). The test \(u=v=1\), \(\psi=(1,-1,0,2)\) annihilates the positive event rows \(\lambda,m,A\) and has signed-core value \(-4d\). Physical rotation excludes each triple and its subsets. Therefore a physical finite noncommuting contact with \(d>0\) has **coercive** support. No limiting argument or denominator inversion is used.

Consequently its logarithmic angular Hessian
\[
\mathcal H_L=\sum_e w_e ee^T\succ0.\tag{8}
\]
With \(P=\left(\begin{smallmatrix}H_3&b\\b^T&c\end{smallmatrix}\right)\), \(H_3\succ0\), the completed square is
\[
F=(z+H_3^{-1}b)^TH_3(z+H_3^{-1}b)+s(u,v),
\qquad s=\frac{\det P}{\det H_3}.
\]
At contact, \(s=\nabla s=0\); its angular Hessian is a positive scalar multiple of a diagonal congruence of (8). Therefore the full Hessian of \(F\) in \((W,X,Y,u,v)\) at \(Z=1\) is positive definite. These calculations establish the invertible Jacobian needed for EN's analytic boundary opening.

## Lemma 3. Noncoercive support means one of twelve exclusive omission pairs

The event cycle is
\[
(\lambda,D,C,A,\alpha,m,\beta,n).\tag{9}
\]
Let a pair be **exclusive** when its two cycle indices have cyclic distance 3 or 4. There are eight pairs at distance 3 and four at distance 4.

The Laurent support is noncoercive exactly when there exists a nonzero direction \(w\) with \(w\cdot e\ge0\) for every active exponent \(e\). Because the eight exponents lie on the boundary of the \(3\times3\) lattice square, a feasible direction may be chosen among \((\pm1,0),(0,\pm1),(\pm1,\pm1)\): the closed cone of feasible directions has a boundary ray perpendicular to some nonzero active exponent. For any such direction, **every** packet in (4) with \(w\cdot e<0\) must be zero; the omitted events must hit every underlying positive product.

For \(w=(1,0)\), the forced zeros \(q_{00}=q_{01}=q_{02}=0\) say
\[
\lambda\alpha=\lambda D=\alpha m=mD=0.
\]
These four edges form a four-cycle, whose only minimal omission sets are \(\{\lambda,m\}\) and \(\{\alpha,D\}\). For the diagonal \(w=(1,1)\), the forced zeros are \(q_{00}=q_{01}=q_{10}=0\); its minimal omission sets are \(\{\lambda,\alpha\}\), \(\{\lambda,m,A\}\), and \(\{D,\alpha,n\}\). The two triples already contain axis-exclusive pairs. The remaining six normals give, analogously, the additional minimal pairs:

| \(w\) | Minimal exclusive omitted pairs |
|---|---|
| \((1,0)\) | \(\{\lambda,m\},\{\alpha,D\}\) |
| \((-1,0)\) | \(\{C,n\},\{A,\beta\}\) |
| \((0,1)\) | \(\{\lambda,A\},\{\alpha,n\}\) |
| \((0,-1)\) | \(\{D,\beta\},\{C,m\}\) |
| \((1,1)\) | \(\{\lambda,\alpha\}\) |
| \((-1,-1)\) | \(\{C,\beta\}\) |
| \((1,-1)\) | \(\{D,m\}\) |
| \((-1,1)\) | \(\{A,n\}\) |

Each diagonal also has two minimal three-omission sets, already containing one of the eight axis pairs; the full algebraic hitting-set rule follows from (4). Thus the minimal noncoercive omissions are **exactly** the twelve exclusive pairs. Every such pair forces a noncoercive exponent support by reversing the same product implications, proving the converse:
\[
\boxed{\ L\text{ coercive}\iff\text{every two omitted events are compatible}.\ }\tag{10}
\]
Physical compatibility means assignments on their common local observables agree; inspection of the eight event definitions gives compatible pairs at cyclic distance 1 or 2. Three pairwise compatible losing events extend to one deterministic assignment, so they number at most three, and the only compatible triples are the eight three-consecutive cycle vertices. In particular, every coercive source has at least five positive event weights.

The physical transformations \(g_A,g_B,g_T\) given explicitly in the attachment act on the cycle index \(i\bmod8\) by \(5-i,3-i,2-i\), respectively. Their compositions generate a rotation by one and a reflection, hence the whole finite dihedral group. They preserve the **entire** Bell parent (see the exact polynomial reflection identities in the checker), not merely its determinant. Consequently distance-3 omitted pairs form one orbit, represented by \(\{\lambda,m\}\), and distance-4 pairs form another, represented by \(\{\lambda,\alpha\}\). These are precisely the SG12/SG13 closed faces. Extra omitted events lie in their intersections and require no third closed-face type. All eight PR orientations are obtained separately by local output flips because any two PR correlation sign matrices have product \(-1\) and differ by row/column signs.

## Lemma 4. Extreme coercive sources have finite contact and five/six positive events

Normalize \(d+\sum c_i=1\). Suppose a valid \(d>0\) source has coercive support. Then \(L(U,V)\) tends to \(+\infty\) as \(\|(\log U,\log V)\|\to\infty\), because the active exponent polygon surrounds zero. It has a global minimum \(\mu\ge0\) at finite positive \(U,V\).

If \(\mu>0\), choose a nonzero direction among the active source coefficients preserving their total sum and all literal zero coefficients. Such a direction exists since (10) leaves at least five positive events plus \(d>0\). Small opposite perturbations retain positivity of the active source coordinates. Each active \(q_e\) varies by at most an arbitrarily small relative factor \(\varepsilon\) and the central \(h\) by absolute \(\varepsilon\). Since the omitted source events remain exactly zero, the absent \(q_e\)'s remain absent. Uniformly in \(U,V>0\),
\[
L'\ge(1-\varepsilon)L+\varepsilon h-\varepsilon
\ge(1-\varepsilon)\mu+\varepsilon h-\varepsilon>0
\]
for sufficiently small \(\varepsilon\). Both perturbed sources are globally physical by Lemma 1's three-positive-eigenvalue test and the direct projective endpoints, contradicting extremality. Thus \(\mu=0\) and a finite stationary zero exists.

At that zero, \(P\psi=0\) gives four equations and the two angle jets give two more, all homogeneous **linear equations in the original source coefficients**. If seven or eight events were positive, the active coefficient space would have dimension at least eight (counting \(d\)). The kernel of six linear equations has dimension at least two. There is a nonzero kernel direction preserving total normalization; small plus/minus shifts retain all coefficient signs and the same stationary zero. By Lemma 1 both shifted sources remain globally physical. This contradicts extremality. Hence at most six positive events. Lemma 3 supplies at least five.

The dihedral group makes any two compatible omissions equivalent to either the adjacent pair \(\{\lambda,n\}\) (C1) or the distance-two pair \(\{\lambda,\beta\}\) (C2). An allowed third omission in C1 must be \(D\) or \(\beta\); in C2 it must be \(n\). In both cases the three omissions are consecutive. This proves (T), **including all zero-event boundary supports**, without assuming their admission into any SOS space.

## Lemma 5. Both closed-face source pencils are exact, singular cases included

**Face \(\lambda=m=0\).** Define the symmetric quadratic matrix polynomial \(K(v)=K_0+vK_1+v^2K_2\) by
\[
K_0=\operatorname{diag}(4C+4D+2d,4A,4\alpha),\quad
K_2=\operatorname{diag}(4C,4A+4\alpha+2d,4D),
\]
\[
K_1=\begin{pmatrix}0&0&-(4D+2d)\\0&0&-(4\alpha+2d)\\-(4D+2d)&-(4\alpha+2d)&0\end{pmatrix}.
\]
Necessity of \(K(v)\succeq0\) follows from the physical parent by setting \(W=a/u\), letting \(u\to0\) *over the real field*, then putting \(R=a-Y\): the surviving parent is \(2t_5(1+v^2)X^2+(R,Y,Z)^TK(v)(R,Y,Z)\). Conversely, \(K_0+vK_1+v^2K_2\) has a positive definite leading two-state block for \(v\ne0\). Its determinant identity is
\[
\det K(v)=\frac{8(1+v^2)}{t_5}\,[q_{10}+h v^2+q_{12}v^4].
\]
Thus \(K(v)\succeq0\) forces \(h+q_{10}/V+q_{12}V\ge0\), and
\[
L=h+q_{10}/V+q_{12}V+U(q_{20}/V+q_{21}+q_{22}V)\ge0.
\]
Lemma 1 makes the source physical. The exact univariate quadratic matrix SOS lemma says \(K(v)\succeq0\) for all real \(v\) iff some skew \(S\) satisfies
\[
\boxed{\begin{pmatrix}K_0&K_1/2+S\\K_1/2-S&K_2\end{pmatrix}\succeq0.}\tag{11}
\]
The reason is finite-dimensional separation: the Gram image is closed because its PSD trace is \(\operatorname{tr}(K_0+K_2)\); its dual Hankel blocks decompose into finitely many moment atoms \((z,tz)(z,tz)^T\) and an endpoint-at-infinity remainder via the Schur complement \(W-VU^+V\succeq0\). This proves the quadratic matrix lemma and both directions of (11). Every statement is a fixed-dimensional first-order formula and transfers to real closed fields; no non-Archimedean limiting argument is required. At the stated strict point \(C=A=2,D=\alpha=d=n=\beta=1\), \(S=0\), the six leading principal minors are \((14,112,448,2576,28000,62500)\), so the homogeneous lift has a Slater point even though some represented sources may be singular.

**Face \(\lambda=\alpha=0\).** Put \(p=4C+4D+2d\), \(q=4m+4\beta+2d\), \(r=4D+2d\). The exact condition is \(\exists\xi\in\mathbb R\) with
\[
B_1=\begin{pmatrix}4m&-4m&\xi\\-4m&q+4n&q\\\xi&q&q\end{pmatrix}\succeq0,
\quad B_2=\begin{pmatrix}p+4A&4D&-p\\4D&4D&-(4m+r)-\xi\\-p&-(4m+r)-\xi&p\end{pmatrix}\succeq0.
\tag{12}
\]
For necessity, use the homogeneous physical limit \(u=\varepsilon s,v=\varepsilon,(W,Z)=(w,z)/\varepsilon\), leaving \(X=x,Y=y\). The limiting state form is exactly
\[
\begin{split}
G_s={}&4nx^2+4m(w-x-sz)^2+4C(sw-y)^2+4Ay^2\\
&+r(sw-y-z)^2+(4\beta+2d)(x+sz)^2-4dyz-2dz^2.
\end{split}
\]
It is PSD for all \(s\). Apply (11) to this quadratic matrix polynomial in \(s\). The zero coefficients of \(s^2x^2\) and \(s^2y^2\) force the corresponding Gram rows to vanish; coefficient matching on the six remaining features \((w,x,y,z,sw,sz)\) forces the cross-group entries to zero and leaves precisely one freedom \(\xi\), because \(w(sz)=z(sw)\). Reordering gives the two PSD blocks (12). Conversely the full algebraic identity
\[
G_s=(w,x,sz)B_1(w,x,sz)^T+(y,z,sw)B_2(y,z,sw)^T
\]
proves \(G_s\succeq0\). The limiting state matrix has three positive transverse eigenvalues and
\[
\det G_s=16(q_{02}+hs^2+q_{20}s^4),
\]
so \(h+q_{02}/t+q_{20}t\ge0\) for every \(t=s^2>0\). Adding the nonnegative Laurent remainder \(q_{12}V+q_{21}U+q_{22}UV\) proves \(L\ge0\); hence Lemma 1 gives physicality. At \(C=A=n=\beta=2,m=D=d=1,\xi=-5\), the respective leading minors are \((4,72,234)\) and \((22,72,234)\), giving strict feasibility. Both equivalences include all zero-coefficient and rank-deficient boundaries.

The full-parent source lifts (11)--(12) provide the already published *separated tensor* certificates by attained finite PSD polar pairing; that further argument is not used to classify the sources.

## Lemma 6. Every stationary C1/C2 source enters the stated strict native chart

**First the uniform real contact graph.** In the strict six-event C1 and C2 support faces, fix \(C=2\) and allow the other five positive event weights to vary through \((0,\infty)^5\). At \(d=0\) the Laurent central coefficient is strictly positive, containing a positive \(mCAD\) Cauchy--Binet term; its outer coefficients are nonnegative. Therefore some \(d>0\) is physical. A fixed CHSH-violating real test bounds \(d\) from above. At the maximal \(d_*(w)>0\), the coercive Laurent minimum must be zero (otherwise the previous relative-perturbation argument allows increasing \(d\)), and uniqueness follows from the positive angular Hessian (8). The state kernel has dimension one and may be normalized \(Z=1\).

Let \(z=(W,X,Y,u,v)\). At that contact, the Jacobian of \((\nabla_zF,F)=0\) in \((z,d)\) is
\[
\begin{pmatrix}\nabla_z^2F&\partial_d\nabla_zF\\0&F_H\end{pmatrix},
\qquad \nabla_z^2F\succ0,\quad F_H<0.
\tag{13}
\]
The core value \(F_H<0\) follows because some retained positive event amplitude must be nonzero: if all vanished, \(P_H\psi=0\), contradicting
\(\det P_H=-16D_u^2D_v^2u^2v^2\ne0\). The positive-state-core Schur form (checked exactly) is
\[
\frac{2D_uD_v}{h_0^2-u^2v^2}
\bigl\{u^2v^2(W^2+Z^2)-2uvh_0WZ\bigr\},\quad h_0=1+v^2+u^2v^2.
\tag{14}
\]
It is nonnegative for \(WZ\le0\); hence \(WZ>0\). After \(Z=1\), \(W>0\). Equation (13) is invertible. The real analytic implicit-function theorem provides a local contact branch, and Lemma 1 makes it physical. Its core value makes \(d\) maximal; uniqueness patches the branches into a single analytic graph over the connected five-dimensional positive source orthant.

**C1 signs.** Set \(\psi=(S,u(\theta-1),uY_0,1)\). For \(B=\operatorname{adj}P\), direct 3-by-3 cofactor expansion gives
\[
u B_{44}+B_{24}=16uv^2D_u^2D_v\,\mathfrak F>0,
\]
\[
\mathfrak F=(d+2D)\{(d+2\alpha)(mD_v+Cu^2v^2)+2AD_v(m+Cu^2)\}
+2mC\{d(u^2v^2+D_v)+2\alpha+2Dv^2D_u\}.
\tag{15}
\]
Since \(B=\gamma\psi\psi^T\) with \(\gamma>0\), this gives \(X+u>0\), i.e. \(\theta>0\). The walls \(Y_0=0\) and \(S-Y_0=0\) would delete the positive \(A\) or \(C\) event while preserving all of (5). Lemma 1 would leave a physical finite-contact source with an exclusive omitted pair, contrary to Lemma 2. Hence both walls are impossible.

Define \(R_m=vS-u(\theta-1)+u^2vY_0-u\), \(R_D=uS+u^2v(\theta-1)-uY_0-v\), and
\[
\mathcal D=Y_0D_u(uY_0+v)(1-uvY_0)
[S+uv(\theta-1)]R_m(-R_D).
\]
The **actual six-column stationary determinant**, for columns \((d,m,A,D,\alpha,\beta)\) after fixing \(C=2\), is exactly
\[
\boxed{\det M_{\rm C1}=-8192\theta u^2vD_u^2D_v^2\mathcal D.}\tag{16}
\]
Replacing the \(m\)-column by the stationary right-hand side \(-2\partial_C\) produces \(\det M_{m}\). Its **fully factored numerator** is
\[
\det M_m=16384\,Y_0\theta u^5v D_u^3D_v^3
(S-Y_0)(\theta-1)(uY_0+v)(uvY_0-1)R_D.
\]
Exact polynomial expansion, *without assuming invertibility*, gives the Cramer identity
\[
\det M_{m}\,[S+uv(\theta-1)]R_m
-\det M_{\rm C1}\,2u^3(S-Y_0)(1-\theta)D_v=0.
\]
Where (16) is nonzero, this is EN7:
\[
\boxed{m[S+uv(\theta-1)]R_m
=2u^3(S-Y_0)(1-\theta)D_v.}\tag{17}
\]
It is not legitimate to cancel (16) where it vanishes. Instead, the analytic contact graph above contains the exact seed
\[
(u,v;S,Y_0,\theta)=(1,1;1,1/4,1/2),\quad
(d,m,C,A,D,\alpha,\beta)=(27/2,4,2,9,5/4,13/4,9/4),
\]
whose kernel and jets vanish and for which (16) is nonzero. Thus the analytic determinant on the connected graph is not identically zero, so its nonzero set is dense; the **polynomial** equation (17), already proved there by the explicitly factored Cramer relation, extends by continuity to every positive support point. If \(\theta=1\), then \(m,S>0\) and (17) force \(R_m=0\); deleting the \(m\)-event gives the same forbidden finite-contact exclusive pair. Hence \(\theta\ne1\). Connectedness and the seed fix \(0<Y_0<S\) and \(0<\theta<1\) throughout.

The exact source/test separated C1 reflection is
\[
u^2F_{\sigma_1 c}(1/u,v;G_v\psi)=D_vF_c(u,v;\psi),
\quad G_v(W,X,Y,Z)=(vW-X,-W-vX,Z-vY,Y+vZ).
\tag{18}
\]
Its transformed last coordinate is \(uY_0+v>0\). Applying the just-proved three C1 signs at the transformed contact yields, respectively,
\(1-uvY_0>0\), \(S+uv(\theta-1)>0\), \(R_m>0\), and \(R_D<0\). Together these are all eight signs EN2, including \(X<0,Y>0\) and \(W+vX>0\). With these signs, every factor of \(\mathcal D\) in (16) is strictly positive, so the selected stationary matrix is invertible.

For the literal third omissions \(D=0\) or \(\beta=0\), open the missing coefficient to \(\varepsilon>0\) and use (13) at the **given finite contact** to solve for \((z,d)\). The branch remains physical by Lemma 1. The inequalities hold on the strict branch, and at \(\varepsilon=0\) all potential equality walls are excluded by (15), (14), deletion, (17), and its reflected equations. Thus all signs stay strict at both boundaries. No fictitious claim about a globally bounded optimizer or an infinitesimal source residue is used.

**C2 signs.** The five positive event-amplitude walls \(Y=0,uW-Y=0,R_m=0,R_D=0,vY-1=0\) delete, respectively, \(A,C,m,D,\alpha\); each would introduce an exclusive omitted pair next to \(\lambda=\beta=0\), so is excluded. At \(n>0\), \(X=0\) would delete \(n\), landing on the already proved C1 triple face and contradicting its sign \(X<0\). Consequently the signs of all six active event amplitudes are constant along the connected strict C2 graph. The exact stationary seed
\[
(d,n,m,C,A,D,\alpha)=(1,1/8,3/8,1/6,5/8,3/8,1/6),
\quad (u,v;\psi)=(1,1;1,-1/4,1/4,1)
\]
fixes the first six EN2 signs. At \(n=0\) use the C1 triple face directly.

For the omitted-row sign, write \(P_0\) for the C2 state matrix with \(n=0\). The term \(4nD_uD_vX^2\) alters only the XX entry. Therefore \(B_{24}\) is unchanged and \(B_{44}\) increases by \(4nD_uD_v\det(P_0|_{W,Y})\), so the exact cofactor equation is
\[
u B_{44}+B_{24}=16uv^2D_u^2D_v\mathfrak F
+4unD_uD_v\det(P_0|_{W,Y})>0.
\tag{19}
\]
The two-state minor is strictly positive by (2). Thus \(X+u>0\). The full C2 reflection
\[
u^2v^2F_{\sigma_2c}(1/v,1/u;Z,X,Y,W)=F_c(u,v;W,X,Y,Z)
\tag{20}
\]
exchanges \((m,D),(C,\alpha),(\lambda,\beta)\) and gives \(W+vX>0\). These establish EN2 also on \(n=0\), where the same cofactor formula remains strict.

**Native coordinates and uniqueness.** From all strict C1 signs define \(k=u_c\), \(S=W_c\), \(r=Y_c/(kS)\), \(\theta=1+X_c/k\), \(h=1-\theta\), \(\eta=v_c/(kS)\), \(\zeta=Sv_c/k\), \(e=k/S\), \(a=(1-r)/\eta\), \(t=\theta/\zeta\), and \(x=k^2\), after the allowed reflection making \(0<k\le1\). Then
\[
\begin{gathered}
0<r,\theta<1,\quad a,t,\eta,\zeta,e>0,\quad 0<a,t<2,\\
M=1+xr-t=R_m/(k\zeta)>0,\quad J=1+xh-a=-R_D/v_c>0,\\
E=1-x\eta h=(S+u_cv_c(\theta-1))/S>0,
\quad Q=1-x\zeta r=1-v_cY_c>0.
\end{gathered}\tag{21}
\]
The original native coefficients SD1 have the rational formulas (with \(C=2\))
\[
\begin{gathered}
d=2a\eta c_d/c_3,\quad m=2ae^2c_2/(Mc_3),\quad A=2a\eta^2c_4/c_3,\\
\alpha=2ae^{-2}k^6\eta^2c_\alpha/c_3,\quad
D+d/2=2ac_5/(Jc_3),\quad
\beta+d/2=2ae^2c_b/(tc_3).
\end{gathered}\tag{22}
\]
Here \(c_2,c_b,c_3,c_5,c_4,c_\alpha,c_d\) are explicit source polynomials in \((k^2,\eta,\zeta,a,t)\); their factored formulas and the **six** polynomial identities SD9--SD10 are reproduced and checked in the separate checker. The exact native contact SD11 satisfies the full kernel and angular jets by those identities. By (16) and (21), this solution is unique, so it equals the original physical source. Formula (22) applies unchanged when \(D=0\) or \(\beta=0\), provided the original source coefficient signs are imposed.

For C2 normalize \(d=1\), and set \(c=u_c,a=v_c/u_c,b=aW_c,r=-X_c/u_c,t=(u_cW_c-Y_c)/v_c,z=c^2\). Then
\[
\begin{split}
h&=b-a^2t=aY_c/u_c>0,\\
g&=(1+z)b-(1-r+a^2zt)=R_m/u_c>0,\\
j&=1+zr-t=-R_D/v_c>0,\\
&0<r<1,\quad 1-zh=1-v_cY_c>0,
\quad b-a^2zr=a(W_c+v_cX_c)>0.
\end{split}\tag{23}
\]
The six native coefficient formulas in the repository `c2_short.tex` give exactly the kernel and both angle jets; the checker verifies those six identities symbolically without assuming \(n>0\). The unweighted selected six-column stationary determinant is
\[
\boxed{-64ac^6rgt\,hj(1-zh)(1+z)^2(1+a^2z)^2\Delta,}\tag{24}
\]
where
\[
\Delta=1+(1+z)b+z(1+z)b^2+2z(1+z)br+z^2r^2+2zr
+a^2z\{(1+z)b+1+zr(2-r)\}>0.
\tag{25}
\]
Every displayed factor in (24) is nonzero on (23), including \(n=0\), and (25) is manifestly positive because \(a,b,z>0\) and \(0<r<1\). The determinant was recomputed directly from the six unweighted event-packet gradients. Consequently the six native stationary coefficients are unique for that contact.

Conversely, for native parameters satisfying (21) or (23), the displayed native coefficient formulas and the **required nonnegative original event weights**, the verified native contact equations give (5), hence global physicality by Lemma 1. If a nonzero such source decomposes into nonnegative canonical physical sources, every summand shares its literal zero event weights and must vanish with the same full kernel and two angle jets at the source's quantum contact. The selected nonzero determinant (16) or (24) makes the space of source coefficients satisfying these six equations one-dimensional (before fixing \(C=2\) or \(d=1\)). Thus every summand is proportional to the original source, proving extremality even with the allowed third omission. This is the promised reverse direction for actual native generators; **no SDP or tensor-SOS assumption is involved.**

## Appendix: native source-coefficient algebra used in the converse

These formulas make the contact-to-coefficient step reproducible without relying on a verbal appeal to "native stationarity." The script checks every displayed polynomial identity. All parameters below are **source** parameters and are unrelated to the independent testing angles.

**C1.** Use \(\kappa=k^2\), \(V_c=\kappa\eta\zeta\), \(L_c=1-a\eta\), \(H_c=1+(1-a)\eta\), \(Q_c=1-\kappa\zeta+aV_c\), \(R=1+\kappa L_c\), \(T=1+\kappa-\kappa\zeta t\). Define
\[
\begin{split}
P_4={}&1+\kappa+\eta\kappa(a-\kappa-1)-\zeta\kappa t
 +V_c[-\zeta\kappa t^2+2\kappa t+a(1-t)+t],\\
P_\alpha={}&1+\kappa-a\eta\kappa+\zeta[(1+\kappa)-(2+\kappa)t]
 +V_c[-\zeta t^2+a(t-1)+t],\\
P_d={}&(1+\kappa)^2-a\eta\kappa^2
 +\zeta\kappa[(\kappa+1)-(\kappa+2)t]
 +V_c[-\zeta\kappa t^2+\kappa a(t-1)+\kappa t+\kappa+1].
\end{split}
\]
Set
\[
\begin{array}{ll}
c_2=2(1+\kappa)(1-t\zeta)L_c(1+V_c)H_cQ_c,&
 c_b=2(1+\kappa)L_c[1+(1-t)\zeta]H_cQ_c,\\
c_3=2(1+\kappa)L_cH_cQ_c(1-\kappa\eta+tV_c),&
 c_5=2(1+\kappa)L_cH_cQ_c(1+V_c),\\
c_4=2Q_cP_4,&c_\alpha=2L_cH_cP_\alpha,\\
c_d=4L_cQ_cP_d.&
\end{array}
\]
The source weights are precisely (22). The six elementary polynomial balances (SD9--SD10), all identically zero, are
\[
\begin{aligned}
c_b-c_2-\zeta c_3&=0,&c_3-c_5+\kappa\eta c_b&=0,&c_2-(1-t\zeta)c_5&=0,\\
2\kappa c_2+2(1+V_c)L_cc_4-2\kappa^2\zeta Q_cc_\alpha-Q_cc_d&=0,\\
\kappa L_cc_2+c_5+(1+V_c)L_c^2c_4+\kappa Q_c^2c_\alpha
 +\tfrac12(\kappa\zeta L_c^2-2L_c-\eta)c_d&=0,\\
\kappa Rc_2+t\kappa V_cc_b+aV_cc_3+Tc_5+(1+\kappa)V_cL_c^2c_4
 -\kappa^2\zeta(1+\kappa)Q_cL_cc_\alpha
 -\tfrac12(1+\kappa)(Q_cL_c+\eta)c_d&=0.
\end{aligned}
\]
They give the complete kernel and all angular first jets at
\[
(u_c,v_c;W_c,X_c,Y_c,Z_c)=
\left(k,\frac{k^2\eta}{e};\frac{k}{e},-k(1-t\zeta),
\frac{k^2}{e}(1-a\eta),1\right),\qquad e=k\sqrt{\eta/\zeta}.
\]
No source-dependent testing equality was imposed. Original signs \(d,m,C,A,\alpha>0\), \(D,\beta\ge0\), \(D+\beta>0\) are part of the native converse hypothesis.

**C2.** With \(z=c^2\), \(A_c=1+zr\), \(h=b-a^2t\), \(g=(1+z)b-(1-r+a^2zt)\), \(j=1+zr-t\), define
\[
\begin{split}
\Delta={}&(1+z)zb^2+2z(1+z)br+(1+z)(1+a^2z)b
 +z^2(1-a^2)r^2+2z(1+a^2z)r+1+a^2z,\\
N_c={}&(1+z)(b+r)(b+a^2(1-t))-a^2(1+a^2z)(1-r)t.
\end{split}
\]
The full normalized source is \(d=1,\lambda=\beta=0\), with
\[
\begin{aligned}
n&=\frac{(b+r)A_c(A_c+a^2zt)-r\Delta}{2r\Delta},&
m&=\frac{(1+a^2z)(1-r)A_c(A_c+a^2zt)}{2g\Delta},\\
C&=\frac{A_cN_c}{2a^2t\Delta},&
A&=\frac{a^2A_c(b+r)(A_c+a^2zt)}{2h\Delta},\\
D&=\frac{(1+z)(b+r)(b-a^2zr)(A_c+a^2zt)}{2a^2j\Delta},&
\alpha&=\frac{z^2(b+r)N_c}{2(1-zh)\Delta}.
\end{aligned}
\]
The kernel and both angle jets vanish **identically** at
\[
(u_c,v_c;W_c,X_c,Y_c,Z_c)
=(c,ac;b/a,-cr,c(b/a-at),1).
\]
The separate short symbolic checker recomputes those six zero expressions. In addition to positive chart gaps, one must retain \(n\ge0\) and all original event nonnegativities. The positive identity
\[
N_c=(1+z)(b+r)h+a^2\{g+j+rt(1+a^2z)\}>0
\]
follows directly by the gap definitions. It establishes \(C,\alpha>0\) on the strict gap domain. The determinant \(\Delta\) here equals the manifestly positive expression (25) exactly.

## Consequence and exact scope

The classification is a complete geometric reduction: physical relabelings, events, the two exact closed-face source cones, and finite native C1/C2 extreme rays (including \(D=0,\beta=0,n=0\)) exhaust the canonical nonnegative source extreme rays. The SG7 stationary-to-global step and the EN positive-gap entrance are proved with exact algebra. The finite-order-field extreme decomposition was audited separately and is not used to manufacture the geometry. The further assertion that *each independent source/test pair has separated tensor SOS* is a different statement; it depends on the remaining RFN and thin-R10 local proof gates and the other printed local certificates. **No assertion of a completed global finite SDP lift or of an explicit SDP pencil follows from this geometric document alone.**
