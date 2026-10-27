# Exact symmetric lift, including homogeneous coefficients

The supplied symmetric lift is correct after making its domain inequalities explicit and supplying the endpoint and singular-case arguments below. A stronger statement holds: the same two Gram blocks represent the full homogeneous support cone for the planar behavior projection with coordinates

\[
M=\langle A_0\rangle+\langle B_0\rangle,\qquad
S=\langle A_0B_0+A_0B_1+A_1B_0-A_1B_1\rangle.
\]

This is an independently reconstructed proof from the physical Bell operator. The exact identity program `check_exact.py` verifies the named polynomial identities used here; its output alone is not the universal proof.

## The theorem

For arbitrary real \(h,a,\lambda\), define

\[
\mathcal C=\{(h,a,\lambda):\lambda\ge
\max_{p\in\mathcal Q}[hS(p)+aM(p)]\}.
\]

Define the real symmetric quartic matrix polynomial

\[
K_{h,a,\lambda}(t)=
\begin{pmatrix}
\lambda-2h+2\lambda t^2+(\lambda+2h)t^4&-2a(1+t^2)^2&0\\
-2a(1+t^2)^2&\lambda-2h+(2\lambda-8h)t^2+(\lambda+2h)t^4&-8ht^3\\
0&-8ht^3&\lambda+2h+(2\lambda+8h)t^2+(\lambda-2h)t^4
\end{pmatrix}.
\]

Then the following are equivalent:

1. \((h,a,\lambda)\in\mathcal C\).
2. \(K_{h,a,\lambda}(t)\succeq0\) for every real \(t\).
3. There are \(U\in\mathbb S^5_+\), \(V\in\mathbb S^4_+\) satisfying the sixteen affine equations displayed below.

In particular, this is one homogeneous affine lift for all three visible variables. Setting \(h=1\) gives the full common-tilt epigraph for every real \(a\). Its restriction to \(0\le a\le1\) is exactly the supplied \(E_{\rm sym}\), and

\[
\operatorname{sxdeg}(E_{\rm sym})=2.
\]

The last equality is an existence theorem for a representation using blocks of order at most two; the explicit Gram lift here has orders five and four.

## 1. Physical block and the local bound

Jordan reduction permits

\[
A_0=B_0=Z,\quad A_1=cZ+rX,\quad B_1=dZ+sX,
\quad c,d\in[-1,1],\quad r=\sqrt{1-c^2},\ s=\sqrt{1-d^2}.
\]

Changing the sign of a transverse component is a local conjugation, so nonnegative \(r,s\) lose no spectra. One-dimensional blocks embed as endpoint observables and suitable states. The Bell matrix is real symmetric, so a real eigenvector attains its largest eigenvalue. Thus support optimization is over a compact set of such two-qubit blocks and states; direct sums cannot increase it.

With \(\Omega=1+c+d-cd\), direct expansion of
\(H=h[A_0(B_0+B_1)+A_1(B_0-B_1)]+a(A_0+B_0)\) gives

\[
H=\begin{pmatrix}
h\Omega+2a&h(1-c)s&h(1-d)r&-hrs\\
h(1-c)s&-h\Omega&-hrs&h(d-1)r\\
h(1-d)r&-hrs&-h\Omega&h(c-1)s\\
-hrs&h(d-1)r&h(c-1)s&h\Omega-2a
\end{pmatrix}.
\]

Write \(H_*=|h|\), \(b=|a|\), and \(L=2H_*+2b\). Deterministic assignments attain \(L\): take \(A_0=B_0=\operatorname{sign}(a)\); for \(h\ge0\) set \(A_1=B_1=\operatorname{sign}(a)\), and for \(h<0\) take their negatives. If \(a=0\), either sign works. Thus a universal upper bound must satisfy \(\lambda\ge L\).

## 2. Determinant positivity and global diagonal reduction

For \(N=\lambda I-H\), write
\(\det(uI+N)=u^4+e_1u^3+e_2u^2+e_3u+e_4\). Exact identities give

\[
e_1=4\lambda,\quad e_2=6\lambda^2-4a^2-8h^2,
\quad e_3=4\lambda^3-8\lambda a^2-16\lambda h^2-8ha^2\Omega.
\]

Here \(-2\le\Omega\le2\). For \(\lambda=L+v\), \(v\ge0\),

\[
e_2=16H_*^2+48H_*b+20b^2+24(H_*+b)v+6v^2\ge0,
\]

and

\[
\begin{aligned}
e_3&\ge4\lambda^3-8\lambda b^2-16\lambda H_*^2-16H_*b^2\\
&=16b(b+2H_*)^2+(32H_*^2+96H_*b+40b^2)v
 +24(H_*+b)v^2+4v^3\ge0.
\end{aligned}
\]

Thus \(N\succeq0\) is equivalent to \(e_4=\det N\ge0\). Indeed, if all coefficients are nonnegative, the displayed polynomial is positive at every \(u>0\), whereas a negative eigenvalue of \(N\) would give a positive root. This also covers zero coefficients, zero determinant and multiple zero eigenvalues.

Set \(x=(1-c)/2\), \(y=(1-d)/2\), so \(x,y\in[0,1]\). The determinant is

\[
\begin{aligned}
D(x,y)&=D_0+64h^2(4h^2-a^2)xy(x+y)\\
&\quad +(32\lambda ha^2+128h^2a^2-256h^4)xy-256h^4x^2y^2,\\
D_0&=[(\lambda-2h)^2-4a^2](\lambda+2h)^2.
\end{aligned}
\]

For \(z=\sqrt{xy}\in[0,1]\),

\[
D(x,y)-D(z,z)=64h^2(4h^2-a^2)xy(\sqrt x-\sqrt y)^2.
\]

When \(|a|\le |h|\), this is nonnegative. Consequently positivity of every diagonal block \(c=d\), together with \(\lambda\ge L\), implies positivity of every physical block. The converse is immediate. No division by \(x,y,h\), or a determinant is involved, so the edges and all singular blocks remain included. The determinant argument in fact works throughout \(|a|\le2|h|\), but the smaller range suffices.

## 3. Critical and supercritical coefficients

The missing range is handled directly at the operator level. Write
\(\widehat S=A_0(B_0+B_1)+A_1(B_0-B_1)\), and put
\(C_B=(B_0+B_1)/2\), \(D_B=(B_0-B_1)/2\). Dichotomic observables satisfy \(C_B^2+D_B^2=I\), \(C_BD_B+D_BC_B=0\), and cross-party commutation. Hence

\[
4I-[\widehat S+2A_0]
=(I-C_B-A_1D_B)^2+
\frac14\sum_{y=0}^1[(I-A_0)(I+B_y)]^2\succeq0.
\]

Party exchange gives \(\widehat S+2B_0\le4I\). Local output flips give
\(\sigma\widehat S+2\tau A_0\le4I\) and
\(\sigma\widehat S+2\tau B_0\le4I\) for either \(\sigma,\tau\in\{\pm1\}\). Their average gives

\[
\sigma\widehat S+\tau(A_0+B_0)\le4I.
\]

For \(b\ge H_*\), multiply this inequality by \(H_*\), with \(\sigma=\operatorname{sign}(h)\), \(\tau=\operatorname{sign}(a)\), and add
\((b-H_*)[2I-\tau(A_0+B_0)]\succeq0\). The result is

\[
h\widehat S+a(A_0+B_0)\le (2H_*+2b)I.
\]

The deterministic assignments above attain it. If \(h=0\), the same conclusion follows immediately from \(-2I\le A_0+B_0\le2I\). Thus the exact support value is \(L\) on the whole critical and supercritical range, including both signs of both coefficients.

## 4. Three-dimensional block, endpoints and quartic polynomial

On \(c=d\), use the ordered orthonormal basis

\[
\phi_+=(|00\rangle+|11\rangle)/\sqrt2,\quad
\phi_-=(|00\rangle-|11\rangle)/\sqrt2,\quad
\psi_+=(|01\rangle+|10\rangle)/\sqrt2,\quad
\psi_-=(|01\rangle-|10\rangle)/\sqrt2.
\]

The first three basis vectors give

\[
H_3(c)=\begin{pmatrix}
2hc&2a&0\\
2a&2h(1+c-c^2)&2h(1-c)r\\
0&2h(1-c)r&-2h(1+c-c^2)
\end{pmatrix},
\]

and \(\psi_-\) has eigenvalue \(-2hc\), bounded above by \(2|h|\le L\).

With \(c=(1-t^2)/(1+t^2)\), \(r=2t/(1+t^2)\), exact substitution gives

\[
K_{h,a,\lambda}(t)=(1+t^2)^2[\lambda I_3-H_3(c(t))].
\]

Finite \(t\ge0\) covers \(-1<c\le1\); \(c=-1\) follows by continuity. Negative \(t\) gives the same PSD test because
\(K(-t)=R K(t)R\), \(R=\operatorname{diag}(1,1,-1)\). Thus no measurement endpoint is discarded and no extra physical condition is introduced.

Conversely, positivity of \(K\) forces \(\lambda\ge L\) without appending a nonlinear constraint. The upper-left two-dimensional block of \(K(0)\) forces \(\lambda\ge2h+2|a|\); the corresponding block of
\(\lim_{|t|\to\infty}t^{-4}K(t)\succeq0\) forces \(\lambda\ge-2h+2|a|\). Their combination is exactly \(\lambda\ge2|h|+2|a|\). The scalar sector is therefore controlled. The low-tilt determinant argument and the supercritical operator inequality now prove the equivalence of statements 1 and 2 for all real \(h,a,\lambda\).

At \(h=0\), explicitly,

\[
K=(1+t^2)^2\begin{pmatrix}\lambda&-2a&0\\-2a&\lambda&0\\0&0&\lambda\end{pmatrix},
\]

so the equivalence includes the completely degenerate origin and the rays \(\lambda=2|a|\).

## 5. Matrix SOS and all sixteen equations

Hanselka and Sinn, Corollary 2.8, applies to every PSD real symmetric univariate polynomial matrix, including an identically singular one. It gives \(K=Q^TQ\) with polynomial \(Q\). Each diagonal entry \(K_{jj}\) is a sum of squares of the entries in column \(j\) of \(Q\), so \(\deg K_{jj}\le4\) forces every factor entry to have degree at most two: the highest-degree sum of squares cannot cancel. This argument also covers vanishing diagonal entries, lower degree and rank-deficient leading coefficients.

For \(\xi=(\xi_1,\xi_2,\xi_3)^T\), every factor is therefore a linear combination of the nine monomials \(t^k\xi_j\), \(0\le k\le2\). Average its Gram matrix under \((t,\xi)\mapsto(-t,R\xi)\). PSD is preserved; opposite parity blocks vanish. The even and odd monomial vectors are

\[
Z_+=(\xi_1,\xi_2,t\xi_3,t^2\xi_1,t^2\xi_2)^T,
\quad Z_-=(\xi_3,t\xi_1,t\xi_2,t^2\xi_3)^T.
\]

Consequently \(\xi^TK(t)\xi=Z_+^TUZ_++Z_-^TVZ_-\), with \(U\succeq0,V\succeq0\). Coefficient matching is precisely

\[
\begin{array}{lll}
u_{11}=\lambda-2h,&u_{12}=-2a,&u_{22}=\lambda-2h,\\
v_{11}=\lambda+2h,&u_{13}+v_{12}=0,&u_{23}+v_{13}=0,\\
2u_{14}+v_{22}=2\lambda,&u_{15}+u_{24}+v_{23}=-4a,&
2u_{25}+v_{33}=2\lambda-8h,\\
u_{33}+2v_{14}=2\lambda+8h,&u_{34}+v_{24}=0,&u_{35}+v_{34}=-8h,\\
u_{44}=\lambda+2h,&u_{45}=-2a,&u_{55}=\lambda+2h,\\
v_{44}=\lambda-2h.&&
\end{array}
\]

There are exactly sixteen scalar monomials compatible with the parity: three degrees for each diagonal quadratic form, three for \(\xi_1\xi_2\), and two for each of \(\xi_1\xi_3\), \(\xi_2\xi_3\). Conversely these equations make the Gram identity exact, so PSD of \(U,V\) implies PSD of \(K(t)\) for every real \(t\). This proves the theorem.

To describe the original \(E_{\rm sym}\), impose \(h=1\) and the scalar affine inequalities \(a\ge0\), \(1-a\ge0\). They must not be dropped: the sixteen Gram equations by themselves describe the all-real-tilt epigraph. These scalar inequalities do not increase the largest PSD block order.

## 6. Extension degree two without importing a boundary expansion

Let \(f(a)=\max_{p\in\mathcal Q}[S(p)+aM(p)]\). Compactness makes \(f\) continuous; it is convex as a supremum of affine functions. The fixed Gram lift makes its epigraph semialgebraic by real algebraic elimination. Hence \(E_{\rm sym}\) is a closed convex semialgebraic subset of \(\mathbb R^2\).

Scheiderer's Theorem 0.1 gives \(\operatorname{sxdeg}(E_{\rm sym})\le2\). Its hypotheses do not require compactness or smoothness of the set.

For completeness, nonpolyhedrality can be proved without the supplied near-critical cubic expansion. Put \(w=2\sqrt2\). For the untilted operator,

\[
\widehat S(c,d)^2=4I+4rs\,Y\otimes Y,
\]

so reaching \(w\) requires \(r=s=1\), hence \(c=d=0\). Its top eigenspace is one-dimensional, spanned by

\[
v=\cos(\pi/8)\phi_-+\sin(\pi/8)\psi_+.
\]

This state has both reduced states \(I/2\), so every CHSH-maximizing Jordan block has \(M=0\). If \(p_a\) maximizes \(f(a)\), compactness and the bound \(|M|\le2\) imply that all limits of \(p_a\) as \(a\downarrow0\) maximize CHSH. Thus \(M(p_a)\to0\), and

\[
0\le\frac{f(a)-w}{a}\le M(p_a)\longrightarrow0.
\]

But \(f(a)>w\) for every \(a>0\). Indeed, keep \(c=d=0\), write \(q=\cos(\pi/8)>0\), and test the normalized vector \((v+\varepsilon\phi_+)/\sqrt{1+\varepsilon^2}\). Its Bell value is

\[
\frac{w+4aq\varepsilon}{1+\varepsilon^2}>w
\quad\text{at}\quad\varepsilon=\frac{2aq}{w}.
\]

If the restricted epigraph were a polyhedron, its lower boundary would be affine on some interval \([0,\delta]\). The zero right derivative would make that segment constant, contradicting the strict inequality. Thus \(E_{\rm sym}\) is nonpolyhedral. A lift with only order-one PSD blocks is a finite linear inequality system, whose projection is a polyhedron. Therefore \(\operatorname{sxdeg}(E_{\rm sym})\ge2\), proving equality. The same argument gives extension degree two for the all-real common-tilt epigraph.

## 7. Optional exact quartic boundary identity

For \(h=1\), the diagonal determinant is

\[
D(z,z)=D_0+32z^2[(\lambda+4)a^2-8]+128(4-a^2)z^3-256z^4.
\]

The exact stationary resultant is

\[
\operatorname{Res}_z\bigl(D,\partial_zD/z\bigr)
=2^{36}(\lambda-2)^2(2\lambda+a^2)^2F_4(\lambda,a),
\]

where

\[
\begin{aligned}
F_4={}&4\lambda^4+(16-4a^2)\lambda^3
 +(11a^4-48a^2-16)\lambda^2\\
&+(8a^6-4a^4-80a^2-128)\lambda
 +20a^6-84a^4+64a^2-128.
\end{aligned}
\]

This identity has been checked exactly. It is not by itself a boundary-root selection rule: endpoint minima and extraneous resultant factors must still be handled. None of the lift or extension-degree proof relies on choosing a quartic root.

## 8. Exact dual membership of the full planar behavior projection

The parent's independently constructed adjoint provides the following completion. Let \(A(U,V)\) be the vector of the sixteen left-hand sides above, and let \(Bx\) be the right-hand side for \(x=(\lambda,a,h)\). Then

\[
\widetilde{\mathcal C}=\{x:\exists U,V\succeq0,\ A(U,V)=Bx\}
\]

is the support cone proved above, with its coordinates reordered. This conic system is strictly feasible: at \(x_0=(1,0,0)\),

\[
U_0=\operatorname{diag}(1,1,2,1,1)\succ0,\qquad
V_0=\operatorname{diag}(1,2,2,1)\succ0,
\quad A(U_0,V_0)=Bx_0.
\]

Conic Slater duality (equivalently, the intersection-dual rule with this strict feasible point) therefore gives the attained representation

\[
\widetilde{\mathcal C}^{\,*}
=\{B^Tz:A^*z\in\mathbb S^5_+\times\mathbb S^4_+\}.
\]

The image on the right is also directly closed. Writing \(A^*z=(\widehat U(z),\widehat V(z))\), its first coordinate is

\[
\ell=z_1+z_3+z_4+2z_7+2z_9+2z_{10}+z_{13}+z_{15}+z_{16}
=\operatorname{Tr}(\widehat U U_0)+\operatorname{Tr}(\widehat V V_0).
\]

A convergent image sequence has bounded \(\ell\); positive definiteness of \(U_0,V_0\) bounds the traces and hence all entries of the PSD matrices \(\widehat U,\widehat V\). The adjoint is injective (rank sixteen), so \(z\) is bounded as well. A convergent subsequence supplies the required limit preimage. The strict-feasibility argument above is the license for conic duality; closedness alone would not justify interchanging a preimage with closure.

Let \(W=\{(M(p),S(p)):p\in\mathcal Q\}\), which is compact and convex by Jordan reduction and convex mixing. The cone generated by \(\{(1,-M,-S):(M,S)\in W\}\) is closed, and its dual is exactly \(\widetilde{\mathcal C}\). The bipolar theorem now gives

\[
(M,S)\in W\quad\Longleftrightarrow\quad
\exists z:\widehat U(z)\succeq0,\ \widehat V(z)\succeq0,
\quad B^Tz=(1,-M,-S).
\]

Thus the homogeneous extension genuinely yields one affine lift of the full planar behavior projection, with the same orders five and four. It is not restricted to a subinterval of support directions. Its semialgebraicity follows from this lift. It is nonpolyhedral because a compact polyhedral \(W\) would make \(f(a)=\max_W(S+aM)\) a finite maximum of affine functions, contradicting Section 6. Scheiderer's planar theorem therefore also gives
\(\operatorname{sxdeg}(W)=2\).

The explicit matrices and the three affine equations are in the parent's `extension/verify_membership_and_asymmetric.py`; its coefficient/adjoint checks have been inspected, without duplicating that program here.

## Primary sources and bibliography

- Claus Scheiderer, *Second-Order Cone Representation for Convex Sets in the Plane*, SIAM Journal on Applied Algebra and Geometry **5**(1), 114–139 (2021), DOI [10.1137/20M133717X](https://epubs.siam.org/doi/10.1137/20M133717X). [Theorem 0.1 in the public author manuscript](https://arxiv.org/html/2004.04196v2) applies exactly to closed convex semialgebraic subsets of the plane. The supplied arXiv reference identifies the right preprint but should preferably be updated to this publication.
- Christoph Hanselka and Rainer Sinn, *Positive semidefinite univariate matrix polynomials*, Mathematische Zeitschrift **292**, 83–101 (2019), DOI 10.1007/s00209-018-2137-7. [Corollary 2.8, p. 9 of the public author manuscript](https://arxiv.org/pdf/1707.08261) explicitly includes singular polynomial matrices. The supplied bibliography entry is correct. The degree bound above is an immediate independent consequence of sums of real squares, so no genericity or nonzero-determinant hypothesis is being imported.

Scope: the original symmetric proof, its homogeneous extension, and the resulting exact full planar projection lift are established here. The NPA nontermination theorem, the general asymmetric reduction, and explicit order-two Gram blocks are not established by this report.
