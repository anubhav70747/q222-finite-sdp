# A partially entangled extreme inequality outside 1+AB

The additional certificates must cover changes in the Schmidt coefficients,
as well as unequal measurement angles. The following rational example makes
this requirement explicit. It is globally quantum valid, generates an extreme
ray of the full Bell polar, and has no `1+AB` certificate. Both its validity
and its exclusion have short exact proofs.

## The inequality and its contact

Use the [canonical slack and event order](../extracted-lifts/c2-equal-angle-cone/)
and write

$$f=H+4\sum_i a_i p_i,\qquad
H=1-\frac{E_{00}-E_{01}+E_{10}+E_{11}}2.$$

Take $\lambda=\beta=0$ and

$$
(d,n,m,C,A,D,\alpha)
=\left(1,\frac7{66},\frac5{33},\frac5{22},
              \frac{10}{33},\frac{16}{33},\frac9{22}\right).\tag{1}
$$

For $A_0=B_0=Z$ and $A_1=B_1=X$, the Bell operator annihilates

$$\psi=(5/4,-1/4,1/2,1)^T.\tag{2}$$

Its normalized concurrence is $22/23$, so the contact is not maximally
entangled. In the native C2 coordinates, (1) corresponds to
$(a,b,c,t,r)=(1,5/4,1,3/4,1/4)$.

Global validity can be proved directly. Parametrize the second measurements by

$$A_1=\frac{(1-u^2)Z+2uX}{1+u^2},\qquad
B_1=\frac{(1-v^2)Z+2vX}{1+v^2},$$

with the tensor factors understood, and let $B_f(u,v)$ be the Bell operator
of (1). Its determinant is

$$
\det\!\left[66(1+u^2)(1+v^2)B_f(u,v)\right]
=2048000(1+u^2)^2(1+v^2)^2 P(u,v),\tag{3}
$$

where

$$
\begin{aligned}
P={}&189u^4v^2+98u^4+216u^2v^4+270u^2\\
   &+152v^4+135v^2-1060u^2v^2.
\end{aligned}
$$

The six positive coefficients sum to $1060$. Their weighted mean exponents
are $(2,2)$. Weighted AM--GM therefore gives $P\ge0$ for all real $u,v$;
zeros of either variable follow by continuity.

To pass from the determinant to positivity, the CHSH core has three strictly
positive eigenvalues whenever the two pairs of measurements are
noncommuting. Indeed its eigenvalues are
$1\pm\sqrt{1+\kappa}$ and $1\pm\sqrt{1-\kappa}$, with
$\kappa=4|uv|/[(1+u^2)(1+v^2)]\in(0,1]$.
Adding the nonnegative event operators preserves those three
positive eigenvalues. Equation (3) forces the remaining eigenvalue to be
nonnegative. When a pair commutes, the core itself is positive semidefinite.
Continuity also includes infinite half-angle coordinates. The two-observable
block decomposition on each party then establishes validity in arbitrary
dimension.

## Six forbidden words exclude an NPA certificate

The two deterministic zeros of (1) compress any nine-word `1+AB` Gram to
the seven roots

$$
R=\bigl(A_1-I,B_1-I,A_0B_0-I,A_1B_1-I,
        B_0-A_0,A_0B_1-A_0,A_1B_0-A_0\bigr)^T.
$$

The [seven-root reduction](EXACT_LIFT_ROUTE.md) imposes the full universal
word identity $f=R^\dagger Q R$, $Q\succeq0$. At the contact (2), let $T$
have columns $R_i\psi$. Its first four columns have determinant $-1/2$,
and its kernel has basis

$$
N=\begin{pmatrix}
-1&3/4&9/8\\
1&-5/4&-7/8\\
3&0&3\\
-3&-3/4&-9/8\\
1&0&0\\0&1&0\\0&0&1
\end{pmatrix}.\tag{4}
$$

Since $f$ vanishes there and $Q\succeq0$, $QT^T=0$. Thus
$Q=NSN^T$ for a symmetric $3$-by-$3$ matrix $S\succeq0$.
In the order $(S_{11},S_{12},S_{13},S_{22},S_{23},S_{33})$, the coefficients
of the six forbidden words

$$B_0B_1,\ A_0A_1,\ A_0A_1B_0B_1,\ A_0A_1B_0,
\ A_0A_1B_1B_0,\ A_0A_1B_1$$

are given by

$$
M=\begin{pmatrix}
1&7/4&-31/8&0&9/4&-9/8\\
1&-11/4&23/8&-3/2&-3&15/8\\
-9&-9/4&-99/8&0&-9/4&-27/8\\
-3&9/4&-5/8&0&5/4&19/8\\
0&0&0&0&1&0\\
3&11/4&33/8&3/2&3&9/8
\end{pmatrix},\qquad \det M=-405.\tag{5}
$$

All these coefficients must vanish. Hence $S=0$ and $Q=0$, contradicting
the nonzero target (1). The products in (5) use only reflection identities
and commutation between the parties; no same-party commutation is imposed.
A complex Hermitian Gram gives the same contradiction by taking its real
part in this real word basis.

## Extremality and the required additional coverage

The six derivatives of the six active event packets, with respect to
$(W,X,Y,Z,u,v)$ at (2), form a matrix of determinant $-297$.
Every canonical valid summand in a decomposition of (1) has
$\lambda=\beta=0$, because these nonnegative coefficients sum to zero.
Each summand is tight at the contact, so its six derivatives vanish there.
Once $d$ is fixed, these equations uniquely determine its six active event
coefficients.

This proves extremality in the nonnegative canonical cone. It also proves
extremality in the full Bell polar. Each canonical event coefficient of a
valid summand is nonnegative, since it is one quarter of the slack's value
at the corresponding deterministic facet vertex. At (2),

$$H=-17/46<0,\qquad p_A=2/23.$$

A tight valid summand cannot have a negative core coefficient: its core
term would be strictly positive and all its event terms nonnegative.
All summands therefore belong to the canonical cone, where the preceding
rank argument forces proportionality.

The equal-angle maximally entangled cone obeys $d=4A-2m-2D$.
For (1), the right side is $-2/33$. Thus this extreme ray belongs neither
to that cone nor to `1+AB`, and cannot be supplied by their sum.
This statement concerns the specified equal-angle cone; membership in all
its physical relabelings is a separate question.

The NPA obstruction also holds generically in the strict five-parameter
native C2 family. The contact, the matrix $T$, and
$N=\begin{pmatrix}-T_0^{-1}T_{\rm last}\\I_3\end{pmatrix}$ depend rationally on $(a,b,c,t,r)$ wherever
$T_0$ is invertible. The same six-word matrix then depends rationally on
these parameters. Its determinant is nonzero at (1), so the two required
minors remain nonzero on a nonempty open neighborhood and on a relative
Zariski-open dense subset of the full-dimensional strict interior. This
includes unequal-angle, nonmaximally entangled contacts. Exceptional loci,
including lower-dimensional faces, require their own calculation.

The [exact checker](checks/c2_nonmaximal_npa_obstruction.py) reconstructs
(1)--(5), the determinant (3), the AM--GM balances and the stationary minor.
No numerical optimization is used. The example identifies an additional
family requiring certificates; it does not provide their uniform admission.
