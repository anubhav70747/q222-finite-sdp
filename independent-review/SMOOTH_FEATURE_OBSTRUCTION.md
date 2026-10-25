# Fixed smooth angular features cannot resolve the cubic quantum scale

The no-finite-NPA argument applies to more than measurement words. It
excludes any fixed finite angular feature space whose functions have
polynomial expansions through cubic order at the commuting origin.
The Gram coefficients may depend arbitrarily on the Bell source.

Let $\Omega=(0,\infty)^2$ and fix real functions $f_1,\ldots,f_N$ on
$\Omega$. Assume

$$
f_i(u,v)=p_i(u,v)+o(\|(u,v)\|^3),\qquad \deg p_i\leq3,
\tag{1}
$$

as the origin is approached within $\Omega$. A $C^3$ extension to a
neighborhood of the origin suffices; a one-sided polynomial Peano
expansion (1) is all the proof uses. Put
$Z_f=(f_i\psi_\ell)_{i,\ell}$, with all four state coordinates retained.
An admitted parent has an exact identity
$F_c=Z_f^TG_cZ_f$, $G_c\succeq0$, for every positive angle pair and every
real state vector.

**Theorem.** Write $\ell_T=4-2T$ for the local bound of the doubly tilted
operator $B_T$ used in the [native multiplier theorem](NATIVE_MULTIPLIER_OBSTRUCTION.md),
namely
$$B_T=A_0B_0-A_0B_1+A_1B_0+A_1B_1+(1-T)(A_1+B_0),$$
with the real qubit observables of [SG2](../latex/chapters/04_geometry.tex).
Put $D_u=1+u^2$, $D_v=1+v^2$, and set
$F_{T,w}=D_uD_v\psi^T(wI-B_T)\psi$.
For every finite $K>0$, there is $T_0>0$ such that no admitted parent
has

$$w=\ell_T+\kappa T^3,\qquad 0<T<T_0,\quad 0\leq\kappa\leq K.$$

If $w_f(T)$ is the infimum of admitted bounds above the quantum maximum
$b_T$, taking $w_f=+\infty$ when there are none, then

$$\frac{w_f(T)-b_T}{T^3}\longrightarrow+\infty.$$

## Bounded sources give bounded Grams

Remove linear dependencies among the $f_i$. There are finitely many
points $a_j\in\Omega$ for which
$M=\sum_j f(a_j)f(a_j)^T\succ0$: otherwise the evaluation vectors would
span a proper subspace, contradicting independence. Evaluate the Gram
identity at each $a_j$ and each state basis vector $e_\ell$. Then

$$
\operatorname{tr}G_c\leq
\frac{\sum_{j,\ell}F_c(a_j,e_\ell)}{\lambda_{\min}(M)}.
\tag{2}
$$

The numerator is linear in the source coefficients. Bounded sources
therefore have bounded positive semidefinite Grams. This also proves
closedness of the admitted source cone by taking a convergent Gram
subsequence. No regularity away from the origin is needed for this
finite evaluation argument.
The zero-dimensional feature space admits only the zero parent and
already satisfies the asserted exclusion.

## The cubic expansion forces polynomial square roots

Suppose contrary to the theorem that certificates exist along
$T_n=t_n^2\downarrow0$ with bounded $\kappa_n$. Pass to a subsequence
with $\kappa_n\to\kappa_0$. The exact polynomial test from the native
multiplier theorem is

$$
\psi_t(x,y)=\left(
1+\frac23t^4x^2y^2,
-\frac13t^3x^2y,
-tx+\frac23t^3xy^2,
-t^2xy-\frac23t^4x^3y
\right)^T.
$$

It gives, locally uniformly in $x,y$,

$$
t_n^{-6}F_{T_n,\ell_{T_n}+\kappa_nT_n^3}
(-t_nx,t_ny;\psi_{t_n})
\longrightarrow\frac43\left(M(x,y)+\frac{3\kappa_0}{4}-1\right),
\tag{3}
$$

where $M=1+x^4y^2+x^2y^4-3x^2y^2$. Choose an open rectangle
with $x<0<y$. Its tests $(-t_nx,t_ny)$ lie inside $\Omega$; no
extension of the feature functions to negative angles is required.

Let $B_n=G_n^{1/2}$. The source coefficients and hence $B_n$ are bounded
by (2). Equation (1) and the displayed state yield

$$
Z_f(-tx,ty;\psi_t)=\sum_{j=0}^3t^jZ_j(x,y)+o(t^3),
\qquad \deg Z_j\leq j.
$$

The normalized roots therefore equal

$$
t_n^{-3}B_n Z_f(-t_nx,t_ny;\psi_{t_n})
=P_n(x,y)+o(1),\qquad
P_n=\sum_{j=0}^3t_n^{j-3}B_nZ_j.
$$

Their squared norm is the left side of (3), so it is bounded on every
compact subrectangle. The remainder tends uniformly to zero. Values
at ten fixed unisolvent points thus bound every coefficient of the
vector polynomial $P_n$, whose degree is at most three. A subsequence
converges coefficientwise to $P$, giving

$$
\frac43\left(M+\frac{3\kappa_0}{4}-1\right)=\sum_jP_j^2.
\tag{4}
$$

The identity initially holds on the rectangle and hence holds everywhere
as a polynomial identity. Along each coordinate axis its right side is
constant. Each cubic root must consequently have the form
$r_j+s_jxy+h_jx^2y+k_jxy^2$. The coefficient of $x^2y^2$ in their square
sum is $\sum_js_j^2\geq0$, whereas (4) requires $-4$. This is a
contradiction. If the axis constant were negative, positivity would
already give the contradiction.

Finally $b_T-\ell_T=\Theta(T^3)$, with
$b_T<\ell_T+4T^3$ for small $T$. Excluding every bounded $\kappa$
therefore proves the stated divergence, whether or not the infimum
defining $w_f$ is attained.

## Consequence for native contacts and a general finite lift

No common feature library satisfying (1) admits all strict native C1/C2
sources and their physical relabelings in the original testing frame.
Indeed, augment the library by the event features $1,u,v,uv$. The
[contact-preserving decomposition](NATIVE_MULTIPLIER_OBSTRUCTION.md)
writes the tight doubly tilted source as a finite sum of native-orbit
sources and events. Third-omission sources are limits of strict contacts,
and closedness from (2) passes their certificates to those limits.
Adding the Grams would admit the tight parent, contradicting the theorem.

The common-frame hypothesis matters. A feature smooth in one native
chart need not remain smooth after a reciprocal physical relabeling.
A finite collection of libraries is excluded when all their pulled-back
features satisfy (1) at the same commuting origin: concatenate the lists.

Evaluated finite NPA words have analytic angular coefficients there and
fall within this obstruction. An arbitrary finite SDP lift need not use
this angular-state factorization. Fractional functions without cubic
polynomial expansions are also outside the theorem. The [scalar
half-power identity](SIXTH_ROOT_FEATURES.md) demonstrates that distinction,
without proving admission of the full parent. Neither $C^{2,1}$ nor
$C^2$ libraries are settled here, and no sharp regularity threshold is
claimed.
