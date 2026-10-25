# Two C1 extreme inequalities outside 1+AB

The fixed nine-word NPA cone also misses native C1 sources. A single
integer positive-definite matrix separates the two examples below,
including one on the boundary $\beta=0$. Their quantum validity follows
from elementary AM--GM identities, and their extremality follows from
eight linear contact equations.

## Sources and global validity

Use $f=dH+4\sum_i a_i p_i$ in the
[canonical event convention](../extracted-lifts/c2-equal-angle-cone/), with
$\lambda=n=0$. In the active order $(d,m,C,A,D,\alpha,\beta)$, take

$$
\begin{aligned}
c_{\rm int}&=(1,20/29,12/29,13/58,51/58,3/29,3/58),\\
c_{\rm bdy}&=(1,5/6,3/10,1/3,2,1/13,0).
\end{aligned}\tag{1}
$$

Both have $A_0=B_0=Z$ and a contact at half-angle coordinates
$u=1,v=1/2$, with state

$$\psi_r=(1,-1/2,r,1)^T,\qquad
r=1/2\ \text{or}\ 3/8,$$

respectively. They are native C1 sources with
$x=1$, $\eta=\zeta=\theta=1/2$ in
[the rational dictionary](../latex/certificates/c1_short.tex).

Let $B_f(u,v)$ be their Bell operator, with
$A_1=[(1-u^2)Z+2uX]/(1+u^2)$ and similarly for $B_1,v$.
Put $P_f=(1+u^2)(1+v^2)B_f$. Its determinant has the exact form

$$\det P_f=\gamma_f(1+u^2)^2(1+v^2)^2G_f,$$

where $\gamma_{\rm int}=18432/707281$ and
$\gamma_{\rm bdy}=64/585$. To display $G_f$, set $U=u^2$, $V=v^2$,
$T=4V$ and define

$$M_1=(UT)^2+U+T-3UT,\qquad
M_2=U^2T+T^2+U-3UT.$$

Each is nonnegative by AM--GM: the product of its three positive terms
is $(UT)^3$. The complete determinant factors are

$$
\begin{aligned}
G_{\rm int}&=670V(U-1)^2+\frac{449}{2}U(T-1)^2
             +\frac{55}{2}M_1+255M_2,\\
G_{\rm bdy}&=235V(U-1)^2+\frac{153}{2}U(T-1)^2
             +\frac{295}{2}M_2.
\end{aligned}\tag{2}
$$

Thus both determinants are nonnegative at every real test. For noncommuting
measurements the CHSH core has three strictly positive eigenvalues;
adding the nonnegative event operators preserves them. The nonnegative
determinant makes the fourth eigenvalue nonnegative. When either pair
commutes, the core itself is positive semidefinite. Continuity includes
the projective angle endpoints. The two-observable block decomposition
then gives validity in arbitrary dimension.

## A universal word functional excludes both certificates

The two deterministic zeros are the assignments
$(A_0,A_1,B_0,B_1)=(1,1,1,1)$ and $(1,-1,-1,-1)$.
Their common kernel compresses any nine-word `1+AB` Gram to

$$R=(A_0-I,A_1B_0-I,A_1B_1-I,B_0-A_1,B_1-A_1,
      A_0B_0-A_1,A_0B_1-A_1)^T.$$

An NPA certificate therefore requires the universal identity
$f=R^\dagger QR$, $Q\succeq0$. Define a real linear word functional
$\ell$ by the following table and $\ell(w^\dagger)=\ell(w)$.
Every other word occurring in these products has value zero.

| $w$ | $\ell(w)$ | $w$ | $\ell(w)$ |
|---|---:|---|---:|
| $I$ | $11/2$ | $A_0A_1$ | $-28$ |
| $A_0$ | $-189$ | $B_0B_1$ | $-91$ |
| $A_1$ | $49/2$ | $A_0A_1B_0$ | $-133$ |
| $B_0$ | $7/2$ | $A_0A_1B_1$ | $-385/2$ |
| $B_1$ | $21$ | $A_0B_0B_1$ | $-259/2$ |
| $A_0B_0$ | $154$ | $A_0A_1B_0B_1$ | $329/2$ |
| $A_0B_1$ | $-105/2$ | $A_0A_1B_1B_0$ | $-315/2$ |
| $A_1B_0$ | $-147/2$ | $A_1B_1$ | $0$ |

Reducing only adjacent reflection squares and commuting opposite parties
gives $\ell(R_i^\dagger R_j)=K_{ij}$, with

$$K=\begin{pmatrix}
389&135&2&203&-21&-98&126\\
135&158&-12&42&0&-161&-84\\
2&-12&11&0&7&14&28\\
203&42&0&158&-12&23&142\\
-21&0&7&-12&11&9&9\\
-98&-161&14&23&9&277&240\\
126&-84&28&142&9&240&396
\end{pmatrix}.\tag{3}$$

Its leading principal minors are
$389,43237,412479,17674904,41152768,1682481152,73148743680$.
They are strictly positive, so $K\succ0$. Consequently every C1 source
admitted by this NPA cone satisfies

$$
0\le\operatorname{tr}(KQ)=\ell(f)
=-61d+9m+58C+44A+2D+121\alpha+100\beta.\tag{4}
$$

For the two sources in (1), the right side is respectively
$-43/29$ and $-3169/390$. Neither has a `1+AB` certificate.
Complex Grams do not change this conclusion: averaging a Gram with its
complex conjugate preserves positivity and the real universal identity.

## Both rays are extreme in the full Bell polar

Let $O=(I,A_0,A_1,B_0,B_1,A_0B_0,A_0B_1,A_1B_0,A_1B_1)$.
Every valid summand of either inequality in (1) vanishes at its two
deterministic zeros and at its quantum contact. Its Bell coefficients
therefore obey eight homogeneous equations: four operator-kernel
equations, two angular derivatives, and the two deterministic evaluations.
Explicitly, column $j$ of this system is

$$
\begin{pmatrix}
O_j\psi_r\\
\psi_r^T(\partial_uO_j)\psi_r\\
\psi_r^T(\partial_vO_j)\psi_r\\
O_j(1,1,1,1)\\
O_j(1,-1,-1,-1)
\end{pmatrix}_{u=1,\,v=1/2}.\tag{5}
$$

The determinants of its first eight columns are $-576/125$ and
$-273/100$, respectively. Thus each system has rank eight. Its kernel
is the line through the corresponding source, and every valid summand
must be a nonnegative multiple of that source. This proves extremality
without restricting the summands to one canonical coefficient orthant.

The examples require additional C1 certificates in any exact hybrid
construction retaining `1+AB`. The strict violation in the interior
example persists in an open native neighborhood. The calculation does
not classify every C1 source or supply their uniform finite lift.

Exact, self-contained checks reconstruct the
[global determinants and contact systems](checks/c1_rational_validity.py)
and the [universal word functional](checks/c1_npa_separator.py).

## Exclusion on an open dense strict C1 parameter set

A contact-kernel calculation gives a stronger generic statement. Consider
the additional source, in the same active order as (1),

$$c_*=(27/2,4,2,9,5/4,13/4,9/4),
\qquad (u_c,v_c;\psi_c)=(1,1;1,-1/2,1/4,1).\tag{6}$$

It is the exact seed EN9 of the native contact construction. Its coordinates
are $k=S=v_c=1$, $r=1/4$, $\theta=1/2$, hence
$\eta=\zeta=1$, $a=3/4$, $t=1/2$. The native gaps are
$E=1/2$, $Q=M=J=3/4$, $B_5=25/32$ and $B_\beta=15/16$;
all six retained event weights are strictly positive.

There is a direct global validity check. With $U=u^2,V=v^2$, the full
cleared state matrix has determinant

$$
\det P_*=1152(1+u^2)^2(1+v^2)^2UV
\left(\frac{494}{U}+\frac{210V}{U}+\frac{845}{V}
       +451V+520U+184UV-2704\right).
$$

The six coefficients sum to $2704$, and their exponent-weighted sums
in both $U$ and $V$ vanish. Weighted AM--GM therefore makes the bracket
nonnegative for all $U,V>0$. The positive three-state block, followed
by continuity at compatible measurements and projective endpoints,
proves global quantum validity as above. The four state-kernel equations
and both scalar angular derivatives vanish at the contact in (6).

Let $T$ have columns $R_i\psi_c$. Its first four columns have determinant
$27/32$, so $\operatorname{rank}T=4$. An exact basis of its kernel is

$$
N=\begin{pmatrix}
1&5/3&2\\
-7/4&-4/3&-17/4\\
-5/3&-16/9&-5\\
-3/4&-7/3&-13/4\\
1&0&0\\
0&1&0\\
0&0&1
\end{pmatrix}.\tag{7}
$$

Positivity and the zero expectation at the quantum contact force every
possible seven-root Gram to have the form $Q=NSN^T$, $S\succeq0$.
Reduce $R^\dagger NSN^TR$ in universal words, and take the coefficients of

$$
(B_0B_1,\ A_0B_0B_1,\ A_0A_1,\ A_0A_1B_0,
  \ A_0A_1B_0B_1,\ A_0A_1B_1).
$$

All six must vanish in a Bell functional. In coordinate order
$(S_{11},S_{22},S_{33},S_{12},S_{13},S_{23})$, their coefficient matrix
$M$ satisfies

$$\det M=\frac{7360}{243}\ne0.\tag{8}$$

Thus $S=Q=0$, contradicting the nonzero source. This is a second
exact exclusion argument, independent of the separating functional (4).

To extend it, use the five independent strict native contact coordinates

$$
p=(k,\ell,S,r,\theta),\qquad
(u_c,v_c)=(k,\ell),\qquad
\psi_c=(S,k(\theta-1),kSr,1).
$$

Here $\eta=\ell/(kS)$ and $\zeta=S\ell/k$; the native coefficient
formulas are rational in $p$. Strict positivity of the native gaps and
retained event weights defines an open five-dimensional parameter domain.
The contact action $T(p)$ is rational. Multiplication by
$(1+k^2)(1+\ell^2)$ gives a polynomial matrix $\widehat T(p)$.
Put $\Delta(p)=\det\widehat T_{[:,1:4]}$. Where $\Delta\ne0$, choose

$$
N(p)=\begin{pmatrix}
-\widehat T_{[:,1:4]}^{-1}\widehat T_{[:,5:7]}\\ I_3
\end{pmatrix}
$$

and form $M(p)$ by the same six universal coefficients. Its entries have
denominators dividing $\Delta^2$, so
$J(p)=\Delta(p)^{12}\det M(p)$ is a polynomial. At (6),
$\Delta=216$ and (8) holds; neither polynomial is identically zero.
Consequently $\Delta J\ne0$ is open and dense in every nonempty open
component of the strict parameter domain. Each such point is excluded
from `1+AB` by the same contact compression. No connectedness assumption
is needed: a nonzero polynomial cannot vanish on an open parameter set.
If the orientation $k\le1$ is imposed, the result holds on its open part
$k<1$ and is relatively dense up to $k=1$; all other seed inequalities
are strict, so the seed is approached from that open part.

This is a generic exclusion in the strict **contact parameter domain**,
not an open subset of the ambient source space. The selected exceptional
sets $\Delta=0$ or $J=0$ are not classified. A vanishing selected minor
does not imply admission, and this argument makes no assertion on the
literal $D=0$ or $\beta=0$ faces or on degenerate contacts. A uniform
additional C1 certificate cone remains unresolved.

The [exact seed checker](checks/c1_generic_seed.py) reconstructs (6)--(8),
the native gaps and the complete determinant identity. It uses the local
universal-word reduction helper and performs no numerical optimization.
