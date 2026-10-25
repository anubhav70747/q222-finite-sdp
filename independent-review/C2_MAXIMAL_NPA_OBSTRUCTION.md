# Maximally entangled C2 contacts require additional certificates

Every finite noncommuting maximally entangled contact in the native C2
family lies outside the fixed `1+AB` certificate cone. This includes unequal
measurement angles and the face $n=0$. The obstruction follows from six
universal word coefficients after the contact reduces the Gram to order three.

## Contact and domain

Use the [canonical slack and event order](../extracted-lifts/c2-equal-angle-cone/):
$f=dH+4\sum_i a_i p_i$, with $d>0$ and $\lambda=\beta=0$.
Put $M(t)=((1-t^2)Z+2tX)/(1+t^2)$ and take

$$
A_0=B_0=Z,\qquad A_1=M(\sqrt z),\qquad B_1=M(a\sqrt z),\qquad
\psi=(1,-r\sqrt z,r\sqrt z,1)^T.\tag{1}
$$

The finite native domain satisfies

$$a,z>0,\qquad 0<r<1,\qquad arz<1,\qquad
K=a+r-1+arz>0.\tag{2}$$

Here and below tensor factors are understood. The native parameters are
$c=\sqrt z$, $b=a$ and $t=(1-r)/a$. The source inequalities further
restrict (2); the obstruction holds throughout (2) whenever the slack has
the contact (1). In particular, it imposes no strict inequality on $n$.

## The contact leaves only three Gram directions

Every C2 slack vanishes at the deterministic assignments
$(A_0,A_1,B_0,B_1)=(1,1,1,1)$ and $(-1,1,-1,1)$.
A positive semidefinite nine-word Gram must annihilate both evaluation
vectors. Their common kernel has the following seven-root basis:

$$
h=(A_1-I,\ B_1-I,\ A_0B_0-I,\ A_1B_1-I,\ B_0-A_0,
       \ A_0B_1-A_0,\ A_1B_0-A_0)^T.
$$

Thus a certificate would have the universal identity $f=h^\dagger Qh$,
$Q\succeq0$. Averaging a complex Gram with its conjugate permits $Q$ to be
real. At (1), define $T=(h_1\psi,\ldots,h_7\psi)$. Zero expectation gives
$TQT^T=0$, hence $TQ=0$.

For an explicit kernel calculation put $U=1+z$, $V=1+a^2z$. Then
$T=\operatorname{diag}(2z,2\sqrt z,2\sqrt z,2)T_0$, where

$$
T_0=\begin{pmatrix}
\frac{r-1}{U}&-\frac{a(a+r)}V&0&-\frac{(a-1)K}{UV}&0&-\frac{a(a+r)}V&\frac{r-1}U\\
\frac{1+rz}U&\frac{a+r}V&r&\frac{(1+az)K}{UV}&r&\frac{a+r}V&\frac{r-1}U\\
\frac{1-r}U&-\frac{a(arz-1)}V&-r&-\frac{(1+az)K}{UV}&r&\frac{a(arz-1)}V&\frac{1+rz}U\\
-\frac{1+rz}U&\frac{arz-1}V&0&-\frac{z(a-1)K}{UV}&0&-\frac{arz-1}V&\frac{1+rz}U
\end{pmatrix}.\tag{3}
$$

The columns $P=(1,2,3,5)$ have determinant

$$\det(T_{0,P})=-\frac{2r^2(a+1)K}{UV}\ne0.\tag{4}$$

Consequently the kernel has dimension three. An explicit basis $N$ is
specified without a choice of normalization by

$$N_J=I_3,\qquad N_P=-T_{0,P}^{-1}T_{0,J},\qquad J=(4,6,7).\tag{5}$$

Every possible contact Gram is therefore $Q=NSN^T$, $S\succeq0$.

## Six coefficients force this Gram to vanish

The Bell slack contains only the identity, single observables and one
observable from each party. Reducing products using $A_x^2=B_y^2=I$ and
$[A_x,B_y]=0$, the following six forbidden words have coefficients

| Word | Coefficient in $h^TQh$ |
|---|---|
| $B_0B_1$ | $Q_{25}+Q_{36}+Q_{47}$ |
| $A_0A_1$ | $-Q_{15}-Q_{16}-Q_{17}+Q_{37}+Q_{46}$ |
| $A_0A_1B_0$ | $Q_{13}-Q_{57}-Q_{67}-Q_{77}$ |
| $A_0A_1B_0B_1$ | $Q_{34}$ |
| $A_0A_1B_1$ | $Q_{16}-Q_{45}-Q_{46}-Q_{47}$ |
| $A_0A_1B_1B_0$ | $Q_{67}$ |

All six must be zero. Substitute (5), and let $M$ be their coefficient
matrix in the order
$q=(S_{11},S_{22},S_{33},S_{12},S_{13},S_{23})^T$.
This completely specifies $M$ as a rational six-by-six matrix. Its determinant is

$$
\det M=
\frac{32a^2(a+r)^2(1-r)(1+rz)(1-arz)^2
       [1+z(a+r-ar)]}
     {r^3(a+1)^3(1+z)(1+a^2z)^3}.\tag{6}
$$

The calculation can be reduced to order four. The last equation is
$S_{23}=0$. The fourth equation is
$N_{31}S_{11}+N_{32}S_{12}+N_{33}S_{13}=0$, with

$$N_{31}=-\frac{K[1+z(a+r-ar)]}{rUV}\ne0.$$

Eliminating these two variables gives $\det M=-N_{31}\det M_4$.
The [exact checker](checks/c2_maximal_npa.py) reconstructs all six
coefficients from the universal words and records the row and column factors
of this four-by-four calculation in [its output](checks/c2_maximal_npa.json).

Every factor in (6) is positive on (2), since
$a+r-ar=a(1-r)+r>0$. Therefore $Mq=0$ forces $S=0$, and hence $Q=0$.
This contradicts $d>0$. The contradiction excludes even an indefinite
$S$ after contact compression; positivity was needed to obtain that
compression from a nine-word certificate.

## Consequence for the finite lift

The additional cones must cover the entire finite maximally entangled C2
contact family. Equal angles $a=1$ and $n=0$ introduce no exception to
(6). The [fixed-fiber pencils](C2_FIXED_FIBERS.md) supply exact certificates
on each of their specified source cones; combining their continuum into
one fixed affine lift remains open.

The conditions (2) exclude degenerations such as $r=0$, $r=1$, $arz=1$,
$z=0$, $a=0$ or $K=0$. Nonmembership does not pass to these limits:
the complement of the closed NPA cone need not be closed. General partially
entangled contacts are covered only by the separate
[generic obstruction](C2_NONMAXIMAL_NPA_OBSTRUCTION.md).
