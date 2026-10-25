# An exact small block beyond NPA 1+AB

A continuous family of tight quantum inequalities has no certificate in the
nine-word NPA `1+AB` space. Nevertheless, all inequalities in this family can
be imposed together with one 2-by-2 PSD matrix and one nonnegative scalar.
The construction works because a positive rescaling makes their coefficients
quadratic in a single parameter.

The resulting block sharpens `1+AB` and is valid on the complete quantum set.
It represents this family exactly; coverage of the remaining C1/C2 sources
is still required for a representation of the full quantum set.

## The inequalities

Let the binary observables be reflections, with Alice's observables commuting
with Bob's. Write

$$H=1-\frac{E_{00}-E_{01}+E_{10}+E_{11}}2,$$

and order the eight event probabilities as

$$
(p_\lambda,p_n,p_m,p_C,p_A,p_D,p_\alpha,p_\beta)
=\bigl(p(++|01),p(+-|00),p(+-|11),p(-+|10),
p(-+|00),p(-+|11),p(--|01),p(+-|10)\bigr).
$$

A source vector $c=(d,\lambda,n,m,C,A,D,\alpha,\beta)$ specifies the affine
functional $f_c=dH+4\sum_i c_i p_i$, where the sum runs over the eight event
weights. Put $b=2/5$. For $0<z\le b$, consider

$$
\begin{gathered}
d=1,\qquad \lambda=\beta=0,\qquad
n=\frac{2-5z}{8(1+2z)},\qquad m=\frac{2+z}{8(1+2z)},\\
C=4m,\qquad A=3m,\qquad
D=\frac{3(2-z)}{8(1+2z)},\qquad
\alpha=\frac{3z^2}{2(2-z)(1+2z)}.
\end{gathered}
$$

These functionals are nonnegative on every quantum behaviour and attain
zero at explicit real two-qubit strategies. With $s=\sqrt z$, choose

$$
A_0=B_0=\sigma_z,\qquad
A_1=B_1=\frac{(1-s^2)\sigma_z+2s\sigma_x}{1+s^2},\qquad
|\psi_s\rangle\ \propto\ (1,-s/2,s/2,1)^T.
$$

The endpoint $z=0$ is included by continuity. The proof of validity and the
obstruction to a `1+AB` certificate appear below.

## The exact source cone

Multiply $c(z)$ by the positive factor $8(2-z)(1+2z)$. Its ray is unchanged,
and its coefficients become

$$
\widetilde c(z)=
\begin{pmatrix}
16+24z-16z^2\\0\\4-12z+5z^2\\4-z^2\\16-4z^2\\
12-3z^2\\12-12z+3z^2\\12z^2\\0
\end{pmatrix}
=T\begin{pmatrix}1\\z\\z^2\end{pmatrix},\qquad
T=\begin{pmatrix}
16&24&-16\\0&0&0\\4&-12&5\\4&0&-1\\16&0&-4\\
12&0&-3\\12&-12&3\\0&0&12\\0&0&0
\end{pmatrix}.
$$

Let $K$ be the closed conic hull of these sources for $0\le z\le b$.
Then

$$
c\in K\quad\Longleftrightarrow\quad
\exists h_0,h_1,h_2:\quad c=T(h_0,h_1,h_2)^T,\quad
\begin{pmatrix}h_0&h_1\\h_1&h_2\end{pmatrix}\succeq0,\quad
bh_1-h_2\ge0.\tag{1}
$$

To see both inclusions, interpret $h_j$ as the first three moments of a
nonnegative measure on $[0,b]$. Such moments satisfy the two constraints.
Conversely, if $h_2>0$, those constraints give

$$
t=\frac{h_2}{h_1}\in(0,b],\qquad
w=\frac{h_1^2}{h_2}\le h_0.
$$

The two atoms $w\delta_t+(h_0-w)\delta_0$ have exactly the prescribed
moments. If $h_2=0$, PSD forces $h_1=0$, and the single atom $h_0\delta_0$
suffices. Thus two source rays already generate every point of $K$.

There is no additional closure hidden in (1). The map $T$ is injective:

$$h_2=\alpha/12,\qquad h_0=m/4+\alpha/48,\qquad
h_1=(d-4m+\alpha)/24.$$

It therefore maps the closed moment cone homeomorphically onto its image.

## One block imposes the whole family

For any candidate behaviour $q$, define the three affine quantities

$$
\begin{aligned}
a_0(q)&=16H+16p_n+16p_m+64p_C+48p_A+48p_D,\\
a_1(q)&=24H-48p_n-48p_D,\\
a_2(q)&=-16H+20p_n-4p_m-16p_C-12p_A+12p_D+48p_\alpha.
\end{aligned}
$$

They are the coefficients of $f_{\widetilde c(z)}(q)=a_0+a_1z+a_2z^2$.
Consequently $q$ obeys every inequality in the family exactly when

$$
\exists\tau\ge0:\qquad
\begin{pmatrix}
a_0&(a_1-b\tau)/2\\
(a_1-b\tau)/2&a_2+\tau
\end{pmatrix}\succeq0.\tag{2}
$$

Here is an elementary proof that also constructs the auxiliary variable.
For any PSD matrix $Q$ in (2),

$$p(z)=(1,z)Q(1,z)^T+\tau z(b-z)$$

is nonnegative on the interval. Conversely, suppose $p\ge0$ there and set
$t_0=(2a_0+ba_1)/b^2$. If $t_0\ge0$, choose $\tau=t_0$. The required matrix is

$$
Q=a_0\begin{pmatrix}1\\-1/b\end{pmatrix}
\begin{pmatrix}1&-1/b\end{pmatrix}
+p(b)\begin{pmatrix}0\\1/b\end{pmatrix}
\begin{pmatrix}0&1/b\end{pmatrix}\succeq0.
$$

If $t_0<0$, then $a_1<0$, $a_2>0$, and the vertex of $p$ lies in $(0,b)$.
Its minimum is therefore nonnegative on the whole real line. Taking
$\tau=0$ gives $Q\succeq0$. This proves (2), including singular endpoints,
without a dual-attainment assumption.

Thus the family uses total PSD order three: a 2-by-2 block and the scalar
block $\tau$. Both its source cone and its constraints on behaviours are
explicit. Finite physical relabellings can be added as further copies.

## Quantum validity from four AM–GM terms

For arbitrary real measurement coordinates $u,v$, replace $s$ by $u$ in
Alice's observable and by $v$ in Bob's. Write $D_u=1+u^2$, $D_v=1+v^2$,
and use $(W,X,Y,Z)$ for an unnormalised real state. The cleared expectation
$D_uD_v\|\psi\|^2f_{\widetilde c(z)}$ is the quadratic form

$$
\begin{aligned}
F={}&4nD_uD_vX^2+4m(vW-X+uvY-uZ)^2\\
&+4CD_v(uW-Y)^2+4AD_uD_vY^2\\
&+(4D+2d)(uW+uvX-Y-vZ)^2\\
&+(4\alpha+2d)D_u(vY-Z)^2+2dD_v(X+uZ)^2-2dD_uD_vZ^2,
\end{aligned}
$$

with the **scaled** weights $\widetilde c(z)$. Let $P$ be its state matrix.
When $uv\ne0$, its $Z=0$ principal block is positive definite: the $d$ terms
alone give

$$2d\bigl[(uW+uvX-Y)^2+D_uv^2Y^2+D_vX^2\bigr].$$

It remains to determine the sign of the one remaining eigenvalue, which is
the sign of the determinant. Set $x=u^2/z$, $y=v^2/z$ and
$G=432(2-z)^3(z+2)^3$. Direct expansion gives

$$\det P=16G D_u^2D_v^2u^2v^2\,R_z(x,y),$$

where

$$
\begin{aligned}
R_z={}&16z(x+x^{-1}-2)
+16z(2-z)(y+y^{-1}-2)\\
&+3(z+2)(2-5z)(x/y+y/x-2)\\
&+16z(z+1)(y/x+x+y^{-1}-3).
\end{aligned}
$$

Every coefficient is nonnegative on $0<z\le2/5$. Each parenthesis is
nonnegative by AM–GM; the three factors in the last one have product one.
Thus $P\succeq0$ at every finite nonzero $u,v$. Continuity includes zero
and infinite coordinates after restoring the bounded observable matrices.
The usual two-reflection block decomposition reduces arbitrary-dimensional
strategies to these qubit blocks. Complex states give the sum of the
expectations on their real and imaginary parts, and mixtures preserve
nonnegativity. This proves global quantum validity. The displayed state
at $u=v=\sqrt z$ has zero expectation, as direct substitution confirms.

This calculation also shows why the endpoint $n=0$ is harmless: only the
third AM–GM term disappears there.

## Why 1+AB misses it

The nine-word list is
$I,A_0,A_1,B_0,B_1,A_0B_0,A_0B_1,A_1B_0,A_1B_1$.
The two deterministic zeros, $(A_0,A_1,B_0,B_1)=(1,1,1,1)$ and
$(-1,1,-1,1)$, force every square root of a `1+AB` certificate into the
span of

$$
h=(A_1-I,\ B_1-I,\ A_0B_0-I,\ A_1B_1-I,\
B_0-A_0,\ A_0B_1-A_0,\ A_1B_0-A_0)^T.
$$

Hence such a certificate would have the form $f=h^TQh$, $Q\succeq0$.
The nonlocal zero further forces the range of $Q$ into the kernel of the
four-by-seven matrix with columns $h_i\psi_s$. Put

$$a=\frac{3(2-z)}{2(1+z)},\quad
b_0=-\frac{z+4}{2(1+z)},\quad d_0=-b_0,\quad
e=-\frac{z+2}{2(1+z)}.$$

This kernel has basis

$$
N=\begin{pmatrix}
0&a&d_0\\0&b_0&e\\-1&-a&-e\\1&0&0\\0&0&-1\\0&1&0\\0&0&1
\end{pmatrix}.
$$

Indeed $(h_1\psi_s,h_2\psi_s,h_3\psi_s,h_5\psi_s)$ has determinant
$-8s^4/(1+s^2)\ne0$, and multiplying the matrix of root actions by $N$
gives zero. Thus $Q=NSN^T$ for a symmetric
3-by-3 matrix $S$.

Coefficient comparison in the universal reflection algebra gives four
homogeneous identities:

$$Q_{34}=Q_{67}=0,\qquad Q_{13}=Q_{57}+Q_{77},\qquad Q_{23}=-Q_{56}.$$

Since $ae-b_0d_0=1$, these identities imply

$$S_{23}=0,\qquad S_{12}=-aS_{22},\qquad
S_{13}=-eS_{33},\qquad S_{11}=a^2S_{22}+e^2S_{33}.$$

Writing $k=N^Th$, the proposed certificate consequently reduces to just
two squares:

$$S_{22}(k_2-ak_1)^2+S_{33}(k_3-ek_1)^2.$$

Their roots factor as

$$
k_2-ak_1=(I-B_1)(aA_1-A_0-b_0I),\qquad
k_3-ek_1=(A_1-I)(B_0-eB_1+d_0I).
$$

The coefficients of the words $A_0A_1$ and $B_0B_1$ are respectively
$-2aS_{22}$ and $-2eS_{33}$. Neither word occurs in a Bell functional.
Because $a>0$ and $e<0$, both coefficients vanish only when
$S_{22}=S_{33}=0$. The whole Gram matrix would then vanish, contradicting
$d=1$. This excludes every $0<z\le2/5$, including $n=0$.

For example, $z=1/4$ gives the rational source

$$(d,\lambda,n,m,C,A,D,\alpha,\beta)
=(1,0,1/16,3/16,3/4,9/16,7/16,1/28,0).$$

The exclusion also implies a strict improvement of the outer relaxation.
In the full word-product coefficient space, the image of PSD nine-word
Grams is closed: the identity coefficient equals the Gram trace and bounds
every entry. Separation supplies a PSD NPA moment matrix on which the
excluded functional is negative. Its identity moment is positive, since
a zero identity moment would make every diagonal, hence the whole moment
matrix, zero. Normalizing it produces a `1+AB` behaviour excluded by (2).

The positive-parameter sources are also extreme rays of the canonical cone
of valid nonnegative sources. Let $B_c(u,v)$ denote the uncleared Bell
operator of $f_c$. At the displayed strategy, use the four entries of
$B_c\psi_s$ followed by $\psi_s^T(\partial_uB_c)\psi_s$ and
$\psi_s^T(\partial_vB_c)\psi_s$, with the unnormalised state above. These
six equations are linear in
$(d,n,m,C,A,D,\alpha)$. The minor obtained by deleting the $\alpha$ column
is $-768s^{10}/(1+s^2)^5$, so their common solution space is one-dimensional.
In any decomposition into valid nonnegative sources, the two omitted
weights remain zero in every summand. Each summand must vanish at the
contact and have zero first derivatives there, so it belongs to that same
one-dimensional space. Every summand is therefore proportional to $c(z)$.
This includes $z=2/5$.

## Reproduction

The supplied symbolic checks use rational arithmetic and the universal
relations $A_x^2=B_y^2=I$ and $[A_x,B_y]=0$. They do not impose commutation
between measurements on the same party. The accompanying validity proof
uses the explicit state matrix, its determinant and the contact equations.

Run the four short checks with Python and SymPy:

```text
python npa_gram.py
python verify_npa_exclusion.py
python verify_physicality.py
python verify_curve_lift.py
```

The JSON files contain the coefficient maps and exact identity results.
Equations (1) and (2) supply the complete matrices; the proofs above establish
their domains and both inclusions.
