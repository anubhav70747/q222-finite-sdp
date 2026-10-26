# Equal-angle maximally entangled C2 sources lie outside W2

Let

$$
W_2=(I,A_0,A_1,B_0,B_1,A_0A_1,A_1A_0,B_0B_1,B_1B_0,
A_0B_0,A_0B_1,A_1B_0,A_1B_1).
$$

The whole positive equal-angle maximally entangled C2 family has no universal SOS in this fixed 13-word space. This includes its $n=0$ boundary. The obstruction below uses the reflection algebra, not equality only in one qubit representation.

## Sources and contact

Take $r,x>0$, $r+x+3rx\le1$, and put $z=x/r$. In the canonical slack $f=dH+4\sum_i a_ip_i$, set $\lambda=\beta=0$ and

$$
\begin{aligned}
d&=1,&
n&=\frac{1-r-x-3rx}{4(r+x+2rx)},\\
m&=\frac{(1-r)(1+x)}{4(r+x+2rx)},&
A&=\frac{(1+r)(1+x)}{4(r+x+2rx)},\\
D&=\frac{(1+r)(1-x)}{4(r+x+2rx)},&
C&=\frac{r^2(1+x)}{(1-r)(r+x+2rx)},\\
\alpha&=\frac{x^2(1+r)}{(1-x)(r+x+2rx)}.
\end{aligned}
$$

These are the native C2 formulas at $a=b=1$, $t_{\rm source}=1-r$, $c=\sqrt z$. All native gaps are strict: $0<r<1$, $0<x<1$, $h=r$, and $g=j=r(1+z)>0$. All retained event weights are positive except that $n$ may vanish.

Their global validity follows directly from the existing [order-seven source cone](../extracted-lifts/c2-equal-angle-cone/). Indeed, with $D_0=r+x+2rx$,

$$
A-m=\frac{r(1+x)}{2D_0}>0,\qquad
A-D=\frac{x(1+r)}{2D_0}>0,\qquad
C=\frac{(A-m)^2}{m},\quad \alpha=\frac{(A-D)^2}{D},
$$

and $n=m+D-A\ge0$, $d=4A-2m-2D=1$. Thus both source blocks $\left(\begin{smallmatrix}m&A-m\\A-m&C\end{smallmatrix}\right)$ and $\left(\begin{smallmatrix}D&A-D\\A-D&\alpha\end{smallmatrix}\right)$ are rank-one PSD, including $n=0$. The native stationary-contact theorem supplies the contact used below.

With

$$
M(t)=\frac{(1-t^2)\sigma_z+2t\sigma_x}{1+t^2},
$$

the actual contact is

$$
A_0=B_0=\sigma_z,\quad A_1=B_1=M(\sqrt z),\quad
\psi_c=(1,-r\sqrt z,r\sqrt z,1)^T.
$$

The state is homogeneous; after normalization it is maximally entangled.

## The entire W2 Gram reduces to seven coordinates

Every source above vanishes at the deterministic assignments $(1,1,1,1)$ and $(-1,1,-1,1)$ in order $(A_0,A_1,B_0,B_1)$. Their independent evaluation vectors have an eleven-dimensional common annihilator, with fixed basis

$$
\begin{aligned}
h=(\;&A_1-I,\ B_1-I,\ A_0B_0-I,\ A_1B_1-I,\\
&B_0-A_0,\ A_0B_1-A_0,\ A_1B_0-A_0,\\
&A_0A_1-A_0,\ A_1A_0-A_0,\ B_0B_1-B_0,\ B_1B_0-B_0\;)^T.
\end{aligned}
$$

Indeed each root has its own nonpivot word with coefficient one. A PSD $W_2$ Gram annihilates both evaluation vectors, so any certificate compresses to $f=h^\dagger Qh$, $Q\succeq0$. Complex Grams can be averaged with their conjugates, giving real symmetric $Q$ for the same real universal target.

At the contact let $T=(h_1\psi_c,\ldots,h_{11}\psi_c)$. Its columns $(1,2,3,5)$ have determinant

$$
\det T_{[:,(1,2,3,5)]}=-\frac{64r^3z^2}{1+z}\ne0.
$$

Put $R=r(1+z)$, $p=(1+r)(1-rz)$, $q=1+r^2z$, and $s=(r-1)(1+rz)$. An exact basis of $\ker T$ is

$$
N=\frac1R\begin{pmatrix}
0&p&q&q&p&p&q\\
0&-q&s&s&-q&-q&s\\
-R&-p&-s&-q&q&q&-q\\
R&0&0&0&0&0&0\\
0&0&-R&-R&0&R&0\\
0&R&0&0&0&0&0\\
0&0&R&0&0&0&0\\
0&0&0&R&0&0&0\\
0&0&0&0&R&0&0\\
0&0&0&0&0&R&0\\
0&0&0&0&0&0&R
\end{pmatrix}.
$$

Direct Pauli action gives $TN=0$; rows $(4,6,7,8,9,10,11)$ form $I_7$. Zero contact expectation forces $TQ=0$, so $Q=NSN^T$, $S\succeq0$. Consequently every possible certificate must satisfy

$$
f=(N^Th)^\dagger S(N^Th)
$$

as a universal identity. Some roots are not Hermitian: the adjoint reverses same-party word order and must be retained.

## Linear consequences force a PSD contradiction

Universal coefficient matching gives

$$
\begin{gathered}
S_{11}=S_{56}-2S_{57}-S_{77},\qquad S_{55}=-S_{56},\\
S_{15}+S_{17}=0,\qquad S_{16}=S_{67}=0,\\
S_{35}=\frac{r(1+r)(1+z)^2}{32(1-r)(1+z+2rz)}>0.
\end{gathered}
$$

Their bounded triangular derivation is given below. For $v=e_1+e_5+e_7$, these equations imply $v^TSv=0$. Positivity forces $Sv=0$, whose sixth component is

$$
(Sv)_6=S_{16}+S_{56}+S_{67}=S_{56}=-S_{55}.
$$

Thus $S_{55}=0$. A PSD matrix with a zero fifth diagonal entry has zero fifth row, contradicting $S_{35}>0$. This excludes every $W_2$ Gram.

At $r=1/4,z=1$ the source is

$$
(d,n,m,C,A,D,\alpha)=(1,1/8,3/8,1/6,5/8,3/8,1/6),
\qquad S_{35}=1/48.
$$

At $r=1/3,z=1$ it is

$$
(d,n,m,C,A,D,\alpha)=(1,0,1/4,1/4,1/2,1/4,1/4),
\qquad S_{35}=1/32.
$$

No step divides by $n$.

## Reproducing the coefficient elimination

Here is a finite prescription for every coefficient and each elimination step. Reduce a universal word to the pair of its alternating Alice and Bob words, cancelling adjacent equal letters within each party and commuting only across parties. For an adjoint, reverse both same-party words.

Put $K=RN$, $\ell=K^Th$, and

$$
\mathcal D=4r(1-r)(1-rz)(1+z+2rz),\qquad U_{ij}=\mathcal D S_{ij}.
$$

For each reduced word $w$, set

$$
A_{w,ii}=[w](\ell_i^\dagger\ell_i),\qquad
A_{w,ij}=[w](\ell_i^\dagger\ell_j+\ell_j^\dagger\ell_i)\quad(i<j).
$$

There are 41 word rows. Their equations are

$$
\sum_{i\le j}A_{w,ij}U_{ij}=R^2\mathcal D\, b_w.
$$

The target coefficients $b_w$ are the following nine, and zero for all other words:

$$
\begin{array}{ll}
b_I=1+n+m+C+A+D+\alpha,&b_{A_0}=n-A-\alpha,\\
b_{A_1}=m-C-D,&b_{B_0}=-n+C+A,\\
b_{B_1}=-m+D-\alpha,&b_{A_0B_0}=-1/2-n-A,\\
b_{A_0B_1}=1/2+\alpha,&b_{A_1B_0}=-1/2-C,\\
b_{A_1B_1}=-1/2-m-D.&
\end{array}
$$

Use upper-triangle variable order $(U_{11},U_{12},\ldots,U_{17},U_{22},\ldots,U_{77})$. In each table row, take the equation for the named word after eliminating the preceding pivot variables, divide by the displayed pivot, and eliminate that variable from every remaining equation. Back substitution after these 23 steps gives the displayed relations; it also gives $S_{14}=S_{45}=0$ and $S_{47}=-S_{77}$. The five free variables are $U_{44},U_{56},U_{57},U_{66},U_{77}$.

| Step | Word equation | Pivot variable | Pivot coefficient |
|---|---|---|---|
| 1 | $B_0B_1B_0B_1$ | $U_{67}$ | $r^{2} \left(z + 1\right)^{2}$ |
| 2 | $A_0A_1A_0A_1$ | $U_{45}$ | $r^{2} \left(z + 1\right)^{2}$ |
| 3 | $A_0B_1B_0B_1$ | $U_{26}$ | $2 r^{2} \left(z + 1\right)^{2}$ |
| 4 | $A_0A_1A_0B_1$ | $U_{25}$ | $2 r^{2} \left(z + 1\right)^{2}$ |
| 5 | $A_1B_0B_1B_0$ | $U_{37}$ | $2 r^{2} \left(z + 1\right)^{2}$ |
| 6 | $A_1B_1B_0B_1$ | $U_{16}$ | $2 r^{2} \left(z + 1\right)^{2}$ |
| 7 | $A_1A_0A_1B_0$ | $U_{34}$ | $2 r^{2} \left(z + 1\right)^{2}$ |
| 8 | $A_1A_0A_1B_1$ | $U_{14}$ | $2 r^{2} \left(z + 1\right)^{2}$ |
| 9 | $B_0B_1B_0$ | $U_{47}$ | $- 2 r^{2} \left(z + 1\right)^{2}$ |
| 10 | $A_0A_1A_0$ | $U_{55}$ | $- 2 r^{2} \left(z + 1\right)^{2}$ |
| 11 | $A_0A_1B_1B_0$ | $U_{23}$ | $r^{2} \left(z + 1\right)^{2}$ |
| 12 | $A_0B_0B_1B_0$ | $U_{17}$ | $- 2 r^{2} \left(z + 1\right)^{2}$ |
| 13 | $A_0A_1A_0B_0$ | $U_{15}$ | $- 2 r^{2} \left(z + 1\right)^{2}$ |
| 14 | $B_1B_0B_1$ | $U_{36}$ | $2 r \left(r - 1\right) \left(z + 1\right) \left(r z + 1\right)$ |
| 15 | $A_1A_0A_1$ | $U_{24}$ | $- 2 r \left(r + 1\right) \left(z + 1\right) \left(r z - 1\right)$ |
| 16 | $A_0B_1$ | $U_{22}$ | $4 r \left(z + 1\right) \left(r^{2} z + 1\right)$ |
| 17 | $A_0A_1B_1$ | $U_{12}$ | $- r^{2} \left(z + 1\right)^{2}$ |
| 18 | $A_0B_0$ | $U_{46}$ | $4 \left(r^{2} z + 1\right)^{2}$ |
| 19 | $A_0B_0B_1$ | $U_{35}$ | $\left(8 r - 8\right) \left(r z + 1\right)$ |
| 20 | $A_1B_0B_1$ | $U_{13}$ | $- 2 r^{2} \left(z + 1\right)^{2}$ |
| 21 | $A_1$ | $U_{33}$ | $- 4 r^{2} \left(z + 1\right)^{2}$ |
| 22 | $A_1B_1$ | $U_{27}$ | $\frac{\left(8 - 8 r\right) \left(r + 1\right) \left(r z - 1\right) \left(r z + 1\right)}{r^{2} z + 1}$ |
| 23 | $A_0A_1B_0B_1$ | $U_{11}$ | $- r^{2} \left(z + 1\right)^{2}$ |

Every pivot is nonzero throughout the entire domain. Its factors lie among

$$
r,\quad1-r,\quad1+r,\quad1-rz,\quad1+rz,\quad1+z,\quad1+r^2z.
$$

Back-substitution denominators additionally use $1+z+2rz>0$. Therefore the displayed consequences hold at every allowed parameter, including $n=0$, rather than only at generic points. The source denominators are nonzero for the same reason.

The [standalone exact checker](checks/c2_w2_exclusion.py) reconstructs the universal coefficient matrix independently from the words, the Pauli contact action and the native source formulas. Its [output](checks/c2_w2_exclusion.json) records all pivot factors and identities. It uses exact symbolic arithmetic with a bounded sparse elimination, no numerical solver or parameter sweep.

The existing [exact order-seven source cone](../extracted-lifts/c2-equal-angle-cone/) contains this whole family, so appending its dual block strictly sharpens $W_2$. Other certificate cones could also supply these inequalities. This result concerns one two-parameter family and does not establish an exact lift for all C2 sources or for $\mathcal Q_{222}$.
