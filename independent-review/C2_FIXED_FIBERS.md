# Exact convexification at fixed angle and contact parameters

The maximally entangled C2 sources admit a simple convexification after two
parameters are fixed. Their remaining contact curve becomes a hyperbola.
Its closed conic hull is one $2$-by-$2$ PSD constraint and two scalar
inequalities, with no auxiliary variables. Every point is a sum of at most
two contact rays or their endpoint limits.

Use the active source order $(d,n,m,C,A,D,\alpha)$, with
$\lambda=\beta=0$, and the native dictionary in
[C2_CONTACT_GEOMETRY.md](C2_CONTACT_GEOMETRY.md). Set

$$x=\frac{a+r}{1-r},\qquad y=\frac{a+s}{1-s},\qquad k=xy.$$

Here $k$ denotes this product throughout. Its inverse is
$r=(x-a)/(x+1)$, $s=(y-a)/(y+1)$. Fix constants

$$a>0,\qquad k>\max(1,a^2),$$

and put

$$w_- =\frac{k(k+a^2-2a)}{a(k-a)},\qquad w_+=a+k/a.$$

The native parameter domain is exactly

$$xy=k,\qquad x,y>a,\qquad x+y\ge w_- .\tag{1}$$

Indeed, the two native sign conditions become

$$g=\frac{(a+1)(k-1)}{(x+1)(y+1)},\qquad
B_0=\frac{(a+1)^3a(k-a)(w_--x-y)}{(x+1)^2(y+1)^2}.$$

In particular, $x+y<w_+$. Equality at $w_+$ is allowed in the closure
and corresponds to $r=0$ or $s=0$. Every stated pair $(a,k)$ has a
nonempty native curve, since
$w_+-w_-=(k-a^2)/(k-a)>0$ and $w_+>2\sqrt{k}$.

## The source pencil

Define the linear source expressions

$$
X=\frac{k[(a-1)A+a(k-1)D]}{k-a},\qquad
Y=\frac{k[(a-1)A+a(k-1)m]}{k-a}.
$$

The closed conic hull of the native sources in (1) is exactly

$$
\begin{gathered}
a(d+2n)=2A,\\
C+d/2=(k-1)D,\qquad \alpha+d/2=(k-1)m,\\
(k-a^2)(A+n)=a^2(k-1)(m+D),\\[1mm]
\begin{pmatrix}X&kA\\kA&kY\end{pmatrix}\succeq0,\qquad
n\ge0,\qquad aA-kn\ge0.\tag{2}
\end{gathered}
$$

All entries are affine in the source coefficients because $a,k$ are fixed.
Joining the three PSD blocks gives an order-four pencil. The other four
constraints are linear equalities. Additional coordinatewise positivity
constraints are unnecessary; (2) already implies them.

To prove the forward inclusion, direct substitution of the native
coefficients gives $X=Ax$, $Y=Ay$ and the four linear equalities. The PSD
block has rank one. If
$c_w=a(k-a)/[k(k-a^2)]>0$, the scalar conditions are

$$n=c_w A(x+y-w_-),\qquad
aA-kn=kc_w A(w_+-x-y).\tag{3}$$

They follow from (1). Continuity includes both zero-opening endpoints.

For the converse, the scalar inequalities first imply $A\ge0$.
If $A=0$, they force $n=d=0$; the fourth equality and PSD diagonals then
force every other coefficient to vanish. Suppose $A>0$ and write
$w=(X+Y)/A$. The linear equalities give (3) with $x+y$ replaced by $w$,
while the PSD block gives

$$X,Y>0,\qquad XY\ge kA^2,\qquad
\max(w_-,2\sqrt{k})\le w\le w_+.$$

Define

$$x_\pm=\frac{w\pm\sqrt{w^2-4k}}2.$$

Then $x_+x_-=k$ and $x_\pm\ge a$. Indeed,
$(a-x_-)(a-x_+)=a(w_+-w)\ge0$, while $x_+>a$.
Moreover,

$$
(X/A-x_-)(x_+-X/A)=XY/A^2-k\ge0.
$$

Thus $(X/A,Y/A)$ lies on the segment joining
$(x_+,x_-)$ to $(x_-,x_+)$. When the endpoints differ, its weight on the
first is $(X/A-x_-)/(x_+-x_-)$. At a double root, one point suffices.
The four linear equalities uniquely recover the entire source from
$A,X,Y$. The same convex combination therefore reconstructs every
source coefficient. Both endpoints satisfy (1) or are limits at $w=w_+$.
This proves (2), including all coefficient signs and finite endpoints.

Every nonzero point has

$$d=\frac{2(k-a^2)A}{ak}+\frac{2(aA-kn)}k>0.$$

The boundary $n=0$ is retained whenever $w_-\ge2\sqrt{k}$.
There is also strict feasibility in the three-dimensional source subspace:
choose $A=1$, $X=Y=w/2$ and
$\max(w_-,2\sqrt{k})<w<w_+$. All three PSD blocks are then positive
definite. Consequently the usual SDP polar representation is attained
after imposing the four source equalities.

## A short converse for individual generators

On the four linear equalities in (2), let $M=m+D$ and define

$$R=(k-a^2)A^2+a(1-a)kAM-a^2k(k-1)mD.$$

The quartic contact residual from the geometry note factors as

$$
\frac{(a+1)^2(k-1)^2M^2}{a^2(k-a^2)}R,
\qquad
XY-kA^2=-\frac{k(k-1)}{(k-a)^2}R.\tag{4}
$$

Thus the quartic is exactly the rank-one boundary of this cone.
For $d>0,n\ge0$ and $m,C,A,D,\alpha>0$, recover

$$a=\frac{2A}{d+2n},\qquad k=1+\frac{C+d/2}{D}.$$

The three contact minors imply the other linear equalities in (2) and
$k>\max(1,a^2)$. If the quartic holds, the single additional branch
condition

$$aA-kn>0\tag{5}$$

is necessary and sufficient for an interior native generator. Equations
(3)--(4) give $XY=kA^2$, a positive sum $X+Y$, and
$w_-\le(X+Y)/A<w_+$. Hence (2) holds at rank one. The unique inverse is
$x=X/A$, $y=Y/A$, followed by (1)'s inverse formulas. In the original
five-parameter native chart, use
$b=a$, $t=(1-r)/a$, $z=s/(ar)$.

This branch condition is essential. At $a=2,k=11$, the positive source
$(1,38,77/4,33/4,77,7/8,192)$ obeys the angle law, minors and quartic,
but $aA-kn=-264$. Its algebraic inverse has $x=3/2<a$.

## The remaining convexification

Let $F_{a,k}$ denote (2). The closed cone generated by all native
maximally entangled C2 contacts is

$$\overline{\operatorname{cone}\!\bigcup_{a>0,\ k>\max(1,a^2)}F_{a,k}}.$$

The unresolved step is a fixed finite affine representation of this
continuum. Allowing $a,k$ to vary inside (2) introduces products with
source coefficients. The individual order-four pencils do not establish
a single finite lift of their combined cone.

The [exact checker](checks/c2_fixed_fibers.py) verifies the native
dictionary, domain identities, inverse reconstruction and both
factorizations in (4). The two-generator proof above supplies the
convex-hull argument.
