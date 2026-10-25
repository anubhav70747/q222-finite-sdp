# The unequal-angle C2 contact equations

Allowing the two measurement angles to differ adds one parameter to the
maximally entangled C2 contacts. Their coefficients still obey a small set
of identities. Three quadratic equations determine ratios of source
coefficients; a fourth-degree equation fixes a direction that all three
quadratics miss. This distinction matters when constructing a lift: the
identities describe individual contact sources, whereas an SDP must describe
their convex cone.

## A symmetric parametrization

Use the canonical slack and event order from the
[equal-angle construction](../extracted-lifts/c2-equal-angle-cone/).
Here $\lambda=\beta=0$, and the active order is
$(d,n,m,C,A,D,\alpha)$. Let

$$a>0,\qquad 0<r,s<1,\qquad g=a+r+s-1>0.$$

Write $e=r+s$, $p=rs$ and define

$$
\begin{aligned}
L&=a(e-p)+p,& Q&=ag+L,\\
\Delta&=a^3p+a^2(e^2-p)+ap(2e-p)+p^2,\\
B_0&=a^3p+a^2(e^2+p-e)+a((3e-2)p-p^2)+p^2.
\end{aligned}
$$

The remaining domain condition is $B_0\le0$. The normalized sources are

$$
\begin{aligned}
d&=1,&
n&=-\frac{aB_0}{2(a+1)\Delta},\\
m&=\frac{(a+s)(1-r)(as+r)L}{2(a+1)g\Delta},&
C&=\frac{r^2(a+s)Q}{2(1-r)\Delta},\\
A&=\frac{a(a+r)(a+s)L}{2(a+1)\Delta},&
D&=\frac{(a+r)(1-s)(ar+s)L}{2(a+1)g\Delta},\\
\alpha&=\frac{s^2(a+r)Q}{2(1-s)\Delta}.
\end{aligned}\tag{1}
$$

Every denominator is positive. In particular, $e-p>0$,
$e^2-p=r^2+rs+s^2>0$ and $2e-p>0$, so $L,Q,\Delta>0$.
The domain also implies $e<1$: if $e\ge1$, every displayed
group in $B_0$ is positive. The face $n=0$ is included by $B_0=0$;
none of these formulas divides by $n$.

In the native C2 parametrization, set
$b_{\rm source}=a$, $t_{\rm source}=(1-r)/a$ and $z=s/(ar)$.
The contact has
$u_c=\sqrt z$, $v_c=a\sqrt z$ and
$\psi_c=(1,-r\sqrt z,r\sqrt z,1)^T$.
The state is maximally entangled. All three parameters in (1) vary
independently over the stated domain.

These are globally valid native sources, not merely positive operators
at the contact. The argument in SG3–SG7 of
[the source geometry](../latex/chapters/04_geometry.tex) applies to the
complete parent: its three-state principal block is positive for nonzero
testing angles, while its determinant is a positive testing factor times
a Laurent polynomial with nonnegative outer coefficients. The contact
equations give a stationary zero of that Laurent polynomial. Convexity
in the logarithms of the squared half-angle coordinates makes this zero
a global minimum. Schur complementation then gives the full
four-state inequality; continuity includes compatible measurements.

## The coefficient ratios

A direct consequence of (1) is

$$a(d+2n)=2A.\tag{2}$$

Thus the angle ratio is recovered from a source ray. Introduce

$$
B=A+n,\quad U=C+d/2,\quad V=\alpha+d/2,\quad
S=U+V,\quad M=m+D,\quad \Gamma=a^2(B+S)-B.
$$

All $B,U,V,S,M$ are positive on the domain of (1). The three relations are

$$
mU=DV,\qquad m\Gamma=BV,\qquad D\Gamma=BU.\tag{3}
$$

Equivalently, every $2$-by-$2$ minor of
$\begin{pmatrix}m&D&B\\V&U&\Gamma\end{pmatrix}$ vanishes.
Summing the last two equations gives

$$a^2M(B+S)=B(M+S).\tag{4}$$

These equations alone cannot determine the contact surface. For fixed $a$,
adding a multiple of

$$
w=\left(1,-\frac{a}{2(a+1)},0,-\frac12,
             \frac{a}{2(a+1)},0,-\frac12\right)\tag{5}
$$

in active source order preserves (2), $B,U,V,M$, and hence all of (3).
It changes $d$.

The missing relation along this direction is

$$
S^2(B+d/2)(M+d/2)
       =(a+1)^2UVM(B+S).\tag{6}
$$

Substitution of (1) proves (6) throughout the domain, including $n=0$.
For fixed $B,U,V,M,a$, its left-hand side is strictly increasing in
$d\ge0$. Thus (6) removes the ambiguity in (5).

## Why a fourth-degree relation appears

At a maximally entangled state the local marginals vanish. The slack becomes

$$
T-\bigl(P E_{00}-V E_{01}+U E_{10}+R E_{11}\bigr),
\qquad
T=B+M+S,\quad P=B+d/2,\quad R=M+d/2.
$$

For Bob's two unit vectors with inner product $\tau$, optimizing Alice's
vectors gives

$$
\sqrt{P^2+V^2-2PV\tau}
+\sqrt{U^2+R^2+2UR\tau}.
$$

At an interior optimum, differentiation fixes the ratio of the two square
roots. Eliminating $\tau$ yields the squared maximum

$$
\frac{(PV+UR)(PU+VR)(PR+UV)}{PRUV}.
$$

The native contact is interior and saturates the valid inequality.
Equating this expression to $T^2$ explains the fourth-degree equation.
To see its two algebraic branches, put

$$
k=B/S,\quad t=M/S,\quad y=d/(2S),\quad
\theta=UV/S^2,\quad x=(k+y)(t+y).
$$

The equality reduces to

$$
\bigl[x-(2kt+k+t)\theta\bigr]^2
       =4kt(k+1)(t+1)\theta^2.\tag{7}
$$

Equation (4) says $a^2t(k+1)=k(t+1)$. Since $a>0$, the roots in (7) are

$$x_\pm=(a\pm1)^2t(k+1)\theta.$$

The native sources take $x_+$, exactly (6); direct substitution verifies
this branch. Squaring the support equation would leave the other branch
unless it were excluded separately. An intrinsic converse to (2), (3)
and (6), with the necessary domain and branch conditions, is still needed
before these equations can replace the parametrization.

## Two exact tests of proposed simplifications

At $a=2,r=1/2,z=1/16$, the native source is

$$
c=\left(1,\frac{43}{299},\frac{77}{1495},\frac{891}{598},
        \frac{385}{299},\frac{119}{598},\frac9{598}\right).
$$

Every active coefficient of $c+w/1000$ is positive. Equations (2)–(3)
remain true, but its normalized slack at
$u=1/4$, $v=1/2$, $\psi=(1,-1/8,1/8,1)^T$ is

$$-\frac{963}{1105000}.$$

The left side minus the right side of (6) is instead

$$\frac{899901143099}{106923596000000}>0.$$

Thus source signs and the quadratic equations admit an invalid inequality;
the fourth-degree equation detects it.

There is a separate obstruction to applying the generator equations
directly to mixtures. Let
$R_c=\begin{pmatrix}m&D\\ \alpha+d/2&C+d/2\end{pmatrix}$.
Each generator has rank one, with second row
$\kappa$ times the first, where
$\kappa=(a+1)g/[(1-r)(1-s)]$. For two generators,

$$
\det(R_1+R_2)=(\kappa_2-\kappa_1)(m_1D_2-D_1m_2).
$$

Choose $a_1=1$, $a_2=2$ and assign
$(r,s)=(1/8,1/16)$ and $(1/16,1/8)$, respectively.
The determinant is $-16557677/50331190$. Reversing the ordered pairs
gives its positive counterpart. Both sums are valid.
Neither determinant sign therefore describes the conic hull.

The [order-six curve cone](../extracted-lifts/c2-unequal-angle-curve/)
gives an exact finite representation when $r=1/4$, $s=1/8$ and
$3/4\le a\le2$. The full three-parameter convex hull and the remaining
general C1/C2 contacts still require uniform certificates.

Run
[the exact identity check](checks/c2_contact_geometry.py)
with Python and SymPy to reproduce (2)–(7), both signed mixtures, and the
invalid-source slack and quartic residual. These checks concern the
displayed identities and witnesses; global validity uses the native-source
argument above.
