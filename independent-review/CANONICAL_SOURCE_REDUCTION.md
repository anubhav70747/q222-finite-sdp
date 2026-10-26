# Reducing canonical Bell sources to two closed faces and two contact families

For $d>0$ and nonnegative canonical event weights, a finite noncommuting stationary zero of the complete Bell state matrix determines a global quantum inequality.
This makes the extreme-source classification a finite support argument. Native entry then requires proving every denominator sign.

Write $H=1-(E_{00}-E_{01}+E_{10}+E_{11})/2$. In a chosen PR orientation, the eight losing events give the unique affine expansion
$$
f=dH+4\sum_i a_ip_i,\qquad
d=-f(\mathrm{PR}),\quad a_i=f(v_i)/4,
\qquad (a_i)=(\lambda,n,m,C,A,D,\alpha,\beta).
$$
Here $v_i$ is the deterministic facet neighbour on which $p_i=1$ and the other seven losing events vanish. Quantum validity therefore gives
$a_i\geq0$. A valid slack nonnegative at every PR vertex is a nonnegative sum of probability inequalities, since the nonsignalling vertices are the
deterministic and PR behaviours. Otherwise choose a PR orientation with $d>0$. The event definitions, full parent and determinant coefficients
are given in [SG1–SG6](../latex/chapters/04_geometry.tex).

**Source reduction.** Up to physical relabellings, every canonical valid source is a finite nonnegative sum of event sources, sources on the two
closed faces $\lambda=m=0$ and $\lambda=\alpha=0$, and contact sources in the following domains:
$$
\begin{array}{ll}
\mathrm{C1}:&d,m,C,A,\alpha>0,\quad \lambda=n=0,\quad
D,\beta\geq0,\quad D+\beta>0,\\
\mathrm{C2}:&d,m,C,A,D,\alpha>0,\quad\lambda=\beta=0,\quad n\geq0.
\end{array} \tag{1}
$$
Every extreme source outside the two closed-face orbits has one of these supports and a finite noncommuting contact. Its native coordinates have
strict original gaps, including at the allowed third omissions $D=0$ or $\beta=0$ in C1 and $n=0$ in C2. Native parameters conversely
give valid extreme rays when these gaps, the stationary formulas and all original source-coefficient signs are imposed.

## A stationary contact proves global validity

Let $F_c(u,v,\psi)=\psi^TP_c(u,v)\psi$ be the complete cleared parent SG2, with $\psi=(W,X,Y,Z)$ and $D_u=1+u^2$, $D_v=1+v^2$.
It equals $D_uD_v\|\psi\|^2 f(q)$. For $d>0$ and $uv\ne0$, the
$Z=0$ principal block is positive definite: the signed core already contains
$$
2d\{(uW+uvX-Y)^2+D_uv^2Y^2+D_vX^2\}.
$$
Consequently $P_c\succeq0$ exactly when $\det P_c\geq0$, and any zero
has a one-dimensional state kernel with $Z\ne0$. The determinant is
$$
\det P_c=16D_u^2D_v^2u^2v^2L(u^2,v^2),\qquad
L(U,V)=h+\sum_{e\in\{-1,0,1\}^2\setminus\{0\}}q_eU^{e_1}V^{e_2},
\quad q_e\geq0. \tag{2}
$$
Suppose $\psi_c\ne0$, $P_c\psi_c=0$ and
$\psi_c^TP_{c,u}\psi_c=\psi_c^TP_{c,v}\psi_c=0$ at positive
$u_c,v_c$. The adjugate is a positive multiple of $\psi_c\psi_c^T$,
so Jacobi's formula gives $L=L_U=L_V=0$ there. Set $U_c=u_c^2$, $V_c=v_c^2$,
and $w_e=q_eU_c^{e_1}V_c^{e_2}$. Then
$$
h=-\sum_e w_e,\qquad \sum_e w_e e=0,\qquad
L(U,V)=\sum_e w_e\{(U/U_c)^{e_1}(V/V_c)^{e_2}-1\}. \tag{3}
$$
Weighted AM–GM makes the last expression nonnegative. Algebraically,
the balanced vector decomposes into opposite pairs and triangles with
integer weights $(1,1)$, $(1,1,1)$ or $(1,1,2)$; each subtraction
eliminates an active exponent. The quadratic, cubic and quartic identities
in SG5–SG7 prove AM–GM for these finitely many cases over a real closed field.
Thus the full state matrix is PSD at every finite nonzero test.
At commuting or projective endpoints the CHSH core is nonnegative and
each event is nonnegative; homogeneous clearing includes these endpoints.
This proves global validity from the six stationary equations.

Every physical finite contact has coercive exponent support: zero lies
in the interior of the active exponent polygon. Otherwise balance in (3)
forces all active exponents onto a line through zero. Horizontal or vertical
opposite pairs in SG6 force an off-line positive event product. The diagonal
pairs leave two alternating four-event matchings, related by a physical
cycle rotation. For the matching $n,m,A,D$, the full parent at
$u=v=1$, $\psi=(1,-\epsilon,\epsilon,1)$ is
$16\epsilon\{(d+n+m+A+D)\epsilon-d\}<0$ for sufficiently small
$\epsilon>0$. If all outer coefficients vanish, the maximal remaining
event sets are rotations of $\{\lambda,m,A\}$; the test
$\psi=(1,-1,0,2)$ at $u=v=1$ gives $F_c=-4d$. Neither case is physical.

Coercivity also gives strict contact curvature:
$\sum_e w_e ee^T\succ0$. Completing the three state coordinates in the
single source field shows that the full Hessian in $(W,X,Y,u,v)$,
with $Z=1$, is positive definite. Deleting a positive event whose
amplitude vanishes preserves all six stationary equations and hence
validity. Such a deletion cannot create noncoercive support.

## Six linear equations restrict the extreme supports

In the event cycle $(\lambda,D,C,A,\alpha,m,\beta,n)$, noncoercivity
is equivalent to omitting a pair at cyclic distance three or four.
This follows directly from SG6: a separating direction may be chosen
among the four axes and four diagonals. For example, one axis requires
$\lambda\alpha=\lambda D=\alpha m=mD=0$, whose minimal omissions
are $\{\lambda,m\}$ and $\{\alpha,D\}$; the diagonal adds
$\{\lambda,\alpha\}$. Its other minimal omission sets already contain
an axial pair. Physical rotations and reflections give the twelve pairs.
Compatible omissions have distance one or two; at most three can occur,
and a compatible triple consists of three consecutive cycle vertices.

Normalize $d+\sum_i a_i=1$. Coercivity gives a finite minimum
$\mu\geq0$ of $L$. If $\mu>0$, small opposite nonradial perturbations
inside the active source face preserve validity. Indeed, positive
$q_e$ change by at most relative $\delta$, absent ones stay zero,
and $h$ changes by at most $\delta$, so uniformly
$L'\geq(1-\delta)\mu+\delta h-\delta>0$ for small $\delta$.
An extreme source therefore has $\mu=0$ and a finite stationary contact.
At that fixed contact the four state equations and two angle equations
are homogeneous linear equations in the source coefficients. Seven positive
events would give at least eight active coordinates including $d$, hence
a nonradial normalized direction preserving all six equations.
Small opposite perturbations would remain valid by (3).
Thus an extreme coercive source has five or six positive events.
The compatible omission patterns give exactly (1) and its stated third omissions.

The normalized source body is compact and eight-dimensional, so at most
nine extreme atoms suffice. Its finite quantified decomposition and the
contact implications transfer to real closed fields; this uses no
topological compactness assertion inside a non-Archimedean field.

## The two noncoercive faces include their singular boundaries

The twelve incompatible pairs form two physical orbits, represented by
$\lambda=m=0$ and $\lambda=\alpha=0$. Their exact source pencils are
[SG11–SG13](../latex/chapters/04_geometry.tex). A homogeneous physical limit
gives a quadratic matrix polynomial in one angle. Such a polynomial
$K_0+tK_1+t^2K_2$ is PSD for all real $t$ exactly when
$$
\begin{pmatrix}K_0&K_1/2+S\\K_1/2-S&K_2\end{pmatrix}\succeq0
\quad\text{for some }S^T=-S.
$$
For the first face this is the six-dimensional pencil SG12.
For the second, coefficient matching leaves the single scalar $\xi$
and two $3\times3$ blocks $B_1,B_2$ in SG13.
Their determinant identities recover the required part of $L$;
its remaining Laurent terms are nonnegative. This proves sufficiency,
while the homogeneous limits prove necessity, including singular blocks.
The SG13 transverse-positive-block argument requires $m+A>0$.
Feasibility guarantees it: if $m=A=0$, $B_1\succeq0$ forces $\xi=0$;
then $q=(1,0,1)^T$ has $q^TB_2q=0$ but
$B_2q=(0,-2d,0)^T\ne0$, impossible for a PSD matrix.

## Native entry follows from an analytic contact graph

Over $\mathbb R$, fix $C=2$ in either strict six-event support and vary the other five
positive event weights. Physical values of $d$ form a closed interval of
positive length: at $d=0$ the central determinant coefficient contains
$16mCAD>0$. A fixed CHSH-violating test bounds $d$ above.
At its maximum $d_*$, coercivity and strict logarithmic convexity give
one squared-angle contact and one state with $Z=1$.
Its core value $F_H$ is negative: otherwise every retained event amplitude
would vanish, forcing $P_H\psi=0$, whereas $\det P_H\ne0$.
The Jacobian in $z=(W,X,Y,u,v)$ and $d$ is
$\left(\begin{smallmatrix}\nabla_z^2F&\partial_d\nabla_zF\\0&F_H\end{smallmatrix}\right)$,
which is invertible. The analytic implicit-function theorem gives local
branches; stationarity makes each physical and maximal, and uniqueness
patches them into a connected analytic graph.

On this graph, the exact cofactors, Cramer identities, nonzero seeds and
full-parent reflections in [EN5–EN12](../latex/chapters/05_entry.tex)
fix every contact sign. A Cramer identity is first proved where its
selected minor is nonzero, then extended as a polynomial identity over
the graph: analyticity and a nonzero seed make that set dense.
Event-deletion contradictions exclude the possible equality walls.
This order avoids assuming the native denominator signs in order to prove them.
At a third omission, the same invertible Jacobian opens the missing
event coefficient at the fixed finite contact. The cofactor, deletion and
reflection identities exclude equality in the limiting inequalities;
strictness is not inferred from continuity alone.
The C2 boundary $n=0$ is already a C1 triple-omission orbit.

The resulting strict signs make every denominator in
[EN13–EN16](../latex/chapters/05_entry.tex) nonzero. The native rational
stationary formulas therefore give the unique source at that contact.
Conversely, those formulas with the original coefficient signs give
stationarity and hence global validity by (3). Their rank-six stationary
system makes any nonnegative valid decomposition proportional, proving
canonical extremality also at the allowed third omissions. The resulting finite first-order implications transfer to real closed fields.

This is a geometric source reduction. The next requirement is a certificate
for the **whole native parent** at every independently labelled source/test
pair as a sum of squares in $K_c\otimes_{\mathbb R}K_t$.
Nonnegativity of its Laurent determinant does not supply that certificate
or complete a finite SDP lift.
