# A fixed polynomial multiplier cannot cover all native contacts

The obstruction from doubly tilted CHSH already occurs among the native
C1/C2 contacts. Tightening the tilted inequality forces every source in
a decomposition to share its noncommuting contact. This removes the two
noncoercive closed faces from the argument.

Write $F_c(u,v,\psi)$ for the complete polynomial parent in
[SG2](../latex/chapters/04_geometry.tex). A native family below includes
every finite contact in its prescribed two-omission chart; its strict
part has six positive event coefficients.

**Theorem.**

1. No fixed nonzero polynomials $p(u),q(v)$, each nonnegative on
   $\mathbb R$, make $p(u)q(v)F_c$ a polynomial sum of squares for every
   strict native C1 and C2 source.
2. No fixed real polynomial $P(u,v)$ with $P(0,0)>0$ makes $PF_c$ a
   polynomial sum of squares for every strict native C1/C2 source and
   all its physical relabelings in the original testing frame.

The roots may depend on the source. No degree bound is assumed: for a
fixed multiplier, the degree of the parent bounds the degree of every
root. The second statement allows nonseparable multipliers, and does
not assume that the multiplier itself is a sum of squares.

## Tighten the tilted source

Put $D_u=1+u^2$, $D_v=1+v^2$, use the real qubit observables of SG2,
and define

$$
B_T=A_0B_0-A_0B_1+A_1B_0+A_1B_1+(1-T)(A_1+B_0).
$$

Its local maximum is $\ell_T=4-2T$. Let $b_T$ be its quantum maximum.
The [fixed-multiplier theorem](../extracted-lifts/canonical-completion/polynomial-obstruction/source/fixed_polynomial_multipliers.tex)
proves that $w_T=4-2T+4T^3$ is a strict upper bound for
$0<T\leq1/100$. Its proof is an explicit determinant and inertia
calculation for both independent measurement angles.

The same theorem supplies a polynomial state test. Set $T=t^2$ and

$$
\psi_t(x,y)=\left(
1+\frac23t^4x^2y^2,
-\frac13t^3x^2y,
-tx+\frac23t^3xy^2,
-t^2xy-\frac23t^4x^3y
\right)^T.
$$

For $F_T=D_uD_v\psi^T(w_TI-B_T)\psi$ it gives

$$
F_{t^2}(-tx,ty;\psi_t)
=\frac43t^6\bigl(M(x,y)+2\bigr)+O(t^8),
\qquad
M=x^4y^2+x^2y^4-3x^2y^2+1.
\tag{1}
$$

The remainder is polynomial and its degree in $(x,y)$ is at most
twelve. At $x=y=1$, the corresponding normalized Bell expectation is
$\ell_T+\frac43T^3+O(T^4)$. Hence, for small positive $T$,

$$
\ell_T<b_T<w_T<4.
\tag{2}
$$

The compact real two-qubit model attains $b_T$. Neither measurement
pair can commute at a maximizing test, since a commuting pair gives a
local behavior. Both half-angles are therefore finite and nonzero;
state sign congruences make them positive.

The tight parent $F_T^*=D_uD_v\psi^T(b_TI-B_T)\psi$ has coefficients

$$
\begin{aligned}
d&=4-b_T,&\lambda=m=A&=(b_T-\ell_T)/4,\\
n=D=\alpha&=(b_T-2T)/4,&C=\beta&=(b_T-2)/4.
\end{aligned}
$$

They are all positive by (2). Thus $F_T^*$ is a canonical source with
a finite noncommuting contact.

## Decompose while keeping that contact

Fix its contact $(u_c,v_c,\psi_c)$. In the nine source coefficients,
the equations

$$
P_c\psi_c=0,\qquad
\psi_c^T(\partial_uP_c)\psi_c=0,\qquad
\psi_c^T(\partial_vP_c)\psi_c=0
\tag{3}
$$

are six homogeneous linear equations. Intersect them with $c\geq0$
and $\sum_i c_i=1$. The result is a compact polytope, and each vertex
has at most seven positive coordinates: otherwise the seven equations,
including normalization, admit a nonzero perturbation on its support.

Every vertex with $d>0$ is globally valid by SG7. SG17 shows that a
physical source with $d>0$ and this finite noncommuting contact has
coercive Laurent support. SG9 then requires at least five positive
event weights; the vertex bound allows at most six. SG10 identifies
exactly C1, C2 and their permitted third omissions, up to the finite
physical symmetries. A vertex with $d=0$ is an event combination.
Consequently

$$
F_T^*=N_T+E_T,
\tag{4}
$$

where $N_T$ is a nonnegative sum of native-orbit parents and $E_T$ is
an event sum of polynomial squares. Each summand is nonnegative at
every independent real test.

The finite third-omission contacts follow from strict ones. The
[opening argument](EXACT_LIFT_ROUTE.md#only-the-strict-contact-families-need-new-admission-proofs)
keeps the two prescribed omissions and solves the five contact
stationarity equations and zero value for the five contact variables
and $d$. Its Jacobian has positive-definite contact Hessian and
nonzero core derivative. The implicit-function theorem therefore opens
the third weight. For any fixed multiplier, coefficient convergence
and closedness of the fixed-degree SOS cone preserve its certificate.

This decomposition uses SG7, SG9, SG10 and SG17 at a finite contact.
It requires no polynomial certificate for the full parents of SG12 or
SG13; their finite source pencils alone would not establish one.

## The limiting contradiction

First prove statement 2. Define
$\kappa(T)=(b_T-\ell_T)/T^3$. Equation (2) and the state test give
$\liminf\kappa\geq4/3$ and $\kappa<4$. Along a sequence $T=t^2\downarrow0$,
choose $\kappa(T)\to\kappa_0\in[4/3,4]$. Subtracting the strict gap
from (1) yields

$$
t^{-6}F_T^*(-tx,ty;\psi_t)
\longrightarrow
\frac43\left(M+\frac{3\kappa_0}{4}-1\right).
\tag{5}
$$

The added constant belongs to $[0,2]$. For every $a\geq0$,
$M+a$ is not polynomial SOS: a root has degree at most three and
must be constant on either coordinate axis, where the sum of squares
is $1+a$. Each root thus has the form
$r+sxy+hx^2y+kxy^2$. Its square sum has a nonnegative $x^2y^2$
coefficient, whereas that coefficient in $M+a$ is $-3$.

Apply the same substitution and division to (4), writing
$f_t=n_t+e_t$. These polynomials have degree at most twelve;
$n_t,e_t\geq0$ and $e_t$ is SOS. On $[-1,1]^2$,

$$
0\leq\int e_t\leq\int f_t.
$$

The right side stays bounded by (5). If $z$ lists monomials through
degree six, the matrix $\int zz^T$ is positive definite. Hence this
integral bound bounds the trace of a positive semidefinite Gram for
$e_t$. A subsequence converges to a polynomial SOS $e$, and
$n_t=f_t-e_t$ converges coefficientwise to $n$.

The assumed common multiplier makes $P(-tx,ty)n_t$ SOS. Its degree
is bounded by $12+\deg P$ and its limit is $P(0,0)n$. The same
integral-Gram argument proves closedness of this SOS cone, so $n$ is
SOS. Then the limit $f=n+e$ in (5) is SOS, contradicting the preceding
coefficient argument. The event terms were bounded before applying
the multiplier; no SOS assumption on $P$ was used.

For statement 1, transport the assumed separable multiplier to the
finite physical orbits. The exact transport is given below. It yields
one fixed separable nonnegative multiplier $M_*$ that certifies every
term in (4), including events. Thus $M_*F_T^*$ is SOS. Adding back

$$
(w_T-b_T)D_uD_v\|\psi\|^2
$$

preserves SOS after multiplication by $M_*$. This would certify
$M_*F_T$ for every small $T$, contrary to the fixed-multiplier theorem.

### Finite relabeling transport

Let $N(w)=\left(\begin{smallmatrix}1&w\\-w&1\end{smallmatrix}\right)$,
$X=\left(\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\right)$, and let
$\Pi$ exchange the two state factors. Three physical symmetries give
the exact parent identities

$$
\begin{aligned}
D_uF_{g_Ac}(u,v;\psi)
 &=v^2F_c(-u,-1/v;(N(u)\otimes I)\psi),\\
D_vF_{g_Bc}(u,v;\psi)
 &=u^2F_c(1/u,-v;(X\otimes N(v))\psi),\\
F_{g_Tc}(u,v;\psi)
 &=u^2v^2F_c(-1/v,-1/u;\Pi\psi).
\end{aligned}
\tag{6}
$$

In the zero-based event cycle $(\lambda,D,C,A,\alpha,m,\beta,n)$,
these act as $i\mapsto5-i$, $3-i$ and $2-i$ modulo eight. They
generate the sixteen dihedral event symmetries. Equations (6) are
polynomial identities after clearing; axes are included by continuity.

Starting from an SOS of $p(u)q(v)F_c$, substitute one identity and
clear reciprocal powers by sufficiently large even powers of $u,v$.
Multiply by the displayed $D_u$ or $D_v$ where required. This gives
a fixed separable nonnegative multiplier $M_g$ for the transformed
family. The degree bound on the original roots fixes the clearing
power independently of the source. Iterate over the finite group.
The product $M_*=\prod_g M_g$ works for every orbit, because each
quotient $M_*/M_g$ is a product of nonnegative univariate polynomials
and hence SOS. The same reasoning covers event parents in SG2.

## Consequence for rational testing functions

A fixed finite rational angular span, linear in the state, cannot
certify all native orbits if its features have a common scalar
denominator $Q(u,v)$ with $Q(0,0)\ne0$. Clearing the denominator of
a source-only positive semidefinite Gram would give a polynomial SOS
of $Q^2F_c$, contradicting statement 2. A finite collection of such
spans is covered by taking the product of their denominators.
Statement 1 also excludes common separable denominators for the two
representative families, even when those denominators vanish at zero.

The argument does not settle a nonseparable multiplier vanishing at
the origin, or a nonseparable multiplier supplied only on the two
representative charts. Fractional-power testing spaces remain eligible.
Their uniform admission problem is still open.

The [exact checker](checks/native_multiplier_obstruction.py) verifies
all three identities (6), their sixteen-element group, and the
constants in (1), (2) and (5). The contact decomposition and limiting
argument above supply the quantifiers that algebraic checks alone do
not establish.
