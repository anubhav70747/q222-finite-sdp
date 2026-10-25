# An exact cone for two independent contact parameters

The equal-angle C2 contact family admits a direct representation by two
2-by-2 PSD matrices and three scalar inequalities. The construction covers
two independently varying parameters, includes the earlier quadratic curve,
and gives an explicit additional constraint for the NPA relaxation.

The main step is a decomposition: every feasible source is one tight contact
source plus nonnegative multiples of two event inequalities. This proves
both inclusions without constructing a large Gram matrix for the state
polynomial.

## The family and its cone

Use

$$H=1-\frac{E_{00}-E_{01}+E_{10}+E_{11}}2,\qquad
f_c=dH+4\sum_i c_i p_i,$$

with source order $(d,\lambda,n,m,C,A,D,\alpha,\beta)$ and event order

$$
(++|01),(+-|00),(+-|11),(-+|10),(-+|00),(-+|11),(--|01),(+-|10).
$$

Let $r,x>0$ satisfy

$$r+x+3rx\le1,$$

and put $S=r+x+2rx$. The sources are

$$
\begin{gathered}
d=1,\quad\lambda=\beta=0,\quad
n=\frac{1-r-x-3rx}{4S},\\
m=\frac{(1-r)(1+x)}{4S},\quad
A=\frac{(1+r)(1+x)}{4S},\quad
D=\frac{(1+r)(1-x)}{4S},\\
C=\frac{r^2(1+x)}{(1-r)S},\qquad
\alpha=\frac{x^2(1+r)}{(1-x)S}.
\end{gathered}\tag{1}
$$

The domain implies $r,x<1$, so every denominator is positive. In the native
C2 formulas these are $a=b_{\rm source}=1$, $t=1-r$, $z=x/r$.
Writing $s=\sqrt{x/r}$, their contact is

$$
A_0=B_0=\sigma_z,\qquad
A_1=B_1=\frac{(1-s^2)\sigma_z+2s\sigma_x}{1+s^2},\qquad
\psi\propto(1,-sr,sr,1)^T.
$$

The state is maximally entangled. The measurement angle and the rotation
of the state vary independently through $x,r$. Global validity of the
inequalities is proved below.

**Exact source cone.** The closed conic hull $K$ of (1) is precisely

$$
\begin{gathered}
\lambda=\beta=0,\qquad d=4A-2m-2D,\qquad n=m+D-A,\\
\begin{pmatrix}m&A-m\\A-m&C\end{pmatrix}\succeq0,\qquad
\begin{pmatrix}D&A-D\\A-D&\alpha\end{pmatrix}\succeq0,\\
A-m\ge0,\qquad A-D\ge0,\qquad m+D-A\ge0.
\end{gathered}\tag{2}
$$

All source coefficients, including $d$, are nonnegative under these
conditions. The total PSD order is seven. The only equalities are the
displayed linear relations; no contact parameter occurs in (2).

## Why the representation is exact

Substitution of (1) gives

$$mC=(A-m)^2,\qquad D\alpha=(A-D)^2,$$

as well as the linear relations and scalar inequalities in (2). Thus every
source ray, and hence its closed conic hull, obeys (2).

Conversely, first suppose $m,D>0$. Keep $A,m,D$ fixed and set

$$C_0=\frac{(A-m)^2}{m},\qquad
\alpha_0=\frac{(A-D)^2}{D}.$$

The PSD constraints say $C\ge C_0$ and $\alpha\ge\alpha_0$. The source
with $C=C_0$, $\alpha=\alpha_0$ is a scaled source from (1), or an endpoint
limit, with

$$r=\frac{A-m}{A+m},\qquad x=\frac{A-D}{A+D}.$$

Indeed, $m/A=(1-r)/(1+r)$ and $D/A=(1-x)/(1+x)$; the remaining scalar
inequality is exactly $r+x+3rx\le1$. These formulas also prove uniqueness
of the contact ray for fixed positive $A,m,D$.

The excesses $(C-C_0)e_C$ and $(\alpha-\alpha_0)e_\alpha$ are event rays.
Both belong to the closure of the contact family. For $e_C$, take
$r=1-\varepsilon$, $x=\varepsilon/8$ in (1) and divide the source by its
$C$ coefficient; the limit as $\varepsilon\downarrow0$ is $e_C$.
Interchanging $r,x$ gives $e_\alpha$.

The remaining endpoints are equally direct. Put $h=m/A$, $k=D/A$ for
$A>0$. Their domain is the closed triangle

$$0\le h,k\le1,\qquad h+k\ge1.$$

For $h,k>0$, the source with $A=1$ and minimal $C,\alpha$ is continuous:

$$
d=4-2h-2k,\quad n=h+k-1,\quad m=h,\quad D=k,\quad
C=(1-h)^2/h,\quad \alpha=(1-k)^2/k.
$$

Approaching the triangle boundary from $0<h,k<1$, $h+k>1$ proves membership
of every finite endpoint, including $h=k=1$, where $d=0$.
If $m=0$ or $D=0$, PSD and the scalar constraints in (2) force
$A=m=D=0$. Only the two nonnegative event axes remain. This proves the
reverse inclusion and all finite limits.

## The explicit constraint on behaviours

For a behaviour $q$, substitute the two linear source relations into its
slack. The coefficients of $(A,m,D,C,\alpha)$ become

$$
\begin{aligned}
L_A&=4(H-p_n+p_A),& L_m&=-2H+4(p_n+p_m),\\
L_D&=-2H+4(p_n+p_D),& L_C&=4p_C,& L_\alpha&=4p_\alpha.
\end{aligned}
$$

Every inequality in $K$ holds exactly when four real auxiliaries
$\xi,\eta,\rho,\sigma$ satisfy

$$
\begin{pmatrix}
L_A+L_m-2\eta-\sigma&\xi\\\xi&L_C
\end{pmatrix}\succeq0,\qquad
\begin{pmatrix}
L_A+L_D-2\xi-\rho&\eta\\\eta&L_\alpha
\end{pmatrix}\succeq0,\tag{3}
$$

$$\rho\ge0,\qquad\sigma\ge0,\qquad
2\xi+2\eta+\rho+\sigma-L_A\ge0.\tag{4}$$

To verify the forward direction, pair the two PSD matrices in (3) with
the source matrices in (2). Add $\rho(A-m)$, $\sigma(A-D)$ and
$(2\xi+2\eta+\rho+\sigma-L_A)(m+D-A)$. The sum is exactly $f_c(q)$.

For the converse, the source system is strictly feasible at
$A=1$, $m=D=3/4$, $C=\alpha=1$. Both 2-by-2 determinants equal $11/16$;
the scalar margins are $1/4,1/4,1/2$. Pairing with this fixed strict point
bounds the trace of any PSD dual preimage of a bounded coefficient vector.
The adjoint image is therefore closed, so finite-dimensional separation
gives an attained dual representation. Eliminating its fifth scalar
auxiliary gives (3)–(4). Thus these are exact finite affine constraints.

## Global quantum validity

Multiply the source (1) by $4S(1-r)(1-x)>0$. Its coefficients are

$$
\begin{gathered}
\widetilde d=4S(1-r)(1-x),\qquad
\widetilde n=(1-r-x-3rx)(1-r)(1-x),\\
\widetilde m=(1-r)^2(1-x^2),\quad
\widetilde C=4r^2(1-x^2),\quad
\widetilde A=(1-r^2)(1-x^2),\\
\widetilde D=(1-r^2)(1-x)^2,\qquad
\widetilde\alpha=4x^2(1-r^2).
\end{gathered}
$$

For arbitrary testing angles $u,v$, let $P(u,v)$ be the state matrix of
$(1+u^2)(1+v^2)\langle\psi|f_{\widetilde c}|\psi\rangle$.
Its $Z=0$ principal block is positive definite for $uv\ne0$, because the
signed core alone restricts to

$$2\widetilde d\bigl[(uW+uvX-Y)^2+(1+u^2)v^2Y^2+(1+v^2)X^2\bigr].$$

Set $z=x/r$, $X_0=u^2/z$, $Y_0=v^2/z$ and
$G=16(1-r^2)^3(1-x^2)^3$. The determinant is

$$\det P=16G(1+u^2)^2(1+v^2)^2u^2v^2 R(X_0,Y_0),$$

where

$$
\begin{aligned}
R={}&8rx(1-r)(X_0+X_0^{-1}-2)\\
&+8rx(1-x)(Y_0+Y_0^{-1}-2)\\
&+(1+r)(1+x)(1-r-x-3rx)(X_0/Y_0+Y_0/X_0-2)\\
&+8rx(r+x)(Y_0/X_0+X_0+Y_0^{-1}-3).
\end{aligned}
$$

Every coefficient and every parenthesis is nonnegative. For the last
parenthesis, the three factors have product one. Hence the determinant is
nonnegative. The positive three-state block then gives $P\succeq0$ by
Schur complementation. Continuity includes zero and infinite testing
angles after restoring the bounded observables. The two-reflection block
decomposition, real and imaginary state parts, and convex mixtures extend
validity to arbitrary quantum strategies. At the stated contact the slack
and its first derivatives vanish exactly.

## Relation to NPA and the remaining construction

Every source (1) with $r,x>0$ is an extreme ray of the canonical valid
nonnegative source cone. To see this, fix its contact and form the six linear
source equations from the four entries of $B_c\psi$, followed by
$\psi^T(\partial_uB_c)\psi$ and $\psi^T(\partial_vB_c)\psi$, where $B_c$
is the uncleared Bell operator and the unnormalised state is the one above.
For columns $(d,n,m,C,A,D,\alpha)$, the minor deleting $d$ is

$$-\frac{32768s^6r^4(r-1)(s^2r-1)(2s^2r+s^2+1)}{(1+s^2)^5}\ne0.$$

Any valid nonnegative decomposition preserves $\lambda=\beta=0$ and the
contact zero. Every summand therefore solves the same six equations and
is proportional to the original source. This also includes $n=0$.

Setting $r=1/2$, $x=z/2$ recovers the
[quadratic curve](../c2-quadratic-family/). Every positive-parameter member
of that curve is an extreme canonical ray outside `1+AB`. Therefore the
larger block (3)–(4) also strictly sharpens `1+AB`, while retaining every
quantum behaviour. It supersedes the curve block when these two contact
parameters are allowed to vary independently.

The exact representation here concerns the closed cone generated by (1).
General C1/C2 contacts have further independent parameters. Covering them
uniformly, together with their limits and the existing closed-face cones,
remains necessary for a full Q222 lift.

## Exact matrices and reproduction

[`source_pencil.json`](source_pencil.json) contains the order-seven matrix
in source coordinates $(A,m,D,C,\alpha)$ and the linear embedding into the
nine canonical coefficients. [`behavior_pencil.json`](behavior_pencil.json)
contains the order-seven matrix in the eight observed expectations and
four auxiliary variables. Every matrix entry is rational. The record
convention is specified in each file; off-diagonal entries are mirrored once.

With Python and SymPy, run:

```text
python verify_surface.py
python verify_dual.py
python build_pencils.py
```

The checks verify the full four-variable determinant, its equality to the
cleared Bell operator, the contact equations, source decomposition, endpoint
limits, both exported pencils and their exact trace pairing. No numerical
optimization is used.
