# A finite testing cover for compact strict native families

A compact set of strict native contact sources admits one finite, globally sound SDP certificate cone. The construction uses a fixed testing grid; the source contact changes Gram coefficients, not testing functions.

**Theorem.** Normalize $d=1$. Let $K$ be a compact subset of either strict six-event native domain, or of their union:

$$
\mathsf C_1:\ \lambda=n=0,\quad m,C,A,D,\alpha,\beta>0;
\qquad
\mathsf C_2:\ \lambda=\beta=0,\quad n,m,C,A,D,\alpha>0.
$$

Require the finite native contact parameters and all native denominators and positive gaps to stay in their strict domain, bounded away from zero on $K$. There is a finite integer $N_K$ and a closed semidefinite-representable cone $\mathcal C_K$ such that

$$\operatorname{cone}(K)\subseteq\mathcal C_K\subseteq\mathcal C_{\mathrm{can}},$$

where $\mathcal C_{\mathrm{can}}$ is the cone of nonnegative canonical sources whose complete Bell slack is quantum-valid. The representation uses at most $4N_K^2$ PSD blocks of order 36. Its size depends on $K$; no uniform bound as $K$ approaches a singular source boundary is asserted.

## The zeros and their local completions

Write $F_c(u,v;\psi)=\psi^TP_c(u,v)\psi$ for the complete parent [SG2](../latex/chapters/04_geometry.tex), with $\psi=(W,X,Y,Z)^T$. Use the four testing charts with coordinates $(x,y)$ equal to $(u,v),(1/u,v),(u,1/v),(1/u,1/v)$ on $[0,1]^2$. For each reciprocal coordinate, multiply by its square: for example, $\widehat F_c(x,y;\psi)=x^2F_c(1/x,y;\psi)$. Each chart matrix $\widehat P_c$ is again polynomial of bidegree at most $(2,2)$.

Include **all four finite coordinate images of the moving contact** in the bounds below, even when an image lies outside a chart's unit square. Compactness of the strict native parameter set keeps these images bounded and away from zero. Thus a contact crossing $u=1$ or $v=1$ is not lost from a neighboring chart's estimates.

The positive projective testing square has one finite noncommuting zero, and the following deterministic zeros:

| Family | Corner | Kernel state |
|---|---|---|
| Both | $(u,v)=(0,0)$ | $W$ |
| C1 | $(u,v)=(\infty,0)$ | $X$ |
| C2 | $(u,v)=(\infty,\infty)$ | $Z$ |

Indeed, [SG6](../latex/chapters/04_geometry.tex) gives $q_{01},q_{10}>0$. The balanced Laurent determinant is strictly convex in logarithmic angles, so its finite positive zero is unique. On compatible edges, SG2 has no common kernel of the core and retained event rows except at the listed corners. The corner scalar forms, with the kernel coordinate set to one, are

$$
\begin{aligned}
g_{00}(u,v)&=(4C+4D+2)u^2+4mv^2+4Cu^2v^2,\\
g_{\infty0}(p,v)&=(4m+4\beta+2)p^2+(4D+2)v^2+(4\beta+2)p^2v^2,\\
g_{\infty\infty}(p,q)&=4Dp^2+(4m+4\alpha+2)q^2+4\alpha p^2q^2.
\end{aligned}
$$

Their scalar Grams on $(s,t,st)$ are positive definite. The complementary three-state blocks are positive definite by SG2. The leading chart determinants are respectively

$$
16(q_{10}u^2+q_{01}v^2),\qquad
16(q_{10}p^2+q_{21}v^2),\qquad
16(q_{12}p^2+q_{21}q^2),
$$

up to terms of total degree at least four. Each displayed coefficient is positive in its applicable strict family. Dividing by the positive complementary determinant shows that the state-normal/angular Hessian is positive definite. At the moving contact this follows from [SG18–SG19](../latex/chapters/04_geometry.tex); its scalar Gram is positive definite by [the contact-scalar result](CONTACT_SCALAR_GRAM.md). Reciprocal clearing acts invertibly on the bilinear scalar-root space, preserving this assertion in every chart.

At any of these zeros use source-only state coordinates $(\eta,z)$, with $\eta\in\mathbb R^3$ and kernel coordinate $z$. Put $s=x-x_c$, $t=y-y_c$ and $q=(s,t,st)^T$. Then

$$
\widehat F=\eta^TH\eta+2z\eta^T\zeta+z^2q^TGq,\qquad
G=\begin{pmatrix}J&k\\k^T&g\end{pmatrix}\succ0.
$$

Write $\zeta=\sum_{i+j>0}\zeta_{ij}s^it^j$, let $A_0=(\zeta_{10},\zeta_{01})$, and choose $w=A_0J^{-1}k$. The following bilinear columns give $Vq=\zeta$ exactly:

$$
\begin{aligned}
V_1&=\zeta_{10}+\zeta_{20}s+(\zeta_{11}-w)t,\\
V_2&=\zeta_{01}+\zeta_{02}t,\\
V_3&=w+\zeta_{21}s+\zeta_{12}t+\zeta_{22}st.
\end{aligned}
$$

Consequently

$$
\widehat F=(V^T\eta+Gqz)^TG^{-1}(V^T\eta+Gqz)+\eta^TR\eta,
\qquad R=H-VG^{-1}V^T,
$$

and

$$R(0,0)=H(0,0)-A_0J^{-1}A_0^T\succ0.$$

The last inequality is the Schur complement of the strict contact Hessian. Every inverse is source-only, and $R$ has bidegree at most $(2,2)$.

## One fixed grid

Over $K$ and the finite zero-chart collection, choose

$$
\Lambda=\min\lambda_{\min}(R_{00})>0,\qquad
M=\max\sum_{i+j>0}\|R_{ij}\|_F<\infty,\qquad
\rho=\min\left\{\frac18,\frac{\Lambda}{4(M+\Lambda)}\right\}.
$$

The coefficients here are Taylor coefficients about the corresponding zero. On any cell contained in $|s|,|t|\le2\rho$, every degree-$(2,2)$ Bernstein coefficient of $s^it^j$ has magnitude at most $(2\rho)^{i+j}$. Hence

$$\|B_{ij}(R)-R_{00}\|\le2\rho M\le\Lambda/2,$$

so all these matrix coefficients are positive definite.

On the complement of the radius-$\rho$ zero squares in each chart, let $\kappa>0$ be the uniform minimum eigenvalue of $\widehat P_c$. This exists because that compact set contains no zero. Let $M_1$ bound the sum of nonconstant Taylor coefficient norms of $\widehat P_c$ about all chart points. Choose

$$
N_K\ge\left\lceil\max\left\{\frac2\rho,\frac{4(M_1+\kappa)}{\kappa}\right\}\right\rceil.
$$

A grid cell with midpoint within $\rho$ of a zero lies in its radius-$2\rho$ square and uses the completed form above. Every other midpoint has $\widehat P_c\succeq\kappa I$; its Bernstein coefficient error is at most $M_1/N_K\le\kappa/2$. This covers every cell, including the transition between the two cases.

For a fixed interval $I=[L,R]$, put $\theta=(x-L)/(R-L)$ and

$$b_0^I=1-\theta,\qquad b_1^I=\sqrt{2\theta(1-\theta)},\qquad b_2^I=\theta.$$

Their squares are the quadratic Bernstein basis. On $I\times J$, use the fixed column

$$\Phi_{I,J}=\bigl(b_i^I(x)b_j^J(y)\psi_\ell\bigr)_{0\le i,j\le2,\ 1\le\ell\le4}\in\mathbb R^{36}.$$

Positive Bernstein coefficient matrices give squares in this span. In a local completion, the other square roots are bilinear polynomials times source-only linear state forms; these use only $b_0,b_2$ in each coordinate and lie in the same span. Thus every $c\in K$ satisfies, on every chart cell,

$$\widehat F_c=\Phi_{I,J}^TQ_{I,J}\Phi_{I,J},\qquad Q_{I,J}\succeq0.$$

The branch points are fixed grid endpoints. Source-dependent contact shifts occur only in Gram coefficients.

## The certificate cone

Define $\mathcal C_K$ by $c\ge0$ and all these cell identities simultaneously. Products of the fixed algebraic features span a finite-dimensional function space, so coefficient matching is linear in $c$ and the Gram entries. This is one finite SDP representation. The construction contains $K$, and homogeneity gives $\operatorname{cone}(K)\subseteq\mathcal C_K$.

Every represented parent is nonnegative on all four projective charts. The exact [angle-sign congruences](CONTACT_SCALAR_GRAM.md) give the other quadrants, and continuity includes endpoints. The two-observable qubit reduction then gives $\mathcal C_K\subseteq\mathcal C_{\mathrm{can}}$.

Each cell cone is closed. Integrate $\Phi\Phi^T$ over the cell interior and the state sphere. Linear independence of the 36 features makes this moment matrix positive definite. For any coefficient-convergent sequence of represented parents, its integrals bound the traces of PSD Grams, so a Gram subsequence converges and retains the identity. A finite intersection of these cell cones is closed.

The constants $\Lambda,M,\kappa,M_1$, hence $N_K$, depend on the chosen compact strict set. Uniform control of these quantities under singular native limits remains a separate obligation.
