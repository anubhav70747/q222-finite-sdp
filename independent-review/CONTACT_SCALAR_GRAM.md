# Three squares for the contact scalar

Fix a globally valid source and its state at a finite noncommuting contact. Evaluating the Bell matrix on that state, while allowing both measurement angles to vary independently, gives a scalar polynomial with a simple exact certificate. It needs at most three squares from the fixed testing span $\{1,u,v,uv\}$. The remaining difficulty is completing its three couplings to the other state coordinates.

## A biquadratic polynomial with a zero

Let $g(u,v)\ge0$ for all real $u,v$, with degree at most two in each variable and $g(u_c,v_c)=0$. Put $s=u-u_c$ and $t=v-v_c$. Its gradient vanishes at the zero, so uniquely

$$
g=A s^2+B t^2+2Cst+2D s^2t+2E st^2+F s^2t^2
=q^TGq,\qquad
q=(s,t,st)^T,\quad
G=\begin{pmatrix}A&C&D\\C&B&E\\D&E&F\end{pmatrix}.
$$

**The matrix $G$ is positive semidefinite.** For any $(x,y,z)$ with $xyz\ne0$, take $s=z/y$, $t=z/x$ and $\lambda=z/(xy)$. Then $q=\lambda(x,y,z)^T$, and nonnegativity of $g$ gives nonnegativity of the quadratic form on this direction. These directions are dense, which proves the claim. Conversely, $G\succeq0$ immediately gives $g\ge0$.

All entries are explicit derivatives at the contact:

$$
G=\begin{pmatrix}
g_{uu}/2&g_{uv}/2&g_{uuv}/4\\
g_{uv}/2&g_{vv}/2&g_{uvv}/4\\
g_{uuv}/4&g_{uvv}/4&g_{uuvv}/4
\end{pmatrix}_{(u_c,v_c)}.
$$

The shifts only change source coefficients. Indeed, for $\phi=(1,u,v,uv)^T$,

$$
q=D_c\phi,\qquad
D_c=\begin{pmatrix}-u_c&1&0&0\\-v_c&0&1&0\\u_cv_c&-v_c&-u_c&1\end{pmatrix},
\qquad g=\phi^T(D_c^TGD_c)\phi.
$$

Factoring $G$ therefore gives at most three squares in one fixed testing space. Singular $G$ is allowed.

## Applying the lemma to the Bell matrix

Write $F_c=\psi^TP_c(u,v)\psi$ for the complete canonical parent with $d>0$ in [SG2](../latex/chapters/04_geometry.tex). Let $\psi_c=\operatorname{vec}(\Psi_c)$ be its contact state in row-major order, at finite $u_c,v_c>0$. The entries of $P_c$ have degree at most two in each angle. Positivity extends to the whole real testing plane by

$$
P_c(-u,v)=S_uP_c(u,v)S_u,\quad S_u=\operatorname{diag}(1,1,-1,-1),
$$
$$
P_c(u,-v)=S_vP_c(u,v)S_v,\quad S_v=\operatorname{diag}(1,-1,1,-1),
$$

and by continuity on the axes. Hence the lemma already applies to $\psi_c^TP_c(u,v)\psi_c$.

To locate this scalar within the full matrix, assume $\delta=\det\Psi_c\ne0$, set

$$
\epsilon=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
L=\epsilon\Psi_c^T,\quad T_c=I\otimes L^{-1},
\quad E=(e_{00},e_{01}+e_{10},e_{11}),\quad h=(0,-1,1,0)^T.
$$

Use the constant basis $B=(E,h)$ to write

$$
T_cP_cT_c^T=B\begin{pmatrix}H&\zeta\\\zeta^T&\gamma\end{pmatrix}B^T.
$$

The identities $\Psi_cL=\delta\epsilon$ and $\operatorname{tr}(P_c(\epsilon\otimes\epsilon))=0$ give

$$
\gamma=\frac{\psi_c^TP_c\psi_c}{4\delta^2},\qquad H_{13}=H_{22}-\gamma.
$$

For the second identity, $E^T(\epsilon\otimes\epsilon)E$ has entries $J_{13}=J_{31}=1$, $J_{22}=-2$, while $h^T(\epsilon\otimes\epsilon)h=2$ and the cross terms vanish. The local congruence preserves the zero trace condition. At the contact, $\gamma=\zeta=0$ and $H\succ0$ is a three-by-three Hankel matrix. The source-only congruence retains all independent testing angles away from that contact.

The scalar lemma gives $\gamma=q^TG_\gamma(c)q$ with $G_\gamma(c)\succeq0$. No multiplier is needed. For native C2, $\delta=(b+c^2r(b-a^2t))/a>0$; for C1, $\delta=S(1+k^2r(1-\theta))>0$. The construction includes the permitted zero-event faces and never divides by their vanishing event weights.

## The native scalar Gram is positive definite

At every finite native C1 or C2 contact, $G_\gamma\succ0$. This includes the permitted $D=0$ or $\beta=0$ faces of C1 and $n=0$ face of C2. The same rank-three statement holds on the maximally entangled subfamily.

In both families, $m,A,\alpha>0$ and all four contact-state coordinates are nonzero. The [C1 state](../latex/certificates/c1_short.tex) is $(S,k(\theta-1),kSr,1)$ with $k,S>0$ and $0<r,\theta<1$; the [C2 state](../latex/certificates/c2_short.tex) is $(b/a,-cr,c(b-a^2t)/a,1)$ with $a,b,c,r>0$ and $b-a^2t>0$. The determinant coefficients [SG5–SG6](../latex/chapters/04_geometry.tex) satisfy

$$
q_{01}=4\alpha m t_5t_8>0,\qquad
q_{10}=4\alpha A t_7t_3>0.
$$

Their Laurent exponents $(-1,0)$ and $(0,-1)$ are independent. Thus the Hessian of $L(e^x,e^y)$, namely $\sum_e q_e e^{e\cdot(x,y)}ee^T$, is positive definite everywhere. Its stationary zero is unique. Consequently, at nonzero finite angles the determinant of $P_c$ vanishes only at $(\pm u_c,\pm v_c)$. The state kernel there is one-dimensional by SG3. The sign congruences above and the four nonzero contact amplitudes imply that the fixed state $\psi_c$ belongs to that kernel only at $(u_c,v_c)$.

Compatible measurements also exclude zeros of this fixed-state scalar. If either local pair commutes, the CHSH core is positive semidefinite, as are all event terms. The positive $A$-event gives

$$
\frac{\gamma(u,v)}{D_uD_v}\ge\frac{A Y_c^2}{\delta^2}>0
$$

on every compatible boundary, including $u=0$, $v=0$ and the projective endpoints $u=\infty$ or $v=\infty$. Hence $\gamma$ has exactly one finite zero and is strictly positive at every projective endpoint after the indicated normalization.

Its angular Hessian at the zero is positive definite too. Let $H_3$ be the $WXY$ principal block of $P_c$, with remaining column $b_3$, and normalize $Z_c=1$. For the fixed state $\psi_c=(\xi_c,1)$, ordinary matrix completion gives

$$
\psi_c^TP_c\psi_c
=(\xi_c+H_3^{-1}b_3)^TH_3(\xi_c+H_3^{-1}b_3)
+\frac{\det P_c}{\det H_3}.
$$

At the contact the first term has a positive semidefinite Hessian; the second has a positive definite Hessian by the strict Laurent Hessian and the positive determinant prefactor. Thus the leading two-by-two block of $G_\gamma$ is positive definite. This ordinary completion proves only a derivative statement at the source contact; no mixed inverse enters the SOS certificate.

Finally, suppose a nonzero $(x,y,z)$ lies in $\ker G_\gamma$. If $xyz\ne0$, the projective-direction argument above gives a second finite zero of $\gamma$. If $z=0$, the positive leading two-by-two block excludes it. If $x=0$, $y\ne0$, $z\ne0$, fix $s=z/y$ and send $t$ to infinity: $\gamma/t^2$ tends to $(0,1,s)G_\gamma(0,1,s)^T=0$, contradicting the projective boundary bound. Interchange $s,t$ for the case $y=0$, $x\ne0$. The remaining direction $(0,0,z)$ would make the joint-projective leading coefficient zero, with the same contradiction. Therefore $\ker G_\gamma=\{0\}$.

An ordinary source-only inverse $G_\gamma^{-1}$ is consequently available throughout these finite native families. No source-uniform lower eigenvalue bound is asserted at their degenerating limits.

## The remaining coupling condition

Each component of $\zeta$ vanishes at the contact. If its shifted coefficients are $a_{jk}$, then it has the exact representation $\zeta_i=V_{i\cdot}q$, where

$$
V_{i1}=a_{10}+a_{20}s,\qquad V_{i2}=a_{01}+a_{02}t,\qquad
V_{i3}=a_{11}+a_{21}s+a_{12}t+a_{22}st.
$$

Other choices of $V$ can represent the same coupling. For the native contacts just treated, the source-only inverse gives

$$
\begin{pmatrix}H&\zeta\\\zeta^T&\gamma\end{pmatrix}
=\begin{pmatrix}H-VG_\gamma^{-1}V^T&0\\0&0\end{pmatrix}
+\begin{pmatrix}V\\q^TG_\gamma\end{pmatrix}
G_\gamma^{-1}
\begin{pmatrix}V^T&G_\gamma q\end{pmatrix}.
$$

The second term is a fixed-feature SOS. A choice of $V$ admitting a uniform certificate for $H-VG_\gamma^{-1}V^T$ remains to be proved. Pointwise positivity of the original matrix only controls the smaller subtraction $\zeta\zeta^T/\gamma$: away from the contact, where $\gamma>0$,

$$
VG_\gamma^{-1}V^T-\frac{\zeta\zeta^T}{\gamma}
=VG_\gamma^{-1/2}\left(I-\frac{aa^T}{a^Ta}\right)G_\gamma^{-1/2}V^T\succeq0,
\qquad a=G_\gamma^{1/2}q.
$$

Thus the scalar certificate isolates the remaining three couplings; it does not complete the uniform state-matrix certificate.

The [native polynomial obstruction](C2_POLYNOMIAL_OBSTRUCTION.md) further
restricts this completion: an ordinary polynomial SOS for the full C2
parent requires $h=a^2r$. For generic sources, polynomial choices of $V$
cannot leave an ordinary polynomial SOS remainder. Fractional-power
certificates remain possible.
