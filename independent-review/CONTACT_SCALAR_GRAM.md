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

## The remaining coupling condition

Each component of $\zeta$ vanishes at the contact. If its shifted coefficients are $a_{jk}$, then it has the exact representation $\zeta_i=V_{i\cdot}q$, where

$$
V_{i1}=a_{10}+a_{20}s,\qquad V_{i2}=a_{01}+a_{02}t,\qquad
V_{i3}=a_{11}+a_{21}s+a_{12}t+a_{22}st.
$$

Other choices of $V$ can represent the same coupling. If $V$ annihilates $\ker G_\gamma$, its source-only generalized inverse gives

$$
\begin{pmatrix}H&\zeta\\\zeta^T&\gamma\end{pmatrix}
=\begin{pmatrix}H-VG_\gamma^\dagger V^T&0\\0&0\end{pmatrix}
+\begin{pmatrix}V\\q^TG_\gamma\end{pmatrix}
G_\gamma^\dagger
\begin{pmatrix}V^T&G_\gamma q\end{pmatrix}.
$$

The second term is a fixed-feature SOS. The range condition and a uniform certificate for $H-VG_\gamma^\dagger V^T$ remain to be proved. Pointwise positivity of the original matrix only controls the smaller subtraction $\zeta\zeta^T/\gamma$: for $G_\gamma\succ0$ and $\gamma>0$,

$$
VG_\gamma^{-1}V^T-\frac{\zeta\zeta^T}{\gamma}
=VG_\gamma^{-1/2}\left(I-\frac{aa^T}{a^Ta}\right)G_\gamma^{-1/2}V^T\succeq0,
\qquad a=G_\gamma^{1/2}q.
$$

Thus the scalar certificate isolates the remaining three couplings; it does not complete the uniform state-matrix certificate.
