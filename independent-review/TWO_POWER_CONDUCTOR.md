# A two-power conductor for two complete normal rows

Two source factors in the scalar coefficients and one in each first-normal
coefficient suffice for both contact divisions. No parity condition on the
testing polynomials is needed.

Let $K_c,K_t$ be separately labelled real closed extensions of $\mathbb R$,
let $\mathcal A\subset B_c$ be an $\mathbb R$-algebra domain of bounded source
elements, and put $R=\mathcal A[u,w,v]$. The independent testing variables
$u,w,v$ have bounded nonnegative values in $K_t$. Every element of $R$
therefore maps into $\mathcal B=B_c\otimes_{\mathbb R}B_t$.
Take positive $\varepsilon,j,E\in\mathcal A$ and define
$$
k=\varepsilon j,\quad v_0=\varepsilon j^2E,\quad
J=(u-k,w-j,v-v_0),\qquad P_c=P(k,j,v_0),
$$
$$
A=\varepsilon wv,\quad B=v,\quad C=Euw,\quad D_0=\varepsilon^2j^3E,\quad
(x,y,z,h)=(B-C,A-kC,uB-kC,kC-D_0).
$$
Write $(r_x,r_y,r_z,r_h)=(x,y,z,h)$.
Set
$q_1=\sqrt\varepsilon\sqrt{uw}(v-v_0)$ and
$q_2=\varepsilon\sqrt{jE}\sqrt{uv}(w-j)$.
All radicals are pure source or pure testing. Direct multiplication gives
$$
r_*=y+z-h-2kx\in J^2,\qquad
kCr_*=q_1^2-yz=q_2^2+zh. \tag{1}
$$

**Theorem.** Let $N,M$ be formal normal variables and suppose
$$
\mathscr F=F_0+2NG_N+2MG_M+aN^2+2bNM+cM^2,
$$
$$
F_0=f_AA^2+f_{BB}B^2+f_{BC}BC+f_{CC}C^2+\varpi D_0^2,\qquad
G_T=g_{TA}A+g_{TB}B+g_{TC}C\quad(T=N,M).
$$
Assume $a,b,c\in R$, $F_0\in J^2$, $G_N,G_M\in J$, and supply exact
factorizations
$$
f_i=\varepsilon^2\widehat f_i,\qquad
g_{T\bullet}=\varepsilon\widehat g_{T\bullet},
\qquad \widehat f_i,\widehat g_{T\bullet}\in R. \tag{2}
$$
Then the construction below gives a symmetric $8\times8$ matrix $H$ over
$\mathcal B$ with
$\mathscr F=\mathscr R^TH\mathscr R$ for
$\mathscr R=(N,M,x,y,z,h,q_1,q_2)^T$.

The factors in (2) must be bounded polynomial coquotients, rather than
merely small values of the raw expressions. Testing denominators are
cleared before applying the theorem. Any source localization used in $R$
must still lie in $B_c$, and every supplied coquotient must remain in $R$.

## Two filtered divisions

For $P\in R$, let $T_u(P),T_w(P),T_v(P)$ be the ordered monic quotients
$$
T_u=\frac{P-P(k,w,v)}{u-k},\quad
T_w=\frac{P(k,w,v)-P(k,j,v)}{w-j},\quad
T_v=\frac{P(k,j,v)-P_c}{v-v_0}.
$$
These are polynomial divided differences, so they preserve every
source factor $\varepsilon^m$. If $P_c=0$ and $T_w(P)=\varepsilon U_P$
with $U_P\in R$, define
$$
\mathcal L_x=-uT_u-\varepsilon wU_P-EkwT_v,\qquad
\mathcal L_y=U_P+EjT_v,\qquad
\mathcal L_z=T_u+EwT_v. \tag{3}
$$
Then $CP=x\mathcal L_x+y\mathcal L_y+z\mathcal L_z$. This follows by
telescoping $P$ and using
$C(u-k)=z-ux$,
$\varepsilon C(w-j)=y-\varepsilon wx$, and
$C(v-v_0)=E(wz+jy-kwx)$.

For the scalar, put
$$
S=k^2f_A+f_{BB}+f_{BC}+f_{CC}+k^2\varpi,\qquad
Q=f_Ay^2+f_{BB}x^2+\varpi h^2.
$$
Substitution of $A=y+kC$, $B=x+C$, $D_0=kC-h$ gives
$$
F_0=Q+C\{x(2f_{BB}+f_{BC})+2kf_Ay-2k\varpi h\}+C^2S. \tag{4}
$$
Since $C_c=\varepsilon j^2E\ne0$, $F_{0,c}=0$ implies $S_c=0$ by
cancellation in the source domain. Write $S=\varepsilon^2S_2$.
The first division uses $U_S=\varepsilon T_w(S_2)$ in (3), giving
$$
(\mathcal L_x(S),\mathcal L_y(S),\mathcal L_z(S))
\in\varepsilon^2R\times\varepsilon R\times\varepsilon^2R.
$$
Define
$$
(H_x,H_y,H_z,H_h)
=(2f_{BB}+f_{BC}+\mathcal L_x(S),\
2kf_A+\mathcal L_y(S),\ \mathcal L_z(S),\ -2k\varpi),
$$
and $\lambda=2k\varpi_c$, $J_i=H_i-\lambda(-2k,1,1,-1)_i$. Then
$$
(J_x,J_y,J_z,J_h)\in
\varepsilon^2R\times\varepsilon R\times\varepsilon^2R\times\varepsilon^3R.
\tag{5}
$$
Every $J_{i,c}$ vanishes. Indeed, differentiating (4) after the first
division and using $F_0,Q,r_*\in J^2$ gives
$\sum_iJ_{i,c}\nabla r_{i,c}=0$, while $J_{h,c}=0$.
The other three gradients have determinant
$-E^2\varepsilon^3j^4\ne0$, so cancellation proves the assertion.
No inverse of that determinant enters a coefficient.

Write $J_i=\varepsilon J_{i,1}$ using (5). The second division uses
$U_{J_i}=T_w(J_{i,1})$ in (3), with no additional divisibility or parity
test. Write $CJ_i=x\ell_{ix}+y\ell_{iy}+z\ell_{iz}$. Equation (1) gives
$$
F_0=Q+2\varpi_c(q_1^2-yz)
+\sum_{i=x,y,z,h}r_i(x\ell_{ix}+y\ell_{iy}+z\ell_{iz}). \tag{6}
$$

## Both normal rows and the matrix

For $T=N,M$, set $S_T=kg_{TA}+g_{TB}+g_{TC}=\varepsilon S_{T,1}$.
The identity $G_T=g_{TB}x+g_{TA}y+CS_T$ and $G_T\in J$ imply $(S_T)_c=0$.
Use $U_{S_T}=T_w(S_{T,1})$ in (3). It yields
$$
G_T=\gamma_T^T(x,y,z,h,q_1,q_2)^T,\quad
\gamma_T=(g_{TB}+\mathcal L_x(S_T),\
g_{TA}+\mathcal L_y(S_T),\ \mathcal L_z(S_T),0,0,0)^T. \tag{7}
$$
In the soft six-dimensional basis let $e_i$ be coordinate vectors and
$\widetilde\ell_i=(\ell_{ix},\ell_{iy},\ell_{iz},0,0,0)^T$. Define
$$
\begin{aligned}
H_s={}&f_{BB}e_xe_x^T+f_Ae_ye_y^T+\varpi e_he_h^T
+\varpi_c(e_{q_1}e_{q_1}^T+e_{q_2}e_{q_2}^T)\\
&-\tfrac{\varpi_c}{2}(e_ye_z^T+e_ze_y^T)
+\tfrac{\varpi_c}{2}(e_ze_h^T+e_he_z^T)
+\tfrac12\sum_{i=x,y,z,h}(e_i\widetilde\ell_i^T+\widetilde\ell_ie_i^T),\\
H={}&\begin{pmatrix}a&b&\gamma_N^T\\b&c&\gamma_M^T\\
\gamma_N&\gamma_M&H_s\end{pmatrix}.
\end{aligned}
$$
The equal circuit split in (1) makes this exactly (6), and (7) supplies
both complete normal terms. All divided differences and supplied
coquotients lie in $R$; thus every coefficient is bounded. This proves
the theorem without introducing a singular source or mixed inverse.

For the CA/SE change of soft normalization, take $\delta>0$. If each CA monomial is
$\delta$ times its SE counterpart, so is each of the six soft roots.
The unchanged hard normals give
$D=\operatorname{diag}(1,1,\delta,\delta,\delta,\delta,\delta,\delta)$ and
$H_{\rm SE}=D^TH_{\rm CA}D$. This preserves boundedness and strict real
residue only when $\delta$ and $\delta^{-1}$ are bounded source units.

The theorem constructs a symmetric identity. Positivity of $H$ and the
hypotheses for a particular native stress require separate proofs.
