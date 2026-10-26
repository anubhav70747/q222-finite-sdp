# An explicit bounded conductor for two state normals

The following polynomial construction writes a scalar remainder and two
complete first-normal terms as one quadratic form in eight specified roots.
It preserves bounded coefficients through both divisions.

## Ring, roots and sufficient hypotheses

Let $\mathcal D=\mathbb Q[\varepsilon,j,E,\xi_1,\ldots,\xi_m]$ and
$\mathcal R=\mathcal D[u,w,v]$. Source specialization sends every generator
of $\mathcal D$ to a bounded element of a real closed extension $K_c$ of $\mathbb R$,
with $\varepsilon,j,E>0$. The independent nonnegative testing variables
$u,w,v$ specialize to bounded elements of a separately labelled real closed
extension $K_t$ of $\mathbb R$. Write $\mathcal B=B_c\otimes_{\mathbb R}B_t$, where $B_c,B_t$
are the rings of elements bounded by real constants. A permitted localization
must keep $\varepsilon$ a nonunit, preserve the formal reductions and quotients,
have bounded inverses of its
pure source or pure testing denominators, and give bounded values to every
formal coquotient below. Polynomial inputs with bounded generators satisfy
these boundedness requirements automatically.

Put
$$
k=\varepsilon j,\quad v_0=\varepsilon j^2E,\quad
J=(u-k,w-j,v-v_0),\quad P_c=P(k,j,v_0),
$$
$$
A=\varepsilon wv,\quad B=v,\quad C=Euw,\quad
D_0=\varepsilon^2j^3E,\qquad
(x,y,z,h)=(B-C,A-kC,uB-kC,kC-D_0).
$$
Choose the nonnegative separated radicals
$$
a_1^2=A,\quad a_2^2=uB,\quad a_3^2=kC,\quad a_4^2=D_0,\qquad
q_1=a_1a_2-a_3a_4,\quad q_2=a_1a_3-a_2a_4.
$$
Equivalently,
$q_1=\sqrt\varepsilon\sqrt{uw}(v-v_0)$ and
$q_2=\varepsilon\sqrt{jE}\sqrt{uv}(w-j)$.
Each factor is pure source or pure testing. Direct multiplication gives
$$
r_*=y+z-h-2kx
=(v-v_0)\{u-k+\varepsilon(w-j)\}-Ek(u-k)(w-j)\in J^2,
\qquad kCr_*=q_1^2-yz=q_2^2+zh. \tag{1}
$$

Let $N,M$ be formal normal variables, with their actual forms substituted
after the polynomial identities are established. Suppose
$$
\mathscr F=F_0+2NG_N+2MG_M+aN^2+2bNM+cM^2,
$$
$$
F_0=f_AA^2+f_{BB}B^2+f_{BC}BC+f_{CC}C^2+\varpi D_0^2,\qquad
G_T=g_{TA}A+g_{TB}B+g_{TC}C\quad(T=N,M), \tag{2}
$$
with all displayed coefficients in $\mathcal R$. Assume:

1. $F_0\in J^2$ and $G_N,G_M\in J$ as formal polynomial identities.
2. For $S=k^2f_A+f_{BB}+f_{BC}+f_{CC}+k^2\varpi$,
   $S-S|_{w=j}\in\varepsilon^2\mathcal R$,
   $\overline S=s(u^2,v^2)$, and
   $\overline{2f_{BB}+f_{BC}}$ is independent of $w$.
3. Each $\overline{k g_{TA}+g_{TB}+g_{TC}}$ is independent of $w$.

Here a bar means reduction modulo $\varepsilon$; the other formal source
letters remain. All raw coefficients and all formal coquotients used below
must have bounded specializations.

**Theorem.** Under these hypotheses there is an explicitly constructed
symmetric $8\times8$ matrix $H$ over $\mathcal B$ such that
$$
\mathscr F=\mathscr R^TH\mathscr R,\qquad
\mathscr R=(N,M,x,y,z,h,q_1,q_2)^T. \tag{3}
$$
The construction uses polynomial arithmetic and monic division.

## The first-order division

For $P\in\mathcal R$, define ordered polynomial quotients
$$
T_u=\frac{P-P(k,w,v)}{u-k},\quad
T_w=\frac{P(k,w,v)-P(k,j,v)}{w-j},\quad
T_v=\frac{P(k,j,v)-P_c}{v-v_0}.
$$
Then $P=P_c+(u-k)T_u+(w-j)T_w+(v-v_0)T_v$.
If $P_c=0$ and $T_w\in(\varepsilon,v)$, form the polynomial quotients
$$
V_P=\frac{\overline{T_w}}v,\qquad
U_P=\frac{T_w-vV_P}{\varepsilon}.
$$
Lift $\overline{T_w}/v$ using its coefficients independent of $\varepsilon$.
The ideal test is equivalently that
$\overline P(0,w,0)$ is independent of $w$; reduction modulo
$\varepsilon$ is essential in this statement.

Define
$$
\begin{aligned}
\mathcal L_x(P)&=-uT_u-\varepsilon wU_P-Ekw^2V_P-EkwT_v,\\
\mathcal L_y(P)&=U_P+EjwV_P+EjT_v,\\
\mathcal L_z(P)&=T_u+Ew(w-j)V_P+EwT_v.
\end{aligned} \tag{4}
$$
These give $CP=x\mathcal L_x+y\mathcal L_y+z\mathcal L_z$.
Indeed, substitute the telescoping identity and use
$$
C(u-k)=z-ux,\quad
\varepsilon C(w-j)=y-\varepsilon wx,\quad
C(v-v_0)=E(wz+jy-kwx),
$$
$$
Cv(w-j)=Ew(w-j)z+Ejwy-Ekw^2x.
$$
If $T_w\in\varepsilon\mathcal R$, one may instead take $V_P=0$ and
$U_P=T_w/\varepsilon$; if $T_w\in\varepsilon^2\mathcal R$ this
$U_P$ has an additional factor $\varepsilon$.

## The scalar requires two divisions

Substituting $A=y+kC$, $B=x+C$, $D_0=kC-h$ into (2) gives
$$
F_0=Q+C\{x(2f_{BB}+f_{BC})+2kf_Ay-2k\varpi h\}+C^2S,
\qquad Q=f_Ay^2+f_{BB}x^2+\varpi h^2. \tag{5}
$$
At contact, $C_c=\varepsilon j^2E$ is a nonzero polynomial, so
$F_{0,c}=C_c^2S_c=0$ implies $S_c=0$ in the integral domain.
Apply (4) to $S$, choosing $V_S=0$ and
$U_S=T_w(S)/\varepsilon\in\varepsilon\mathcal R$. Write $L_i=\mathcal L_i(S)$.
Then $F_0=Q+C\sum_i r_iH_i$, where
$$
(r_x,r_y,r_z,r_h)=(x,y,z,h),\quad
(H_x,H_y,H_z,H_h)
=(2f_{BB}+f_{BC}+L_x,\ 2kf_A+L_y,\ L_z,\ -2k\varpi).
$$
Set $\lambda=2k\varpi_c$ and
$J_i=H_i-\lambda(-2k,1,1,-1)_i$.
The gradient of (5) at contact, together with (1), gives
$\sum_iJ_{i,c}\nabla r_{i,c}=0$, and $J_{h,c}=0$.
The remaining three gradients have determinant
$$
\det\frac{\partial(x,y,z)}{\partial(u,w,v)}\bigg|_c
=-E^2\varepsilon^3j^4\ne0.
$$
Cancellation in the formal domain therefore gives $J_{i,c}=0$ for all $i$;
no inverse of this determinant enters a coefficient.

It remains to justify the second division. Since $S_c=0$, $s(0,0)=0$.
Put
$$
p=\frac{s(u^2,v^2)-s(0,v^2)}{u^2},\qquad
q=\frac{s(0,v^2)}{v^2}.
$$
These are polynomials. Reduction modulo $\varepsilon$ gives
$T_u(S)\equiv up$, $T_v(S)\equiv vq$; separately,
$T_w(S)\in\varepsilon^2\mathcal R$. Thus (4) implies
$$
(J_x,J_y,J_z,J_h)\equiv
(\overline{2f_{BB}+f_{BC}}-u^2p,\ Ejvq,\ up+Ewvq,\ 0)
\pmod\varepsilon.
$$
The ordered $w$ quotient first sets $u=k$, whose reduction is zero.
Consequently
$$
(T_w(J_x),T_w(J_y),T_w(J_z),T_w(J_h))
\equiv(0,0,Evq,0)\pmod\varepsilon.
$$
Every $T_w(J_i)$ therefore lies in $(\varepsilon,v)$. Apply (4) again:
$CJ_i=x\ell_{ix}+y\ell_{iy}+z\ell_{iz}$.
Using (1) now yields the complete scalar identity
$$
F_0=Q+2\varpi_c(q_1^2-yz)
+\sum_{i=x,y,z,h}r_i(x\ell_{ix}+y\ell_{iy}+z\ell_{iz}). \tag{6}
$$

## The normal rows and the eight-root matrix

For $T=N,M$, put $S_T=kg_{TA}+g_{TB}+g_{TC}$.
Then $G_T=g_{TB}x+g_{TA}y+CS_T$, and $G_{T,c}=0$ implies $(S_T)_c=0$.
Hypothesis 3 gives $T_w(S_T)\in\varepsilon\mathcal R$, so (4), with
$V_{S_T}=0$, gives
$$
G_T=\gamma_{Tx}x+\gamma_{Ty}y+\gamma_{Tz}z,\quad
\gamma_T=(g_{TB}+\mathcal L_x(S_T),\
g_{TA}+\mathcal L_y(S_T),\ \mathcal L_z(S_T),0,0,0)^T. \tag{7}
$$
Let $e_x,e_y,e_z,e_h,e_{q_1},e_{q_2}$ be the standard basis of the soft
six-dimensional space and set
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
\end{aligned} \tag{8}
$$
Equation (1) identifies the two circuit pairings, so their equal split in
$H_s$ reproduces (6). Equation (7) supplies both complete first normals.
This proves (3). Every coefficient is a finite sum of products of bounded
generator images and the stipulated bounded coquotients. Monic division is
performed formally before specialization, with no inversion of a contact
difference or of $\varepsilon,j,E$. Hence $H$ has entries in $\mathcal B$.

The hypotheses concern the restricted support in (2); they do not cover
arbitrary additional terms involving $D_0$. Positive semidefiniteness of
$H$, and verification of these hypotheses for a particular native stress,
require separate arguments.
