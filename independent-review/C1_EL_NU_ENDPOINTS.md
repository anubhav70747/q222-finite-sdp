# Native C1 endpoint EL/NU: independently checked separated tensor certificates

Scope (26 October 2026). Pinned original source: c1_short.tex at 8960b3fead3ce6b185bb1387d148fb87d05d3d23, SD9–SD11, ES1–ES3/SE2–SE3 in component p8c30238a5f7d, CE1–CE10, EB4–EB9, FU6–FU23, NU1–NU26, and EL0–EL21. The claim here is only exact tensor SOS on the two endpoint source domains, not global existence or an evaluated lift. The proof uses the common exact source physicality and strict positive-residue tensor-Gram lemma as stated in the foundational chapter. All sources are in a real-closed field K_c, all independent test and state values in a distinct labelled real-closed K_t. Bounded Gram entries lie in B_c ⊗_R B_t; roots and positive pure-field prefactors may be unbounded. Every sign-restricted circuit is used only on its corresponding testing sign chart.

## Theorem

Let the physical native C1 source have d=1, lambda_event=n=0, m,C,A,alpha>0, D,beta>=0, D+beta>0, a noncommuting finite stationary contact, all original strict native gaps, and infinitesimal total-scale S=sigma+L. Put

$$
x=\kappa S,\quad\nu=\sigma/S,\quad\ell=L/S,\quad \nu+\ell=1,\quad
 \theta=1-h,\quad e=\kappa h/\nu,\quad E=1-e\ge1/18,
$$

with 0<r,h<1, positive native J,M,Q_s, and zeta=ell/kappa=L/x noninfinitesimal. If either

- NU: kappa,ell are positive bounded units and nu>0 is infinitesimal; or
- EL: nu is a positive unit and ell,kappa>0 are infinitesimal, with kappa/ell bounded,

then the entire original four-state parent at every independent projective testing pair and state is a finite SOS in K_c ⊗_R K_t, including either exact D=0 or beta=0 where physical, independent infinitesimal rates, both unrestricted state normals, axes, reciprocals, and complex states. No inverse of a mixed source/test element is used.

## I. Native source and contact, without assumed first jets

Write R_s=1+xr,

$$
\begin{aligned}
 J&=1+xh-\nu S(1-r),&M&=1+xr-\kappa\theta/\ell,& Q_s&=1-\ell Sr,\\
 P_0&=1+S(\nu\ell+2\nu\kappa h+\kappa\ell+\kappa^2h\theta)
 +S^2(\nu r\kappa\ell+\nu r\kappa^2h),\\
 U_s&=(1+x)(1+\nu Sr).
\end{aligned}
$$

The actual normalized source coefficients are

$$
\begin{aligned}
 c_0:=SC(1-r)&=\frac{U_sE}{2P_0},\\
 a_0:=SAr&=\frac{(1+xh)E+(1-r)(\ell S+xh)}{2P_0},\\
 D/S&=\frac{U_s-JP_0}{2JSP_0},\\
 a_1:=\alpha/S&=\frac{(1+\nu Sr)[\ell R_s e+\ell^2ME]}{2Q_sP_0},\\
 \beta\theta&=\frac{U_s(\ell+\kappa h)-\theta P_0}{2P_0},\\
 m/\kappa^2&=\frac{U_sh}{2M\ell\nu P_0}.
\end{aligned}\tag{1}
$$

Reconstruction, rather than reliance on a label. In original SD variables (k^2=x, eta=1/sigma, zeta=L/x, a_SD=sigma(1-r), t_SD=x theta/L), the identities r=1-a_SD eta, h=1-t_SD zeta, and y=x eta zeta=L/sigma are literal. The original half-angle reciprocal multiplier is Delta=sqrt(sigma L)/x. SD11 consequently gives

$$
u_c=\sqrt{x},\quad v_c=\sqrt{L/\sigma},\quad
 (W_c,X_c,Y_c,Z_c)=\left(\sqrt{\sigma L/x},-h\sqrt{x},r\sqrt{\sigma L},1\right).
\tag{2}
$$

The original NC4 source formulas transform into (1): P=sigma P_d=SP_0, P_4=(1+xh)E+(1-r)(L+xh), the numerator of D is S U_s-JP, the beta numerator is U_s(L+xh)-theta P_0 S, and x P_alpha=S R_s e+LM E. The two SE2/SE3 source-sphere equations are verified directly from these polynomial formulas:

$$
\begin{aligned}
 2\theta P\beta&=h B-\sigma H_0,\\
 2P(JD+x\theta\beta)&=\sigma^2\mathcal H,\\
 H_0&=1+L(1-r),\quad B=P+x(J-\sigma),\\
 \mathcal H&=R_sH_0+(1-r)xh(2+xr)-e[2+x(1+h+r\theta)]+(1+x)e^2.
\end{aligned}\tag{3}
$$

U_s-JP_0 has an exact factor S because its constant term is zero; its quotient has unit denominators. The six SD9/SD10 native source balances have been checked with exact polynomial expansion (script verify_sd_source.py). Differentiating the actual complete SD2 parent produces SD11's state and angular jets via these six identities, so all source contact derivatives used below are exact. At the contact (2) we have

$$
Y_c=r u_cW_c,\quad X_c=-hu_c,\quad W_c=(L/x)u_c/v_c,
\quad {X_c+u_cZ_c\over v_c(W_c+u_cY_c)}=p:=\frac{\kappa\theta}{\ell R_s}.
\tag{4}
$$

In particular both NU normals below and the exact centered NU event row vanish at contact, along with the first derivatives of that row's square. The EL circuit/normal roots vanish for the same (2).
The sphere is specifically SE2/SE3 in component p8c30238a5f7d, not the later similarly numbered SE conductor. Writing nonnegative sphere coordinates d_sph^2+b_sph^2=1, equation (3) yields

$$
D={\sigma^2\mathcal H\over2PJ}d_{\rm sph}^2,\qquad
\beta={\sigma^2\mathcal H\over2Px\theta}b_{\rm sph}^2.
\tag{5}
$$

No sphere coordinate is divided by.

## II. Correct finite-helper Gram lemma

If a bounded Gram on the main roots has positive-definite real residue, and all additional main-main terms have zero residue, then finitely many helper roots with bounded arbitrary main-helper and helper-helper couplings can be supplied with sufficiently large fixed real positive diagonals, provided explicit SOS budget identities purchase those diagonals with zero-residue cost on main roots and with noninfinitesimal costs only on already-acquired helpers. Purchase in triangular order, choosing earlier helper pivots sufficiently large for later costs. The resulting real block Schur complement is strictly positive. This does not hold for an arbitrary bounded main-main error: such an error could destroy its positive residue. The NU and EL checks below explicitly establish the required main-block property before any helper budgets are used.
For each budget C_i-H_i^2=\sum_j t_{ij}^2, insert M_i(E_{H_iH_i}-Q(C_i)) into the Gram and retain M_i\sum_jt_{ij}^2. This preserves the exact whole polynomial. Matrix inverses used to choose the real finite M_i are real-only, not tensor coefficients.

## III. NU full local parent and bounded coquotients

Here nu is infinitesimal, kappa,ell units and q=nu/ell infinitesimal. Dividing (5) gives

$$
d_s={D\over Sq^2}={\ell^2\mathcal H\over2P_0J}d_{\rm sph}^2,
 \qquad n_s={\beta\over q^2}={\ell^2\mathcal H\over2\kappa P_0\theta}b_{\rm sph}^2,
\quad Jd_s+\kappa\theta n_s={\ell^2\mathcal H\over2P_0}.
\tag{6}
$$

Every denominator is a positive source unit. The last right-hand side is a positive unit; hence d_s,n_s are bounded nonnegative and at least one has positive unit residue even if the other is literally zero. The first equation in (3), divided by S, gives e_0=kappa_0/(1+kappa_0), E_0=1/(1+kappa_0), h infinitesimal, theta unit, and c_0,a_0,a_1 all have residue E_0/2. Since M=R_s(1-p),

$$
\mu_c=m(1-p)={\kappa e U_s\over2\ell P_0R_s}
\longrightarrow {\kappa_0^2\over2(1+\kappa_0)}>0.
\tag{7}
$$

Choose real 0<epsilon<mu_c. In the homogeneous reciprocal quarter-parent P_w=w^2P(u,1/w) retain the original positive, contact-zero square
$\Pi_m=(m-\epsilon)[w(X+uZ)-p(W+uY)]^2.$ The exact remainder of the original m square has bounded coefficients

$$
m_0B_c^2-2m_1B_cV+m_2V^2,\quad
 m_0=\epsilon,\ m_1=\mu_c+\epsilon p,\ m_2=\mu_c(1+p)+\epsilon p^2,
\tag{8}
$$

where B_c=w(X+uZ), V=W+uY. The retained row vanishes by (4), hence has zero value and all first jets at the source contact. It is used even if m diverges; no unbounded m is left in the remainder.
On independent infinitesimal testing u,w and arbitrary states, put a=u^2,c=w^2, delta=c-q, t=sqrt(S), and use all the roots

$$
X_e=X+ewW,\ N_y={Y-ruW\over\sqrt{Sr(1-r)}},\quad
 N=W-\zeta uwZ,\quad G={uW\over t}-\ell twZ,\quad
 B=t(c-q)Z,\ B_u=u(c-q)Z,
\quad\zeta=\ell/\kappa.
\tag{9}
$$

The original quarter-parent minus Pi_m is, exactly in both unrestricted normals,

$$
F_0+H_xX_e^2+H_yN_y^2+2H_{xy}X_eN_y+2L_xX_e+2L_yN_y.
\tag{10}
$$

For a complete finite coefficient specification define

$$
\begin{aligned}
 A_N={}&m_2(1+ra)^2+2m_1ec(1+ra)+m_0e^2c^2+(\beta+\tfrac12)e^2c(1+c)-de[2(1-r)-e]ac,\\
 A_P={}&(1+c)[c_0(1-r)+a_0r(1+a)] + S[d(1-r)^2c+(\alpha+\tfrac12)(1+a)r^2],\\
 K={}&m_1(1+ra)+m_0ec+d(1-r-e)+(\alpha+\tfrac12)(1+a)r+(\beta+\tfrac12)e(1+c),\\
 F_Z={}&D+\alpha c+\beta a+(m_0+\alpha+\beta+\tfrac12)ac,\\
 J_1={}&K-\ell A_P,\quad L_1=J_1-\zeta A_N,\quad
 H=F_Z-\ell^2ScA_P+ac(\zeta^2A_N-2\zeta J_1),\quad d=D+\tfrac12.
\end{aligned}\tag{11}
$$

Here and below A_P has precisely the factor (1+a) in its alpha contribution as in NU11. Then

$$
F_0=A_NN^2+A_PG^2-2uwZL_1N+HZ^2.
\tag{12}
$$

The complete first-normal and hard arrays, obtained by differentiation of the original P_w (verified symbolically for unrestricted W,Z,X_e,N_y), are

$$
\begin{aligned}
 A_x&=-m_0ec-m_1(1+ra)+da(1-r-e)-e(\beta+\tfrac12)(1+c),\\
 T_x&=\zeta cA_x+\beta-D+c(m_0+\beta+\tfrac12),\\
 Y_P&=(1+c)[-c_0+a_0(1+a)]+S[m_2(1+ra)+m_1ec-dc(1-r-e)+(\alpha+\tfrac12)(1+a)r],\\
 Y_C&=\ell Y_P-m_1a+d-(\alpha+\tfrac12)(1+a),\\
 L_x&=wA_xN+uZT_x,\qquad L_y=\sqrt{r(1-r)}[Y_PG+Y_CtwZ],\\
 H_x&=m_0c+da+(\beta+\tfrac12)(1+c),\\
 H_y&=rc_0(1+c)+(1-r)a_0(1+a)(1+c)
   +Sr(1-r)[m_2a+dc+(\alpha+\tfrac12)(1+a)],\\
 H_{xy}&=-uw\sqrt{Sr(1-r)}(m_1+d).
\end{aligned}\tag{13}
$$

Every source coefficient is bounded because D=Sq^2d_s, beta=q^2n_s, alpha=Sa_1, c_0,a_0,a_1,m_i bounded. No division by r,1-r is left in (11)–(13).
At the native contact (2), all four graph roots and the retained event row vanish. State and two angle first jets, proved by SD9–SD11, force L_1(x,q)=0, H(x,q)=H_a(x,q)=H_c(x,q)=0, T_x(x,q)=Y_C(x,q)=0. The crucial finite monic divisions are in the order a-x, then c-q:

$$
\begin{aligned}
 L_1&=(a-x)L_a+(c-q)L_c,\\
 H&=H(x,c)+(a-x)H_a(x,c)+c(a-x)^2H_{20},\\
 H_a(x,c)&=(c-q)[-\beta/q+cH_{11}],\\
 H(x,c)&=S(c-q)^2[d_s+\kappa n_s+cH_{02}],\\
 T_x&=c(a-x)X_a+(c-q)[(D-\beta)/q+cX_c],\\
 Y_C&=(a-x)Y_a+(c-q)Y_c.
\end{aligned}\tag{14}
$$

Why every coquotient is bounded: H(a,0)=D+beta a exactly, and T_x(a,0)=beta-D. Thus the only apparent quotients are beta/q=q n_s and (D-beta)/q=q(Sd_s-n_s), both bounded, including both literal zero faces. At a=x=kappa S, every monomial of H(x,c) contains S: F_Z does because D,alpha=O(S) and beta x=O(S); the other terms of (11) visibly contain S or a=x. Consequently H(x,c)/S is a bounded coefficient polynomial before its monic (c-q)^2 division. The factor c in the a-quadratic follows from H(a,0) affine in a. Similarly T_x has an explicit factor c in every a-dependent term. All other divisions use monic polynomials over the unit-denominator source ring, so boundedness is coefficientwise, not inferred from evaluation sizes.
Use the exact separated identity

$$
c(a-x)(c-q)Z^2=\zeta^{-1}[G(wB)-N(wB_u)],
 \qquad R=w(a-x)Z=(\sqrt S/\zeta)G-(u/\zeta)N.
\tag{15}
$$

Substitution of (14) into (12) and (13) gives, with all terms retained,

$$
\begin{aligned}
 F_0={}&A_NN^2+A_PG^2+d_sB^2+n_sB_u^2+H_{02}(wB)^2+H_{20}R^2\\
 &+(H_{11}-n_s)\zeta^{-1}[G(wB)-N(wB_u)]-2uL_aNR-2L_cN(wB_u),\\
 L_x={}&wA_xN+uwX_aR+d_s\sqrt S\,quB-n_sqB_u+wX_c(wB_u),\\
 L_y={}&\sqrt{r(1-r)}[Y_PG+\sqrt S Y_aR+Y_c(wB)].
\end{aligned}\tag{16}
$$

The coefficient Gram's main block, in (N,G,X_e,N_y,sqrt(d_s)B,sqrt(n_s)B_u), has already positive-definite residue

$$
\operatorname{diag}(\kappa_0^2(\epsilon+1/2), E_0/2,1/2,E_0/2,1,1).
\tag{17}
$$

All main-main cross terms have zero residue: R has factors sqrt S or u; Y_P has zero residue from c_0=a_0=E_0/2; the face-safe d_s,n_s normal terms have extra factors sqrt S q u or q; and H_xy contains uw sqrt S. This remains true if a scaled angular root vanishes identically.
The two exact SOS budgets, delta=c-q, are

$$
\begin{aligned}
\delta^2N^2+2\zeta\delta wGB+2\zeta\ell w^2B^2-\zeta^2(wB_u)^2&=(\delta W)^2,\\
\delta^2G^2+2\ell\delta wNB_u+2\ell\zeta w^2B_u^2-\ell^2(wB)^2&=(\delta uW/\sqrt S)^2.
\end{aligned}\tag{18}
$$

They were re-expanded exactly. If d_s is unit, purchase wB_u first: the cost on sqrt(d_s)B is O(w^2/d_s) with zero residue; purchase wB second, charging an already acquired wB_u pivot. If n_s is unit, reverse. Every other cost on the main roots has delta or w and zero residue; ell,zeta are units in NU. Hence arbitrarily large real helper diagonals absorb the bounded helper couplings in (16). By the corrected lemma, the eight-row residue is strictly positive. Restore the right sides of (18) and Pi_m. Literal D=0 or beta=0 simply removes an original nonnegative square; the alternate sphere pivot still buys both helpers.

## IV. EL full local parent and its exceptional divisions

Here ell,kappa are infinitesimal, nu unit and g=kappa/ell bounded. beta theta>=0 in (1) gives theta=O(ell) and beta theta=O(ell); hence e=O(kappa), M,J,E~1, c_0,a_0~1/2, D/S~1/2, m=O(kappa^2/ell)=O(kappa), and

$$
\frac{\alpha}{S\ell^2}-\frac{1+g}{2}\in\operatorname{Inf}.
\tag{19}
$$

The face D=0 is therefore impossible here, but beta=0 remains allowed. Set q=ell/nu, b=sqrt(q), j=nu b, p=kappa theta/(ell R_s), and define the unrestricted normals and main angular roots

$$
N=X+uZ-pv(W+uY),\quad N_y={Y-ruW\over\sqrt{Sr(1-r)}},\quad
 P=uW/\sqrt S,\ Q=\sqrt S Z,\quad
 \Gamma=P-\nu vQ,\ \Lambda=vP-\ell Q.
\tag{20}
$$

The original finite-chart P=F/4 is exactly quadratic in both normals. Its scalar graph is

$$
F_0=c(1+ra)^2(L_0+L_1c)W^2+A_PP^2-2vKPQ+F_ZZ^2,
\quad a=u^2,c=v^2,
\tag{21}
$$

where L_0=mu(1-p)+(beta+1/2)p^2, L_1=(beta+1/2)p^2, mu=m(1-p), d=D+1/2, B_c=1-r+pc(1+ra),

$$
A_P=(1+c)[c_0(1-r)+a_0r(1+a)]+S[dB_c^2+(\alpha+1/2)(1+a)r^2c],\quad
 K=(1+a)[dB_c+(\alpha+1/2)r],\quad F_Z=(1+a)[\alpha+Dc+dac].
$$

Set lambda(a,c)=(lambda_0+lambda_1 c)(1+ra)^2, lambda_i=qL_i/kappa^2. The exact cancellations of the apparent 1/ell poles are

$$
\lambda_1={\theta(2\beta\theta+\theta)\over2\ell\nu R_s^2}=O(\ell),\qquad
\lambda_0={U_shM\over2\nu^2P_0R_s^2}+\lambda_1\sim1/2.
\tag{22}
$$

Define widehat A=A_P-kappa lambda, widehat K=K-2kappa nu lambda, widehat F=F_Z-lambda(nu^2 ac+kappa nu ell S), J=widehat K-nu widehat A, and H=widehat F+nu^2 Sc widehat A-2nu Sc widehat K. Then exact expansion of the whole scalar gives

$$
F_0=\lambda\Psi+\widehat A\Gamma^2-2vJ\Gamma Q+HZ^2,
\quad\bar{\widehat A}=1/2,
\tag{23}
$$

where the genuine separated four-monomial circuit is

$$
\Psi=\sum_{i=1}^4 A_i^2-4\kappa\nu uvWZ,
\quad(A_1,A_2,A_3,A_4)=((\kappa/b)vW,\nu uvZ,\nu b\sqrt x Z,\sqrt\kappa P).
$$

On the positive-product sign chart its two paired AM–GM decompositions reserve positive real 1/16 pivots on E12=A1-A2 and E23=A2-A3, leaving an exact SOS remainder. On the opposite-product chart the signed cross term itself is a separated positive square; the remaining constant real quadratic is positive. No mixed radical of lambda is taken, as lambda-1/4 is a nonnegative sum of separated monomial squares by (22).
At (a,c)=(x,q) the native contact (2) makes both normals, both circuit differences, Gamma, and the auxiliary contact rows zero. Their squares have zero contact value and first jets. The six native first derivatives imply J(x,q)=0 and H(x,q)=H_a(x,q)=H_c(x,q)=0. Finite monic division in a-x, then c-q gives

$$
J=(a-x)J_a+(c-q)J_c,\quad
 H=(a-x)^2H_{20}+(a-x)(c-q)H_{11}+(c-q)^2H_{02}.
\tag{24}
$$

The four critical coquotients, proved by evaluating the exact formula for H(a,0) and the coefficient of c^3, are

$$
\begin{aligned}
 H_{20}(a,0)&=-\kappa\nu\ell S\lambda_0r^2,\\
 H_{11}(0)&=-\alpha/q+2\kappa\nu^2S\lambda_0r(1+rx),\\
 H_{02}(0)/S&=\nu^2[(1+x)\alpha/(S\ell^2)-(\kappa/\ell)\nu\lambda_0(1+rx)^2]\sim1/2,\\
 H_{02}(c)&=H_{02}(0)+\nu^2S^2dp^2(1+rx)^2c.
\end{aligned}\tag{25}
$$

Thus H11(0)/(S ell) is bounded because alpha/(S ell^2) and kappa/ell are bounded; H02(c)/S has positive-unit residue. Do not infer any of these from the undivided weight of H alone. The remaining coefficients tilde H20, tilde H11 are bounded monic polynomial coquotients over the source unit-denominator ring.
The complete first-normal and hard arrays, checked directly from the full original parent, are

$$
\begin{aligned}
 B_0(c)&=[-\mu+(\beta+1/2)p(1+c)]/\kappa,\\
 F_a&=r\kappa B_0(c)+d(1-r)+dpc(1+ra),\quad D_t=d(1+a)-\nu SF_a,\\
 L_N&=bB_0(c)E_{12}+uv\sqrt S F_a\Gamma-uvD_aR_v
 +[-D_t(x,q)+\nu b B_{01}(v+b)-v(v+b)D_c]V_h,\\
 L_y&=\sqrt{r(1-r)}[Y_P\Gamma+\sqrt S Y_aR_v+vY_cB],\\
 H_N&=m+dac+(\beta+1/2)(1+c),\\
 H_y&=rc_0(1+c)+(1-r)a_0(1+a)(1+c)
  +Sr(1-r)[d(-1+pac)^2+(\alpha+1/2)(1+a)c+ac(L_0+L_1c)],\\
 H_{Ny}&=\sqrt{Sr(1-r)}uv[\kappa B_0(c)-d+dpac].
\end{aligned}\tag{26}
$$

Here Y_P=(1+c)[-c_0+a_0(1+a)]+S[c(1+ra)(L_0+L_1c)+d(-1+pac)B_c+(alpha+1/2)(1+a)rc], Y_K=nu Y_P+(1+a)(D-alpha-dpac)=(a-x)Y_a+(c-q)Y_c, D_t(a,c)-D_t(x,q)=(a-x)D_a+(c-q)D_c and B_0(c)-B_0(q)=(c-q)B_{01}. All divided arrays are bounded: mu/kappa=O(kappa/ell), (beta+1/2)p/kappa=(beta theta+theta/2)/(ell R_s)=O(1), and the source ratios M,J,P_0,R_s are units.
Independent readback of the delicate L_N identity. Direct differentiation before any source simplification yields
$L_N=v[\kappa B_0(c)W+aF_aW-duD_uvZ].$ Using E12=(v/b)(kappa W-nu b uZ) and Gamma=P-nu vQ, the difference between this derivative and (26)'s first-normal formula is identically

$$
uvZ\,b[\nu B_0(q)-D_t(x,q)].
\tag{27}
$$

At the source contact all right-hand roots in (26) vanish; the left derivative vanishes by SD9–SD11, so the pure-source bracket equals zero. No coefficient is bounded by dividing a mixed element. This proves the entire first-normal identity analytically, not by claiming a scalar Schur consequence. The L_y, H_N,H_y,H_Ny equalities are separately verified by unrestricted symbolic state expansion.
With B=sqrt S(c-q)Z, R_v=v(a-x)Z, R_0=ell sqrt S(a-x)Z, V_h=v(v-b)uZ, the final scalar identity, obtained by applying (24), is

$$
\begin{aligned}
 F_0={}&\lambda\Psi+\widehat A\Gamma^2+H_{02}(c)B^2/S
  -2\sqrt S J_a\Gamma R_v-2vJ_c\Gamma B\\
 &-(\kappa/\ell)\nu\lambda_0r^2R_0^2+\widetilde H_{20}R_v^2
  +[H_{11}(0)/(S\ell)]R_0B\\
 &+\widetilde H_{11}[ (c-q)E_{23}^2/\nu^2+(2\sqrt{\kappa q}/\nu)E_{23}B-\kappa B^2].
\end{aligned}\tag{28}
$$

The last bracket is exactly c(a-x)(c-q)Z^2, proved by expanding E23=nu(uv-sqrt(xq))Z; it is not a dropped mixed-radial term. The main six roots (E12,E23,Gamma,B,sqrt(1+beta)N,N_y) already have strictly positive full residue, with diagonal bounded below by (1/16,1/16,1/2,1/2,1/2,1/2) and every main-main off-diagonal having zero residue. This follows from (25)–(28), Y_P having zero residue, and the source-only hard normalization. The auxiliary signed terms act only on (V_h,R_v,R_0) and their bounded couplings to the main list.
The exact SOS budget identities EL8–EL11 buy real helper pivots in order V_h, R_v, R_0. Their main costs are infinitesimal because

$$
\kappa/j=(\kappa/\ell)b\in\operatorname{Inf},\qquad
 x^2/(qS)=\kappa S\cdot(\kappa\nu/\ell)\in\operatorname{Inf}.
$$

The independent exact EL8 SOS remainder and all EL9–EL11 and EL18 polynomial identities are reproduced by verify_budgets.py. The two later budgets charge only earlier helpers and zero-residue main coefficients. By the corrected helper lemma a bounded strictly positive-residue nine-root Gram results. Restoring the budget SOS remainders proves the entire original parent, retaining both normals, arbitrary state ratios and either sign chart. All unbounded beta-dependence is absorbed into sqrt(1+beta)N with diagonal residue at least 1/2, not treated as a bounded error.

## V. Independent projective testing joins, with scopes checked

CE (both endpoints). c_0,a_0 are positive units and S infinitesimal, hence C=c_0/[S(1-r)] and A=a_0/[Sr] are unbounded. For a noninfinitesimal Alice half-angle the original positive C,A events may be capped at a finite real M. Exact four-state CE2 and homogeneous CE6 identities were re-expanded. On the common kernel of the two capped events the real signed core is strictly positive: for bounded unit u, its (B=X+uZ,V=vZ) determinant is D_u u^2>0; the u-reciprocal and Bob-reciprocal forms CE4–CE8 have positive definite full-state residues. Choose M by a real finite-dimensional penalty Schur complement. Retain all excess original events whole. CE uses no EL/NU source theorem, so the join is noncircular.
EB/SE unmatched Bob (EL). EL has e,m infinitesimal, D/S-\nu^2/2, alpha/S-ell^2/2 infinitesimal, and c_0,a_0\sim1/2. These are exactly the source bounds required by the homogeneous all-state EB4–EB7 parent, even when zeta=ell/kappa is a unit (the original EB subsection's narrower zeta-unbounded source premise is not imported). The EB5 core and all three EB6 event amplitude identities were re-expanded independently; their harmonic C/A completion is exact. After retaining the positive ASu^2y^2 and beta(X+uZ)^2, the six-root Gram residue is diag(1,1,1/2,1/2) plus

$$
\frac12\begin{pmatrix}1&-s_0c_0^{\rm test}\\-s_0c_0^{\rm test}&\nu_0^2s_0^2+\ell_0^2(c_0^{\rm test})^2\end{pmatrix},\quad
 \det=\frac14(s_0^2-\ell_0)^2.
\tag{29}
$$

Here (c^{test},s) are normalized pure-testing homogeneous Bob coordinates. In EL ell_0=0, so all Bob tests with nonzero sine-square residue, including Bob infinity, are strict. The remaining Bob chart has v infinitesimal or literally zero and is the EL local chart already proved. No mixed inverse is used by EB.
FU finite Bob (NU). The exact FU8 scalar, complete FU9 first normal and FU9 hard coefficient have each been re-expanded without using the FU subsection's nu-unit source assumptions. For NU, D/S=q^2d_s and beta=q^2n_s tend to zero, a_1,c_0,a_0\to E_0/2>0, and m has a positive lower bound. For v infinitesimal (including zero), FU22–FU23's six roots have strict residue

$$
m_c(T_v/\kappa-B_0)^2+\tfrac12 B_0^2+\tfrac12 U_v^2
     +\tfrac E2(P^2+N_y^2)+a_1Q^2.
$$

For finite positive-unit v, the real capped FU10–FU15 reduction is valid with these NU source bounds. If kappa_0<1, substitute m_0=kappa_0^2/[2(1-kappa_0^2)] and E_0=1/(1+kappa_0) into its full Schur polynomial; exact rational simplification gives

$$
\mathcal P(y,c)=E_0(m_0+\tfrac12)c(y-\kappa_0)^2
 +(m_0+\tfrac12)E_0^2y>0 \quad (y,c>0).
\tag{30}
$$

Thus the real gauge TU-kappa PQ=0 separates two strictly positive 2x2 Gram blocks (FU18), and the complete N_y row has residue E_0(1+v_0^2)/2>0. If kappa_0=1, m is unbounded but
$\mathcal P/m\to\frac12 c(y-1)^2+\frac14y>0;$ choose a sufficiently large finite real cap of the original m event using the real kernel-penalty argument, and restore its whole positive excess. This is an application of the FU full-parent identities plus new NU source estimates, not a circular application of the FU source theorem. Finally Bob infinity is NU's w=0 reciprocal local chart. CE covers every noninfinitesimal Alice test, and exact sign congruences, reciprocal homogenizations and real/imaginary state splitting finish the independent all-tests join.
Status and limits. Conditional only on the common strict bounded-Gram lemma and standard native source-entry assumptions that were separately audited through SD9–SD11 here, NU and EL both have explicit finite coefficient/monic-division tensor-SOS certificates for the complete original parent. Nothing here proves other native C1/C2 sectors, Scheiderer's universal tensor criterion, or a finite SDP matrix lift.

## Reproducibility

Run `python independent-review/checks/el-nu-endpoints/run_exact_checks.py` from the repository root (Python and SymPy). The included exact SymPy checkers verify the original SD9/SD10 balances, the SE2/SE3 sphere, the entire EL3/EL11/EL12 readback and EL13 modulo its exact SD contact jet, NU12/NU19 entire full-parent readback, the NU10/NU15 and EL7–EL11/EL18 helper identities, full CE2/CE6 signed-core identities, EB5 and EB6 complete-event amplitudes, FU8/FU9 complete finite-Bob parent, and the precise NU25/NU26 factorization. The checker does not replace the analytical boundedness, helper-purchase order, main-positive-residue, or source/stationarity arguments above.

The complete proof was recovered from the independent review and read against the native source. All eight symbolic modules were replayed locally on 26 October 2026. Their identities support the source and Gram arguments above; the other C1/C2 source domains remain separate.
