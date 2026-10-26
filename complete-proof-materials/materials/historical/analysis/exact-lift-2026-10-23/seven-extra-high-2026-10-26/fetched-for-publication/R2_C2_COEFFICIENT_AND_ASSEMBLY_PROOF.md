# Q222 native C2: R2 bounded-current and bounded-shape closure, with noncircular source ledger

**26 October 2026. Pinned source:** `anubhav70747/q222-finite-sdp`, commit `8960b3fead3ce6b185bb1387d148fb87d05d3d23`, `latex/certificates/c2_short.tex` (RER1–RER8, RSR2–RSR6, RFN4–RFN22, RL1–RL24, SH1–SH8, TU1–TU23, R2O1–R2O24), and `latex/certificates/c1_short.tex` (NC2–NC4, SC3–SC5, SC16–SC20, QI2–QI10). The independent UL/R2O theorem is an imported result precisely on unbounded relative current. **This note establishes no general Q222 spectrahedral-shadow theorem and no fixed matrix size.**

## 1. Quantified R2 statement and imported lemmas

Let the native **physical** C2 source lie in a real closed source field \(K_c\supset\mathbb R\), testing parameters and real state in an independently labelled real closed field \(K_t\), and the separated coefficient ring be \(K_c\otimes_{\mathbb R}K_t\). Assume the original complete stationary coefficients RER1, every strict source gap, \(d=1,\lambda=\beta=0,n\ge0\), and the R2 source domain
\[
 t_s\text{ positive unbounded},\qquad B_s,\ a_sB_s,\ t_sB_s\text{ positive infinitesimals}.
\tag{1}
\]
No positive lower bound is assumed on \(r,h,1-h,e=1-r,g_s\), or on a source/testing rate ratio. The quarter-parent \(P=F/4\) is the **whole** polynomial
\[
\begin{aligned}
P={}&m(vW-X+uvY-uZ)^2+nD_uD_vX^2+CD_v(uW-Y)^2+AD_uD_vY^2\\
&+(D+\tfrac12)(uW+uvX-Y-vZ)^2+(\alpha+\tfrac12)D_u(vY-Z)^2\\
&+\tfrac12D_v(X+uZ)^2-\tfrac12D_uD_vZ^2,
\end{aligned}\qquad D_u=1+u^2,\ D_v=1+v^2.
\tag{2}
\]
**R2 theorem (relative to explicit prior algebraic lemmas).** The full parent is a finite separated tensor SOS for every independent real/projective testing pair and arbitrary real or complex state under (1). The proof uses only: the native coefficient formulas and actual contact, the generic four-state polynomial identities RFN7–RFN22 and RSR12–RSR21; the explicitly displayed real PSD matrices RFN19 and SH1–SH8; the purely algebraic SC3–SC5, SC16–SC20 and QI2–QI10 coefficient/conductor identities, without asserting physicality of a formal negative-beta C1 tuple; the SJ3–SJ10 sign/endpoint identities; and the full-state R2O10–R2O24 outer square identities. The independently checked UL theorem enters **only** when \(a_s\) and the relative current are both unbounded. Each import is a polynomial identity or a finite bounded-Gram lemma on its stated positive-chart domain, not another all-C2 theorem.

### Source split (exhaustive, without overlap assumption)

Write
\[
 p=\frac{a_se}{1+ha_sB_s},\quad 0<p<1,
 \qquad\lambda_n=\frac{n}{t_sB_s}\ge0.
\tag{3}
\]
There are exactly four alternatives: (i) \(a_s=O(1)\), handled by RL; (ii) \(a_s\) unbounded and \(p\) infinitesimal, also RL; (iii) \(a_s,p\) unbounded/unit respectively, \(\lambda_n=O(1)\), handled by TU below; (iv) \(a_s\) unbounded, \(p\) a unit, \(\lambda_n\) unbounded, handled by the independent UL proof. If \(n=0\), it remains permitted in RL and TU; in particular the exact OP4 same-contact map identifies the physical C1 boundary. **Neither RL nor TU invokes C1 physical-source validity at a formal negative \(\beta\).**

## 2. A common native source bound before the two local constructions

In native weak coordinates, put \(z=a_sB_s\), \(y=t_sB_s\), \(R_s=1+hz\), \(g_s=R_s-a_se>0\). The exact stationary numerator identity is
\[
2r\Delta n=eQ-B_sK_s,\qquad K_s=(1+a_s)(h+t_s)(1+z),
\tag{4}
\]
\[
 Q=1+B_s+y+z(h+t_s+2r)+z^2\{t_s(1+er)+h(1+r)\}+B_sz(h+t_s).
\]
The inequality \(n\ge0\) gives \(eQ\ge B_sK_s\). Initially \(Q\le C(1+a_sy)\) because the only potentially large monomials are bounded by \(a_sy=a_st_sB_s\). The gap \(a_se<R_s\sim1\), multiplied into the physical inequality, forces \((1+a_s)a_sy=O(1)\), hence
\[
 a_sy\in\mathrm I_c,\quad a_s(1+a_s)y=O(1),\quad Q\sim1,\quad
 \frac{(1+a_s)y}{e}=O(1).
\tag{5}
\]
The exact coefficients then imply, writing \(S=1/t_s\), \(l=1-h\) and \(k=(1+a_sr)/(1+a_s)\in(1/2,1)\),
\[
 SC l\sim SAh\sim D/S\sim k/2,
\quad\frac{\alpha}{Sy^2(1+a_s)k^2}\sim\frac12,
\quad mg=\frac{a_s^2e(1+y)(1+zr)s_c}{2\Delta}=O(1).
\tag{6}
\]
These are *original-source*, not independently stipulated, bounds. The R2O source estimates and contact derivations depend on no additional rate relation.

## 3. RL: bounded native \(a_s\), or infinitesimal paired parameter \(p\)

Use RSR2 exactly:
\[
 j=\frac{a_s(1+a_sB_sr)}{R_s},\quad B=\frac{B_sR_s}{1+a_sB_sr},\quad
 T=\frac{t_s(1+a_sB_sr)}{R_s},\quad H=\frac{l}{t_s(1+a_sB_sr)},
\quad x=jB,\ y=TB.
\]
Write \(S=T^{-1}\), \(B=Sy\), \(H=S\ell\), \(x=jSy\), \(K=j+1-p\), \(\bar h=1-\ell\). The exact positive domain gives \(0<p<\min(1,j),0<\ell<1\); the physical estimates imply \(S,y,x,jy\in\mathrm I_c\). The compact RSR polynomials \(\mathcal P,\mathcal L,\mathcal N,\mathcal F,\mathcal K\) are those printed in RSR5/15/18, not free parameters. Direct termwise division (each nonleading term contains \(S,y,jy\), or \(p/j\le1\)) gives
\[
\frac{\mathcal P}{j+1}\sim1,\quad\frac{S\mathcal L}{K}\sim1,\quad
\frac{S\mathcal F}{F_b}\sim1,\quad
 F_b=j+1+\frac{p(1-p)}j,\quad \frac12<\frac K{j+1}<1.
\tag{7}
\]
There is no inversion of \(\ell,\bar h,j-p\). Crucially, the **exact finite** zero-\(p\) source numerator is
\[
\boxed{\mathcal N\big|_{p=0}=-yj(j+1)\{1+S(1-\ell)\}
 \left\{1+\frac{S\ell yj}{j+1}\right\}.}
\tag{8}
\]
Since \(\mathcal N-\mathcal N|_{p=0}\) is divisible by \(p\), with quotient \(1+\mathrm I_c\) after normalization, \(n=\mathcal N/[2(j-p)\mathcal P]\ge0\) forces \(\eta_2=yj(j+1)/p=O(1)\). If \(\eta_2\) is infinitesimal, \(n\sim p/[2(j-p)(j+1)]>0\). The source masses therefore satisfy
\[
\begin{split}
 &\mu=m(1-p)\sim\frac{jp}{2(j+1)},\quad h_0=\mu+\tfrac12\asymp1,\quad
 CH\sim SA\bar h\sim D/S\sim\frac{k}{2},\\
 &\alpha\sim\frac{Sy^2K^2}{2(j+1)},\quad
 f_c\sim\frac{j^2y^2F_b}{2(j+1)},\quad
 r_c\sim\frac{jy^2(1-p)K}{2(j+1)}.
\end{split}\tag{9}
\]
With \(\kappa=-k_c>0\), the exact RSR19 positive numerator implies \(\kappa/r_c=O(p+Sy)\) if \(p\) is infinitesimal, and otherwise (necessarily bounded-unit \(p,j\))
\[
 \frac\kappa{f_c}\sim\frac{pK}{jF_b}\le\frac12,
 \qquad jF_b-2pK=(j-p)K\ge0.\tag{10}
\]

### Entire state-normal identity and bounded masses

Normalize the *testing* coordinate \(w=(W+uY)/D_u\) to 1 only after the signed-core state-entry lemma, and leave \(z_t=(Z-uX)/D_u>0\) **arbitrary**. Put \(N_0=(X+uZ-pD_uv)/D_u\), \(M_0=(uW-Y-HD_uvz_t)/D_u\). The complete inverse is
\[
 X=N_0+pv-uz_t,\ Z=z_t+uN_0+puv,\quad
 W=1+u(M_0+Hvz_t),\ Y=u-M_0-Hvz_t.
\]
After retaining \(mpD_u^2N_0^2\), the full identity is **exactly**
\[
\begin{aligned}
G={}&\mu D_u[N_0^2-2vN_0+(1-p)v^2]+nD_v(N_0+pv-uz_t)^2\\
&+CD_uD_v(M_0+Hvz_t)^2+AD_v(u-M_0-Hvz_t)^2\\
&+(D+\tfrac12)D_u(M_0-(1-H)vz_t)^2\\
&+(\alpha+\tfrac12)[u(1-p)v-(1+Hv^2)z_t-uN_0-vM_0]^2\\
&+\tfrac12D_uD_v(N_0+pv)^2-\tfrac12D_v(z_t+uN_0+puv)^2.
\end{aligned}\tag{11}
\]
No normal has been suppressed. Its complete hard block is
\[
 H_N=h_0D_u+v^2/2+\alpha u^2+nD_v,\qquad
 H_M=D_u(CD_v+D+\tfrac12)+AD_v+(\alpha+\tfrac12)v^2,\quad
 H_{NM}=(\alpha+\tfrac12)uv.
\tag{12}
\]
The scalar is the *exact* RFN12 mass/circuit polynomial, and the two first normals are RFN21–RFN22. Their complete source balances can be verified without guessing a source: let \(\mu_N=2\mathcal P\mu\), \(\alpha_N=2E\mathcal P\alpha\), \(L_*=S\mathcal L\),
\(A_N=(1+x)(1+Hy)(K+Hjy)\), \(D_N=(1+Hy)(K+Hjy)(1-y[(1-H)j-p])\), \(E=1+Sy(\ell+p-1)\). Then the four **identically zero polynomial numerators** are
\[
\begin{aligned}
\mathsf S_0&=L_*-SD_N+H\mathcal P-A_N,\\
\mathsf S_1&=j(1+y)L_*-jSD_N-j(1-H)\mathcal P+\mu_N,\\
\mathsf S_2&=-\mu_NE Sy+(L_*-A_N)Ey
 +(\alpha_N+E\mathcal P)Sy(\ell-1+p)+\alpha_N,\\
\mathsf S_N&=-E\mathcal N-E\mu_N(1+x)+E\mathcal Pp
 -x(1-p)(\alpha_N+E\mathcal P)+j\alpha_N\\
 &\qquad+y[E\mathcal P(p-j)+jH(\alpha_N+E\mathcal P)-E\mathcal N].
\end{aligned}\tag{13}
\]
Their vanishing is checked by finite expansion of RSR5/6 in the exact checker. \(\mathsf S_0,\mathsf S_1,\mathsf S_2\) are precisely the three source identities required in RFN22; \(\mathsf S_N\) is the complete RFN21 normal identity. These equations give full normal equality, not just contact evaluation.

Set the harmonic current \(q_t=vz_t-u/S\), and \(M_0=L-Hq_t\) with \(L=\sqrt{S\ell\bar h}\,\widehat L\). A direct expansion of (11), with (13), gives the exact second-normal coefficient
\[
 C_q-HH_M=-D-H(C+D+\tfrac12)u^2;\tag{14}
\]
the first-normal row simplifies to RL22, with \(K_s=-\mu+(\alpha+\tfrac12)(\ell+p-1)\), \(\kappa_N=-(n+\tfrac12)(1-p/j)\). In particular the dangerous factor \(n(j-p)/j=O(p/[j(j+1)])\) cancels before a bounded-Gram coefficient is formed.

The shifted scalar retains a strictly positive \((\sqrt S q_t)^2\) mass with its **exact** coefficient
\[
D_q=S\{A\bar h^2+\ell^2(C+D+\tfrac12)D_u\}
 -\frac{f_c}{Tx}-TH^2Cu^2v^2-\frac{k_cv^2}{yTx}.
\tag{15}
\]
The first two terms have residue \(k/2\): \(SA\bar h\sim CH\sim k/2\), and \(\bar h+\ell=1\). The discarded/reallocated ratios are bounded with vanishing effects:
\[
 f_c/(Tx)=O(jy),\quad TH^2C=\ell CH=O(1),\quad
 -k_c/(yTx)=\kappa/(jy^2)=O(1).
\tag{16}
\]
This is a **full scalar** coefficient result: (11), RFN12 and the exact finite harmonic shift determine every remaining term, specifically the opposite row and its two cross terms RL15. The signed pair is purchased with SH2 and the exact SH1, \(a=3/2\), plus Young's square:
\[
\begin{aligned}
\kappa Z_*^2(V^2-1)^2\ \preceq_{\rm SOS}\
\frac74\kappa\Psi+\frac{21}{2}\kappa V^2Q_*^2+\frac72\kappa V^2B_*^2.
\end{aligned}\tag{17}
\]
Indeed \(Z_*(V^2-1)=-VQ_*-R_*\), \((R_*+VQ_*)^2\le\frac76R_*^2+7V^2Q_*^2\), and SH1 says \(R_*^2\preceq\frac32\Psi+3V^2(Q_*^2+B_*^2)\). In the positive-unit \(p\) branch, (10) bounds its circuit cost by \(7/8+o(1)\); in the infinitesimal \(p\) branch the direct \(r_c\) budget costs \(O(p+Sy)\). The other signed \(e_0\) cost is \(O(S(y+v^2))\). Thus positive **real** fractions of the strict scalar rows remain.

Both complete normals use the printed **finite explicit**, non-pointwise square identities RFN19 and SH1–SH8. The free rows are \(R_u=v(u^2-x)\), \(E_u=u(uvz_t-jy)\), \(E_v=v(uvz_t-jy)\), \(pv^2(v-u z_t/j)\), \(\sqrt S\mu uv(v-uz_t/j)\), and the whole SH6 opposite row. Each has its displayed **SOS budget**, with costs RL18–RL21, not merely a pointwise inequality. The correlated \(n\) completion must precede allocation:
\[
 h_0N_0^2+n(N_0+jb)^2=(n+h_0)(N_0+\tfrac{nj}{n+h_0}b)^2
 +\frac{nh_0j^2}{n+h_0}b^2,
\tag{18}
\]
\(b=v-uz_t/j\); when \(n=0\) the second row is literally absent and no \(n^{-1}\) is formed. The positive diagonal on \(\widehat N=\sqrt{n+h_0}(N_0+c_bb)\) has residue 1, on \(\widehat L\) and \(\sqrt S q_t\) residue \(k/2\), and the RFN19/SH strict scalar/helper pivots have positive real reserve. Equation (12) controls **all** hard cross entries, the remaining one being \((\alpha+\tfrac12)uv\sqrt{S\ell\bar h}/\sqrt{n+h_0}\in\mathrm I\). All normalized normal couplings have bounded source/testing arrays and vanish or are absorbed by fixed real helper pivots. This gives the full local bounded positive-residue Gram, with no source/testing rate condition.

## 4. TU: unbounded \(a_s\), unit \(p\), bounded relative current

Write
\[
S=t_s^{-1},\ z=a_s^{-1},\ q=a_s(t_sB_s),\ w=a_se,
\quad B_s=Szq,\quad e=zw,\quad \lambda_n=\frac{n}{qz}=O(1),\quad
 \tau=\frac q{zw}.
\tag{19}
\]
Here \(S,z,q\) are infinitesimal and \(0<w<R=1+Sqh\sim1\); in this branch \(p=w/R\) is a **positive unit**. The exact RER numerator identity divided by \(e\) is
\[
 Q_{\rm RER}=\tau\big[(1+z)(1+Sh)(1+Sq)+2r(z\Delta)\lambda_n\big].
\tag{20}
\]
Because \(Q_{\rm RER},z\Delta, r\) are units, the bracket is a positive unit, so
\[
 \tau\sim\frac1{1+2\overline\lambda_n}>0,
\qquad q=zw\tau,\quad n=\lambda_nqz,
\quad H_0\asymp q+z\asymp z.
\tag{21}
\]
In particular **neither \(q/z\) nor \(1/\tau\) diverges**; the opening of \(w\) is a unit in this branch.

For **algebraic comparison only**, use the formal same-contact C1 parameters \(k^2=Sq,\eta=1/S,\zeta=z,r_1=h,t_1=w,a_1=S(1-h)\). Their \(\beta^0\) is not presumed nonnegative. NC2–NC4 reconstruct their coefficient ratios exactly. Define the six actual-minus-formal differences using \(c_b=SC\), \(a_b=SA\), \(D_b=D/S\), \(\alpha_b=\alpha/(Sz)\), \(h_2=mz(R-w)\), \(h_b=zw(\beta+\tfrac12)\), and use \(\beta=0\) in the actual source. Set \(e=zw\) and
\[
\begin{aligned}
P={}&1+qz+Sq(2+z-2e)+Sq^2(z+e-e^2)+S^2hq^2(1+z-e),\\
P_C={}&1+qz+S\{1+q(2+2z-e)+q^2z\}\
 &+S^2\{q[h(1+z)+1-e]+q^2[hez-e(1+z)+2z+1]\}+S^3hq^2(1+z-e),\\
P_A={}&1+qz+S\{q(2+z-2e)+q^2z(1-e)\}+S^2q^2(1+z-e)(1+z-e-hz),\\
P_D={}&1+q(z+e)+S\{q[h(1+z)+1-e]+q^2(z+e-e^2)\}+S^2hq^2(1+z-e),\\
P_m={}&1+qz+S\{q(2+z-2e)+q^2z\}+S^2q^2\{(1-e)^2+z(1-e)+hez\},\\
P_b={}&(1+z)(1+qz)+Sq(1+z)(2+z-2e)+Sq^2z(1+z-e^2)+S^2q^2(1+z-e)^2.
\end{aligned}\tag{22}
\]
**Six exact factored numerator identities**, independently reconstructed from RER1 and NC2–NC4:
\[
\boxed{\begin{aligned}
(1-h)\delta c_b&=-\frac{n(1-e)P_C}{zP},&h\delta a_b&=-\frac{n(1-e)P_A}{zP},\\
\delta D_b&=-\frac{n(1-e)(1+qz)P_D}{zJP},&\delta h_2&=\frac{n(1-e)(1+qz)P_m}{zP},\\
\delta h_b&=\frac{n(1-e)P_b}{zP},&
\delta\alpha_b&=\frac{n(1-e)q^2(1+qz)[(1+S)e-S(1-h)(1+z)]}{zQP}.
\end{aligned}}\tag{23}
\]
Here \(J=1+Sq(1-e)-S(1-h)\), \(Q=1-Sqzh\) are positive units, and \(P\sim1\). Each identity is a cross-multiplied polynomial statement. The exact checker reconstructs RER1 and NC4, clears denominators, and verifies all six identically. Also \(P_b-(1+qz)P_m\) is divisible by \(z\) **after \(e=zw\)**, so
\[
 \delta(h_b-h_2)=n(1-e)\frac{P_b-(1+qz)P_m}{zP}
\tag{24}
\]
has no \(z^{-1}\) pole. Substituting \(n=\lambda_nqz\), and cancelling \(h\) or \(1-h\) **before** using (23), yields \(q\lambda_n O(1)\) on each of \((1-h)\delta c_b,h\delta a_b,\delta D_b,\delta h_2,\delta h_b\), \(\delta\alpha_b/q^2=\lambda_n qO(S+e)\), and \(\delta(h_b-h_2)/(qz)=\lambda_n O(1)\). Every \(O(1)\) is a rational polynomial array with denominators only \(P,J,Q,R\), never a hidden rate comparison.

### Complete additional current, both normal rows, hard matrix

Retain the entire original SC3 full-state polynomial (an algebraic identity irrespective of formal \(\beta^0\)), with actual coefficients; set \(\widehat N=\sqrt{H_0}N\), \(N_y=\sqrt\sigma(D_uY_b-\rho_0H)\), \(\sigma=c_b+a_b\), \(\rho_0=(1+Sq)h/R\), \(H_0=4h_2/R+4h_b/w\). Both opening cancellations are **exact**:
\[
1-\rho_0=(1-h)/R,
\qquad \delta(\sigma\rho_0-c_b)=\rho_0\delta a_b-(1-\rho_0)\delta c_b=q\lambda_nO(1).
\tag{25}
\]
The separate original \(n\)-packet is handled by a whole positive square, never expanded as an unbounded \(u/b\) remainder. With the four separated positive toric monomials \(A_*,B_*,C_*,D_*\) of TU15 and \(\Gamma^2=U_*B_*\), the exact full-normal identity is
\[
\boxed{\frac{4n}{S}D_uD_vX^2=
 4\lambda_nD_uD_v(U_*-B_*-czN)^2
 +4\lambda_nD_uD_v[2(1-e)\Gamma^2-(1-e^2)B_*^2-2(1-e)czB_*N].}\tag{26}
\]
This follows by expanding \((U_*-eB_*-czN)^2\); the identity \(cgX/b=czN+eB_*-U_*\) proves its equality to the original packet. The first square vanishes at the common contact. Its coefficient is literal zero when \(n=0\), and no division by \(n\) is used.

The **entire** scalar perturbation and the **two complete first-normal perturbations** are the coefficient-extraction identities TU21–TU23 (or (10)–(12) in the independently audited UL record), with \(\delta\)-arrays (23). Their exact normalized coefficient tests are
\[
\begin{gathered}
(a_1-c_1,b_1+2c_1,R^2P_L/q,4\alpha_bD_u/q^2)\big|_{\rm precurrent}
\longrightarrow(2,-8-8\lambda_0,2+4\lambda_0,2),\\
\text{additional current (26)}:\quad(0,+8\lambda_0,-4\lambda_0,0),
\qquad\lambda_0=\overline\lambda_n,
\end{gathered}\tag{27}
\]
so the combined full scalar has **exactly** the strict original residues \((2,-8,2,2)\). For example the differences of the two essential coefficient quotients at \(q,S,z=0\) follow immediately from (23): \(\delta(a_1-c_1)/q=0\), \(\delta(b_1+2c_1)/q=-8\lambda_0\), and \(\delta(R^2P_L)/q=4\lambda_0\). The bounded finite coquotients are obtained by monic numerator divisions in \(q,U=u^2,V=v^2\), after the factors in (25) are removed.

For the **first full normal**, the exact TU11 polynomial \(\mathsf F_N=[P_bD_v-(1+qz)P_m]/P-zD_uD_v\) has the finite form \(V-zUD_v+qF_{Nq}\), where \(F_{Nq}\) has unit denominators. In particular its normalized coefficients are bounded because \(H_0\asymp z\), \(c/\sqrt{H_0}=O(1)\), \(g/\sqrt{H_0}=O(1)\), and
\[
 (h_bD_v-h_2)/\sqrt{qH_0}=O(z+v^2).
\tag{28}
\]
The whole extra \(B_*\widehat N\) current has coefficient \(\lambda_ncz/\sqrt{H_0}=O(\lambda_n\sqrt q\,\sqrt z)\in\mathrm I\). For the **second full harmonic normal**, TU20 and the QI10 algebraic numerator allow a finite monic division \(R\mathsf F_H+\mathsf F_v=qT_q+UT_U+VT_V\) with bounded arrays, despite \(h,1-h\) having no lower bounds. All coefficients of the subsequent contact-root products have zero residue. Neither normal has been set to zero in obtaining these equations.

After subtracting the whole current square in (26) and the paired nonnegative current packets, the **full hard matrix** on \((\widehat N,N_y)\), in symmetric-matrix convention, is
\[
\boxed{\mathcal H=\begin{pmatrix}
1+\frac{z u^2v^2(4D_bS+2)}{H_0}
&-\frac{2\sqrt{Sz}\,uv(1+2D_bS)}{\sqrt{H_0\sigma}}\\
* &\displaystyle\frac{4D_v(c_b+a_b/D_u)+4D_bS^2+2S(1+v^2/D_u)+4\alpha_bzS^2v^2/D_u}{\sigma}
\end{pmatrix}.}\tag{29}
\]
This is obtained by differentiating the *complete* SC3 parent. The cross term contains both the \(D_bJ_b^2\) and the signed core \(2K_b\). In particular \(\bar {\mathcal H}=\operatorname{diag}(1,4)\) with all source/test infinitesimal rates independent. There are no additional hard-state directions.

The finite QI2–QI10/SC16–SC20 contact-root algorithm now applies to the *actual* bounded coquotient arrays: retain positive fractions of \(E_*^2,F_*^2,R_E^2,R_F^2\), and the full strict 16-root toric \(2T_4=2(\sum A_i^2-4\Gamma^2)\) representation. The scalar differences have zero source residue after (27), the complete first normals are zero-residue contact-root combinations by (28) and the harmonic monic division, and the hard block is (29). The contact value and three first jets vanish because (2) and every retained square share the native contact. Each remaining bare-\(\Gamma\) coefficient is monically divided at that contact, converted by SC18/20 to fixed toric roots, and bounded with a zero-residue coefficient. These are finitely many exact numerator extractions, not an approximation or an invocation of a C1 *physicality* theorem. The strict bounded-Gram lemma yields local tensor SOS even for arbitrarily small positive \(\lambda_n\). Literal \(n=0\) may alternatively use the actual physical C1 equality OP4.

## 5. Independent tests and projective endpoints

The source-independent SJ3 sign identities reduce finite-angle signs to positive angles, with a direct original-core/event SOS on the nonnegative-core branch; SJ9–SJ10 give literal homogeneous projective axes. On the negative-core branch, \(w=(W+uY)/D_u>0\) and \(z_t=(Z-uX)/D_u>0\); dividing by the nonzero pure-testing \(w\) is legal, and RL estimates never restrict \(z_t\). TU uses \(Z=1\) only on its actual state chart; at \(Z=0\), use the original signed-core square restriction. Complex states split into real/imaginary parts.

The exact full-parent reflection R2O4–R2O6 has \(G^TG=D_uD_vI\), sends the source to \((\lambda',n',m',C',A',D',\alpha',\beta')=(\alpha,m,n,C,D,A,0,0)\), and satisfies \(F_{c'}(v,u;G\psi)=D_uD_vF_c(u,v;\psi)\). Under (1), \(g_*=(CA)/(C+A)\) is unbounded and \(Sg_*\sim k/2\). Therefore R2O10 is strict when the **reflected Bob** angle is noninfinitesimal. When reflected Bob is infinitesimal and reflected Alice noninfinitesimal, R2O11 has a full four-state residue \(\operatorname{diag}(1,2,2,k_0^2-c_{u,0}^2)\).

**New simplification of the join.** If \(a_s\) is unbounded (both TU and UL, and RL's unbounded-\(a_s\) infinitesimal-\(p\) case), then
\[
 1-k=\frac{a_se}{1+a_s}<\frac{R_s}{1+a_s}\in\mathrm I_c,
 \quad\therefore\quad k_0=1.
\tag{30}
\]
For noninfinitesimal reflected Alice angle, \(c_{u,0}^2<1\), so **R2O10–R2O11 alone** cover every remote testing pair in these three subdomains. They require no finite caps.

For bounded \(a_s\), the only difficult remote case has positive-unit reflected Alice \(u\), infinitesimal reflected Bob \(v\), and \(k_0<1\). Then (R2O13–R2O24) force \(a_s,e\) positive units and use the original actual event coefficients, not a physical capped Bell inequality. The exact selected three-state matrix on \((X,Z,V)\) is
\[
\mathcal Q_3=m_c\{(D_u/u)V-X-uZ\}^2+n_cD_uX^2+\tfrac12D_uV^2-D_uVZ
 +\tfrac12(X+uZ)^2-2bV(uX-Z).
\tag{31}
\]
For \(p_0<1,j_0>p_0\), choose \(b=k_0/2\), \(m_0=j_0p_0/[2(1-p_0)(j_0+1)]\), \(n_0=p_0/[2(j_0-p_0)(j_0+1)]\): the last two Sylvester minors are
\[
\frac{p_0u_0^2D_u(j_0+1-p_0)}{4(j_0-p_0)(j_0+1)^2(1-p_0)},\qquad
\frac{p_0u_0^2D_u(j_0+1-p_0)^2}{8(j_0-p_0)(j_0+1)^3}>0.
\tag{32}
\]
At \(j_0=p_0<1\), use the large \(n_c\) event and the exact \((Z,V)\) determinant \(u_0^2/[4(j_0+1)^2]>0\). At \(p_0=1<j_0\), use the large \(m_c\) event; the \((X,V)\) determinant is
\[-\frac{D_u^2[2b(j_0-1)-j_0][2b(j_0+1)-j_0]}{4u_0^2(j_0^2-1)}>0\]
for \(j_0/[2(j_0+1)]<b<\min(k_0\sqrt{D_u}/2,j_0/[2(j_0-1)])\). At \(p_0=j_0=1\) choose \(1/4<b<\sqrt{D_u}/4\): the double event kernel has \(D_u(4b-1)V^2/(2u_0^2)>0\). All these are actual **real cap** positive-definiteness proofs, checked by exact symbolic minors; a positive-unit claim is not substituted for the requisite real lower bound. This finishes independent R2 tests without using a hypothetical global source theorem.

## 6. Noncircular C2 source partition: what this R2 theorem does and does not close

The exact involution
\((u,v;W,X,Y,Z;n,m,C,A,D,\alpha)\mapsto(1/v,1/u;Z,X,Y,W;n,D,\alpha,A,m,C)\)
preserves the whole parent after the homogeneous factor \(u^2v^2\). Together with the original event-cycle tests, it reduces the nine independent *source-angle* residue configurations \((\mathrm I,\mathrm U,\mathrm O)^2\) to the directed weak-II and the one-angle/strong patterns printed in the first 365 lines of the pinned C2 source. It is not permissible to replace an arbitrary strong-C2 contact with an R2 source: R2 specifically assumes \(t_s\) unbounded and \(B_s,a_sB_s,t_sB_s\) infinitesimal.

Within the **directed weak-II** chart, the source variables \(t_s\) infinitesimal/unit/unbounded give disjoint R1, middle, and R2 branches. The whole middle source domain, if native \(a_s\) is bounded, belongs to the original bounded-\(a_s\) theorem. If \(a_s\) is unbounded, the actual source gaps give \(e\) infinitesimal and the following mutually exhaustive leaf choices:

| Physical opening | Leaf | Required independent result |
|---|---|---|
| \(h\) unit, \(1-h\) infinitesimal | R3 | EN contact plus all-tests join |
| \(h,1-h\) units, \(j_{\rm gap}\) infinitesimal | R4 | EN contact plus all-tests join |
| \(h,1-h,j_{\rm gap}\) units, \(g_s\) infinitesimal | R6 | RB/RI, both \(a_s^2B_s\) scales |
| \(h,g_s\) infinitesimal | R7 | RU, both \(a_s^2B_s\) scales |
| \(h,j_{\rm gap}\) infinitesimal, \(g_s\) unit | **R8** | **Independent joint-D gate, not implied by R7 or R10** |
| \(h\) infinitesimal, \(g_s,j_{\rm gap}\) units | R10 | Weighted/thin/canonical-zero alternatives |

The source ledger also requires the separately checked R1, bounded-\(a_s\), and now **complete R2** branches; literal \(n=0\) is physical C1 via OP4. The original UU/five-event entry may invoke the real five-event lemma only when \(m,C,A,D,\alpha\) each have a **positive real lower bound**; original 'strong' plus the source estimates supplies this, mere positivity would not.

**R8 is provably nonempty** and therefore cannot be omitted from a global C2 source-coverage implication. Over a positive infinitesimal \(s\), choose
\[
 a_s=s^{-1},\ B_s=s^2/4,\ h=s,\ r=1-s/2,\ e=s/2,\quad
 t_s=\frac{1-s}{1+(s/4)r-s}.
\tag{33}
\]
Then \(z=a_sB_s=s/4\), the exact original \(j_{\rm gap}=1+zr-(1-h)/t_s=s\) is infinitesimal, \(g_s=1-a_se+hz=1/2+s^2/4\) is a unit, and \(t_s\to1\). The physical \(n\)-numerator (4) satisfies
\[
 \lim_{s\downarrow0}(eQ-B_sK_s)/s=1/4>0.
\tag{34}
\]
All remaining native gap inequalities hold for sufficiently small positive \(s\). Thus this is an actual positive-native R8 contact with \(n>0\), disjoint from both R7 and R10. The proposition **'R7 plus R10 covers all infinitesimal-\(h\) sources' is false**. R8 needs its own checked whole-state bounded coefficient proof (or an independently proved theorem matching precisely its source domain), not a recitation of section titles.

The separate one-angle IU (bounded/unbounded native \(b\) and their \(g,r\) openings), strong IH/HI, radial two-weak complementary-event and directed weak OO reflected-II leaves also require their **own stated entire-parent and all-testing theorems**. The source table is exhaustive algebraically; the checked external list in the current assignment does **not**, by itself, certify every one of these leaves. The correct noncircular assembly implication is therefore: **if** each listed original domain has an independently proved bounded entire-state tensor-SOS lemma, including R8, these separate one-angle/strong and two-weak joins, the literal omitted-event faces, and the positively-capped UU case, **then** the exact source partition and reflection yield the all-tests C2 tensor-SOS assertion. This R2 note verifies the R2 premise of that implication, not its other premises, and does not assert qualitative/extracted global Q222 lift existence.

## 7. Exact reproductions and provenance

The accompanying symbolic scripts check, without solver calls or sampling: (i) the six **original** TU9 C2-versus-NC4 C1 rational source identities and TU10 divisibility; (ii) the generic full four-state RFN7 parent, both hard diagonals and the hard cross entry; (iii) RL compact RSR5 physical \(\mathcal N|_{p=0}\) numerator and **all four** polynomial native source balances RFN21/RFN22 after clearing unit denominators; (iv) SH1/SH2, RL17 and TU17 signed-square expansions, plus SH4; (v) R2O21–R2O24 exact finite-cap minors; and (vi) the positive Puiseux R8 branch (33)–(34). Positive-residue and source-boundedness conclusions in the text are proofs from the displayed ratio estimates and positive square budgets, not claims inferred from the number of checks. The imported SC17–SC20/QI finite conductor is used as the pinned explicit polynomial lemma, not as a C1 physical-source result.
