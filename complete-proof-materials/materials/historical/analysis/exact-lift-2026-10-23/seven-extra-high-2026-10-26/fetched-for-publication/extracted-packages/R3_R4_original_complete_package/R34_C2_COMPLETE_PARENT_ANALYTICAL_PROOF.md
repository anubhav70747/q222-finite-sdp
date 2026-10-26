# Native C2 R3/R4: full separated SOS, all independent tests

**Source:** `latex/certificates/c2_short.tex`, commit `171a7e94b8f52f2bb601562eda362badfb687a9a`, RSR2–RSR20, RFN4–RFN24, EN1–EN25. **Scope:** only physical C2 R3/R4, not R1, R10, all C2 or global Q222.

## Theorem

Let \(K_c,K_t\supset\mathbb R\) be independently labelled real-closed extensions. Use the **actual** stationary C2 source with \(d=1,\lambda_{\rm event}=\beta=0\), \(m,C,A,D,\alpha>0,n\ge0\), all its strict native gaps and an actual noncommuting contact. Assume its two contact half-angles are infinitesimal; its total source scale is unbounded; and the original weak-II source variables \(a_s,t_s,h\) satisfy \(a_s\) unbounded, \(t_s,h\) positive bounded units and \(e=1-r\) positive infinitesimal. Assume either

- **R3:** \(1-h\) infinitesimal; or
- **R4:** \(1-h\) a positive bounded unit and \(1+a_sB_sr-(1-h)/t_s\) infinitesimal.

Then the complete four-state parent \(F\) is a sum of finitely many squares in \(K_c\otimes_{\mathbb R}K_t\) for every independent projective testing pair and every homogeneous real or complex state, including all independent source/testing rates, source \(n=0\), axes, signs and projective infinities. The proof below is constructive **per source/test pair**; no fixed feature degree or source-independent Gram is claimed.

The only general algebraic input is the **strict bounded-tensor-Gram lemma**: a finite symmetric matrix over \(\mathcal B=B_c\otimes_{\mathbb R}B_t\) with strictly positive-definite real residue induces a separated tensor SOS on arbitrary (possibly unbounded) root columns. It follows by real congruence and the exact diagonal/\((e_i\pm e_j)\) decomposition with positive-residue scalar coefficients; the scalar positive-residue tensor lemma expands each coefficient as a sum of pure-factor squares. No mixed inverses are used.

## 1. Exact native source chart and complete bounds

Use \(c_s=1+a_sB_sr\), \(R_s=1+ha_sB_s\), and put
\[
 j=\frac{a_sc_s}{R_s},\quad B=\frac{B_sR_s}{c_s},\quad
 T=\frac{t_sc_s}{R_s},\quad p=\frac{a_se}{R_s},\quad
 H=\frac{1-h}{t_sc_s},\quad x=jB,\ y=TB,\ q_z=jT.
\tag{1}
\]
The *actual* source inverse is
\[
 r=\frac{j-p}{j+px},\quad h=\frac{1-HT}{1+HTx},\quad
 a_s=\frac{j+px}{1+HTx},\quad
 B_s=\frac{B(1+HTx)}{1+pB},\quad
 t_s=\frac{T(1+pB)}{1+HTx}.
\tag{2}
\]
In the present regimes \(j\) is unbounded, \(B,x,y,e,p/j\) infinitesimal, \(T\) and \(1-HT\) are positive units; \(0<p<1\), \(0<H<1\). In R3, \(H\) is infinitesimal and \(\bar\sigma=1+\bar T\); in R4, \(1-H\) is infinitesimal, \(0<\bar T<1\), and \(\bar\sigma=1\), where \(\sigma=1+(1-H)T\).

The exact original polynomials \(\mathcal P,\mathcal L,\mathcal N\) and six event weights are RSR5–RSR6. To reproduce the scale calculation without imposing a relation between infinitesimals, substitute **formally** \(j=d^{-1}, B=dx\) and reduce the canceled rational arrays at \(d=x=0\). This gives, exactly in the bounded source coefficient ring,
\[
 \frac{\mathcal P}{j+1}\sim1,\quad
 \frac{\mathcal L}{j+1}\sim\sigma,\quad
 \frac{\mathcal F}{j+1}\sim\sigma,\quad
 \frac\mu p\sim\tfrac12,\quad A\sim\frac{T}{2(1-HT)},\quad
 CH\sim\frac{\sigma}{2T},\quad D(1-H)\sim\frac1{2T},\quad
 \frac\alpha{B^2(j+1)}\sim\frac\sigma2,
\tag{3}
\]
where \(\mu=m(1-p)\), \(\mathcal F\) is the finite RSR15 polynomial. Here \(\sim\) means the ratio has real residue 1, with the displayed bounded positive expression used as denominator. All displayed denominators are units. In particular \(m\) itself may diverge at \(1-p\to0\) and is *not* inserted into a bounded Gram; \(\mu\) is bounded.

The exact RSR8 numerator grouping is
\[
 \mathcal N=p(1+\epsilon_p)-Bj(j+1)\sigma(1+\epsilon_B),
 \qquad \epsilon_p,\epsilon_B\in\operatorname{Inf}(K_c),
\]
with finite source coquotients from RSR5 (the ungrouped remainder is bounded by \(C BT\,Bj(j+1)\)). Since \(n=\mathcal N/[2(j-p)\mathcal P]\ge0\),
\[
 \eta_*:=\frac{Bj(j+1)}p=O(1),\qquad
 n=O(p/j^2)\in\operatorname{Inf}(K_c).
\tag{4}
\]
If \(\eta_*\) is infinitesimal, then \(n\sim p/(2j^2)>0\). If \(n=0\), necessarily \(\eta_*\) is a positive unit. The correlated current
\[
 h_0=\mu+\tfrac12,\qquad
 c_b=\frac{nj}{n+h_0},\qquad
 M_e=\frac{nh_0j^2}{n+h_0}
\tag{5}
\]
obeys \(h_0\asymp1\) and, when \(n>0\),
\[
 \frac{c_b^2}{M_e}=\frac{n}{h_0(n+h_0)}\in\operatorname{Inf};\qquad
 \eta_*\in\operatorname{Inf}\Longrightarrow M_e\sim p/2,
 \ \frac{p^2}{M_e}=O(1).
\tag{6}
\]
The canceled residues of (3) give \((C,A,D)=(\infty,T/2,1/(2T))\) in R3 and \((1/(2T),T/[2(1-T)],\infty)\) in R4. The *positive-unit* finite residues, not merely positivity, authorize the real event caps used later.

## 2. Exact complete parent and the two full normal rows

Let \(D_u=1+u^2,D_v=1+v^2\). The unmodified quarter-parent is
\[
 \begin{aligned}
 P={}&m(vW-X+uvY-uZ)^2+nD_uD_vX^2+CD_v(uW-Y)^2+AD_uD_vY^2\\
 &+(D+\tfrac12)(uW+uvX-Y-vZ)^2
 +(\alpha+\tfrac12)D_u(vY-Z)^2
 +\tfrac12D_v(X+uZ)^2-\tfrac12D_uD_vZ^2.
 \end{aligned}\tag{7}
\]
Literal axes have direct original-parent square identities. In the only remaining double-infinitesimal negative-core chart, normalize the *pure-testing* rotated coordinate \((W+uY)/D_u=1\) and let \(z=(Z-uX)/D_u>0\), arbitrary and not presumed bounded. Put
\[
 X=N_0+pv-uz,\quad Z=z+uN_0+puv,\quad
 W=1+u(M_0+Hvz),\quad Y=u-M_0-Hvz.
\tag{8}
\]
Retain the exact positive source square \(mpD_u^2N_0^2\) and divide the remainder by the positive **pure-testing** \(D_u\). Then \(G=(P-mpD_u^2N_0^2)/D_u\) is the exact entire polynomial
\[
\begin{aligned}
G={}&\mu D_u[N_0^2-2vN_0+(1-p)v^2]+nD_v(N_0+pv-uz)^2\\
&+CD_uD_v(M_0+Hvz)^2+AD_v(u-M_0-Hvz)^2\\
&+(D+\tfrac12)D_u(M_0-(1-H)vz)^2\\
&+(\alpha+\tfrac12)[u(1-p)v-(1+Hv^2)z-uN_0-vM_0]^2\\
&+\tfrac12D_uD_v(N_0+pv)^2-\tfrac12D_v(z+uN_0+puv)^2.
\end{aligned}\tag{9}
\]
In particular its hard block, before any state-normal completion, is exactly
\[
 H_N=h_0D_u+\tfrac12v^2+\alpha u^2+nD_v,\quad
 H_M=D_u(CD_v+D+\tfrac12)+AD_v+(\alpha+\tfrac12)v^2,
 \quad H_{NM}=(\alpha+\tfrac12)uv,
\tag{10}
\]
where the mixed term in \(G\) is \(2H_{NM}N_0M_0\).

Put
\[
 b_0=v-uz/j,\quad q_0=vz-Tu,\quad R_0=z-uv/B,\quad R_u=v(u^2-x).
\]
Define the bounded source combinations
\[
\begin{aligned}
K_c&=p/2+j[H(\alpha+\tfrac12)-\tfrac12]+n(p-j),\\
K_R&=-\mu-(\alpha+\tfrac12)(1-p)+(T/j)K_c,\\
K_d&=CH-(D+\tfrac12)(1-H),\\
C_q&=K_d+AH+\tfrac12+v^2[CH(1+u^2)+AH+(\alpha+\tfrac12)H].
\end{aligned}\tag{11}
\]
The **complete** half-derivatives of (9), on \(N_0=M_0=0\), are exactly
\[
\boxed{\begin{aligned}
L_N={}&[j(n-\alpha)+p(n+\tfrac12)v^2]b_0+K_RR_u+(K_c/j)uvq_0,\\
L_M={}&C_qq_0+\mu uvb_0+v(\alpha-BTCHu^2)R_0.
\end{aligned}}\tag{12}
\]
Indeed direct expansion of (9) minus (12) equals respectively \(vS_N\) and \(uS_0+u^2vzS_1/j+uv^2(S_2-S_1)\), with
\[
\begin{aligned}
S_N&=np-\mu+p/2+[-\mu-(\alpha+\tfrac12)(1-p)]x+j(\alpha-n)+yK_c,\\
S_0&=T(K_d+AH+\tfrac12)-A,\\
S_1&=jK_d+\mu+xTCH,\\
S_2&=jK_d+xTCH+T[CH+AH+(\alpha+\tfrac12)H]
 -A-(\alpha+\tfrac12)(1-p)+\alpha/B.
\end{aligned}\tag{13}
\]
All four vanish identically on the *actual rational source RSR6*. The checker clears their source denominators and verifies four zero polynomial numerators. In particular
\[
 K_d=-\mu/j-BTCH,
 \qquad C_q=A/T+v^2[CH(1+u^2)+AH+(\alpha+\tfrac12)H],
\tag{14}
\]
so \(K_R,K_c/j,C_q\) are bounded even at the small \(H\), \(1-H\), \(1-p\) openings; no term requires a mixed inverse or a separate positive bound on \(m\).

## 3. Exact scalar circuit and signed-mass spending

Set \(U=u/\sqrt x,V=v/\sqrt y,Z_*=z/\sqrt{jT}\), and
\[
 Q_*=U-VZ_*=-q_0/(T\sqrt x),\quad B_*=V-UZ_*=b_0/\sqrt y,\quad R_*=Z_*-UV=R_0/\sqrt{jT},
\]
\[
\Psi=U^2+V^2+Z_*^2+U^2V^2Z_*^2-4UVZ_*
 =(U-V)^2+Z_*^2(1-UV)^2+2UV(1-Z_*)^2.
\tag{15}
\]
The last term is a separated square because \(\sqrt{UV}=(xy)^{-1/4}\sqrt{uv}\). With the ten explicit monomial masses RFN12–RFN13 (computed directly from (9)), the *entire* scalar is exactly
\[
\begin{aligned}
G|_{N_0=M_0=0}={}&q_cQ_*^2+nj^2b_0^2+r_cR_*^2+d_cV^2Q_*^2
 +(np^2y^2+e_0)V^2B_*^2\\
 &+f_c\Psi+g_cV^2\Psi+k_cZ_*^2(V^2-1)^2,
 \quad e_0=p^2y^2/2-g_c.
\end{aligned}\tag{16}
\]
An exact coefficient check substitutes the RSR6 weights in both (9) and (16) and verifies every monomial in \((u,v,z)\), not only the contact restriction.

Write \(\kappa=-k_c>0\). The source polynomials RSR15–RSR19 give
\[
\begin{gathered}
f_c\sim Tx^2\sigma/2>0,\quad
q_c\sim Tx/[2(1-HT)]>0,\quad
r_c\sim TjB^2(1-p)\sigma/2>0,\\
g_c/f_c=Hy\,\mathcal L/\mathcal F\in\operatorname{Inf},\quad
 d_c=M_{b^2c}-k_c>0,\\
\frac\kappa{f_c}=O((p+B)/j)\in\operatorname{Inf},\quad
\frac\kappa{q_cy}=O(p+B),\quad
\frac\kappa{r_c}=O((p+B)/(1-p)).
\end{gathered}\tag{17}
\]
All source coquotients follow by cancelling \(x=jB\) **before** specialization; (3) makes every remaining denominator a unit. There is no bound involving a testing rate.

The signed \(e_0\) packet is paid **exactly** by
\[
 g_cV^2\Psi+e_0V^2B_*^2
 =g_c[(1+V^2)\Psi-V^2B_*^2]-g_c\Psi
       +\tfrac12p^2y^2V^2B_*^2.
\tag{18}
\]
The bracket is a separated SOS by the exact RFN18b identity, and \(g_c/f_c\) is infinitesimal. For the other signed packet use
\[
 Z_*(V^2-1)=-VQ_*-R_*,\quad
 (VQ_*+R_*)^2\le_{\rm SOS}2V^2Q_*^2+2R_*^2.
\]
Purchase \(2\kappa R_*^2\) in the following exhaustive cases:

| Source case | Exact SOS budget | Relative costs against the retained \(f_c\Psi,q_cQ_*^2,M_eb_0^2,r_cR_*^2\) |
|---|---|---|
| \(p\) infinitesimal | direct \(r_cR_*^2\) | \(O(p+B)\) |
| \(p\) unit, \(\eta_*\) unit | RFN18a: \(R_*^2\preceq(1+U^2)\Psi\) | \(O((p+B)/j)+O(u^2)\), since \(j^2B\asymp\eta_*p\) is a unit |
| \(p\) unit, \(\eta_*\) infinitesimal | RFN18c: \(R_*^2\preceq2(\Psi+V^2Q_*^2+U^2B_*^2)\) | \(O((p+B)/j),\ O(v^2),\ O(u^2)\), since \(M_e\sim p/2\) is a unit |

The other \(2\kappa V^2Q_*^2\) costs \(O((p+B)v^2)\) relative to \(q_cQ_*^2\). Every quoted cost is in the maximal ideal of the bounded tensor ring for **every independent** \(u,v\) infinitesimal. Thus fixed positive real fractions of \(q_cQ_*^2,r_cR_*^2,f_c\Psi\), and \(M_eb_0^2\) when present survive. The other positive packets \(d_cV^2Q_*^2\), \(\frac12p^2y^2V^2B_*^2\) and the explicit SOS brackets remain whole.

## 4. Free rows and a finite Gram algorithm for **both normals**

Complete the *whole* correlated \(n\)-current with the exact source identity
\[
 h_0N_0^2+n(N_0+jb_0)^2
 =(n+h_0)(N_0+c_bb_0)^2+M_eb_0^2.
\tag{19}
\]
Keep \(nv^2(N_0+pb_0)^2\), \(h_0(uN_0)^2\), \(\frac12(vN_0)^2\), \(\alpha(uN_0)^2\) and the original \(mpD_u^2N_0^2\) as **whole nonnegative squares**, never expanding them into an uncontrolled bounded matrix. At \(n=0\), set \(c_b=M_e=0\) and omit the effective-current root. Normalize
\[
 \widehat N=\sqrt{n+h_0}(N_0+c_bb_0),\qquad
 \widehat M=M_0/\sqrt{TH(1-H)}.
\tag{20}
\]
The coefficient of \(\widehat N^2\) is 1. From (3), (10) and (14), its other diagonal obeys
\[
 d_M=TH(1-H)H_M\in\mathcal B,\qquad
 \overline{d_M}=\begin{cases}(1+\overline T)/2,&\mathrm{R3},\\1/2,&\mathrm{R4}.
 \end{cases}\tag{21}
\]
This is the exact complete hard diagonal, not a determinant test. The mixed \(\widehat N\widehat M\) coefficient is infinitesimal because it contains \(uv\sqrt{TH(1-H)}\). The \(q_0/\sqrt T\) scalar pivot is \(q_c/(Tx)\sim[2(1-HT)]^{-1}\), a positive unit. Its full couplings to the normals are \((K_c/j)uv\sqrt T\) and \(T\sqrt{H(1-H)}C_q\), both infinitesimal. The shift \(c_b/\sqrt{M_e}\) is infinitesimal by (6), when defined.

We now give a **finite row-purchasing and matrix-extraction algorithm**; it is often shorter than writing the many Gram entries. A row \(h\) is called *free* if an exact SOS identity gives \(h^2\preceq_{\rm SOS}\sum a_iS_i^2\), where each \(a_i\) divided by the appropriate surviving positive scalar mass has zero residue. Any prescribed **finite real** coefficient \(L>0\) on \(h^2\) can then be purchased, with the old scalar diagonals changed only infinitesimally, by the exact identity
\[
 \sum M_iS_i^2=Lh^2+\sum(M_i-La_i)S_i^2
        +L\big(\sum a_iS_i^2-h^2\big).
\tag{22}
\]
The last summand is an explicit separated SOS, not a positivity inequality in a single field.

The following five rows are free, with the **actual** source bounds:

| Free row \(h\) | Exact spending proof and cost |
|---|---|
| \(R_u=v(u^2-x)\) | RFN19 (7- and 3-row exact positive matrices), or RFN20. For \(1-p\) a unit, relative costs \(O(B),O(Bu^2),O(u^2/(1-p))\). For \(1-p\) infinitesimal and \(\eta_*\) a unit, costs \(O(B+u^2/j+u^4/(j^2B)),O(Bu^2)\), and \(j^2B\) is a unit. For both \(1-p,\eta_*\) infinitesimal, costs \(O(u^2/j),O(u^4/j),[8u^4+2(u^2-x)^2]/M_e\); here \(M_e\sim p/2\sim1/2\) is a unit. |
| \(j\alpha b_0\) | RFN18a: \(B_*^2\preceq(1+U^2)\Psi\). Its costs relative to \(f_c\Psi\) are \(O(B^3(j+1)^2)+O(B^2(j+1)^2u^2/j)\), both infinitesimal because \(B(j+1)=B+x\) is infinitesimal. This works for all R3/R4 openings, with no split at \(p=1\). |
| \(pv^2b_0\) | If \(\eta_*\) infinitesimal, charge \(M_eb_0^2\) at \(O(v^4p^2/M_e)=O(v^4)\). Otherwise use RFN18b \(V^2B_*^2\preceq(1+V^2)\Psi\) at cost \(O(v^2(y+v^2))\). |
| \(\sqrt T\,\mu uvb_0\) | If \(\eta_*\) infinitesimal, charge \(M_eb_0^2\) at \(O(u^2v^2p^2/M_e)\). Otherwise use RFN18b at \(O(Tu^2(y+v^2))\). |
| \(v(\alpha-BTCHu^2)R_0\) | Put \(J_s=B(j+1)+u^2\). The exact bound \(\alpha-BTCHu^2=O(BJ_s)\) gives squared M-normal cost \(a_R=T(jT)v^2(\alpha-BTCHu^2)^2\). Its ratios to \(f_c,q_cy,r_c,M_eyx\) are \(O((T/j)v^2J_s^2),O(v^2J_s^2),O(Tv^2J_s^2/(1-p)),O(T(1+1/n)v^2J_s^2/j^2)\). Use the direct \(r_c\) mass when \(1-p\) is a unit, RFN18a when \(1-p\) infinitesimal and \(\eta_*\) a unit (then \(j^2B\) is a unit), and RFN18c when both are infinitesimal (then \(M_e\sim1/2\)). Each actual relative cost has an extra infinitesimal testing factor. |

The exact RFN19 identity, with \(D_*=V(U^2-1)\), is
\[
 4\Psi+4U^2Q_*^2+4U^2R_*^2-D_*^2
 =t_7^TG_7t_7+4VZ_*(1-U)^2+UZ_*t_3^TG_3t_3,
\tag{23}
\]
where
\[
\begin{split}
t_7={}&(Z_*-U^2V,V-U^2V,U-U^2V,UZ_*-U^2V,UV-U^2V,UVZ_*-U^2V,U^2-U^2V),\\
t_3={}&(1-UV,V-UV,U-UV),\\
G_7={}&\begin{pmatrix}
4&-2&-3/2&0&-1/2&0&1/2\\
-2&3&0&-3/2&0&-1&2\\
-3/2&0&4&3/2&-2&-3&0\\
0&-3/2&3/2&4&-3/2&0&-3/2\\
-1/2&0&-2&-3/2&3&2&-1\\
0&-1&-3&0&2&8&-4\\
1/2&2&0&-3/2&-1&-4&4
\end{pmatrix},\quad
G_3=\begin{pmatrix}3&-2&-2\\-2&2&2\\-2&2&3\end{pmatrix}.
\end{split}
\]
Its leading minors are respectively \(4,8,101/4,1049/16,907/16,3279/16,743/16\), and \(3,2,2\). RFN20 is the exact identity
\[
\begin{split}
&8U^2\Psi+8U^4Q_*^2+[8U^4+2(U^2-1)^2]B_*^2-[V(U^2-1)]^2\\
&=2U^2\{4\Psi+4U^2Q_*^2+4U^2B_*^2-[Z_*(U^2-1)]^2\}
  +\{UZ_*(U^2-1)-(U^2-1)B_*\}^2.
\end{split}\tag{24}
\]
The braced expression is (23) after swapping \(V\) and \(Z_*\) and is SOS. The checker verifies both displayed identities exactly. RFN18a–c are likewise literal finite SOS polynomial identities (including the five nonnegative terms in RFN18c), not estimates.

**Matrix recipe.** Begin with (9). Substitute (12), the exact scalar (16), (18), (19), and (20). Pay the \(k_c\) term by the three-case procedure in §3. Buy a positive real pivot on *each* of the five rows in the table using (22) and the displayed exact SOS budgets. Represent \(\Psi\) by its three roots in (15). Keep the strictly positive normalized rows
\[
 \widehat N,\ \widehat M,\ q_0/\sqrt T,\ \sqrt{r_c}R_*,\
 \sqrt{f_c}(U-V),\ \sqrt{f_c}Z_*(1-UV),\
 \sqrt{2f_cUV}(1-Z_*),\ \sqrt{M_e}b_0\ (n>0),
\tag{25}
\]
and the five purchased helpers. Collect each remaining term \(a_ir_i^2+b_{ij}r_ir_j\) with the finite rule \(H_{ii}\mathrel{+}=a_i,H_{ij}=H_{ji}\mathrel{+}=b_{ij}/2\); all other terms are the explicitly listed retained nonnegative squares. This is a finite coefficient-extraction algorithm producing **every entry** of the full symmetric matrix; its input contains no higher normal terms, by (9). All coefficient arrays are bounded, as proved above, and every division is by a pure-source or pure-testing positive unit or a pure-source positive scale already absorbed in a row. Neither \(z\), \(U\), \(V\), \(R_u\) nor a test/source difference is inverted in a matrix coefficient.

All off-diagonal residues except the bounded couplings of \(\widehat N\) to purchased helpers vanish: \((K_c/j)uv\sqrt T\to0\), \(T\sqrt{H(1-H)}C_q\to0\), the complete \(\widehat N\widehat M\) coupling vanishes, and the current shift is infinitesimal by (6). Source \(K_R\) is bounded; indeed \(\bar K_R=-(1+\bar T)/2\) in R3 and \(\bar K_R=-1/2\) in R4. Choose the five finite **real** helper pivots successively large by the elementary real Schur complement. This makes the entire real residue Gram strictly positive definite, while (22) spends only an infinitesimal fraction of the old scalar masses. The bounded-tensor-Gram lemma gives the exact separated SOS for the full \(G\), hence for \(P\) after restoring the original \(mpD_u^2N_0^2\) and the pure-testing factor \(D_u\). The original \(F=4P\) follows by the real square 2. This completes the double-infinitesimal local contact chart **including both unrestricted normals and the weak-\(m\) limit**.

## 5. Exact independent testing/state join

The signed CHSH core has a pure-testing Schur decomposition. At \(Z=0\), it is a nonnegative sum of squares. When \(Z\ne0\), scale \(Z=1\) by a pure-testing factor. For \(u,v>0\), its Schur complement in the rotated coordinate \(w_t=(W+uY)/D_u\) is a positive pure-testing multiple of \(D_uuv w_t^2-2D_uD_vw_t+uv\), which is nonnegative when \(w_t\le0\). Therefore a negative-core evaluation necessarily has \(w_t>0\) and can be normalized to \(w_t=1\) by another pure-testing factor; \(z\) remains arbitrary. This justifies **every state chart** used in §2, without constraining the state ratio. Literal axes have the direct original-parent square forms. The sign congruences of the full parent permit \(u,v\ge0\).

For every other testing residue, use the *actual full parent*, cap the unique unbounded source coefficient \(C\) (R3) or \(D\) (R4) at a sufficiently large positive real, and retain its excess as a whole event square. Keep \(n,\alpha\) whole and replace \(m\) only by \(\widetilde m=m/(1+m)\le m\), retaining the excess. The bounded finite nonzero source weights are given exactly in (3): \(A_0=T_0/2,D_0=1/(2T_0)\) in R3 and \(C_0=1/(2T_0),A_0=T_0/[2(1-T_0)]\) in R4.

- **Finite Bob with nonzero residue, bounded Alice:** use full rows \((uW,Y,X,Z,\sqrt{\widetilde m}W)\); the unbounded-row kernel restriction has the two explicit three-state forms EN22. Their three leading principal minors in R3 are
  \[
  \frac{D_u(TD_v+v^2)}2,\quad
  \frac{D_u(TD_v+v^2)(Tu^2v^2+Tv^2+T+u^2v^2)}{4T},\quad
  \frac{v^2(1+T)D_u^2D_v(Tu^2+v^2)}{8T};
  \]
  in R4 they are
  \[
  \frac{D_u(T+v^2)}{2(1-T)},\quad
  \frac{D_uD_v(T+v^2)(T+u^2v^2)}{4T(1-T)},\quad
  \frac{v^2D_u^2D_v^2(Tu^2+v^2)}{8T(1-T)}.
  \]
  All are strictly positive because \(v\) has a nonzero real residue. A real penalty cap on the strong row gives a strictly positive full five-row residue, including \(m\to0\).
- **Infinitesimal Bob, Alice a bounded unit:** use \((uW-Y,X+uZ,vZ,Y)\). The exact complete residue is EN23. On the R3 strong-row kernel its \(2\times2\) determinant, after the positive \(D_u\) factor is removed, is \((1+T)u^2/4>0\); on the R4 strong-row kernel it is \(u^2/[4(1-T)]>0\). The positive normal pivot and a real cap make the entire form strict. Literal \(v=0\) is its polynomial value.
- **Bob infinite, Alice bounded:** the exact homogeneous five-root EN24 residue \([S+\sqrt{\widetilde m}uY]^2+c_0(uW-Y)^2+(a_0+1/2)D_uY^2+d_0(uX-Z)^2+D_uX^2/2\) is strictly positive on its full root list.
- **Alice infinite:** use the exact homogeneous EN25 frames: finite unit \(v\) gives \(c_0D_vW^2+a_0D_vY^2+(d_0+1/2)(W+vX)^2+(vY-Z)^2/2\); infinitesimal \(v\) uses \((sX+Z,W,Y,vX)\), giving \(c_0W^2+a_0Y^2+(d_0+1/2)(W+vX)^2+(sX+Z)^2/2\); both infinite use \((W,X,Y,sZ,(1/v)Z)\), giving \(c_0W^2+a_0Y^2+d_0(X-sZ)^2+X^2/2+(Y-(1/v)Z)^2/2\). All these residue square forms have trivial common kernel. The exact original homogeneous polynomials define the literal endpoints; no zero reciprocal is inverted.

The positive-real caps are chosen *after* the individual real testing residues, justified by the real Schur penalty lemma; no cap must work uniformly across different independently chosen tests. Every retained excess is an original globally positive event packet, not a modified physical source. Signs, the zero state and complex states follow by exact congruences and real/imaginary splitting. This exhausts the independent projective testing charts and completes R3/R4.

## External dependencies and reproducibility

Used, but not re-proved in this R3/R4 sector: (i) the native source-entry theorem identifying the exact stationary physical C2 source and its positive gaps; (ii) the elementary bounded-positive-residue tensor-SOS lemma; (iii) the pure-testing signed-core Schur formula and exact sign/homogeneous congruences. **Not** used: any complete R1/R10 theorem, physical C1 admission, overall C2 theorem, Scheiderer conclusion, or fixed sixth-root candidate.

The companion `R34_FULL_PARENT_EXACT_CHECK.py` checks without sampling: RFN7 hard block; RFN21/22 complete derivative readback; all four original-source zero numerators \(S_N,S_0,S_1,S_2\); the *entire* RFN12 scalar against (9) after RSR6 substitution; RFN11, RFN18a/b/c, RFN19 with exact \(G_7,G_3\) and minors, RFN20; EN14; all six EN22 principal minors; and eight canceled source-coquotient limits (3). It does **not** claim to certify the order-based bounds by finite sampling: those are derived from the exact positive source polynomials and displayed ratios in this proof.

**Conclusion:** This closes the R3/R4 assigned local/all-tests leaf, conditional only on the three independently needed structural entry/tensor/core lemmas stated above. It supplies no conclusion about unreviewed C2 leaves or the global Q222 finite SDP lift.
