# Q222 native C2 R6/R7: whole-state bounded coefficient proof

**Scope (26 October 2026).** This is a coefficient-level proof for the R6/R7 portions of the *actual* native C2 contact family in `latex/certificates/c2_short.tex`, pinned at commit `171a7e94b8f52f2bb601562eda362badfb687a9a`. It uses the quarter-parent \(P=F/4\), keeps the source and testing real-closed fields independently labelled, and admits arbitrary relative infinitesimal rates. It does **not** give the global C1/C2 theorem, a qualitative full-Q222 existence theorem, or a preassigned finite affine lift. At the end, the exact imported, noncircular linear-algebra lemmas and remote identities are identified. The full scalar, both normal rows, and hard block below are identities, not evaluations at a contact.

## 1. Actual-source dictionary and uniform regularity

On weak-C2 R6/R7 put
\[
 s=1/a,\quad b=a^2B,\quad x=sb,\quad B=s^2b,\quad
 l=1-h,\quad R=1+hx,\quad w_s=ae=R-g,\quad
 e=sw_s,\quad r=1-e,
\]
\[
 K=1+x,\qquad T=1+xr,\qquad j_R=T-l/t,\qquad
 S_c=1+B(l+r/s),\qquad \Delta_*=s\Delta.
\]
Here the letter \(j_R\) is the original R10/R6/R7 gap, **not** the contact scale \(\jmath\) introduced below. The original BA1 denominator is
\[
\begin{split}
 \Delta={}&s^{-1}+K+KB+2xKr+s^{-3}B^2r^2
  +2s^{-2}Br+tB\{1+s^{-1}+x+s^{-2}Br(2-r)\}.
\end{split}
\]
Introduce the *polynomials with cancelled poles*
\[
\begin{split}
 M_*&=K(s+r)(h+t)-s w_s(1+tB)l,\\
 Q_*&=1+B(1+t+h/s+2r/s+t/s)\\
 &\quad+B^2\{(h+t)/s+[h(r+1)+t(1+r-r^2)]/s^2\},\\
 K_*&=(1+s)(h+t)K,\qquad \mathcal N=w_sQ_*-bK_*.
\end{split}\tag{1.1}
\]
These expressions are polynomials in \((s,b,h,t,g)\) after \(B=s^2b\) (there are no genuine negative powers). Substituting the **original native stationary C2 coefficients**, not independent formal weights, yields the following exact table:
\[
\boxed{\begin{array}{c|l}
 n/s^2&\mathcal N/(2r\Delta_*)\\
 mg&(1+tB)w_sT S_c/(2\Delta_*)\\
 C&T M_*/(2l\Delta_*)\\
 A_h:=Ah&tT(s+r)S_c/(2\Delta_*)\\
 \beta_D:=Dj_R/T&K(s+r)(1-txr)S_c/(2tT\Delta_*)\\
 \alpha/(s^3b^2)&(s+r)M_*/[2(1-hB)\Delta_*]\\
 n_1:=n\gamma&sK\mathcal N/(2R\Delta_*)\\
 n_2:=n\gamma^2&rK^2\mathcal N/(2R^2\Delta_*)
\end{array}}\tag{1.2}
\]
where \(\gamma=rK/(sR)\). The exact identity
\[
 2r\Delta_*\frac{n}{s^2}=w_sQ_*-bK_*\tag{1.3}
\]
follows already from the first row. Before bounding \(b\), all nonconstant terms in \(Q_*\) and \(\Delta_*\) contain \(s\) or \(x=sb\). As \(x,s\) are infinitesimal and \(w_s=R-g\) is bounded, \(Q_*\sim\Delta_*\sim1\), while \(K_*\sim h+t\) (R6) or \(\sim t\) (R7), a positive unit. The physical sign \(n\ge0\) therefore gives
\[
 0<b\le w_sQ_*/K_*,\qquad b\in O(1),\qquad n=O(s^2).\tag{1.4}
\]
Consequently \(s,b,x,B,g\) are bounded and \(R,K,T,S_c,\Delta_*,1-hB,t,1+s,r\) have positive unit denominators where used. In R6 additionally \(h,l,j_R\) are positive units; in R7 \(h,g\) are infinitesimal and \(l\sim1\), while **\(j_R\) is never inverted**. There is no comparison with a testing infinitesimal.

The joint residues of (1.2) are
\[
\begin{array}{c|ccccc}
 &mg&C&A_h&\beta_D&\alpha/(s^3b^2)\\\hline
 \mathrm{R6}&\tfrac12&(h+t)/(2l)&t/2&1/(2t)&(h+t)/2\\
 \mathrm{R7}&\tfrac12&t/2&t/2&1/(2t)&t/2
\end{array}.\tag{1.5}
\]
For R6 \(A=A_h/h\) and \(D=\beta_D/(j_R/T)\) are units with residues \(A_*=t/(2h)\) and \(D_*=1/[2(t-l)]\); \(C_*=(h+t)/(2l)\). In R7, \(A\) is unbounded, \(D\) is either a positive unit (\(t_0>1\)) or unbounded (\(t_0=1\)), and \(C_* =t/2\) remains a positive unit. The real caps in the remote joins are available because \(mg\sim1/2\) with infinitesimal \(g\), \(Ah\sim t/2\) with infinitesimal R7 \(h\), and \(\beta_D\sim1/(2t)\) with infinitesimal R7 \(j_R\). Weak \(n,\alpha\) are **always kept whole and positive** (they may literally vanish for \(n\)).

For \(b\) infinitesimal, take the two source-only scales
\[
 \varepsilon=K\sqrt{bt}/R,\quad
 \jmath=\sqrt{s/t}\,R/K,\quad k=\varepsilon\jmath=\sqrt x,\quad E=tK/R.\tag{1.6}
\]
Exactly \(b=\varepsilon^2R^2/(K^2t)\), \(s=t\jmath^2K^2/R^2\), \(\varepsilon,\jmath\in\operatorname{Inf}\), and the contact is \((u,\omega,v)=(k,\jmath,\varepsilon\jmath^2E)\). If \(b\) is a unit, \(\varepsilon\) is also a unit and \(\jmath\) is infinitesimal. The contact \(\omega=(W+uY)/D_u\) is a *testing coordinate*; it is not the source \(w_s\).

## 2. One exact whole-parent normal lemma for both R6 and R7

Let \(D_u=1+u^2,D_v=1+v^2\), and on the genuine state chart \(Z=1\) define
\[
 p=w_s/R,\quad c=l/(tT),\quad q=1-c=j_R/T,\quad
 A_0=v\omega,\quad V_0=1-puv\omega,
\]
\[
 N=X+u-pD_uA_0,\quad q_B=v(1-uX),\quad
 M_B=uW-Y-cq_B.\tag{2.1}
\]
These coordinates have the full inverse
\[
 X=N+pD_uA_0-u,\quad q_B=v(D_u-uN-puD_uA_0),\quad
 W=\omega+u(M_B+cq_B)/D_u,\quad Y=u\omega-(M_B+cq_B)/D_u.\tag{2.2}
\]
Set \(T_B=M_B+cq_B\), \(J=M_B-qq_B\), \(H=\alpha+1/2\), \(\mu_1=mg/R\), \(\mu_2=(mg)g/R^2\), \(n_1=n\gamma\), \(n_2=n\gamma^2\) and
\[
 F_0=\frac{p^2}{2}-n_2-2n_1p+\mu_2,\qquad
 A_N=p-2\mu_1-2n_1.\tag{2.3}
\]
The two source numerators \(F_0,A_N\) have **no poles at \(b=0\)**. More precisely, direct substitution of (1.2) gives the *polynomial numerators*
\[
\boxed{\begin{aligned}
 \Delta_*R A_N
  &=U_N:=w_s[\Delta_*-(1+tB)TS_c]-sK\mathcal N,\\
 2\Delta_*R^2F_0
  &=U_F:=w_s^2\Delta_*-(rK^2+2sKw_s)\mathcal N
                    +(1+tB)w_sTS_c\,g.
\end{aligned}}\tag{2.3a}
\]
They are polynomials in \(s,b,h,t,g\): the apparent \(s^{-1}\) in (1.1) has cancelled. The exact checker enumerates monomials and proves \(U_N\in(sb)\), \(U_F\in(b)\), with coefficient \([sb]U_N|_{g=0}=h+t\) and \([b]U_F|_{g=0}=h+t\). The divisible numerator is therefore established without taking a limit through a pole. To see the key cancellation by hand, at \(b=0\) one has \(K=R=T=S_c=1,\Delta_*=1+s,\mathcal N=w_s=1-g,r=1-sw_s\), so both numerators vanish identically. At \(s=0\), \(U_N=w_s(1-1)=0\) for arbitrary \(b,h,t,g\), proving the additional factor \(s\). This is a polynomial-factor argument, not division by either opening.

Equivalently, at \(b=0\), with \(w_s=1-g\),
\[
 p=w_s,\quad mg=\frac{w_s}{2(1+s)},\quad
 n_1=\frac{s w_s}{2(1+s)},\quad
 n_2=\frac{w_s(1-sw_s)}{2(1+s)}.\tag{2.4}
\]
Insertion gives \(F_0=A_N=0\) identically. All denominators are units there, so the *polynomial numerators* have a factor \(b\). Moreover the polynomial numerator of \(A_N\) also vanishes at \(s=0\) for arbitrary \(b,h,t,g\); hence \(A_N/(bs)\) has a unit denominator. The exact quotients satisfy
\[
 F_0/b\sim(h+t)/2,\qquad A_N/(bs)\sim h+t\quad(g=0),\tag{2.5}
\]
with R7 specialization \(t/2,t\) when \(h\) is infinitesimal. This is **not** an assumption of independent coefficient cancellation: (1.2)--(2.4) derive it from the native source.

Use the four separated raw monomials
\[
 \mathsf A=\varepsilon\omega v,\quad \mathsf B=v,\quad
 \mathsf C=Eu\omega,\quad \mathsf D=\varepsilon^2\jmath^3E,
 \qquad \mathsf D^2=s^3b^2t.\tag{2.6}
\]
Writing \(\mathcal G\) for the *whole* parent remainder after retaining the positive squares below, put
\[
 \mathcal G=G_0+N L_N+M_B L_M+
  H_{NN}N^2+2H_{NM}NM_B+H_{MM}M_B^2.\tag{2.7}
\]
The **same finite expansion** gives the exact scalar restricted support
\[
 G_0=f_A\mathsf A^2+f_{BB}\mathsf B^2+
 f_{BC}\mathsf B\mathsf C+f_{CC}\mathsf C^2+
 \varrho\mathsf D^2,\quad
 f_A=D_u^2D_vF_0/\varepsilon^2,\quad
 \varrho=\alpha D_u/\mathsf D^2.\tag{2.8}
\]

**R6:** retain \((m-K_m)N^2\) for a fixed positive real \(K_m<m\), together with the whole \(nD_uD_v(X+\gamma v\omega)^2\). The remaining \(\mathcal G\) is exactly RB11. Put \(L=AD_v+Hv^2\). Then
\[
\begin{aligned}
 f_{BB}&=[Cc^2D_u^2D_v+(D+1/2)q^2D_u^2+D_uLc^2]V_0^2
       +2HD_ucV_0-D_u/2-D_u^2\mu_2v^2\omega^2,\\
 f_{BC}&=[-2D_uLcV_0-2HD_u+2n_1D_uD_v]/E,\\
 f_{CC}&=[D_uL+D_uD_vn_2v^2]/E^2.
\end{aligned}\tag{2.9}
\]
Both **complete** normal rows are
\[
 L_N=R_A\mathsf A+R_B\mathsf B+R_C\mathsf C,\quad
 L_M=S_B\mathsf B+S_C\mathsf C,\tag{2.10}
\]
\[
\begin{aligned}
 R_A&=D_u A_N/\varepsilon,\qquad R_C=2cuvL/E,\qquad S_C=-2L/E,\\
 R_B&=D_u(p-2n_1)v^2\omega
  -2uv\,[Cc^2D_uD_v+Ac^2D_v+(D+1/2)q^2D_u+Hc^2v^2]V_0
  -2Hcuv,\\
 S_B&=[2CcD_uD_v-2(D+1/2)qD_u+2AcD_v+2Hcv^2]V_0+2H.
\end{aligned}\tag{2.11}
\]
The hard matrix entries are
\[
\begin{aligned}
 H_{NN}&=K_m+D_v/2+Cc^2D_vu^2v^2+(Ac^2D_v/D_u)u^2v^2
 +(D+1/2)q^2u^2v^2+(Hc^2/D_u)u^2v^4,\\
 H_{NM}&=-CcD_vuv-(AcD_v/D_u)uv+(D+1/2)quv-(Hc/D_u)uv^3,\\
 H_{MM}&=CD_v+AD_v/D_u+D+1/2+Hv^2/D_u.
\end{aligned}\tag{2.12}
\]
Here \(n_1,n_2,\alpha\) are present in the *full parent*, not dropped. They survive in scalar and first-normal entries after the selected positive \(n\) square is removed. At \(u=v=0\) the hard block is exactly \(\operatorname{diag}(K_m+1/2,A+C+D+1/2)\).

**R7:** retain the four whole squares
\[
 mN^2,\quad nD_uD_v(X+\gamma v\omega)^2,\quad
 AD_uD_v(Y-huW)^2,\quad D M_B^2.\tag{2.13}
\]
All vanish at the actual contact, including when \(n=0\). The exact subtraction (RU6)--(RU7) eliminates every unbounded \(m,A,D\) without dividing by the small gaps. Now
\[
\begin{aligned}
 f_{BB}&=[Cc^2D_u^2D_v+(\beta_Dq+q^2/2)D_u^2
 +H D_uv^2c^2-A_hD_uD_v(2u^2+hu^4)c^2]V_0^2\\
 &\quad+2HD_ucV_0-D_u/2-D_u^2\mu_2v^2\omega^2,\\
 f_{BC}&=[2A_hD_uD_v(lu^2-1)cV_0-2HD_uv^2cV_0
    -2HD_u+2n_1D_uD_v]/E,\\
 f_{CC}&=[A_h(2-h)D_uD_v+HD_uv^2+D_uD_vn_2v^2]/E^2.
\end{aligned}\tag{2.14}
\]
The same row form (2.10) holds, now with
\[
\begin{aligned}
 R_A&=D_uA_N/\varepsilon,\\
 R_C&=2cuv\{A_hD_v(1-lu^2)+Hv^2\}/E,\\
 R_B&=D_u(p-2n_1)v^2\omega
 -2uv\,[Cc^2D_uD_v+(\beta_Dq+q^2/2)D_u+Hc^2v^2
       -A_hc^2D_v(2u^2+hu^4)]V_0-2Hcuv,\\
 S_C&=2\{A_hD_v(lu^2-1)-Hv^2\}/E,\\
 S_B&=[2CcD_uD_v-(2\beta_D+q)D_u
    -2A_hcD_v(2u^2+hu^4)+2Hcv^2]V_0+2H.
\end{aligned}\tag{2.15}
\]
Finally
\[
\begin{aligned}
 H_{NN}&=D_v/2+Cc^2D_vu^2v^2
   -(A_hc^2D_v/D_u)(2u^2+hu^4)u^2v^2
   +(\beta_Dq+q^2/2)u^2v^2+(Hc^2/D_u)u^2v^4,\\
 H_{NM}&=-CcD_vuv+(A_hcD_v/D_u)(2u^2+hu^4)uv
  +(\beta_D+q/2)uv-(Hc/D_u)uv^3,\\
 H_{MM}&=CD_v-(A_hD_v/D_u)(2u^2+hu^4)+1/2+Hv^2/D_u.
\end{aligned}\tag{2.16}
\]
At \(u=v=0\) this is \(\operatorname{diag}(1/2,C+1/2)\).

**Exact verification and boundedness.** The checker `verify_arrays.py` derives (2.8)--(2.16) by coefficient extraction from the *whole parent* with free independent normals and establishes twelve zero polynomial differences, including all three hard entries of each family. This is stronger than checking stationary substitution. By (1.2), (2.4)--(2.5) the only potentially singular factors in (2.8)--(2.16) are
\[
 \frac{F_0}{\varepsilon^2}=\frac{R^2}{K^2t}\frac{F_0}{b},\quad
 \frac{A_N}{\varepsilon}=\frac{\varepsilon R^2s}{K^2t}\frac{A_N}{bs},\quad
 \frac{\alpha}{\mathsf D^2}=\frac1t\frac{\alpha}{s^3b^2}.\tag{2.17}
\]
All denominators on the right are units; the quotients have unit denominators, and there is no inversion of \(\varepsilon,\jmath,h,j_R,g,n\), nor of a mixed testing expression. The remaining coefficients depend only on the bounded collection
\(p,c,q,E,C,A\) (R6) or \(p,c,q,E,C,A_h,\beta_D\) (R7), \(\mu_1,\mu_2,n_1,n_2,H\) and \(D_u^{-1}\) (a pure-testing unit). Hence the **scalar, BOTH full normals and ALL hard entries** are bounded simultaneously, regardless of independent source/test rates.

## 3. Division-free conductor and actual coefficient filtration

Put \(\mathsf x=\mathsf B-\mathsf C\), \(\mathsf y=\mathsf A-k\mathsf C\),
\(\mathsf z=u\mathsf B-k\mathsf C\), \(\mathsf h=k\mathsf C-\mathsf D\). Then the following are literal identities:
\[
\begin{split}
 \mathsf C(u-k)&=\mathsf z-u\mathsf x,\\
 \varepsilon\mathsf C(\omega-\jmath)&=\mathsf y-\varepsilon\omega\mathsf x,\\
 \mathsf C(v-\varepsilon\jmath^2E)&=E(\omega\mathsf z+\jmath\mathsf y-k\omega\mathsf x),\\
 \mathsf C v(\omega-\jmath)&=E\omega(\omega-\jmath)\mathsf z+E\jmath\omega\mathsf y-Ek\omega^2\mathsf x.
\end{split}\tag{3.1}
\]
For any polynomial coefficient \(S\) with \(S_c=0\), do **monic division** in order \(u,\omega,v\):
\[
 S=(u-k)T_u+(\omega-\jmath)T_\omega+(v-\varepsilon\jmath^2E)T_v.
\]
If \(T_\omega=\varepsilon U+vV\), with polynomial \(U,V\), (3.1) implies exactly
\[
\begin{aligned}
\mathsf CS={}&\mathsf x[-uT_u-\varepsilon\omega U-Ek\omega^2V-Ek\omega T_v]\\
 &+\mathsf y[U+E\jmath\omega V+E\jmath T_v]\\
 &+\mathsf z[T_u+E\omega(\omega-\jmath)V+E\omega T_v].
\end{aligned}\tag{3.2}
\]
These transformations use only addition, multiplication and exact polynomial coefficient division; they never divide by a source opening or by the contact Jacobian.

The hypotheses for the **actual arrays** are now explicit. After clearing a fixed finite product of source unit denominators and powers of \(D_u\), every occurrence of \(\omega\) in the **five scalar coefficient functions** (2.8)--(2.9)/(2.14) is through \(uv\omega\) or \(v^2\omega^2\); the other angular monomials are even in \(u,v\). The *normal coefficients* do additionally contain \(v^2\omega\) in \(R_B\), and the monomial \(\mathsf A=\varepsilon\omega v\) itself, so the preceding sentence is **not** asserted for an unsorted full first-normal polynomial. At \(u=k=\varepsilon\jmath\), the first scalar \(\omega\)-divided difference is in \((\varepsilon v,v^2)\), so choose \(U,V\in(v)\). The next \(v\)-quotient belongs to \((\varepsilon,v)\), and the first \(u\)-quotient vanishes at \(\varepsilon=u=v=0\). This proves the first application of (3.2) to the coefficient of \(\mathsf C^2\) after the substitutions \(\mathsf A=\mathsf y+k\mathsf C\), \(\mathsf B=\mathsf x+\mathsf C\), \(\mathsf D=k\mathsf C-\mathsf h\). Its contact value vanishes by the *whole parent* and the retained contact-zero squares.

Its post-division \(\mathsf h\) coefficient is \(-2k\varrho\). The full first angular jet gives the **single gradient relation** \((-2k,1,1,-1)\) on \((\mathsf x,\mathsf y,\mathsf z,\mathsf h)\); the first three gradients have determinant \(-\varepsilon^3\jmath^4E^2\ne0\), but its inverse is **not used**. Subtract the contact multiple \(2k\varrho_c(-2k,1,1,-1)\) and use the exact syzygy
\[
 k\mathsf C(\mathsf y+\mathsf z-\mathsf h-2k\mathsf x)
       =Q_1^2-\mathsf y\mathsf z,\tag{3.3}
\]
with \(a_i^2=(\mathsf A,u\mathsf B,k\mathsf C,\mathsf D)\), \(Q_1=a_1a_2-a_3a_4\), \(Q_2=a_1a_3-a_2a_4\). The corrected second-stage coefficients have, at \(\varepsilon=u=v=0\), \(\omega\)-independent residues (the \(\mathsf x\) one restricts to \(2f_{BB}+f_{BC}\), the other three to zero); their next \(\omega\)-divided differences are in \((\varepsilon,v)\). Apply (3.2) again, choosing \(V=(T_\omega-T_\omega|_{v=0})/v\), \(U=T_\omega|_{v=0}/\varepsilon\), both **polynomial factor extractions**. The raw scalar support remains *only* \((\mathsf A,\mathsf B,\mathsf C)^2+(\mathsf D^2)\).

For the two full normal rows (2.10), substitution produces \(S_N=kR_A+R_B+R_C\) and \(S_M=S_B+S_C\) multiplying \(\mathsf C\); their values at contact vanish from the independent full \(N,M_B\) first jets. Their residues at \(\varepsilon=u=v=0\) are independent of \(\omega\) directly from (2.11)/(2.15), so their divided differences are in \((\varepsilon,v)\). One application of (3.2) lifts each normal with bounded coefficients. The hard entries already belong to the bounded coefficient ring by (2.17). This completes the **entire finite eight-root symmetric matrix**, without evaluating or discarding either normal.

For reference, the *strict circuit* identity is
\[
 T_4=\mathsf A^2+(u\mathsf B)^2+(k\mathsf C)^2+\mathsf D^2
       -4k^2\mathsf B\mathsf C
 =\frac{(\mathsf y-\mathsf z)^2+\mathsf h^2+\mathsf y^2+(\mathsf z+\mathsf h)^2}{2}+Q_1^2+Q_2^2.\tag{3.4}
\]
Its \((\mathsf y,\mathsf z,\mathsf h)\) contrast matrix has determinant \(1/2>0\). All \(a_i\) are products of a positive source root and a positive testing root on the physical positive chart. The bounded conductor (3.2) preserves the source ideal \((h,g)\) because all its coefficients and the contact substitutions have unit denominators and no inverses of \(h,g\).

### The finite source jets and coefficient filtrations

At \(g=0\), the **native** source gives the following exact Taylor coefficients (meaning polynomial divisibility after clearing unit denominators, not a numerical asymptotic substitution), with total degree at most two in \((s,b)\):
\[
\begin{aligned}
 C&=\frac{h+t}{2l}-\frac{t+1}{2l}s+\frac{t+1}{2l}s^2+O_{\deg(s,b)}(3),\\
 A_h&=\frac t2(1-s+s^2)+O_{\deg(s,b)}(3),\\
 \beta_D&=\frac{1-s+s^2-(t+1)sb}{2t}+O_{\deg(s,b)}(3),\\
 n_1&=\frac{s-s^2-(h+t)sb}{2}+O_{\deg(s,b)}(3).
\end{aligned}\tag{3.5}
\]
RI uses \(h,l,j_R\) units and \(A=A_h/h,D=\beta_D/q\); RU specializes \(h=0\) *only in this filtered leading model*, not in the actual positive source. These jets follow by expanding (1.2) and dividing by the unit \(\Delta_*\); `check_jets.py` checks each coefficient exactly.

For **R6, unit \(b\)** the expansion is in \(s\), with \(b\) a bounded source unit, and the actual first coefficients are
\[
\begin{aligned}
 mg&=1/2+(bh-1)s/2+O(s^2),\quad C=C_*-(t+1)s/(2l)+O(s^2),\\
 A&=A_*(1-s)+O(s^2),\\
 D&=D_*-[bt(h+t)+t-l]s/[2(t-l)^2]+O(s^2),\\
 \alpha&=s^3b^2(h+t)/2+O(s^4).
\end{aligned}\tag{3.6}
\]
Use \(A_0=\omega v,B_0=\delta v,C_0=\delta Eu\omega,D_0=\delta^2k^3E\) with \(\delta=\varepsilon^{-1}\) a source unit and \(\mathsf x_0=B_0-C_0\). Weight \((k,u,\omega,v)=(1,1,1,2)\). The complete leading collection of the scalar and both normals, modulo \(g\), is
\[
\begin{aligned}
 G_0&=\frac b{2h}(1+u^2)\mathsf x_0^2
      +\frac{b(h+t)}2T_4
      -\frac{k^2[1+bh(h+t)]}{2h}\mathsf x_0^2+O_{\rm wt}(8),\\
 L_N&=O_{\rm wt}(5),\quad
 L_M=\mathsf x_0/(h\delta)+O_{\rm wt}(4).
\end{aligned}\tag{3.7}
\]
The hard block has residue \(\operatorname{diag}(K_m+1/2,A_*+C_*+D_*+1/2)\). Restoring the contact circuit, the strict coupled block on \((\mathsf x_0,\delta M_B)\) is
\[
 \begin{pmatrix}b/(2h)&bt/(2h)\\bt/(2h)&bt(A_*+C_*+D_*+1/2)\end{pmatrix},\quad
 \det=\frac{b^2t}{2h}(C_*+D_*+1/2)>0.\tag{3.8}
\]
Here the \(N\) pivot is \(K_m+1/2>0\) and the five circuit roots have the positive coefficient \(b(h+t)/2\). The correction lifted by (3.2) has only strictly positive weights or a factor \(g\), hence all Gram corrections have infinitesimal coefficients.

For **R6, infinitesimal \(b\)**, weight \((\varepsilon,\jmath,u,\omega,v)=(1,1,2,1,3)\), \(\kappa=(h+t)/(2t)\),
\[
 \Lambda=\frac{1-s+s^2+u^2-[2l+h(h+t)]k^2}{2ht}.
\]
The *raw* coefficients from (2.9) obey the exact source-jet/weight statements
\[
 \begin{gathered}
 f_A-\kappa,\varrho-\kappa\text{ have weight }\ge2,\\
 f_{BB}-\Lambda-\kappa u^2,\quad f_{BC}+2\Lambda+4\kappa k^2,\quad
 f_{CC}-\Lambda-\kappa k^2\text{ have weight }\ge6,\\
 S_B-1/[h(1+s)],\quad S_C+1/[h(1+s)]\text{ have weight }\ge4,\\
 (R_A,R_B,R_C)\text{ have weights }\ge(3,5,5).
 \end{gathered}\tag{3.9}
\]
The exact scalar difference from \(\Lambda\mathsf x^2+\kappa T_4\) starts at weight 12, the \(L_M-\mathsf x/[h(1+s)]\) difference at 7, and \(L_N\) at 8. The source jets (3.5) yield every low-weight term; higher weights have no real residue and are lifted with bounded coefficients by (3.2), independently of source/test rates. The strict coupled block on \((\mathsf x,M_B)\) is
\[
 \begin{pmatrix}1/(2ht)&1/(2h)\\1/(2h)&A_*+C_*+D_*+1/2\end{pmatrix},\quad
 \det=\frac{C_*+D_*+1/2}{2ht}>0.\tag{3.10}
\]
The circuit coefficient \(\kappa>0\) and \(N\)-pivot \(K_m+1/2>0\) complete the strict bounded eight-root Gram. If \(b\) is infinitesimal, the source estimate \(n/s^2\sim1/2\) proves \(n>0\); \(n=0\) can only occur in the unit-\(b\) branch.

For **R7**, both unit and infinitesimal \(b\), the actual centered differences use (2.14)--(2.16). Let \(\Lambda=1/[t(1+s)]\) and
\[
 G_{\rm mod}=\Lambda\mathsf x^2+\frac12T_4-
       (1/t+1/2)k^2\mathsf x^2+(k^2-u^2)\mathsf C\mathsf x/t.\tag{3.11}
\]
The last term is *exactly* \((u+k)u\mathsf x^2/t-(u+k)\mathsf x\mathsf z/t\), since \(\mathsf C(u-k)=\mathsf z-u\mathsf x\). This is an infinitesimal Gram perturbation on the contact-root frame; it cannot be discarded before the full identity is assembled.

In the unit-\(b\) case give \((k,u,\omega,v)\) weights \((1,1,1,2)\); in the infinitesimal-\(b\) case use \((\varepsilon,\jmath,u,\omega,v)=(1,1,2,1,3)\). In both cases
\[
\begin{array}{c|cc}
 &b\text{ unit}&b\text{ infinitesimal}\\\hline
 f_A-1/2,\varrho-1/2&\ge2&\ge2\\
 (f_{BB}-F_{BB},f_{BC}-F_{BC},f_{CC}-F_{CC})&\ge4&\ge6\\
 S_B-1/(1+s),S_C+1/(1+s)&\ge2&\ge4\\
 (R_A,R_B,R_C)&\ge(2,3,3)&\ge(3,5,5)
\end{array},\tag{3.12}
\]
where
\[
 F_{BB}=\Lambda+u^2/2-(1/t+1/2)k^2,\quad
 F_{BC}=-2\Lambda-u^2/t+(3/t-1)k^2,\quad
 F_{CC}=\Lambda+u^2/t-2k^2/t.
\]
Thus the full raw scalar correction begins at weight \(8\) or \(12\), the \(M_B\)-normal correction at \(4\) or \(7\), and the \(N\)-normal correction at \(5\) or \(8\). Every \((h,g)\)-correction is lifted boundedly by (3.2) and is infinitesimal, with *no restriction on \(h/g\), \(j_R/h\), \(s/b\)*. The hard-block residue is \(\operatorname{diag}(1/2,(t+1)/2)\). The complete eight-root residue has the strict circuit \(T_4/2\), the \(N\)-pivot \(1/2\), and
\[
 \boxed{\begin{pmatrix}1/t&1/2\\1/2&(t+1)/2\end{pmatrix}\succ0,
 \quad\det=(t+2)/(4t)>0.}\tag{3.13}
\]
The exact jet checker `verify_source_jets_table.py` establishes all six R6/R7 infinitesimal-\(b\) scalar coefficient comparisons and the corresponding \(u^2\) coefficients at residue; `verify_unit_b_RU.py` checks all six R7 unit-\(b\) constant and first-\(s\) comparisons. The formulas for each *actual* complete normal row were checked as polynomial identities independently in `verify_arrays.py`. These are coefficientwise checks, not a scalar sample or a determinant-only argument.

## 4. Remaining state and independent projective testing charts

On the negative-core chart, the pure-testing sign/Schur identity SD5 permits \(u,v>0\). If \(Z\) is a unit, normalize \(Z=1\); \(\omega>0\) is bounded. Infinitesimal \(\omega\) is treated by §3. When \(\omega\) is a testing unit, retain the *original* positive \(nD_uD_vX^2\) instead of its contact-centered square, so the \(n_1,n_2\) negative lines of (2.9)/(2.14) disappear. On the exact frame \((A_o,\mathsf x,\mathsf D,N,M_B)\), \(A_o=\omega v\), use \(\mathsf A=\varepsilon A_o\), \(\mathsf B=A_o/\omega\), \(\mathsf C=A_o/\omega-\mathsf x\). All coefficient quotients cancel using only the **bounded pure-testing** inverse \(1/\omega\); the residual \(N\) pivot and \((\mathsf x,M_B)\) block are respectively those of (3.10) for R6-inf-\(b\), and (3.13) for R7. In R6 unit-\(b\) the \(\delta\)-scaled block (3.8) applies. These are entire-state identities with no frozen first normal. If the state instead has infinitesimal \(Z\) after normalization by a maximal state component, \(W\) is a unit on the negative-core branch; the algebraic BE4 selected frame \((X,Y,uW,vW)\) has real residue \(2d[X^2+(uW-Y)^2]+4m_c(vW-X)^2+4A_cY^2\succ0\). The coefficients \(m_c,A_c\) are **positive real caps** below unbounded actual \(m,A\); retain all other nonnegative event packets, including \(n,\alpha\). This exact imported selected-frame identity requires only those two caps, independent infinitesimal \(u,v\), and \(W\) a unit/\(Z\) infinitesimal. It imports no conclusion about R6/R7 or the full C1 theorem.

For R6 with \(u\) a testing unit and \(v\) infinitesimal, the entire remainder in \((N,R_t,Y,V_Z)\), with \(R_t=uW-Y,V_Z=vZ\), has residue RB22; its three-state determinant is
\[
 (1+u^2)^2A_*C_*(D_*+1/2)u^2>0.\tag{4.1}
\]
For finite remote \(v\) a testing unit, removing the unbounded centered \(m\)-row and restricting to its kernel yields the explicit three-state determinant
\[
 \frac{(1+U)^2(1+V)^2(h+t)}{8hl(t-l)}
  [V^2+t(h+t)U^2V+hUV^2]>0\quad(V>0).\tag{4.2}
\]
This was checked exactly from the original whole parent after eliminating \(M=0\) (`check_remote.py`). A positive real cap below \(m\) restores strictness off the kernel by a finite-dimensional Schur complement. The \(Z=0\) state restriction has positive \(A,C\) event squares and a nonnegative signed core. Literal projective endpoints follow from the homogeneous R6.5--R6.7 four-state matrices, specialized to \(n_0=\alpha_0=0\): their nontrivial corner block determinants are \(D_*/2\) and \(A_*(K_m+1/2)\), both positive. The \((\infty,0)\) exceptional homogeneous BE frame retains positive real caps of \(C,A\) and has strict roots \((p_uX+Z,W,Y,vX)\). These imported exact identities require precisely positive real lower bounds on \(C,A,D\) and an unbounded \(m\), supplied by (1.5); they are **not** conclusions imported from another R6 theorem.

For R7, cap only the unbounded \(m,A\), and, when \(t_0=1\), the unbounded \(D\), all by positive real constants. Retain original \(n,\alpha\). For finite remote \(V=v^2>0\), on the common \(Y=M=0\) kernel the remaining \((W,Z)\) form is
\[
 CD_vU W^2+(D_c+1/2)(uD_vW-vD_uZ)^2
       +VD_vW^2/2-D_uVZ^2/2.
\]
Its determinant is identically
\[
 VD_uD_v\left\{D_cV/2+U[CD_c-(D_c+1/2)/2]
                     +C(D_c+1/2)U^2\right\}.\tag{4.3}
\]
If \(t_0>1\), use the *actual* residues \(C=t/2,D_c=1/[2(t-1)]\) to obtain \(D_uD_v[V^2+t^2U^2V]/[4(t-1)]>0\). If \(t_0=1\), use \(C=1/2\) and a sufficiently large positive real \(D_c<D\): the braces are \(D_c(V+U^2)/2+(U^2-U)/4>0\) for every fixed real testing \(U\ge0,V>0\). A real cap lemma then gives positive real \(m_c,A_c\) to make the full state residue strict.

If \(u\) is a positive testing unit and \(v\) infinitesimal (including literally zero), the full four-root remainder in \((N,R_t,Y,V_Z)\) is RU27, with restriction \(Y=0\) of determinant \(t^2D_u(D_u-1)/[4(t-1)]>0\) when \(t_0>1\). At \(t_0=1\), the common kernel \(Y=0,R_t=D_uV_Z\) has strictly positive value \(D_uu^2V_Z^2/2\). Thus the same cap lemma yields strict entire-state PSD residues. On reciprocal charts the **full homogenized parent**, not a pointwise state restriction, has respectively the following strict square rows (RU28--RU30):
\[
\begin{array}{ll}
 (u=\infty,v\text{ unit}):& CD_vW^2+A_cD_vY^2+(D_c+1/2)(W+vX)^2+(m_c+1/2)(vY-Z)^2,\\
 (v=\infty,u\text{ bounded}):&m_c(W+uY)^2+C(uW-Y)^2+(A_c+1/2)D_uY^2+D_c(uX-Z)^2+D_uX^2/2,\\
 (u,v=\infty):& CW^2+A_cY^2+D_c(X-P_Z)^2+X^2/2+(m_c+1/2)(Y-Q_Z)^2,
\end{array}\tag{4.4}
\]
with the appropriate homogeneous state coordinates \(P_Z=u^{-1}Z,Q_Z=v^{-1}Z\). They have trivial common kernels, with last two nontrivial block determinants \(D_c/2\) and \(A_c(m_c+1/2)\). At \((u,v)=(\infty,0)\), use the same exact BE selected frame as in R6, with positive real caps \(C,A_c\); all other events are retained as nonnegative squares. Literal zero/infinite angles are polynomial evaluations on their respective homogeneous charts; no inverse is taken at zero. Nonnegative-core tests are direct original event/pure-testing-core squares. Both angle signs follow from exact sign congruences, and real/imaginary state components add.

**Strict bounded-Gram step.** In every nontrivial local or remote chart the complete remaining source/testing coefficient matrix is bounded and has a positive-definite *real* residue, by (3.8), (3.10), (3.13), or (4.1)--(4.4). Apply the established strict bounded-Gram tensor lemma (WA10): any finite bounded symmetric matrix over the separated bounded coefficient ring with positive-definite real residue admits a finite tensor-square decomposition. Restore each retained nonnegative event/centered row and each positive pure-source or pure-testing square factor. This uses no fixed uniform degree across source strata and no inversion of a mixed-source/testing multiplier.

## 5. Exact dependencies, proved checks and status

**External inputs, with precise noncircular hypotheses:** (i) SG2 native quarter-parent, native stationary source and physical gaps, plus exact sign congruences; (ii) the previously checked CA/SE monic contact-conductor lemma CC1--CC6 / SE5--SE7, applicable exactly to a contact-zero scalar with raw support \((\mathsf A,\mathsf B,\mathsf C)^2+(\mathsf D^2)\) and both contact-zero normal rows, whose source unit denominators are bounded and whose two stages satisfy \(T_\omega\in(\varepsilon v,v^2)\) and \(T_\omega\in(\varepsilon,v)\); (iii) WA10 strict bounded-Gram lemma; (iv) SD5 sign/core reduction; (v) BE4 and the homogeneous BE/R6.5--R6.7 selected-frame *polynomial identities* under the positive-real-cap hypotheses supplied here, not any C1/C2 global admission theorem. The coefficient verification in §§1--3 proves every new premise needed by (ii), including both normals. The remote calculations in §4 establish the explicit determinants and state-kernel tests under (iv)--(v).

**Reproduction.** `verify_arrays.py` checks the exact full scalar, both normal rows and the hard block (12 symbolic identities). `verify_fast_divisibility.py` checks the complete polynomial numerators \(U_F\in(b)\), \(U_N\in(bs)\) and both quotient residues without large expression expansion. `check_jets.py`, `check_b_unit_jets.py`, `verify_source_jets_table.py`, and `verify_unit_b_RU.py` check actual-source residues and all low-weight raw coefficient comparisons. `check_conductor.py` checks every algebraic identity (3.1)--(3.2) and the full contact Jacobian. `check_RB17_coeffs.py` verifies all nine RB17 low-weight scalar coefficient equations for arbitrary \(b,h,t\). `check_remote.py` verifies the R6 kernel determinant (4.2). The scripts perform exact rational/polynomial algebra, not numerical SDP optimization. Their available outputs are logged in `verification.log`.

**Result:** Under the listed, independently established/noncircular algebraic lemmas, the actual complete C2 parent is admitted on R6 (\(b\) unit and infinitesimal, including permitted \(n=0\)) and R7 (both \(b\) alternatives and the physical \(n=0\) boundary), for all independent tests and state normals. This **does not close the other C2 strata or the global Q222 theorem**. It also does not prove an explicit fixed SDP lift, because the local contact frames and degrees are not assembled into a single fixed affine pencil.
