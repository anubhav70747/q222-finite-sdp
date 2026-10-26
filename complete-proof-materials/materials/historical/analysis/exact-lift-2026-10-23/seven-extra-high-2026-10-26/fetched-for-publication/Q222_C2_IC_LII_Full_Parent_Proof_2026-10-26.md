# Q222: weak-IH IC and strong-II LII, full-parent tensor-SOS admission

**Scope and verdict (26 October 2026).** This is an analytical certificate for two **separate native physical C2 source cells** of the normalized full parent. Weak-IH means actual contact squared angles `z` infinitesimal and `y` unbounded, `C` weak, `B=H+L` unbounded. Strong-II means `z,y` infinitesimal and `B=H+L` a positive bounded **unit**, with strong complementary events. Both impose the original native gaps and `n>=0`; `n=0` is literal and allowed. Source variables are in a real closed `K_c`, testing variables in an independent real closed `K_t`, and all identities are in `K_c tensor_R K_t`. No comparison of source and testing infinitesimal rates is used. The parent is certified for **all four homogeneous real state amplitudes**, both unrestricted normals, all signed projective measurement tests, and permitted event zeros. No mixed inverse or mixed square root is used.

The result is **conditional only on these external facts with the stated scope**: (i) the already checked canonical physical native C2 entry and its exact full state-kernel/angular stationarity identities for the original weights; (ii) the algebraic *strict bounded-residue Gram lemma* for separated fields; and (iii) the separate arbitrary-source, noncommuting-residue theorem only where explicitly mentioned (it is not invoked in the two delicate local charts). In fact the complementary charts below are direct parent calculations, so they do not invoke the whole C2 theorem. This document is not a complete proof of a Q222 SDP lift.

## I. Source formulas and the shared finite contact principle

Set `D_u=1+u^2`, `D_v=1+v^2`, and use the complete original C2 parent, divided by four,

\[
\begin{split}
P={}&nD_uD_vX^2+m(vW-X+uvY-uZ)^2+CD_v(uW-Y)^2+AD_uD_vY^2\\
&+(D+\tfrac12)(uW+uvX-Y-vZ)^2
 +(\alpha+\tfrac12)D_u(vY-Z)^2\\
&+\tfrac12D_v(X+uZ)^2-\tfrac12D_uD_vZ^2.
\end{split}\tag{1}
\]

For a reciprocal Bob chart `w=1/v`, the unnormalized homogeneous parent `\widehat F=4w^2P(u,1/w)` is instead

\[
\begin{split}
\widehat F={}&4nD_uD_wX^2+4m(W-wX+uY-uwZ)^2
+4CD_w(uW-Y)^2+4AD_uD_wY^2\\
&+4D(uwW+uX-wY-Z)^2+4\alpha D_u(Y-wZ)^2\\
&+2\{(uwW+uX-wY-Z)^2+D_u(Y-wZ)^2
  +D_w(X+uZ)^2-D_uD_wZ^2\}. 
\end{split}\tag{2}
\]

Equation (2) defines `w=0` literally and is an ordinary polynomial in testing variables. The other projective endpoints use analogous original homogeneous forms, **never** a raised-degree contact graph evaluated at infinity.

The two source cells use the same exact rational contact dictionary (`WI1` in the cited source). Define

\[
R=zr,\quad L=B-H,\quad \mathcal S=1+R+L,\quad
G=B+R-z(1-H),
\]
\[
\Delta=(1+z)B^2+2(1+z)BR+(1+z)(1+y)B
 +(z-y)R^2+2z(1+y)R+z(1+y),
\]
\[
\mathcal N=(1+z)(B+R)(H+y)-(1+y)(z-R)L.
\]

The original six coefficients, with no asymptotic replacement, are

\[
\begin{array}{ll}
C=(1+R)\mathcal N/(2L\Delta), &
\alpha=(B+R)\mathcal N/[2(1-H)\Delta],\\
A=y(1+R)(B+R)\mathcal S/(2H\Delta),&
D=(1+z)(B+R)(B-yR)\mathcal S/[2\{y(1+R)-L\}\Delta],\\
m=z(1+y)(z-R)(1+R)\mathcal S/(2G\Delta),&
n+\tfrac12=z(B+R)(1+R)\mathcal S/(2R\Delta).
\end{array}\tag{3}
\]

All denominators here are positive on the **actual native source domain**. Importantly, coefficients which are infinitesimal need not be set to zero. The contact in the reciprocal graph has `u_c^2=z, w_c^2=1/y, W_c/Z_c=B/(sqrt(zy))` and `X_c=-r u_c Z_c`, `Y_c=H w_c Z_c`.

**Contact-division lemma.** Suppose a quadratic parent and a proposed sum of contact-centered root products agree modulo
\[
\mathcal E=A\,uwWZ+(B_1u^2+B_2w^2+B_3)Z^2
\tag{4}
\]
with *source-only* `A,B_i`. Suppose the parent has zero value and full first derivatives at the above finite contact, and all proposed root products have zero value and first derivatives there. Because `u_c,w_c,Z_c` are nonzero, `\partial_W\mathcal E=0` implies `A=0`; then `\partial_u\mathcal E=0` and `\partial_w\mathcal E=0` give `B_1=B_2=0`, and the value gives `B_3=0`. This is a **literal polynomial identity**, requiring no genericity, density, or evaluation of source/testing variables in the same field. It also works over real closed coefficient fields. This is how IC7 is proved below without a large source substitution. For LII, the analogous assertion is monic division at a double contact, proved explicitly in III.

## II. Weak-IH: interior source shape, both normals and endpoint shapes

Assume `z` positive infinitesimal, `y` positive unbounded, `B` unbounded, and put `rho=B/y`; all native inequalities, including `0<H<1` and `n>=0`, remain in force. From (3), after **positive-unit rational cancellation**,
\[
BC,\quad Bm/z^2,\quad AH/B,\quad \alpha(1-H)/B
\tag{5}
\]
are bounded positive units, `m,C` are infinitesimal, and `nr` is bounded. In the interior `0<rho_0<1`,
\[
\begin{gathered}
(BC)_0=(AH/B)_0=(\alpha(1-H)/B)_0={1\over2(1+\rho_0)},\\
(Bm/z^2)_0={1-r_0\over2(1+\rho_0)},\quad
\nu_0:=[(4n+2)r]_0={2\rho_0\over1+\rho_0},\quad
D_0={\rho_0^2\over2(1-\rho_0^2)},\quad
0\le r_0\le {\rho_0\over1+\rho_0}.
\end{gathered}\tag{6}
\]
No inverse of `r`, `H`, or `1-H` is introduced in a bounded Gram coefficient.

At independent `u,w` infinitesimal or zero, put **without restricting the state**
\[
X=-ruZ+X',\quad Y=HwZ+Y',\quad
E_X=\sqrt{4n+2}\,X',\quad
E_Y=\sqrt{B/[H(1-H)]}\,Y'.
\tag{7}
\]

On `X'=Y'=0` use the five contact roots
\[
\begin{gathered}
R_0=(uW-BwZ)/\sqrt B,\quad R_1=(zW-BuwZ)/\sqrt B,\\
R_2=uwW-\rho Z,\quad T=(u^2-z)Z,\quad
E=\sqrt{BH}(w^2-1/y)Z.
\end{gathered}\tag{8}
\]

Set
\[
K=4AH^2+8DHr+4\alpha(1-H)^2+4m(H+r-1)^2
 +2(H+r-1)^2+4nr^2,
\]
\[
k=K-4B^2m/z^2+8znr^2+8AH^2/y.\tag{9}
\]

The exact graph certificate is **IC7**: its diagonal root coefficients are
\[
4BC-kz/B,\quad4Bm/z^2,\quad2(2C+2D+1),\quad
r^2\{4nw^2+2(2D+2n+1)\},\quad
(H/B)\{4Au^2+2(2A+2C+2D+1)\},
\]
with cross products
\[
(k u/B)R_0R_1+[-4r(2D+1)+k/B]R_2T
  -4\sqrt{H/B}(2C+2D+1)R_2E.
\tag{10}
\]
The checker independently expands (2) on the graph and proves that its difference from (10) has **only** the four monomials (4). Each root in (8) vanishes at the original contact. The contact-division lemma therefore proves the full identity using the *original* contact value and angle/state jets. This is not determinant-only reasoning. All its coefficient combinations are bounded by (5), `nr`, and `D` in the interior.

The exact zero identities
\[
B(u^2-z)(w^2-1/y)Z^2
=(u^2-z)ZR_2-zR_0^2+uR_0R_1,\tag{11}
\]
\[
\rho T^2+(u^2-z)TR_2-u(zR_0-uR_1)(uR_0-R_1)=0\tag{12}
\]
are direct polynomial products. Add (12) to (10). All newly introduced coefficients except the `+rho T^2` diagonal have zero residue, independent of source/testing rates. The scalar residue has three strictly positive diagonal pivots on `(R0,R1,E)`:
\[
{2\over1+\rho_0},\quad{2(1-r_0)\over1+\rho_0},\quad{2\over1+\rho_0},
\]
plus the `(R2,T)` block
\[
\begin{pmatrix}
 2/(1-\rho_0^2)&-r_0/(1-\rho_0)\\
 -r_0/(1-\rho_0)&\rho_0+r_0\nu_0+4D_0r_0^2
\end{pmatrix}.\tag{13}
\]
Its determinant is exactly
\[
{2\rho_0\over1-\rho_0^2}
+{r_0\{4\rho_0-r_0(3\rho_0+1)\}\over
 (1-\rho_0)(1+\rho_0)^2}>0
\tag{14}
\]
since `r0<=rho0/(1+rho0)`. These facts also cover literal `n=0`.

**Both normal readbacks.** Write `U=u^2,V=w^2`, and let `L_X=(1/4)\partial_X\widehat F|graph`, `L_Y=(1/4)\partial_Y\widehat F|graph`. Direct expansion, retaining every state term, gives
\[
\begin{aligned}
L_X={}&(2D+1)uR_2-(2m\sqrt B/z)wR_1+uZP_X(U,V),\\
L_Y={}&\sqrt B\,[2m-2C(1+V)-(2D+1)V]R_0+wZP_Y(U,V).
\end{aligned}\tag{15}
\]
The explicit finite polynomials are
\[
\begin{split}
P_X={}&(2D+1)\rho-2mBV/z-2nr(1+U)(1+V)
 -r(2D+1)U-2D-r\\
&+V\{1-H-r-2DH-2m(H+r-1)\},\\
P_Y={}&B[2m-2C(1+V)-(2D+1)V]
 +2AH(1+U)(1+V)+2CH(1+V)\\
&+2DHV+2DrU+2D-2\alpha(1-H)(1+U)
 +(2m+1)(H+r-1)U+H(1+V).
\end{split}\tag{16}
\]
Contact stationarity implies `P_X(z,1/y)=P_Y(z,1/y)=0`. Define the monic quotients `P_X=(U-z)f_X+(V-1/y)g_X`, `P_Y=(U-z)f_Y+(V-1/y)g_Y`. **Without dividing by a test difference**, they are
\[
\begin{array}{ll}
f_X=-2nr(1+V)-r(2D+1),&
 g_X=-2mB/z-2nr(1+z)+1-H-r-2DH-2m(H+r-1),\\
f_Y=2AH(1+V)+2Dr-2\alpha(1-H)+(2m+1)(H+r-1),&
 g_Y=-B(2C+2D+1)+2AH(1+z)+2CH+2DH+H.
\end{array}\tag{17}
\]
Inserting (11) and `B(V-1/y)Z=R2-sqrt(B) wR0` gives exactly
\[
\begin{split}
L_X={}&(2D+1+g_X/B)uR_2-(2m\sqrt B/z)wR_1
       +uf_XT-(uwg_X/\sqrt B)R_0,\\
L_Y/\sqrt{B/[H(1-H)]}={}&\sqrt{H(1-H)}(f+zf_Y/B)R_0\\
&-\sqrt{H(1-H)}(uf_Y/B)R_1
 +(\sqrt{1-H}wg_Y/B)E,
\end{split}\tag{18}
\]
where `f=2m-2C(1+V)-(2D+1)V`. Every coefficient in (18) is bounded **with zero residue**: `m sqrt(B)/z = (Bm/z^2) z/sqrt(B)`; `f_Y/B,g_Y/B` are bounded by the normalized source packets (5); and the other terms carry `u,w`, or `1/sqrt(B)`. The normal quadratic Gram has residue `diag(1,2/(1+rho0))` on `(E_X,E_Y)`; its cross-normal entry is infinitesimal. Therefore (13) plus these two normal pivots is one **strict bounded seven-row coefficient Gram** for the *entire* original parent. The positive-residue lemma yields separated tensor SOS.

**Shape endpoint `rho0=0`.** Put `s=r/rho` bounded, `tau=1/B` infinitesimal. Substitute `B=1/tau`, `y=1/(tau rho)`, `r=s rho` in (3). The normalized rational numerators `Delta/[B(B+y)]`, `N/[yB]`, `S/B`, `G/B` have denominators which remain positive units; the source residues give
\[
\nu/\rho\to2,\quad nr/\rho\to(1-s_0)/2,\quad
D/\rho^2\to1/2,\quad
BC,Bm/z^2,AH/B,\alpha(1-H)/B\to1/2,
\]
\[
(k/B\rho)\to2s_0,\qquad 0\le s_0\le1.\tag{19}
\]
For the last and most delicate assertion, insert (3) in (9) *before estimating terms*: the numerator of `k/B`, after clearing only unit denominators, vanishes **identically at `rho=0`** for arbitrary `z,tau,H,s`; its derivative in `rho` at `rho=z=tau=0` is `2s`. The checker performs these two exact polynomial tests. Therefore no `z/rho` or `tau/rho` comparison is hidden. Replace `T` by `T_rho=rho T`. After multiplying (12) by `rho`, the `(R2,T_rho)` residue is
\[
\begin{pmatrix}2&-s_0\\-s_0&1+2s_0\end{pmatrix},
\quad\det=2+4s_0-s_0^2\ge2.
\tag{20}
\]
The other scalar pivots are `2,2,2`. In (18) the sole affected coefficient is `f_X/rho=-2(nr/rho)(1+V)-s(2D+1)`, bounded, multiplied by testing `u`. The two normal pivots are `1,2`. This gives a strict full-parent Gram even if `r` is infinitesimal and `n` unbounded.

**Shape endpoint `rho0=1`.** Let `tau=1/B` and `e=1+zr-rho+rho H tau>0` infinitesimal. The same *simultaneously normalized* exact rational source formulas give
\[
De\to1/4,\quad\nu\to1,\quad nr\to(1-2r_0)/4,
\quad BC,AH/B,\alpha(1-H)/B\to1/4,
\quad Bm/z^2\to(1-r_0)/4,
\tag{21}
\]
`0<=r0<=1/2`, and `(k-8DHr)/B ->r0`. Thus one cannot put the unbounded `D` into the bounded Gram. Retain the **full-state** positive row
\[
F_D=\sqrt D\,(J+eZ),\qquad J=uwW+uX-wY-Z.
\]
On the graph `Q=R2-rT-sqrt(H/B)E=J+eZ`. The identity valid for every state is
\[
R_2={F_D\over\sqrt D}+rT+\sqrt{H/B}E
  -u{E_X\over\sqrt{4n+2}}
  +w{E_Y\over\sqrt{B/[H(1-H)]}}.
\tag{22}
\]
The graph subtraction has the following exact readback, with `a=sqrt(H/B)` and `kap=k-8DHr`:
\[
\begin{split}
G-4DQ^2={}&(4BC-\mathrm{kap}\,z/B)R_0^2+4BmR_1^2/z^2
 +2(2C+1)R_2^2\\
&+r^2(4nw^2+4n+2)T^2+{H\over B}\{4A(1+u^2)+4C+2\}E^2\\
&+(\mathrm{kap}\,u/B)R_0R_1+(-4r+\mathrm{kap}/B)R_2T
 -4\sqrt{H/B}(2C+1)R_2E.
\end{split}\tag{22a}
\]
It follows by expanding `Q=R2-rT-aE` and applying the *already exact* toric zero (11); the checker verifies it with independent formal weights. This subtraction is essential when `D` is unbounded. 

Subtracting `4F_D^2`, using (22) in the remainder, and adding (12) leaves a *bounded* full-state Gram on `(R0,R1,T,E,F_D,E_X,E_Y)` with residue
\[
\operatorname{diag}(1,1-r_0,1+r_0-r_0^2,1,4,1,1)\succ0.
\tag{23}
\]
The normal readback is exactly (15)-(18) with the `D` event removed and `De` retained. Its only new constant coefficient is `De`, bounded by (21); `1/sqrt(D)` in (22) is an infinitesimal pure-source factor. No bare `D` enters the selected Gram. This covers `H` and `1-H` infinitesimal openings and literal `n=0`.

## III. Strong-II: exact four-packet subtraction and bounded double-contact division

Now `z,y` are positive infinitesimals, `B=H+L` is a positive bounded unit, `L=yt>0`, `R=zr`, `j=1+zr-t>0`; `t+j=1+zr` is a unit. The normalized exact weights (3), with no order restriction between `z` and `y`, yield
\[
\begin{gathered}
(CL)_0=(Dyj)_0={h\over2(1+h)},\quad
(m/z^2)_0={1-r_0\over2h^2(1+h)},\quad
(A/y)_0={1\over2h(1+h)},\\
(nr)_0={1\over2(1+h)}-{r_0\over2},\quad
\{\alpha(1-H)\}_0={h^2\over2(1+h)},
\quad 0<h=B_0\le1.
\end{gathered}\tag{24}
\]
Here `n>=0` gives `r0<=1/(1+h)`, so `1-r` is a unit. Set
\[
c=CL,\quad d=Dyj,\quad \ell=y(C+D)=c/t+d/j,\quad
q=D/(C+D),\quad g=yCD/(C+D)=cd/(cj+dt),\quad s=(1+zr)/B.
\tag{25}
\]
Then `c,d,g,s` are bounded positive units, `0<q<1` bounded, and `ell` has a positive real lower bound but may be unbounded. Choose **positive real caps** `0<ell_*<=ell`, `0<alpha_*<=alpha` after fixing this source-field point. They need not be uniform over different sources. Their excesses remain in whole positive packets. These caps are not required to work uniformly across sources.

Use the full invertible state graph
\[
X=-zrvW/B+\xi,\quad Y=HuW/B+\eta,\quad Z=uvW/B+\zeta.
\]
Define `k=sqrt(y)`, `Q_n=u xi-zeta`, `K=ukW`,
\[
R_h={vQ_n-su(V-y)W\over k},\qquad
N=\eta/k-qR_h,\qquad U=u^2,\ V=v^2.\tag{26}
\]
The **exact** two-event identity, before any asymptotic estimate, is
\[
\begin{split}
4C(uW-Y)^2+4DJ^2
={}&4\ell N^2+4gR_h^2
+8(d-c)KN/B-8gsKR_h+4(ct+dj)K^2/B^2,\\
4CV(uW-Y)^2={}&4CV\eta^2+4ctVK^2/B^2
      -8cVK(N+qR_h)/B.
\end{split}\tag{27}
\]
It is verified by expansion using `J/k=(1-q)R_h-N-(j/B)K`, `\eta=k(N+qR_h)`, and `uW-Y=k(tK/B-N-qR_h)`.

Retain **exactly**
\[
\Pi=4(\ell-\ell_*)N^2+4CV\eta^2
 +4n(1+U)(1+V)\xi^2
 +4(\alpha-\alpha_*)(1+U)(v\eta-\zeta)^2.\tag{28}
\]
These packets are nonnegative for every independent test, vanish to first order at the original contact, and include literal `n=0`. Define `G=F-\Pi`, with `F=4P` after the graph. Equations (1), (26)-(28) are a **complete finite coefficient-extraction algorithm** for every scalar, both full first normals, and the hard block. The independent checker expands them identically to the printed LII6, for *arbitrary formal source weights*, before imposing (3). In particular no `C` or `D` event-normal term is silently dropped.

**One dependent-root gauge and its bounded coefficient.** Treat `N,xi,zeta,R_h` as temporary formal roots. The coefficient of `R_h W` in `G` has the exact form `uk\,\mathcal R(U,V)`, where
\[
\begin{split}
\mathcal R={}&-8gs-{8cq\over B}V
-{8\alpha(1-H)q\over B}(1+U)V
+{8mq\over B}V\{B+zr+U(H-1)\}\\
&+{8AHq\over B}(1+U)(1+V)
-{4yj(1-q)\over B}+{4(H-1)q\over B}(1+U)V.
\end{split}\tag{29}
\]
All coefficients in (29) are bounded by (24)-(25); especially, neither `C`, `D`, `alpha`, `n` nor `ell` appears naked. This identity is checked exactly from (27) and the complete signed core. Replace `K R_h` using
\[
K R_h=uvWQ_n-sU(V-y)W^2.\tag{30}
\]
This is a polynomial zero identity on the actual dependent roots. The resulting **formal coefficient polynomial** has exactly
\[
\widetilde G=Q_0(N,\xi,\zeta,R_h)
 +W\{vP_\xi(U,V)\xi+uvP_\zeta(U,V)\zeta+ukP_N(U,V)N\}
 +\phi(U,V)W^2.\tag{31}
\]
The coefficient-extraction rules (27)-(31), rather than any enormous expansion, specify each entry *uniquely*.

Every coefficient of `Q0,Pxi,Pzeta,PN,phi` is bounded. This follows term by term from (27)-(29) after noting `m=z^2\widehat m` with `\widehat m` a unit, `nr` bounded, `A=y\widehat A` with `\widehat A` a unit, `c,d,g,q,s` bounded, and `alpha` occurring only through `alpha(1-H)` or the capped `alpha_*`. Specifically, the `n` remainder after retaining `4nDuDv xi^2` is `4nDuDv[z^2r^2VW^2/B^2-2zrvWxi/B]`, and the `alpha` remainder is `4alpha Du[(1-H)^2UVW^2/B^2-2(1-H)uvW(v eta-zeta)/B]`: both are bounded. The initial `ell` and bare `C,D` terms have been removed by (27)-(28). The checker finds `deg_U(phi)<=2`, `phi(0,V)=V(d_s+e_sV)`, and `Pxi` in the bounded-coefficient ideal `(U,z)` (at `U=z=0` its only otherwise surviving formal term is `-8m=-8z^2\widehat m`). Here
\[
d_s={4\over B^2}\{m(B+zr)^2+(n+\tfrac12)z^2r^2\},\qquad
 e_s={4\over B^2}(n+\tfrac12)z^2r^2,\tag{32}
\]
and both `d_s/z^2` and `e_s/z^2` are bounded.

**Complete contact-division lemma.** At the physical contact `(U,V)=(z,y)`, every retained packet (28) is contact-centered, and the original full kernel and angle jets vanish. The gauge (30) removes the `R_h W` ambiguity. Hence
\[
P_\xi(z,y)=P_\zeta(z,y)=P_N(z,y)=0,\quad
\phi(z,y)=\phi_U(z,y)=\phi_V(z,y)=0.
\tag{33}
\]
With no inverse of a test difference, finite *monic polynomial divisions* yield
\[
\begin{aligned}
P_\xi&=(U-z)A_\xi+z(V-y)B_\xi,\\
P_\zeta&=(U-z)A_\zeta+(V-y)B_\zeta,\\
P_N&=(U-z)A_N+(V-y)B_N,\\
\phi&={d_s+e_sV\over z^2}\,V(U-z)^2
  +U(U-z)(V-y)A_s(V)+U(V-y)^2C_s(V).
\end{aligned}\tag{34}
\]
**Why the bounded division is exact.** Subtract the first term in `phi`. The remainder is `U[a_1(V)U+b_1(V)]`, because the constant coefficient `phi(0,V)` agrees by (32). The double zero at `U=z,V=y` makes both `a_1(y),b_1(y)` zero. Divide by the monic `V-y` to obtain `a_2,b_2`; the `V`-jet implies `b_2(y)+za_2(y)=0`. Set `A_s=a_2` and `C_s=(b_2+za_2)/(V-y)`. All coefficients are bounded, since `d_s/z^2,e_s/z^2` are bounded. For `Pxi`, first divide by `U-z`; the remainder `Pxi(z,V)` has a **factor z with bounded quotient**, by the `(U,z)` support proved above, then divide that quotient by `V-y`. The other two first normals require only ordinary monic divisions. The differences `(U-z)` and `(V-y)` are *never inverted inside the tensor ring*. This proves (34) for every independent rate and opening.

Use the two contact roots `R=v(U-z)W`, `T=uk(U-z)W`. The five **universal allocation identities** are
\[
\begin{aligned}
zv(V-y)W&=[uVQ_n-uvkR_h]/s-(V-y)R,\\
uv(V-y)W&=[VQ_n-vkR_h]/s,\\
uk(V-y)W&=[vkQ_n-yR_h]/s,\\
U(U-z)(V-y)W^2&=[-R_hT+uRQ_n]/s,\\
U(V-y)^2W^2&=[VQ_n^2-2vkR_hQ_n+yR_h^2]/s^2.
\end{aligned}\tag{35}
\]
All five are exact polynomial identities (checked). Since `s` is a positive bounded unit and `u,v,k,z,y` are infinitesimal, (34)-(35) give one bounded six-row coefficient Gram on `(N,xi,zeta,R_h,R,T)`, with **every state-normal term** included. The five-by-five residue excluding `T` has isolated positive pivots
\[
4\ell_*,\quad4\alpha_*,\quad4g_0
\]
and the sole nontrivial `(xi,R)` block
\[
\begin{pmatrix}
2&2/[h(1+h)]\\2/[h(1+h)]&2/[h^2(1+h)]
\end{pmatrix},\quad
\text{Schur pivot}={2\over h(1+h)^2}>0.\tag{36}
\]
The two exact source entries are `st(Axi)=4/[h(1+h)]`, `st(d_s/z^2)=2/[h^2(1+h)]`. The first follows from the unscaled `xi-W` coefficient `-4nzr(1+y)+4m(1-H)-(4D+2)yj+2(1+y) ->2/(1+h)` and its multiplier `2/B`; the second from `(m/z^2)_0` and `((n+1/2)r^2)_0=r_0/[2(1+h)]` in (24). Neither calculation divides by `n,t,j,1-H`.

Finally add a positive real `Lambda` times the exact **zero**
\[
T^2-UR^2-{U-z\over s}R_hT+{u(U-z)\over s}RQ_n=0.
\tag{37}
\]
Its only nonzero residue coefficient is `+Lambda T^2`. The five-by-five matrix (36) is strictly positive, and every `T` cross entry is bounded, so a sufficiently large *real* `Lambda` makes the entire six-by-six residue positive by the ordinary Schur complement. Its entries are in the **bounded separated tensor ring** and its roots need not be bounded. The strict bounded-Gram lemma therefore gives an exact finite separated SOS for `G`. Adding (28) gives the *entire original parent*, including `n=0`, both normals, and all retained openings.

## IV. Exhaustive complementary testing and state joins

**Weak-IH:** the source shapes `rho0` interior, zero, one are exhaustive. The double-corner `u,w` infinitesimal/zero is handled in II. For `u` a positive testing unit, `w` infinitesimal/zero, the independent original-parent calculation OA3--OA9 (no contact graph specialization at infinity) gives: interior `rho0` a strictly positive `(P,E)` block provided `a<s/(1+rho0)` and a strictly positive `(V,Z,X)` block provided `s/(s+rho0)<a<s/(s-rho0)`, with `s=sqrt(1+u_0^2)>1`. These intervals intersect. At `rho0=1`, cap the unbounded original `D` event by a sufficiently large real constant; at `rho0=0`, retain it whole and use the explicit two-square sums OA9. The `Y` normal is restored at every opening by a bounded coefficient with zero residue, since `(AH/B)_0=(alpha(1-H)/B)_0`. For testing `v` finite (including `v=0`), the exact original parent and the original `A,alpha` positive rows give the full-state diagonal strict residues IR1--IR3 for testing `u` infinitesimal, unit, unbounded respectively. For both tests unbounded, IR4--IR8 yields a bounded strict Gram; its only nontrivial `(P,E_X)` block has determinant `8` after retaining the original `n,m` packets and any `D` excess. These cases are exhaustive of all nine projective angle types. The formulas use real caps **only after** the particular positive testing residue is fixed; no uniform cap is claimed. Signs are carried by the exact homogeneous reflection symmetries; complex states split into real and imaginary parts.

**Strong-II:** if testing `v` is noninfinitesimal, including literal infinity, the original homogeneous parent RR3 and RR4 yield a bounded strict Gram on `(sqrt(C)E,sqrt(D)J,bL)` with residue `diag(1,1,1/2)`, since `C,D` are unbounded. If testing `v` is infinitesimal/zero but `u` is noninfinitesimal, RR5 decomposes `CD_v E^2+D(E+vK)^2` into a retained positive square plus `gamma V^2` with `gamma=yCD/(C+D)` a unit. The exact toric identity `ax ell-b^2 z_0 ell=bVq` permits the finite real cross-shift `k`. For testing `b_0=1/sqrt(1+u_0^2)>0`, choose the real shift in RR10: its interval is nonempty because `sqrt(1+u_0^2)>1`; the `(x,z_0,ell)` determinant RR11 and `(V,q)` determinant RR12 are strictly positive. For `b_0=0` (testing Alice infinity), take `k=-1` and use the explicit three squares RR13. All potentially negative `v^2` corrections have zero-residue coefficients. The remaining original `n,m,A,alpha` packets are retained as positive squares. The local two-infinitesimal chart is III. Thus **every** projective testing pair, including literal `u=0,v=0` and either infinity, is included. No testing state normalization was imposed; the full four-state identities already cover it.

**Conclusion.** Subject to the independently established original native contact equations (which are used only for the stationary jets and the exact source formulas), both the weak-IH and strong-II C2 sectors have exact full-parent separated tensor-SOS certificates at all independent tests. The decisive bounded-admission ingredients are the weak-IH four-monomial *contact-division lemma*, its two normalized endpoint source cancellations, and the strong-II *four-packet subtraction + bounded dependent-root gauge + monic division + toric `T` pivot*. Scalar determinant nonnegativity alone is nowhere used as a substitute for the unrestricted two-normal Gram.

**Reproducibility:** Run the bundled exact SymPy scripts; they check the generic IC four-monomial residual support, both IC normal expansions, WIE5, exact simultaneous endpoint source limits, LII4/LII6 full-parent identities, the five LII15 conductors and LII20, the LII dependent-root coefficient (29), `phi(0,V)`, and the two residue determinants. These are formal identities; the ordered-field inequalities are proved above. No numerical scans or SDP solver are used.
