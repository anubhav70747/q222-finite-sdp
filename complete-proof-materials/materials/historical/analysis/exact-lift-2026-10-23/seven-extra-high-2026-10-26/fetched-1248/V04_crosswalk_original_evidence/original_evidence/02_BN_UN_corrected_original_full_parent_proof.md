# Q222 C1 C14.4, corrected BN/UN full-parent proof

26 October 2026. Mathematical scope: **the physical C14.4 source cell** with `S,kappa,nu` positive infinitesimals, `ell=1-nu` a positive unit, `q=nu/ell`, the original strictly positive native gaps, `D,beta>=0`, and all permitted limits of `r,1-r,theta,D,beta`. BN has bounded `tau=h/q`; UN has unbounded `tau`, `e=kappa*tau/ell` infinitesimal (as forced by the exact native sphere). Source and test lie in **separately labelled** real closed fields `K_c,K_t`; roots are in `T=K_c tensor_R K_t`. We prove the local matched-Bob SOS with two unrestricted state normals, then join independently scoped existing CE/SE exits. This is **not** a proof of the full Q222 finite lift or a common finite angular feature basis.

Authoritative unchanged input: [pinned c1_short.tex](https://github.com/anubhav70747/q222-finite-sdp/blob/171a7e94b8f52f2bb601562eda362badfb687a9a/latex/certificates/c1_short.tex), specifically ES1–ES10, FM2–FM3, the native-contact identities preceding SN, SN5–SN18, BN, UN, the CE1–CE10 adjacent-current exit and SE10–SE12 unmatched-Bob exit; TF1–TF2 from the tensor chapter. Native positivity and stationary contact are **hypotheses**, verified there from the original coquotients, not from a global C1 SOS assertion. Everything needed for the present assembly is displayed or specified below.

## 1. Source normalizations and independent-rate facts

Write `x=kappa*S`, `sigma=nu*S`, `L=ell*S`, `theta=1-h`, `0<r,h<1`, `zeta=ell/kappa`, and `e=kappa*h/nu=kappa*tau/ell`. Set

\[
\begin{aligned}
J&=1+xh-\nu S(1-r),&M&=1+xr-\kappa\theta/\ell,\\
P_0&=1+S\{\nu\ell+2\nu\kappa h+\kappa\ell+\kappa^2h\theta\}
   +S^2\{\nu r\kappa\ell+\nu r\kappa^2h\},\\
U&=(1+x)(1+\nu Sr),&H_0&=1+\ell S(1-r),\\
B/S&=P_0+\kappa(J-\nu S).
\end{aligned}
\]

All of `J,M,P_0,B/S,H_0` are positive bounded units, independently of the rates of `S,kappa,nu,r,1-r,theta`; `e` is infinitesimal in both BN and UN. Original coquotient substitution gives

\[
2\theta P_0\beta=h(B/S)-\nu H_0,
\qquad
\mu\ell^2=\frac{\tau U}{2MP_0},\quad \mu=m/\kappa^2,
\]

\[
Jd_s+\kappa\theta n_s=\frac{\ell^2\mathcal H}{2P_0},\quad
 d_s=\frac{D}{Sq^2},\quad n_s=\frac\beta{q^2},\quad
\mathcal H\sim1.
\tag{S1}
\]

These are polynomial identities after the positive source denominators are cleared; no inverse of an opening is used in the **final bounded coefficients**. The sphere coordinates `d_0,b_0 >=0`, `d_0^2+b_0^2=1`, give
\[
d_s=\frac{\ell^2\mathcal H d_0^2}{2P_0J},\qquad
\kappa\theta n_s=\frac{\ell^2\mathcal H b_0^2}{2P_0}.
\tag{S2}
\]

**BN:** Put `b=beta/q`, `z=kappa*beta/q^2`. The exact beta identity is
\[
 b=\frac{\tau(B/S)-\ell H_0}{2\theta P_0}.
\]
Since `beta>=0`, `tau >= ell H_0/(B/S)` is bounded below by a positive real; since BN assumes it bounded, `tau,mu` are positive bounded units. Moreover `b,z,d_s` are bounded, `Jd_s+theta*z~1/2`, so either `d_s` or `z` is a unit. A cancellation, to be performed **before** taking residues, gives
\[
\mu\ell^2-b=
\frac{\ell H_0}{2\theta P_0}+
\frac\tau{2P_0}\left(\frac U M-\frac{B/S}{\theta}\right)\sim\tfrac12.
\tag{S3}
\]
No restriction on `beta/q^2` is used. Also `c_0=SC(1-r), a_0=SAr, a_1=alpha/S~1/2`, whereas `m,D,alpha,beta` are infinitesimal.

**UN:** Define `n=beta*theta^2+2Dh+h^2/2>0`, `z_n=kappa*n/q^2`, `A_00=m*zeta^2-n/q`, and `M_*=ell*M`. Substituting ES5–ES7 into this definition gives the **exact** bounded expression
\[
 z_n=\frac{\ell^2\mathcal H\theta b_0^2}{2P_0}
 +\frac{\ell^2(eH_0+h\mathcal H b_0^2)}{2(B/S)}
 +2\kappa S d_s h.
\tag{S4}
\]
Thus `d_s,z_n` are bounded and `d_s+z_n~1/2`, even if standalone `beta` diverges. Write
\[
T=\ell JU+JM_*P_0-JM_*\kappa\theta(J-\nu S)-2M_*U,
\quad T_\nu=r-1-\nu-\nu S(1-r)(r+\ell).
\]
Direct substitution of the printed short source polynomials shows `T|_{kappa=0}=ell*nu*S*T_nu`. Hence the **monic source quotient** `T_kappa=(T-ell*nu*S*T_nu)/kappa` is a bounded polynomial, and
\[
 A_{00}=
 \frac{\theta\ell H_0}{2P_0}
 +\frac{eT_\kappa+h\ell ST_\nu}{2JMP_0},
 \qquad A_{00}-\theta/2\in\mathfrak m_c.
 \tag{S5}
\]
This explicitly cancels the independently unbounded `m*zeta^2` and `n/q`. In particular `n/q~tau/2`, so `q/n` is infinitesimal. For `theta` infinitesimal, `h~1`, `n~1/2`, and `kappa/q=e*ell/h` is infinitesimal. Finally
\[
 t_0=\frac{D+h/2-\beta\theta}{q}
 =\frac Dq+\frac{\ell H_0-\kappa\tau(J-\nu S)}{2P_0}\sim\tfrac12.
 \tag{S6}
\]
In UN literal `beta=0` cannot occur: the beta equality would force `tau=ell H_0/(B/S)` bounded. Literal `D=0` is permitted.

## 2. Complete original parent, including both normals

With `D_u=1+u^2`, `D_w=1+w^2`, and `d=D+1/2`, the **actual reciprocal-Bob quarter-parent** is
\[
\begin{aligned}
\mathcal F={}&m(W-wX+uY-uwZ)^2+C D_w(uW-Y)^2+A D_uD_wY^2\\
 &+d(uwW+uX-wY-Z)^2+(\alpha+1/2)D_u(Y-wZ)^2\\
 &+(\beta+1/2)D_w(X+uZ)^2-\tfrac12D_uD_wZ^2.
\end{aligned}\tag{P1}
\]
Hereafter `a=u^2`, `c=w^2`, `delta=c-q`. In BN use `X_e=X+ewW`, in UN use `N_x=X+huZ`, and in both cases `N_y=(Y-ruW)/sqrt(Sr(1-r))`. These are independent unrestricted state coordinates. The full expansion, for either X-normal, is
\[
\mathcal F=F_0+H_X X_*^2+H_YN_y^2+2H_{XY}X_*N_y+2L_X X_*+2L_YN_y,
\]
\[
\begin{aligned}
H_X&=mc+da+(\beta+1/2)(1+c),\\
H_Y&=rc_0(1+c)+(1-r)a_0(1+a)(1+c)
   +Sr(1-r)[ma+dc+(\alpha+1/2)(1+a)],\\
H_{XY}&=-uw\sqrt{Sr(1-r)}(m+d).
\end{aligned}\tag{P2}
\]
Thus `H_Y~1/2` uniformly across both `r` openings, and `H_XY` has zero residue. BN has `H_X~1/2`; UN subtracts the actual positive packets `da N_x^2` and `beta(1+c)N_x^2`, leaving `(mc+(1+c)/2)N_x^2~N_x^2/2` (see Section 5). Nothing in (P1) or (P2) specializes a normal or state ratio.

The corrected homogeneous radial rows are
\[
\begin{gathered}
 \boxed{N_0=N/\ell=W/\zeta-uwZ},\quad N=\kappa W-\ell uwZ,\\
 G=uW/\sqrt S-\ell\sqrt S wZ,
 \quad B=\sqrt S\,\delta Z,\quad B_u=u\delta Z,\\
 R=w(a-x)Z=\frac{\sqrt S}{\zeta}G-uN_0.
\end{gathered}\tag{P3}
\]
The recovered `N_0=N/kappa` is a **typo**; it contradicts both the definition and the subsequent identities. All identities here and below use `N_0=N/ell`.

## 3. Exact monic quotients and **complete** scalar and normal identities

BN graph (`X_e=N_y=0`) is
\[
 F_0=A_NW^2+(a/S)A_PW^2-2KuwWZ+F_ZZ^2,
\]
\[
\begin{aligned}
 A_N&=m(1+ra+ec)^2+(\beta+1/2)e^2c(1+c),\\
 A_P&=(1+c)[c_0(1-r)+a_0r(1+a)]
   +S[d(1-r-e)^2c+(\alpha+1/2)(1+a)r^2],\\
 K&=m(1+ra+ec)+d(1-r-e)+(\alpha+1/2)(1+a)r+(\beta+1/2)e(1+c),\\
 F_Z&=D+\alpha c+\beta a+(m+\alpha+\beta+1/2)ac.
\end{aligned}\tag{P4}
\]
Define `J_1=K-ell A_P`, `L_1=J_1-zeta A_N` and
\[
H=F_Z-\ell^2ScA_P+ac(\zeta^2A_N-2\zeta J_1).
\]
The original native stationary contact gives `L_1(x,q)=H(x,q)=H_a(x,q)=H_c(x,q)=0`, and the **full** parent gives `H(a,0)=D+beta*a`. Finite monic subtraction first in `a-x`, then `c-q` therefore yields
\[
 L_1=(a-x)L_a+\delta L_c,
\quad H=S\delta^2(d_s+\kappa n_s+cH_{02})+(a-x)\delta(-qn_s+cH_{11})+c(a-x)^2H_{20}.
\tag{P5}
\]
Proof of the specialized factors, without a huge source expansion: `H(x,c)` is divisible by `(c-q)^2` by its two contact jets, and its value at `c=0` is `D+beta*x`, fixing `(d_s+kappa*n_s)`. Likewise `H_a(x,c)` is divisible by `c-q` and its value at `c=0` is `beta`, fixing `-q*n_s`. Finally `H(a,0)` is affine in `a`, so the second divided difference in `a` has a factor `c`. The coefficient of `c^3` in `H(x,c)` gives exactly `H_02=kappa*tau^2*(m+beta+1/2)`. All divisions are monic; no evaluation of `(a-x)^{-1}` or `(c-q)^{-1}` is made.

The source-safe combinations are
\[
\begin{aligned}
 C_R&=H_{20}+2\zeta L_a=-\ell\mu r(a\ell r+2c\kappa\tau+2\ell),\\
 C_G&=\kappa(H_{11}-n_s)/\ell,\\
 C_N&=(n_s-H_{11})/\ell-2L_c/\kappa.
\end{aligned}
\]
The exact scalar identity (not just contact evaluation) is
\[
F_0=(A_N/\kappa^2)N^2+A_PG^2+d_sB^2+n_sB_u^2+H_{02}(wB)^2
 +C_RR^2-2\sqrt S L_aGR+C_GG(wB)+C_NN(wB_u).
\tag{P6}
\]
The conversion uses only `uN=kappa sqrt(S) G-ell R` and
\[
 c(a-x)\delta Z^2=
 \frac{\kappa G(wB)-N(wB_u)}\ell.
\tag{P7}
\]
In particular the original `n_s B_u^2` is not omitted.

For BN subtract the true four-term circuit `Psi` below. Its polynomial effect is exactly
\[
\begin{gathered}
 A=\zeta^2 A_N-b,\quad A'_P=A_P-zc/\ell^2,\quad
 H'_{11}=H_{11}-2n_s,\quad L'_c=L_c+\kappa n_s/\ell,
\end{gathered}
\]
with `b=beta/q`, `z=kappa*beta/q^2`, and `L'_a=L_a`, `H'_{20}=H_{20}`, `H'_{02}=H_{02}`. The whole scalar becomes
\[
\begin{aligned}
F_0={}&\Psi+A N_0^2+A'_PG^2+d_sB^2+H_{02}(wB)^2+C_RR^2
 -2\sqrt S L_aGR\\
 &+C'_GG(wB)+C'_NN_0(wB_u),\\
 C'_G&=\kappa(H_{11}-2n_s)/\ell,\\
 C'_N&=-Sa_0\ell^2r+(\beta+m+1/2)\tau^2(c+1/\ell)
            +2\ell\kappa\mu\tau-m\tau^2.
\end{aligned}\tag{P8}
\]
The cancellation in `C'_N` follows by extracting coefficients of `c^2,c^3` in `H_a(x,c)` and using `ell(1+q)=1`, not by bounding `H_11` and `L_c/kappa` separately. `A~A'_P~1/2` by (S3). Moreover `C_R,L_a,C'_G,C'_N` are bounded; all non-helper cross coefficients have zero residue after the exact `R` identity in (P3).

The **complete BN first normals**, not graph-only terms, are
\[
\begin{aligned}
L_X={}&\ell w\widehat A_XN_0+\frac{uw\sqrt S}\ell\widehat X_aG
  +d_s\sqrt S quB-bB_u+wX_c(wB_u),\\
L_Y={}&\sqrt{r(1-r)}[Y_PG+\sqrt S Y_aR+Y_c(wB)],
\end{aligned}\tag{P9}
\]
where `widehat A_X=-kappa*mu-(tau/ell)[(beta+1/2)(1+c)+mc]`, `widehat X_a=ell[d(1-r-e)-mr]`, `X_c=-tau(m+beta+1/2)`, and `Y_P=(1+c)(-c_0+a_0(1+a))+S[m(1+ra+ec)-dc(1-r-e)+(alpha+1/2)(1+a)r]`. The monic quotients `Y_a,Y_c` of `Y_C=ell Y_P-ma+d-(alpha+1/2)(1+a)` are bounded, as is `Y_P` (indeed infinitesimal). To check the first row, the raw X half-derivative is `w A_X W+u T_X Z` with
\[
 T_X=(\ell/\kappa)c A_X+\beta-D+c(m+\beta+1/2).
\]
Native contact and `T_X(a,0)=beta-D` give by monic division
`T_X=c(a-x)X_a+delta[(D-beta)/q+c X_c]`. Use (P3), and the *exact*, face-safe identity
\[
  (D-\beta) B_u/q=d_s\sqrt S quB-n_s qB_u.
\]
The Y row follows by dividing `Y_C(x,q)=0` monically. No independently unbounded `n_s` is admitted into a Gram coefficient.

**UN:** Retain the different actual normal `N_x=X+huZ`. Its exact scalar is (P4) with
\[
\begin{aligned}
 A_N&\rightsquigarrow m(1+ra)^2,\\
 A_P&=(1+c)[c_0(1-r)+a_0r(1+a)]
      +S[d(1-r)^2c+(\alpha+1/2)(1+a)r^2],\\
 K&=m\theta(1+ra)+d(1-r)(1+ha)+(\alpha+1/2)(1+a)r,\\
 F_Z&=D+\alpha c+na+[m\theta^2+\alpha+(\beta+1/2)\theta^2]ac+dh^2a^2.
\end{aligned}
\]
Subtract the circuit `Psi_n`; the resulting `H'` has `H'(a,0)=D+dh^2a^2`, has a second-order stationary zero at `(x,q)`, and is quadratic in `c`. Monic extraction fixes every coefficient in
\[
H'=(D+dh^2x^2)\delta^2/q^2
 +(a-x)\delta(-2dh^2x/q+cH_{11})+(a-x)^2(dh^2+cH_{20}).
\tag{P10}
\]
The coefficient of `c^2` in `\partial_aH'(x,c)` is `H_11`. Extracting it and the secant of `L'` gives directly, *without native source substitution*,
\[
\begin{aligned}
C_R&=H_{20}+2\zeta L_a=-m\zeta^2r(2+ra),\\
C_G&=H_{11}/\zeta=2\ell P_1(x)+\ell xa_0r-2z_n/\ell,\\
C_N&=-H_{11}-2\zeta L_c=-Sa_0\ell^2r,
\end{aligned}\tag{P11}
\]
where `P_1(a)=[c]A_P`. **Do not bound `C_R` separately**: use `R=(sqrt(S)/zeta)G-uN_0` and the *exact* `A+a C_R=A_00`. The complete scalar is
\[
\begin{aligned}
F_0={}&\Psi_n+A_{00}N_0^2+\widetilde A_PG^2
  +2u\sqrt S(L_a-C_R/\zeta)N_0G\\
 &+d_sB^2+dG_h^2+C_GG(wB)+C_NN_0(wB_u),\\
\widetilde A_P&=A'_P+SC_R/\zeta^2-2SL_a/\zeta\sim1/2,\\
G_h&=h(a-xc/q)Z.
\end{aligned}\tag{P12}
\]
`C_R/zeta=-m*zeta*r(2+ra)` is infinitesimal, `C_R/zeta^2=-mr(2+ra)` bounded; `C_G` bounded and `C_N` infinitesimal. All terms in (P12) have bounded coefficient arrays for any independent source/test rates.

The **complete UN X-normal** has, before its exact positive exchange,
\[
 L_X=-m\zeta wN_0+d\tau a(wN_0+B_u)
 +uw\sqrt S[d(1-r-e)-mr]G+t_0B_u.
\]
Set `Q=wN_0+B_u=wW/zeta-uqZ`. Then
\[
 G_h=-\tau uQ+ew\sqrt S G,\qquad N_x+\tau Q=X+ewW,
\]
\[
 dG_h^2+daN_x^2+2d\tau a QN_x
 =da(X+ewW)^2-2de\tau uw\sqrt S\,QG+de^2Sw^2G^2.
\tag{P13}
\]
Retain the complete positive squares `da(X+ewW)^2`, `beta(1+c)N_x^2`, and optionally the last nonnegative term. The remaining X hard diagonal is `mc+(1+c)/2~1/2`; the remaining signed QG coefficient has zero residue since `e*tau` is bounded. The only surviving non-small X coupling is bounded `t_0 B_u`. The UN Y row is (P9) with
\[
Y_P=(1+c)(-c_0+a_0(1+a))
+S[m(1+ra)-d(1-r)c+(\alpha+1/2)(1+a)r]
\]
and `Y_C=ell Y_P-m*theta*a+d(1+h*a)-(alpha+1/2)(1+a)`, divided monically as above. This includes the *entire* mixed `H_{XY}` of (P2).

## 4. BN: positive packets and finite ordered Gram completion

On the `u,w>=0` chart, put
\[
(A_1,A_2,A_3,A_4)=\sqrt\beta\,(uZ,\ W/(\zeta\sqrt q),\ \sqrt{x/q}\,wZ,\ uwW/(\zeta q\sqrt x)),
\]
\[
\Psi=\sum_iA_i^2-4\beta uwWZ/(\zeta q),\quad E_1=A_1-A_2,\ E_2=A_3-A_4.
\]
For `WZ>=0`, two geometric-mean squares give `Psi=E_1²+E_2²+2(sqrt(A_1A_2)-sqrt(A_3A_4))²`; for `WZ<0`, `Psi-1/2(E_1²+E_2²)=1/2(A_1+A_2)²+1/2(A_3+A_4)²+4beta|uwWZ|/(zeta q)`. Every positive radical is a product of pure-source and pure-testing roots. Therefore `Psi-1/2(E_1²+E_2²)` is tensor SOS, also when `beta=0` (omit the zero rows).

For `d_s` a unit, `wB=(w/sqrt(d_s))sqrt(d_s)B`, a zero-residue-cost helper. Otherwise `z=kappa*beta/q²` is a unit and
\[
 wB=-\frac{w+\sqrt q}{\sqrt z}E_2-\frac{w(w+\sqrt q)}\ell G.
\tag{B1}
\]
Both coefficients are infinitesimal. **Only after buying `wB`**, use the exact homogeneous zero identity
\[
 \delta^2N_0^2+\frac{2\kappa\delta}\ell G(wB)
 +2\kappa(wB)^2-(wB_u)^2=(\delta W/\zeta)^2
\tag{B2}
\]
to buy any prescribed real pivot on `wB_u`, at zero-residue cost on `N_0,G,wB`. Finally the remaining BN X coupling is paid by
\[
 bB_u=-(w+\sqrt q)\sqrt b\,E_1-b(w+\sqrt q)N_0,
\tag{B3}
\]
which uses no inverse of `beta` and whose coefficients are infinitesimal. The only bounded non-small off-diagonal terms are incident to the purchased helpers (`C'_G G(wB)`, `C'_N N_0(wB_u)` and `Y_c N_y(wB)`).

For each fixed pair of independently labelled fields and each actual source/test, retain a positive fixed fraction of each available circuit square, the `d_sB²` packet if nonzero, and the four strict main rows `N_0,G,X_e,N_y`. Their leading diagonal is `(1/2,1/2,1/2,1/2)` and all mutual off-diagonals have zero residue. Buy `wB` **first**, then `wB_u`; add sufficiently large **real** diagonal budgets to dominate all bounded helper couplings by an ordinary real Schur complement. The cost on preceding Gram entries always has residue zero. The resulting *full* coefficient Gram has positive-definite real residue and bounded entries. TF2 gives an exact tensor SOS for (P1), retaining the circuit remainder and every omitted positive original packet. This is independent of all relative infinitesimal rates.

## 5. UN: source-safe cancellation, positive exchange and ordered helpers

Take the same four `A_i` with `sqrt(beta)` replaced by `sqrt(n)` and define `E_3=A_1-A_3`, `Psi_n=sum A_i²-4n uwWZ/(zeta q)`. Averaging the pairings `12|34` and `13|24` gives
\[
 \Psi_n-\tfrac14(E_1^2+E_2^2+E_3^2)\in\Sigma T^2.
\tag{U1}
\]
For `WZ<0` the residual quadratic matrix on `(A_1,...,A_4)` contains `I-(1/4)L_{2-1-3-4}\succ0`, since the path Laplacian has largest eigenvalue `2+sqrt(2)<4`; its positive signed cross packet is a separated square. There is no sign restriction on the state.

1. Buy `wB` using `sqrt(d_s)B` if `d_s` is a unit, otherwise `z_n` is a unit and use (B1) with `z_n,E_2` in place of `z,E_2`.
2. If `theta` is a unit, `A_00~theta/2` already gives a strict `N_0` pivot. If `theta` is infinitesimal, `n~1/2`, `q/n` and `kappa/q` are infinitesimal, and the exact relation
\[
N_0=\sqrt{q/n}(E_2-E_1)-\frac{w-\sqrt q}{\sqrt n}E_3
       +\frac{\sqrt{\kappa/q}}\ell wG
\tag{U2}
\]
buys any fixed real `N_0²` pivot at zero-residue cost. The formerly possibly signed bounded `A_00` has residue zero and is then harmless.
3. With `N_0` strict, use
\[
 B_u=-(w+\sqrt q)\sqrt{q/n}\,E_1-(w+\sqrt q)N_0
\tag{U3}
\]
to buy any fixed real `B_u²` at zero-residue cost. Finally `wB_u=w B_u` has infinitesimal cost. **The order is noncircular**: the purchases rely only on already strict circuit/main roots.

After (P13), the hard normal residues on `N_x,N_y` are `(1/2,1/2)`, and the strict scalar residue on `G` is `1/2`. `N_0` is strict either from `theta/2` or (U2). The terms with possibly nonzero bounded cross coefficients, `C_G G(wB)` and `t_0 N_xB_u`, are absorbed by finite real helper pivots. Every remaining cross coefficient has zero residue (including `e*tau*u*w*sqrt(S)`, `m*zeta*w`, `C_N`, and the entire `Y` row). The retained `Psi_n` remainder, `da(X+ewW)²`, `beta(1+c)N_x²` and positive exchanges are independent separated SOS. TF2 proves the **complete** (P1) is tensor SOS, including `theta` infinitesimal and literal `D=0`, with every opening rate allowed.

## 6. All tests and the exact dependency boundary

The literal local axes `u=0,w=0` require no division. For negative angle signs use the exact physical sign gauges, then separate the two signs of `WZ` as above. For arbitrary complex states apply the real quadratic identities separately to real and imaginary components. No state coordinate has been normalized or bounded.

Since `S` is infinitesimal and `c_0,a_0~1/2`, the *actual original positive events* satisfy `C=c_0/[S(1-r)]` and `A=a_0/(Sr)`, hence both are **unbounded**, not just positive. The independent adjacent-C,A theorem CE1–CE10 applies to **noninfinitesimal Alice tests** with the other events unselected and nonnegative; its required positive-real lower bounds follow from unboundedness. In this C14.4 sector `e` and `m` are infinitesimal and the source equalities imply SE8. The independently scoped unmatched-Bob theorem SE10–SE12 therefore covers every normalized homogeneous Bob sine-square `s²` with `res(s²)!=res(ell)=1`, including Bob zero and infinity. If `res(s²)=1`, the complementary cosine is infinitesimal and Bob's reciprocal coordinate `w=cos/sin` is infinitesimal. When Alice is also infinitesimal, Sections 4–5 cover this remaining *entire* chart. These two exits are imported **only** with their stated source and angle hypotheses; neither imports global C1 SOS. Thus the BN/UN source cells are exhaustively certified at all independent projective tests **provided those separately established CE/SE theorems are available**.

This closes BN/UN as a component of the qualitative tensor proof, but it makes **no** claim that the other native source sectors are completed or that an explicit, source-affine finite SDP lift has been extracted.

## Reproducibility record

The three accompanying self-contained Python/SymPy checkers are exact, no random evaluations or optimizers:
- `verify_bn_un_independent.py`: original ES/FM source-sphere identities, UN3/UN4 source cancellations and UN15 pole, corrected N0 and radial/circuit/helper identities, UN16 exchange.
- `verify_bn_un_parent_full.py`: twelve coefficient identities comparing the original seven-event parent with BN and UN scalar, both full normal diagonals, mixed normal block, and both linear normals.
- `verify_bn_un_structural_coeffs.py`: SN12 radial cancellation, SN13 coefficient formula, H02, BN10 shift, BN12 pole-free coefficient, and UN13 three exact cancellation identities. It uses coefficient extraction before source substitution, and the identity `ell(1+q)=1` when required.

**Careful status:** The source-rate inequalities are derived above from the explicit native source/sphere identities, whose exact polynomial formulas and independent native-contact verification are inputs. They were not replaced by an approximate one-parameter expansion. The fixed test feature library and the complete Q222 lift are outside this component's scope.
