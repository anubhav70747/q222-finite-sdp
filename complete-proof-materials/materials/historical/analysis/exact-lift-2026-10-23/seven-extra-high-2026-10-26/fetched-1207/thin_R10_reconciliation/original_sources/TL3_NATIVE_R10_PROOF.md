# Exact TL3 reconstruction and the original R10 bounded-\(b\) gate

**Scope.** This note proves the seven coefficient identities in TL3, the exact formal same-contact C1 dictionary, the signs and scale of \(\pi,\beta_1,\kappa\), and the bounded-\(b\) premise, directly from native C2 and C1 source formulas. All Bell parents use **\(F/4\)**. It proves an algebraic dependency of TL6, **not** the TL6 conductor theorem or a global Q222 finite SDP lift. Source and testing parameters remain independent.

**Source locations.** The supplied original excerpts give BA1, R10.2, OP1–OP5 and TL1–TL3. The pinned repository contains the original C2 formula in `latex/certificates/c2_short.tex`, section *Stationary coefficients and contact*, and native C1 identities NC1–NC4 in `latex/certificates/c1_short.tex` (roughly lines 195–240): [C2 native](https://github.com/anubhav70747/q222-finite-sdp/blob/ada9c64d79dc4e50f862b154e1743374c5f9bbe8/latex/certificates/c2_short.tex) and [C1 native](https://github.com/anubhav70747/q222-finite-sdp/blob/ada9c64d79dc4e50f862b154e1743374c5f9bbe8/latex/certificates/c1_short.tex).

## 1. Independent source dictionary

In the actual weak-C2 R10 source, let

\[
\begin{gathered}
 a>0,\quad B>0,\quad t>0,\quad 0<h<1,\quad 0<r<1,\\
 e=1-r,\quad s=a^{-1},\quad w=ae,\quad b=a^2B,\quad x=aB,\quad y=tB,\quad\ell=1-h,\\
 K=1+x,\quad c_a=1+xr,\quad S_c=1+B\ell+xr,\\
 Q_c=1-hB,\quad E=1-txr,\quad R=1+xh,\\
 g=R-w=1-ae+hx,\qquad j_R=c_a-\ell/t.
\end{gathered}\tag{1}
\]

The original strict source gaps ensure \(E,Q_c,g,j_R>0\). For the thin-R10 valuation, \(s,B,h,e,x\) are positive infinitesimals, \(a\) unbounded, and \(t,g,j_R\) positive bounded units; \(\rho=n/e^2>0\) is infinitesimal. **No bound on \(b\) is added.**

### Deriving BA1 from the original native C2 formula

Let \((a_0,z_0,b_0,t_0,r_0)\) denote the parameters of the *original* native C2 stationary formulas (where \(z_0=c_0^2\)). Their exact dictionary is

\[
\boxed{a_0^2=t/a,\quad z_0=x,\quad b_0=a^{-1},\quad t_0=\ell/t,\quad r_0=r.}\tag{2}
\]

The original gaps, determinant denominator, and positive C2 numerator become

\[
\begin{array}{c|c}
\text{original C2 quantity}&\text{weak-R10 expression}\\\hline
 b_0-a_0^2t_0&h/a\\
 g_0=(1+z_0)b_0-(1-r+a_0^2z_0t_0)&g/a\\
 j_0=1+z_0r-t_0&j_R\\
 1-z_0(b_0-a_0^2t_0)&Q_c\\
 b_0-a_0^2z_0r&E/a\\
 1+z_0r+a_0^2z_0t_0&S_c\\
 \Delta_0&\Delta/a\\
 N_0^{\mathrm{C2}}&M/a^2.
\end{array}\tag{3}
\]

For example, the last equality follows without expansion:

\[
\begin{aligned}
 a^2 N_0^{\mathrm{C2}}
 &=a^2\left[(1+x)\frac{1+ar}{a}\frac{h+t}{a}
 -\frac ta(1+tB)e\frac\ell t\right]\\
 &=K(1+ar)(h+t)-a(1+y)e\ell=:M.
\end{aligned}
\]

The denominator identity \(a\Delta_0=\Delta\) follows by inserting (2) into its six printed terms and collecting powers of \(B\). Explicitly,

\[
\boxed{
\Delta=a+(1+x)+(1+x)B+2x(1+x)r+a^3B^2r^2
       +2a^2Br+tB[1+a+x+a^2Br(2-r)].
}\tag{4}
\]

Thus, in the common normalization \(d_2=1\), \(\lambda_2=\beta_2=0\), the original native C2 coefficients become, with \(N:=2r\Delta n\),

\[
\boxed{\begin{aligned}
 N&=(1+ar)c_aS_c-r\Delta,\qquad n=N/(2r\Delta),\\
 m_2&=\frac{a^2(1+y)e c_aS_c}{2g\Delta},&
 C_2&=\frac{c_a M}{2\ell\Delta},\\
 (Ah)_2&=\frac{t c_a(1+ar)S_c}{2\Delta},&
 D_2&=\frac{K(1+ar)ES_c}{2t j_R\Delta},\\
 \alpha_2&=\frac{B^2(1+ar)M}{2Q_c\Delta}.
\end{aligned}}\tag{5}
\]

The entire source contact is
\[
(u_c,v_c;W_c,X_c,Y_c,Z_c)
=\left(\sqrt{x},\sqrt{tB};\frac{B}{\sqrt{x tB}},
-r\sqrt{x},\frac{hB}{\sqrt{tB}},1\right).
\tag{6}
\]
No testing quantity has been set equal to this contact in any coefficient identity below.

### Reconstructing the formal same-contact C1 formulas

The map OP4 is
\[
\boxed{x_1=x,\quad \eta_1=t,\quad \zeta_1=s,\quad r_1=h,
\quad\theta_1=e,\quad t_1=\theta_1/\zeta_1=w,
\quad a_1=\ell/t.}\tag{7}
\]
It gives the same full contact (6). In the original C1 NC1 notation,
\[
\begin{gathered}
 y_1=x\eta_1\zeta_1=y,\quad K_1=K,\quad R_1=R,\quad
 T_1=c_a,\quad M_1=g,\quad J_1=j_R,\\
 E_1=E,\quad Q_1=Q_c,\quad G_1=r+s,\quad
 h_{\mathrm{C1,NC1}}=1-\theta_1=r.
\end{gathered}
\]
**Warning:** the local centered-state \(h_1=a_1/T_1=\ell/(t c_a)\) is neither the R10 infinitesimal \(h\) nor the C1 NC1 letter \(h_{\mathrm{C1,NC1}}=r\). Likewise \(j_R\) is not the CA/SE source coordinate \(j_s=\delta_{\rm CA}k\). The operator \(\delta L=L_2-\kappa L_1\) below is unrelated to the original positive C2 chart denominator \(\delta_{\rm C2}=B\Delta\) or the positive CA source scale \(\delta_{\rm CA}\).

Substitution in C1 NC2 gives the **reconstructed** short polynomials

\[
\boxed{\begin{aligned}
 P_d&=1+B(1+t)+2xr+xB(h+t)+x^2r(h+te),\\
 P_4&=c_a E+x\ell(r+s),\\
 P_{\alpha,1}&=Rr(1+y)+s g E.
\end{aligned}}\tag{8}
\]

For clarity, NC2 originally reads \(P_d=p_0+\eta_1p_1\), with
\(p_0=1+2x(1-\theta_1)+r_1x^2G_1+x\zeta_1\) and
\(p_1=x[K\zeta_1+x(1-\theta_1)\theta_1]\). Inserting (7) gives the first line of (8) by \(xs=B\). The other two are exactly NC2's factorizations \(P_4=T_1E_1+x(1-r_1)G_1\) and \(P_{\alpha,1}=R_1(1-\theta_1)(1+y_1)+\zeta_1M_1E_1\).

NC4, in the normalization \(d_1=1\), now gives

\[
\boxed{\begin{aligned}
 m_1&=\frac{Kx(h+t)r(1+y)}{2gsP_d},&
 C_1&=\frac{K(h+t)E}{2\ell P_d},\\
 (Ah)_1&=\frac{tP_4}{2P_d},&
 D_1+\frac12&=\frac{K(1+y)(h+t)}{2t j_RP_d},\\
 \alpha_1&=\frac{x^2s(h+t)P_{\alpha,1}}{2Q_cP_d},&
 \beta_1+\frac12&=\frac{Kx(h+t)(r+s)}{2eP_d},\\
 n_1&=\lambda_1=0.
\end{aligned}}\tag{9}
\]

In particular the D row uses \(B_5=K(1+y)(h+t)-t j_RP_d\), so the \(-1/2\) in \(D_1\) has **not** been discarded. Formula (9) is an algebraic *formal* C1 source formula, not an invocation of C1 physicality when \(\beta_1<0\).

## 2. Three factored native numerator bridges

The following three independently reconstructed numerators all equal the *same* \(N\):

\[
\boxed{\begin{aligned}
 N&=(1+ar)c_aS_c-r\Delta,\\
 &=eP_d-KB(h+t)(1+ar),\\
 &=P_4-rP_d,\\
 &=eS_cE-rBM.
\end{aligned}}\tag{10}
\]

This equality is elementary enough to check without a computer. Substitute \(x=aB\), \(e=1-r\), \(\ell=1-h\) into **each** of the four expressions in (10). All four are quadratic in \(B\), with precisely the same factored coefficient table:

| coefficient | common value |
|---|---|
| \([B^0]\) | \(e\) |
| \([B^1]\) | \(e-h-rt+ar(2e-h-t)\) |
| \([B^2]\) | \(-ar\big[(h+t)+ar(h+t(2-r))\big]\) |

Therefore their difference is identically zero, over any coefficient field, without a division by \(h,B,e\) or by an infinitesimal. This simultaneously proves the OP2 source equality \(eP_d-KB(h+t)(1+ar)=2r\Delta n\) from **BA1**, rather than importing it as a new assumption.

Three other printed TL3 polynomials have unexpectedly short forms:

\[
\boxed{\begin{aligned}
 P_\kappa&=c_aS_c+KB(h+t),\\
 P_D&=c_a t j_R-K(h+t),\\
 P_\alpha&=(r+s)c_a-eQ_c.
\end{aligned}}\tag{11}
\]

To verify (11) against TL3's displayed definitions, use \(B=sx\):

\[
\begin{aligned}
 c_aS_c+KB(h+t)
 &=c_a^2+sx\big[1+t+x(r+he+t)\big],\\
 c_a t j_R-K(h+t)
 &=-c_aE-xe(h+t)\\
 &=tx^2r^2+x\{t-1+e(1-h-2t)\}-1,\\
 (r+s)c_a-eQ_c
 &=1+s-2e+x\{r^2+s(r+he)\}.
\end{aligned}\tag{12}
\]

For the final two rows of TL3, two small bridges suffice:

\[
\boxed{\begin{aligned}
 P_d&=c_a(1+y)+B(1+ar)(KE+xtj_R),\\
 c_aP_{\alpha,1}&=EP_\alpha+rBM.
\end{aligned}}\tag{13}
\]

The first follows immediately from (8) using
\(KE+xtj_R=1+x(h+te)\) and \(tr+te=t\). For the second, note the explicitly divisible identities

\[
\begin{gathered}
 R(1+y)-E=B\big[ah+t+atr+atBh\big],\\
 eQ_c+sc_a(g-1)=rB(hK-ae),\\
 M=c_a\big[ah+t+atr+atBh\big]+E(hK-ae).
\end{gathered}\tag{14}
\]

The last equality is obtained by expanding the printed \(M=K(1+ar)(h+t)-a(1+y)e\ell\) and collecting \(1,B,B^2\). Multiplying the first line by \(rc_a\), the second by \(E\), and adding proves
\[
 c_aP_{\alpha,1}-EP_\alpha
 =rc_a[R(1+y)-E]+E[sc_a(g-1)+eQ_c]
 =rBM.
\]
This gives a short numerator derivation, not a bare CAS assertion.

## 3. Exact slack, sign, and normalization

The original C2 OP source denominator is \(\delta_{\rm C2}=B\Delta\) and \(P_n=BP_d\). Its OP3 positive weight becomes

\[
 \pi=q=\frac{K(h+t)(B+xr)}{eP_d}
 =\frac{KB(h+t)(1+ar)}{eP_d}.
\tag{15}
\]

The last entry of (9) agrees because \(x(r+s)=B(1+ar)\). By the second equality in (10),

\[
\boxed{\begin{gathered}
 s_n=1-\pi=\frac{N}{eP_d}=\frac{2r\Delta n}{eP_d},\\
 \beta_1=-\frac{s_n}{2}=-\frac{r\Delta n}{eP_d}<0,\\
 1+2n=\frac{(1+ar)c_aS_c}{r\Delta},\\
 \kappa:=\frac{1+2n}{\pi}
 =\frac{e c_aS_cP_d}{r\Delta KB(h+t)}.
\end{gathered}}\tag{16}
\]

The sign is strict in the thin layer \(n>0\), while \(\pi>0\) follows from the positive product in (15). The \(n=0\) boundary gives \(\beta_1=0\), but is not silently substituted for \(n>0\).

Write \(\Theta=n/B=N/(2r\Delta B)\). From (10)–(11):

\[
\begin{aligned}
 e c_aS_cP_d-r\Delta KB(h+t)
 &=c_aS_c\{KB(h+t)(1+ar)+N\}\\
 &\quad -\{(1+ar)c_aS_c-N\}KB(h+t)\\
 &=N[c_aS_c+KB(h+t)]=NP_\kappa.
\end{aligned}
\]

Hence the exact scaling difference is

\[
\boxed{\kappa-1=\frac{NP_\kappa}{r\Delta KB(h+t)}
 =\frac{2\Theta P_\kappa}{K(h+t)}.}\tag{17}
\]

## 4. All seven TL3 rows with numerator derivations

Set \(\delta L=L_2-\kappa L_1\), with \(d_2=d_1=1\). The following cancellations use only (5), (9), and (10)–(17).

**(i) \(\delta m=0\).**
\[
 \kappa m_1
 =\frac{ec_aS_c}{r\Delta KB(h+t)}\frac{Kx(h+t)r(1+y)}{2gs}
 =\frac{a^2(1+y)ec_aS_c}{2g\Delta}=m_2,
\]
where \(x/(Bs)=a^2\).

**(ii) \(\delta C\).** By the fourth equality of (10),
\[
 \delta C=\frac{c_a}{2\ell r\Delta B}
 (rBM-eS_cE)
 =-\frac{c_aN}{2\ell r\Delta B}
 =\boxed{-\Theta c_a/\ell}.
\]

**(iii) \(\delta(Ah)\).** The second and third equalities in (10) give
\(eP_4=rKB(h+t)(1+ar)+N\). Therefore
\[
\begin{aligned}
\delta(Ah)
 &=\frac{tc_aS_c}{2r\Delta KB(h+t)}
   [rKB(h+t)(1+ar)-eP_4]\\
 &=-\frac{tc_aS_cN}{2r\Delta KB(h+t)}
 =\boxed{-\Theta\frac{tc_aS_c}{K(h+t)}}.
\end{aligned}
\]

**(iv) \(\delta D\).** An efficient route proves the *paired* formula first:

\[
 \boxed{D_2+\tfrac12-\kappa(D_1+\tfrac12)
 =-\Theta\left(B+\frac{S_c}{tj_R}\right).}\tag{18}
\]

Multiply the difference between the two sides of (18) by the positive pure-source product \(2r\Delta Btj_R\). After (5) and (9) the result is

\[
\begin{aligned}
&S_c\{rBK(1+ar)E-ec_a(1+y)+N\}
 +Btj_R(r\Delta+N)\\
&=S_c\Big\{e[P_d-c_a(1+y)]
   +B(1+ar)[rKE+c_atj_R-K(h+t)]\Big\}\\
&=S_c B(1+ar)\Big\{e(KE+xtj_R)-e(KE+xtj_R)\Big\}=0,
\end{aligned}
\]

where we used (10), (13), \(r\Delta+N=(1+ar)c_aS_c\), \(E+tj_R=h+t\), and \(K-c_a=xe\). Thus (18) is an exact numerator factorization. Finally add \((\kappa-1)/2\) to (18) and insert (11),(17):

\[
\begin{aligned}
 \delta D
 &=\Theta\left[\frac{P_\kappa}{K(h+t)}-B-\frac{S_c}{tj_R}\right]\\
 &=\Theta S_c\left[\frac{c_a}{K(h+t)}-\frac1{tj_R}\right]
 =\boxed{\Theta\frac{S_cP_D}{K(h+t)tj_R}}.
\end{aligned}
\]

**(v) \(\delta\alpha\).** Use \(S_c=Q_c+x(r+s)\), \(K-c_a=xe\), the second bridge (13), and \(P_\alpha=(r+s)c_a-eQ_c\):

\[
\begin{aligned}
 ec_aS_cP_{\alpha,1}-rKB(r+s)M
 &=eS_cE P_\alpha+rBM\{eS_c-K(r+s)\}\\
 &=eS_cE P_\alpha-rBM P_\alpha
 =NP_\alpha.
\end{aligned}\tag{19}
\]

Consequently
\[
\begin{aligned}
 \delta\alpha
 &=\frac{rKB^2(1+ar)M-ec_aS_c xP_{\alpha,1}}{2r\Delta KQ_c}\\
 &=-\frac{xNP_\alpha}{2r\Delta KQ_c}
 =\boxed{-\frac{nxP_\alpha}{KQ_c}},
\end{aligned}
\]
where \(B(1+ar)=x(r+s)\) supplies the factor \(x\).

**(vi) \(\delta d\).** Both normalizations have \(d=1\). Equation (17) gives
\[
 \boxed{\delta d=1-\kappa
 =-\frac{2\Theta P_\kappa}{K(h+t)}}.
\]

**(vii) The paired \(\beta,d\) row.** Here \(\beta_2=0\) and \(\beta_1+1/2=\pi/2\). Hence without expansion,
\[
\boxed{\delta\beta+\frac{\delta d}{2}
 =\frac12-\kappa\left(\beta_1+\frac12\right)
 =\frac{1-\kappa\pi}{2}=-n.}
\]

Collecting the results gives precisely the original TL3 table, and in the stipulated \(F/4\) normalization:

\[
\boxed{\begin{array}{c|l}
\delta m&0\\
\delta C&-\Theta c_a/(1-h)\\
\delta(Ah)&-\Theta tc_aS_c/[K(h+t)]\\
\delta D&\Theta S_cP_D/[K(h+t)tj_R]\\
\delta\alpha&-nxP_\alpha/(KQ_c)\\
\delta d&-2\Theta P_\kappa/[K(h+t)]\\
\delta\beta+\delta d/2&-n.
\end{array}}\tag{TL3}
\]

## 5. Original R10.2 forces bounded \(b\) (no additional premise)

In R10, \(r=1-sw\), \(B=s^2b\), \(x=sb\), and
\[
\boxed{2r\Delta_*\,\frac{n}{s^2}=wQ_*-bK_*,\qquad
\Delta_*=s\Delta,\quad K_*=(1+s)(h+t)(1+x).}\tag{20}
\]

The exact same identity is obtained by inserting these substitutions into the BA1 numerator \(N=2r\Delta n\):
\(N=s(wQ_*-bK_*)\). It is not an extra hypothesis. Although \(b\) was initially unrestricted, \(Q_*\) and \(\Delta_*\) admit the **exact** expressions

\[
\begin{aligned}
Q_*={}&1+sx(1+t)+x(h+2r+t)+sx^2(h+t)\\
 &+x^2\{h(r+1)+t(1+r-r^2)\},\\[1mm]
\Delta_*={}&1+s(1+x)+s^2x(1+x)+2sx(1+x)r+x^2r^2+2xr\\
 &+t\{s^2x(1+x)+sx+sx^2r(2-r)\}.
\end{aligned}\tag{21}
\]

These contain only the *already assumed* small variables \(s,x,h\), bounded \(r\), and unit \(t\): there is no standalone \(b\) and no \(b^{-1}\). Consequently \(Q_*\sim\Delta_*\sim1\), and \(K_*\sim t\) is a positive bounded unit. In addition, \(g=1-w+hx>0\) gives \(0<w<1+hx\), so \(w\) is bounded; since \(g\) is a positive unit and \(hx\) infinitesimal, \(1-w\) has positive unit residue.

The original physical condition \(n\ge0\), applied **before** the thin specialization, gives from (20)

\[
\boxed{0<b\le w\frac{Q_*}{K_*},\qquad b=O(w),\quad b\text{ bounded}.}\tag{22}
\]

For the thin case \(n=\rho e^2=\rho s^2w^2\), with positive infinitesimal \(\rho\), dividing (20) by the positive source \(w\) gives exactly

\[
\boxed{\frac bw=\frac{Q_*-2r\Delta_*\rho w}{K_*}\sim\frac1t,
\qquad\frac wb\sim t.}\tag{23}
\]

Thus the claimed alternatives \(tb\) a positive bounded unit or infinitesimal are **exhaustive**. Indeed \(tb\) has the same valuation type as \(w\); no unbounded \(tb\) case remains. Using \(e=sw\), \(\Delta=\Delta_*/s\), and \(P_d\sim1\), (16) becomes

\[
\boxed{s_n=\frac{2r\Delta_*\rho w}{P_d}>0,\quad
\beta_1=-\frac{r\Delta_*\rho w}{P_d}<0,\quad
\Theta=\frac{n}{B}=\frac{\rho w^2}{b}=O(\rho w),
\quad\kappa\sim1.}\tag{24}
\]

In particular \(s_n/w=O(\rho)\), \(\beta_1\) is negative infinitesimal, \(\pi=1-s_n\) is a positive unit, and \(\kappa-1=O(\rho w)\) with no comparison to testing scales. All TL3 denominators \(K,\ell,h+t,Q_c,t,j_R,P_d\) are source-field positive units in R10; the large \(A_1\) itself is never divided by \(h\) after centering \(Ah\).

**Boundary qualification.** Formula (22) also holds at \(n=0\), where OP4 yields \(\beta_1=0\). Formula (23) uses \(\rho>0\) infinitesimal and is not asserted for a general nonthin positive \(n\). None of these calculations treats the formal negative-\(\beta_1\) C1 tuple as a physical C1 source.

## 6. Verification and scope

The companion `verify_TL3_exact.py` reconstructs the original C2 stationary formulas through the dictionary (2), reconstructs C1 NC2/NC4 through (7), verifies all factored numerator identities, all seven TL3 rows, and the literal R10.2 identity. It checks **46 exact symbolic equalities** and uses no solver or sampling. The human proof is the coefficient table in section 2, the short bridges (13)–(14), and the numerator chains in section 4. Positivity/valuation in section 5 is a separate consequence of the stated strict source gaps and original R10.2.

**Proved here:** the entire TL3 algebraic dependency; exact \(\beta_1\), \(\pi\), \(\kappa\); the formal same-contact dictionary; and the missing bounded-\(b\) premise. **Not proved here:** TL6 bounded conductor, positivity of the inherited CA/SE source matrices, remote test joins, universal C1/C2 certificate admission, or the global SDP-lift claim. Those are distinct dependencies.
