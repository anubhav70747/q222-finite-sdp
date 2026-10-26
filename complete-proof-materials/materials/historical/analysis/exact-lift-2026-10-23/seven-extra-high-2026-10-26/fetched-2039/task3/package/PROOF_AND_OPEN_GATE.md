# Direct K364 C2 strict admission: all-domain fractional opening and the unclosed coupled transfer

**Scope and verdict (27 October 2026).** This is a verified **positive partial Gram construction on every genuinely physical strict native C2 source**. It does **not** prove universal H2, any PSD completion of the residual, or that a feasible K364 Gram must dominate this chosen carrier. It does not reuse the frozen tau=2 syzygy or the n=0 wall as an input. The qualitative SDP-existence theorem, verified open-end theorem, and global conditional 364 SDP are left unchanged.

Let `theta=(a,b,c,t,r)` satisfy the *actual* strict native C2 conditions: `a,b,c,t>0`, `0<r<1`, all six native weights `n,m,C,A,D,alpha>0`, and positive gaps

\[
h=b-a^2t>0,\quad g=(1+c^2)b-(1-r+a^2c^2t)>0,\quad j=1+c^2r-t>0,\quad k=1-c^2h>0,\quad B=b-a^2c^2r>0,
\]

with the native source numerator/denominator conditions `Delta>0,N_*>0,N_n>0`. Its coefficient order is `(d,lambda,n,m,C,A,D,alpha,beta)=(1,0,n,m,C,A,D,alpha,0)`. The genuine full parent, independent in testing variables `u,v` and state `(W,X,Y,Z)`, is

\[
\begin{split}
F_\theta={}&4n(1+u^2)(1+v^2)X^2+4m(vW-X+uvY-uZ)^2\\
&+4C(1+v^2)(uW-Y)^2+4A(1+u^2)(1+v^2)Y^2\\
&+(4D+2)(uW+uvX-Y-vZ)^2+(4\alpha+2)(1+u^2)(vY-Z)^2\\
&+2(1+v^2)(X+uZ)^2-2(1+u^2)(1+v^2)Z^2.
\end{split}\tag{1}
\]

The stationary contact is
\[
(u_c,v_c;W_c,X_c,Y_c,Z_c)=\left(c,ac;\frac ba,-cr,\frac{ch}{a},1\right).
\tag{2}
\]

## Theorem A. Exact manifestly positive extreme-u packets

Writing `Dv=1+v²`, the *entire unrestricted four-state* parent is
\[
 F_\theta=F_0(v;\psi)+uF_1(v;\psi)+u^2F_2(v;\psi),\tag{3}
\]
where
\[
\begin{aligned}
F_0={}&(4n+2)D_vX^2+4m(vW-X)^2+(4C+4A+2)D_vY^2\\
&+4D(Y+vZ)^2+4\alpha(vY-Z)^2,\\
F_1={}&8m(vW-X)(vY-Z)-8CD_vWY\\
&-2(4D+2)(W+vX)(Y+vZ)+4D_vXZ,\\
F_2={}&4nD_vX^2+4CD_vW^2+4AD_vY^2\\
&+(4m+4\alpha+2)(vY-Z)^2+(4D+2)(W+vX)^2.
\end{aligned}\tag{4}
\]

All coefficients of `F0` and `F2` in (4) are positive throughout `Omega`; each is an explicit sum of squares of polynomial state-linear roots of `v`-degree <=1. After multiplication by `Dv` these packets have fixed-grid roots of `v`-degree <=2. The mixed term `uF1` must be treated jointly and is *not* asserted SOS. Identity (3) is coefficientwise, including every state normal.

## Theorem B. Two-square, source-only, moving-contact fractional opening

For every pure-source gauge `0<=q<=1`, define
\[
\begin{split}
L_q&=(u-qv/a)X+(1-q)(cr/b)vW,\\
H_q&=\sqrt{2q/a}\,\sqrt{uv}\bigl(X+(acr/b)W\bigr).
\end{split}\tag{5}
\]

Both roots vanish **exactly** at (2): `u_c-v_c/a=0` and `X_c+(acr/b)W_c=0`. The source-only radicals exist in its real closed coefficient field. Their angular monomials are `uX,vX,vW,sqrt(uv)X,sqrt(uv)W`; multiplying either root by `v` still lies in the *same fixed* 364-feature dictionary
\[
Z_{6,12}=\{u^{i/6}v^{j/6}\psi_\ell:0\le i\le6,0\le j\le12,\ell=W,X,Y,Z\}.\tag{6}
\]

Consequently the *four-root* construction
\[
\mathcal C_{n,q}=4n\big(L_q^2+H_q^2+(vL_q)^2+(vH_q)^2\big)
=4nD_v(L_q^2+H_q^2)\tag{7}
\]
is a genuine globally sound PSD Gram in K364's polynomial feature space for every native strict `theta`. It uses no testing-dependent denominator or coefficient-field/test-field inverse.

The complete coefficient identity is
\[
\boxed{
L_q^2+H_q^2 =
\left(u^2+q^2\frac{v^2}{a^2}\right)X^2
+2\frac{cr}{b}(1+q)uvWX
+(1-q)^2\frac{c^2r^2}{b^2}v^2W^2
+2qa\frac{c^2r^2}{b^2}uvW^2
-2q(1-q)\frac{cr}{ab}v^2XW.}
\tag{8}
\]

In particular, **both** the mandatory coefficient `[u²X²]F=4n` and the zero coefficient `[uvX²]F=0` are matched *exactly* by (7). The first coefficient is a unique diagonal within K364: `uX` is the only feature pair yielding `u²X²` with the upper u exponent 1, so any PSD Gram uses diagonal `4n` there. Equation (7) gives a contact-compatible realization, not the invalid independent n-event packet.

The fixed source-algebraic choice
\[
 q_*(a)=\frac a{1+a}\in(0,1)
\tag{9}
\]
retains the rigorously positive unused raw diagonal budget
\[
[v^2X^2]\big(F_\theta-4n(L_{q_*}^2+H_{q_*}^2)\big)
=2+4n\frac{a(a+2)}{(1+a)^2}>0.
\tag{10}
\]
No analogous positivity is asserted for every other residual coefficient, much less for the residual *Gram*.

## The single remaining implication for H2

For `q=q_*(a)`, the residual `R_theta=Dv F_theta-C_{n,q_*}` has compulsory couplings
\[
[u v W X]R_\theta=[u v^3 W X]R_\theta=-8n\frac{cr}{b}(1+q_*)<0,\tag{11}
\]
while these coefficients are identically zero in `Dv F_theta`. Matching (11), all other coefficients and the moving-contact kernel by a **simultaneously PSD** Gram is not accomplished by (7). One *sufficient but not necessary* route is proving `R_theta` has a PSD Gram in the remaining fixed features; an arbitrary admitting Gram need not dominate the selected carrier. The true all-source H2 remains
\[
\boxed{\forall\theta\in\Omega\quad\exists P_\theta\succeq0:\quad
\operatorname{coeff}(Z_{6,12}^TP_\theta Z_{6,12})
=\operatorname{coeff}(D_vF_\theta).}\tag{12}
\]

This proof supplies **no such P_theta**, nor an exact separating functional for K364. It does not assert that adding just one uX row to the 18-feature n=0 wall certificate is sufficient; its new fractional contact-zero root and the remaining coupled transfers are essential distinctions.

## Reproduction

`verify_fractional_opening.py` independently expands (1)--(4), proves the exact source-contact cancellation and (8) with formal symbols, verifies the two compulsory coefficient keys and (10)--(11), reconstructs all six strict event weights at the physically genuine rational seed `(a,b,c,t,r)=(1,1,1,3/4,1/4)`, and verifies every root is a fixed-grid feature. It uses only SymPy, with no solver, sweep, external Git data, prior pickle, or downloaded package. `EXACT_DATA.json` records the constructed roots and the nonadmission status. The original reviewed source material is not modified.
