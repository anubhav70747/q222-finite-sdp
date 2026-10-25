# A uniform multiplied SOS at native C2 contacts

The complete C2 state polynomial has a sum-of-squares identity after multiplication by its positive three-state principal determinant. The identity uses at most 183 squares from one fixed 676-function span. The multiplier depends on both the source and the testing angles; the identity therefore does not give an affine SDP lift for the unmultiplied source cone.

Let $D_u=1+u^2$, $D_v=1+v^2$, and let $d>0$, $n,m,C,A,D,\alpha\ge0$, with $\lambda=\beta=0$. The full polynomial is

$$
\begin{aligned}
F={}&4nD_uD_vX^2+4m(vW-X+uvY-uZ)^2
+4CD_v(uW-Y)^2+4AD_uD_vY^2\\
&+(4D+2d)(uW+uvX-Y-vZ)^2
+(4\alpha+2d)D_u(vY-Z)^2\\
&+2dD_v(X+uZ)^2-2dD_uD_vZ^2.
\end{aligned}
$$

Assume it has an actual stationary zero at fixed source-side $u_c,v_c>0$ and a nonzero state: $P(u_c,v_c)\psi_c=0$ and $\psi_c^TP_u\psi_c=\psi_c^TP_v\psi_c=0$, where $F=\psi^TP\psi$. The native C2 formulas satisfy these equations, including $n=0$. All testing variables below are independent of this contact.

## The scalar determinant

The determinant formula SG5–SG6 gives

$$
\det P=16D_u^2D_v^2u^2v^2L(u^2,v^2),\qquad
L(U,V)=h+\sum_{i=1}^6q_iU^{(e_i)_1}V^{(e_i)_2},
$$

where

$$
(e_1,\ldots,e_6)=((-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0)).
$$

For completeness, the nonnegative coefficients are

$$
(q_1,\ldots,q_6)=4(\alpha mt_5t_8,mDt_5t_6,
\alpha At_7t_3,CDt_5t_2,nAt_3t_4,CAt_1t_4),
$$

with

$$
\begin{array}{ll}
t_1=d+2(n+D),&t_2=d+2(\alpha+m+A),\\
t_3=d+2(C+D),&t_4=d+2(\alpha+m),\\
t_5=d+2n,&t_6=d+2(\alpha+C+A),\\
t_7=d+2(n+m),&t_8=d+2(C+A+D).
\end{array}
$$

Set $U_c=u_c^2,V_c=v_c^2$, $s_i=q_iU_c^{(e_i)_1}V_c^{(e_i)_2}$, and $x=(u^2/U_c,v^2/V_c)$. The contact equations imply

$$
h=-\sum_i s_i,\qquad
s_1-s_6=s_5-s_2=s_4-s_3=\varepsilon.
$$

Let $\tau_+=\max(\varepsilon,0)$, $\tau_-=\max(-\varepsilon,0)$, and

$$
\rho_1=s_1-\tau_+=s_6-\tau_-,\quad
\rho_2=s_5-\tau_+=s_2-\tau_-,\quad
\rho_3=s_4-\tau_+=s_3-\tau_-.
$$

These quantities are nonnegative. With $T_+=(e_1,e_4,e_5)$ and $T_-=(e_2,e_3,e_6)$,

$$
L=\sum_{i=1}^3\rho_i(x^{e_i}+x^{-e_i}-2)
+\tau_+\left(\sum_{e\in T_+}x^e-3\right)
+\tau_-\left(\sum_{e\in T_-}x^e-3\right).
$$

Each opposite pair gives the square

$$
u^2v^2(x^e+x^{-e}-2)=\big[uv(x^{e/2}-x^{-e/2})\big]^2.
$$

For either triangle $(f_1,f_2,f_3)$, whose sum is zero, the cubic identity gives

$$
u^2v^2\left(\sum_i x^{f_i}-3\right)
=\frac12\sum_{h=1}^3\sum_{i<j}
\big[uvx^{f_h/6}(x^{f_i/3}-x^{f_j/3})\big]^2.
$$

Indeed this is $a^3+b^3+c^3-3abc=(a+b+c)\sum_{i<j}(a_i-a_j)^2/2$, with $a_i=x^{f_i/3}$ and $abc=1$. Only one triangle has nonzero weight. Hence

$$
u^2v^2L=\sum_{\ell=1}^{12}\ell_\ell^2,
$$

where zero terms may be omitted. Every root is a linear combination of $u^{i/3}v^{j/3}$, $0\le i,j\le6$, with source-only coefficients. Its exponents are nonnegative: a triangle monomial has angular exponent $1+(f_h+2f_i)/3\in[0,2]^2$.

## Completing all four state coordinates

Write $\xi=(W,X,Y)^T$ and

$$
P=\begin{pmatrix}H&k\\ k^T&\gamma\end{pmatrix},\qquad
k=\begin{pmatrix}
-uv(4m+4D+2d)\\
u(4m+2d-4Dv^2)\\
v[4(D-\alpha)-(4m+4\alpha+2d)u^2]
\end{pmatrix}.
$$

The 19 rows $\varrho_i$ of a Gram factor for $H$ are specified by the following table. Multiply each row in the right column by the factor on its left.

| Factor | Rows |
|---|---|
| $\sqrt{2d}$ | $(u,uv,-1),(0,0,v),(0,0,uv),(0,1,0),(0,v,0)$ |
| $2\sqrt n$ | $(0,1,0),(0,u,0),(0,v,0),(0,uv,0)$ |
| $2\sqrt m$ | $(v,-1,uv)$ |
| $2\sqrt C$ | $(u,0,-1),(uv,0,-v)$ |
| $2\sqrt A$ | $(0,0,1),(0,0,u),(0,0,v),(0,0,uv)$ |
| $2\sqrt D$ | $(u,uv,-1)$ |
| $2\sqrt\alpha$ | $(0,0,v),(0,0,uv)$ |

Thus $H=\sum_i\varrho_i\varrho_i^T$, and the three-dimensional Cauchy–Binet identity is

$$
\operatorname{adj}H=\sum_{i<j}
(\varrho_i\times\varrho_j)(\varrho_i\times\varrho_j)^T.
$$

The core gives $H\succ0$ at all $u,v>0$, and

$$
\delta_H:=\det H\ge8d^3u^2v^2D_uD_v>0.
$$

The fraction-free Schur identity, with $\eta=H\xi+kZ$, now gives

$$
\boxed{\delta_H F
=\sum_{i<j}\big[(\varrho_i\times\varrho_j)\cdot\eta\big]^2
+\sum_{\ell=1}^{12}\big[4D_uD_v\ell_\ell Z\big]^2.}
$$

To verify it, $\eta^T\operatorname{adj}(H)\eta+\det(P)Z^2=\det(H)F$; then substitute the two SOS identities above. No inverse of a mixed source/testing expression occurs in these squares.

There are at most $\binom{19}{2}+12=183$ squares. Each cross product has bidegree at most $(2,2)$, while $H\xi+kZ$ has bidegree at most $(2,2)$. Each scalar root multiplied by $D_uD_v$ has degree at most four in either angle. All roots therefore belong to

$$
\mathcal Z=\{u^{i/3}v^{j/3}\psi_l:
0\le i,j\le12,\ 1\le l\le4\},\qquad |\mathcal Z|=676.
$$

This is one fixed span. Square coefficients, the contact constants and circuit weights are source-only. At $n=0$, four Gram rows vanish and the same identity holds literally. The construction concerns independent $u,v>0$; commuting and projective endpoints of the original slack are handled directly by positivity of its core and event operators, rather than by dividing by a vanishing multiplier.

## Two limitations verified exactly

The strict native source

$$
(d,n,m,C,A,D,\alpha)
=(1,7/66,5/33,5/22,10/33,16/33,9/22)
$$

has parameters $(a,b,c,t,r)=(1,5/4,1,3/4,1/4)$, with $(h,g,j)=(1/2,1,1/2)$. Define $R_m=vW-X+uvY-uZ$, $R_D=uW+uvX-Y-vZ$, and its contact normals

$$
H_n=X+ruZ,\quad H_m=R_m-guZ,\quad
H_A=Y-(h/a)uZ,\quad H_D=R_D+jvZ.
$$

At the independent test $u=2,v=1,\psi=(0,0,0,1)$,

$$
F-4nD_uD_vH_n^2-4mH_m^2-4AD_uD_vH_A^2-4DH_D^2
=-31/11.
$$

Thus subtracting these four entire normal squares need not leave a positive remainder. This does not rule out a coupled Gram for the full parent.

For this source, $\delta_H(2,1)/\delta_H(1,1)=22125/1408$. For the strict native seed

$$
(1,59/320,73/230,23/100,73/200,73/230,23/100),
$$

the ratio is $460405/31648$. They differ, so this particular multiplier cannot factor as one source scalar times one fixed testing polynomial across native C2. Since $\delta_H$ is cubic in source coefficients, the target $\delta_HF$ is quartic in them. Its Gram identity is not an affine coefficient identity for $F$. Neither this ratio test nor the multiplied SOS settles existence of another fixed multiplier or removal of the present one.

Run `python checks/c2_multiplied_certificate.py` from this directory with SymPy to reproduce [the exact checks](checks/c2_multiplied_certificate.json). The checker reconstructs the full parent, Gram factor, adjugate and Schur identities, both circuit signs and exponent bounds, and both rational examples. It uses no numerical optimization.
