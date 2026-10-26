# Direct strict-C2 admission in K364: exact normal syzygies and the remaining scalar SOS gate

**Scope:** genuine six-event native C2 source; full four-state parent; source coefficients independent of testing variables. **Status:** a new exact division-free, fixed-feature full-parent normal completion and a positive universal highest-Bob-degree pivot are proved. **Universal K364 PSD admission remains open.** The independently checked strict C2 open-end theorem, the older scalar determinant circuits and the qualitative finite-SDP existence proof are not altered.

## 1. Native hypotheses and fixed library

Take the complete native strict C2 source `(a,b,c,t,r)` in
\[
\Omega=\{a,b,c,t>0,\ 0<r<1,\ h=b-a^2t>0,\ g=(1+c^2)b-(1-r+a^2c^2t)>0,\ j=1+c^2r-t>0,\ 1-c^2h>0,\ b-a^2c^2r>0,\ N_*>0,\ N_n>0\}.
\]
Here
\[
\begin{split}
\Delta={}&1+(1+c^2)b+c^2(1+c^2)b^2+2c^2(1+c^2)br+c^4r^2+2c^2r\\
&+a^2c^2[(1+c^2)b+1+c^2r(2-r)],\\
N_*={}&(1+c^2)(b+r)[b+a^2(1-t)]-a^2(1+a^2c^2)(1-r)t,\\
N_n={}&(b+r)(1+c^2r)(1+c^2r+a^2c^2t)-r\Delta.
\end{split}
\]
The six weights are precisely the positive rational native values
\[
\begin{split}
n&=N_n/(2r\Delta),\\
m&=(1+a^2c^2)(1-r)(1+c^2r)(1+c^2r+a^2c^2t)/(2g\Delta),\\
C&=(1+c^2r)N_*/(2a^2t\Delta),\\
A&=a^2(1+c^2r)(b+r)(1+c^2r+a^2c^2t)/(2h\Delta),\\
D&=(1+c^2)(b+r)(b-a^2c^2r)(1+c^2r+a^2c^2t)/(2a^2j\Delta),\\
\alpha&=c^4(b+r)N_*/(2(1-c^2h)\Delta).
\end{split}
\]
Set \(d=1,\lambda=\beta=0\), \(D_u=1+u^2\), \(D_v=1+v^2\). The complete cleared parent is
\[
\begin{split}
F={}&4nD_uD_vX^2+4m(vW-X+uvY-uZ)^2+4CD_v(uW-Y)^2+4AD_uD_vY^2\\
&+(4D+2)(uW+uvX-Y-vZ)^2+(4\alpha+2)D_u(vY-Z)^2\\
&+2D_v(X+uZ)^2-2D_uD_vZ^2.
\end{split}\tag{1}
\]
The fixed candidate cone is \(K_{364}=\{s\ge0: D_vF_s=Z_{(6,12)}^TPZ_{(6,12)},\ P\succeq0\}\), with \(Z_{(6,12)}=(u^{i/6}v^{j/6}\psi_l)_{0\le i\le6,0\le j\le12,l=W,X,Y,Z}\). The native contact is
\[
(u_c,v_c;W_c,X_c,Y_c,Z_c)=(c,ac;b/a,-cr,c(b/a-at),1).
\tag{2}
\]
The exact full-state and angular contact equations \(F=\partial_{W,X,Y,Z,u,v}F=0\) at (2) are the verified native-source inputs.

## 2. A division-free full-state identity in exactly 32 fixed polynomial rows

Let \(\xi=(W,X,Y)^T\). Use the **16 normalized row shapes** \(r_i=(\mathbf r_i,q_i)\) in the following order; each actual full row is \(k_i r_i\), with its positive **pure-source** scaling \(k_i\):

| Index | Normalized shape `(W,X,Y,Z)` | `k_i` |
|---:|---|---|
| 0–3 | `(0,fg,0,0)`, `f=1,u`, `g=1,v` lexically | `2 sqrt(n)` |
| 4 | `(v,-1,uv,-u)` | `2 sqrt(m)` |
| 5–6 | `(ug,0,-g,0)`, `g=1,v` | `2 sqrt(C)` |
| 7–10 | `(0,0,fg,0)`, `f=1,u`, `g=1,v` lexically | `2 sqrt(A)` |
| 11 | `(u,uv,-1,-v)` | `sqrt(4D+2)` |
| 12–13 | `(0,0,fv,-f)`, `f=1,u` | `sqrt(4alpha+2)` |
| 14–15 | `(0,g,0,ug)`, `g=1,v` | `sqrt(2)` |

Write \(R_0\in\mathbb R[u,v]^{16\times3}\) for the first three columns of these normalized rows, \(q_0\) for the last column, and \(K_0=\operatorname{diag}(k_0,\ldots,k_{15})\). Set
\[
\widehat R=\begin{pmatrix}K_0R_0\\vK_0R_0\end{pmatrix},\quad
\widehat q=\begin{pmatrix}K_0q_0\\vK_0q_0\end{pmatrix},\quad
\mathsf K=\operatorname{diag}(K_0,K_0).
\]
By expanding precisely the *original* nine terms in (1),
\[
\boxed{D_vF=\|\widehat R\xi+\widehat qZ\|^2-2D_uD_v^2Z^2.}\tag{3}
\]
Every row of \((\widehat R,\widehat q)\) has angular bidegree at most \((1,2)\). It is therefore a root in the **fixed** 364-feature dictionary. No state normal or angle is specialized.

## 3. Twenty-nine exact source-independent polynomial syzygies

Number the original *normalized* rows \(e_i\), and the 16 Bob-shifted rows \(e_i'=v e_i\). All the following are identities **in the three-dimensional \((W,X,Y)\) part**, not identities of the fourth-state coefficient:

1. For the five rows `i=1,2,3,14,15` proportional to \(X\): \(e_i-f_i g_i e_0=0\), where \(f_i g_i=(R_0)_{i,X}\).
2. For the five rows `i=8,9,10,12,13` proportional to \(Y\): \(e_i-f_i g_i e_7=0\), where \(f_i g_i=(R_0)_{i,Y}\).
3. \(e_6-ve_5=0\).
4. \(e_{11}-e_5-uv e_0=0\).
5. \(u e_4-ve_5+u e_0-ve_7-uv e_9=0\).
6. For each \(0\le i<16\), \(e_i'-ve_i=0\).

These yield **29 linearly independent kernel vectors** \(S(u,v)\in\mathbb R[u,v]^{32\times29}\), with every entry of bidegree at most \((1,1)\), such that
\[
\boxed{(R_0;vR_0)^TS(u,v)=0.}\tag{4}
\]
Their full sparse entries, row order and exact coefficients are machine-readable in `K364_SYZYGY_DATA.json`. Multiplying each vector by \(\mathsf K^{-1}\) produces actual source-scaled syzygies \(t\) satisfying \(\widehat R^Tt=0\) with the *same* testing degrees.

**Contact-kernel interpolation without any testing inverse.** At arbitrary \(u_c,v_c>0\), rows \(e_0,e_7,e_5\) form an invertible \(3\times3\) matrix of determinant \(u_c\). Consequently the 29 syzygies evaluated at the contact form a basis of \(\ker(R_0(u_c,v_c);v_cR_0(u_c,v_c))^T\). There is an explicit triangular interpolation procedure: for an arbitrary normalized contact-kernel vector \(b=(b_0,\ldots,b_{31})\), put \(\eta_i=b_i+v_cb_{16+i}\), set the 16 shift coefficients to \(b_{16+i}\), set every simple 13-column coefficient to \(\eta_i\) for its distinguished row except \(m\), set \(\lambda_m=\eta_4/u_c\), and adjust the \(Y_u\) coefficient to \(\lambda_{Y_u}=\eta_9+v_c\eta_4\). The three unused pivot-row equations then hold automatically from the kernel equation. **All divisions are by the positive pure-source \(u_c\)**; none is by \(u\), \(v\), or a mixed source/testing sum.

## 4. An explicit Bob-degree-two transfer and a positive all-domain high pivot

The following **additional** normalized syzygy has bidegree at most \((1,2)\):
\[
\boxed{T_{\rm high}=v^2e_{11}'-uv e_4'-uv^2e_{15}'-uv e_0'+v e_8'+uv e_{10}'.}\tag{5}
\]
Indeed its \(W,X,Y\) components cancel as
\[
 v^2(uv,uv^2,-v)-uv(v^2,-v,uv^2)-uv^2(0,v^2,0)
 -uv(0,v,0)+v(0,0,v^2)+uv(0,0,uv^2)=0.
\]
Unlike an unrestricted Koszul syzygy, this uses no \(u^2\) or \(v^3\) coefficient and **stays inside the fixed K364 degrees**.

For any source contact, put \(e_c=\widehat R_c\xi_c+\widehat q_c\), \(\xi_c=(W_c,X_c,Y_c)\). Native stationarity gives \(\widehat R_c^Te_c=0\). Let \(\Lambda(b)\) be the triangular interpolation operator of Section 3. Define source-only normalized contact vectors \(b_c=-\mathsf K e_c\), \(b_h=T_{\rm high}(u_c,v_c)\) and the **explicit polynomial syzygy**
\[
\boxed{t_\theta(u,v)=\mathsf K^{-1}\left(S(u,v)\Lambda(b_c)
  +2\left[T_{\rm high}(u,v)-S(u,v)\Lambda(b_h)\right]\right).}\tag{6}
\]
Every parameter in (6) depends on the source only; its entries have bidegree at most \((1,2)\), and
\[
\widehat R^Tt_\theta=0,\qquad t_\theta(u_c,v_c)=-e_c.\tag{7}
\]
In particular the rows of \(\widehat R\xi+(\widehat q+t_\theta)Z\) **all vanish at the actual four-state contact**.

For **every** polynomial syzygy \(t\) with \(\widehat R^Tt=0\), we have the exact full-parent normal-completion identity
\[
\boxed{D_vF=\sum_{i=1}^{32}\left[\widehat R_i\cdot\xi+(\widehat q_i+t_i)Z\right]^2
  +\mathscr S_t(u,v)Z^2,\quad
\mathscr S_t=-2D_uD_v^2-2\widehat q^Tt-t^Tt.}\tag{8}
\]
This is a literal polynomial identity (not merely at contact or in residue); its first term already belongs to the fixed K364 cone for all \(\theta\in\Omega\). With \(t=t_\theta\), the scalar remainder \(\mathscr S_t\) vanishes with both angular first derivatives at the native contact, because every retained square vanishes there and the entire original parent is stationary.

The degree-four Bob pivot is **universally positive, with no restriction on the sign of \(h-a^2r\)**. The \(v^2\) coefficient of the vector \(\widehat q\) has precisely two entries \(-\sqrt{4D+2}\) in row 27 and \(u\sqrt2\) in row 31. The contact interpolation \(S\Lambda\) in (6) has degree at most one in \(v\); the high syzygy contributes \(2/\sqrt{4D+2}\) to row 27 and \(-2u/\sqrt2\) to row 31 at degree two. Therefore **exactly**
\[
\begin{split}
[v^4]\mathscr S_{t_\theta}
&=-2(1+u^2)+4(1+u^2)
  -4\left[\frac1{4D+2}+\frac{u^2}2\right]\\
&=\boxed{\frac{4D}{2D+1}>0.}
\end{split}\tag{9}
\]
This sign is valid for every strict native source, including arbitrarily small source gaps (pointwise, not uniform at their zero limits). No numerical comparison of source and testing rates is involved.

## 5. The single unresolved PSD implication for this constructive route

Equation (8) certifies the complete four-state parent in K364 **if** the scalar remainder from (6), or from an admissible source-only syzygy deformation of it, has a finite SOS using scalar roots in the fixed grid
\[
\big\{u^{i/6}v^{j/6}Z:0\le i\le6,\ 0\le j\le12\big\}.
\tag{10}
\]
This is a specific *sufficient normal-assembly ansatz*, not a necessary characterization of all possible K364 Grams. The positive leading coefficient (9), exact contact tangency, and full state-normal matching **do not establish** that scalar SOS, nor positivity of the entire scalar polynomial. Lower Bob degrees, angular contact circuits, and their Schur couplings are the first unproved implication. One must find source-algebraic coefficients and an exact positive Gram for (10), possibly after varying the admissible syzygy with zero contact evaluation. The original K364 cone might admit a different Gram even if this particular scalar-completion ansatz fails.

For comparison, the original 141-square identity involves the test-dependent multiplier \(\det B\). The present construction **never introduces that multiplier in (8)**. Independently, exact expansion of the three-state principal block gives
\[
\det B=8(1+u^2)(1+v^2)H(u^2,v^2),\qquad
H=\sum_{(i,j)\in\{0,1,2\}^2\setminus\{(0,0)\}}h_{ij}\,u^{2i}v^{2j},
\tag{11}
\]
with all eight \(h_{ij}\) strictly positive for every \(n,m,C,A,D,\alpha>0\); the checker writes all eight polynomials to `K364_SYZYGY_DATA.json`. Notice \(H(0,0)=0\), while \(1+v^2=1\) there. Thus the old multiplied-parent certificate cannot be divided by \(8(1+u^2)H\) inside this fixed Gram library; (8) is the direct denominator-free alternative, whose scalar PSD step is genuinely still open.

**No theorem of all-strict-C2 K364 admission, nor of the global explicit SDP lift, is claimed.** The orientation-restricted C1 problem is being solved independently. The verified open-end theorem, source-density conditional, reflection, and 128-relabeling assembly retain their previous stated scopes.

## 6. Reproduction

Run `python verify_syzygy.py`. It uses SymPy only and checks, without optimization or samples: original four-state parent row expansion; all 29 syzygies; the admissible extra v² syzygy; the full rank-29 contact interpolation with an explicit pure-source inverse; exact \((1+v^2)F\) identity; contact compatibility of the v² correction; the strict source-uniform pivot (9); and the 8 positive coefficient formulas (11). The output is `K364_SYZYGY_DATA.json`, including every sparse syzygy and positive determinant coefficient. None of this is an all-source PSD check of (10).
