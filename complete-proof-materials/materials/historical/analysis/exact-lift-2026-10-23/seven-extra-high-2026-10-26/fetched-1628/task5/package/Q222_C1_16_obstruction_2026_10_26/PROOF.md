# Strict C1: exact negative resolution of OPEN-C1-16 and the fixed-library boundary

**26 October 2026.** This proof is about the **entire unrestricted four-state parent and independently variable positive measurement half-angle coordinates**. It works with the genuine native C1 dictionary SD1, not a formal source or a scalar contact restriction. There is no SDP solve, search over witness points, numerical sweep, or publication operation. It retains the already-proved 42-row, rank-36 exact contact coefficient map as an input and does not reprove that map.

## Theorem (OPEN-C1-16 is false on A1)

Let A1 be the **exact strict domain printed in the recovered C1 16-feature proof**: the native positive parameters `k,f,g,a,t` and the coquotient/denominator and strict source-sign inequalities, with **no additional inequality `k<=1` appended**. Let `Z16` contain the fourteen integer rows

\[
(uW,vW,X,uvX,Y,uvY,uZ,vZ,uvW,vX,uY,vY,Z,uvZ)
\]

and the two rows \(\sqrt{uv}W,\sqrt{uv}Z\). On this domain, the assertion that every genuine strict native C1 source has a PSD Gram in `Z16` is **false**. Specifically

\[
\boxed{(k,f,g,a,t)=(4,1/1000,1/8,1,1)}\tag{1}
\]

satisfies all strict native source inequalities, has a genuine globally quantum-valid stationary contact, but **has no PSD Gram in these sixteen features, regardless of the eight gauges**. In fact an open subset of strict A1 sources, including a sufficiently small positive-`f` part of the family `(k,g,a,t)=(4,1/8,1,1)`, is excluded.

This theorem neither excludes a PSD Gram using additional sixth-root features nor refutes the complete qualitative finite-SDP existence theorem. The physical Alice reflection may orient to `0<k<=1`, but that reflection is **not proved to preserve this exact 16-feature subspace**. Thus the stated obstruction is for A1 as recovered, which permits `k>1`; it does not purport to refute an extra, separately imposed `k<=1` representative-only claim.

## 1. Exact strict native source and full-state validity

Retain the native SD1 notation

\[
\kappa=k^2,\ \eta=f^2,\ \zeta=g^2,\quad
L_c=1-a\eta,\ H_c=1+(1-a)\eta,\ Q_c=1-\kappa\zeta+a\kappa\eta\zeta,
\]
\[
R=1+\kappa L_c,\quad T=1+\kappa-\kappa\zeta t,\quad M=R-t,\quad J=T-a.
\]

The seven source coquotients `c2, cb, c3, c5, c4, calpha, cd` and `P4, Palpha, Pd` are those printed in `c1_native.py` and the recovered original C1 TeX; this is an exact, fully algebraic source dictionary. The complete physical strict domain A1 requires

\[
\begin{gathered}
 k,f,g,a,t,L_c,H_c,Q_c,M,J>0,\quad 1-t\zeta>0,\quad
 1-\kappa\eta+t\kappa\eta\zeta>0,\\
 P_4,P_\alpha,P_d,c_2,c_b,c_3,c_5,c_4,c_\alpha,c_d>0,\\
 D_{\rm e}=2ac_5/(Jc_3)-d_{\rm B}/2>0,\qquad
 \beta=2ae^2c_b/(tc_3)-d_{\rm B}/2>0,
\end{gathered}\tag{2}
\]

with `e=kf/g`, `d_B=2a eta cd/c3`; the other three physical weights `m,A,alpha` are strictly positive by these formulas. These signs make every native denominator legal and give the actual stationary contact

\[
(u_c,v_c;W_c,X_c,Y_c,Z_c)=
(k,kfg;g/f,-k(1-tg^2),kgL_c/f,1).\tag{3}
\]

At (1), the contact is

\[
\boxed{(u_c,v_c;\psi_c)=(4,1/2000;125,-63/16,999999/2000,1).}\tag{4}
\]

For exact reproducibility, the entire original source is

\[
\begin{split}
(d,\lambda,n,m,C,A,D,\alpha,\beta)=(&923999207/13599785800000,0,0,\\
&8000002/63491000001,2,\\
&53599207/27199544400428400000,\\
&11449570655261/90188053200000,\\
&228841893/21473353157782,\\
&2883242147/1431556400000).
\end{split}\tag{5}
\]

All nonomitted weights are **strictly positive**. The exact positive gaps are, among others,

\[
 L_c=999999/1000000,\quad Q_c=3000001/4000000,\quad
 M=999999/62500,\quad J=63/4,\quad 1-tg^2=63/64.
\tag{6}
\]

The exact checker evaluates the whole SG2 parent with the seven weights (5), verifies `F=0` at (4), and verifies **all four unrestricted state gradients and both angular gradients vanish identically**. The previously proved SG7 stationary-to-validity theorem therefore yields genuine global quantum validity. In particular, the signed-core and all independent test/state directions, not just the contact graph, are covered.

## 2. Source-wide necessary **lower** bound on the sole fractional Gram gauge

Make the invertible source-only normalization

\[
 u=u_cU,\quad v=v_cV,\quad
 (W,X,Y,Z)=\operatorname{diag}(W_c,X_c,Y_c,Z_c)(\phi_W,\phi_X,\phi_Y,\phi_Z).
\]

At the positive contact `U=V=1`, `phi=(1,1,1,1)`. The established contact-kernel coefficient map shows that every PSD `Z16` Gram would have the form

\[
\widehat F=\sum_j p_j(U,V,\phi)^2+\gamma UV(\phi_W-\phi_Z)^2,
\qquad\gamma\ge0,\tag{7}
\]

where the `p_j` are **ordinary real polynomials** in the fourteen integer features. We do not assume that the integer polynomial alone inherits the Bell sign invariance. Instead use the four exact sign-related contact zeros of the **original** parent:

\[
 (U,V;\phi_W,\phi_X,\phi_Y,\phi_Z)
 =(\sigma,\tau;1,\tau,\sigma,\sigma\tau),
 \qquad\sigma,\tau\in\{\pm1\}.\tag{8}
\]

These zeros are exact by the physical sign congruences of SG2. At `sigma*tau=1`, the fractional term in (7) is zero, so the entire integer-root vector `p=(p_j)` is zero. At `sigma*tau=-1`, the fractional term is `-4 gamma` and `F=0`, whence

\[
\|p_{+-}\|^2=\|p_{-+}\|^2=4\gamma.\tag{9}
\]

The average of every integer feature over (8) is zero **except** `V phi_X`, `U phi_Y`, and `UV phi_Z`. Let `a,b,e` denote their coefficient vectors in the polynomial-square representation. Averaging the four root vectors gives

\[
 a+b+e=\tfrac14(p_{+-}+p_{-+}),\qquad
 \|a+b+e\|\le\sqrt\gamma.\tag{10}
\]

Because the monomials `V^2 phi_X^2`, `U^2 phi_Y^2`, and `U^2V^2 phi_Z^2` each have a **unique integer feature square** and are absent from the fractional square (7), their squared coefficient norms are exactly

\[
\boxed{S_c=(4\beta+2d)v_c^2X_c^2,\quad
S_b=4A u_c^2Y_c^2,\quad
S_e=4\beta u_c^2v_c^2Z_c^2.}\tag{11}
\]

Thus, with `L=max(0,sqrt(S_b)-sqrt(S_c)-sqrt(S_e))`, the reverse triangle inequality proves the **source-wide necessary bound**

\[
\boxed{\gamma\ge L^2.}\tag{12}
\]

There is no square root of a mixed tensor element in this argument: (12) is a real inequality for the **source-only** Gram coefficient, derived from four ordinary real polynomial evaluations and Euclidean norms.

## 3. Source-wide forced **upper** bound from the strict Schur pivot

The previously verified 42-row, rank-36 coefficient map, with eight source-only gauge directions, is retained. It fixes the following reduced entries for **every** exact contact-matching Gram in the integer odd-parity `7x7` block (zero-based indices):

\[
\begin{aligned}
H_{00}&=2Bu_c^2W_c^2,&H_{01}&=-\gamma/2,&H_{02}&=0,\\
H_{04}&=-2Bu_cW_cY_c,&H_{14}&=H_{24}&=0,\\
H_{11}&=4mv_c^2W_c^2,&H_{12}&=-4mv_cW_cX_c,\\
H_{22}&=2(2\beta+d+2m)X_c^2,&H_{44}&=2(2A+B)Y_c^2,
\end{aligned}\tag{13}
\]

where `B=2C+2D+d`. The identity `H01=-gamma/2` is the exact cancellation of the `UV phi_W^2` coefficient between the two integer features and the fractional square. These entries and the displayed zeros are unaffected by the other seven gauges. The separate `H5[11,11]=4 alpha Z_c^2>0` is also fixed.

The forced block `P=H7[{1,2}]` satisfies

\[
\delta=\det P=8m(2\beta+d)v_c^2W_c^2X_c^2>0,\tag{14}
\]

while the forced `H7[{0,4}]` has determinant

\[
8ABu_c^2W_c^2Y_c^2>0.\tag{15}
\]

Schur-complement `P` in the *four-by-four* principal matrix on `{0,1,2,4}`. Since `H02=H14=H24=0`, its remaining `(0,4)` matrix is

\[
\begin{pmatrix}
H_{00}-\gamma^2H_{22}/(4\delta)&H_{04}\\
H_{04}&H_{44}
\end{pmatrix}\succeq0.\tag{16}
\]

Its determinant, including all source/contact factors, yields the necessary **source-wide upper bound**

\[
\boxed{\gamma^2\le U_\theta:=
\frac{64mA(2\beta+d)(2C+2D+d)u_c^2v_c^2W_c^4}
{(2\beta+d+2m)(2A+2C+2D+d)}.}\tag{17}
\]

This is forced even if the *full* Gram is singular. No positive-definiteness of the Gram, gauge continuity, or source-uniform cap is used.

## 4. Exact rational infeasibility certificate

At (1), all three squared lengths (11) are **rational**:

\[
\begin{aligned}
S_c&=\frac{3969}{124998031250},\\
S_b&=\frac{53599153400793}{1699973225000000000},\\
S_e&=\frac{2883242147}{89472275000000000}.
\end{aligned}\tag{18}
\]

Direct cross-multiplication in the integers proves

\[
S_b>\frac1{200^2},\qquad
S_c<\frac1{4000^2},\qquad
S_e<\frac1{4000^2}.\tag{19}
\]

The Fourier inequality (12) therefore forces

\[
\boxed{\gamma>\left(\tfrac1{200}-\tfrac2{4000}\right)^2
 =\frac{81}{4\,000\,000}.}\tag{20}
\]

On the other hand the exact rational value of (17) is

\[
U_\theta=
\frac{1097712033787939840000000}
 {75086759814055431645960567199993023}
 <\frac1{100000^2}.\tag{21}
\]

Consequently any PSD Gram would require

\[
\boxed{0\le\gamma<\frac1{100000}<\frac{81}{4\,000\,000}<\gamma,}\tag{22}
\]

an explicit **algebraic impossibility certificate**. Equations (19) and (21) are checked by exact rational ordering, not floating-point approximations. Moreover, the two PSD consequences can be combined into **one exact linear separator on the exposed contact-zero face**, as follows.


### A single rational separating functional on the exposed contact face

The Fourier argument admits a **linear** Cauchy variant. Write `S=a+b+e` in the notation of (10). Then `b=S-a-e` and the Euclidean inequality `||r+s+t||^2 <= 3(||r||^2+||s||^2+||t||^2)` imply, using `||S||^2<=gamma`,

\[
\boxed{\gamma\ge L_{\rm lin}:=\frac{S_b}{3}-S_c-S_e.}\tag{22a}
\]

At the exact source (1),

\[
\boxed{L_{\rm lin}=\frac{8878812233069}{849986612500000000}>\gamma_0:=\frac1{100000}.}\tag{22b}
\]

Choose the following **fixed rational Gram test vector**, supported entirely on the forced coordinates `(0,1,2,4)` of `H7`:

\[
\begin{split}
q_0&=e_0+a_0e_1+b_0e_2+c_0e_4,\\
a_0&=\frac{\gamma_0H_{22}}{2\delta}
=\frac{276058888004033}{102400025600000}>0,\\
b_0&=-\frac{H_{12}}{H_{22}}a_0
=-\frac{3999937}{1612800000},\\
c_0&=-\frac{H_{04}}{H_{44}}
=\frac{68000000000000}{67999932000063}.
\end{split}\tag{22c}
\]

For **any** contact-reduced exact `Z16` Gram, all entries touched by this vector are forced and

\[
q_0^TH_7(\gamma)q_0
=Q_0+a_0(\gamma_0-\gamma),\qquad
Q_0=q_0^TH_7(\gamma_0)q_0
=-\frac{64109639476176033245960567199993023}
{5570468085829570501939179678720000000000}<0.
\tag{22d}
\]

Thus **both** `q_0^TH7q_0>=0` and `gamma-L_lin>=0` are mandatory PSD consequences. Multiply the second by the positive rational `a0` and add to the first. The `gamma` coefficient cancels **exactly**, leaving the source-only coefficient functional

\[
\boxed{\begin{aligned}
\mathcal L(\widehat F)&:=q_0^TH_7(\gamma)q_0+a_0(\gamma-L_{\rm lin})\\
&=Q_0+a_0(\gamma_0-L_{\rm lin})\\
&=-\frac{75230062831940553240656946536904079000951}
 {5918622341193918658310378408640000000000000000}<0.
\end{aligned}}\tag{22e}
\]

Every term defining `\mathcal L` is **linear in the normalized parent coefficients**; the joint-half gauge disappears. For every `Z16`-represented parent **on the exposed contact-zero face** `F(1,1;1,1,1,1)=0`, its PSD Gram yields `\mathcal L>=0`. The displayed strict negative rational therefore constitutes a reproducible **linear separating certificate for that face**. This stronger certificate is not claimed to be positive on the entire unexposed ambient polynomial cone or on the full `Z6` cone. No numerical eigensolver or gauge selection is involved.

## 5. An entire analytic strict C1 branch, not an isolated witness

Keep `(k,g,a,t)=(4,1/8,1,1)` and let `f>0` tend to zero. The native rational coquotients give the following **exact leading coefficients**:

\[
\frac d{f^2}\longrightarrow\frac{1155}{17},\quad
\frac m{f^2}\longrightarrow126,\quad
\frac A{f^4}\longrightarrow\frac{67}{34},\quad
\frac\beta{f^2}\longrightarrow\frac{68477}{34},\quad
D\longrightarrow\frac8{63}.\tag{23}
\]

All other native gaps and denominator limits remain strictly positive; `m,A,alpha,beta,d>0` for sufficiently small positive `f`. In this family

\[
(u_c,v_c,W_c,X_c,Y_c,Z_c)=
\left(4,\frac f2,\frac1{8f},-\frac{63}{16},\frac{1-f^2}{2f},1\right).
\]

Substitution into the source-wide necessary estimates gives

\[
\boxed{
\lim_{f\downarrow0}\frac{(\sqrt{S_b}-\sqrt{S_c}-\sqrt{S_e})^2}{f^2}
 =\frac{536}{17},\qquad
\lim_{f\downarrow0}\frac{U_\theta}{f^4}=
 \frac{270144}{18479}.}\tag{24}
\]

The exact comparison

\[
\left(\frac{536}{17}\right)^2-\frac{270144}{18479}
 =\frac{307698304}{314143}>0\tag{25}
\]

shows `gamma_min>sqrt(U_theta)` for **every sufficiently small** positive `f`. By continuity and strict native signs, every such source has an open five-parameter neighborhood of strict sources also outside the 16-feature cone. This proof is analytical and makes no grid, sweep, or connectedness-to-PSD argument.

## 6. Complete forced Schur reduction, including **singular** certificates

For completeness, the promised actual elimination of the two universally strict inputs works at **every** strict A1 source, not just at (1). Keep the fully reconstructed affine blocks `H7(theta;s),H5(theta;s)` of the original eight-gauge constructor. Put

\[
P=H_7[\{1,2\},\{1,2\}]\succ0,
\quad a_4=H_{44}>0,
\quad a_{11}=H_{11,11}^{(5)}=4\alpha Z_c^2>0.
\]

For `I=(0,3,4,5,6)` define the **exact** five-by-five Schur matrix

\[
J_{ij}=H_{ij}-
 \frac{H_{22}H_{i1}H_{j1}-H_{12}(H_{i1}H_{j2}+H_{i2}H_{j1})+H_{11}H_{i2}H_{j2}}
 {\delta},\quad i,j\in I.\tag{26}
\]

Since `H14=H24=0`, eliminating the strictly positive node `4` commutes with (26). For `I'=(0,3,5,6)` put

\[
T_{ij}=J_{ij}-\frac{J_{i4}J_{4j}}{H_{44}},\qquad i,j\in I'.\tag{27}
\]

Then **exactly**

\[
H_7\succeq0\ \Longleftrightarrow\ T\succeq0.
\]

The first scalar pivot of `T` simplifies to

\[
\boxed{d_0=T_{00}=
\frac{4A(2C+2D+d)u_c^2W_c^2}{2A+2C+2D+d}
-\frac{(2\beta+d+2m)\gamma^2}
{16m(2\beta+d)v_c^2W_c^2}.}\tag{28}
\]

Thus `d0>=0` is precisely (17). For `R=(3,5,6)`, let `b=(T03,T05,T06)^T` and `B3=T[R,R]`. If `d0>0`, the **complete** condition is the one three-by-three inequality

\[
\boxed{B_3-bb^T/d_0\succeq0.}\tag{29}
\]

If `d0=0`, there is no division by a vanishing Schur pivot: the **necessary and sufficient** condition is

\[
\boxed{b=0,\quad B_3\succeq0.}\tag{30}
\]

For the other block, with `J5=(7,8,9,10)`, define

\[
S^{(5)}_{ij}=H_{ij}-\frac{H_{i,11}H_{11,j}}{4\alpha Z_c^2},\qquad i,j\in J5.\tag{31}
\]

Then `H5>=0` if and only if `S^(5)>=0`. If desired one can further eliminate `i=10` with scalar pivot

\[
 e_{10}=2v_c^2Y_c^2(2A+2C+2\alpha+d)-\frac{s_6^2}{4\alpha Z_c^2}.\tag{32}
\]

For `e10>0`, a final three-by-three Schur inequality on indices `(7,8,9)` is exact; for `e10=0`, its complete off-diagonal vector must vanish and the remaining three-by-three principal matrix must be PSD. Therefore no unnecessary strict-Gram feasibility assumption is present anywhere in the reduction.

The **only remaining admission test after the forced eliminations**, for sources not already excluded by (12)–(17), is the simultaneous PSD feasibility of (27), (31) and `gamma>=0` in the seven other real gauges. It is not supplied by the merely positive forced masses. Crucially the source (1) fails **before** this unresolved residual test, so OPEN-C1-16 itself is decisively settled negatively.

## 7. Necessary finite-library enlargement and the unproved sufficiency

All ordinary integer-angle polynomial squares representing a bidegree-`(2,2)` parent have testing degree at most `(1,1)` in each root: an uncancelled highest-degree coefficient of a sum of real squares cannot disappear. The C1 Newton support forces the constant `W` and `uX` rows to vanish. Thus the fourteen integer rows of `Z16` exhaust **all permitted integer-angle features** of this degree. The theorem therefore proves:

> **Necessity.** Any single fixed monomial sublibrary of the full sixth-root grid that admits **all** strict A1 C1 sources must activate at least one genuinely new **noninteger-exponent** row outside the existing shared half-power `sqrt(uv)W,Z` pair. Adding more integer-degree roots or changing only source-dependent Gram gauges cannot fix the exhibited failure.

A natural explicit finite enlargement inside the full original `Z6` grid is

\[
\boxed{Z_{18}=Z_{16}\cup\{\sqrt{uv}X,\sqrt{uv}Y\}\subset Z_6.}\tag{33}
\]

Its certificate cone is finite affine-SDP and globally sound for the **canonical Bell source polynomials**, because a positive Gram proves the original parent nonnegative for all positive test angles and all unrestricted states, and the original Bell sign congruences restore the other quadrants and homogeneous endpoints. Its ambient angular polynomial SOS cone is **strictly** larger: `uv(X-Y)^2` is the single added square `sqrt(uv)(X-Y)` and is not representable in `Z16`, since setting `W=Z=0` kills all old fractional features, while any SOS of the remaining integer polynomials is nonnegative also at `u=-1,v=1,X=1,Y=0`, where `uv(X-Y)^2=-1`.

This strictly larger **ambient** cone does **not** constitute a proof that the native witness (1), let alone the entire strict C1 family, belongs to the enlarged **physical-source** cone. The eight-gauge contradiction does not transfer unchanged to (33), because the fractional block now has additional contact-root directions. Neither sufficiency of `Z18` nor sufficiency of the full 196-feature `Z6` has been proved by this task. That is the next actual admission implication, not another audit of the qualitative theorem.

## 8. Exact reproducibility and provenance

- `verify_exact_no16.py` checks every A1 sign and original weight, all six physical stationarity equations, the exact four-sign character support, the rational Fourier lower bound, the factorized Schur upper bound, and their strict contradiction.
- `verify_forced_schur.py` independently reconstructs the needed forced universal coefficient entries from the already retained integer map and verifies the complete strict-block Schur expression in symbolic source/angle letters, without a solver.
- `verify_analytic_family.py` verifies all five native rational leading coefficients, the two asymptotic gamma limits, and their rational strict comparison.
- `original_baseline/` retains unchanged original SD1/SD2 source, contact rank, gauge-constructor and verification files. The source map and prior 42/36 result are inputs, not rediscovered verdicts.

The evaluated global Q222 fixed lift, the positive C1 `Z18` or `Z6` admission and any claim of an orientation-restricted `k<=1` counterexample are **not** established here. The independently proved qualitative existence and C2 strict open-end results are unaffected.
