# Native strict C2: exact sixth-root interior reduction and the remaining PSD inequality

**Status (26 October 2026).** The strict open-`c` end admission from the independently verified `C2_FiniteC_End_Theorem` is preserved and **not rederived** here. The complete strict five-parameter interior does **not** receive a positive semidefinite Gram theorem from the present argument. This is an exact full-parent coefficient, contact and dual-cone reduction with a sharp, explicit unresolved semialgebraic inequality. All assertions labelled proved are independent of SDP solvers, sample searches, or missing pickle files.

## 1. Actual native source, parameter domain, and full parent

Work over positive real parameters \(\theta=(a,b,c,t,r)\). Set \(z=c^2\), \(s=1+zr\) and
\[
\begin{split}
h&=b-a^2t,\quad g=(1+z)b-(1-r+a^2zt),\quad j=1+zr-t,\\
\Delta&=1+(1+z)b+z(1+z)b^2+2z(1+z)br+z^2r^2+2zr
  +a^2z[(1+z)b+1+zr(2-r)],\\
N_*&=(1+z)(b+r)[b+a^2(1-t)]-a^2(1+a^2z)(1-r)t,\\
N_n&=(b+r)s(s+a^2zt)-r\Delta.
\end{split}
\]
The precise **strict native domain** \(\Omega\) used here is
\[
 a,b,c,t>0,\quad 0<r<1,\quad h,g,j,1-zh,b-a^2zr,N_*,N_n>0.
\tag{1}
\]
The polynomial \(\Delta>0\) on this domain. Its six strictly positive C2 event weights, with \(d=1,\lambda=\beta=0\), are the original rational family
\[
\begin{split}
n&=\frac{N_n}{2r\Delta},\\
m&=\frac{(1+a^2z)(1-r)s(s+a^2zt)}{2g\Delta},\\
C&=\frac{sN_*}{2a^2t\Delta},\\
A&=\frac{a^2s(b+r)(s+a^2zt)}{2h\Delta},\\
D&=\frac{(1+z)(b+r)(b-a^2zr)(s+a^2zt)}{2a^2j\Delta},\\
\alpha&=\frac{z^2(b+r)N_*}{2(1-zh)\Delta}.
\end{split}\tag{2}
\]
The *entire* parent, retaining every real state coordinate and every mixed state coefficient, is
\[
\begin{split}
F_\theta(u,v;W,X,Y,Z)={}&4nD_uD_vX^2+4m(vW-X+uvY-uZ)^2\\
&+4CD_v(uW-Y)^2+4AD_uD_vY^2\\
&+(4D+2)(uW+uvX-Y-vZ)^2\\
&+(4\alpha+2)D_u(vY-Z)^2+2D_v(X+uZ)^2-2D_uD_vZ^2,
\end{split}\tag{3}
\]
with \(D_u=1+u^2,D_v=1+v^2\). The native contact is
\[
(u_c,v_c;\psi_c)=\left(c,ac;\frac ba,-cr,c(\frac ba-at),1\right).
\tag{4}
\]
Every coordinate in (4) is nonzero on (1). All six source-field partial derivatives \(\partial_{W,X,Y,Z,u,v}F_\theta\) vanish at (4). **Independent verification:** `verify_native_jets.py` constructs the six original polynomial derivatives, multiplies by the explicitly factored positive native denominator, and verifies that all resulting five-variable polynomial numerators are identically zero. This was executed successfully in the current run. No state normal was frozen.

## 2. Full fixed global library and exact constant coefficient map

Put \(u=cU,v=acV\), \(\psi=\operatorname{diag}(b/a,-cr,c(b/a-at),1)\phi\), so that the entire contact becomes \((U,V;\phi)=(1,1;\mathbf1)\). This is invertible over the *source* field; the testing angles \(U,V>0\) and four real state coordinates \(\phi\) remain independent.

Let \(z_{ij\ell}=U^{i/6}V^{j/6}\phi_\ell\), \(0\le i,j\le6\), \(\ell\in\{W,X,Y,Z\}\). For the **strict C2 support face**, the Newton polygons force (and, for PSD Grams, necessitate)
\[
 i+j\ge6\quad\text{for }W,\qquad i+j\le6\quad\text{for }Z.
\tag{5}
\]
There are 154 retained features, split by angular sign parity into block orders \((52,36,36,30)\). Within each block choose its first lexicographic feature \(z_{q,0}\) and form all differences \(d_{q,k}=z_{q,k}-z_{q,0}\). They vanish identically at the contact for any testing/state rate. Their orders are \((51,35,35,29)\) and they contain \(3021\) independent symmetric Gram entries.

Let \(\mathscr S\) be the 412 ordered monomial keys \((I,J,\ell,m)\), \(0\le I,J\le12,\ell\le m\), arising as products of features *within a parity block*. Define the constant integer map
\[
\mathcal A:\bigoplus_{q=1}^{4}\mathbb S^{n_q-1}\to\mathbb R^{412},\qquad
\mathcal A(H)=\operatorname{coeff}_{\mathscr S}\sum_{q=1}^4 d_q^TH_qd_q.
\tag{6}
\]
For each off-diagonal Gram variable \(H_{ij},i<j\), its coefficient column is twice the polynomial product \(d_id_j\); for a diagonal variable it is \(d_i^2\). This fixes every factor of two. `c2_fixed_map.py` regenerates the complete integer columns, their ordering, and all 412 rows. `C2_154_CONTACT_MAP.json` stores every nonzero sparse coefficient; it is not a numerical approximation.

Let the six rows of \(\mathcal N\in\mathbb Z^{6\times412}\) be the four state derivatives \(\delta_{\ell,k}+\delta_{m,k}\) and two scaled angular derivatives \(I,J\), evaluated at the standard contact. Every column of (6) is killed by \(\mathcal N\). The six rows are independent. Exact modular Gaussian elimination at prime 65521 produces a **406-column, 406-row nonzero minor**, whose original row and column indices and nonzero modular pivot residues are saved in `C2_RANK406_MODULAR_PIVOTS.json` and reproducible by `verify_rank_mod65521.py`.

Therefore
\[
\boxed{\operatorname{rank}_{\mathbb Q}\mathcal A=406,
\qquad \operatorname{im}\mathcal A=\ker\mathcal N.}\tag{7}
\]
The first inequality follows from the nonzero 406-minor mod 65521; the second from dimension counting and the six exact annihilators. The contact equations verified in §1 imply that the complete source-normalized coefficient vector \(f(\theta)\) belongs to \(\ker\mathcal N\) **for every strict native source**, and hence admits at least one **source-rational coefficient-matching Gram** \(\mathcal A(H)=f(\theta)\). One deterministic rational construction takes the saved nonsingular minor \(B=\mathcal A_{I,J}\), sets all unselected Gram variables to zero, and solves \(H_J=B^{-1}f_I(\theta)\). It is exact but may be indefinite; this is not a PSD admission theorem.

## 3. A single positive source denominator and 24-term polynomial target

Define
\[
 \mathscr D(\theta)=2a^4t\,r g h j(1-c^2h)\Delta>0.\tag{8}
\]
After the contact normalization, \(f(\theta)=P(\theta)/\mathscr D(\theta)\), where **each of the 412 components of \(P\)** is an integer polynomial in \(a,b,c,t,r\). Exactly **24 components are nonzero**. `generate_cleared_targets.py` derives their explicit polynomial numerators by multiplying each original event's native numerator by the product of *all other* denominators; no potentially small gap is inverted after clearing. `C2_CLEARED_POLYNOMIAL_TARGET.json` contains these 24 complete polynomials and the zero rows. Its construction passed an exact rational cross-check against the original complete parent. This is an exact 24-term polynomial presentation, not a contact-state scalar approximation.

## 4. Exact compact dual positivity test, and the *first unresolved inequality*

For \(y\in\mathbb R^{412}\), define a linear functional \(L_y\) on the 412 monomials by their coordinates \(y_\kappa\). Its four **contact-truncated moment matrices** are
\[
 M_q(y)_{ij}=L_y(d_{q,i}d_{q,j}),\qquad
 q=1,\ldots,4.\tag{9}
\]
These are explicit integer-linear functions of \(y\), read directly from (6). Define the **source-independent compact spectrahedron**
\[
\mathcal Y=\Big\{y\in\mathbb R^{412}:
 M_q(y)\succeq0\ (q=1,\dots,4),\quad
 \mathcal N y=0,\quad \sum_q\operatorname{tr}M_q(y)=1\Big\}.
\tag{10}
\]
The gauge \(\mathcal Ny=0\) eliminates the six-dimensional kernel of \(\mathcal A^*\). If \(y\in\mathcal Y\), all four PSD matrices have trace sum one and hence bounded entries. Since \(y\mapsto\oplus_q M_q(y)\) is injective on \(\ker\mathcal N\), \(\mathcal Y\) is bounded, hence compact. It is nonempty, for example by taking a noncontact real evaluation functional, removing its \(\operatorname{im}\mathcal N^T\) component, and normalizing its positive trace.

Let
\[
\delta(\theta)=\min_{y\in\mathcal Y} L_y(f(\theta))
=\frac1{\mathscr D(\theta)}\min_{y\in\mathcal Y}
\underbrace{\sum_{\kappa\in\mathscr S}P_\kappa(\theta)y_\kappa}_{\textbf{24 explicit native polynomial terms}}.
\tag{11}
\]
The source-normalized Gram cone \(\mathcal K=\mathcal A(\bigoplus_q\mathbb S_+^{n_q-1})\) is **closed**. Indeed, integrate each difference-root Gram on a fixed product of a positive-angle compact rectangle and a real state sphere. Each difference-column is linearly independent within its parity block; therefore its integral moment matrix is positive definite. A convergent sequence of coefficient vectors bounds the traces of all PSD Gram blocks, and has a convergent Gram subsequence. Thus the limiting coefficient vector remains in \(\mathcal K\).

Finite-dimensional separation and (7) now give the exact equivalences
\[
\boxed{\begin{aligned}
 f(\theta)\in\mathcal K &\iff \delta(\theta)\ge0,\\
 \exists\ H_q\succ0\ \forall q:\ \mathcal A(H)=f(\theta)&\iff \delta(\theta)>0,\\
 f(\theta)\notin\mathcal K&\iff\delta(\theta)<0.
\end{aligned}}\tag{12}
\]
For the strict part: if \(f=\mathcal A(H)\) and \(H_q\succ0\), then \(L_y(f)=\sum_q\operatorname{tr}(H_qM_q(y))\ge\min_q\lambda_{\min}(H_q)>0\). Conversely, let \(e=\mathcal A(I,I,I,I)\). Every \(y\in\mathcal Y\) satisfies \(L_y(e)=1\). If \(\delta(f)>0\), choose \(0<\xi<\delta(f)\); then \(L_y(f-\xi e)>0\) for all \(y\in\mathcal Y\), hence \(f-\xi e\in\mathcal K\), so \(f=\mathcal A(H'+\xi I)\) with every block positive definite. The nonstrict equivalence uses the bipolar theorem and closedness proved above. The dual normalization is therefore exact, with no hidden Slater or image-closedness assumption.

**The exact unresolved global strict-C2 admission statement is**
\[
\boxed{\forall\theta\in\Omega\quad
\forall y\in\mathcal Y:\quad
\sum_{\kappa\in\mathscr S}P_\kappa(\theta)y_\kappa\ge0.}\tag{13}
\]
This is one completely explicit finite polynomial inequality involving only the fixed moment matrices and the 24 printed source polynomials. It is **not proved** here. Proving the stronger **strict** inequality in (13) for every \(y\in\mathcal Y\) would yield a nondegenerate Gram for all strict sources, and then compact-uniform source margins. Merely proving (13) nonstrictly would suffice for their fixed-\(L_6\) admission. An exact \(y\in\mathcal Y\) with negative pairing would instead be an exclusion certificate, directly checkable by `verify_exact_certificates.py`. No such y or positive universal Gram is claimed.

This is the first unproved PSD implication. There is **no remaining affine coefficient, missing test-normal, or kernel-matching condition** in the fixed sixth-root support face.

## 5. Strict physical angular curvature does not prove (13)

For strict C2 sources, the native determinant Laurent polynomial (SG5–SG6) has exactly six positive outer coefficients \(q_{01},q_{02},q_{10},q_{12},q_{20},q_{21}\). Their exponent vectors are
\[
E=\{(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0)\}.
\]
Let \(w_e=q_e U_c^{e_1}V_c^{e_2}>0\) with \(U_c=c^2,V_c=a^2c^2\), and put \(z=(\log(U/U_c),\log(V/V_c))\). The contact value and log-gradient vanish, so
\[
 L(U,V)=\sum_{e\in E}w_e(e^{e\cdot z}-1-e\cdot z).
\tag{14}
\]
Write \(A_*=w_{(-1,0)}+w_{(1,0)}\), \(B_*=w_{(0,-1)}+w_{(0,1)}\) and \(C_*=w_{(-1,1)}+w_{(1,-1)}\), all strictly positive. Then
\[
 \nabla^2 L(0)=\begin{pmatrix}A_*+C_*&-C_*\\-C_*&B_*+C_*\end{pmatrix},
 \qquad\det\nabla^2 L(0)=A_*B_*+C_*(A_*+B_*)>0.
\tag{15}
\]
For \(\|z\|\le R\) the elementary integral remainder inequality yields
\[
 L(U,V)\ge\tfrac12e^{-\sqrt2R}z^T\nabla^2L(0)z.
\tag{16}
\]
These inequalities are source-uniform on compact subsets of \(\Omega\). Together with SG3's strictly positive three-state block, they prove a unique positive-angle physical contact ray and nondegenerate transverse curvature. In particular real evaluation functionals on any compact strictly positive testing-angle window have a positive lower ratio between the parent and the contact-difference norm (extended across the contact by the positive transverse Hessian).

**But (15)–(16) control physical test/state evaluations, not every \(y\in\mathcal Y\).** The latter are positive *truncated pseudo-moments* for the 154 sixth-root features and need not possess a representing measure. Their separation is exactly the extra algebraic issue in (13). It would be incorrect to substitute the scalar determinant Hessian, a physical positivity bound, or rank 406 for that inequality.

## 6. End preservation, symmetry and precise scope

The previously independently verified 16-feature finite-\(c\) theorem proves **PSD** membership for sufficiently large \(c\) near each strict \(c\to\infty\) face point, uniformly on compact interior face subsets, and for the corresponding \(c\to0\) faces by the exact rational source inversion. The original files are copied byte-for-byte under `preserved-open-end/` in this package and SHA256 checksums are recorded.

The verified source inversion is
\[
 c'=1/(ac),\quad a'=a,\quad b'=a^2/b,\quad
 r'=a^2c^2r/b,\quad t'=[1-c^2(b-a^2t)]/b.
\tag{17}
\]
It permutes the exact native weights \((n,m,C,A,D,\alpha)'=(n,D,\alpha,A,m,C)\) and the complete four-state parent by
\(u^2v^2F_{\sigma\theta}(1/v,1/u;Z,X,Y,W)=F_\theta(u,v;W,X,Y,Z)\). Since this permutation maps sixth-root indices \((i,j)\to(6-j,6-i)\), it preserves the **whole fixed** \(Z_6\) library. For interior admission one may therefore restrict to \(ac^2\ge1\), applying inversion to any source with \(ac^2<1\). This reduces the parameter range but does not prove (13).

**Warning about nondegenerate continuation.** The verified 16-feature open-end Grams, embedded by zero rows in the 154-feature support basis, are generally singular after removing the four contact vectors. Thus the open-end theorem proves \(\delta(\theta)\ge0\) there, **not automatically** \(\delta(\theta)>0\). It cannot be turned into an all-interior nondegeneracy continuation argument by openness and closedness alone.

The `L_2` strict rational exclusion in `independent-review/C2_HALF_POWER_OBSTRUCTION.md` additionally rules out a universal 16-feature or half-power certificate. No such claim is made. The entire strict native C2 family, its source walls, and the global hybrid lift remain open for **explicit** finite-degree admission until (13), or another concrete globally sound certificate, is proved. The separately reviewed qualitative spectrahedral-shadow existence assertion is not used as an explicit extraction.

## 7. Original deterministic files

- `c2_fixed_map.py`: actual native dictionary (2), the full original parent (3), exact 154 support features, all 412 source-independent monomials, 3021 integer coefficient columns and six contact rows. No external files.
- `verify_native_jets.py`: directly checks six full-state/angle native source identities, with a monomial-cleared denominator. Verified **six exact zero polynomials** in this run.
- `verify_rank_mod65521.py` and `C2_RANK406_MODULAR_PIVOTS.json`: reproduces the exact rank406 modular minor and records every selected original row/column (sufficient for a rational right inverse).
- `C2_154_CONTACT_MAP.json`: every original fixed monomial support row and sparse coefficient column. These are complete, not pseudocode.
- `generate_cleared_targets.py` and `C2_CLEARED_POLYNOMIAL_TARGET.json`: exact 24 integer-polynomial source coefficients and the positive common denominator. A deterministic rational-source cross-check confirmed all 24 active entries.
- `verify_exact_certificates.py`: validates any proposed rational PSD primal Gram or dual pseudo-moment witness, with exact full-parent coefficient matching and exact zero/positive rational Schur pivots. It **does not search** for a certificate.
- `generate_affine_preimage.py`: finite source-rational affine solution via the saved 406-minor. It is marked explicitly **AFFINE ONLY / NOT PSD**.
- `verify_angular_curvature.py`: checks the exact strict hexagon Hessian and its positive determinant, not a Gram admission.

### Bottom line

**Proved here:** exact full-parent affine rational reachability for every strict native C2 parameter, contact kernels/normal rows and their six native equations, source-polynomial target clearing, rank406 data, exact fixed compact dual moment criterion, strict physical angular curvature, and preservation of the independently verified open-end result. **Not proved:** inequality (13) for all strict C2 parameters, an all-interior 154-feature PSD Gram, or global \(L_6\) equality. The first genuine unresolved implication is the nonnegative 24-term native polynomial pairing over the compact dual moment spectrahedron \(\mathcal Y\).
