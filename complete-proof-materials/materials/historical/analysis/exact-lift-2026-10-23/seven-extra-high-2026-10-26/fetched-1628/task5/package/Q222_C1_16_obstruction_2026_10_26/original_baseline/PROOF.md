# Native C1: exact fixed-feature coefficient map, source gauges and PSD frontier

**26 October 2026.** This is an original analytical construction from the genuine physical C1 coquotients SD1--SD11 and the complete SG2/SD2 four-state parent. It proves an exact *universal coefficient-matching theorem* in 16 fixed sixth-root-compatible features and identifies an eight-scalar, source-only PSD completion. It **does not prove** that the eight-scalar completion is positive for every strict C1 native source. Consequently neither complete strict C1 admission to `Z_6` nor an evaluated Q222 global lift follows here. No SDP solver, numerical sampling, parameter sweep or GitHub write was used.

## Theorem A. Exact native strict domain and stationary parent

Use five positive real parameters `k,f,g,a,t`. Put

\[
\begin{gathered}
 \kappa=k^2,\quad\eta=f^2,\quad\zeta=g^2,\quad e=kf/g,\quad\Delta=e^{-1},
 \quad V_c=\kappa\eta\zeta,\\
 L_c=1-a\eta,\quad H_c=1+(1-a)\eta,\quad Q_c=1-\kappa\zeta+aV_c,\\
 R=1+\kappa(1-a\eta),\quad T=1+\kappa-\kappa\zeta t,\quad
 M=R-t,\quad J=T-a.
\end{gathered}
\]

The exact polynomial definitions of `P4`, `Palpha`, `Pd` and the seven native coquotients `c2, cb, c3, c5, c4, calpha, cd` are in `c1_native.py`; they are reproduced verbatim from the source. Work on the **physical strict C1** locus where

\[
\begin{gathered}
L_c,H_c,Q_c>0,\quad 1-t\zeta>0,\quad 1-\kappa\eta+tV_c>0,\\
P_4,P_\alpha,P_d>0,\quad M,J>0,\\
 2c_5/J-\eta c_d>0,\qquad 2e^2c_b/t-\eta c_d>0.
\end{gathered}\tag{A1}
\]

In particular `R>t>0`, `T>a>0`, `c2,cb,c3,c5,c4,calpha,cd>0`, and none of the source denominators displayed below vanishes. No arbitrary positive choices of these coefficients are allowed. The original Bell tuple is

\[
\begin{gathered}
(d,\lambda,n,m,C,A,D,\alpha,\beta)
=\bigl(d_B,0,0,m,2,A,D_e,\alpha,\beta\bigr),\\
d_B=2a\eta c_d/c_3,
\qquad m=2ae^2c_2/(Mc_3),\quad A=2a\eta^2c_4/c_3,\\
\alpha=2a\Delta^2k^6\eta^2c_\alpha/c_3,\quad
D_e=2ac_5/(Jc_3)-d_B/2,\quad
\beta=2ae^2c_b/(tc_3)-d_B/2.
\end{gathered}\tag{A2}
\]

Every non-omitted coefficient is strictly positive by (A1). The contact is

\[
\boxed{u_c=k,\quad v_c=kfg,\qquad
(W_c,X_c,Y_c,Z_c)=\bigl(g/f,-k(1-tg^2),\,kgL_c/f,1\bigr).}\tag{A3}
\]

All six entries are nonzero, with `Xc<0` and the other state entries positive. The complete SD2 source-normalized parent is **identically**

\[
\Phi=\Delta^2\Lambda F_c,\qquad \Lambda=RTc_3/(2a)>0,
\tag{A4}
\]

without restricting a state coordinate or testing angle. `verify_parent_scale.py` proves (A4) term-by-term in all seven original event/core coefficients. `verify_c1_native_stationarity.py` differentiates the original SD2 four-state parent and proves, as six exact multivariate identities using the seven native coquotients,

\[
 \partial_W\Phi=\partial_X\Phi=\partial_Y\Phi=\partial_Z\Phi
 =\partial_u\Phi=\partial_v\Phi=0\quad\text{at (A3)}.\tag{A5}
\]

This verifies the state and angular contact of the *physical* parent. The established stationary-to-global-validity theorem (SG) makes this source quantum-valid on the native strict locus. The present problem is strictly stronger: a fixed-source Gram admission.

## Theorem B. Universal 16-feature **exact coefficient matching**

Use the **original fixed** testing/state library

\[
\begin{aligned}
 z_o&=(uW,vW,X,uvX,Y,uvY,uZ,vZ),\\
 z_e&=(uvW,vX,uY,vY,Z,uvZ),\\
 z_h&=(\sqrt{uv}\,W,\sqrt{uv}\,Z).
\end{aligned}\tag{B1}
\]

Every row is a member of the full `Z_6` feature library, with testing exponents in `{0,1/2,1}`. Its definition does **not** depend on a source or state contact. The labels `o,e` refer to parity under `(u,v;X,Y) -> (-u,-v;-X,-Y)` (leaving `W,Z` fixed), an exact symmetry of the whole parent.

To derive an explicit source-rational contact Gram, make only **invertible source-side** substitutions

\[
u=u_cU,\quad v=v_cV,\qquad
(W,X,Y,Z)=\operatorname{diag}(W_c,X_c,Y_c,Z_c)\phi.\tag{B2}
\]

Normalize each feature of (B1) by its nonzero contact evaluation. In the normalized variables `U,V,phi`, both integer feature vectors and both half-power features take the value `1` at the standard contact `(1,1;1,1,1,1)`.

Introduce exactly twelve difference roots, indexed from zero:

\[
\begin{aligned}
 r_0,\ldots,r_6&=(UW-VZ,\ VW-VZ,\ X-VZ,\ UVX-VZ,\\
 &\hspace{20mm}Y-VZ,\ UVY-VZ,\ UZ-VZ),\\
 r_7,\ldots,r_{11}&=(UVW-UVZ,\ VX-UVZ,\ UY-UVZ,\\
 &\hspace{26mm}VY-UVZ,\ Z-UVZ),\\
 h&=\sqrt{UV}(W-Z).
\end{aligned}\tag{B3}
\]

Here all `W,X,Y,Z` in (B3) denote the normalized state variables `phi`. Define the fixed linear integer coefficient map

\[
\boxed{\mathcal A(H_7,H_5,g)
  =\operatorname{coeff}_{42}\left(
  r_{0:6}^T H_7r_{0:6}+r_{7:11}^TH_5r_{7:11}+gh^2\right).}
\tag{B4}
\]

The `42` coefficients are the invariant products `U^iV^j phi_a phi_b` in the full bidegree `(2,2)` support, with parity
`i+j+1_{a in {X,Y}}+1_{b in {X,Y}}` even. All 42 keys are regenerated in `c1_small_map.py`. The map has `28+15+1=44` entries.

**Exact integer rank:**

\[
\boxed{\operatorname{rank}_{\mathbb Q}\mathcal A=36,
\qquad \dim\ker\mathcal A=8.}\tag{B5}
\]

A fully specified reproducible minor uses the lexicographically first 36 independent columns and, among them, the first 36 independent rows. Their zero-based indices appear in `verify_c1_map_exact.py` and `verification.log`; its **integer determinant is exactly `-2^24`**. Hence the rank is at least 36.

Every difference root in (B3) vanishes at the standard contact, so its square and all mixed products annihilate the four state derivatives and the two testing derivatives. These are six linearly independent functionals on the 42-dimensional coefficient space: their lexicographic `6x6` minor indexed by `(0,1,3,4,6,14)` has determinant `4`. Thus the rank is at most `42-6=36`, proving (B5) without a genericity assumption.

It follows that

\[
 \boxed{\operatorname{im}\mathcal A
 =\{b\in\mathbb R^{42}:\partial_{\phi_W,\phi_X,\phi_Y,\phi_Z,U,V} b(1,1;\mathbf1)=0\}.}\tag{B6}
\]

The native parent (A2) obeys these six equations by (A5); it has the required parity and 42-term support. Therefore **for every strict native C1 parameter tuple** there are real symmetric (not yet necessarily PSD) matrices with the universal *whole-parent* identity

\[
\boxed{
 F_{c(\theta)}(u_cU,v_cV;\operatorname{diag}(\psi_c)\phi)
 =r_o^TH_7(\theta)r_o+r_e^TH_5(\theta)r_e+g(\theta)h^2.}
\tag{B7}
\]

The full four-state, all-testing coefficient identity is exact, not a Taylor equality or a contact evaluation.

### A deterministic, all-parameter, source-rational solution

Choose the fixed `36x36` invertible integer minor above. Set the 8 nonpivot Gram coordinates to zero, solve those 36 equations by its fixed **rational inverse**, and write the 44-vector solution as `h_part(theta)`. Every entry is a rational function of the five native source parameters, since the only operations are (A2), (B2), coefficient extraction and multiplication by this one constant rational matrix. It is defined throughout the strict domain (A1); no testing value is inverted. The full image equation follows automatically from (B6), not a source-specific rank computation.

All exact solutions are

\[
\boxed{h(\theta;s)=h_{\rm part}(\theta)+\sum_{j=0}^7s_jK_j,
\qquad s_j\in\mathbb R,}\tag{B8}
\]

where the following source-independent sparse kernel vectors `K_j` are stated by their **only nonzero symmetric Gram coordinates** (write `Hij` for a matrix entry on the 12 difference roots and `g` for the odd square):

| Gauge | Nonzero entries |
|---|---|
| 0 | `H05=-1, H13=1, H16=-1, H56=1, H78=-1, H79=1` |
| 1 | `H13=1, H15=-1, H78=-1, H7,10=1` |
| 2 | `H13=1, H16=-1, H78=-1, H7,11=1` |
| 3 | `H25=1, H56=-1, H89=-1, H8,11=1` |
| 4 | `H34=1, H45=-1, H89=-1, H9,10=1` |
| 5 | `H34=1, H46=-1, H89=-1, H9,11=1` |
| 6 | `H34=1, H56=-1, H89=-1, H10,11=1` |
| 7 | `H01=-1/2, H16=1/2, g=1` |

All matrix indices are zero-based, and symmetric off-diagonal positions mean both `Hij` and `Hji`. The exact checker regenerates the kernel, so no one must trust a displayed nullspace vector or stored binary. In the lexicographic particular solution the odd `g` is zero; **the eighth gauge is exactly the C2-style odd source-only Gram gauge**. Seven other gauges change off-diagonal connections within the two reduced parity blocks. No gauge changes any of the twelve reduced diagonal entries.

## Theorem C. Forced positive masses and three genuinely strict Schur blocks

Write the original strict source coefficients as `(d, m,C,A,D,alpha,beta)` and the native contact as `(u_c,v_c;W_c,X_c,Y_c,Z_c)` from (A3). Every PSD or indefinite exact reduced Gram (B7) has the following **forced**, positive diagonal coefficients:

| Difference root index | Fixed diagonal |
|---|---|
| 0 | `2 uc^2 Wc^2(2C+2D+d)` |
| 1 | `4m vc^2 Wc^2` |
| 2 | `2Xc^2(2beta+d+2m)` |
| 3 | `2uc^2 vc^2 Xc^2(2D+d)` |
| 4 | `2Yc^2(2A+2C+2D+d)` |
| 5 | `2uc^2 vc^2 Yc^2(2A+2m+2alpha+d)` |
| 6 | `2uc^2 Zc^2(2alpha+2beta+d+2m)` |
| 7 | `4C uc^2 vc^2 Wc^2` |
| 8 | `2vc^2 Xc^2(2beta+d)` |
| 9 | `4A uc^2 Yc^2` |
| 10 | `2vc^2 Yc^2(2A+2C+2alpha+d)` |
| 11 | `4alpha Zc^2` |

The proof is exact coefficient comparison: for each node `i`, its **own monomial square is the only contribution** to the listed angular/state monomial of the full parent. Thus its diagonal is uniquely fixed. Every entry is strictly positive throughout (A1), including simultaneous approaches toward a wall **at a still-strict source**. This does not assert a uniform positive lower bound as source gaps vanish.

Three off-diagonal entries are also uniquely fixed by single mixed monomials:

\[
\begin{split}
H_{12}&=-4m v_cW_cX_c,\\
H_{04}&=-2(2C+2D+d)u_cW_cY_c,\\
H_{03}&=2(2D+d)u_c^2v_cW_cX_c.
\end{split}\tag{C1}
\]

Their entire forced `2x2` principal determinants are, respectively,

\[
\boxed{\begin{aligned}
\det H_{\{1,2\}}&=8m(2\beta+d)v_c^2W_c^2X_c^2>0,\\
\det H_{\{0,4\}}&=8A(2C+2D+d)u_c^2W_c^2Y_c^2>0,\\
\det H_{\{0,3\}}&=8C(2D+d)u_c^4v_c^2W_c^2X_c^2>0.
\end{aligned}}\tag{C2}
\]

`verify_forced_positivity.py` establishes these as rational polynomial identities for arbitrary source weights and contacts, **not** at chosen witnesses. None of the eight gauges modifies these entries. The three strict inequalities survive every strict native opening; they do not control the other principal minors.

## Exact PSD criterion, soundness, and first unproved implication

Define the **one fixed affine SDP cone** from the original source coefficients and original fixed 16 features (B1), without contact-dependent coefficients:

\[
\mathcal L_{16}^{\rm C1}=\bigl\{c:\ F_c(u,v,\psi)=Z_{16}(u,v,\psi)^TGZ_{16}(u,v,\psi),\quad G\succeq0\bigr\}.
\tag{D1}
\]

The identity means finite coefficient matching of all integer and fractional powers; it is affine in the original Bell source coefficients. It is globally sound: a PSD Gram gives nonnegativity for independent `u,v>0` and all four real states, the exact original parent sign congruences give other quadrants, homogeneous limits cover endpoints, and real/imaginary splitting covers complex states. There is no source-dependent testing feature.

For a strict native C1 source, (D1) is **equivalent** to the existence of eight real source-only gauge values in (B8) satisfying

\[
\boxed{H_7(\theta;s)\succeq0,\qquad
 H_5(\theta;s)\succeq0,\qquad g(\theta;s)=s_7\ge0.}
\tag{D2}
\]

**Necessity proof:** fractional-versus-integer angular monomial parity lets one delete integer/half-power Gram cross blocks while retaining their positive principal Grams. On the integer features, averaging the exact simultaneous sign symmetry splits odd and even blocks without changing the whole parent or positivity. Since the contact is a common zero, every PSD block kills its respective contact-evaluation vector. Source-only normalization (B2) and subtraction of anchor features then give the 7- and 5-dimensional contact differences, while the odd `2x2` block becomes `s_7 h^2`. The complete coefficient solution must have the form (B8) because the fixed map has exactly eight kernel directions. Sufficiency is immediate from the block PSD identity (B7).

**The first unproved implication for this explicit 16-feature construction is precisely**

\[
\boxed{\forall\theta\text{ satisfying the complete strict native inequalities (A1)},\quad
 \exists s\in\mathbb R^8\text{ satisfying all of (D2)}.}\tag{OPEN-C1-16}
\]

The twelve positive diagonals (Theorem C), three positive `2x2` minors, six source stationary equations and exact rational coefficient matching do **not** imply (OPEN-C1-16). There is no proven universal choice of gauge functions making both full blocks positive semidefinite. Our analytic construction therefore **does not prove** all strict C1 sources belong to this cone or even to the larger `Z_6` cone. Failure of (OPEN-C1-16), if eventually demonstrated, would disprove only this 16-feature *subcone*, not sixth-root admission in the full 196-feature library. A different globally positive fixed feature alternative remains possible.

The established closedness of the fixed `Z_6` cone can pass *proved membership* to limiting native faces, but cannot manufacture missing strict-source PSD certificates. Constant contact-map rank plus local positivity establishes openness only. Neither source connectedness nor compact-subset grids supplies the required closed-nondegeneracy implication for (OPEN-C1-16).

## Reproduction and original evidence

`c1_native.py` reproduces SD1 and SD2; `verify_parent_scale.py` proves SD2=positive source multiple of the original canonical full parent; `verify_c1_native_stationarity.py` verifies the original native contact through all six exact derivatives. `c1_small_map.py` lists the 14+2 fixed features and constructs the integer Gram coefficient map. `verify_c1_map_exact.py` proves the ranks, two exact nonzero minors and the universal all-contact coefficient identity. `verify_forced_positivity.py` verifies all twelve positive forced masses and three Schur determinants. `eight_gauge_constructor.py` supplies a deterministic exact formula for any input native source and any eight rational gauges, with no positivity/optimization claims. `verify_all.py` replays every completed exact test and writes `verification.log`.

The uploaded source text, original C1 TeX, and the separately verified C2 strict-open-end mechanism are evidentiary inputs; none is taken as proof of (OPEN-C1-16). The full qualitative tensor-SOS/global-conjunction proof and the other chat's C2 interior work are outside this assignment.
