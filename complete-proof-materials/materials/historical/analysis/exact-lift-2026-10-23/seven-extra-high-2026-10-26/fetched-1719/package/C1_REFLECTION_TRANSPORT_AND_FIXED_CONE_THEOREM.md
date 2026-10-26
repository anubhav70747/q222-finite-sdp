# Q222: exact C1 orientation, Bell reflection, and one common fixed-cone entrance

**Original targeted theorem, 26 October 2026.** Baseline: the attached `C1_REFLECTION_SOURCE_PACKET.txt`, notably SG1–SG2, EN1–EN16, the independently corrected C1 obstruction, the exposed-contact-ray theorem, the published singular-density/closed-sum result, and the rational strict-feasibility construction. **No C1 or C2 strict-family PSD assertion is established here.** The proofs below are finite coefficient-space identities, not local contact certificate searches.

## Theorem (conditional explicit finite global cone)

Let `Z16` consist, in the indicated order, of

```
uW, vW, X, uvX, Y, uvY, uZ, vZ,
uvW, vX, uY, vY, Z, uvZ, sqrt(uv)W, sqrt(uv)Z.
```

Assume only these **two pending construction inputs**:

* **Oriented C1 admission:** every genuinely physical **strict** C1 native source obeying the original strict A1 conditions and `0 < k = u_c <= 1` has a PSD Gram in this fixed `Z16` library (after its positive overall C=2 normalization). The claim for unrestricted `k` is false, as proved by the published `k=4` witness. No use of that false claim is made.
* **Strict C2 admission:** every genuinely physical strict six-event C2 source obeying the complete native domain has a PSD Gram in the original full 196-row `Z6` sixth-root library. This is the unproved universal strict C2 dual inequality, not the already established affine rank result.

Set

\[
\begin{split}
 Z_{6,12}(u,v;\psi)&=\bigl(u^{i/6}v^{j/6}\psi_l\bigr)_{0\le i\le6,\,0\le j\le12,\,l=W,X,Y,Z},\\
 \mathcal K_{6,12}&=\left\{c\in\mathbb R^9_+:\ (1+v^2)F_c(u,v;\psi)
      =Z_{6,12}^{T}PZ_{6,12}\text{ identically, for some }P\succeq0\right\}.
\end{split} \tag{1}
\]

`F_c` is the complete original nine-coefficient SG2 parent with four unrestricted state amplitudes. Every identity in (1) means **coefficient matching** as a polynomial in `u^(1/6),v^(1/6)` and a quadratic polynomial in all four state variables. No contact or source normalization enters this definition. It is a fixed, globally sound, **closed** rational-data affine-SDP cone, with **364** monomial rows and at most **3250** coefficient equations.

**Conclusion:** Under the two stated conditional inputs, **every real strict C1 source in either orientation, and every real strict C2 source, belongs individually to the same cone** `K_{6,12}`. Consequently, by the already published five-event-contact density theorem and the original SG9–SG14 source classification, its finite physical symmetry orbit together with the existing closed-face cones SG12/SG13 and event rays gives exactly the complete quantum-valid Bell source cone. The finite sum is closed. The corresponding dual intersection is a single finite rational SDP lift of Q222, with the target coefficient adjoint and source orthant both retained. None of these conclusions assumes membership of the reflected C1 sources in the *unchanged* `Z6` cone.

## 1. Bell/source reflection is a fixed integral linear map

Use source coordinates `c=(d,lambda,n,m,C,A,D,alpha,beta)`. In the exact SG1 basis the EN10 reflection is

\[
 \boxed{\sigma c=(d,n,\lambda,C,m,\alpha,\beta,A,D),\qquad \sigma^2=I.} \tag{2}
\]

Thus the mandatory C1 omissions `lambda=n=0` are preserved, and `D,beta` are interchanged. Positive retained events remain positive. In conventional affine Bell coefficient order

`(constant,A0,A1,B0,B1,E00,E01,E10,E11)`, it acts by the following signed permutation:

\[
 (t,a_0,a_1,b_0,b_1,e_{00},e_{01},e_{10},e_{11})
 \longmapsto(t,a_0,-a_1,-b_1,-b_0,-e_{01},-e_{00},e_{11},e_{10}). \tag{3}
\]

Here `E_SG1` is the fixed 9-by-9 coefficient map from `f=dH+4 sum a_i p_i`; (3) is exactly `E_SG1 sigma E_SG1^{-1}`, whose determinant and involution are checked in the included script. The dual action on behaviours is the transpose of (3), equal to the same signed permutation in this basis. No nonlinear source transformation occurs in a source cone.

Define the complete state congruence

\[
 G(v)=\begin{pmatrix}v&-1&0&0\\-1&-v&0&0\\0&0&-v&1\\0&0&1&v\end{pmatrix},
 \quad G(v)^T=G(v),\quad G(v)^2=(1+v^2)I.\tag{4}
\]

A direct independent expansion of all nine SG2 coefficients proves the **entire four-state polynomial identity**

\[
 \boxed{u^2F_{\sigma c}(1/u,v;G(v)\psi)=(1+v^2)F_c(u,v;\psi).}\tag{5}
\]

It holds before imposing physicality, contact, omitted weights or source scaling. Because `G(v)` is invertible for real finite `v`, the inverse identity follows by applying it twice. For `u>0`, positive-angle canonical physicality is preserved; the homogeneous/projective endpoints follow by the polynomial identity as described in §4.

## 2. Exact strict-native parameter transport, including C=2 restoration

Use the original C1 parameters `(k,f,g,a,t)`, with positive 

\[
\begin{gathered}
\kappa=k^2,\quad \eta=f^2,\quad\zeta=g^2,\quad r=1-a\eta,\quad h=1-t\zeta,\\
H=r+\eta=1+(1-a)\eta,\quad G_0=h+\zeta=1+(1-t)\zeta,\\
R=1+\kappa r,\quad T=1+\kappa h,\quad M=R-t,\quad J=T-a,\\
E=1-\kappa\eta h,\quad Q=1-\kappa\zeta r,\quad V_0=\kappa\eta\zeta.
\end{gathered}\tag{6}
\]

The actual A3 contact is

\[
(u_c,v_c;W_c,X_c,Y_c,Z_c)=(k,kfg;g/f,-kh,kg r/f,1).
\]

Its reflected state has the **exact** nonzero positive last coordinate

\[
 G(kfg)\psi_c=(kG_0,-(g/f)E,Q,kgH/f),\quad
 s_c=kgH/f>0. \tag{7}
\]

Divide the reflected state by `s_c`; this is a positive **source-only** state normalization, not a testing-dependent Gram coefficient. Its four normalized entries are

\[
 \boxed{\psi'_c=\left(\frac{fG_0}{gH},-\frac{E}{kH},\frac{fQ}{kgH},1\right).}\tag{8}
\]

Applying the exact inverse native map EN13 to this contact gives

\[
\boxed{\begin{aligned}
 k'&=1/k,& \eta'=f'^2&=\kappa\zeta H/G_0,&\zeta'=g'^2&=\kappa\eta G_0/H,\\
 a'&=M/(\kappa H),& t'&=J/(\kappa G_0),&v'_c&=v_c.
\end{aligned}}\tag{9}
\]

Only positive square roots `f',g'` are chosen. Direct multiplication verifies `k'f'g'=kfg`. The **complete opening-safe native sign dictionary** includes

\[
\begin{array}{ll}
 r'=1-a'\eta'=Q/G_0>0,&h'=1-t'\zeta'=E/H>0,\\
 H'=r'+\eta'=(1+V_0)/G_0>0,& G'_0=h'+\zeta'=(1+V_0)/H>0,\\
 E'=1-\kappa'\eta'h'=h(1+V_0)/G_0>0,&
 Q'=1-\kappa'\zeta'r'=r(1+V_0)/H>0,\\
 M'=1+\kappa'r'-t'=a(1+V_0)/(\kappa G_0)>0,&
 J'=1+\kappa'h'-a'=t(1+V_0)/(\kappa H)>0,\\
 V'_0=\kappa'\eta'\zeta'=V_0,&
 1-\kappa'\eta'+t'V'_0=h(1+V_0)/G_0>0.
\end{array}\tag{10}
\]

Here `kappa'=1/kappa`. More explicitly, `G_0-Q=zeta*M>0` and `H-E=eta*J>0` show the **strict upper bounds** `0<r'<1` and `0<h'<1`, not merely positivity. Every denominator in (9)–(10) is strictly positive in the original **physical** native domain. The remaining native coquotient signs (`P4',Palpha',Pd'` and both event-weight numerator signs) follow from the strictly positive *reflected original weights*, the positive gaps (10), and the stationary-uniqueness theorem EN8/SD1. Thus (9) lies in **the entire actual oriented A1 domain**, not an artificial k-truncated source chart. No sign is inferred merely from a finite limit.

**Crucial scale convention.** The native formula SD1 sets `C=2`, whereas (2) changes the new C weight to the old `m`. If the original representative has `C=2`, the **correct** C-normalized reflected representative is

\[
 \boxed{\widehat{\sigma c}=(2/m)\,\sigma c.} \tag{11}
\]

The factor is strictly positive and pure source. It is used **only to invoke the oriented native theorem**; source cones and their semidefinite pencils always use the **linear** involution (2), never the nonlinear normalization (11). The checker independently compares all nine exact SD1 source coefficients under (9),(11) at the published strict rational `k=4` source, confirming this easily missed factor.

If `k<=1`, use the source itself; if `k>1`, (9) gives `0<k'<1`. Therefore just the identity orientation and the single reflection cover **every strict C1 ray**, including `k=1`. The `k=4` Z16 obstruction consequently does *not* obstruct the conditional oriented statement: it makes no assertion about `sigma c` or the larger cone below.

## 3. Fixed transport of the PSD Grams, without an illegitimate Z6 inference

Let `z16` be the displayed ordered sixteen-row column. Write

\[
\begin{gathered}
 W'=vW-X,\quad X'=-W-vX,\quad Y'=Z-vY,\quad Z'=Y+vZ,\\
 h_{16}(u,v;\psi):=u\,z_{16}(1/u,v;G(v)\psi).
\end{gathered}\tag{12}
\]

All sixteen fixed rows are, in original testing/state coordinates,

\[
\boxed{\begin{array}{rrrr}
 W',&uvW',&uX',&vX',\\
 uY',&vY',&Z',&uvZ',\\
 vW',&uvX',&Y',&uvY',\\
 uZ',&vZ',&\sqrt{uv}\,W',&\sqrt{uv}\,Z'.
\end{array}}\tag{13}
\]

The display is a *fixed* feature column, not an additional parametrization. Its integer rows use only `u^0,u^1` and `v^0,v^1,v^2`; its two fractional rows use `u^(1/2)v^(1/2)` and `u^(1/2)v^(3/2)`. Hence each is an **integer-coefficient linear combination** of the 364 fixed `Z_{6,12}` monomials. The included JSON gives the full sparse matrix `L in Z^{364 x 16}` with exactly 32 nonzero signed entries, so `h16=L^T Z_{6,12}`.

The reflection identity proves the equivalence for every canonical source:

\[
 F_{\sigma c}=z_{16}^TQz_{16}
 \quad\Longleftrightarrow\quad
 (1+v^2)F_c=h_{16}^TQh_{16}, \qquad Q\succeq0. \tag{14}
\]

The reverse implication follows by replacing `u` by `1/u` and `psi` by `G(v)psi/(1+v²)`, using `G(v)^2=(1+v²)I`. Thus the reflection of the oriented 16-feature cone is exactly the **closed fixed** cone cut out by `(1+v²)F_c=h16^TQh16`. It is **not claimed** to be contained in the original 196-row `Z6` cone: `v²` and `v^(3/2)` occur in (13), and the factor `(1+v²)` must not be discarded.

For any full `Z6` representation `F_c=z6^TQz6`, there is the exact identity

\[
 (1+v^2)F_c=z6^TQz6+(v z6)^TQ(v z6).
 \tag{15}
\]

Both columns `z6` and `v z6` also lie in the fixed `Z_{6,12}` span. Formula (15) embeds **every** `Z6` PSD source into `K_{6,12}`. Formula (14) embeds every reflected oriented C1 PSD source into `K_{6,12}`. This proves the conditional strict-source coverage asserted in the theorem, with *one cone, one testing multiplier and one fixed finite rational coefficient map*. No mixed source/testing square root is taken. Over independently labelled real-closed source and test fields, the fixed Gram equality transfers by first-order real-closed-field transfer. Factoring the source-only PSD Gram and dividing each resulting separated square by the **pure-testing** positive root `sqrt(1+v²)` gives a literal separated tensor SOS for `F_c`; the multiplier is not inverted in a mixed field.

## 4. Global soundness, literal axes, projective endpoints, and closedness

Let `R=u^(1/6)`, `S=v^(1/6)` only as formal independent positive test variables. Equality in (1) is a finite polynomial identity in `(R,S)` and `psi` with rational coefficients linear in `(c,P)`. It imposes at most
`13 * 25 * 10 = 3250` real linear coefficient equations, one for each exponent/state-symmetric quadratic key. Gram off-diagonal entries contribute their correct factor **two**. This is a fixed *affine linear* projection of `S_+^364 x R_+^9`. The native parameters in (9) do not appear in its SDP variables.

For any feasible Gram, `(1+v²)F_c >= 0` for all positive finite `u,v` and every unrestricted real state. The multiplier is positive and **pure testing**, hence `F_c>=0`. The two exact SG2 sign congruences transport this to every finite real testing quadrant. At either literal axis, the fractional rows have well-defined zero values and polynomial continuity gives `F_c>=0`. At the projective endpoints, homogenize the original parent by dividing by `u²`, `v²`, or `u²v²` before taking the corresponding reciprocal limit. Because every row in (13) and every `Z_{6,12}` monomial has `u` exponent <=1 and `v` exponent <=2, the certificates also have exact finite leading-square limits. In particular at Bob infinity

\[
 \lim_{v\to\infty}\frac{F_c(u,v)}{v^2}
 =\lim_{v\to\infty}\frac{Z_{6,12}^TPZ_{6,12}}{v^4}
 =\left(\lim_{v\to\infty}v^{-2}Z_{6,12}\right)^TP
 \left(\lim_{v\to\infty}v^{-2}Z_{6,12}\right)\ge0.
\]

The corresponding `u` and double-endpoint limits are identical. Real/imaginary parts of a complex state split because SG2 has real coefficients. Thus `K_{6,12}` is globally sound, including all boundary tests.

For closedness, integrate `Z_{6,12} Z_{6,12}^T` over a fixed positive testing rectangle `[1,2]^2` and the real unit state sphere `S³`. The resulting `364 x 364` moment matrix `M_ext` is strictly positive definite, because its 364 analytic monomial/state functions are linearly independent. If `c_j->c` and `P_j>=0` represent `(1+v²)F_{c_j}`, then

\[
 \lambda_{\min}(M_\mathrm{ext})\operatorname{tr}P_j
 \le \langle M_\mathrm{ext},P_j\rangle
 =\int(1+v²)F_{c_j},
\]

whose right side is bounded and linear in `c_j`. A Gram subsequence converges and preserves every identity, as do the orthant constraints. Thus **`K_{6,12}` is closed in original Bell coefficient space**, even when the *native contact chart and a chart-normalized Gram diverge*. Source-dependent contact limits are not used in this argument.

## 5. Exact hybrid assembly and the source-dual adjoint, including the orthant

The separately reviewed density theorem now applies **without change** to the single fixed closed `K_{6,12}`: under the two pending strict PSD assumptions, the permitted five-event C1/C2 contact rays lie in this same cone by coefficient-space closure. The SG9/SG10 classification and EN1–EN16 source entry leave no extra canonical extreme-contact type. The published contact-ray exposure theorem requires each extreme ray to be individually admitted by an actual retained cone; this is satisfied by §3 and density, *not* inferred from a decomposition into other contact rays.

Let `E_event` denote the finite event-ray cone, `F12,F13` the two SG12/SG13 closed-face cones, and `Gamma` the finite physical source-relabeling group of SG8. The exact conditional coefficient-space assertion is

\[
 \boxed{\mathcal P_{\rm Bell}=
 E_{\rm event}+\sum_{g\in\Gamma}
        g\bigl(\mathcal K_{6,12}+F_{12}+F_{13}\bigr).}\tag{16}
\]

All summands are globally sound and separately **closed**. Their finite sum is also closed *for this Bell problem*: every sound affine functional `f=t+sum_i b_i chi_i` has `t>=0` and `|b_i|<=t`, by averaging over the 16 deterministic assignments. A bounded sum therefore bounds each sound summand, and termwise subsequences converge. Equation (16) uses no false general principle that sums of closed cones are automatically closed. Literal coefficient-zero faces are preserved under each nonnegative decomposition.

For explicit dual assembly, index the `3250` coefficient/state keys and define **fixed rational linear maps**

\[
 \mathsf T c=\operatorname{coeff}[(1+v²)F_c],\qquad
 \mathsf A P=\operatorname{coeff}[Z_{6,12}^T P Z_{6,12}].
\]

The actual canonical source cone is **not** merely `A P = T c, P>=0`; its source-orthant constraints must be included:

\[
 \mathcal K_{6,12}=
 \{c\in\mathbb R_+^9:\exists P\succeq0,\ \mathsf T c=\mathsf A P\}. \tag{17}
\]

Its exact *attained* dual is

\[
 \boxed{\mathcal K_{6,12}^*=
 \bigl\{b\in\mathbb R^9:\exists y\in\mathbb R^{3250},\ r\in\mathbb R_+^9,
 \quad b=\mathsf T^Ty+r,\quad \mathsf A^*y\succeq0\bigr\}.}\tag{18}
\]

**Proof of equality and attainment, not a shorthand.** One inclusion follows from
`b.c = <A* y,P>+r.c >=0` whenever `T c=A P`. Conversely, the Gram image `A(S_+^364)` is closed by the same positive-feature-moment trace argument. If a source `c` is nonnegative on the right side of (18), choosing `y=0` forces `c>=0`; varying `y` with `A*y>=0` forces `T c` into `A(S_+)` by finite-dimensional separation. Thus the dual of the right side of (18) is exactly `K_{6,12}`. It remains to show the right side is closed.

The **previously published exact rational strict point** `c*=(1/21,1,...,1)>0` has a *strictly positive definite* full-`Z6` Gram `P6*>0`. Embed `Z6` into the extended grid twice: `J0` uses Bob indices `0..6`, and `J1` uses indices `6..12`. Their row ranges cover **all 364** feature coordinates. Hence

\[
 P_*:=J_0P_6^*J_0^T+J_1P_6^*J_1^T\succ0,
 \qquad \mathsf T c_* =\mathsf A P_*.
\]

To see closedness in (18), take a convergent sequence `b_j=T^T y_j+r_j` with `A*y_j>=0,r_j>=0`. Pair with `c*`:

\[
 b_j\cdot c_* =\langle A^*y_j,P_*\rangle+r_j\cdot c_*.
\]

Both nonnegative quantities are bounded. Strict positivity of `P*` bounds `A*y_j` and strict positivity of `c*` bounds `r_j`; therefore `T^Ty_j=b_j-r_j` is bounded. Subtracting the joint kernel of the fixed linear map `y -> (A*y,T^Ty)` gives bounded representatives `y_j`. A subsequence converges to a feasible `y,r`, proving that the right side of (18) is **closed** and hence equals the full dual by bipolarity.

For an actual Bell behaviour `q`, the target **adjoint** is the nine-component vector

\[
 \boxed{b(q)=\bigl(H(q),4p_\lambda(q),4p_n(q),4p_m(q),4p_C(q),
 4p_A(q),4p_D(q),4p_\alpha(q),4p_\beta(q)\bigr),}\tag{19}
\]

because `f_c(q)=c.b(q)` by SG1. Membership `b(q) in K_{6,12}*` is the actual affine constraint `b(q)=T^T y+r` together with `A*y>=0,r>=0`. For each group element, apply the fixed transpose source map to `b(q)`; impose the corresponding SG12/SG13 source-dual pencils (including their source-coordinate adjoints and orthant constraints), and the event inequalities. The variables for different cones remain independent and the physical eight Bell coordinates are shared. Taking one block-diagonal product of these finitely many SDP systems gives a **conditional exact finite rational SDP lift of Q222**. Its only new mathematical assumptions are the oriented strict C1 PSD assertion and the strict C2 PSD assertion at the start of this theorem.

**No new qualitative-existence inference is used**, and (18) explicitly prevents the known error of writing only `A*y>=0` as a Bell behaviour pencil.

## Source and reproducibility boundaries

The exact source identities are the original EN10 reflection and EN13 native inverse in `latex/chapters/05_entry.tex`, and SG1–SG2 in `latex/chapters/04_geometry.tex`. The closure/density result and strict 196-feature rational Slater point are previously independently verified imports, not reproved here. The `k=4` obstruction is used only to emphasize **why the orientation is essential**, not as input to the transport algebra. The general strict C1 and C2 PSD feasibility statements remain **conditional**. The included `verify_reflection_transport.py` performs exact SG2 expansion, all native/contact gap transformations, 9-by-9 Bell coefficient map, 364-grid support with a 32-entry fixed signed transport, exact source normalization at the obstruction witness, and both 196-row embeddings. It makes no solver calls or numerical sweeps. `transport_sparse_maps.json` lists every original indexing and nonzero mapping coefficient.
