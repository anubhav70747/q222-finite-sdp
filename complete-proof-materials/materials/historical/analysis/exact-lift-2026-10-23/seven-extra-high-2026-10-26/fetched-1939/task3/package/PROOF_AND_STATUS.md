# Direct C2 K364: full polynomial jet repair and exact frozen-top no-go

**Status (27 October 2026).** This note does **not** prove strict C2 admission H2 to `K364` and does **not** disprove `K364`. It proves a new exact necessary jet condition for the entire polynomial 32-row positive-packet normal completion at one genuine strict source, and an exact, stronger obstruction to retaining the earlier `tau=2` **full top-degree packet** while varying all lower polynomial syzygies. The missing freedom is an independently variable degree-`v^2` charge vector or a genuinely different PSD Gram, including fractional features. All tests and all four state components remain independent in the parent identity.

## 1. Physical source and complete parent

Use the native strict C2 contact source

\[
(a,b,c,t,r)=(1,1,1,3/4,1/4),\qquad
(d,\lambda,n,m,C,A,D,\alpha,\beta)
=(1,0,1/8,3/8,1/6,5/8,3/8,1/6,0).
\]

All strict native gaps are positive: `h=1/4`, `g=j=1/2`, `1-c²h=3/4`, `b-a²c²r=3/4`, `N_*=2`, `N_n=5/8`, and `Delta=10`. The genuine stationary contact is `(u,v;W,X,Y,Z)=(1,1;1,-1/4,1/4,1)`.

Put `Du=1+u²`, `Dv=1+v²`. With the sixteen complete positive event/core roots `r_i` and their squared weights `k_i²` explicitly generated in `verify_exact.py`, the *entire four-state parent*, not a contact restriction, is the exact identity

\[
(1+v^2)F=\sum_{i=0}^{31}k_i^2(r_i\cdot\psi)^2
                -2(1+u^2)(1+v^2)^2Z^2,
\quad r_{16+i}=vr_i,\quad k_{16+i}=k_i.
\tag{1}
\]

The 32 rows have angular bidegree at most `(1,2)`, hence are in `Z_(6,12)`. Their first three state columns form `R`, and their fourth column is `q`. For any polynomial syzygy `t` of bidegree at most `(1,2)` satisfying `R^Tt=0`, the *full state* identity is

\[
(1+v^2)F=\sum_{i=0}^{31}\frac{[k_i^2(r_{i,W}W+r_{i,X}X+r_{i,Y}Y)+(k_i^2r_{i,Z}+t_i)Z]^2}{k_i^2}
   +S_t(u,v)Z^2,
\]
\[
S_t=-2D_uD_v^2-2\sum_i r_{i,Z}t_i-\sum_i t_i^2/k_i^2.
\tag{2}
\]

The positive sum in (2) is a genuine fixed-library Gram. Whether `S_t` is scalar SOS is a separate assertion, not a consequence of the identity.

## 2. Complete contact-angular-jet freedom and optimal budget

Write each `t_i=Σ_{p=0}¹Σ_{q=0}² x_{ipq}u^pv^q` and form the exact 45-row polynomial coefficient map `A` for `R^Tt=0`. Let `E` evaluate the 32 `t_i` at the contact, and let `U` take their 32 derivatives in `u`. The exact integer map `[A;E;U]` has rank 99; after contact-zero reduction, its attainable `u`-jet has dimension 26 and **six** independent equations. They are reconstructed directly as the `u` portions of the rational left kernel of `[A;E;U]`.

Take the old contact-interpolated `tau=2` syzygy `t_0`, using the **correct** `Y_u` interpolation correction `eta_9+v_c eta_4`. Its 32 completed-root numerator derivatives `g_u` have weighted energy

\[
\|g_u\|_{K^{-2}}^2=1523/40>289/24.
\]

The original parent on the *fixed physical contact state* and the independent testing line `v=1` is exactly

\[
(1+v^2)F(u,1;1,-1/4,1/4,1)=\frac{289}{24}(u-1)^2.
\tag{3}
\]

For any nonnegative scalar remainder in (2), its `u`-curvature at contact cannot be negative. Let `L` be the exact `6×32` matrix of the attainable variation-jet annihilators, and `K²=diag(k_i²)`. The unique sharp minimum is given by the **pure-source rational projection**

\[
g_u^*=K^2L^T(LK^2L^T)^{-1}Lg_u,
\qquad
\|g_u^*\|_{K^{-2}}^2=\frac{289}{24}.
\tag{4}
\]

Therefore every nonnegative completion of this entire 32-root normal-packet ansatz must use the **same optimum contact `u`-jet**. Equality is permitted. Formula (4) and all six constraints are checked by rational arithmetic; the inverse is of a strictly positive definite **source-only** six-by-six matrix. The completed `u`-jet is attained by a polynomial syzygy, not only by a formal tangent vector.

## 3. Exact no-go for freezing the old complete high-degree allocation

In addition to the necessary contact and optimum-jet conditions above, impose the extra ansatz restriction that every `v²` and `uv²` coefficient of the polynomial syzygy remain equal to those of the old `tau=2` top-degree transfer. This restriction is **not necessary for `K364` membership**. It preserves the earlier positive highest-Bob pivot but fixes too much charge.

All such polynomial corrections form the exact rational affine space

\[
x=x_0+N\zeta,\quad \zeta\in\mathbb R^{38},
\tag{5}
\]

obtained by equating (i) all 45 state-normal polynomial coefficients, (ii) all 32 contact values, (iii) all 32 optimal first `u`-jet entries, and (iv) both 32-vectors of `v²`/`uv²` coefficients. Every equality is required and generated from (1), not postulated from a scalar graph.

Two necessary conditions for *any* scalar nonnegative remainder are its constant state-`Z` value at `(u,v)=(0,0)` and its highest `u²` coefficient along `v=0`, evaluated at the physical contact state:

\[
M_0(\zeta)=\frac23-\|O\zeta+o\|_{K^{-2}}^2\ge0,
\qquad
M_u(\zeta)=\frac{409}{48}-\|B\zeta+b\|_{K^{-2}}^2\ge0.
\tag{6}
\]

`B,O` are the two **explicit** `32×38` selector matrices computed from (5); their offset vectors `b,o` are the original full-parent coefficient vectors. They are deterministically reconstructed from the sixteen roots in `verify_exact.py`. Continuity of an SOS certificate on `u,v>0` makes the literal zero tests legitimate necessary conditions.

The `64×38` stacked matrix `E_*=[B;O]` has full column rank. Exact rational least squares, with `W=diag(K^{-2},K^{-2})`, gives

\[
\begin{aligned}
\inf_{\zeta}\big(\|B\zeta+b\|_{K^{-2}}^2+\|O\zeta+o\|_{K^{-2}}^2\big)
&=d^TWd-d^TWE_*(E_*^TWE_*)^{-1}E_*^TWd\\
&=\boxed{\frac{57818093029}{6288328800}},\qquad d=(b,o).
\end{aligned}
\tag{7}
\]

Since the sum of the capacities is `409/48+2/3=147/16`, every member of this frozen-top packet family satisfies the **exact universal inequality**

\[
\boxed{M_u(\zeta)+M_0(\zeta)
\le-\frac{44072179}{6288328800}<0.}
\tag{8}
\]

Thus even **all** lower-degree bidegree-`(1,2)` syzygy gauges cannot repair the old fixed top-degree packet at this genuine strict source. This is a no-go theorem for the specified *normal assembly*, not for `K364` or for all 32-root assemblies with variable top coefficients.

## 4. Why this does not exclude K364

Releasing the unnecessarily frozen `v²` coefficients, the *same two necessary budgets* have positive maximum combined margins:

\[
\begin{array}{c|c|c}
\text{Allowed polynomial syzygies after optimal contact jet}&\text{minimum combined norm}&\text{max combined margin}\\\hline
\text{All bidegree }(1,2)&1019/144&19/9\\
\text{Only the compelled }uv²\text{ zero retained}&233201/29952&41983/29952
\end{array}
\tag{9}
\]

These are exact dual projections, not certificates of simultaneous PSD at *all* testing angles. They show precisely that (8) is caused by freezing the old full `v²` charge vector. The full `K364` library additionally permits fractional rows and non-packet Gram allocations.

## 5. Remaining actual H2 implication

For **every** genuine strict physical C2 source `θ` in its full native gap domain `Ω`, the unresolved requirement remains

\[
\boxed{\exists P_\theta\succeq0:\quad
\mathrm{coeff}(Z_{(6,12)}^TP_\theta Z_{(6,12)})
=\mathrm{coeff}[(1+v^2)F_{c(\theta)}].}
\tag{H2}
\]

Equivalently, the exact `K364` moment dual must evaluate the complete parent nonnegatively for every admissible source. This note does not prove `(H2)` or produce a dual pseudo-moment separating an actual source from `K364`. The independent qualitative existence proof, strict open-end result, and conditional global SDP assembly are unchanged.

## 6. Reproduction

Run `python verify_exact.py` (SymPy only). It builds the *whole four-state parent* and 32 event rows from scratch, verifies all six contact derivatives, 29 polynomial syzygies and the correct `Y_u` interpolation, computes the full bidegree-`(1,2)` syzygy map, its exact six jet constraints and optimal contact curvature, reconstructs the **entire** frozen-top affine space, and checks (7)–(9) in rational arithmetic. The output `EXACT_RESULTS.json` records the exact rational conclusions. No unavailable `.pkl`, network, GitHub write, numerical sweep, or SDP solver is involved in this final proof replay.
