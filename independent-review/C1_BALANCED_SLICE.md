# An exact C1 source slice and a quantified admission tube

27 October 2026. The source cone below has an exact rational SDP representation of total PSD order five. Its entire boundary curve also has a verified global Gram in the fixed Z18 library. That Gram gives an explicit five-native-parameter neighborhood of strict C1 admission. These are scoped global constructions. They do not prove H1 on the complete native domain or prove H2. All new exact checks here are same-assistant checks; independent review is not claimed.

## 1. Exact source cone

Use the canonical source order (d,lambda,n,m,C,A,D,alpha,beta), the complete four-state SG2 parent F, and the physical source functional dH+4 sum_e a_e p_e. Consider

\[
s(m,A,\alpha)=(4m,0,0,m,m,A,A-m,\alpha,\alpha-m).
\tag{1}
\]

The canonical source orthant on this slice is exactly

\[m\ge0,\qquad A\ge m,\qquad\alpha\ge m.\tag{2}\]

**Theorem 1.** Within (2), the source (1) is quantum-valid if and only if

\[A\alpha\ge2m^2.\tag{3}\]

Consequently the cone is represented exactly by

\[
\begin{pmatrix}A&2m\\2m&2\alpha\end{pmatrix}\succeq0,
\qquad m\ge0,\ A-m\ge0,\ \alpha-m\ge0.
\tag{4}
\]

The PSD block orders are (2,1,1,1), total five. No minimality claim is made.

### Necessity and a physical separating strategy

At the orthogonal measurement pair u=v=1, restrict the state to

\[\psi=(z_0,-z_1,z_1,z_0).
\]

Direct substitution in the full parent, before any state restriction is used for sufficiency, gives

\[
F_s(1,1;z_0,-z_1,z_1,z_0)
=16\begin{pmatrix}z_0&z_1\end{pmatrix}
\begin{pmatrix}
\alpha&-\alpha-2m\\
-\alpha-2m&2A+\alpha+4m
\end{pmatrix}
\binom{z_0}{z_1}.
\tag{5}
\]

The determinant of this 2-by-2 matrix is 2(A alpha-2m²). When m>0, (2) makes its lower-right entry positive. If (3) fails, the explicit choice

\[z_0=1,\qquad z_1=\frac{\alpha+2m}{2A+\alpha+4m}
\]

gives cleared slack 32(A alpha-2m²)/(2A+alpha+4m)<0. After normalization this is a genuine two-qubit witness to invalidity. Every state of this displayed form is maximally entangled: its amplitude matrix times its transpose is (z0²+z1²)I. Thus orthogonal measurement pairs and maximally entangled states detect every failure of (3) inside the canonical orthant. For m=0, (3) follows directly from A,alpha>=0.

This is a separating strategy, not a claim about the maximum violation of every invalid source. The ambient scenario remains Q222 with two binary choices per party and all eight marginal/correlator coordinates; the established two-qubit reduction licenses the dimension-unrestricted validity conclusion.

### Sufficiency and the entire boundary curve

For 1/3<=r<=1/2 define

\[
c_r=\left(8,0,0,2,2,
\frac{2(1-r)}r,\frac{2(1-2r)}r,
\frac{4r}{1-r},\frac{2(3r-1)}{1-r}\right).
\tag{6}
\]

For the open interval, this is the actual strict native C1 source with

\[(k,f,g,a,t)=(1,1,1,1-r,1-r),
\quad (u_c,v_c;\psi_c)=(1,1;1,-r,r,1).
\tag{7}
\]

All native signs hold on the whole interval; the checker verifies them after the positive rational parametrization r=(3z+2)/(6(z+1)), z>0. SG7 therefore gives global physical validity. The endpoints follow by coefficient-space continuity, and the direct fixed-feature certificate below also includes them.

Suppose m>0 and (2)-(3) hold. Put

\[
A_* =\min(A,2m),\quad
\alpha_* =2m^2/A_*,\quad r=m/(A_*+m).
\]

Then m<=A_*<=2m, m<=alpha_*<=alpha and 1/3<=r<=1/2. With e_(A+D) denoting one unit in each of A,D, and e_(alpha+beta) one unit in each of alpha,beta,

\[
s(m,A,\alpha)=\frac m2c_r
 +(A-A_*)e_{A+D}+(\alpha-\alpha_*)e_{\alpha+\beta}.
\tag{8}
\]

Both additional terms are nonnegative sums of event sources. Equation (8) proves sufficiency. If m=0, the source is exactly A e_(A+D)+alpha e_(alpha+beta), so no division by m is needed.

This proof uses the established native stationary-source theorem under its actual hypotheses. It does not infer global validity from the positive matrix at one testing angle. The necessity calculation explains why that particular test already captures all quantum restrictions on this slice.

The full-parent identity in Section3 supplies a direct algebraic proof of the same restricted sufficiency statement. Once that fixed rational certificate is checked, no source-wide Gram-admission hypothesis is needed for this slice.

## 2. Evaluated behavior dual

Write b=(b_m,b_A,b_alpha)=S^T E^T(1,q), where S is the explicit source injection in (1), E is the canonical source-to-Bell map and q has all eight physical coordinates. The [evaluated matrix data](checks/c1-balanced-slice/source_and_behaviour_pencil.json) give the complete 3-by-9 target.

The behavior constraint corresponding to every source in this slice is exactly

\[
\begin{split}
b_m&=4Y_{01}+z_m-z_D-z_\beta,\\
b_A&=Y_{00}+z_D,\\
b_\alpha&=2Y_{11}+z_\beta,\\
Y&\succeq0,\qquad z_m,z_D,z_\beta\ge0.
\end{split}\tag{9}
\]

The source point (m,A,alpha)=(1,2,2) is strict for every block of (4), so the displayed source dual is attained. It has the same total PSD order five and three affine equalities. It imposes exactly the inequalities of this source slice; it is not an exact representation of full Q222.

## 3. A fixed global Z18 certificate on the whole interval

The basis is the existing C1 Z18 basis. Its four parity blocks are

\[
\begin{split}
z_A={}&(uW,vW,X,uvX,Y,uvY,uZ,vZ),\\
z_B={}&(uvW,vX,uY,vY,Z,uvZ),\\
z_W={}&\sqrt{uv}(W,Z),\qquad
z_X=\sqrt{uv}(X,Y).
\end{split}
\]

The following definitions use original state coordinates, not contact-normalized coordinates:

\[
\begin{split}
h={}&(vW-uW,\ uZ-uW,\ X+Y,\ uv(X+Y),\ X-uvX),\\
g={}&(uvW-Z,\ uY-vY,\ vX+uY),\\
h_W={}&\sqrt{uv}(W-Z),\qquad h_X=\sqrt{uv}(X+Y).
\end{split}
\]

Define eight root columns

\[
\begin{array}{ll}
R_A=(h,\ X+ruW,\ (r-1/2)(vZ-uW),\ rh),& |R_A|=12,\\
R_B=(g,\ uY-rZ,\ (r-1/3)(uvZ-Z),\ rg),& |R_B|=8,\\
R_W=(h_W,rh_W),\quad R_X=(h_X,rh_X),& |R_W|=|R_X|=2,\\
T_A=(h,vZ-uW),\quad T_B=(g,uvZ-Z),& |T_A|=6,\ |T_B|=4,\\
T_W=(h_W),\quad T_X=(h_X).&
\end{array}
\]

There are eight fixed rational matrices Q_A,Q_B,Q_W,Q_X,P_A,P_B,P_W,P_X, supplied completely in [certificate.json](checks/c1-balanced-slice/certificate.json), of these respective orders. Every matrix is strictly greater than I/4. With omega=(r-1/3)(1/2-r), the complete polynomial identity is

\[
\boxed{r(1-r)F_{c_r}
=\sum_{b=A,B,W,X}R_b^TQ_bR_b
 +\omega\sum_{b=A,B,W,X}T_b^TP_bT_b.}
\tag{10}
\]

It holds for all independent testing u,v and all four state coordinates, interpreted through x=sqrt(u), y=sqrt(v). The checker reconstructs every root, expands every coefficient and verifies the identity in Q[r,x,y,W,X,Y,Z]. It proves all eight bounds Q_b,P_b>I/4 by exact rational LDL of the shifted matrices. Source endpoints cause no inversion of D or beta: omega simply vanishes there. Dividing by the positive pure-source number r(1-r) is legitimate throughout the closed interval.

The rational matrices were found by numerical discovery followed by exact rational correction. The proof consists of the delivered exact identity and rational PSD certificates, with no assumption about the solver's accuracy or optimality. Reproduction needs SymPy only.

The event-pair sources in (8) also have polynomial squares in Z18. Thus the entire order-five source cone, not only its tight boundary, is admitted by this same fixed feature library. The usual SG2 sign congruences and homogeneous closure give all physical testing quadrants and projective endpoints.

## 4. Strictness and a uniform five-parameter tube

The interval certificate has more content than pointwise membership. Normalize the contact by

\[u=U,\quad v=V,\quad\psi=(W,-rX,rY,Z).
\]

The first seven roots of R_A, after deleting the vZ anchor row, have coefficient determinant 1/2-r. The first five roots of R_B, after deleting the uvZ anchor, have coefficient determinant r-1/3. Their inverse entries belong to

\[\{0,\pm1,\ \pm1/(r-r_*),\ \pm r/(r-r_*)\},
\]

where r_* is the corresponding endpoint. This is checked exactly. For

\[0<\eta<1/12,\qquad 1/3+\eta\le r\le1/2-\eta,
\]

every inverse entry has magnitude at most 1/eta. Using Q_A,Q_B>=I/4, r>=1/3 and r(1-r)<=1/4 gives lower bounds on the two contact-reduced integer Grams:

\[H_7^0(r)\succeq\eta^2I/441,\qquad
H_5^0(r)\succeq\eta^2I/225.
\tag{11}
\]

For example, ||B_7^-1||_2<=7/eta, the state normalization has smallest absolute diagonal at least 1/3, and 1/[r(1-r)]>=4; their product with 1/4 gives the first bound. The two half-power charges satisfy gamma>=1 and delta>=1/9.

For an arbitrary actual strict native C1 source theta, use its actual finite nonzero contact to form

\[\widehat F_\theta(U,V;\Phi)
=F_{c(\theta)}(u_cU,v_cV;\operatorname{diag}(\psi_c)\Phi).
\]

Let b(theta) be its 42-key coefficient vector, ordered as in [contact_map.json](checks/c1-balanced-slice/contact_map.json). The fixed C1 map M has shape 42-by-45 and rank36. Its image is exactly the kernel of the six contact derivatives. A fixed rational right inverse E on that kernel changes only the two integer blocks, with entrywise row-sum bounds 40 and17/2. These are evaluated, not unspecified generic maps.

Let x0(r) be the exact interval Gram vector in the normalized contact frame. Set

\[x(\theta)=x_0(r)+E\,[b(\theta)-b(\theta_0(r))],
\qquad\theta_0(r)=(1,1,1,1-r,1-r).
\tag{12}
\]

All coefficients match because actual native stationarity supplies the six compatibility equations. If

\[
\boxed{\|b(\theta)-b(\theta_0(r))\|_\infty
<\eta^2/246960,}
\tag{13}
\]

the spectral perturbations of H7 and H5 are at most 280 times that distance. Equations (11)-(13) preserve positive definiteness, and E leaves gamma,delta unchanged. Thus every such source has an exact global Z18 Gram. The formula is rational in the source parameters and r, with pure-source denominators only.

This is an explicit relative open region in all five native parameters, not only the one-dimensional curve (7). For example the checker takes r=2/5, eta=1/24 and a nonzero rational epsilon=2^-40, perturbs

\[(k,f,g,a,t)=(1-\epsilon,1+\epsilon,1-\epsilon,3/5+\epsilon,3/5-\epsilon),
\]

and verifies every native sign, the exact bound (13), all coefficient equations and rational PSD. The theorem itself covers every actual native point satisfying (13), without a sampling argument.

## 5. Scope and reproduction

Run

```bash
python independent-review/checks/c1-balanced-slice/verify_c1_slice.py
```

The checker writes the [source/behavior pencil](checks/c1-balanced-slice/source_and_behaviour_pencil.json) and the [exact result file](checks/c1-balanced-slice/exact_results.json), including tube constants and the off-family Gram. The copied native dictionary is the unchanged original SD1 dictionary; no revised free-weight surrogate is used.

The order-five slice is a complete restricted representation with both directions proved. The tube is a new sufficient global all-testing admission region. It does not cover the complete oriented native domain, arbitrary coupled source limits or the full K364 H1 claim. The hyperbolic boundary has a maximally entangled noncommuting contact, but these Bell sources also have deterministic zeros; no self-testing claim follows.
