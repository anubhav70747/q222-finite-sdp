# Native C2 c-end finite-c continuation: recovered exact frontier

**Recovery only (2026-10-26).** This is a reconstruction of the previous run's saved *results and displayed source/output*. It is not an independent mathematical run. Original runtime binary `.pkl` files and some scripts were not mounted into the current container session; see `MISSING_DEPENDENCIES.txt`. Do not mistake the recovered code for a fully executable standalone checker.

## 1. Target and exact strict parameter domain

Let
\[
\epsilon=c^{-2}>0,\qquad b=a^2t+\eta\epsilon,\qquad r=\rho\epsilon,
\]
with
\[
 a,t,\eta,\rho>0,\qquad 0<\rho<t,\quad 0<t-\rho<1,\quad 0<\eta<1,\quad a^2t+\eta>1.
\]
For each fixed strict face point, native C2 source coefficients are rational in \(\epsilon\), with all original gaps and all six event weights strictly positive for sufficiently small \(\epsilon\). The fixed certificate target is \(Z_6=\{U^{i/6}V^{j/6}\psi_\ell:0\le i,j\le6,\ \ell\in\{W,X,Y,Z\}\}\), without new testing denominators. All limit-face results were obtained in the original pinned exact checker, included here as `check_c2_c_limits.py`.

The exact **complete four-state finite-epsilon parent** is
\[
\begin{aligned}
\widetilde F_\epsilon={}&4n(U^2+\epsilon)(V^2+\epsilon)X^2+4A(U^2+\epsilon)(V^2+\epsilon)Y^2\\
&+4m(VW+UVY-UZ-\epsilon X)^2+(4D+2)(UW+UVX-VZ-\epsilon Y)^2\\
&+(4\alpha+2)(U^2+\epsilon)(VY-Z)^2
+\frac{4C}{\epsilon}(V^2+\epsilon)(UW-\epsilon Y)^2\\
&+2(V^2+\epsilon)(2UXZ+\epsilon X^2-Z^2).
\end{aligned}
\]
The **actual native coefficients** and their exact rational parameter dictionary are in `check_c2_c_limits.py`. At the source contact, the rescaled coordinates are
\[
(U_c,V_c;W_c,X_c,Y_c,Z_c)=(1,a; b/a,-\rho,\eta/a,1).
\]
An invertible pure-source normalization fixes the contact to \((U,V,\psi)=(1,1,\mathbf1)\), using \(V_{\rm old}=aV,\ W_{\rm old}=(b/a)W,\ X_{\rm old}=-\rho X,\ Y_{\rm old}=(\eta/a)Y,\ Z_{\rm old}=Z\). All subsequent algebraic coefficient maps refer to this normalized form.

## 2. Preserved exact low-grade coefficient system

The saved computation used the integer part of the sixth-root grid with fourteen state-linear rows:

- **Old endpoint:** \((UW,UVX,VZ,VW,UZ,UVY)\).
- **First-order charged:** \((UVW,Z,UX,VX,UY,VY)\).
- **Second-order:** \((X,Y)\).

The separate odd sector contains \((\sqrt{UV}W,\sqrt{UV}Z)\). This *subgrid* is not the whole 196-feature \(Z_6\), and its coefficient reachability is only a sufficient route to admission.

For retained rows \(z_i=U^{p_i}V^{q_i}\psi_{s_i}\), the affine coefficient map is
\[
  [U^pV^q\psi_s\psi_t]\;(z^T Gz)
   = \sum_{\substack{i,j:\ p_i+p_j=p,\ q_i+q_j=q\\\{s_i,s_j\}=\{s,t\}}}
         (2-\delta_{ij})G_{ij},
\]
where the sum is over unordered pairs \(i\le j\), and the odd 2-by-2 sector contributes separately. The saved map for the 14 integer rows had **78 distinct polynomial coefficient rows**, **106 symmetric Gram variables including the odd-sector gauge**, and rank **78** before contact constraints. The contact-augmented map had rank **86**, with all source contact and full first-jet equations enforced. A narrower twelve-row version was ruled out because it cannot produce the independent \(U\,WY\) and \(V\,WX\) first-order terms; admitting \(X,Y\) at second grade restores exact first-order affine compatibility.

The generic first-order charge matrix had rank 72 and left nullity **34**. The saved exact source substitution gave zero for every relation indexed 0 through 33, including all tied rate parameters. A second two-kernel extension was a \(106\times106\) rational map of rank **93** and left nullity **13**. All thirteen generic source-compatibility identities vanished exactly. At the rational point below, **all thirteen full finite-epsilon coefficient compatibility relations vanished identically in epsilon**; this means exact affine reachability on the prescribed two-kernel face, not PSD.

## 3. Forced vertex charges from the FULL parent

The unique extreme diagonal coefficients force
\[
\begin{array}{c|c}
\text{feature}&\text{Gram diagonal}\\\hline
UVW&4C/\epsilon\\ Z&4\alpha\epsilon\\ UX&4n\epsilon\\ VX&(4n+2)\epsilon\\
UY&4A\epsilon\\ VY&(4A+4C+4\alpha+2)\epsilon\\
X&(4n+4m+2)\epsilon^2\\ Y&(4A+4C+4D+2)\epsilon^2
\end{array}
\]
and the exact zeros \(G_{UX,UY}=G_{VX,VY}=0\). Three other original full-parent coefficients that must be matched simultaneously are
\[
 [UXZ]\widetilde F=4\epsilon(2m+1),\quad
 [VYZ]\widetilde F=8\epsilon(D-\alpha),\quad
 [UVXY]\widetilde F=-4\epsilon(2D+2m+1).
\]
The all-face first-order charge graph uses six vertices \((P,Z,UX,VX,UY,VY)\), where \(P=UVW\), and the exact 13 allowed edges enumerated in `charge_graph_rational.csv`. Its connected Laplacian has the prescribed six positive degrees, preserves the contact kernel, and omits the two forbidden \(UX-UY\), \(VX-VY\) edges. The generic vertex-degree proof and torus-orbit normalization were already saved in the earlier charge audit but the original checker file is not mounted here.

## 4. Rational face point and exact graded positivity data

At
\[
  (a,t,\eta,\rho)=\left(1,\frac32,\frac14,1\right),
\]
all inequalities above are strict. The exact endpoint coefficients are
\[
 (n_0,m_0,C_0,A_0,D_0,\alpha_0)
 =\left(\frac1{10},\frac8{15},0,\frac{12}{5},\frac3{10},\frac56\right),
\]
with \(C_2>0\) at second source order. The previous run's first-order source derivatives were \((-33/50,-48/25,0,-66/25,3/25,-10/9)\). A previously chosen 13-edge strictly positive rational charge graph is saved as `charge_graph_rational.csv` (the edge numbers are literal exact outputs from `charge_good_rational.py`). Its six charged rows have rank five after imposing the contact kernel.

The endpoint 6-by-6 Gram has rank five and kernel exactly the contact. In the natural five-dimensional complement, leading principal minors of the recovered rational \(G_0\) restriction were
\[
 \frac{36}{5},\ \frac{72}{25},\ \frac{175392}{78125},\
 \frac{3701376}{390625},\ \frac{677376}{78125}>0.
\]
The order-two projected corrected matrix on the old-contact, new-contact, X, Y directions was recorded as
\[
 H_2=\begin{pmatrix}
16/3&0&-68/15&-4/5\\
0&0&0&0\\
-68/15&0&68/15&0\\
-4/5&0&0&4/5
\end{pmatrix}.
\]
After varying the allowed thirteen charge edges, the saved first-order old/second-state Schur costs were
\[
 E_{XX}=\frac{329207}{131250},\qquad
 E_{YY}=\frac{501079}{2362500},\qquad
 E_{XY}=-\frac{2469}{21875}.
\]
Hence the exact recovered residual two-by-two Schur block is
\[
 S_2=\begin{pmatrix}
\frac{265793}{131250}&\frac{2469}{21875}\\
\frac{2469}{21875}&\frac{1388921}{2362500}
\end{pmatrix},\qquad
\det S_2=\frac{2921722373}{2480625000}>0.
\]
These are strict exact rational margins. They are a second-grade condition for this chosen charge allocation, not a certificate for the entire higher-degree Gram cone. They establish local positivity of the displayed graded block at the rational point, but **do not prove** global or compact-uniform finite-\(c\) admission.

## 5. Symmetry and remaining step

The exact native C2 symmetry is
\[
 u^2v^2F_{\sigma c}(1/v,1/u;J\psi)=F_c(u,v;\psi),\quad
 J(W,X,Y,Z)=(Z,X,Y,W).
\]
Its source parameter map is
\[
 c'=\frac1{ac},\quad a'=a,\quad b'=\frac{a^2}{b},\quad
 r'=\frac{a^2c^2r}{b},\quad t'=\frac{1-c^2(b-a^2t)}{b}.
\]
It transforms the open \(c\to\infty\) face onto the open \(c\to0\) face, and sends every fixed sixth-root angular index \((i,j)\) to \((6-j,6-i)\). Hence proving the infinity end would settle the zero end without enlarging \(Z_6\).

**First exact unresolved mathematical implication:** Construct, with the *full 196-feature coefficient map* or with a rigorously sufficient graded subgrid, a rational \(G(\epsilon)\succeq0\) satisfying all finite-\(\epsilon\) coefficient equations and every required contact kernel. The saved affine compatibility checks permit a rational right inverse on the selected two-kernel face; to prove a finite-\(\epsilon\) PSD theorem one still needs a verified complete graded Schur expansion, positive bounds on all surviving higher-order eigenvalue/pivots, and uniform control over every strict open-face parameter (or an explicitly stated smaller neighborhood). No such complete proof or full certificate file was delivered by the terminated run.

## 6. Source and script provenance

- User-uploaded pinned full parent and both c-end exact rational source/endpoint verifier: `check_c2_c_limits.py` copied **verbatim** from the uploaded `Pasted text(20261026-082840).txt` (not rerun here).
- Earlier charge skeleton checker and 13-kB audit: historically saved as `verify_charge_identities.py`, `C2_C_End_Charge_Audit.md`, `exact_check_output.txt`, in `/mnt/data/c2_c_continuation/`; their old runtime-backed bytes are **not currently mounted**.
- Newer previous-run checkers historically saved as `first_order_rank_probe.py`, `second_grade_rational.py`, `charge_edge_family.py`, `charge_good_rational.py`, `graded_schur_exact.py`, `construct_two_kernel_grade2.py`, `two_kernel_rank.py`, `generic_first_order_charge_verify.py`, `generic_two_kernel_verify.py`, `exact_two_kernel_compat.py`, `new_contact_exact.py`, `new_contact_formal.py`, and others. Their exact dependency graph is listed separately.
- Previous-run binary intermediates `first_order_matrices.pkl`, `second_grade_data.pkl`, `charge_edge_family.pkl`, `graded_schur_data.pkl`, `two_kernel_grade2.pkl`, and `exact_full_coeffs.pkl` are **not mounted here**, so the old checker chain is not directly replayable without recovering or regenerating these inputs.

This reconstructed record uses only prior source material and saved outputs, without new calculations or optimization.
