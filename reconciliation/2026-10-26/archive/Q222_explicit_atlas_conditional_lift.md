# A finite affine source atlas: explicit conditional Q222 lift

**Research note, 24 October 2026.** This note accepts the finite-spectrahedral-shadow theorem and the exact source calculations in [anubhav70747/q222-finite-sdp](https://github.com/anubhav70747/q222-finite-sdp). It constructs explicit *sound source sub-cones*, supplies exact block-order bounds for their dual outer approximations, and identifies a necessary-and-sufficient *extreme-source coverage test* for the atlas to equal the entire canonical source cone. The coverage test is **not yet proved**, so these are not asserted to be exact lifts of the full quantum set.

## 1. Source cones and a 12-face atlas

Let \(C\subset\mathbb R_+^9\) be the canonical Bell-validity cone in the fixed source order
\(c=(d,\lambda,n,m,C,A,D,\alpha,\beta)\). Let \(R\subsetneq C\) be the repository's exact rational order-49 *ramified quartic completion* cone, with 148 auxiliary coordinates (see `extracted-lifts/canonical-completion/ramified-completion/data/ramified_order49_pencil.json`). The source-label involution

\[
 \pi(d,\lambda,n,m,C,A,D,\alpha,\beta)
 =(d,C,A,m,\lambda,n,D,\beta,\alpha)
\]

preserves \(C\), and the repository proves that \(R\) and \(\pi R\) are incomparable. Both are sound cones. The strict source \(\bar c\in R\), excluded image \(\pi\bar c\notin R\), and rational dual separator are in `canonical-completion/ramified-obstruction/`.

Write the cycle of event labels as \(\Omega=(\lambda,D,C,A,\alpha,m,\beta,n)\), with indices modulo eight. SG8–SG13 classify all noncoercive sources as lying in exactly the following twelve *possible* source faces:

\[
 F_{i,3}=C\cap\{c_{\Omega_i}=c_{\Omega_{i+3}}=0\}
 \quad(0\le i<8),\qquad
 F_{i,4}=C\cap\{c_{\Omega_i}=c_{\Omega_{i+4}}=0\}
 \quad(0\le i<4).
\]

Every distance-three face is a canonical relabeling of \(\lambda=m=0\), with explicit degree-two matrix-positive Gram \(\begin{psmallmatrix}K_0&K_1/2+S\\K_1/2-S&K_2\end{psmallmatrix}\succeq0\) of order six, \(S\in\mathfrak{so}(3)\) (SG12). Every distance-four face is a relabeling of \(\lambda=\alpha=0\), with two positive \(3\times3\) matrices depending on a scalar \(\xi\) (SG13). In each face, only **seven source coordinates are free**; adding seven diagonal nonnegative-scalar blocks gives an order-13 pencil. The two omitted source coordinates are linear equalities, not zero-size Slater blocks. These relative-face pencils have strict feasible points in their seven-dimensional spans (SG12/SG13), and hence admit attained SDP dual descriptions.

Define the **two-chart plus faces** cone

\[
 K_2=R+\pi R+\sum_{i=0}^7F_{i,3}+\sum_{i=0}^3F_{i,4}\subseteq C.
 \tag{K2}
\]

Introduce separate source vectors \(c^{(0)},\dots,c^{(13)}\) and impose \(c=\sum_jc^{(j)}\), with an order-49 source pencil for each of the first two components and the twelve order-13 face pencils for the others. This is a fully specified *construction algorithm* for a rational affine PSD lift of \(K_2\), of total block order

\[
 \boxed{2\cdot49+12\cdot13=254}.\tag{S254}
\]

The matrices required are the repository's published order-49 JSON pencil, SG12 and SG13, and the eight-event dihedral permutations; no coefficient fitting or degree choice is needed. The order count does not count scalar linear equalities or auxiliary variables, but does count every actual PSD block.

## 2. A symmetry-invariant alternative

Let \(\Gamma\cong D_8\) be the sixteen-element CHSH-preserving relabeling group acting on the cycle \(\Omega\), with \(d\) fixed. Set

\[
 K_{16}=\sum_{\gamma\in\Gamma}\gamma R
 +\sum_{i=0}^7F_{i,3}+\sum_{i=0}^3F_{i,4}\subseteq C.
 \tag{K16}
\]

This is \(\Gamma\)-invariant by construction, and has an explicit rational source PSD lift of total block order at most

\[
 \boxed{16\cdot49+12\cdot13=940}.\tag{S940}
\]

Some ramified images could coincide, in which case a smaller order suffices; no such identification is assumed here. In contrast, \(K_2\) requires no symmetry closure and is more economical.

The cones are closed. Each summand is closed (for \(R\), see the repo's positive-definite Gram-moment integral in `proof/ramified_regularity.tex`; for the two faces, SG11–SG13), and each is contained in the nonnegative orthant. If \(c^{(n)}=\sum_jc_j^{(n)}\to c\), then coordinatewise \(0\le c_j^{(n)}\le c^{(n)}\). All component sequences are bounded, and a simultaneous convergent subsequence proves that their limits lie in the corresponding closed cones and sum to \(c\). The relabeling images inherit the same property.

## 3. Explicit outer Q222 SDP from any such atlas

The repository constructs eight rational invertible source-to-Bell-functional maps \(M_g\), indexed by the eight CHSH/PR orientations, and establishes the exact full polar assembly \(Q_{222}^{\vee}=\sum_g M_g C\), including deterministic and PR sectors. Given \(K\in\{K_2,K_{16}\}\), form

\[
 Q_K=\left\{p\ \middle|\ M_g^Tp\in K^*
 \text{ for all eight }g\right\}\cap\{\text{normalization}\}.
 \tag{QK}
\]

Since \(K\subseteq C\), this is a **certified outer approximation**: \(Q_{222}\subseteq Q_K\subseteq\mathrm{NS}_{222}\). Here the repository's *no-enlargement* theorem (`canonical-completion/full-polar/full_polar_assembly.tex`) applies because \(K\) is closed and includes the complete \(d=0\) event cone (already contained in \(R\)). It gives the stronger **if-and-only-if** conclusion
\[
 \boxed{Q_K=Q_{222}\quad\Longleftrightarrow\quad K=C}.\tag{NE}
\]
Thus there is no loophole in which the eight orientations make a deficient canonical source cone accidentally exact.

Here is an evaluated-matrix prescription for \(Q_K\). For each source summand \(K_j\), let its homogeneous published pencil be
\(\Lambda_j(x,z)=\sum_{r\in I_j}x_rA_{jr}+\sum_kz_kB_{jk}\succeq0\), where \(I_j\) is the set of its seven or nine active source coordinates. For every \((g,j)\), introduce \(Y_{gj}\succeq0\) of the same block order as \(\Lambda_j\), with linear equations

\[
 \langle Y_{gj},B_{jk}\rangle=0\quad(\forall k),\qquad
 (M_g^Tp)_r=\langle Y_{gj},A_{jr}\rangle\quad(\forall r\in I_j).
 \tag{D}
\]

On a source face, components outside \(I_j\) are unconstrained (the omitted sources are fixed to zero). Slater in each source-face span, and in the ramified cone, makes these attained dual systems exact. Thus \(Q_{K_2}\) is an explicit rational affine outer pencil of order at most **\(8\cdot254=2032\)**, and \(Q_{K_{16}}\) one of order at most **\(8\cdot940=7520\)**. The construction uses matrices already given in the repo, not unspecified positive polynomials. These are *outer*, not proven exact.

## 4. Every genuine C1/C2 contact ray is exposed

Here is a useful converse to the extreme-source reduction in SG10. **Every physical C1 or C2 source with an actual finite contact is itself an exposed extreme ray of \(C\)**, including the literal third-omission boundaries. This follows from the *nonvanishing contact minors* that are already computed in Chapter 5, without invoking the difficult tensor-SOS atlas.

Fix a finite contact \(\tau=(u_c,v_c,\psi_c)\) of \(c\), with \(u_c,v_c>0\), \(Z_c=1\), and define the six linear source equations \(\mathcal L_\tau(c')=0\): the four components of \(P_{c'}(u_c,v_c)\psi_c=0\) and the two angular derivatives \(\psi_c^T\partial_uP_{c'}\psi_c=\psi_c^T\partial_vP_{c'}\psi_c=0\).

For C1, fix the two absent coordinates \(\lambda=n=0\). The selected \(6\times6\) minor on columns \((d,m,A,D,\alpha,\beta)\), with \(C\) the remaining coefficient, is explicitly nonzero by EN8: its factors are the positive native parameter \(\vartheta=1+X_c/u_c\), \(u_c,v_c,D_u,D_v,Y_0,(u_cY_0+v_c),(1-u_cv_cY_0),(S+u_cv_c(\vartheta-1)),R_m,-R_D\), all strictly positive on the C1 contact graph by EN2/EN10. Thus \(\mathcal L_\tau\) restricted to the seven-dimensional active source-coordinate space has rank six and **one-dimensional kernel**, necessarily spanned by \(c\). For C2, fix \(\lambda=\beta=0\). EN16 gives a nonzero selected minor on columns \((n,m,C,A,D,\alpha)\), the remaining coefficient being \(d\); all factors are strictly positive by EN15–EN16, including at \(n=0\). It gives the same one-dimensional-kernel conclusion. The repo verifies strictness also when a third event weight vanishes.

Now suppose \(c'\in C\) satisfies both the same two omitted coordinates and \(F_{c'}(u_c,v_c;\psi_c)=0\). Because \(c'\) is quantum-valid, \(F_{c'}\ge0\) on all state and angle tests. Its zero at an interior test forces all first state and angle derivatives to vanish, i.e. \(\mathcal L_\tau(c')=0\). By the rank-six result, \(c'\) lies on \(\mathbb R_+c\). Consequently the source-dual functional

\[
 \ell_c(c')=F_{c'}(u_c,v_c;\psi_c)+
 \begin{cases}c'_{\lambda}+c'_n,&\text{C1},\\
 c'_{\lambda}+c'_{\beta},&\text{C2},\end{cases}
\]

is nonnegative on all of \(C\) and vanishes **only** on the ray \(\mathbb R_+c\). This explicitly exposes that source ray. In physical terms, a single saturating two-qubit test plus two nonnegative omitted-event coordinates identifies its Bell inequality uniquely up to scale.

**Direct reconstruction from one experiment.** Order the six jet rows as \((F,\partial_WF,\partial_XF,\partial_YF,\partial_uF,\partial_vF)\). Let \(J_\tau\) be the \(6\times9\) matrix whose column for each source label is the corresponding jet of the elementary SG2 packet. For C1 let \(A_1(\tau)\) consist of columns \((d,m,A,D,\alpha,\beta)\), and let \(b_1(\tau)\) be the column for \(C\). For C2 take columns \((n,m,C,A,D,\alpha)\) and the remaining column \(b_2(\tau)\) for \(d\). Then, with normalizations \(C=2\) and \(d=1\), respectively,

\[
\begin{aligned}
(d,m,A,D,\alpha,\beta)^T&=-2A_1(\tau)^{-1}b_1(\tau),\\
(n,m,C,A,D,\alpha)^T&=-A_2(\tau)^{-1}b_2(\tau).
\end{aligned}
\tag{REC}
\]

All entries are rational functions of the **five real contact parameters** \((u_c,v_c,W_c,X_c,Y_c)\), and the two determinants are nonzero throughout their physical domains. Thus a measured contact (in an agreed gauge) determines the corresponding extremal Bell-validity inequality up to an overall positive scale. The exact-rational independent seed check in `verify_q222_exposed_contacts.py` gives \(\det A_1=69120\) and \(\det A_2=184320\) in the displayed jet convention; all six stationary equalities hold identically at those seeds. These seed checks do not replace the universal sign proofs EN8 and EN16.

**A nonlocal exposing behaviour in the full eight-coordinate quantum set.** The source-ray exposure above has a direct operational form. Denote by \(q_\tau\in\mathcal Q_{222}\) the normalized correlation from the contact experiment, and by \(v_i,v_j\) the deterministic CHSH-saturating neighbours associated with its two omitted losing-event weights. The polar source map satisfies \(f_{c'}(v_k)=4c'_k\), and SG2 is \(F_{c'}(u_c,v_c;\psi_c)=w\,f_{c'}(q_\tau)\), where \(w=(1+u_c^2)(1+v_c^2)\|\psi_c\|^2>0\). Therefore the explicit mixture
\[
\boxed{
q_{\mathrm{exp}}=\frac{wq_\tau+\tfrac14v_i+\tfrac14v_j}{w+\tfrac12}
\ \in\mathcal Q_{222}
}
\tag{EXP}
\]
has the property that \(f_{c'}(q_{\mathrm{exp}})=\ell_c(c')/(w+\tfrac12)\). It exposes \(f_c\) uniquely **within the canonical source cone** by the rank-six argument. It is also strictly CHSH-nonlocal: each deterministic neighbour has canonical CHSH value two, while \(q_\tau\) has CHSH value strictly larger than two because its source contact satisfies \(f_c(q_\tau)=0\) with \(d>0\), positive event coefficients and strictly positive event probabilities.

In fact, \(q_{\mathrm{exp}}\) exposes the ray \(\mathbb R_+f_c\) in the **entire Bell-validity polar**, not only in the canonical sector. Here is the additional global argument. EN2 makes each of the eight canonical losing-event amplitudes
\[
W+vX,\ X,\ R_m,\ uW-Y,\ Y,\ R_D,\ vY-1,\ X+u
\]
strictly nonzero at the C1/C2 contact. Thus all eight canonical losing probabilities of \(q_\tau\) are strictly positive. Every *winning* probability is also strictly positive at any nonsignalling behaviour with canonical CHSH value \(S>2\): in an extremal decomposition of the nonsignalling polytope, the canonical PR vertex must have weight \(t\ge(S-2)/2>0\), since every other vertex has CHSH value at most two. That PR vertex assigns probability \(1/2\) to each winning outcome, so \(p_{\mathrm{win}}\ge(S-2)/4>0\). Therefore \(q_\tau\), and hence \(q_{\mathrm{exp}}\), have **full sixteen-event support**. Every other CHSH orientation is strictly subcritical: its sign pattern differs from the violated one on two or four edges, so either \(S+S_g\le4\) (the two shared terms each contribute at most one) or \(S_g=-S\); in both cases \(S_g<2\).

If a full quantum-valid affine functional \(f\) vanishes at \(q_{\mathrm{exp}}\), it vanishes separately on the three positively weighted quantum behaviours \(q_\tau,v_i,v_j\). By the full-polar theorem, \(f\) is either a nonnegative event combination or lies in one of the eight oriented source sectors. The event combination cannot vanish on the full-support \(q_\tau\) unless it is zero. A sector other than the canonical one cannot vanish there because its CHSH term has strictly positive slack and all its event terms are nonnegative. Hence \(f=f_{c'}\) for a canonical valid source \(c'\); the two deterministic equalities force \(c'_i=c'_j=0\), and the contact-minor argument forces \(c'\in\mathbb R_+c\). Consequently
\[
\boxed{\{f\in Q_{222}^{\vee}:f(q_{\mathrm{exp}})=0\}=\mathbb R_+f_c.}
\tag{FULL-EXPOSE}
\]
Thus every coercive C1/C2 contact ray is **exposed in the full Bell-functional polar by a concrete nonlocal quantum correlation**, which has a unique supporting Bell inequality up to positive scale. This is a local smoothness statement at the exposing correlation, not a claim of self-testing or finite-size cryptographic security.

For the repository's two exact rational seeds, the exposing mixtures and their CHSH values are especially simple. The C1 seed \((u,v;\psi)=(1,1;1,-1/2,1/4,1)\) has \(w=37/4\), contact CHSH value \(103/37\), and
\[
q_{\mathrm{exp}}^{\mathrm{C1}}=\frac{37}{39}q_\tau
 +\frac1{39}v_\lambda+\frac1{39}v_n,
\qquad S(q_{\mathrm{exp}}^{\mathrm{C1}})=\frac{107}{39}>2.
\]
The C2 seed \((u,v;\psi)=(1,1;1,-1/4,1/4,1)\) has \(w=17/2\), contact CHSH value \(46/17\), and
\[
q_{\mathrm{exp}}^{\mathrm{C2}}=\frac{17}{18}q_\tau
 +\frac1{36}v_\lambda+\frac1{36}v_\beta,
\qquad S(q_{\mathrm{exp}}^{\mathrm{C2}})=\frac83>2.
\]
Here the deterministic neighbours are completely elementary: \(v_\lambda\) uses outcome assignments \(A=(+,+)\), \(B=(+,+)\); \(v_n\) uses \(A=(+,-)\), \(B=(-,-)\); and \(v_\beta\) uses \(A=(-,+)\), \(B=(-,+)\). The contact CHSH values follow directly from the signed-core coefficient \(F_d(\tau)/w=H(q_\tau)\) in SG2. All these ratios were evaluated by exact fractional arithmetic.

The exposed face in the primal quantum set has a simple operational geometry. It contains \(q_\tau\) and the two *distinct deterministic* neighbours \(v_i,v_j\) that saturate the omitted source events. Both neighbours have CHSH value two, whereas \(q_\tau\) violates CHSH, so these three behaviours are affinely independent. Thus **every coercive extreme canonical Bell functional supports an exposed quantum face of dimension at least two, containing a nonlocal point and two local vertices**. If a third event weight is omitted, the face contains a third distinct deterministic vertex as well, so its dimension is at least three: three distinct deterministic assignments are affinely independent, all lie in the CHSH-equals-two hyperplane, and the nonlocal contact point lies outside that hyperplane. The dimension bound is a lower bound, not an assertion that the full exposed face has exactly that dimension.

**Quantitative local witness stability.** Normalize the source by \(\|c\|_1=1\). At its fixed contact \(\tau\), there is an explicit finite constant \(C_\tau\) such that every canonical valid source \(c'\ge0\), also normalized by \(\|c'\|_1=1\), satisfies
\[
\|c'-c\|_1\ \le\ C_\tau\,\sqrt{\ell_c(c')}
\qquad\text{whenever }0\le\ell_c(c')\le\varepsilon_\tau.
\tag{ROB}
\]
To see this without an appeal to compactness, write \(e=\ell_c(c')\). The omitted source weights sum to at most \(e\), and \(0\le F_{c'}(\tau)\le e\). Let \(M_\tau=\max_j\|P_{e_j}(u_c,v_c)\|_{\mathrm{op}}\), where \(e_j\) are the nine unit source directions. Since \(P_{c'}(\tau)\succeq0\), its state gradient satisfies \(\|P_{c'}(\tau)\psi_c\|\le\sqrt{M_\tau e}\). Fix a small closed interval around each positive contact angle, and let \(K_u,K_v\) bound there the absolute second angular derivatives of the nine single-source packets at the fixed state \(\psi_c\). Positivity of \(F_{c'}\) throughout these intervals and the two-sided Taylor inequality imply \(|\partial_u F_{c'}(\tau)|\le\sqrt{2K_u e}\), \(|\partial_v F_{c'}(\tau)|\le\sqrt{2K_v e}\) for sufficiently small \(e\). Thus the six-vector of jets \(J_\tau c'\) has Euclidean norm at most \(\kappa_\tau\sqrt e\), where one may take \(\kappa_\tau=\sqrt{1+12M_\tau+2K_u+2K_v}\) for \(e\le1\).

Now eliminate the six selected source coefficients using (REC). With \(A_j(\tau)\) the corresponding invertible selected minor and \(B_\tau\) the largest norm of the two omitted source-jet columns, this gives
\[
\left\|c'-\frac{c'_{\mathrm{scale}}}{c_{\mathrm{scale}}}c\right\|_1
\le\left[\sqrt6\|A_j(\tau)^{-1}\|_2(\kappa_\tau+B_\tau)+1\right]\sqrt e.
\]
The common \(\ell_1\) normalization then doubles the right side at most, proving (ROB) with
\[
C_\tau=2\left[\sqrt6\|A_j(\tau)^{-1}\|_2(\kappa_\tau+B_\tau)+1\right].
\]
Here the scale coordinate is \(C\) in C1 and \(d\) in C2. The threshold \(\varepsilon_\tau>0\) is chosen so the Taylor steps remain within the fixed positive-angle intervals. The bound is **local and nonuniform** near singular support boundaries, as the inverse contact minor can become large. It certifies stability of the *Bell witness*, not device self-testing.

## 5. Exact if-and-only-if contact coverage, and a finite algebraic test

Normalize the pointed cone \(C\subset\mathbb R_+^9\) by \(\sum_rc_r=1\). Its normalized base \(B\) is compact and convex. In finite dimensions, every point of \(B\) is a finite convex combination of its extreme points; the repository in fact gives a nine-atom bound. For any finite family of sound source sub-cones \(K_j\), an extreme nonzero \(c\in C\) can belong to their sum only if it belongs to at least one of them: if \(c=\sum_jc_j\), extremality forces each nonzero \(c_j\) to lie on the same ray. Conversely, if every extreme ray lies in their union, their sum contains \(C\) by convexity. Therefore

\[
\boxed{Q_{K_{16}}=Q_{222}\quad\Longleftrightarrow\quad K_{16}=C\quad\Longleftrightarrow\quad
 \text{every coercive extreme source is in }\bigcup_{\gamma\in\Gamma}\gamma R.}
 \tag{E16}
\]

There is no missing noncoercive case: the twelve \(F_{i,k}\) already cover every noncoercive source. The repo's SG10 reduces every coercive extreme up to \(\Gamma\) to exactly two source support families, C1 and C2; the third omission boundaries are included literally. The exposed-ray lemma above supplies the converse: **every physical C1/C2 contact source is extreme**. Hence E16 is exactly the problem of covering the *entire* C1/C2 contact graphs with the symmetry images of \(R\). No separate extremality quantifier is needed.

For the smaller \(K_2\), the identical criterion holds with \(R\cup\pi R\) in place of the orbit union, but it must be checked on *all* symmetry images of C1/C2, since \(K_2\) is not known to be \(\Gamma\)-invariant.

The resulting **necessary-and-sufficient**, concrete, first-order condition for the order-7520 pencil to become exact is

\[
 \forall c,u,v,\psi:\quad
 [\,d=1,\ c\ge0,\ \mathrm{C1}(c)\vee\mathrm{C2}(c),\
 u,v>0,\ Z=1,\ P_c(u,v)\psi=0,\
 \psi^TP_{c,u}(u,v)\psi=\psi^TP_{c,v}(u,v)\psi=0\,]
 \Longrightarrow
 \bigvee_{\gamma\in\Gamma}\exists z_\gamma:\
 \Lambda_R(\gamma^{-1}c,z_\gamma)\succeq0.
 \tag{T}
\]

All source equations, positivity conditions, and entries of the 49-pencil are rational polynomials; PSD is equivalent to finitely many principal-minor inequalities. Thus (T) is a single finite first-order sentence over the real closed field, with no unspecified polynomial degree. **Its truth is currently undetermined.** Finite decidability is not a claim of computational tractability, and the existing theorem does not imply (T).

If (T) holds, (E16) proves \(K_{16}=C\) and equations (QK)–(D) become a completely explicit rational affine order-7520 SDP lift of \(Q_{222}\). If (T) fails, the failing **contact ray is exposed**, hence is outside the sum \(K_{16}\) and gives a genuine counterexample to the atlas. By the no-enlargement theorem, the proposed global order-7520 outer approximation is then **strict**. Thus (T), (E16), and global exactness of (QK) for \(K_{16}\) are all equivalent.

Using (REC), (T) can be reduced further to two explicit universal tests over **five contact parameters each**. For each family, compute \(c_j(\tau)\) by its rational reconstruction, restrict to positive angles, nonzero contact minor, and the literal coefficient signs of EN1, and require membership in at least one of the 16 symmetry-transformed order-49 pencils. The contact lemma SG7 makes every reconstructed nonnegative source physical; EN8/EN16 show every original C1/C2 source appears. Conversely any failure is an exposed missing extreme ray. This removes the source variables, the stationarity quantifier, and the abstract extreme-ray predicate from the coverage test. Clearing denominators by a *strictly positive square* converts it into a polynomial first-order sentence.

## 6. Interpretation and the remaining obstruction

The angular determinant now has an independent exact order-13 PSD description, reducing to order nine and at most six squares on C1/C2 contacts (`Q222_angular_order13_note.md`). Its exceptional feature is a square root of the positive Bob half-angle, arising from the term \(U\sqrt V(1+\sqrt V)^2\) in the elementary balanced-circuit identities. The known ramified order-49 completion uses exactly such Bob-root testing functions; party exchange moves them to Alice, explaining why a single ramification is not symmetry stable.

The **unresolved mathematical statement** is coverage of the genuinely curved C1/C2 extreme boundary by finitely many *affine source* chart cones. Angular certificates are affine in \(h,q_{ij}\), which are quartic in \(c\), so they do not prove E16 or T. Conversely, the new source-atlas candidate is affine in \(c\) by construction, but its reverse inclusion is not established. This is the precise point at which an explicit outer SDP may or may not become an exact lift.