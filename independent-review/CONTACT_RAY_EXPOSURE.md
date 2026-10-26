# A contact experiment exposes its Bell source ray

This is a deduction from the source classification SG1–SG10 and the nonvanishing contact minors EN1–EN16. It does not use the long tensor-SOS atlas.

Let \(C\subset\mathbb R_+^9\) be the canonical quantum-valid source cone in the order
\((d,\lambda,n,m,C,A,D,\alpha,\beta)\).
Let \(c\in C\) be a physical C1 or C2 source with a finite noncommuting contact
\(\tau=(u_c,v_c,\psi_c)\), gauged by \(u_c,v_c>0\), \(Z_c=1\).
The two mandatory omitted coordinates are \(I=\{\lambda,n\}\) in C1 and \(I=\{\lambda,\beta\}\) in C2. Literal third omissions are allowed.

The six equations consisting of the four components of \(P_{c'}(\tau)\psi_c=0\) and the two angular stationary equations are linear in \(c'\). On the seven-dimensional subspace \(c'_i=0\) for \(i\in I\), they have rank six. For C1 use the selected columns \((d,m,A,D,\alpha,\beta)\) and EN8; for C2 use \((n,m,C,A,D,\alpha)\) and EN16. Their factors are nonzero throughout the stated contact domains, including the permitted third omissions. Hence their kernel is precisely \(\mathbb Rc\).

For \(c'\in C\), define
\[
\ell_c(c')=F_{c'}(\tau)+\sum_{i\in I}c'_i\ge0.
\]
If it vanishes, each summand vanishes. Positivity of the full state matrix gives \(P_{c'}(\tau)\psi_c=0\); a zero at an interior angular test gives both angular derivatives. The selected minor therefore forces \(c'\in\mathbb R_+c\). Thus every such contact ray is exposed in the canonical cone.

## A physical exposing correlation for the full Bell-validity cone

Let \(q_\tau\) be the normalized contact behaviour, let \(v_i,v_j\) be the deterministic CHSH-facet neighbours corresponding to the omitted events, and put
\[
w=(1+u_c^2)(1+v_c^2)\|\psi_c\|^2,\qquad
q_{\mathrm{exp}}=\frac{wq_\tau+\frac14v_i+\frac14v_j}{w+\frac12}.
\]
The identities \(F_{c'}(\tau)=wf_{c'}(q_\tau)\) and \(f_{c'}(v_k)=4c'_k\) give
\[
f_{c'}(q_{\mathrm{exp}})=\frac{\ell_c(c')}{w+\frac12}.
\]

All eight canonical losing probabilities of \(q_\tau\) are strictly positive by EN2. Its CHSH value \(S>2\): the positive retained event contribution at contact must balance \(d(1-S/2)\), with \(d>0\).

Every winning probability is also positive. In a nonsignalling extremal decomposition, the canonical PR weight is at least \((S-2)/2\), because every other vertex has this CHSH value at most two. Its winning probabilities are \(1/2\); hence \(p_{\mathrm{win}}\ge(S-2)/4>0\). Thus \(q_\tau\) has full event support.

For every other CHSH orientation \(g\), either \(S_g=-S\), or \(S+S_g\le4\); therefore \(S_g<2\). A nonzero quantum-valid affine functional vanishing at \(q_{\mathrm{exp}}\) vanishes at each of its three positively weighted constituent behaviours. It cannot be a nonnegative event combination, because \(q_\tau\) has full support, or belong to another PR sector, whose CHSH slack is strictly positive at \(q_\tau\). It belongs to the canonical cone. The two deterministic equalities impose its omitted coordinates, and the contact-minor argument makes it a positive multiple of \(f_c\). Consequently
\[
\{f\in Q_{222}^{\vee}:f(q_{\mathrm{exp}})=0\}=\mathbb R_+f_c.
\]
This proves exposure in the full cone of Bell-valid affine functionals.

The primal exposed face contains \(q_\tau,v_i,v_j\), which are affinely independent: the deterministic points have CHSH value two and the contact has a larger value. Its dimension is at least two. A third omitted event gives a third distinct deterministic vertex, hence dimension at least three. These are lower bounds. The mixture has a unique supporting Bell ray; this is not a self-testing assertion.

## Reconstruction and exact examples

With jet rows \(F,F_W,F_X,F_Y,F_u,F_v\), the fourth state-kernel equation follows by Euler homogeneity at \(Z=1\). Fixing \(C=2\) in C1 or \(d=1\) in C2, the invertible selected minor recovers the six remaining source coefficients rationally from the five contact parameters.

The two repository seeds give:
\[
\begin{array}{c|c|c|c}
&\det J_{\mathrm{selected}}&S(q_\tau)&S(q_{\mathrm{exp}})\\\hline
\mathrm{C1}&69120&103/37&107/39\\
\mathrm{C2}&184320&46/17&8/3
\end{array}
\]
The selected-minor numbers use this jet convention, rather than the raw-kernel convention in EN8/EN16.
The C1 exposing mixture has weights \(37/39,1/39,1/39\); the C2 weights are \(17/18,1/36,1/36\).
Run `node checks/contact-ray-exposure/verify_contact_seeds.js`. The executed exact rational checks test both six-jet systems, both selected minors, all eight losing probabilities, all three deterministic neighbours, and both CHSH values. They corroborate the examples; the universal exposure proof uses EN1–EN16.

## Consequence for finite source atlases

For finitely many sound source cones \(K_j\subseteq C\), an extreme ray can lie in \(\sum_jK_j\) only if one summand contains it. Each nonzero summand in a decomposition of an extreme ray is on that ray, and cones contain all its positive multiples. Conversely, covering every extreme ray covers the compact normalized base and hence \(C\).

The exposed-ray result makes the C1/C2 contact graphs the exact test objects for a proposed finite atlas, after the closed faces and symmetry images have been included. It does not prove that the order-940 source atlas, its order-7520 outer behaviour pencil, or the sixth-root candidate has this coverage.

References: [source geometry](../latex/chapters/04_geometry.tex), [contact signs and minors](../latex/chapters/05_entry.tex), [current construction route](EXACT_LIFT_ROUTE.md).
