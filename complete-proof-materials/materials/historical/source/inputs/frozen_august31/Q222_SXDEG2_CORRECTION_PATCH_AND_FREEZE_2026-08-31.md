# Q222 sxdeg = 2 — surgical correction patch and freeze overlay

**Freeze date:** 31 August 2026  
**Status:** `LOCK_ZERO_COMPUTER_THEOREM_LEVEL_CORRECTED`  
**Claim:**
\[
\boxed{\operatorname{sxdeg}(\mathcal Q_{2,2,2})=2.}
\]

This file is a **surgical overlay** for the exact 30 August parchment bytes. It changes only the six presentation/self-containment defects found in the 31 August hostile line-by-line audit. No theorem module, child classification, valuation argument, tensor-rank argument, or lower-bound argument is otherwise altered.

## Base file guarded by hash

Apply this overlay only to

`Q222_SXDEG2_CONSOLIDATED_PARCHMENT_FINAL_2026-08-30.md`

with SHA-256

`e268a6ae1c11dd02b5819064fa7e8187173074050d92f3204e703dd46f574731`.

The corrected theorem-level proof is defined canonically as **those exact base bytes plus the six replacements below**. This composite definition avoids silently reserializing or perturbing any of the 2,415 unaffected lines.

---

# Repair 1 — Section 1.2: restore the unique-visible-facet proof

Immediately after

> The PR behavior associated with \(s\) exposes exactly the eight positivity normals with \(ab=-s_{xy}\). Their convex hull is a facet \(F_s\) of \(G\),

insert:

For completeness, this facet is the **unique** facet of \(G\) visible from \(C_s\). Facets of \(G\) are indexed, by polarity, by the vertices of the binary \(2\times2\) no-signalling polytope \(G^\circ\): the sixteen deterministic behaviors and the eight PR behaviors. For a deterministic behavior \(q\), the CHSH bound gives
\[
\langle q,C_s\rangle\le1.
\]
For a PR vertex \(q_t\),
\[
\langle q_t,C_s\rangle
=\frac12\sum_{x,y}t_{xy}s_{xy}
=
\begin{cases}
2,&t=s,\\
-2,&t=-s,\\
0,&t\ne\pm s.
\end{cases}
\]
Thus only \(q_s\) defines a supporting hyperplane of \(G\) from which \(C_s\) lies strictly outside. Hence \(C_s\) sees exactly \(F_s\).

The subsequent average identity, pairwise-midpoint identity, and visible-star lemma are unchanged.

---

# Repair 2 — Section 4: restore the real half-angle measurement chart

At the start of Section 4, before

\[
D_u=1+u^2,\qquad D_v=1+v^2,
\]

insert:

Use the real half-angle qubit chart
\[
A_0=B_0=Z,
\]
\[
A_1(u)=\frac{(1-u^2)Z+2uX}{1+u^2},
\qquad
B_1(v)=\frac{(1-v^2)Z+2vX}{1+v^2}.
\]
For a real homogeneous two-qubit state vector \(\psi\), the rows below are the exact denominator-cleared outcome amplitudes in this chart.

This makes the phrase “direct substitution” in (4.3) self-contained. Equations (4.1)–(4.3) are otherwise unchanged.

---

# Repair 3 — Section 7: correct the F1 Bell gauge definition

Replace the F1 display

\[
\phi_1=(u,s,as/p,1/p),
\qquad
G_{F_1}=p^2G_E\circ\phi_1,
\qquad
\det D\phi_1=\frac{s^3}{p^3},
\]

by the fully typed formula

\[
\boxed{p:=Bs^2-u,}
\]
\[
\boxed{
\phi_1(s,u,a,B)
=\left(u,s,\frac{as}{p},\frac1p\right),
}
\]
\[
G_{F_1}=p^2G_E\circ\phi_1,
\qquad
\boxed{\det D\phi_1=\frac{s^3}{p^3}}.
\]

Consequently the cofactor multiplier remains
\[
\boxed{\lambda_{F_1}=s^3p^7>0}
\]
on the open physical F1 gauge. This repair restores the variable dependence required for the printed Jacobian; it does not change the collar theorem.

---

# Repair 4 — Section 7 / Section 14: make the zero-denominator routing explicit

Replace

> Exact zero of a gauge denominator is never treated by division; it is routed before the quotient is taken.

by

> Exact zero of a gauge denominator is never treated by division. Finite affine numerator/denominator factors are routed first through the explicit strict children in (7.9)–(7.10). If the denominator vanishes because the measurement leaves the affine chart, no quotient is taken: the contact is routed to the homogeneous one-projective or simultaneous-projective analysis of Sections 9 and 11–13. Thus every reciprocal or birational division is performed only on its declared open gauge.

In Lemma 14.1, read the sentence

> A zero gauge denominator is, by definition, one of the printed strict boundary children; no quotient is taken there.

as

> A zero gauge denominator is processed before division: finite algebraic boundary factors route to the printed strict children, while an affine-chart denominator zero routes to the homogeneous projective analysis of Sections 9 and 11–13. No quotient is taken on either boundary.

---

# Repair 5 — Section 12: correct the factor-typing sentence after the B0 four-square

Replace

> Every coefficient is nonnegative and lies wholly in the point factor. Hence (12.1) is already a tensor SOS, with rank-one or rank-two roots after the trivial support-function split.

by

> Every coefficient in (12.1) is nonnegative and **separates factorwise** into a nonnegative pure-point factor times a nonnegative pure-support square or elementary support factor. Explicitly, the only support dependence in the scalar prefactors is through \(1+U^2\), \(1+V^2\), and \(U^2\), while \(A,T,p,w,h,\varepsilon_0\) are point-side. Absorbing the pure-support factors into the displayed support vectors therefore gives literal tensor-square roots of separation rank at most two. Hence (12.1) is already a valid tensor SOS.

No change is made to (12.1) or to the \(h=0\) boundary argument.

---

# Repair 6 — reader-route A17 R23: delete the stale K3-dependent pivots

The companion reader route must **not** use
\[
q_\beta=
\frac4{T^3(1+V^2)}
+\frac{3K_3}{2T^4V^2(1+V^2)^2},
\]
\[
q_{e_1}=
\frac4{T^3V^2(1+V^2)}
+\frac{K_3}{4T^4V^2(1+V^2)^2}.
\]
Those expressions are a stale historical regression and are superseded by the centered calculation in the authoritative parchment.

The frozen reader-route equation (R23) is
\[
\boxed{
q_\beta=\frac4{T^3(1+V^2)}>0,
\qquad
q_{e_1}=\frac4{T^3V^2(1+V^2)}>0.
}
\tag{R23}
\]
Equivalently, with \(S=1+V^2\),
\[
q_\beta=\frac4{T^3S},
\qquad
q_{e_1}=\frac4{T^3V^2S}.
\]

The critical value \(m=u^3K_3\ge0\) is still retained exactly and restored separately by the pure-side square in (10.33). The correction therefore removes a stale pivot formula without deleting the positive infinitesimal critical value.

---

# Frozen theorem statement after the repairs

With the six replacements above incorporated, the proof chain is unchanged:

\[
\operatorname{sxdeg}(K_\star)\le2
\Longrightarrow
\operatorname{sxdeg}(\mathcal B_Q)\le2
\Longrightarrow
\operatorname{sxdeg}(\mathcal Q_{2,2,2})\le2,
\]
while the attained cubic correlator boundary
\[
C_{11}=4t^3-3t,
\qquad
1/\sqrt2\le t\le1,
\]
shows that \(\mathcal Q_{2,2,2}\) is nonpolyhedral and therefore
\[
\operatorname{sxdeg}(\mathcal Q_{2,2,2})\ge2.
\]
Hence
\[
\boxed{\operatorname{sxdeg}(\mathcal Q_{2,2,2})=2.}
\]

## Dependency firewall retained

The corrected frozen theorem does not use as a premise:

- numerical SDP output;
- CAS output or sparse machine Gram replay;
- Newton/Pareto enumeration;
- a machine-generated contact atlas or production-cell table;
- the invalid generic filtered-collar/R2.4 implication;
- a weighted/mixed Bell generator as an SOS root;
- identification of independent point and support tensor copies;
- “different specializations imply positive residue”;
- positivity in a valued field as a substitute for positive residue;
- deletion of the A17 positive infinitesimal critical value.

## Scope not claimed

This freeze proves only the theorem-level extension-degree statement. It does **not** freeze:

1. a literal solver-ready affine \(L_3\) formulation;
2. the exact minimum number \(M_{\rm SOC}\) of nondegenerate Lorentz blocks.

Those remain separate constructive/minimality problems.
