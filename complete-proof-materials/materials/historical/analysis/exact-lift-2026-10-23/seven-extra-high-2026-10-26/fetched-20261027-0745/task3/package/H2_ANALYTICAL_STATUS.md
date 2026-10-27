# Universal strict native C2 admission into K364

**Scope: exact mathematical audit, 27 October 2026. Status: H2 NOT PROVED and NOT DISPROVED.**

This document does not claim a new all-source PSD factorization. It preserves the complete native four-state problem and its exact certificate conventions, and isolates the single remaining inequality. In particular, the exact counterexample to the earlier single-`tau` scalar completion and the frozen-top obstruction are *not* counterexamples to `K364`. The independently verified open C2 ends and separate qualitative finite-lift proof are preserved without alteration.

## 1. The precise assertion

Use the native variables `theta=(a,b,c,t,r)`, with `a,b,c,t>0`, `0<r<1`, and abbreviations

\[
\begin{gathered}
z=c^2,\quad s=1+zr,\quad h=b-a^2t,\quad
 g=(1+z)b-(1-r+a^2zt),\quad j=1+zr-t,\\
 k=1-zh,\quad B=b-a^2zr,\\
\Delta=1+(1+z)b+z(1+z)b^2+2z(1+z)br+z^2r^2+2zr
+a^2z\big[(1+z)b+1+zr(2-r)\big],\\
N_*=(1+z)(b+r)(b+a^2(1-t))-a^2(1+a^2z)(1-r)t,\\
N_n=(b+r)s(s+a^2zt)-r\Delta.
\end{gathered}
\]

The actual **strict** native domain `Omega` retains
\[
h,g,j,k,B,N_*,N_n,\Delta>0.
\]
The source weights, in order `(d,lambda,n,m,C,A,D,alpha,beta)`, are `(1,0,n,m,C,A,D,alpha,0)` with
\[
\begin{aligned}
n&=N_n/(2r\Delta),\\
m&=(1+a^2z)(1-r)s(s+a^2zt)/(2g\Delta),\\
C&=sN_* /(2a^2t\Delta),\\
A&=a^2s(b+r)(s+a^2zt)/(2h\Delta),\\
D&=(1+z)(b+r)B(s+a^2zt)/(2a^2j\Delta),\\
\alpha&=z^2(b+r)N_* /(2k\Delta).
\end{aligned}
\]
All six are strictly positive on `Omega`. The physical contact is
\[
(u_c,v_c;\psi_c)=(c,ac;b/a,-cr,ch/a,1),
\]
where `psi=(W,X,Y,Z)` and the full *original* four-state parent is
\[
\begin{aligned}
F_\theta={}&4n(1+u^2)(1+v^2)X^2+4m(vW-X+uvY-uZ)^2
 +4C(1+v^2)(uW-Y)^2\\
&+4A(1+u^2)(1+v^2)Y^2+(4D+2)(uW+uvX-Y-vZ)^2\\
&+(4\alpha+2)(1+u^2)(vY-Z)^2+2(1+v^2)(X+uZ)^2
 -2(1+u^2)(1+v^2)Z^2.
\end{aligned}
\]

Let \(Z_{6,12}\) be the fixed column of **364** monomials
\[
u^{i/6}v^{j/6}\psi_\ell,\quad0\le i\le6,\quad0\le j\le12,\quad\ell\in\{W,X,Y,Z\}.
\]
The requested H2 theorem is exactly
\[
\boxed{\forall\theta\in\Omega\ \exists P(\theta)\succeq0:\quad
(1+v^2)F_\theta=Z_{6,12}^{T}P(\theta)Z_{6,12}}
\]
coefficientwise in **independent** `u,v,W,X,Y,Z`, including every missing coefficient being zero. The source-only Gram may be rational or algebraic. No source-dependent testing denominator or degree is allowed.

## 2. Exact linear structure and what is still missing

The compact source map `COMPACT_K364_SOURCE_MAP.json` contains the exact source-affine coefficients of `(1+v^2)F` in all 3,250 allowable monomial keys. It is the previously reviewed 120-nonzero-entry map, copied unchanged from the user-supplied packet, **not** a recompilation. For Gram features \(e=(i,j,\ell)\) and \(f=(i',j',\ell')\), the exact coefficient key is
\[
\kappa(e,f)=(i+i',j+j',\min(\ell,\ell'),\max(\ell,\ell')).
\]
The polynomial identity of a positive Gram is all of the coefficient equalities
\[
T\,c(\theta)=\mathcal A(P),\qquad P\succeq0.
\]
If any equality at a zero target key is violated, an apparent candidate Gram is invalid even if it agrees at the contact and passes isolated positivity checks.

The full 32-positive-packet syzygy construction has a verified native-contact interpolation, but the `tau`-fixed scalar remainder is negative at the genuine strict seed `(a,b,c,t,r)=(1,1,1,3/4,1/4)`. Its old entire top Bob-degree packet is likewise excluded by an exact finite-dimensional projection. The four-root fractional `n` opening is PSD but introduces `uvWX` coefficients absent from the target, so it cannot be subtracted as a necessity from an arbitrary admitting Gram. No currently verified identity supplies a positive full-364 Gram for every source.

**Single outstanding implication**: prove the nonnegative pairing
\[
\boxed{\langle Tc(\theta),y\rangle\ge0\quad
\forall\theta\in\Omega,\ \forall y:\ \mathcal A^*y\succeq0.}
\]
This is equivalent to the requested Gram membership because the finite-dimensional SOS coefficient cone is closed (the already-reviewed positive feature-moment trace bound). The moment matrices are **not** assumed to have representing measures. An exact negative pairing would disprove this particular fixed cone; a negative scalar remainder from a chosen decomposition does not.

A universally sufficient primal resolution instead consists of a finite source-domain partition and algebraic matrices \(P_k(\theta)=R_k(\theta)^TR_k(\theta)\), or exact source-algebraic LDL factorizations with all diagonal signs established, such that *every* coefficient key agrees identically after the native source relations. Such a construction was not found in this run.

## 3. Reproducibility: exact checker, without old large-matrix recompilation

`verify_primal_factor.py` is a newly written, standalone **verification-only** program. It parses a proposed finite list of exact factor roots in the fixed 364 features, evaluates the actual rational native six-event source, accumulates coefficient products sparsely, and checks the entire polynomial identity including all zero target keys. It additionally checks that each factor vanishes at the genuine source contact. Thus any result of `EXACT SOS WITNESS VERIFIED` is a correct *pointwise* certificate in `K364`. It does not prove all-source H2, and the supplied package contains **no all-source Gram formula or strict-native witness**.

Execute `python verify_primal_factor.py --selftest` for a 6-root, single-event test of the coefficient convention. This self-test is not a C2 admission test. To validate a proposed strict C2 certificate, execute `python verify_primal_factor.py candidate.json`, where `candidate.json` has `theta` as five exact rational strings and `factors` as lists of `[i,j,l,coefficient]`. No missing `.pkl`, dense Gram, numerical sweep, SDP solver, or testing-dependent inverse is used.

## 4. Independent factorization avenue, carefully separated from H2

The existing global physicality theorem makes the 4-by-4 matrix \(F_\theta/[(1+u^2)(1+v^2)]\) a nonnegative **matrix-valued two-torus trigonometric polynomial of degree (1,1)** under the Cayley coordinates `u=tan(phi/2)`, `v=tan(chi/2)`. The operator-valued two-variable Fejér–Riesz theorem of Michael A. Dritschel, *Mathematische Annalen* 391 (2025), Theorem 4.1, implies for each **fixed** real source an analytic SOS spectral factorization with one angular degree bounded by 1 and the other a finite, potentially source-dependent integer `N(theta)`; see DOI 10.1007/s00208-024-02895-9. Clearing the Cayley denominators yields a real polynomial SOS of `(1+u^2)^(N(theta)-1) F_theta` with Bob polynomial degree at most 1.

This is a **per-source rational-factorization fact**, not a fixed finite library: the cited theorem does **not** bound `N(theta)` uniformly over `Omega`, and it does not supply the specific fixed `Z_(6,12)` Gram. It is not substituted for H2 or for any source-affine finite SDP. The proof of H2 remains open exactly as in Section 2.

## Conclusion

The strongest *target-relevant* statement established here is an exact certificate specification and verifier for the full K364 map. **No universal strict C2 PSD construction or dual separator was obtained.** The source-verified strict open-end theorem, the independently reviewed qualitative existence theorem, and the conditional materialized outer SDP retain their previously stated scopes. No repo files were changed.
