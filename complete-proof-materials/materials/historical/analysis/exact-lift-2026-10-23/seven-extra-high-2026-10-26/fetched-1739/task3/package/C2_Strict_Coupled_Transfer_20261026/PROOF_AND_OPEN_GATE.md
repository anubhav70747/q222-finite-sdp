# Strict native C2: positive signed-diagonal carriers, and a complete multiplied-parent SOS

**26 October 2026. Strict scope: genuine physical six-event C2 contact sources only.**

This is an exact mathematical construction and a precise noncompletion record. The fixed `Z6` all-interior admission assertion and the equivalent nonstrict 24-term dual inequality **remain unproved**. No alleged counterexample to `Z6`, new global lift, or invalid cancellation of a testing-dependent multiplier is claimed. The qualitative finite-SDP existence theorem and verified strict open-end `Z2` theorem are unchanged.

## 1. The native domain and exact full parent

Let `theta=(a,b,c,t,r)` satisfy the complete strict native domain `Omega`:

- `a,b,c,t>0`, `0<r<1`, `h=b-a^2 t>0`;
- `g=(1+c^2)b-(1-r+a^2 c^2 t)>0`, `j=1+c^2 r-t>0`;
- `1-c^2 h>0`, `b-a^2 c^2 r>0`, `N_*>0`, and the native `n` numerator `N_n>0`;
- `Delta>0`, with all six exact native coefficients `n,m,C,A,D,alpha>0` as specified by `c2_fixed_map.py`.

Set `d=1`, `lambda=beta=0`, `D_u=1+u^2`, `D_v=1+v^2`. The *entire* four-state parent, with unrestricted real state `(W,X,Y,Z)`, is

```
F = 4n D_u D_v X^2
  + 4m (vW-X+uvY-uZ)^2
  + 4C D_v (uW-Y)^2 + 4A D_u D_v Y^2
  + (4D+2)(uW+uvX-Y-vZ)^2
  + (4alpha+2) D_u(vY-Z)^2 + 2D_v(X+uZ)^2
  - 2D_u D_v Z^2.
```

The source contact is `(u_c,v_c;W_c,X_c,Y_c,Z_c)=(c,ac;b/a,-cr,c(b/a-at),1)`. For each fixed *real* source, normalize `u=c U`, `v=ac V` and `psi=diag(psi_c) phi`. This transformation is a source-dependent change of formal real polynomial coordinates, not an assertion that `U=u/c` belongs to the independent testing field of the qualitative tensor criterion. Every Gram expression below is ultimately a polynomial identity for all independent real `u,v>0` and all four states; no testing state normal is fixed.

The original 154-feature support-face map, 24 active target coefficients, rank-406 proof and compact dual domain are retained as previously verified inputs, not rederived here.

## 2. A positive all-domain transfer across the signed corner defect

Let `K=2c^4 r(2n+1)>0` and `q=a^2 r>0`, and denote the normalized state variables by `phi_l`. Exact extraction of the **complete parent**, followed by the native identity `2Ah=a^2r(2n+1)`, gives

```
[U^2 phi_Y^2] Fhat = 4A c^4 h^2/a^2 = K h,
[V^2 phi_X^2] Fhat = (4n+2)a^2 c^4 r^2 = K q,
[U^2 phi_Y^2] Fhat - [V^2 phi_X^2] Fhat = K(h-q).
```

Both coefficients arise from **unique diagonal entries** of the full `Z6` Gram, namely the features `(i,j,l)=(6,0,Y)` and `(0,6,X)`; no other pair of allowed features has their product exponent/state key. Thus both entries are forced positive, even when their difference changes sign.

A fixed-library, positive, contact-zero carrier for the two forced diagonals is

```
h>=q:  H+ = K[(h-q)(U phi_Y - phi_X)^2
              + q (U phi_Y - V phi_X)^2],
h<=q:  H- = K[(q-h)(V phi_X - phi_X)^2
              + h (U phi_Y - V phi_X)^2].
```

The comparison is **pure-source**, and both expressions are PSD and annihilate the normalized contact (all displayed monomial roots equal one there). At `h=q` they agree and have zero additional source-gap term. The carrier coefficients on `(U phi_Y,V phi_X,phi_X)` have diagonals `(Kh,Kq)` for both signs. This repairs the *sign ambiguity of the two extreme diagonal allocations*, not the complete 14-corner indefinite bridge.

The target has *zero* coefficients at keys `U phi_X phi_Y` and `V phi_X^2`. The `H+` carrier has `-2K(h-q) U phi_X phi_Y`; the `H-` carrier has `-2K(q-h) V phi_X^2`. Thus a PSD residual after subtracting this **particular** carrier would need to supply respectively

```
h>q:   +2K(h-q) [U phi_X phi_Y],
h<q:   +2K(q-h) [V phi_X^2],
```

while its two saturated extreme diagonal rows would necessarily vanish. In the complete sixth-root support, the remaining coefficient products have, respectively, **six** and **three** available feature-pair classes. They are enumerated *exactly* by `verify_c2_transfer.py`. The second requires fractional rows; none can be furnished by the already saturated `V phi_X` feature. This is a precise, finite residual PSD-transfer obligation, **not a proof that it is soluble or that every feasible Gram must dominate this chosen carrier**. The universally equivalent unproved condition remains the 24-term dual nonstrict inequality.

## 3. Five exact positive angular circuits on all of Omega

For C2, the SG5--SG6 determinant has exactly six strictly positive outer Laurent coefficients:

```
q(-1,0)=4 alpha m t5 t8,    q(-1,1)=4 m D t5 t6,
q(0,-1)=4 alpha A t7 t3,    q(0,1)=4 C D t5 t2,
q(1,-1)=4 n A t3 t4,        q(1,0)=4 C A t1 t4,
```

with `t1=1+2(n+D)`, `t2=1+2(alpha+m+A)`, `t3=1+2(C+D)`, `t4=1+2(alpha+m)`, `t5=1+2n`, `t6=1+2(alpha+C+A)`, `t7=1+2(n+m)`, `t8=1+2(C+A+D)`. These are strictly positive throughout the six-event domain. The missing `q(-1,-1),q(1,1)` vanish because `lambda=beta=0`.

Write `U=u^2`, `V=v^2`, `U_c=c^2`, `V_c=a^2c^2`, `s=U/U_c`, `t'=V/V_c`. Put

```
w_A=q(-1,0)/U_c,              w_B=q(-1,1)*V_c/U_c,
w_C=q(0,-1)/V_c,              w_D=q(0,1)*V_c,
w_E=q(1,-1)*U_c/V_c,          w_F=q(1,0)*U_c.
```

Stationarity of the *actual native full contact* gives the exact source balances

```
w_F+w_E=w_A+w_B,
 w_B+w_D=w_C+w_E,
 h_L=-w_A-w_B-w_C-w_D-w_E-w_F.
```

(`h_L` denotes the central Laurent coefficient, not the source gap `h`.) Define the following **source-only nonnegative** numbers:

```
kappa=max(0,w_A-w_F),    p=w_A-kappa=min(w_A,w_F),
ell=w_F-p=max(0,w_F-w_A),
q=w_E-kappa>=0,          z=w_D-kappa>=0.
```

The last two inequalities follow from `w_E-w_B=w_A-w_F=w_D-w_C`. The exact identity is

```
L(U,V) = p (s+s^-1-2)
       + q (s/t'+t'/s-2)
       + z (t'+t'^-1-2)
       + kappa (s^-1+t'+s/t'-3)
       + ell (t'/s+t'^-1+s-3).
```

The three opposite-pair terms become **one polynomial square each** after multiplication by `UV`:

```
UV(s+s^-1-2)         = [v(U-U_c)/sqrt(U_c)]^2,
UV(s/t'+t'/s-2)      = [sqrt(U_c/V_c)*V - sqrt(V_c/U_c)*U]^2,
UV(t'+t'^-1-2)       = [u(V-V_c)/sqrt(V_c)]^2.
```

For the first triangle, let

```
A1=(U_c V)^(1/3),
A2=(U V^2/V_c)^(1/3),
A3=(U^2 V_c/U_c)^(1/3).
```

Their cubes are the three positive terms in `UV(s^-1+t'+s/t')`, and `A1 A2 A3=UV`. For the second triangle set

```
B1=(U_c V^2/V_c)^(1/3),
B2=(U V_c)^(1/3),
B3=(U^2 V/U_c)^(1/3),
```

again with product `UV`. Each triangle has the exact **nine-square** identity

```
A1^3+A2^3+A3^3-3A1 A2 A3
  = (1/2)*sum_{i<j}sum_{k=1}^3 [sqrt(Ak)*(Ai-Aj)]^2,
```

and likewise for `B`. Its root monomials are fractional powers of the test variables with positive source-only coefficients, without a mixed/testing inverse. Thus `UV L(U,V)` is an explicit separated **21-square** polynomial identity, valid also when any `p,q,z,kappa,ell` vanishes. It is finite and covers the entire actual `Omega` (not a compact sampled subset). This is the source-independent *angular circuit list*, with piecewise source-only coefficients.

## 4. All four state normals: a complete multiplied-parent SOS

Write the original four-by-four state matrix as

```
P(u,v) = [[B(u,v), z(u,v)],
          [z(u,v)^T, d0(u,v)]],
xi=(W,X,Y)^T,
F=xi^T B xi+2 xi^T z Z+d0 Z^2.
```

The three-by-three block `B` is positive definite at every `u,v>0`, because its `Z=0` part contains the strictly positive signed core `2[(uW+uvX-Y)^2+(1+u^2)v^2Y^2+(1+v^2)X^2]` and the remaining native events are positive.

Here is a completely explicit **16-row polynomial factorization** `B=R^T R`. All entries below are three-vectors in the `W,X,Y` order; include one row for each listed choice of `f,g`:

```
(0, 2 sqrt(n) f g, 0),              f in {1,u}, g in {1,v},  4 rows;
2 sqrt(m) (v,-1,uv),                                                   1 row;
2 sqrt(C) g (u,0,-1),              g in {1,v},               2 rows;
(0,0,2 sqrt(A) f g),               f in {1,u}, g in {1,v},  4 rows;
sqrt(4D+2) (u,uv,-1),                                               1 row;
(0,0,sqrt(4alpha+2) f v),           f in {1,u},               2 rows;
(0,sqrt(2) g,0),                   g in {1,v},               2 rows.
```

Let `k=B xi+z Z`. The three-dimensional Cauchy--Binet identity says

```
adj(B) = sum_{1<=i<j<=16} (Ri cross Rj)(Ri cross Rj)^T.
```

The **exact entire four-state Schur identity** is

```
det(B) F = k^T adj(B) k + det(P) Z^2.
```

By SG5, `det(P)=16 (1+u^2)^2(1+v^2)^2 UV L(U,V)`. Combining the finite 120 wedge squares and the 21 circuit squares above gives the explicit 141-square identity

```
 det(B) F
  = sum_{i<j} [ (Ri cross Rj) dot (B xi+z Z) ]^2
    + 16 (1+u^2)^2 (1+v^2)^2 Z^2
       * [p*(opposite square 1) + q*(opposite square 2)
          + z*(opposite square 3)
          + kappa*(9-square triangle 1) + ell*(9-square triangle 2)].
```

This holds **identically in all four state coordinates**, for all independent positive test angles and all native strict C2 sources. It does not discard any `W/X/Y/Z` coupling. By `det(B)>0` it also proves physical nonnegativity, but that was already known by the native-source theorem.

**Crucial limitation:** `det(B(u,v;theta))` depends jointly on the source and the independent test angles. Dividing this identity by `det(B)` is **not a polynomial SOS cancellation**, and using such a quotient in the roots would introduce a forbidden testing-dependent denominator. The 141-square certificate is **not** a Gram of the original parent `F` in `Z6`, nor an affine, fixed globally sound SDP source cone. It is a stronger full-state algebraic adjunct to the angular determinant pencil, not a solution of the requested fixed-library admission.

## 5. The exact remaining implication

The still-unproved claim is that there exist **four PSD matrices in the fixed 154-feature C2 contact basis** whose 412 coefficient equations reproduce `F_theta` for every native theta. Equivalently the verified compact dual spectrum must satisfy the **nonstrict 24-term inequality**

```
forall theta in Omega, forall y in Y:
    sum_{24 active k} P_k(theta)*y[k] >= 0.
```

In particular the explicit source-positive carrier in section 2 needs an **all-row, PSD residual**; the two identified fractional/off-support channels are exact first residual equations but are *not* the only Schur constraints. The 141-square multiplied identity cannot be used to assert this missing residual positivity. No exact negative dual witness for the complete `Z6` cone was constructed, and no enlargement is justified.

## 6. Evidence and reproducibility

`verify_c2_transfer.py` deterministically checks: native identity; both uniquely fixed diagonal capacities; the sign-adaptive contact-square carrier and its off-support defects; admissible complete-Z6 coefficient pairings; exact full-parent `B=R^T R` with all state normals; an independent generic three-state Schur identity; and the five-ray Laurent decomposition and cubical AM--GM square identity. Eight exact checks passed. No numerical points, optimizer or solver were used.

`c2_fixed_map.py`, `C2_154_CONTACT_MAP.json`, and `C2_CLEARED_POLYNOMIAL_TARGET.json` are copied **unchanged** from the previously independently checked strict-interior package. SG5--SG6's full determinant coefficient identity is an input from the established native geometry rather than being expanded a second time in this script; the algebraic adjugate and angular-circuit conversions are explicitly proved here. The bundle does **not** claim to independently verify the SG5--SG6 determinant expansion afresh.
