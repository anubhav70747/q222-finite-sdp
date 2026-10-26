# Strict oriented C1: exact charge-capacity theorem and the unresolved simultaneous Gram completion

**Date:** 26 October 2026. **Immutable evidence:** main `3d193cc25ae2fcdfd3817c0bf5ebd54dfada7b55`, the recovered Z18 constructor and its independently checked source identities. **Status:** universal strict-C1 admission to Z18, Z6, or K364 is **not proved** here. No optimizer, numeric sweep, witness fitting, new coefficient-map rank, or Git write is involved.

## 1. The universal target and the only local Gram gate

For the true, oriented C1 source write `theta=(k,f,g,a,t)`, `0<k<=1`, satisfying all original strict source gaps, seven genuine SD1 coquotients, positive retained event weights, and the nonzero native contact. The full parent is the original SG2 polynomial `F_theta(u,v;W,X,Y,Z)` for **arbitrary independent positive u,v** and **all four state coordinates**. The original feature-gauge source and complete coefficient identity are in `original_z18/PROOF.md` and `original_z18/c1_z18_extension.py`; they are inputs, not rederived here.

In the fixed 18-root cone, the exact matching equation is

```
Fhat_theta = r7^T H7(theta;s,gamma,delta) r7
           + r5^T H5(theta;s) r5
           + gamma*uv*(W-Z)^2+delta*uv*(X-Y)^2,
```

where `s=(s0,...,s6)`, `gamma,delta>=0`; state and angles in this display have been rescaled by **source-only nonzero contact values**. The forced positive H7[1,2], H7[4,4] and H5[4,4] have already been exactly eliminated. Let `S7(theta;s,gamma,delta)` and `S5(theta;s)` be the two exact **4x4** Schur matrices in the original constructor. Then the precise unsettled universal assertion is

```
FOR EVERY strict oriented native theta (0<k<=1),
THERE EXIST real s0,...,s6 and gamma,delta>=0
SUCH THAT  S7(theta;s,gamma,delta)>=0 AND S5(theta;s)>=0.
```

This is *necessary and sufficient for this particular Z18 construction*, including singular feasible Grams. It is not known on the full oriented domain. Matching the parent and satisfying scalar charge inequalities do not imply either PSD condition. The earlier `k=4` Z16 separator does not settle the oriented Z18 problem.

The weaker alternative `K364` requires, for every strict oriented theta, a **different** PSD matrix `P` satisfying

```
(1+v^2)*F_theta(u,v;psi)=Z_(6,12)(u,v;psi)^T P Z_(6,12)(u,v;psi)
```

coefficientwise. This fixed cone is globally sound and admits reflected sources as in the verified transport theorem, but no source-wide P or exact separating witness is supplied. The old Z18 gauges do not parameterize *all* K364 Grams. A failure of Z18 therefore does not imply a failure of K364.

## 2. A new exact source-rational charge selection in a genuine oriented singular regime

This theorem addresses the charge-compatibility step, **not** the seven-gauge Gram completion. It is useful precisely because it proves source signs before using the charge caps, rather than varying event weights independently.

**Theorem (oriented small-f charge capacity).** Fix `x,z` satisfying

```
0 < x <= 1,    0 < z < 1,
B(x,z):=x*z^2-(1+x)^2*z+x*(1+x) > 0.
```

Take the **actual native** C1 subfamily

```
k=sqrt(x),  g=sqrt(z),  a=t=1,  f>0.
```

For every such fixed pair `x,z`, all original strict native conditions and original nonzero denominators hold for sufficiently small `f>0`, with original normalized `C=2`, `lambda=n=0`, and `d,m,A,D,alpha,beta>0`. Put

```
S_X=(4*beta+2*d)*v_c^2*X_c^2,
S_Y=4*A*u_c^2*Y_c^2,
S_Z=4*beta*u_c^2*v_c^2*Z_c^2,
L_triangle=max(0,sqrt(S_X)-sqrt(S_Y)-sqrt(S_Z),
             sqrt(S_Y)-sqrt(S_X)-sqrt(S_Z),
             sqrt(S_Z)-sqrt(S_X)-sqrt(S_Y)),
U_delta=16*(2*beta+d+2*m)*(2*D+d)*u_c^2*v_c^2*X_c^4.
```

Then for all sufficiently small positive `f`, the **explicit positive source-rational charges**

```
gamma = 0,        delta = S_Y = 4*A*u_c^2*Y_c^2
```

satisfy **both charge-only inequalities** `gamma+delta>=L_triangle^2` and `delta^2<=U_delta`. The inherited nonnegative-gamma cap is automatic at `gamma=0`. These statements impose no comparison with an independent testing angle. **They do not assert S7>=0, S5>=0 or K364 admission.**

### Proof and exact factorization

The actual seven SD1 coquotients, with `K=1+x`, give (a limit at fixed `x,z`):

```
d/f^2    -> 4*(K^2-x*z)/K,
m/f^2    -> 2*(1-z)/z,
A/f^4    -> 2*(K-x*z)/K,
D        -> 2/[x*(1-z)],
beta/f^2 -> 2*B(x,z)/(z*K).
```

Each printed leading term has positive sign on the theorem domain, as do the remaining coquotients. The original native gaps tend to `Lc=1, Hc=1, Qc=1-x*z, M=x, J=x*(1-z)`; the other strict factors tend to `1-z,1` and the positive polynomials `P4=K-x*z`, `Palpha=K-z`, `Pd=K^2-x*z`. Thus the source is genuinely strict for small `f>0`, not an independent tuple of weights. The contact is precisely

```
u_c=sqrt(x), v_c=sqrt(x)*f*sqrt(z),
W_c=sqrt(z)/f, X_c=-sqrt(x)*(1-z),
Y_c=sqrt(x)*sqrt(z)*(1-f^2)/f, Z_c=1.
```

The exact source expansion yields

```
S_Y/f^2 -> Y2 = 8*z*x^2*(K-x*z)/K >0,
S_X/f^2 -> 0,   S_Z/f^2 -> 0,
U_delta/f^4 -> D4 = 256*x^3*(1-z)^3*(K-z)>0.
```

The source-sign inequality `B(x,z)>0` implies

```
z < zmax := x*K/(x^2+x+1) <= 2/3.
```

Indeed `B = xK-z*(K^2-x)+x*z*(z-1)`, whose last summand is strictly negative. The second inequality is equivalent to `(x-1)*(x+2)<=0`. Define

```
Dpoly=4*K^2*(1-z)^3*(K-z)-x*z^2*(K-x*z)^2.
```

On `0<=z<=zmax`, its first term strictly decreases; its second term increases because `K-2*x*z >= 1-x/3 >= 2/3`. At the maximal permitted `z`, exact polynomial division gives

```
Dpoly(x,zmax)
 = (1-x)*K^3*(x^5+4*x^4+7*x^3+8*x^2+4*x+4)/(x^2+x+1)^4 >=0.
```

Because the native `z<zmax` is strict, `Dpoly>0`, including `x=1`. Finally the **exact** identity

```
D4 - Y2^2 = (64*x^3/K^2)*Dpoly > 0
```

proves `sqrt(D4)>Y2`. Thus `sqrt(U_delta)>S_Y` for all sufficiently small `f>0`. Meanwhile `sqrt(S_Y)=Theta(f)` and `sqrt(S_X),sqrt(S_Z)=O(f^2)`, so `L_triangle=sqrt(S_Y)-sqrt(S_X)-sqrt(S_Z)` for sufficiently small `f`, and `L_triangle^2<S_Y`. Therefore `L_triangle^2<delta=S_Y<sqrt(U_delta)`, as required. All inequalities are exact symbolic or follow from strict positive limits over fixed real parameters. No uniform lower cutoff in f for noncompact x,z is claimed. On a compact strict set of x,z, standard strict margins give a common small-f cutoff for this **scalar** statement only.

The accompanying `verify_source_charge_margin.py` reconstructs all source limits and all displayed polynomial identities from the original native dictionary. No SDP/solver/sample computations enter.

## 3. What this does **not** finish

The source-rational choice `gamma=0,delta=S_Y` establishes the **charge-only** inequalities on the stated oriented subfamily. It may fail either full Schur matrix; no gauge selection `s0,...,s6` with exact positive 4x4 factorizations has been derived, not even universally on this subfamily. In the full oriented strict C1 domain the charge selection itself is not established. The first unproved Z18 implication is exactly the quantified simultaneous `S7,S5` PSD condition in §1, not stationarity, coefficient matching, any one principal minor, or charge compatibility. The broader sufficient K364 coefficientwise-PSD condition in §1 is likewise open; its compiler is assigned elsewhere, and this note does not purport to materialize its matrices.

The independently checked strict C2 open-end result, the 215-component qualitative result, source density, reflection transport and the full physical orbit assembly retain their **previous conditional scopes**. No global fixed SDP or full strict-C1 admission is claimed.
