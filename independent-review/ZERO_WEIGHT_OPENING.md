# Contact-compatible opening of the three zero-weight walls

27 October 2026. The earlier [C2 n=0 opening theorem](C2_N0_OPENING.md) extends to the C1 beta=0 and D=0 walls. The result is conditional on a Gram positive definite on **all four larger wall contact complements**. It proves a full all-testing PSD neighborhood from that premise, not whole-wall PSD selection or universal H1/H2.

## Fixed wall spaces

Normalize a finite nonzero contact by u=u_c U, v=v_c V and psi=diag(psi_c)Phi, so its coordinates are U=V=1 and Phi=(1,1,1,1). All scale factors belong to the source field. The source has lambda=0.

For the common C2 n=0 / C1 beta=0 wall, n=beta=0. Its larger 27-feature half-power support is the support in the C2 opening note. Equivalently it consists of U^(i/2)V^(j/2)Phi_l with i,j=0,1,2 and

\[i+j\ge2\text{ for }W,\quad j\ge i\text{ for }X,
\quad i+j\le2\text{ for }Z,
\]

with every Y feature retained.

For the C1 D=0 wall, n=D=0. Use the same conditions for W,X and all Y features, and replace the Z condition by i>=j. Both supports have parity block orders (13,5,5,4). Taking differences from the first feature in each block gives orders (12,4,4,3).

For **each** wall support, the full Gram coefficient map has shape72-by-104, rank66, and image exactly equal to the kernel of the six full contact derivatives. Specified 66-by-66 minors have determinant +/-2^54. The right inverse may be supported entirely on the 78 integer-block variables. Its maximum absolute row sum is141 for the n/beta wall and106 for the D wall. All keys, maps and right inverses are supplied in [fixed_maps.json](checks/zero-wall-opening/fixed_maps.json).

## Three positive contact-zero carriers

Write x_c=(psi_c)_X and z_c=(psi_c)_Z. Define positive source scalars

\[N=4n u_c^2x_c^2,\qquad
B=4\beta u_c^2v_c^2z_c^2,\qquad
T=4D v_c^2z_c^2.
\]

The following are sums of squares on U,V>0, for all four unrestricted state coordinates:

\[
\begin{split}
C_n=N\big[&(UX-\tfrac12V(X+W))^2
                 +UV(X-W)^2\big],\\
C_\beta=B\big[&(UVZ-\tfrac12(U+V)Z)^2
                 +U^2V(Z-W)^2+UV^2(Z-W)^2\big],\\
C_D=T\big[&(VZ-\tfrac12U(Z+W))^2
                 +UV(Z-W)^2\big].
\end{split}\tag{1}
\]

Every root vanishes at the normalized moving contact, so the value and all six first derivatives vanish. The coefficient N of U²X², B of U²V²Z², and T of V²Z² match the corresponding compulsory new coefficients of the complete SG2 parent exactly. The apparent extra U²VZ² and UV²Z² terms in the first beta square cancel against the two fractional squares. The mixed WX/WZ terms remain and must be matched; they are not discarded.

With all source weights formal, direct expansion proves:

* Fhat-C_n-C_beta has **only** keys in the n/beta wall map.
* Fhat-C_n-C_D has **only** keys in the D wall map.

Actual stationary sources make both residual coefficient vectors satisfy the six contact-jet equations. Thus every residual coefficient, including all state normals and mixed couplings, is reached by the corresponding constant rational wall map. These identities do not rely on a scalar test, positive diagonal budgets or an assumed domination by the carrier.

## PSD continuation theorem

Suppose a wall source has four reduced Grams Q_b^0 positive definite on the full contact complements above. Let x0 be their independent entries and r0 their coefficient vector. For a nearby actual stationary source, subtract the appropriate carrier sum, let r be its complete residual coefficient vector, and put

\[x=x_0+E(r-r_0).
\tag{2}\]

The fixed-map identities give exact coefficient matching. Only the 12-by-12 integer block changes. If mu is any proved positive lower bound for that block, the explicit condition

\[\|r-r_0\|_\infty<\mu/(24C),
\quad C=141\text{ or }106,
\tag{3}\]

keeps the corrected integer block at least mu I/2; the other three blocks remain unchanged and positive. The carrier is PSD because the opened weights are nonnegative. Hence (1)-(3) give a full PSD Gram for the entire parent. Continuity in the source/contact data proves a genuine source neighborhood; no continuity argument at a singular new row is needed because that row is supplied explicitly by (1).

Individual openings use 28 fixed features for n, 30 for beta, and28 for D. Joint n,beta opening uses31 features, with53 after multiplication by1+v². Joint n,D opening uses29 features, with50 after multiplication. Every weighted feature belongs to K364. On the strict C1 branches n=0; on the strict C2 branch beta=0. Joint openings are stated for nearby stationary seven-event sources and are not additional extreme-source classifications.

These are ambient PSD carriers, not individually Bell-source cone elements. The complete matched sum is the source certificate. The existing sign congruences and homogeneous extension give all physical tests and projective endpoints.

For a genuine permitted five-event wall contact, nonemptiness of its allowed strict native sides follows from the previously verified event-opening and native-entry argument L15. The construction here certifies the resulting actual stationary sources. It does not assert that an arbitrary perturbation of the Bell coefficients preserves a stationary contact or physical validity.

## What is completed and what remains

The same explicit mechanism now completes the local first-normal opening task at each of the three zero-weight walls whenever its full relative-PD wall hypothesis holds. It applies in particular to the stored larger-support wall certificates: eight C2 n=0, nine C1 beta=0 and eight C1 D=0 instances, with their actual source domains retained. The smaller18-feature wall witnesses do not automatically imply the larger-support PD premise.

Whole-wall relative-PD selection remains open. Admission away from the proved source neighborhoods, and the complete H1/H2 statements, also remain open. This theorem does not repair the refuted frozen-top C2 polynomial ansatz or prove a global residual for arbitrary strict sources.

Run the standalone exact map/carrier check with SymPy:

```bash
python independent-review/checks/zero-wall-opening/verify_zero_wall_opening.py
```

It reconstructs both constant coefficient maps, proves their ranks and six-jet characterization, expands the entire source parent and all carriers with formal source/contact quantities, and checks fixed K364 support. It does not assume whole-wall PSD. The proof above supplies the quantitative PSD continuation implication. This is a same-assistant derivation and exact replay, not an independent review.
