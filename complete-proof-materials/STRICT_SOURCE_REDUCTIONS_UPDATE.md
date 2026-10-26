# Strict-source reductions and the remaining positivity statements

The new C1 and C2 packages supply exact reductions. Neither proves universal admission to a fixed Gram cone. Their original files are preserved unchanged; the qualifications below govern their use.

## C1: two coupled Schur conditions

Adding sqrt(uv)X and sqrt(uv)Y gives a fixed 18-feature library contained in Z6. The generated 46-by-49 coefficient map has rank 40. Matching its four forbidden parity coefficients leaves the 42-by-45 invariant map, rank 36, with nine free gauges. Its invariant linear image is unchanged. The enlargement adds PSD freedom through a second nonnegative half-power charge delta alongside gamma.

Exact coefficient matching holds for the complete four-state parent and independent tests. Strictly positive forced pivots reduce PSD to two coupled 4-by-4 Schur residuals, together with gamma,delta>=0. Only those strictly positive pivots are inverted, so singular residual certificates remain included. The unresolved assertion is the existence of a simultaneous choice of all nine gauges for every physical strict native source with 0<k<=1. L17 already transports any such Z6 certificate to the opposite orientation.

The four-contact Fourier constraint now bounds gamma+delta from below. Delta does not remove every constraint: H7[2,3]=-delta/2 gives the necessary bound

    delta^2 <= 4 H22 H33
             = 16(2 beta+d+2m)(2D+d) u_c^2 v_c^2 X_c^4.

The old upper bound on gamma also remains. Scalar compatibility is therefore necessary; it is not sufficient for the simultaneous Schur residuals to be PSD. Saturation imposes the corresponding PSD kernel conditions. The added charge cannot be treated as an unlimited correction.

## C2: a positive carrier and a polynomial multiple

The sign-dependent carrier is a PSD sum of contact squares using the two forced endpoint capacities. A residual after subtracting this particular carrier must satisfy the stated zero-row and off-support restoration equations. This is a condition for that selected decomposition. An arbitrary feasible Gram need not dominate the chosen carrier, so the residual problem is not established as an equivalent characterization of all PSD certificates. Its PSD completion remains unproved.

A separate exact identity gives a full-state sum of at most 141 squares for det(B)F, where B is the positive 3-by-3 state block. The 120 adjugate squares and 21 angular squares have the required nonnegative weights on the stated physical domain. The determinant multiplier depends jointly on the source and the tests. Cancelling this positive multiplier does not produce a certificate for F in a fixed polynomial SOS space or the separated tensor ring. Its source dependence also prevents treating the displayed multiplied identity as an affine source-cone constraint. A fixed multiplier or denominator reduction would require an additional proof.

Consequently neither the carrier nor the multiplied-parent identity proves the full strict C2 dual inequality. Their exact identities and support calculations are verified within these scopes.

## C2 direct syzygies: the selected scalar completion fails

The new construction gives a full four-state identity for (1+v^2)F using 32 fixed polynomial roots and a scalar remainder. Its 29 contact-kernel directions, one additional high-degree direction and positive highest-Bob-degree pivot are verified. The generic deformed roots remain in K364's fixed feature span. None of these facts establishes nonnegativity of the scalar remainder.

The particular contact interpolation, even with an arbitrary real coefficient tau on its single additional direction, fails at the genuine strict native source

    (a,b,c,t,r)=(1,1,1,3/4,1/4),
    (n,m,C,A,D,alpha)=(1/8,3/8,1/6,5/8,3/8,1/6).

For this source the exact remainder on v=1 is

    S_tau(u,1)=-(u-1)^2 [379(tau-738/379)^2+295599/379]/30.

It is strictly negative for every real tau and every positive u other than 1. In particular the proposed tau=2 gives S_2(1/2,1)=-781/120, despite its positive v^4 coefficient 6/7. Thus universal positivity/SOS for this selected one-direction family is refuted, not merely unproved. This is a failure of that sufficient scalar-completion ansatz, not a counterexample to K364 admission.

A correction is also required when using the original JSON's interpolation shorthand: besides lambda_m=eta4/u_c, the exceptional coefficient is lambda_Yu=eta9+v_c*eta4. The original proof and executable code contain this adjustment; the abbreviated JSON prose omits it. The original files remain unchanged.

A repair must change the contact roots' angular derivatives through additional admissible zero-contact syzygies, or supply a different full Gram. At the exact seed, the current completed-root squared u-derivative norm is 1523/40, whereas the parent permits at most 289/24. The necessary local filter is Hessian((1+v^2)F at the fixed contact state)-2(dG)^T(dG)>=0. Passing that filter would still not establish global scalar SOS or H2. Retuning tau alone cannot repair this source.

- [Original direct C2 syzygy construction](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1839/task3/package/PROOF_AND_OPEN_GATE.md)
- [Original syzygy proof and checker package](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1839/task3/original-package.zip)
- [Independent analytical review and strict-source obstruction](materials/historical/analysis/c2_direct_syzygy_analytic_verification_2026_10_26.txt)
- [Independent exact replay and coverage](materials/historical/analysis/c2_direct_syzygy_exact_verification_2026_10_26.txt)
- [Independent exact obstruction checker](materials/historical/analysis/c2_direct_syzygy_exact_obstruction_2026_10_26.py)
- [Exact strict-source obstruction data](materials/historical/analysis/c2_direct_syzygy_exact_obstruction_2026_10_26.json)

## Complete polynomial u-jet at the strict seed

The follow-up package exhausts the 192 coefficients of all 32 polynomial charge rows of bidegree (1,2), at the same strict source. Six linear constraints describe the attainable 32-component first u derivative. The exact weighted projection reaches the sharp parent energy 289/24, so every nonnegative completion in this ansatz must attain that derivative. This result concerns the first u derivative at one source; both angular derivatives and global positivity remain to be controlled.

Keeping the old complete v^2 and uv^2 allocations fixed, the full system leaves 38 lower-degree gauges. Exact least squares proves that their two necessary budgets have sum at most -44072179/6288328800. Thus all lower polynomial syzygies together cannot repair that frozen top allocation. Releasing v^2 gives positive optimal combined margins for the tested relaxed systems, which only removes this sum obstruction: it does not establish simultaneous nonnegativity of the two budgets or scalar SOS. A sufficient construction must allow the high-Bob charge vector to vary or use another full K364 Gram. H2 remains open; this seed calculation does not exclude K364.

- [Original complete polynomial u-jet and frozen-top result](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1939/task3/package/PROOF_AND_STATUS.md)
- [Original complete-jet package](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1939/task3/original-package.zip)
- [Independent first-u-jet and frozen-top analytical review](materials/historical/analysis/c2_complete_jet_analytic_verification_2026_10_26.txt)
- [Independent exact replay and direct-parent capacities](materials/historical/analysis/c2_complete_jet_exact_verification_2026_10_26.txt)

## Next implication

The next C1 task is a simultaneous native-domain PSD gauge construction; the next C2 task is full Gram admission or its equivalent dual positivity statement. Another coefficient-rank calculation or isolated carrier is insufficient. A different fixed globally sound library remains an admissible route if it gives a complete proof. A separating witness against Z6 is not a prerequisite for choosing one.

The two strict-family admission inputs remain open. L15/L17 retain their conditional density, reflection and assembly results; the qualitative V01-V06 conjunction is unchanged. No new completion box is closed by this batch.

## Evidence

- [Original Z18 and Schur reduction](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1739/task5/package/PROOF.md)
- [Original complete C1 package](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1739/task5/original-package.zip)
- [Independent Z18 analytical review](materials/historical/analysis/c1_z18_schur_analytic_verification_2026_10_26.txt)
- [Targeted exact replay](materials/historical/analysis/latest_transfer_exact_verification_2026_10_26.txt)
- [Original C2 coupled-transfer and multiplied-parent identities](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1739/task3/package/C2_Strict_Coupled_Transfer_20261026/PROOF_AND_OPEN_GATE.md)
- [Original complete C2 package](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1739/task3/original-package.zip)
- [Independent C2 analytical review](materials/historical/analysis/c2_coupled_transfer_analytic_verification_2026_10_26.txt)
- [Targeted exact replay](materials/historical/analysis/latest_transfer_exact_verification_2026_10_26.txt)
