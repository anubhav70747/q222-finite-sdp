# Reflection transport and a common fixed cone

The original theorem and checker package are preserved unchanged. This note states the verified extension and the symmetry correction required when using that theorem.

The source reflection is the fixed linear involution sigma(d,lambda,n,m,C,A,D,alpha,beta) = (d,n,lambda,C,m,alpha,beta,A,D). Its full four-state identity is

    u^2 F_(sigma c)(1/u,v;G(v)psi) = (1+v^2) F_c(u,v;psi).

The testing-dependent state matrix satisfies G(v)^2=(1+v^2)I. The native contact transformation sends k to 1/k, preserves the complete strict physical C1 domain and restores the C=2 representative by the positive source factor 2/m. Thus the identity orientation and this reflection cover all strict C1 rays from 0<k<=1. The source-cone map remains linear; the representative normalization is used only when invoking the native admission theorem.

## Full sixth-root extension

Let z6 contain u^(i/6)v^(j/6)psi_l for 0<=i,j<=6 and four state coordinates. Let zext contain the same monomials for 0<=i<=6, 0<=j<=12. For every one of the 196 original features,

    h6 = u z6(1/u,v;G(v)psi) = L6^T zext.

The exact integer matrix L6 has 364 rows, 196 columns and 392 signed nonzero entries. Its columns and their polynomial identities are independently checked. Therefore F_(sigma c)=z6^T Q z6 with Q>=0 implies

    (1+v^2)F_c = zext^T (L6 Q L6^T) zext.

For an unreflected z6 certificate, multiplication by 1+v^2 is represented by the sum of its unshifted and v-shifted Gram embeddings. Both use the same 364-feature space. The theorem can consequently assume oriented C1 admission in the full Z6 space, rather than only in Z16. Any smaller fixed library can be used once its embedding in Z6 is established. This extension proves transport; it supplies no universal oriented C1 or strict C2 PSD certificate.

Define K_ext by c>=0 and coefficientwise equality (1+v^2)F_c=zext^T P zext with P>=0. This is a fixed rational-data cone, globally sound for all independent tests and unrestricted states, including axes and projective endpoints. Positive feature moments prove its closedness. The imported strictly positive Z6 Gram extends to a strictly positive 364-row Gram because the unshifted and v-shifted row ranges cover the entire extended grid.

## Conditional global consequence and required symmetry

The weaker sufficient hypotheses are direct admission to K_ext=K364:

- H1: every physical strict C1 source in the actual oriented domain 0<k<=1 belongs to K364.
- H2: every physical strict C2 source in its complete native domain belongs to K364.

Neither requires an unweighted Z6 representation of F. Set Kbar=K364+sigma(K364). Both summands are fixed, closed and globally sound; the deterministic-coefficient bound closes their finite sum. The native reflection sends every opposite-orientation strict C1 ray to H1, while H2 supplies C2. Thus all strict sources enter Kbar, and L15 gives the required boundary sources there. This does not assert sigma(K364)=K364 or individual membership of every unoriented C1 source in K364.

For the full physical group Gamma, Gamma*sigma=Gamma, so

    sum_g g Kbar = sum_g g K364 + sum_g g sigma(K364)
                = sum_g g K364.

The last equality uses conicity and addition closure of the orbit sum. Hence the same 128-orbit global candidate needs no extra blocks. Retain its event and SG12/SG13 face-orbit blocks. [Direct common-cone reduction](materials/historical/analysis/direct_common_cone_admission_reduction_2026_10_26.txt) proves this implication. Both H1 and H2 remain unproved. The Z6 transport above remains one sufficient route to direct K364 admission; it is not a necessary requirement.

Equation (16) of the original report refers to a group Gamma "of SG8". The stabilizer of a single CHSH orientation is insufficient for the complete Bell cone. Use the full 128-element physical relabeling group, including party exchange, or eight PR-sector representatives together with the canonical stabilizer. Apply those fixed linear maps to every required cone before taking the global sum. This correction is part of the conditional global statement.

For coefficient maps T c=A P, the attained source-dual representation is

    K_ext* = { T^T y+r : A* y>=0, r>=0 }.

The strict source and Gram point bounds convergent dual representatives after quotienting the joint kernel. The Bell target vector and source orthant must be retained. The theorem specifies this map and its coefficient indexing, but the original JSON does not materialize the complete sparse T, A or final symmetry-assembled pencil. The 364-feature count and at most 3250 coefficient keys therefore describe a conditional construction, not a delivered exact global SDP.

The two strict-family admission theorems remain open. So do evaluation of every global matrix, verification of the resulting two inclusions, and the short complete proof. The existing qualitative-existence conjunction V01-V06 is unchanged.

## Evidence

- [Original reflection and conditional common-cone theorem](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1719/package/C1_REFLECTION_TRANSPORT_AND_FIXED_CONE_THEOREM.md)
- [Original complete reflection package](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1719/original-reflection-package.zip)
- [Independent analytical verification](materials/historical/analysis/c1_reflection_transport_analytic_verification_2026_10_26.txt)
- [Independent exact replay and full-Z6 extension](materials/historical/analysis/c1_reflection_transport_exact_verification_2026_10_26.txt)
- [Full-Z6 sparse reflection matrix](materials/historical/analysis/c1_reflection_transport_exact_verification_2026_10_26_artifacts/full_Z6_reflection_sparse_map.json)


The original strict Z6 witness still proves the unconditional K364 Slater point and attained duality. Its use requires no universal Z6 admission theorem.
