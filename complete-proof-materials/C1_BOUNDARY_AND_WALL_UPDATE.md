# C1 endpoint construction and exact boundary certificates

The C1 constructions now control the full state matrix near both Bob endpoints and in a two-sided neighborhood of its contact. Joining these regions with one globally positive polynomial remains open. The merged boundary calculations add exact Gram witnesses and a smaller rational family on the C2 n=0 wall. These results sharpen the two remaining admission problems; H1 and H2 remain open.

## C1: one polynomial matrix, two positive ends and a contact anchor

For each genuine strict oriented C1 source, the complete parent can be written as F=psi^T(P0+uP1+u^2P2)psi. A source-only quartic skew matrix V(v) leaves the identity

    (psi,u psi)^T H(v)(psi,u psi)=(1+v^2)F

unchanged, where H has diagonal blocks (1+v^2)P0 and (1+v^2)P2 and cross block (1+v^2)P1/2+V. The construction chooses its endpoint coefficients and first derivative explicitly. Its two-dimensional soft Schur matrix at v=0 is diagonal, with entries

    4m(2beta+d)/(2beta+d+2m),
    4C(2D+d)/(2C+2D+d).

Both are positive on the strict source domain. The complementary six-dimensional block is positive, and the leading matrix at v=infinity is positive definite. Thus the same full 8-by-8 polynomial matrix is positive definite for all sufficiently small positive v and all sufficiently large v, with source-dependent cutoffs.

At the physical contact, Calderon's nonnegative 2-by-n biquadratic theorem supplies a PSD skew completion. A source-only algebraic interpolation matches this contact matrix while retaining both endpoint results. The explicit rational kernel completion alone is not claimed PSD. Physicality of the genuine native source, including the exact Alice sign congruence, is essential to the contact argument.

If this H(v) can be made PSD for every v>=0, Hanselka--Sinn Corollary 2.8 and the degree bound give a Gram on 40 fixed features: 24 integer-power rows and 16 half-integer Bob rows. All belong to K364. The degree bound follows because leading terms in a diagonal sum of real polynomial squares cannot cancel; even/odd averaging then gives the half-line factorization, including singular matrices. Positivity on the intervening v interval is unproved. This 40-feature route is sufficient for H1, not necessary.

- [Original C1 proof](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1959/task5/package/PROOF.md)
- [Original complete C1 package](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1959/task5/original-package.zip)
- [Independent analytical verification](materials/historical/analysis/c1_bob_harmonic_analytic_verification_2026_10_26.txt)
- [Exact replay and additive sign checks](materials/historical/analysis/c1_bob_harmonic_exact_verification_2026_10_26.txt)
- [Primary-source import check](materials/historical/analysis/c1_bob_harmonic_external_imports_2026_10_26.txt)

## C1: a nondegenerate contact and two-sided local positivity

For every genuine real strict oriented C1 source, including k=1, the fixed-Bob form has exactly two projective zeros related by the Alice sign congruence. The strict native Hessian and compactness at that one source allow a small multiple of the contact-annihilator square to be subtracted. Calderon's theorem then supplies a rank-six contact Gram, positive on the complement of the two contact directions.

Its cross-parity skew lies on an explicit source-rational affine line. The strict complement inequalities admit a source-algebraic choice of its one free scalar. A positive 2-by-2 source-only solve adjusts the first Bob derivative so that the full Schur minimum agrees with the physical tangent minimum. A quadratic skew correction divides the positive second-order curvature between the two soft modes. This gives one quadratic skew polynomial whose full 8-by-8 matrix is PSD at contact and PD throughout a two-sided punctured neighborhood.

The neighborhood radius depends on the source. The constructor requires a free scalar satisfying the strict contact inequalities; it does not implement the analytical selection argument. The endpoint polynomial and this contact polynomial have not been joined within one globally PSD degree-four construction. The remaining sufficient target is one source-only polynomial making both complete 2-by-2 Schur residuals PSD for every positive Bob angle, or a direct full K364 Gram. H1 remains open.

- [Original C1 contact-jet proof](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-20261027-0705/task5/package/Q222_C1_K364_CONTACT_JET_2026_10_27/PROOF.md)
- [Original complete contact-jet package](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-20261027-0705/task5/original-package.zip)
- [Independent contact-local theorem review](materials/historical/analysis/c1_contact_jet_analytic_verification_2026_10_27.txt)
- [Targeted exact replay and independent moment checks](materials/historical/analysis/contact_jet_and_dual_exact_verification_2026_10_27.txt)

## Boundary identities, finite witnesses and the generic C2 wall family

The C2 r->0 three-state limit has no integer-power polynomial SOS anywhere on its stated physical face: a forced positive diagonal contradicts the contact kernel. Half-power PSD certificates are exact at nine stored rational face points. They do not establish half-power sufficiency over the whole face or control the full fourth-state normal when r>0.

The graph substitutions for C2 n=0 and C1 beta=0 or D=0 are exact on their actual native domains. The C1 vanishing weights force the corresponding uvZ or vZ Gram row to vanish. Exact half-power PSD certificates cover eight C2 n=0 points, nine C1 beta=0 points and eight C1 D=0 points on the larger supports.

The smaller 18-feature C2 n=0 family has a universal coefficient identity with 28 free entries. Its contact equations have rank 12 over the rational function field, leaving a generic 16-parameter affine family. This does not establish rank or consistency at every physical specialization: exceptional denominators or vanishing pivot minors remain to be controlled. Exact PSD witnesses are supplied at 22 rational wall points. A contact direction occurs in each parity block, so the total nullity is three on this support and four on the larger four-block supports.

A source-algebraic PSD selection on the dense physical open set, followed by the established coefficient-space closedness at exceptional limits, would suffice to prove whole-wall admission. Such a selection has not been proved. This wall also does not replace the strict C2 target: opening n>0 restores the forced uX diagonal and requires a coupled moving-contact completion. Adding positive event squares to a fixed wall Gram alone cannot supply that completion.

- [Independent review of upstream 3b1450d](materials/historical/analysis/upstream_3b1450d_analytic_review_2026_10_26.txt)
- [Four exact replays and witness counts](materials/historical/analysis/upstream_3b1450d_exact_verification_2026_10_26.txt)
- [Original upstream checkers and data](../independent-review/checks/gram-lift/ends/)

The four checkers certify 56 stored witness instances, including repeated sources on different supports, and 202 reduced-block positivity checks. These are finite exact certificates. Universal graph and coefficient identities have broader scope; numerical samples do not extend the PSD conclusions.

## Remaining implication

Construct a source-only PSD completion throughout each actual strict native domain, either by solving the C1 polynomial-matrix problem and a coupled C2 completion or directly in K364. Once H1/H2 hold, the verified density, reflection and full-group assembly results apply to the already evaluated conditional global pencil. L13/L14 receive this evidence without a new completion milestone; L08-L10 and H1/H2 remain open.
