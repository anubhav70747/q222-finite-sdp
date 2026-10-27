# Bounded review of the Sun–Wang–Yu contact matching

27 October 2026. Reviewed [the matching note](SUN_WANG_YU_CONTACT_MATCH.md), [its checker](checks/sun-wang-yu-match/verify_match.py), and the relevant SG1–SG7, SG17–SG19 and EN1 source passages. This is a separate-agent, targeted mathematical review. The reviewer had previously derived the coefficient/gauge map and checked the proposed Hessian deduction; it is not a cold reading or a new audit of the full finite-lift proof.

**Finding:** the coefficient matching, tied-bound application, Hessian scaling and finite-contact coercivity deduction are correct under the stated canonical-source hypotheses. No mathematical correction to the note is required. This review does not establish H1, H2 or unconditional exactness of the evaluated SDP.

## Hypotheses and source matching

The argument requires d>0, all eight event weights nonnegative, quantum validity for every physical test, and an actual finite noncommuting zero. The native C1/C2 domains satisfy these coefficient restrictions, including D=0 or beta=0 in C1 and n=0 in C2. Positivity of every individual event probability is not an assumption.

Equation (1) includes the affine constant and the correct Bob coefficient b1=m+alpha-lambda-D. Its four correlator coefficients have signs (+,-,+,+), with strictness supplied by d/2 even when event weights vanish. Thus frustration and nonzero correlator coefficients follow directly from the canonical orthant. Both native families have lambda=0, and the all-positive deterministic assignment has slack 4lambda=0. Consequently their local and quantum bounds equal kappa. The strict-gap uniqueness theorem would not be an adequate native input; the note correctly uses the frustrated-case determinant and no-axis-line results.

## Gauge, determinant and full Hessian

The displayed Pauli gauge agrees with every SG2 event row. The exact identity is P=Delta M, with Delta=(1+u^2)(1+v^2), and the overlaps and Jacobian in (3) are correct. Positive finite u,v give interior overlaps and an invertible Jacobian, including u=1. The sign congruences correctly cover other nonzero finite angle signs.

Since M has order four, det(P)=Delta^4 det(M). SG3 supplies the positive three-dimensional principal block and the Z=1 state chart; no external state-simplicity theorem is needed. At contact, delta and its first derivatives vanish, so derivatives of the positive prefactor in (6) contribute nothing to the Hessian. Equation (7), including its factor Delta_c^4/det(H3_c), is exact. The state shift in (6) has invertible Jacobian and makes the full Hessian congruent to diag(2H3,D^2R). The conclusion concerns the five-variable projective chart; the homogeneous state variables retain a radial null direction.

The polynomial lemma is valid: a Hessian-null direction eliminates the quadratic term, two-sided nonnegativity eliminates the cubic, and either the negative mixed quartic gives a contradiction or the restriction is an excluded axis-parallel zero line. Its hypotheses do not require a strict local/quantum gap.

## Dependencies and permitted simplification

The proof is noncircular as written. Global validity can be assumed directly or supplied by SG7, whose balanced-Laurent argument precedes and does not use the coercivity or Hessian conclusions. SG3 and the determinant identities are direct algebraic inputs. EN1 is used to identify the source domains, not to import the later EN2–EN16 sign/coordinate conclusions that themselves use contact curvature. That distinction should remain if the alternative proof is integrated into the manuscript.

The coercivity corollary is also correct. The SG5 multiplier is positive at contact; the relevant coordinate changes are invertible. Positive definiteness of the logarithmic Laurent Hessian makes the active exponents span the plane. Their strictly positive balance puts zero in the interior of their convex hull. This proves finite-contact coercivity and strict convexity in logarithmic squared angles, hence uniqueness of the positive contact. It can replace the local support exclusions and contact-curvature derivation, but not the remaining extreme-source classification, deletion arguments or singular tensor certificates.

The real-closed-field transfer is confined to a first-order, single-field algebraic implication after positive denominators are cleared. It supplies neither a uniform positive real residue nor separated source/testing tensor squares. The note preserves this distinction and makes no fixed-degree, full-testing Gram or K364 admission claim.

## Imported results and verification limits

The external inputs are Proposition 4 and Lemma 6 of [Sun–Wang–Yu, arXiv:2610.34364v1](https://arxiv.org/html/2610.34364v1). Their hypotheses were checked against the native map. The relevant no-axis proof was inspected; this review does not claim an independent certification of the entire paper, its cited inputs, face example or membership algorithms. The new Hessian consequence is a deduction from those inputs, not a quoted theorem.

The tiny checker was rerun without `--write-receipt`: Python 3.14, SymPy 1.14.0, all 13 named symbolic checks completed successfully. In particular it checks the full SG2 matrix identity and the stated coefficient/Jacobian identities. Its determinant-coefficient check expands the imported Gram trace formula; it is not a separate derivation of every trace identity or a verification of Lemma 6. The analytic positivity and dependency arguments above remain necessary. No broad baseline tests or L20/L21 audits were performed.

No substantive correction is requested. H1/H2 remain open and the global SDP reverse inclusion remains conditional.
