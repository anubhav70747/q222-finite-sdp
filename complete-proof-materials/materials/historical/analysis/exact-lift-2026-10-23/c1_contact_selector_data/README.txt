EXACT C1 CONTACT CERTIFICATES

Scope: the canonical C1 coefficient cone used by the order-13744 candidate.
This is not an exact full-Q222 lift and not a universal C1 selector.

Original rational scientific data, retained byte for byte:

1. C1_ENLARGED_PRIMAL.coo.json: actual order-59 source certificate pencil,
   four Gram blocks (14,15,14,15), one [d] block, 302 auxiliaries.
2. C1_REPAIRED_TWO_FRAME_COUNTEREXAMPLE.json: the source c*, its complete
   auxiliary vector, four Gram matrices and rational positive factorizations.
3. C1_REPAIRED_CHAMBER_AUXILIARY_MAP.json: exact 302-by-7 linear selector
   for d>=0 and c_i>=d(c*)_i, plus six event-ray certificates.
4. C1_ENLARGED_REGULAR_CONTACT_SYSTEM.json: all constant coefficient maps,
   features, 148 monomials, 392 compressed Gram entries, 142 pivot columns,
   inverse pivot matrix, and 250 free columns.
5. C1_REGULAR_CONTACT_DOMAIN.json: the full 6-by-7 contact matrix from the
   quarter-parent, its regular source chart and additional source data.

New exact data:

6. RATIONAL_OPEN_CONTACT_SELECTOR.json: compressed reference Grams, fixed
   free entries and target, exact inverse traces, quantitative domain and
   restored eight-coordinate support/contact values. Formula:
      h_free=h*_free,
      h_pivot=h*_pivot+P (b(theta)_I-b*_I).
   It is proved for positive native regular source contacts satisfying
      max |b(theta)-b*| < 1/3284140032.
   The four reference compressed Grams have margin at least 2^(-22).
7. ONE_PARAMETER_SOURCE_FAMILY.json: six short rational coefficients for
   contact (1/13,16,A,-1/32,1/50). The selector is proved for
      |A-140|<=2^(-42).
   Larger exact values A=141 and A=70 show this particular selector can
   be indefinite. Those failures do not prove failure of the source cone.

The homogeneous source vector has final coordinate one and is not a unit
state. The stored quantum behavior divides by its squared Euclidean norm.

Run ../check_c1_contact_selector.py with Python and SymPy. The copied data
alone suffice; the historical archive is optional. When a sibling
../outer-candidates/C1_enlarged_certificate_source_pencil.coo.json is
available, the checker compares its literal bytes as well.
