Reviewed pointwise K364 factor verifier

Run python verify_primal_factor_checked.py --selftest for the elementary n-event check. Run python verify_primal_factor_checked.py candidate.json for a proposed exact certificate. The compact source map matches all 120 entries of the previously evaluated L18 map.

This variant preserves the original coefficient and native-source calculations and adds checks for provably finite real coefficients, nonnegative direct-source weights and known source names. It does not generate a Gram or prove universal H2. Inputs are trusted exact SymPy mathematical expressions.

false_positive_example.json records a negative n-event polynomial accepted as an SOS by the original program because it allowed imaginary roots. The reviewed variant rejects it; it also independently rejects a single I factor. The original program and downloaded archive are preserved separately under fetched-20261027-0745/task3/.

review_results.json records the targeted replay, source-map comparison, archive/manifest integrity and corrected checks. No universal native Gram is included.
