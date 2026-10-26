# Q222 C2 fixed Delta angular library: exact audit

* `FIXED_DELTA_C2_PROOF_STATUS.md` - fully self-contained mathematical proofs of soundness, fixed-cone closedness, exact Newton support, and a sharp weighted-Motzkin non-SOS lemma. **No admission/exclusion proof for full native C2.**
* `check_fixed_delta.py` - bounded exact checker using `fractions.Fraction` and SymPy; no solver or sweep.
* `verification.json` - actual results of the checker.

Run: `python check_fixed_delta.py`

The exact native formulas are reconstructed in the checker, including all four state coordinates and both independent testing angles. No archive extraction is needed.
