# Q222 fixed-cone density entrance (26 October 2026)

Files:
- `COORDINATOR_SHORT_THEOREM.md`: concise conditional density and finite-cone assembly theorem
- `FIXED_CONE_SINGULAR_SOURCE_DENSITY_THEOREM.md`: complete analytic argument, exact dependencies, native C1/C2 entry, field-transfer and sum-closedness proof
- `verify_native_entrances.py`: deterministic SymPy symbolic checker for C1/C2 native inverse and gap maps, and 16-point sign orthogonality (no random tests)
- `EXACT_CHECK_OUTPUT.txt`: output from successfully running this checker
- `INPUT_07_FIXED_CONE_DENSITY_sources.txt`: unchanged user-provided frozen source material
- `density_entries.json`: machine-readable map of 3 literal boundary entrances and surviving source conditions

Execute: `python verify_native_entrances.py` (Python 3, SymPy).

The IFT and cone-closedness arguments are analytical in the proof, not claims based on the checker. No strict-source PSD claim is established here. No GitHub write was made.
