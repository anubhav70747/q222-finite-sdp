# Q222 R10 TL3 exact source dependency

This package is a self-contained proof of the exact native C2 minus formal same-contact C1 coefficient cancellation **TL3** and of the **bounded-b premise from original R10.2**. It does not claim the TL6 conductor or full Q222 lift.

* `TL3_NATIVE_R10_PROOF.md`: pen-and-paper source dictionary, formulas, common factored numerator and coefficient table, derivations for all seven TL3 rows, exact nonphysical formal beta, and bounded-b verification.
* `verify_TL3_exact.py`: 46 exact symbolic identity checks from the original native source formulas; requires SymPy but no numerical solver.
* `verification_receipt.txt`: captured successful checker output.

Run from this directory:

```bash
python verify_TL3_exact.py
```

All Bell parents use the **F/4** normalization. The R10 source parameter `h` differs from both the native C1 letter `h_C1=1-theta` and the centered-state coefficient `h1=a1/T1`. The R10 gap `j_R` is not the SE scaling variable `j_s`, and the coefficient difference `delta L` is not a source-chart denominator.
