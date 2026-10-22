# An explicit local C1 pencil of order 669

The seven-coordinate source cone
`K669 = (K_POS intersect K_NEG) + K_SIX` has the exact rational pencil
in [root/C1_UNIFORM_LOCAL_ORDER669.coo.json](root/C1_UNIFORM_LOCAL_ORDER669.coo.json).
It has 120 PSD blocks, maximum block order 12, and 2694 auxiliaries.
Both inclusions refer to this specified local certificate cone.

The [complete derivation](proof/Fixed_Family_Extraction.pdf) gives the
fixed algebraic functions and proves admission of the native stationary
sources on `0<x,sigma<=2^-30`, `0<z<=1`, with the printed gap conditions.
The independent testing box is `|u|<=2^-24`, `|w|<=2^-12`; all state
coordinates and both state signs are included. The ratio `h/sigma` is
unrestricted. The full eight-coordinate quantum lift still requires the
other source and testing regions.

From this directory, run:

```text
python root/verify_final_pencils.py
```

This standard-library checker rebuilds the product maps and all 682
coefficient equations with exact rational arithmetic. All 2687 component
free directions and the seven source-splitting coordinates are checked.
[README.txt](README.txt) gives the full regeneration sequence.

One intermediate row bound in the recovered proof has been corrected:
the `w^2` row is bounded by 5, rather than 4. Its explicit estimate is now
printed next to the bound. The uniform bound 128, final positive reserve,
feature functions and all matrix entries are unchanged. The supplied PDF
is rebuilt from the corrected source. Input and publication hashes are
recorded in `PROVENANCE.json`; exact checks and analytic inequalities have
distinct roles.
