# Explicit scientific data: exact projection, local certificates and outer candidates

**No file in this directory is an evaluated exact lift of the full eight-coordinate Q222 set.**

`tilted-projection/` gives an exact lift of

```
t = a0,  s = E00+E01+E10-E11,
|t| <= 1,  t^2+s^2/4 <= 2.
```

The matrix orders are 5 with one auxiliary, or 4 with no auxiliary. Coefficients are exact pairs in Q(sqrt(2)); sqrt(2) is positive. `TILTED_SCOPE_COUNTEREXAMPLE.json` is an exact nonquantum behavior passing the projection and every relabelled single-marginal tilted-CHSH inequality. Run:

```
python tilted-projection/build_and_check_tilted_extraction.py
```

`regular-c2-local/` contains the full parent certificate on the explicit normalized source/contact box and independent testing box of radius `510101005/350806159660744`. It has 785 PSD(2) blocks, 1570 listed functions, 237 coefficient equations and 2355 Gram variables. The complete parent, coefficient extraction and sign proof are in the accompanying `Explicit_Lifts.tex`. Run:

```
python regular-c2-local/regular_c2_check.py
```

`c1d-local/` gives a different source-derived LOCAL library, radius `1/2788018`: 82928 separated squares, 1212 testing radicals and 3844 testing spanning functions. Original source fixtures, polynomial allocation multipliers, full five-root Gram, defining algebra and restoration formulas are included. It asserts no global cover. SymPy is required for:

```
python c1d-local/code/verify_source_allocations.py
python c1d-local/code/verify_local_square_library.py
```

`outer-candidates/` contains four full-coordinate OUTER constructions:

| Total PSD order | Auxiliaries | Status |
|---:|---:|---|
| 8368 | 16480 | Reverse inclusion disproved: strictly larger than Q222 |
| 10736 | 21472 | Reverse inclusion disproved: strictly larger than Q222 |
| 13744 | 29920 | Forward inclusion proved; reverse inclusion unproved |
| 63888 (C160) | 31111 | Reverse inclusion unproved; full forward inclusion uses the stated extremal-coverage theorem |

For the first three, `[k,i,j,value]` is an upper-triangular matrix entry: k=0 constant, k=1..8 the visible coordinates `(a0,a1,b0,b1,E00,E01,E10,E11)`, k=9+j auxiliary j. Mirror each off-diagonal exactly once. Blocks use global offsets. For C160, `[block,i,j,k,numerator,denominator]` uses local block indices and the same coefficient indexing. All coefficient arrays are exact rationals. The status metadata of the two smaller candidates has been corrected to reflect their later strictness theorem; the matrix and block arrays are unchanged.

The exact physical source, source maps and eight positive dual separator matrices establish why the two smaller candidates fail. They provide a separation proof of strictness, rather than one printed nonquantum behavior with all auxiliary coordinates. They do not refute the larger two constructions. Run the standard-library checks:

```
python verify_outer_data.py
```

This checks every matrix record's indexing/rational encoding and all coefficients of both exact trace-separator identities, including every zero auxiliary coefficient and all positive LDL pivots. It does not solve an SDP or prove the unproved reverse inclusions. `C160_construction.tex` gives the canonical moment construction and identifies its required Barizien–Bancal extremal-coverage input. `compile_c160_candidate.py` uses the included source cache and writes to `generated_c160/`; its two filesystem paths were relocated, with algebra unchanged. Running it is optional and is not needed to read the evaluated matrices.

`conditional-global-search/` provides the exact rational nine-atom membership formula, real-algebraic coefficient-descent argument and an UNSOLVED three-order-two-block quantifier-elimination input. The candidate has 162 coefficient unknowns and nine auxiliaries. It supplies no solved coefficient tuple. General finite-SDP existence does not imply feasibility at that order-two pattern.

`PROVENANCE.json` records the scientific input members and hashes. No chat exports, private agent notes, numerical sweeps or release histories are included. Exact reconstruction checks support their stated identities and domains; they are not substitutes for either missing full-set inclusion.
