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

`contact-selector/` supplies the repaired C1 support in the actual enlarged
source pencil, its four rational PSD factorizations, and a new rational
selector on an explicit open five-parameter contact region. The compressed
Gram orders are 13, 14, 13 and 14. The proof restores the complete Bell
operator and all quantum realizations. The same fixed-free selector fails
outside the certified region; universal adequacy of order 13744 remains open.
The seven JSON files contain the source arrays, inverse coefficient map,
open-domain selector and rational one-parameter family. Reproduce with:

```text
python contact-selector/check_c1_contact_selector.py
```

The script uses SymPy and exact rational arithmetic. It checks the complete
parent, all state and angular contact equations, the inverse-minor identity,
strict Gram margins, restored Bell operator and the displayed selector failures.

`fixed-features/fixed_box_c1.pdf` derives complete local C1 certificates
on `0<x,sigma<=2^-160`, `0<z<=1`, the original source gaps and
`D,beta>=0`, with independent testing angles of modulus at most `2^-100`.
The ratio `h/sigma` is unrestricted. Positive angular coefficient uses
220 labelled testing functions; nonpositive coefficient uses 157. Their
union therefore contains at most 377 functions. Both complete parent
identities include every retained square and state-normal term. The
construction treats all real state signs, testing axes and the zero event
faces; `D=0` occurs only in the positive branch on this source box.

`c1_positive_lambda_local220.json` contains 220 labelled functions,
144 fixed directions, the main Gram and explicit prescriptions for at
most 2196 square roots. `fixed_box_c1_nonpositive_data.json` contains
157 functions, 36 directions, the full six-root Gram and all 14 retained
squares. Source formulas and the separated scalar factorizations are
included. These are evaluated local function lists and symbolic source
coefficients, not the matrices of a full-set affine pencil.

The three portable SymPy checkers reconstruct the displayed polynomial
identities, source-coordinate maps, exact coefficient norms and rational
margins. Run:

```text
python fixed-features/fixed_box_c1_checks.py
python fixed-features/c1_retained_union_checks.py
python fixed-features/fixed_box_c1_nonpositive_checks.py
```

The first checker also retains the larger 348-function main-Gram count
obtained by padding the monomial support; the complete local construction
uses its reduced 180-function main list. The tensor-rank note explains
the further common-cover and product-coefficient steps required globally.

The C160 event note derives all sixteen probability certificates directly
from two-qubit amplitudes. Its JSON contains 64 integer root vectors in the
existing 32-dimensional block. These prove probability positivity, normalization
and no-signalling for the candidate. Reproduce with:

```text
python outer-candidates/check_c160_event_grams.py
```

`local-order669/` contains an evaluated rational pencil for the seven-coordinate
local source cone `(K_POS intersect K_NEG) + K_SIX`. Its 120 blocks have total
order 669 and maximum order 12, with 2694 auxiliary coordinates. The two state
charts and six-root construction use explicit fixed algebraic testing
functions. Both inclusions for this specified certificate cone follow from
exact coefficient elimination. Source admission is proved on
`0<x,sigma<=2^-30`, `0<z<=1`, with the stated source gaps and `D,beta>=0`.
The independent testing box is `|u|<=2^-24`, `|w|<=2^-12`; state coordinates
are unrestricted. These are local source coordinates, not the eight quantum
correlation coordinates. Run from `local-order669/`:

```text
python root/verify_final_pencils.py
```

The checker rebuilds the full-parent coefficient maps with exact fractions,
verifies ranks 281, 221 and 180, all 2687 component free directions, and all
682 equations of the final combined pencil. The corrected 17-page proof,
matrix data, product maps and symbolic regeneration programs are included.


`inverse-alice/` contains a seven-page completion of squares for every
positive source ratio, with 51 fixed testing functions. The actual local
source-cone pencil has total order 62 (19 blocks, 97 auxiliaries); its
intersection with the inner-angle cone has order 731 (139 blocks, 2791
auxiliaries). On the stated bounded-z source domain, the latter combines
certificates for every Alice angle and `|w|<=2^-36`. The directory README
gives the exact domains, matrix-entry convention and reproduction commands.
Both pencils act on seven Bell-source coordinates. The full quantum-set
pencil still requires the other source and Bob-angle regions.

`canonical-completion/` eliminates the Alice angle exactly. For every fixed
Bob angle, the resulting canonical validity condition is equivalent to two
four-by-four PSD inequalities with four free completion entries. The exported
integer arrays retain all nine source coefficients. The angle quantifier
remains: these are pointwise matrices. A rational valid source rules out
every globally polynomial completion in this form. Run from that directory:

```text
python check_selector_certificate.py
python check_source_validity.py
python check_simple_source.py
python check_alice_quadratic_reduction.py
```

The first two checks use exact standard-library fractions; the last two use SymPy.
The PDF derives the two-matrix characterization and the obstruction.

`size-bound/` gives a probability-face chain proving total PSD order at least
nine for any full-set lift, with its exact Boolean checker. It excludes the
historical three-block order-two input. The bound concerns total order,
not the largest block in a product of PSD cones.
