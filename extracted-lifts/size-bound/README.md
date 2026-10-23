# Necessary PSD order for the full Q222 set

Every exact lift of Q222 through a single cone of positive semidefinite
matrices has order at least **9**. For order `N`, redundant auxiliary
directions compress to at most `N(N+1)/2 - 9`.

The complete proof is the proposition in
[`../../latex/Explicit_Lifts.tex`](../../latex/Explicit_Lifts.tex),
Section “Size bounds and exact algebraic search”. It constructs nine
nonempty exposed faces from zero-event conditions, proves that their
inverse images have strictly decreasing PSD supports, and uses compactness
to keep the final support nonzero. It does not assume strict feasibility.

Run the exact standard-library checker from this directory:

```text
python verify_size_bound.py
```

It writes `face_chain_data.json`, enumerating all 65,536 subsets of the
16 zero-event conditions on the 16 deterministic assignments. There are
1,721 distinct nonempty survivor sets. The longest strict nonempty chain
has nine faces. The data list the eight added events, removed witnesses,
remaining assignments, and exact rational affine dimensions.

This bound excludes the stored three-block order-two sentence in
`../conditional-global-search/Q222_SOC_M3_UNSOLVED_INPUT.json` and its
Wolfram Language counterpart: their total order is six. Their historical
input bytes and filenames are preserved; their mathematical status is
**EXCLUDED**, as recorded in `face_chain_data.json`.

The first unexcluded single-cone synthesis template has order 9, at most
36 auxiliary variables, and 2,025 unknown matrix coefficients. This is
a necessary size restriction, not an existence claim at order 9. For
several PSD blocks the bound applies to the **sum** of their orders;
it does not bound the largest block or rule out second-order-cone lifts
with more blocks.
