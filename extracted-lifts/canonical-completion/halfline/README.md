# Quadratic completion on a half-line

The [proof](Q222_Halfline_Completion.pdf) gives an exact affine representation
of the sources admitting a quadratic completion on v >= 0. Its four PSD
blocks have orders 8, 8, 4, 4, with 44 auxiliaries. Nine scalar source-sign
blocks give total order 33. All matrix coefficients are integers.

Every represented source is quantum valid. Whether every valid canonical
source belongs to this cone remains open.

The auxiliaries satisfy |y_i| <= 16 sum(c), so the projected cone is closed.
The source that excludes all global polynomial completions has the explicit
piecewise quadratic completion printed in the proof. Its four rational Gram
matrices have ranks 5, 7, 1, 3. The JSON contains literal matrices and LDL
factors; no numerical solver is needed to verify them.

Run with Python's standard library:

```text
python halfline_gram.py
python check_simple_source_exact_witness.py
python check_strict_trace_certificate.py
```

The first command regenerates both pencils and checks their coefficient
identities. The second replays the fixed rational witness independently.
COO records contain both off-diagonal copies; do not mirror or double them.
Compile `halfline_combined.tex` twice with pdfLaTeX.
