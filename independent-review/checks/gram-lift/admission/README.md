# Exact pointwise admission certificates (sextic Puiseux-Gram lift)

Standalone check: for a rational source `theta` of family C1 or C2 (five chart parameters), produce and verify in exact arithmetic a
positive definite Gram certificate that the source is strictly admitted by the sixth-root (k = 6) lift.

## What is proved

For each listed source, taken at its exact rational value (the printed decimals are rounded to a rational with denominator at most 10^6,
reported as `rational_point` in the JSON):

1. the rational point lies in the native chart of its family (`verify_boundary.c1_source` / `c2_source`, exact `Fraction` arithmetic,
   all chart inequalities, source weights `>= 0`; `d = 1` by construction);
2. there is an explicit rational vector `x` (the reduced Gram blocks `G'_b`) satisfying **all 412 coefficient equations** `A_0 x = b(theta)`
   of the Gram identity `F_c(s^6, t^6, psi) = Z^T G Z` exactly: the run reports `0` rows with nonzero exact residual;
3. each of the 4 parity blocks `G'_b` (sizes 51, 35, 35, 29 after removing the contact vector from the 154 features of sizes 52, 36, 36, 30)
   is **positive definite**, proved by ball Cholesky at 400 bits (python-flint `arb`: every pivot is certified positive for the exact rational matrix).

Together these give `G = D^-1 V G' V^T D^-1 >= 0` with `F_c(s^6, t^6, psi) = Z^T G Z` as polynomials, with `G'` strictly positive on the
complement of the contact vector: sixth-root admission holds **strictly (mu > 0) at that source**. The floating-point SDP (Clarabel) is
only used to propose a candidate `x`; nothing from it is trusted, because the candidate is projected exactly onto the equations in
rational arithmetic and then re-checked from scratch.

## What is not proved

* Nothing about neighbouring sources: no open neighbourhood, no uniform margin, no statement for any other point of the chart,
  its walls or its corners. Each certificate is for one rational point.
* The float margins `t`, `trel` are informational (they are not part of the proof). They are small (`trel` as low as 9e-11) at four of the five points,
  so those certificates rest on the exact check, not on a comfortable float margin.
* The negative control is a sanity check on the pipeline, not a proof of inadmissibility: failing to certify one source says nothing
  about whether that source is admitted.
* The certificate is relative to the shared inputs, which are used but not re-derived here: the native source formulas
  (`../boundary/verify_boundary.py`), the coefficient tensor of the parent form (`../gram_lift.py`, `tensors()`), and the Newton-pruned
  parity-block feature basis (`gramtools.make_blocks`).

## Run

Dependencies: Python 3.9 or newer, `numpy`, `scipy`, `sympy`, `mpmath` (imported by `../gram_lift.py`), `clarabel`, and **`python-flint`**
(`pip install python-flint`; it supplies `arb`, `fmpq`, `fmpq_mat` for the exact solve and the ball Cholesky). Tested with
Python 3.9.6, numpy 2.0.2, scipy 1.13.1, sympy 1.14.0, mpmath 1.3.0, clarabel 0.11.1, python-flint 0.6.0.
No environment variables, no absolute paths, no caches: all imports are resolved relative to this folder, which must stay at
`independent-review/checks/gram-lift/admission/` so that `../gram_lift.py` and `../boundary/verify_boundary.py` are found.

```
cd independent-review/checks/gram-lift/admission
python -I run_points.py                 # the five sources + negative control; writes admission_points.json; exit status 0 iff all as expected
python -I run_points.py --quiet         # same, without the per-point trace
python -I exact_primal.py C1 0.111 0.0967 0.0543 0.114 1.0                    # one source
python -I exact_primal.py C2 1.0 0.5 0.65 0.4 0.6 --allow-outside             # the negative control alone
```

Without `--allow-outside` a source whose rationalisation is outside the native chart is refused immediately; the flag lets the
pipeline run on it anyway (used only for the negative control). One source takes about 0.5 to 1.5 minutes (single-source time is
dominated by the exact rational projection and varies with machine load and BLAS threads); the full run takes a few minutes.
The float stage is deterministic on a given machine but can shift in its last digits with the BLAS library and thread count; the exact
checks decide the outcome, and a bad float candidate cannot produce a false `True`.

Files: `exact_primal.py` (the certificate), `run_points.py` (driver), `hcore.py` (native charts, right-hand side, float SDP),
`rs_core.py` (coefficient matrix `A_0`), `centre.py`, `pcsolve.py` (float helpers), `exact_check.py` (exact right-hand side),
`nd_core.py` (exact coefficient tensor, chart check), `gramtools.py` (feature basis).

## Results

Columns: `t`, `trel` = float margin of the SDP stage and its relative version; `res rows` = rows of the 412 with nonzero exact residual;
`exact PD` = all four blocks positive definite and `res rows = 0`. (Full data: `admission_points.json`.)

| source (family, decimals as given) | rational point | t | trel | solver | res rows | exact PD per block | exact PD | time |
|---|---|---|---|---|---|---|---|---|
| C1 0.111 0.0967 0.0543 0.114 1.0 | 111/1000, 967/10000, 543/10000, 57/500, 1 | 6.970e-10 | 1.754e-10 | AlmostSolved | 0 | T T T T | **True** | 85 s |
| C1 0.3645 0.01901 0.01046 0.8042 1.086 | 729/2000, 1901/100000, 523/50000, 4021/5000, 543/500 | 3.167e-10 | 8.957e-11 | AlmostSolved | 0 | T T T T | **True** | 64 s |
| C1 0.0945 0.146 0.0544 0.335 0.999 | 189/2000, 73/500, 34/625, 67/200, 999/1000 | 2.714e-09 | 1.022e-09 | Solved | 0 | T T T T | **True** | 56 s |
| C2 0.02828 0.0009625 0.001251 0.47 0.9988 | 707/25000, 77/80000, 1251/1000000, 47/100, 2497/2500 | 6.932e-09 | 4.032e-09 | AlmostSolved | 0 | T T T T | **True** | 62 s |
| C2 1 0.5 0.65 0.4 0.3 | 1, 1/2, 13/20, 2/5, 3/10 | 3.171e-02 | 1.277e-02 | Solved | 0 | T T T T | **True** | 50 s |
| C2 1.0 0.5 0.65 0.4 0.6 (control, n = -83/846 < 0) | 1, 1/2, 13/20, 2/5, 3/5 | refused | | | | | False | 0 s |
| same, with `--allow-outside` | same | -1.441e-01 | -6.63e-03 | Solved | 0 | F F F T | **False** | 30 s |

The control satisfies the chart inequalities but its source weight `n = -83/846` is negative, so it is not a chart source: it is refused by
default, and with `--allow-outside` the exact projection still satisfies all 412 equations but three of the four blocks are not positive
definite, so it does not certify.
Total wall time of the run above: 347 s (shared machine; an identical earlier run took 200 s with the same digits).
