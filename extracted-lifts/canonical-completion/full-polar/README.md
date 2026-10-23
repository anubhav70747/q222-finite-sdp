# From canonical certificates to all eight correlations

[The short proof](Q222_Polar_Assembly.pdf) shows that eight CHSH output
relabelings recover the full quantum polar. For any closed sound canonical
cone containing all positive event terms, the assembled correlation set
equals Q222 if and only if the canonical cone is exact. Other orientations
cannot repair a missing canonical certificate.

The proof applies to any canonical SDP. A strictly feasible canonical pencil
gives eight adjoint copies with both full-set inclusions. The supplied
order-33 specialization yields an outer approximation of total order 264
with 384 auxiliaries. The exact counterexample in `../ramified-completion/`, together with
the assembly theorem, proves that this outer approximation strictly contains
Q222.

## Exact data

- `full_polar_maps.json`: the canonical map and inverse, eight relabelings,
  deterministic neighbours, all 28 local PR-midpoint decompositions and
  a rational 92-by-48 kernel eliminating the adjoint equations.
- `slater_witness.json`: a literal rational strict feasible point, all 13
  positive blocks, and their exact leading principal minors.
- `slater_proof.txt`: the coordinate reconstruction proving explicit
  positive lower bounds without a numerical solver.

Run from this directory:

```text
python build_full_polar_maps.py --check
python slater_checker.py
```

The map checker uses SymPy; the Slater checker uses only the Python standard
library. Both read the published order-33 pencil in `../halfline/`.
The map has determinant -32; the full and auxiliary adjoints have ranks 53
and 44. The eight-copy formula has no remaining coefficient equations.
The data check their stated finite identities; universal canonical coverage
is a separate mathematical requirement.
