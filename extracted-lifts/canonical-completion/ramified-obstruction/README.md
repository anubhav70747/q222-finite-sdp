# Obstruction to the ramified order-49 completion

The all-positive valid source

```text
(d,lambda,n,m,C,A,D,alpha,beta)
= (1,1/576,161,2^-24,2^-34,2^-14,1/2,1/250,76)
```

lies outside the ramified completion cone. Four positive-definite integer
Hankel matrices give a separating trace identity: all 148 auxiliary
coefficients vanish, the nine source coefficients equal `12*ell`, and

```text
ell(cdagger) = -32607070960447356138172589/9663676416000.
```

The same source excludes **every polynomial selector in `s=sqrt(v)/4`**:
positive semidefiniteness forces its degree to be at most four. Rational
selectors, finer ramifications, and enlarged testing spaces remain
separate questions. Combined with the polar assembly theorem, the
obstruction makes the eight-copy order-392 outer representation strict.

Read `Q222_Ramified_Quartic_Obstruction.pdf` or
`proof/ramified_obstruction.tex`. The appendix prints all ten integer
four-by-four moments. `data/certificate.json`, both files under `reused/`,
and `verify.py` retain the received bytes.

## Exact checks

From this directory, using Python 3 and only its standard library:

```text
python verify.py
python check_independent.py
python check_symmetry_transport.py
```

The first verifier checks the homogeneous physical transport, the integer
Hankel matrices, all source coefficients, auxiliary cancellations, and
the transported contact. The independent checker uses fresh rational
Schur elimination and contracts the certificate against the sibling
published pencil at `../ramified-completion/data/ramified_order49_pencil.json`.
It checks all nine scalar source-sign blocks as well. Its exact results
are in `independent_checks.json`.
The symmetry checker independently expands all nine homogeneous source
packets, the transported contact and all six first derivatives, and the
fractional feature exponents. It writes `symmetry_transport_checks.json`.

The reused pencil has SHA256
`5428ca18e35f2d2468c2560eba3109cf26b6c98edc1dcca9c7a108aef80d7421`.
The source's physical validity and the admission of its preimage `cbar`
are proved in the referenced canonical-completion notes and matrices.
The source transport is a homogeneous polynomial identity, including
zero angles and projective endpoints.

To compile the note, run `pdflatex` twice on `ramified_obstruction.tex`
from `proof/`. `SHA256SUMS.txt` lists the hashes of the files in this folder.
