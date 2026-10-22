# One certificate for the reciprocal Alice angle

The proof completes the leading three-state form and rescales its two
source-dependent directions. The remaining coupling vector has squared norm
at most one for every positive source ratio. Three small correction matrices
then account for the complete parent, with all other event squares retained.

`inverse_alice_all_z_certificate.pdf` gives the derivation and source bounds.
The native source domain is

```
d=1, z>0, 0<x,sigma<=2^-78/(1+z),
0<b,h<1, bz<1, the printed strict gaps, D,beta>=0.
```

The independent tests satisfy `0<=p=1/u<=2^24`, `0<=w<=2^-36`;
`p=0` includes the projective endpoint. Fixed sign gauges give signed
angles, and the state coordinates are unrestricted. The ratio `h/sigma`
is unrestricted. The explicit source matrices act on 51 fixed testing
functions.

`SOURCE_PENCIL.coo.json` is an exact rational pencil for this specified
weighted certificate cone: 19 blocks, total order 62, largest block 6,
and 97 auxiliary variables. Its seven visible coordinates are
`(d,m,C,A,D,alpha,beta)`. They are canonical Bell-source coefficients.
The file stores actual symmetric matrix entries; mirror each off-diagonal
entry once.

`ALL_ALICE_ORDER731.coo.json` represents the intersection with the existing
inner-angle certificate cone. It has 139 blocks, total order 731 and 2791
auxiliaries. For `0<z<=1` on the above source domain, its two constructions
cover every projective Alice angle and `|w|<=2^-36`.

Both pencils represent their specified certificate cones exactly. Extending
these certificates to all sources and Bob angles remains necessary for a
full eight-coordinate quantum-set pencil.

## Reproduce

Python 3 and SymPy are sufficient. From this directory run:

```text
python inverse_alice_all_z_blocks.py
python check_inverse_alice_all_z_completion.py
python compile_weighted_pencil.py inverse_alice_all_z_blocks.json .
python verify_weighted_pencil.py inverse_alice_all_z_blocks.json .
python combine_alice_pencils.py ../local-order669/root/C1_UNIFORM_LOCAL_ORDER669.coo.json SOURCE_PENCIL.coo.json ALL_ALICE_ORDER731.coo.json
```

The independent verifier rebuilds the complete canonical polynomial and every
weighted feature product. It verifies all 80 coefficient equations, rank 79
and all 97 free Gram directions. Native source inequalities are proved in
the accompanying text. No numerical optimizer is used.
