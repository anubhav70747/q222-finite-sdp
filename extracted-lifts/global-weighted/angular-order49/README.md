# The angular determinant cone

The [three-page proof](angular_cone.pdf) gives an exact rational affine lift of

`L = {(h,q): q_e >= 0, h + sum q_e U^e1 V^e2 >= 0 for all U,V > 0}`,

where `e` runs through the eight nonzero points of `{-1,0,1}^2`.
Twelve balanced circuits give 24 order-two PSD blocks and one scalar block.
The pencil has total order **49**, **9 inputs**, **48 auxiliaries** and
**109 nonzero rational coefficient records**. It has no separate coefficient
equations. Its diagonal entries enforce every allocation's nonnegativity.

The JSON file `ANGULAR_CONE_ORDER49.coo.json` contains the complete evaluated
pencil. Its entries are `[row, column, variable_index, numerator, denominator]`,
with zero-based indices. Mirror each off-diagonal entry without a factor of
one half. The variable names and PSD-block boundaries are included.

Run with standard Python:

```text
python check_angular_pencil.py
```

This enumerates all twelve primitive circuits, reconstructs every input and
auxiliary coefficient matrix independently from the printed formulas, and
checks 18 exact feasible witnesses. The witnesses include zero weights,
proper faces with an unattained infimum, signed thresholds and simultaneous
active circuits. The proof establishes the universal cone equality.

To rebuild the data:

```text
python compile_angular_pencil.py
python check_angular_pencil.py
```

For a specified nonnegative canonical Bell source, the existing determinant
coefficients give the nine inputs. For example:

```text
python source_to_angular.py 1 1 1 1 1 1 1 1 1
```

The arguments are `d lambda n m C A D alpha beta`; fractions such as `1/2`
are accepted. The script uses the previously verified coefficient table in
`../weighted/DETERMINANT_SOURCE_COEFFICIENTS.json` and exact rational arithmetic.
It returns the pencil inputs, not a numerical feasibility verdict.

The resulting test decides canonical Bell validity exactly on the stated
nonnegative source orthant. Its input map is quartic in the Bell coefficients.
This is an affine lift of the Laurent coefficient cone, **not an affine lift
of Q222 or of the Bell-source cone**. The order-133 construction remains the
explicit fixed square certificate for the determinant-weighted full state
polynomial.
