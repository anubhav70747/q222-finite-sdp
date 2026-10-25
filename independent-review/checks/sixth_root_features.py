"""Exact scalar identities and obstruction to the thirds-plus-halves span."""
from fractions import Fraction
from itertools import combinations
import json


def add(*polys):
    result = {}
    for poly in polys:
        for e, value in poly.items():
            result[e] = result.get(e, 0) + value
    return {e: value for e, value in result.items() if value}


def scale(poly, value):
    return {e: value * c for e, c in poly.items() if value * c}


def mul(left, right):
    result = {}
    for (i, j), a in left.items():
        for (k, l), b in right.items():
            e = (i + k, j + l)
            result[e] = result.get(e, 0) + a * b
    return {e: value for e, value in result.items() if value}


def square(poly):
    return mul(poly, poly)


# Substitute u=x^2,v=y^2 to verify the half-power identity polynomially.
motzkin = {(0, 0): 1, (8, 4): 1, (4, 8): 1, (4, 4): -3}
half = add(square({(4, 2): 1, (2, 4): -1}),
           scale(square({(3, 3): 1, (1, 1): -1}), 2),
           square({(2, 2): 1, (0, 0): -1}))
assert half == motzkin

# P(r^6,s^6), with A=r^4s^2, B=r^2s^4, C=1.
abc = [{(4, 2): 1}, {(2, 4): 1}, {(0, 0): 1}]
sqrt_abc = [{(2, 1): 1}, {(1, 2): 1}, {(0, 0): 1}]
roots = [mul(h, add(a, scale(b, -1)))
         for h in sqrt_abc for a, b in combinations(abc, 2)]
assert all(0 <= i <= 6 and 0 <= j <= 6 for root in roots for i, j in root)
target = {(0, 0): 1, (12, 6): 1, (6, 12): 1, (6, 6): -3}
assert add(*(scale(square(root), Fraction(1, 2)) for root in roots)) == target

# Half the Newton triangle: j>=i/2, i>=j/2, i+j<=9.
J = [0, 2, 3, 4, 6]
support = {(i, j) for i in J for j in J
           if 2*j >= i and 2*i >= j and i+j <= 9}
expected = {(0,0),(2,2),(2,3),(2,4),(3,2),(3,3),(3,4),
            (3,6),(4,2),(4,3),(4,4),(6,3)}
assert support == expected


def pairs_for(exponent, allowed):
    return {(e, f) for e in allowed for f in allowed if e <= f
            and (e[0]+f[0], e[1]+f[1]) == exponent
            and (e[0] % 2, e[1] % 2) == (f[0] % 2, f[1] % 2)}


steps = [((4,8),(2,4)), ((4,6),(2,3)), ((8,4),(4,2)),
         ((6,4),(3,2)), ((6,8),(3,4)), ((8,6),(4,3)),
         ((8,8),(4,4)), ((4,4),(2,2))]
for exponent, root in steps:
    assert target.get(exponent, 0) == 0
    assert pairs_for(exponent, support) == {(root, root)}
    support.remove(root)
assert pairs_for((6,6), support) == {((3,3),(3,3))}
assert target[(6,6)] < 0

print(json.dumps({
    "half_power_Motzkin_identity": True,
    "nine_square_full_sixth_grid_identity": True,
    "newton_support_exact": True,
    "successive_forced_zero_steps": len(steps),
    "negative_coefficient_contradiction": True,
    "native_source_claim": False,
}, indent=2))
