"""Exact Boolean check of an event-positivity flag for Q222.

Only the 16 deterministic local behaviours and the 16 event probabilities
are used.  This is not a numerical SDP test or a quantum-membership test.
"""

from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path
import json


VERTICES = list(product((0, 1), repeat=4))  # (a0,a1,b0,b1)
EVENTS = list(product((0, 1), repeat=4))  # (x,y,a,b)
FULL = (1 << 16) - 1
EVENT_MASKS = [
    sum(1 << i for i, v in enumerate(VERTICES)
        if v[x] == a and v[2 + y] == b)
    for x, y, a, b in EVENTS
]


def surviving_vertices(zero_events):
    removed = 0
    for i, mask in enumerate(EVENT_MASKS):
        if (zero_events >> i) & 1:
            removed |= mask
    return FULL & ~removed


@lru_cache(None)
def longest_chain(survivors):
    """Maximum number of strict nonempty transitions from this face."""
    best = (0, ())
    for e, removed in enumerate(EVENT_MASKS):
        child = survivors & ~removed
        if child and child != survivors:
            length, suffix = longest_chain(child)
            candidate = (1 + length, (e,) + suffix)
            if candidate[0] > best[0]:
                best = candidate
    return best


def affine_rank(points):
    rows = []
    for a0, a1, b0, b1 in points:
        a0, a1, b0, b1 = (1 - 2 * a0, 1 - 2 * a1,
                          1 - 2 * b0, 1 - 2 * b1)
        rows.append(list(map(Fraction, (1, a0, a1, b0, b1,
                         a0*b0, a0*b1, a1*b0, a1*b1))))
    rank = 0
    for col in range(9):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        value = rows[rank][col]
        rows[rank] = [x / value for x in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][col]:
                value = rows[i][col]
                rows[i] = [x - value*y for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank - 1


# Enumerate all 2^16 subsets, including the empty survivor set.  The
# recursion above independently explores every nonempty survivor set.
all_sets = {surviving_vertices(s) for s in range(1 << 16)}
length, path = longest_chain(FULL)
assert length == 8
assert len(all_sets) == 1722
assert 0 in all_sets
assert longest_chain.cache_info().currsize == len(all_sets) - 1
assert affine_rank(VERTICES) == 8

# An outcome relabelling of a longest chain makes its terminal vertex 0000.
# Every event is written (x,y,a,b), meaning P(a,b|x,y)=0.
chosen_events = [
    (0, 0, 1, 0),
    (0, 0, 0, 1),
    (0, 1, 1, 1),
    (1, 0, 1, 1),
    (0, 0, 1, 1),
    (1, 1, 1, 1),
    (0, 1, 0, 1),
    (1, 0, 1, 0),
]
survivors = FULL
flag = []
for e in chosen_events:
    child = survivors & ~EVENT_MASKS[EVENTS.index(e)]
    assert child and child != survivors
    witness_index = next(i for i in range(16) if ((survivors ^ child) >> i) & 1)
    witness = VERTICES[witness_index]
    x, y, a, b = e
    assert witness[x] == a and witness[2 + y] == b
    remaining = [VERTICES[i] for i in range(16) if (child >> i) & 1]
    assert all(v[x] != a or v[2 + y] != b for v in remaining)
    flag.append({
        "added_zero_event_xyab": e,
        "removed_deterministic_witness_a0a1b0b1": witness,
        "remaining_deterministic_vertices": remaining,
        "remaining_count": len(remaining),
        "classical_face_affine_dimension": affine_rank(remaining),
    })
    survivors = child

assert [row["remaining_count"] for row in flag] == [12, 8, 6, 5, 4, 3, 2, 1]
assert [row["classical_face_affine_dimension"] for row in flag] == list(range(7, -1, -1))
assert flag[-1]["remaining_deterministic_vertices"] == [(0, 0, 0, 0)]
assert len(flag) == length

result = {
    "method": "Exact Boolean enumeration and rational Gaussian elimination",
    "event_convention": "event=(x,y,a,b), vertex=(a0,a1,b0,b1), bits in {0,1}",
    "event_subsets_enumerated": 1 << 16,
    "distinct_survivor_sets_including_empty": len(all_sets),
    "distinct_nonempty_survivor_sets": len(all_sets) - 1,
    "maximum_strict_nonempty_inclusions": length,
    "maximum_nonempty_faces_in_chain": length + 1,
    "deterministic_hull_affine_dimension": affine_rank(VERTICES),
    "single_psd_cone_order_lower_bound": length + 1,
    "first_unexcluded_order": 9,
    "first_unexcluded_template": {"order": 9, "symmetric_dimension": 45,
                                   "auxiliary_directions_at_most": 36,
                                   "matrix_coefficients_at_most": 2025},
    "flag": flag,
    "scope": "The chain gives a necessary order bound, not an upper bound or a bound on maximum block size in a product lift.",
}
out = Path(__file__).with_name("positivity_face_chain_lower_bound_checks.json")
out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps({key: value for key, value in result.items() if key != "flag"}, indent=2))
