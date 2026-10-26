"""gramtools: Newton-pruned feature basis (parity blocks) of the sextic Puiseux-Gram lift at a rational source (k = 6).
Only make_blocks (and its helper in_hull) is used by the admission check."""
import numpy as np
from scipy.optimize import linprog

K = 6


def in_hull(pts, t):
    pts = np.array(pts, float); A = np.vstack([pts.T, np.ones(len(pts))])
    return linprog(np.zeros(len(pts)), A_eq=A, b_eq=[t[0], t[1], 1.0], bounds=(0, None), method='highs').status == 0


def make_blocks(fc):
    """Newton-pruned basis {(i,j,l)} split into the 4 parity blocks (k = 6 is even)."""
    cls = {}
    for l in range(4):
        S = [(a, b) for (a, b, ll, mm) in fc if ll == l and mm == l and fc[(a, b, ll, mm)] != 0]
        for i in range(K + 1):
            for j in range(K + 1):
                if in_hull(S, (2*i/K, 2*j/K)): cls.setdefault((i % 2, j % 2), []).append((i, j, l))
    return [cls[c] for c in sorted(cls)]


