"""Exact universal-algebra elimination; no PSD solver or parameter sampling.

Requires SymPy. Run this file to print/replay both affine Gram families.
The universal algebra is R[Z2*Z2] tensor R[Z2*Z2].
"""
import json
from pathlib import Path
import sympy as s


def reduce_word(xs):
    out = []
    for x in xs:
        if out and out[-1] == x:
            out.pop()
        else:
            out.append(x)
    return tuple(out)


def multiply(a, b):
    return reduce_word(a[0] + b[0]), reduce_word(a[1] + b[1])


words = [
    ((), ()), ((0,), ()), ((1,), ()), ((), (0,)), ((), (1,)),
    ((0,), (0,)), ((0,), (1,)), ((1,), (0,)), ((1,), (1,)),
]
patterns = {
    "C1": [(1, 0), (7, 0), (8, 0), (3, 2), (4, 2), (5, 2), (6, 2)],
    "C2": [(2, 0), (4, 0), (5, 0), (8, 0), (3, 1), (6, 1), (7, 1)],
}
d, la, n, m, C, A, D, al, be = s.symbols("d lambda n m C A D alpha beta")
sources = (d, la, n, m, C, A, D, al, be)
target = [sum(sources), la+n-A-al, m-C-D+be, -n+C+A-be,
          la-m+D-al, -d/2-n-A, d/2+la+al, -d/2-C-be, -d/2-m-D]
pairs = [(i, j) for i in range(7) for j in range(i, 7)]
unknowns = s.symbols(" ".join(f"q{i+1}{j+1}" for i, j in pairs))


def calculate(name, basis, return_symbolic=False):
    universal = sorted({multiply(a, b) for a in words for b in words})
    rows = {w: i for i, w in enumerate(universal)}
    L = s.zeros(len(universal), 28)
    for col, (a, b) in enumerate(pairs):
        oriented = [(a, b)] if a == b else [(a, b), (b, a)]
        for i, j in oriented:
            for wi, ci in [(basis[i][0], 1), (basis[i][1], -1)]:
                for wj, cj in [(basis[j][0], 1), (basis[j][1], -1)]:
                    L[rows[multiply(words[wi], words[wj])], col] += ci*cj
    rhs = s.zeros(len(universal), 1)
    sub = {la: 0, n: 0} if name == "C1" else {la: 0, be: 0}
    for w, value in zip(words, target):
        rhs[rows[w]] = value.subs(sub)
    solution = list(s.linsolve((L, rhs), unknowns))
    assert len(solution) == 1
    vector = s.Matrix(solution[0])
    assert L*vector == rhs
    gram = s.zeros(7)
    for (i, j), value in zip(pairs, vector):
        gram[i, j] = gram[j, i] = s.factor(value)
    free = sorted(set().union(*(v.free_symbols for v in vector)) & set(unknowns), key=str)
    # Independent reconstruction of the polynomial, without the matching matrix.
    coefficient = {}
    for i in range(7):
        for j in range(7):
            for wi, ci in [(basis[i][0], 1), (basis[i][1], -1)]:
                for wj, cj in [(basis[j][0], 1), (basis[j][1], -1)]:
                    w = multiply(words[wi], words[wj])
                    coefficient[w] = coefficient.get(w, 0)+ci*cj*gram[i, j]
    for w in universal:
        assert s.expand(coefficient.get(w, 0)-rhs[rows[w]]) == 0
    if return_symbolic:
        return gram, free
    return {
        "rank": L.rank(), "free_dimension": len(free), "universal_word_rows": len(universal),
        "free_parameters": [str(x) for x in free],
        "gram": [[str(v) for v in gram.row(i)] for i in range(7)],
        "all_universal_coefficients_reconstructed": True,
    }


if __name__ == "__main__":
    data = {name: calculate(name, basis) for name, basis in patterns.items()}
    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps(data, indent=2))
