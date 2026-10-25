"""Exact universal W2 certificate at one unequal-angle C2 contact.

Requires Python and SymPy only. All rational data are embedded; no external
files, numerical optimizer, parameter sweep, or coefficient-rank computation.
Run: python checks/c2_w2_unequal.py
Writes the matching c2_w2_unequal.json beside this script.
"""
import json
from pathlib import Path
import sympy as S


def reduce_letters(seq):
    out = []
    for letter in seq:
        if out and out[-1] == letter:
            out.pop()
        else:
            out.append(letter)
    return tuple(out)


def multiply(x, y):
    return (reduce_letters(x[0]+y[0]), reduce_letters(x[1]+y[1]))


def adjoint(x):
    return (tuple(reversed(x[0])), tuple(reversed(x[1])))


def product(p, q):
    out = {}
    for x, a in p.items():
        for y, b in q.items():
            w = multiply(adjoint(x), y)
            out[w] = out.get(w, 0)+a*b
    return {w: a for w, a in out.items() if a}


def main():
    N_data = [['1/6', '17/10', '6/5', '6/5', '4/3', '17/10', '221/150'], ['-16/51', '-6/5', '-22/85', '-22/85', '-26/51', '-6/5', '-58/75'], ['-175/68', '-15/4', '33/34', '-65/34', '65/34', '13/4', '-13/4'], ['1', '0', '0', '0', '0', '0', '0'], ['175/204', '5/4', '-45/34', '-45/34', '-65/102', '9/4', '13/12'], ['0', '1', '0', '0', '0', '0', '0'], ['0', '0', '1', '0', '0', '0', '0'], ['0', '0', '0', '1', '0', '0', '0'], ['0', '0', '0', '0', '1', '0', '0'], ['0', '0', '0', '0', '0', '1', '0'], ['0', '0', '0', '0', '0', '0', '1']]
    S_data = [['49146259/1004640000', '-1935923/80371200', '288229/4018560', '0', '-136427/10046400', '0', '-552231/26790400'], ['-1935923/80371200', '12689303/345644000', '-11939/229632', '-422431/975936000', '0', '0', '-98483/11481600'], ['288229/4018560', '-11939/229632', '3015449/12629760', '0', '6503/114816', '-614509/12629760', '0'], ['0', '-422431/975936000', '0', '1/16', '0', '-1309/229632', '-561/12800'], ['-136427/10046400', '0', '6503/114816', '0', '5/48', '-1/32', '3/52'], ['0', '0', '-614509/12629760', '-1309/229632', '-1/32', '1/40', '0'], ['-552231/26790400', '-98483/11481600', '0', '-561/12800', '3/52', '0', '11/128']]
    N = S.Matrix([[S.Rational(q) for q in row] for row in N_data])
    Gram = S.Matrix([[S.Rational(q) for q in row] for row in S_data])
    assert N.shape == (11, 7) and Gram.shape == (7, 7) and Gram == Gram.T
    I, A0, A1 = ((), ()), ((0,), ()), ((1,), ())
    B0, B1 = ((), (0,)), ((), (1,))
    A0B0, A0B1 = multiply(A0, B0), multiply(A0, B1)
    A1B0, A1B1 = multiply(A1, B0), multiply(A1, B1)
    roots = [{A1: 1, I: -1}, {B1: 1, I: -1}, {A0B0: 1, I: -1},
             {A1B1: 1, I: -1}, {B0: 1, A0: -1}, {A0B1: 1, A0: -1},
             {A1B0: 1, A0: -1}, {multiply(A0, A1): 1, A0: -1},
             {multiply(A1, A0): 1, A0: -1}, {multiply(B0, B1): 1, B0: -1},
             {multiply(B1, B0): 1, B0: -1}]
    kroot = []
    for j in range(7):
        p = {}
        for i in range(11):
            for w, q in roots[i].items():
                p[w] = p.get(w, 0)+N[i, j]*q
        kroot.append({w: q for w, q in p.items() if q})
    upper = [(i, j) for i in range(7) for j in range(i, 7)]
    columns = []
    for i, j in upper:
        p = product(kroot[i], kroot[j])
        if i != j:
            for w, q in product(kroot[j], kroot[i]).items():
                p[w] = p.get(w, 0)+q
        columns.append({w: q for w, q in p.items() if q})
    words = sorted(set().union(*(p.keys() for p in columns)))
    coefficient_map = S.Matrix([[p.get(w, 0) for p in columns] for w in words])
    coefficients = coefficient_map*S.Matrix([Gram[i, j] for i, j in upper])

    source = list(map(S.Rational, [1, "43/299", "77/1495", "891/598",
                                  "385/299", "119/598", "9/598"]))
    d, n, m, C, AA, D, al = source
    event_spec = [(1, 1, 0, 1), (1, -1, 0, 0), (1, -1, 1, 1),
                  (-1, 1, 1, 0), (-1, 1, 0, 0), (-1, 1, 1, 1),
                  (-1, -1, 0, 1), (1, -1, 1, 0)]
    weights = [0, n, m, C, AA, D, al, 0]
    target = {I: d, A0B0: -d/2, A0B1: d/2, A1B0: -d/2, A1B1: -d/2}
    for q, (sa, sb, i, j) in zip(weights, event_spec):
        for w, co in ((I, 1), ([A0, A1][i], sa), ([B0, B1][j], sb),
                      (multiply([A0, A1][i], [B0, B1][j]), sa*sb)):
            target[w] = target.get(w, 0)+q*co
    target = {w: q for w, q in target.items() if q}
    assert set(target) <= set(words)
    assert all(coefficients[i] == target.get(w, 0) for i, w in enumerate(words))
    print("Exact universal word identity passed on", len(words), "word rows.", flush=True)

    # Exact LDL, including the seven bounds actually printed in the source.
    LDL = S.eye(7)
    pivots = []
    bounds = list(map(S.Rational, ["1/25", "1/50", "1/10", "1/20", "1/20", "1/250", "1/1000"]))
    for i in range(7):
        pivot = Gram[i, i]-sum(LDL[i, j]**2*pivots[j] for j in range(i))
        assert pivot > bounds[i]
        pivots.append(pivot)
        for j in range(i+1, 7):
            LDL[j, i] = (Gram[j, i]-sum(LDL[j, h]*LDL[i, h]*pivots[h]
                                      for h in range(i)))/pivot
    assert LDL*S.diag(*pivots)*LDL.T == Gram
    print("Exact LDL and all seven published lower bounds passed.", flush=True)

    # Independent Pauli evaluation is for contact/native checks only;
    # it is NOT used to match the universal coefficient identity above.
    u, v = S.symbols("u v", real=True)
    eye, zz, xx = S.eye(2), S.diag(1, -1), S.Matrix([[0, 1], [1, 0]])
    As = [zz, ((1-u*u)*zz+2*u*xx)/(1+u*u)]
    Bs = [zz, ((1-v*v)*zz+2*v*xx)/(1+v*v)]
    def op(w):
        left, right = eye, eye
        for i in w[0]:
            left = left*As[i]
        for j in w[1]:
            right = right*Bs[j]
        return S.kronecker_product(left, right)
    psi = S.Matrix([1, -S.Rational(1, 8), S.Rational(1, 8), 1])
    contact = {u: S.Rational(1, 4), v: S.Rational(1, 2)}
    root_ops = [sum((q*op(w) for w, q in p.items()), S.zeros(4)) for p in roots]
    T = S.Matrix.hstack(*[p.subs(contact)*psi for p in root_ops])
    pivot_columns = [0, 1, 2, 4]
    contact_pivot = T[:, pivot_columns].det()
    assert contact_pivot != 0 and T*N == S.zeros(4, 7)
    free = [3, 5, 6, 7, 8, 9, 10]
    assert N[free, :] == S.eye(7)
    physical = sum((q*op(w) for w, q in target.items()), S.zeros(4))
    assert physical.subs(contact)*psi == S.zeros(4, 1)
    assert (psi.T*physical.diff(u).subs(contact)*psi)[0] == 0
    assert (psi.T*physical.diff(v).subs(contact)*psi)[0] == 0
    a, b, c, tt, r = map(S.Rational, [2, 2, "1/4", "1/4", "1/2"])
    z = c*c
    ac, h, g, j = 1+z*r, b-a*a*tt, (1+z)*b-(1-r+a*a*z*tt), 1+z*r-tt
    Delta = ((1+z)*z*b*b+2*z*(1+z)*b*r+(1+z)*(1+a*a*z)*b
             +z*z*(1-a*a)*r*r+2*z*(1+a*a*z)*r+1+a*a*z)
    nn = (1+z)*(b+r)*h+a*a*(g+j+r*tt*(1+a*a*z))
    native = [1, ((b+r)*ac*(ac+a*a*z*tt)-r*Delta)/(2*r*Delta),
              (1+a*a*z)*(1-r)*ac*(ac+a*a*z*tt)/(2*g*Delta),
              ac*nn/(2*a*a*tt*Delta),
              a*a*ac*(b+r)*(ac+a*a*z*tt)/(2*h*Delta),
              (1+z)*(b+r)*(b-a*a*z*r)*(ac+a*a*z*tt)/(2*a*a*j*Delta),
              z*z*(b+r)*nn/(2*(1-z*h)*Delta)]
    assert native == source
    gaps = [h, g, j, 1-z*h, b-a*a*z*r, n]
    assert all(q > 0 for q in gaps)

    # Full-polar extremality, not just extremality inside one source face.
    bell_words = [I, A0, A1, B0, B1, A0B0, A0B1, A1B0, A1B1]
    bell_ops = [op(w) for w in bell_words]
    kernel = S.Matrix.hstack(*[p.subs(contact)*psi for p in bell_ops])
    jets = S.Matrix([[(psi.T*p.diff(t).subs(contact)*psi)[0] for p in bell_ops] for t in (u, v)])
    def evaluate(w, point):
        value = 1
        for i in w[0]:
            value *= point[i]
        for j in w[1]:
            value *= point[2+j]
        return value
    points = [(1, 1, 1, 1), (-1, 1, -1, 1)]
    evaluations = S.Matrix([[evaluate(w, p) for w in bell_words] for p in points])
    extremality = kernel.col_join(jets).col_join(evaluations)
    pivcols = list(extremality.rref()[1])
    assert len(pivcols) == 8
    extreme_minor = extremality[:, pivcols].det()
    assert extreme_minor != 0
    bell_target = S.Matrix([target.get(w, 0) for w in bell_words])
    assert extremality*bell_target == S.zeros(8, 1)

    def label(w):
        return "".join("A"+str(i) for i in w[0])+"".join("B"+str(i) for i in w[1]) or "I"
    result = {
        "input": "Literal N_star and S_star embedded in this standalone checker",
        "N_literal": [[str(q) for q in N.row(i)] for i in range(11)],
        "S_literal": [[str(q) for q in Gram.row(i)] for i in range(7)],
        "universal_word_count": len(words), "forbidden_word_count": len(words)-len(target),
        "all_universal_coefficients_match": True,
        "matched_word_coefficients": {label(w): str(coefficients[i]) for i, w in enumerate(words)},
        "coefficient_map_rank": "not attempted; no neighborhood claim",
        "LDL_positive_pivots": list(map(str, pivots)),
        "LDL_claimed_strict_lower_bounds": list(map(str, bounds)),
        "all_claimed_LDL_bounds_verified": True,
        "contact_action_pivot_columns_one_based": [i+1 for i in pivot_columns],
        "contact_action_pivot": str(contact_pivot), "N_is_full_contact_kernel": True,
        "native_parameters_a_b_c_t_r": list(map(str, (a, b, c, tt, r))),
        "native_gaps_h_g_j_1zh_ba2zr_n": list(map(str, gaps)),
        "native_source_matches": True, "full_contact_and_scalar_jets_zero": True,
        "full_Bell_extremality_rank": 8,
        "extremality_minor_columns_one_based": [i+1 for i in pivcols],
        "extremality_minor": str(extreme_minor),
        "scope": "One globally valid unequal-angle maximally entangled C2 full-polar extreme ray has an exact universal W2 certificate. No neighborhood or uniform-family admission inferred."
    }
    Path(__file__).with_suffix(".json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k not in
                      ["N_literal", "S_literal", "matched_word_coefficients", "LDL_positive_pivots"]}, indent=2))


if __name__ == "__main__":
    main()
