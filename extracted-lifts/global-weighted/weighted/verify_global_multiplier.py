#!/usr/bin/env python3
"""Exact checks and evaluated source-independent data for the multiplier lemma.

Python 3.10+ and SymPy 1.14.0 suffice.  No SDP solver or numerical tolerance
is used.  Run with --write-schema to regenerate multiplier_schema.json.
This script does not read or recompute the existing order-669 matrices.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
E = tuple(e for e in itertools.product((-1, 0, 1), repeat=2) if e != (0, 0))
N = tuple(sorted(
    {(3 * a, 3 * b) for a, b in E}
    | set(itertools.product((-1, 1), repeat=2))
    | set(itertools.product((-1, 1), (-2, 2)))
    | set(itertools.product((-2, 2), (-1, 1)))
    | {(0, 0)}
))
NI = {n: i for i, n in enumerate(N)}


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def scale(c, x):
    return tuple(c * a for a in x)


def circuits():
    out = []
    for e in E:
        ne = scale(-1, e)
        if e < ne:
            out.append({"kind": "opposite", "exponents": [e, ne], "weights": [1, 1]})
    for a, b in itertools.product((-1, 1), repeat=2):
        out.append({"kind": "triangle111", "exponents": [(a, b), (-a, 0), (0, -b)],
                    "weights": [1, 1, 1]})
    for a in (-1, 1):
        out.append({"kind": "triangle112", "exponents": [(a, -1), (a, 1), (-a, 0)],
                    "weights": [1, 1, 2]})
    for b in (-1, 1):
        out.append({"kind": "triangle112", "exponents": [(-1, b), (1, b), (0, -b)],
                    "weights": [1, 1, 2]})
    assert len(out) == 12
    return out


CIRCUITS = circuits()


def angular_dyads(C):
    """Return (positive rational weight, plus index, minus index)."""
    ee = C["exponents"]
    if C["kind"] == "opposite":
        return [(sp.Integer(1), scale(3, ee[0]), scale(3, ee[1]))]
    if C["kind"] == "triangle111":
        return [(sp.Rational(1, 2), add(h, scale(2, l)), add(h, scale(2, m)))
                for h in ee for l, m in itertools.combinations(ee, 2)]
    e, f, g = ee
    return [(sp.Integer(1), scale(3, e), scale(3, f)),
            (sp.Integer(2), scale(-3, g), scale(3, g))]


def angular_gram(C):
    M = sp.zeros(len(N))
    for w, p, m in angular_dyads(C):
        assert p in NI and m in NI
        v = sp.zeros(len(N), 1)
        v[NI[p]] += 1
        v[NI[m]] -= 1
        M += w * v * v.T
    return M


ANGULAR_GRAMS = [angular_gram(C) for C in CIRCUITS]


def sparse_matrix(M):
    return {
        "shape": list(M.shape),
        "entries": [[i, j, int(M[i, j].p), int(M[i, j].q)]
                    for i in range(M.rows) for j in range(M.cols) if M[i, j] != 0],
    }


def sparse_polynomial(p, variables):
    return [[list(exp), int(coef.p), int(coef.q)]
            for exp, coef in sp.Poly(sp.expand(p), *variables).terms()]


def row_data():
    a0, a1, b0, b1 = sp.symbols("a0 a1 b0 b1")
    W, X, Y, Z = sp.symbols("W X Y Z")
    variables = (a0, a1, b0, b1, W, X, Y, Z)
    R = sp.Matrix([
        [b0, b1, 0, 0],
        [0, 1, 0, 0],
        [a0*b1, -a0*b0, a1*b1, -a1*b0],
        [a1, 0, -a0, 0],
        [0, 0, 1, 0],
        [a1*b0, a1*b1, -a0*b0, -a0*b1],
        [0, 0, b1, -b0],
        [0, a0, 0, a1],
    ])
    return variables, R


def fixed_features():
    variables, R = row_data()
    psi = sp.Matrix(variables[4:])
    feature_data = []
    for i, j in itertools.combinations(range(8), 2):
        for k in range(8):
            if k in (i, j):
                continue
            minor = sp.Matrix.vstack(R[i, :3], R[j, :3], R[k, :3]).det()
            expr = sp.expand(minor * (R[k, :] * psi)[0])
            feature_data.append({
                "ijk_1_based": [i + 1, j + 1, k + 1],
                "terms": sparse_polynomial(expr, variables),
            })
    return variables, R, feature_data


def test_all_positive_circuits():
    """Enumerate every minimal positive dependence of the eight exponents."""
    found = set()
    for count in (2, 3):
        for S in itertools.combinations(E, count):
            A = sp.Matrix(2, count, lambda i, j: S[j][i])
            ker = A.nullspace()
            if len(ker) != 1:
                continue
            v = ker[0]
            if all(t < 0 for t in v):
                v = -v
            if not all(t > 0 for t in v):
                continue
            den = sp.ilcm(*(t.q for t in v))
            iv = [int(t * den) for t in v]
            gcd = sp.igcd(*iv)
            iv = [t // gcd for t in iv]
            found.add(tuple(sorted(zip(S, iv))))
    expected = {tuple(sorted(zip(C["exponents"], C["weights"]))) for C in CIRCUITS}
    assert found == expected


def test_angular_identities():
    s, t = sp.symbols("s t", positive=True)
    mon = sp.Matrix([s**a * t**b for a, b in N])
    for C, G in zip(CIRCUITS, ANGULAR_GRAMS):
        rhs = sum(w * s**(6*e[0]) * t**(6*e[1])
                  for e, w in zip(C["exponents"], C["weights"])) - sum(C["weights"])
        assert sp.cancel((mon.T * G * mon)[0] - rhs) == 0
        assert G == G.T
        # The positive rational dyads themselves are an exact PSD certificate.
        assert all(w > 0 for w, _, _ in angular_dyads(C))
    assert len(N) == 21
    assert all(0 <= 3 + sign*n[i] <= 6
               for n in N for i in (0, 1) for sign in (-1, 1))


def test_source_face_cases():
    """Check the nine fixed source-face branches against all 256 supports."""
    directions = ((1, 0), (0, 1), (1, 1), (1, -1))
    for mask in range(1 << len(E)):
        I = {e for i, e in enumerate(E) if (mask >> i) & 1}
        supported = [set(C["exponents"]) for C in CIRCUITS
                     if set(C["exponents"]) <= I]
        if not supported:
            # There is no positive dependence, hence 0 is outside conv(I).
            continue
        balanced_support = set().union(*supported)
        selected = None
        for e, sign in itertools.product(directions, (1, -1)):
            if e not in I or scale(-1, e) not in I:
                continue
            if all(sign*(e[0]*f[1]-e[1]*f[0]) >= 0 for f in I):
                selected = {e, scale(-1, e)}
                break
        if selected is not None:
            assert balanced_support == selected
        else:
            assert balanced_support == I
            assert sp.Matrix(list(I)).rank() == 2


def test_all_root_signs():
    """Every angular dyad, hence its square, works on the algebraic cover.

    With unrestricted real roots x_i^6=a_i^2, changing a root sign may
    change each dyad by a common sign; its square is unchanged.  This
    checks the exact parity condition on all 48 dyads and every nonzero
    matrix entry, so no inequality on the testing cover is needed.
    """
    for C, G in zip(CIRCUITS, ANGULAR_GRAMS):
        for _, p, m in angular_dyads(C):
            assert all((pi-mi) % 2 == 0 for pi, mi in zip(p, m))
        for i, p in enumerate(N):
            for j, m in enumerate(N):
                if G[i, j] != 0:
                    exponents = (6-p[0]-m[0], 6+p[0]+m[0],
                                 6-p[1]-m[1], 6+p[1]+m[1])
                    assert all(k >= 0 and k % 2 == 0 for k in exponents)


def test_general_adjugate_completion():
    h00, h01, h02, h11, h12, h22 = sp.symbols("h00 h01 h02 h11 h12 h22")
    b0, b1, b2, c = sp.symbols("b0 b1 b2 c")
    H = sp.Matrix([[h00, h01, h02], [h01, h11, h12], [h02, h12, h22]])
    b = sp.Matrix([b0, b1, b2])
    P = H.row_join(b).col_join(sp.Matrix([[b0, b1, b2, c]]))
    M = H.row_join(b)
    rhs = M.T * H.adjugate() * M
    rhs[3, 3] += P.det()
    assert all(sp.expand(x) == 0 for x in H.det()*P-rhs)


def test_compound_gram():
    """Full symbolic Cauchy--Binet identity for a generic 4-by-3 Gram matrix.

    The proof for any number of rows is the same degree-four coefficient
    identity; the exact eight-row specialization is also checked at rational
    values, including zero source weights and zero test coordinates.
    """
    entries = sp.symbols("r0:12")
    R = sp.Matrix(4, 3, entries)
    H = R.T * R
    rhs = sp.zeros(3)
    for i, j in itertools.combinations(range(4), 2):
        cross = R[i, :].T.cross(R[j, :].T)
        rhs += cross * cross.T
    assert all(sp.expand(x) == 0 for x in H.adjugate() - rhs)

    variables, full_R = row_data()
    a0, a1, b0, b1 = variables[:4]
    points = [
        (sp.Rational(3, 5), sp.Rational(4, 5), sp.Rational(5, 13), sp.Rational(12, 13)),
        (0, 1, 1, 0),
        (-1, 0, 0, -1),
        (-sp.Rational(3, 5), sp.Rational(4, 5), -sp.Rational(5, 13), sp.Rational(12, 13)),
    ]
    weights = [(1, 0, 2, 0, 3, 4, 5, 6), (0, 0, 0, 0, 0, 2, 2, 2)]
    for point, kk in itertools.product(points, weights):
        RR = full_R.subs(dict(zip((a0, a1, b0, b1), point)))
        P = sum((k * RR[i, :].T * RR[i, :] for i, k in enumerate(kk)), sp.zeros(4))
        P[3, 3] -= 2
        H = P[:3, :3]
        rhs = sp.zeros(4)
        for i, j in itertools.combinations(range(8), 2):
            cross = RR[i, :3].T.cross(RR[j, :3].T)
            row = sp.zeros(1, 4)
            for k in range(8):
                row += kk[k] * (cross.dot(RR[k, :3].T)) * RR[k, :]
            rhs += kk[i] * kk[j] * row.T * row
        rhs[3, 3] += P.det()
        assert H.det() * P == rhs


def test_whitening_obstruction():
    u, v, a = sp.symbols("u v a")
    W, X, Y, Z = sp.symbols("W X Y Z")
    Du, Dv = 1+u*u, 1+v*v
    F = (4*Du*(W+v*X)**2 + 4*Du*Dv*X**2
         + 4*(v*W-X+u*v*Y-u*Z)**2 + 4*Dv*(u*W-Y)**2
         + 4*a*Du*Dv*Y**2 + 10*(u*W+u*v*X-Y-v*Z)**2
         + 6*Du*(v*Y-Z)**2 + 6*Dv*(X+u*Z)**2 - 2*Du*Dv*Z**2)
    P = sp.hessian(F, (W, X, Y, Z)) / 2
    H = P[:3, :3].subs(u, 1)
    expected = sp.Matrix([[8*v*v+22, 14*v, -14],
                          [14*v, 32*v*v+18, -14*v],
                          [-14, -14*v, 8*a*(v*v+1)+20*v*v+14]])
    assert all(sp.expand(x) == 0 for x in H-expected)
    r = (320*v**4+621*v*v+126) / (2*(v*v+1)*(64*v*v+99))
    inv_yy = H[:2, :2].det() / H.det()
    assert sp.cancel(inv_yy - 1/(8*(v*v+1)*(a+r))) == 0
    assert r.subs(v, sp.Rational(1, 2)) == sp.Rational(241, 230)
    assert r.subs(v, 1) == sp.Rational(1067, 652)


def make_schema():
    variables, R, features = fixed_features()
    Cmat = sp.Matrix(8, 12, lambda i, j: (
        CIRCUITS[j]["weights"][CIRCUITS[j]["exponents"].index(E[i])]
        if E[i] in CIRCUITS[j]["exponents"] else 0))
    return {
        "scope": "Fixed Gram family for (det H) F, not a Gram family for F or an SDP lift of Q222.",
        "source_commit": "fb58c232e374dadf3cbb5f23e7e7413d429de884",
        "indexed_family_size": 189,
        "polynomial_feature_count": 168,
        "angular_feature_count": 21,
        "rank_claim": "Indexed family only; zero functions and dependencies are retained.",
        "cover_variables": ["a0", "a1", "b0", "b1", "W", "X", "Y", "Z", "x0", "x1", "y0", "y1"],
        "cover_equations": ["a0^2+a1^2=1", "b0^2+b1^2=1", "W^2+X^2+Y^2+Z^2=1",
                            "x0^6=a0^2", "x1^6=a1^2", "y0^6=b0^2", "y1^6=b1^2"],
        "cover_inequalities": [],
        "cover_root_convention": "All real solutions of the four sixth-power equations are allowed; no root sign is selected.",
        "root_sign_certificate": "In every circuit dyad, its two N-indices agree coordinatewise modulo2; each squared dyad is invariant under all four independent root-sign flips.",
        "polynomial_variables": [str(x) for x in variables],
        "polynomial_format": "Each term is [exponent_vector, rational_numerator, positive_denominator].",
        "rows_R": [[sparse_polynomial(R[i, j], variables) for j in range(4)] for i in range(8)],
        "kappa": ["4 lambda", "4 n", "4 m", "4 C", "4 A", "4 D+2 d", "4 alpha+2 d", "4 beta+2 d"],
        "polynomial_features": features,
        "E": E,
        "N": N,
        "angular_features": [
            {"nu": n, "coefficient": 4,
             "exponents_x0_x1_y0_y1_Z": [3-n[0], 3+n[0], 3-n[1], 3+n[1], 1]}
            for n in N],
        "circuits": CIRCUITS,
        "circuit_ray_matrix": sparse_matrix(Cmat),
        "angular_gram_matrices": [sparse_matrix(M) for M in ANGULAR_GRAMS],
        "angular_psd_dyads": [
            [[int(w.p), int(w.q), NI[p], NI[m]] for w, p, m in angular_dyads(C)]
            for C in CIRCUITS],
        "state_gram_rule": "For each i<j use the 6x6 block kappa_i*kappa_j*kappa_without_ij*kappa_without_ij^T.",
        "angular_gram_rule": "See multiplier_lemma.tex equations (13)-(16); all choices are source-only.",
        "source_face_selection": {
            "h_nonnegative": "Use diagonal Gram h at index0 and q_e at index3e.",
            "h_negative_pair_directions": [[1, 0], [0, 1], [1, 1], [1, -1]],
            "pair_signs_in_order": [1, -1],
            "pair_test": "q_e>0 and q_-e>0 and q_f=0 whenever sigma*det(e,f)<0; choose first match.",
            "pair_minimizer": "t=sqrt(q_-e/q_e); (Uc,Vc)=(t,1) if e1=1, else (1,t).",
            "no_pair_match": "Use all active exponents; solve the two bidegree-(2,2) stationary equations for the unique positive pair.",
            "balanced_decomposition": "Repeatedly choose first supported circuit in the saved order; subtract min(w_e/b_e)*b; at most8 steps.",
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-schema", action="store_true")
    args = parser.parse_args()
    checks = [test_all_positive_circuits, test_angular_identities, test_source_face_cases,
              test_all_root_signs,
              test_general_adjugate_completion, test_compound_gram,
              test_whitening_obstruction]
    for check in checks:
        check()
        print(f"PASS {check.__name__}")
    schema = make_schema()
    zero_features = sum(all(n == 0 for _, n, _ in f["terms"])
                        for f in schema["polynomial_features"])
    print(f"Indexed features: {schema['indexed_family_size']}; identically zero polynomial indices: {zero_features}")
    print("Evaluated circuit ray matrix: 8 x 12; angular Gram matrices: 12 x (21 x 21)")
    if args.write_schema:
        path = HERE / "multiplier_schema.json"
        path.write_text(json.dumps(schema, indent=2) + "\n", encoding="utf-8")
        print(f"WROTE {path.name}")


if __name__ == "__main__":
    main()
