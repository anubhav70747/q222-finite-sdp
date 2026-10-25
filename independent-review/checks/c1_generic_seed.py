"""Exact C1 seed and generic seven-root obstruction checks.

No solver, sweep, or symbolic five-parameter expansion. Writes its matching JSON.
Uses the public separator's universal-word helpers, without executing its
top-level calculations or JSON writer. Requires SymPy.
"""
import json
from pathlib import Path
import sympy as S
from sympy.polys.matrices import DomainMatrix


def main():
    checks = {}
    folder = Path(__file__).resolve().parent
    helper_path = Path(__file__).with_name("c1_npa_separator.py")
    text = helper_path.read_text(encoding="utf-8")
    prefix = text.split("upper=[", 1)[0]
    assert "write_text" not in prefix
    helper = {}
    exec(compile(prefix, str(helper_path), "exec"), helper)
    roots, product = helper["roots"], helper["product"]
    u, v = S.symbols("u v", positive=True)
    Du, Dv = 1+u*u, 1+v*v
    I, Z, X = S.eye(2), S.diag(1, -1), S.Matrix([[0, 1], [1, 0]])
    A = [Z, ((1-u*u)*Z+2*u*X)/Du]
    B = [Z, ((1-v*v)*Z+2*v*X)/Dv]

    def op(word):
        left, right = I, I
        for i in word[0]:
            left = left*A[i]
        for j in word[1]:
            right = right*B[j]
        return S.kronecker_product(left, right)

    root_ops = [sum((q*op(w) for w, q in root.items()), S.zeros(4))
                for root in roots]
    psi = S.Matrix([1, -S.Rational(1, 2), S.Rational(1, 4), 1])
    contact = {u: 1, v: 1}
    T = S.Matrix.hstack(*[r.subs(contact)*psi for r in root_ops])
    N = S.Matrix([[1, S.Rational(5, 3), 2],
                  [-S.Rational(7, 4), -S.Rational(4, 3), -S.Rational(17, 4)],
                  [-S.Rational(5, 3), -S.Rational(16, 9), -5],
                  [-S.Rational(3, 4), -S.Rational(7, 3), -S.Rational(13, 4)],
                  [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert T*N == S.zeros(4, 3) and N[4:, :] == S.eye(3)
    pivot = T[:, :4].det()
    assert pivot != 0 and T.rank() == 4
    checks["supplied_kernel_basis_and_first_four_column_pivot"] = True

    forbidden = [((), (0, 1)), ((0,), (0, 1)), ((0, 1), ()),
                 ((0, 1), (0,)), ((0, 1), (0, 1)), ((0, 1), (1,))]
    sbasis = []
    for i, j in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]:
        b = S.zeros(3)
        b[i, j] = b[j, i] = 1
        sbasis.append(b)
    coeffs = []
    for b in sbasis:
        Q = N*b*N.T
        polynomial = {}
        for i in range(7):
            for j in range(7):
                for w, value in product(roots[i], roots[j]).items():
                    polynomial[w] = polynomial.get(w, 0)+Q[i, j]*value
        coeffs.append([polynomial.get(w, 0) for w in forbidden])
    M = S.Matrix.hstack(*[S.Matrix(c) for c in coeffs])
    determinant = S.factor(M.det())
    assert determinant == S.Rational(7360, 243)
    checks["six_actual_universal_forbidden_word_coefficients_det_7360_over_243"] = True

    # The source is the exact EN9 seed; check NC4 independently at EN13 values.
    active = list(map(S.Rational, ["27/2", 4, 2, 9, "5/4", "13/4", "9/4"]))
    d, m, C, AA, D, al, beta = active
    k = SS = ell = eta = zeta = xx = S.Integer(1)
    rr, theta = S.Rational(1, 4), S.Rational(1, 2)
    an, tn = (1-rr)/eta, theta/zeta
    h, y = 1-theta, xx*eta*zeta
    K, R, TT, GG = 1+xx, 1+xx*rr, 1+xx*h, h+zeta
    MM, J, E, QQ = R-tn, TT-an, 1-xx*eta*h, 1-xx*zeta*rr
    P4 = TT*E+xx*(1-rr)*GG
    Pal = R*h*(1+y)+zeta*MM*E
    Pd = (1+2*xx*h+rr*xx*xx*GG+xx*zeta
          +eta*xx*(K*zeta+xx*h*theta))
    B5 = K*(1+y)*(rr+eta)-(eta*TT-1+rr)*Pd
    Bb = K*xx*(rr+eta)*GG-theta*Pd
    normalized = [1, K*xx*(rr+eta)*h*(1+y)/(2*MM*zeta*Pd),
                  K*(rr+eta)*E/(2*(1-rr)*Pd), eta*P4/(2*rr*Pd),
                  B5/(2*J*eta*Pd), xx*xx*zeta*(rr+eta)*Pal/(2*QQ*Pd),
                  Bb/(2*theta*Pd)]
    assert [d*q for q in normalized] == active
    gaps = [rr, 1-rr, theta, h, E, QQ, MM, J, B5, Bb]
    assert all(q > 0 for q in gaps) and 0 < an < 2 and 0 < tn < 2
    checks["strict_native_NC4_membership"] = True

    # Full raw parent, independently from event projectors.
    AB = [[S.kronecker_product(ai, bj) for bj in B] for ai in A]
    H = S.eye(4)-(AB[0][0]-AB[0][1]+AB[1][0]+AB[1][1])/2
    events = [(1, 1, 0, 1), (1, -1, 0, 0), (1, -1, 1, 1),
              (-1, 1, 1, 0), (-1, 1, 0, 0), (-1, 1, 1, 1),
              (-1, -1, 0, 1), (1, -1, 1, 0)]
    weights = [0, 0, m, C, AA, D, al, beta]
    physical = d*H+sum((q*S.kronecker_product(I+sa*A[i], I+sb*B[j])
                              for q, (sa, sb, i, j) in zip(weights, events)), S.zeros(4))
    P = (Du*Dv*physical).applyfunc(S.cancel)
    assert P.subs(contact)*psi == S.zeros(4, 1)
    assert (psi.T*P.diff(u).subs(contact)*psi)[0] == 0
    assert (psi.T*P.diff(v).subs(contact)*psi)[0] == 0
    checks["full_physical_contact_and_six_stationary_equations"] = True

    # Complete determinant, with six nonnegative Laurent outer coefficients.
    lam = n = S.Integer(0)
    t = [d+2*(lam+n+D), d+2*(al+m+AA), d+2*(lam+C+D),
         d+2*(al+beta+m), d+2*(lam+beta+n), d+2*(al+C+AA),
         d+2*(beta+n+m), d+2*(C+AA+D)]
    q = {(-1, -1): 4*lam*al*t[6]*t[7],
         (-1, 0): 4*(lam*D*t[6]*t[5]+al*m*t[4]*t[7]),
         (-1, 1): 4*m*D*t[4]*t[5],
         (0, -1): 4*(al*AA*t[6]*t[2]+lam*n*t[7]*t[3]),
         (0, 1): 4*(beta*m*t[0]*t[5]+C*D*t[4]*t[1]),
         (1, -1): 4*n*AA*t[2]*t[3],
         (1, 0): 4*(beta*n*t[2]*t[1]+C*AA*t[0]*t[3]),
         (1, 1): 4*beta*C*t[0]*t[1]}
    assert all(value >= 0 for value in q.values())
    assert sum(e[0]*value for e, value in q.items()) == 0
    assert sum(e[1]*value for e, value in q.items()) == 0
    center = -sum(q.values())
    Laurent = center+sum(value*u**(2*e[0])*v**(2*e[1]) for e, value in q.items())
    determinant_parent = DomainMatrix.from_Matrix(P).det().as_expr()
    assert S.cancel(determinant_parent-16*Du**2*Dv**2*u*u*v*v*Laurent) == 0
    checks["full_determinant_balanced_AMGM_certificate"] = True
    integers = [int(value) for value in q.values() if value]
    assert all(value.is_Integer for value in q.values())
    from math import gcd
    from functools import reduce
    common = reduce(gcd, integers)
    reduced_q = {str(e): str(value/common) for e, value in q.items() if value}
    result = {
        "checks": checks, "number_of_checks": len(checks),
        "active_source": list(map(str, active)), "contact": ["1", "1", *map(str, psi)],
        "T": [[str(x) for x in T.row(i)] for i in range(4)],
        "first_four_columns_contact_pivot": str(pivot),
        "N": [[str(x) for x in N.row(i)] for i in range(7)],
        "forbidden_word_order": ["B0B1", "A0B0B1", "A0A1", "A0A1B0", "A0A1B0B1", "A0A1B1"],
        "symmetric_parameter_order": ["S11", "S22", "S33", "S12", "S13", "S23"],
        "M": [[str(x) for x in M.row(i)] for i in range(6)],
        "det_M": str(determinant),
        "native_coordinates": {"k": "1", "S": "1", "v_contact": "1", "r": "1/4", "theta": "1/2",
                               "eta": "1", "zeta": "1", "a": str(an), "t": str(tn)},
        "native_gaps": dict(zip(["E", "Q", "M", "J", "B5", "Bbeta"], map(str, (E, QQ, MM, J, B5, Bb)))),
        "Laurent_coefficient_common_factor": common, "reduced_outer_coefficients": reduced_q,
        "reduced_central_coefficient": str(center/common),
        "generic_result": "Open dense exclusion in the strict five-dimensional rational C1 contact parameter domain, away from the two selected rational-minor zero sets; no complete exceptional-locus classification.",
        "public_helper_executed": "definitions and root list only; no public writer or top-level separator computation"
    }
    Path(__file__).with_suffix(".json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
