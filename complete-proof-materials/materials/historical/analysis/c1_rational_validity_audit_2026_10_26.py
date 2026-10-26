"""Physical validity/extremality of the two task-4 C1 sources.

Self-contained exact arithmetic; requires SymPy. No NPA separator check.
"""
from pathlib import Path
import json
import sympy as S
from sympy.polys.matrices import DomainMatrix


def main():
    checks = {}
    def zero(expr, label):
        assert S.cancel(expr) == 0, label
        checks[label] = True

    u, v = S.symbols("u v", real=True)
    Du, Dv = 1+u*u, 1+v*v
    I2, I4 = S.eye(2), S.eye(4)
    Z = S.diag(1, -1)
    X = S.Matrix([[0, 1], [1, 0]])
    A = [Z, ((1-u*u)*Z+2*u*X)/Du]
    B = [Z, ((1-v*v)*Z+2*v*X)/Dv]
    kron = S.kronecker_product
    AB = [[kron(a, b) for b in B] for a in A]
    H = I4-(AB[0][0]-AB[0][1]+AB[1][0]+AB[1][1])/2
    # Repository event order lambda,n,m,C,A,D,alpha,beta.
    events = [(1, 1, 0, 1), (1, -1, 0, 0), (1, -1, 1, 1),
              (-1, 1, 1, 0), (-1, 1, 0, 0), (-1, 1, 1, 1),
              (-1, -1, 0, 1), (1, -1, 1, 0)]
    columns = [H]+[kron(I2+aa*A[i], I2+bb*B[j])
                   for aa, bb, i, j in events]
    raw_columns = [(Du*Dv*column).applyfunc(S.cancel) for column in columns]
    e = [S.Matrix(q) for q in ((1, v, 0, 0), (0, 1, 0, 0),
          (v, -1, u*v, -u), (u, 0, -1, 0), (0, 0, 1, 0),
          (u, u*v, -1, -v), (0, 0, v, -1), (0, 1, 0, u))]
    packet_columns = [2*(e[5]*e[5].T+Du*e[6]*e[6].T
                         +Dv*e[7]*e[7].T-Du*Dv*S.diag(0, 0, 0, 1))]
    packet_columns += [4*f*q*q.T for f, q in zip(
        (Du, Du*Dv, 1, Dv, Du*Dv, 1, Du, Dv), e)]
    for i, (actual, packet) in enumerate(zip(raw_columns, packet_columns)):
        assert all(S.cancel(t) == 0 for t in actual-packet)
        checks[f"Pauli_vs_full_parent_column_{i}"] = True

    bell_columns = [I4, kron(A[0], I2), kron(A[1], I2),
                    kron(I2, B[0]), kron(I2, B[1]),
                    AB[0][0], AB[0][1], AB[1][0], AB[1][1]]
    deterministic = [(1, 1, 1, 1), (1, -1, -1, -1)]
    deterministic_rows = [[1, aa, a1, bb, b1, aa*bb, aa*b1, a1*bb, a1*b1]
                          for aa, a1, bb, b1 in deterministic]
    cases = [
        ("interior", S.Rational(1, 2),
         list(map(S.Rational, [1, "20/29", "12/29", "13/58",
                               "51/58", "3/29", "3/58"]))),
        ("boundary", S.Rational(3, 8),
         list(map(S.Rational, [1, "5/6", "3/10", "1/3", 2, "1/13", 0])))
    ]
    U, V, T = S.symbols("U V T", nonnegative=True)
    M1 = (U*T)**2+U+T-3*U*T
    M2 = U*U*T+T*T+U-3*U*T
    zero((U*T)**2*U*T-(U*T)**3, "AMGM_product_M1")
    zero(U*U*T*T*T*U-(U*T)**3, "AMGM_product_M2")
    output = []
    for label, rr, active in cases:
        d, m, C, aa, D, al, beta = active
        c = [d, 0, 0, m, C, aa, D, al, beta]
        # Reconstruct NC1--NC4 at the specified source parameters.
        xx, eta, zeta, theta = map(S.Rational, [1, "1/2", "1/2", "1/2"])
        h = 1-theta
        an, tn = (1-rr)/eta, theta/zeta
        K, R, TT, G = 1+xx, 1+xx*rr, 1+xx*h, h+zeta
        y = xx*eta*zeta
        MM, J, E, Q = R-tn, TT-an, 1-xx*eta*h, 1-xx*zeta*rr
        P4 = TT*E+xx*(1-rr)*G
        Pal = R*h*(1+y)+zeta*MM*E
        p0 = 1+2*xx*h+rr*xx*xx*G+xx*zeta
        p1 = xx*(K*zeta+xx*h*theta)
        Pd = p0+eta*p1
        B5 = K*(1+y)*(rr+eta)-(eta*TT-1+rr)*Pd
        Bb = K*xx*(rr+eta)*G-theta*Pd
        native = [1,
            K*xx*(rr+eta)*h*(1+y)/(2*MM*zeta*Pd),
            K*(rr+eta)*E/(2*(1-rr)*Pd),
            eta*P4/(2*rr*Pd), B5/(2*J*eta*Pd),
            xx*xx*zeta*(rr+eta)*Pal/(2*Q*Pd), Bb/(2*theta*Pd)]
        assert native == active
        assert all(q > 0 for q in (rr, theta, 1-rr, h, E, Q, MM, J))
        assert 0 < an < 2 and 0 < tn < 2 and B5 >= 0 and Bb >= 0
        checks[f"{label}_native_dictionary_and_all_gaps"] = True

        P = sum((ci*pi for ci, pi in zip(c, raw_columns)), S.zeros(4))
        psi = S.Matrix([1, -S.Rational(1, 2), rr, 1])
        contact = {u: 1, v: S.Rational(1, 2)}
        assert P.subs(contact)*psi == S.zeros(4, 1)
        assert (psi.T*P.diff(u).subs(contact)*psi)[0] == 0
        assert (psi.T*P.diff(v).subs(contact)*psi)[0] == 0
        checks[f"{label}_four_kernel_and_two_stationary_identities"] = True
        det = DomainMatrix.from_Matrix(P).det().as_expr()
        if label == "interior":
            gp = (670*V*(U-1)**2+S.Rational(449, 2)*U*(T-1)**2
                  +S.Rational(55, 2)*M1+255*M2)
        else:
            gp = (235*V*(U-1)**2+S.Rational(153, 2)*U*(T-1)**2
                  +S.Rational(295, 2)*M2)
        Gtest = gp.subs(T, 4*V).subs({U: u*u, V: v*v})
        prefactor = S.cancel(det/(Du**2*Dv**2*Gtest))
        assert prefactor.is_Rational and prefactor > 0
        zero(det-prefactor*Du**2*Dv**2*Gtest, f"{label}_full_determinant_AMGM_identity")
        origin = P.subs({u: 0, v: 0})
        assert origin == S.diag(0, 4*m+4*beta+2*d,
                                4*C+4*aa+4*D+2*d, 4*al)
        checks[f"{label}_origin_PSD_diagonal"] = True

        # Full Bell coefficient space, not merely the canonical orthant.
        kernel_rows = S.Matrix.hstack(*[o.subs(contact)*psi for o in bell_columns])
        jets = S.Matrix([[(psi.T*o.diff(t).subs(contact)*psi)[0]
                          for o in bell_columns] for t in (u, v)])
        equations = kernel_rows.col_join(jets).col_join(S.Matrix(deterministic_rows))
        minor = S.factor(equations[:, :8].det())
        assert minor != 0 and equations.rank() == 8
        expected_minor = (-S.Rational(576, 125) if label == "interior"
                          else -S.Rational(273, 100))
        assert minor == expected_minor
        coeff = S.Matrix([d, 0, 0, 0, 0, -d/2, d/2, -d/2, -d/2])
        for weight, (sa, sb, i, j) in zip(c[1:], events):
            coeff[0] += weight
            coeff[1+i] += sa*weight
            coeff[3+j] += sb*weight
            coeff[5+2*i+j] += sa*sb*weight
        assert equations*coeff == S.zeros(8, 1)
        checks[f"{label}_full_Bell_contact_rank_eight"] = True
        checks[f"{label}_source_in_full_contact_nullspace"] = True
        normalized = sum((ci*pi for ci, pi in zip(c, columns)), S.zeros(4))
        assert normalized.subs(contact)*psi == S.zeros(4, 1)
        norm2 = (psi.T*psi)[0]
        hcontact = S.factor((psi.T*H.subs(contact)*psi)[0]/norm2)
        assert hcontact < 0
        output.append({
            "case": label, "active_order": ["d", "m", "C", "A", "D", "alpha", "beta"],
            "active_source": list(map(str, active)), "contact_state": list(map(str, psi)),
            "native_gaps": dict(zip(["E", "Q", "M", "J", "B5", "Bbeta"],
                                     map(str, (E, Q, MM, J, B5, Bb)))),
            "determinant_prefactor": str(prefactor),
            "determinant_identity": "det P = prefactor*(1+u^2)^2*(1+v^2)^2*G(u^2,v^2,4*v^2)",
            "full_Bell_column_order": ["1", "A0", "A1", "B0", "B1", "E00", "E01", "E10", "E11"],
            "row_order": ["four normalized operator-kernel rows", "u scalar jet", "v scalar jet", "v_lambda", "v_n"],
            "first_eight_columns_minor": str(minor),
            "claimed_minor_matches": minor == expected_minor,
            "full_Bell_source": list(map(str, coeff)),
            "H_at_normalized_contact": str(hcontact),
            "origin_matrix": [[str(x) for x in origin.row(i)] for i in range(4)]
        })
    result = {"checks": checks, "number_of_checks": len(checks), "sources": output,
              "scope": "Direct global quantum validity and full nine-coordinate Bell-polar extremality of two sources; universal seven-root separator excluded from this audit."}
    Path(__file__).with_suffix(".json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
