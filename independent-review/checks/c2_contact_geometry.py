"""Exact checks for ../C2_CONTACT_GEOMETRY.md. Requires Python and SymPy.

Self-contained: no private analysis files, network, numerical solver, or
source archive. This checks formulas and witnesses, not a full-hull theorem.
Run: python independent-review/checks/c2_contact_geometry.py
"""
import json
from pathlib import Path
import sympy as sp


def main():
    checks = {}

    def zero(expression, name):
        assert sp.cancel(expression) == 0, name
        checks[name] = True

    a, r, s, z = sp.symbols("a r s z", positive=True)
    e, p = r+s, r*s
    g = a+e-1
    L = a*(e-p)+p
    Q = a*g+L
    Delta = a**3*p+a*a*(e*e-p)+a*p*(2*e-p)+p*p
    B0 = a**3*p+a*a*(e*e+p-e)+a*((3*e-2)*p-p*p)+p*p
    # Active order: d,n,m,C,A,D,alpha.
    source = [sp.Integer(1), -a*B0/(2*(a+1)*Delta),
              (a+s)*(1-r)*(a*s+r)*L/(2*(a+1)*g*Delta),
              r*r*(a+s)*Q/(2*(1-r)*Delta),
              a*(a+r)*(a+s)*L/(2*(a+1)*Delta),
              (a+r)*(1-s)*(a*r+s)*L/(2*(a+1)*g*Delta),
              s*s*(a+r)*Q/(2*(1-s)*Delta)]

    # Independent native chart, z=c^2. This is rebuilt here, not imported.
    q = 1+(a*(1-r)+r)*z
    ell = a+r-1+a*r*z
    j = 1-a*r*z
    rr = ell+r*q
    w = 1+(a*a+a+2*r)*z+(a*a+a*r*(2-r)+r*r)*z*z
    nn = ((1-r)+(a-r*(a*a+3*a)+r*(2-3*r))*z
          -r*(a*a+a*r*(3-r)+r*r)*z*z)
    native = [sp.Integer(1), a*nn/(2*r*(a+1)*w),
              (1-r)*(1+a*a*z)*(1+r*z)*q/(2*(a+1)*ell*w),
              (1+r*z)*rr/(2*(1-r)*w),
              a*(a+r)*(1+r*z)*q/(2*r*(a+1)*w),
              (a+r)*(1+z)*j*q/(2*(a+1)*ell*w),
              a*z*z*(a+r)*rr/(2*j*w)]
    for i, (left, right) in enumerate(zip(source, native)):
        zero(left-right.subs(z, s/(a*r)), f"eq1_native_source_{i}")
    zero(nn.subs(z, 1/(a*r))+(a+1)**2*(a+r)/a**2,
         "native_positive_root_excludes_s_ge_1")

    # Positive common denominator; homogeneous identities become polynomials.
    scale = 2*r*(a+1)*(1-r)*ell*j*w
    poly = [sp.Poly(sp.cancel(value*scale), a, r, z, domain=sp.QQ)
            for value in native]
    d, n, m, C, A, D, alpha = poly
    half = d.mul_ground(sp.Rational(1, 2))
    B, U, V, M = A+n, C+half, alpha+half, m+D
    S = U+V
    Gamma = a*a*(B+S)-B
    identities = {
        "eq2_angle_recovery": a*(d+2*n)-2*A,
        "eq3_minor_mU_DV": m*U-D*V,
        "eq3_minor_mGamma_BV": m*Gamma-B*V,
        "eq3_minor_DGamma_BU": D*Gamma-B*U,
        "eq4_summed_minors": a*a*M*(B+S)-B*(M+S),
        "eq6_full_contact_quartic": S*S*(B+half)*(M+half)
                                  -(a+1)**2*U*V*M*(B+S),
    }
    for name, value in identities.items():
        assert value.is_zero, name
        checks[name] = True

    # Eq(5): prove invariance symbolically before choosing a source.
    d0, n0, m0, C0, A0, D0, alpha0 = sp.symbols("d n m C A D alpha")
    coordinates = sp.Matrix([d0, n0, m0, C0, A0, D0, alpha0])
    direction = sp.Matrix([1, -a/(2*(a+1)), 0, -sp.Rational(1, 2),
                           a/(2*(a+1)), 0, -sp.Rational(1, 2)])
    epsilon = sp.symbols("epsilon", real=True)

    def invariants(values):
        di, ni, mi, Ci, Ai, Di, ai = values
        bi, ui, vi, Mi = Ai+ni, Ci+di/2, ai+di/2, mi+Di
        si = ui+vi
        gi = a*a*(bi+si)-bi
        return [bi, ui, vi, Mi, a*(di+2*ni)-2*Ai,
                mi*ui-Di*vi, mi*gi-bi*vi, Di*gi-bi*ui]

    for i, (old, new) in enumerate(zip(invariants(coordinates),
                                      invariants(coordinates+epsilon*direction))):
        zero(new-old, f"eq5_invariant_{i}")

    # Eq(7), including the distinction between its two algebraic branches.
    k, t, x, theta = sp.symbols("k t x theta")
    weighted = ((k+t+1)**2*x*theta
                -(x+(k-t)**2*theta)*(x+theta))
    quadratic = (x-(2*k*t+k+t)*theta)**2-4*k*t*(k+1)*(t+1)*theta**2
    zero(weighted+quadratic, "eq7_weighted_CHSH_factorization")
    pa, qa = t*(k+1), k*(t+1)
    angle = a*a*pa-qa
    factored = (x-(a+1)**2*pa*theta)*(x-(a-1)**2*pa*theta)
    remainder = sp.rem(sp.Poly(sp.expand(factored-quadratic), a),
                       sp.Poly(angle, a))
    assert remainder.is_zero
    checks["eq7_two_root_factorization_modulo_eq4"] = True

    # Two signed mixtures are evaluated as exact rational source sums.
    witnesses = []
    pairs = [((sp.Rational(1, 8), sp.Rational(1, 16)),
              (sp.Rational(1, 16), sp.Rational(1, 8))),
             ((sp.Rational(1, 16), sp.Rational(1, 8)),
              (sp.Rational(1, 8), sp.Rational(1, 16)))]
    for rs1, rs2 in pairs:
        c1 = [sp.cancel(v.subs({a: 1, r: rs1[0], s: rs1[1]})) for v in source]
        c2 = [sp.cancel(v.subs({a: 2, r: rs2[0], s: rs2[1]})) for v in source]
        assert all(v > 0 for v in c1+c2)
        R1 = sp.Matrix([[c1[2], c1[5]], [c1[6]+c1[0]/2, c1[3]+c1[0]/2]])
        R2 = sp.Matrix([[c2[2], c2[5]], [c2[6]+c2[0]/2, c2[3]+c2[0]/2]])
        determinant = (R1+R2).det()
        kappa1 = sp.cancel(R1[1, 0]/R1[0, 0])
        kappa2 = sp.cancel(R2[1, 0]/R2[0, 0])
        assert determinant == (kappa2-kappa1)*(c1[2]*c2[5]-c1[5]*c2[2])
        witnesses.append({"a1": 1, "a2": 2, "rs1": list(map(str, rs1)),
                          "rs2": list(map(str, rs2)),
                          "source1": list(map(str, c1)), "source2": list(map(str, c2)),
                          "determinant": str(determinant)})
    expected = sp.Rational(16557677, 50331190)
    assert [sp.Rational(v["determinant"]) for v in witnesses] == [-expected, expected]
    checks["both_signed_mixtures"] = True

    # Reconstruct the invalid source directly, rather than loading a witness.
    base = sp.Matrix([sp.cancel(v.subs({a: 2, r: sp.Rational(1, 2),
                                        s: sp.Rational(1, 16)})) for v in source])
    literal = sp.Matrix([1, sp.Rational(43, 299), sp.Rational(77, 1495),
                         sp.Rational(891, 598), sp.Rational(385, 299),
                         sp.Rational(119, 598), sp.Rational(9, 598)])
    assert base == literal
    shifted = base+sp.Rational(1, 1000)*direction.subs(a, 2)
    assert all(v > 0 for v in shifted)
    di, ni, mi, Ci, Ai, Di, al = shifted
    bi, ui, vi, Mi = Ai+ni, Ci+di/2, al+di/2, mi+Di
    si = ui+vi
    quartic_residual = sp.factor(si**2*(bi+di/2)*(Mi+di/2)
                                 -9*ui*vi*Mi*(bi+si))
    assert quartic_residual == sp.Rational(899901143099, 106923596000000)

    # Full uncleared Bell operator in Pauli half-angle coordinates.
    u, v = sp.symbols("u v", real=True)
    I = sp.eye(2)
    Z = sp.diag(1, -1)
    Xp = sp.Matrix([[0, 1], [1, 0]])
    opA0 = sp.kronecker_product(Z, I)
    opA1 = sp.kronecker_product(((1-u*u)*Z+2*u*Xp)/(1+u*u), I)
    opB0 = sp.kronecker_product(I, Z)
    opB1 = sp.kronecker_product(I, ((1-v*v)*Z+2*v*Xp)/(1+v*v))
    operators = [sp.eye(4), opA0, opA1, opB0, opB1,
                 opA0*opB0, opA0*opB1, opA1*opB0, opA1*opB1]

    def bell(values):
        di, ni, mi, Ci, Ai, Di, al = values
        coefficient = [di+ni+mi+Ci+Ai+Di+al, ni-Ai-al, mi-Ci-Di,
                       -ni+Ci+Ai, -mi+Di-al, -di/2-ni-Ai,
                       di/2+al, -di/2-Ci, -di/2-mi-Di]
        return sum((q*op for q, op in zip(coefficient, operators)), sp.zeros(4))

    contact = {u: sp.Rational(1, 4), v: sp.Rational(1, 2)}
    psi = sp.Matrix([1, -sp.Rational(1, 8), sp.Rational(1, 8), 1])
    original = bell(base)
    assert original.subs(contact)*psi == sp.zeros(4, 1)
    assert all(sp.factor((psi.T*sp.diff(original, angle).subs(contact)*psi)[0]) == 0
               for angle in (u, v))
    slack = sp.factor((psi.T*bell(shifted).subs(contact)*psi)[0]/psi.dot(psi))
    assert slack == -sp.Rational(963, 1105000)
    checks["base_full_state_kernel_and_two_angular_derivatives"] = True
    checks["invalid_source_positive_and_preserves_eq2_eq3"] = all(
        sp.factor(q.subs(a, 2)) == 0 for q in invariants(shifted)[4:])
    assert checks["invalid_source_positive_and_preserves_eq2_eq3"]
    checks["invalid_source_normalized_slack_and_quartic_residual"] = True

    result = {
        "checks": checks,
        "signed_mixtures": witnesses,
        "invalid_source": {"a": 2, "base": list(map(str, base)),
                           "shifted": list(map(str, shifted)),
                           "normalized_slack": str(slack),
                           "quartic_residual": str(quartic_residual)},
        "branch": "Native sources select x_plus=(a+1)^2*t*(k+1)*theta; the squared equation also has x_minus.",
        "scope": "Exact displayed identities and witnesses. No converse, full conic hull, or standalone global-validity theorem.",
        "dependencies": "Python standard library and SymPy only; no private files or recovered source ZIP."
    }
    target = Path(__file__).with_suffix(".json")
    target.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(f"Verified {len(checks)} exact checks; wrote {target.name}.")
    print("Mixture determinants:", [v["determinant"] for v in witnesses])
    print("Invalid-source slack:", slack, "; quartic residual:", quartic_residual)


if __name__ == "__main__":
    main()
