#!/usr/bin/env python3
"""Regenerate exact, local diagonal/B0/event-contact records.

Requires SymPy 1.14.0. All records are deterministic. The certificates are
local or module-level; this script does not certify a global Q222 lift.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as s

HERE = Path(__file__).resolve().parent
GROUPS: list[dict] = []


def require(ok: bool, name: str) -> None:
    if not ok:
        raise AssertionError(name)


def zero(expr) -> bool:
    return s.expand(expr) == 0


def rational(value) -> str:
    value = s.Rational(value)
    return str(value)


def polynomial(expr, variables) -> dict:
    p = s.Poly(s.expand(expr), *variables, domain=s.QQ)
    return {
        "variables": [str(v) for v in variables],
        "total_degree": int(p.total_degree()),
        "terms": [
            {"powers": list(powers), "coefficient": rational(coefficient)}
            for powers, coefficient in p.terms() if coefficient != 0
        ],
    }


def write_json(name: str, value: dict) -> None:
    (HERE / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def separated_identity() -> dict:
    # Three arbitrary separated summands verify the fully symbolic formula;
    # the accompanying proof gives the identical termwise argument for all n.
    n = 3
    alpha = s.symbols("alpha0:3")
    beta = s.symbols("beta0:3")
    aa = s.symbols("A0:3")
    bb = s.symbols("B0:3")
    delta = s.Symbol("delta")
    c = sum(alpha[i] * beta[i] for i in range(n))
    gamma = sum((alpha[i] + aa[i]) * (beta[i] + bb[i]) for i in range(n))
    p = (c - n * delta**2) / 2 + sum((beta[i] - delta) * aa[i] for i in range(n))
    q = (c - n * delta**2) / 2 + sum((alpha[i] - delta) * bb[i] for i in range(n))
    rhs = p + q + sum((aa[i] + delta) * (bb[i] + delta) for i in range(n))
    require(zero(gamma - rhs), "separated coefficient identity")
    x, y, l, k, a, b = s.symbols("x y l k a b")
    toy = (x + s.Rational(1, 4)) * (y + s.Rational(1, 4)) + s.Rational(15, 32) - x / 4 + s.Rational(15, 32) - y / 4
    require(zero(1 + x*y - toy), "toy separated coefficient identity")
    v = s.Matrix([1, k])
    gram = s.Matrix([[l*l, -l], [-l, 1]])
    require(zero(a*b * (v.T * gram * v)[0] - a*b*(k-l)**2), "PSD2 square identity")
    require(gram.det() == 0, "rank-one Gram determinant")
    GROUPS.append({"name": "separated positive coefficient identities", "status": "PASS", "checks": 4})
    return {
        "scope": "General theorem proved in diagonal_proofs.tex; symbolic corroboration and exact toy data only.",
        "generic_symbolic_summands_checked": 3,
        "positivity_condition": "3*n*delta^2 + 2*delta*max(sum(abs(alpha_i)),sum(abs(beta_i))) < c",
        "decomposition": "gamma = p(x) + q(y) + sum_i (A_i(x)+delta)*(B_i(y)+delta)",
        "toy": {
            "source_interval": ["-1/4", "1/4"],
            "test_interval": ["-1/4", "1/4"],
            "source_factors": ["x+1/4", "15/32-x/4", "1"],
            "testing_factors": ["y+1/4", "1", "15/32-y/4"],
            "pure_factor_lower_bound": "13/32",
            "psd_block_orders": [2, 2, 2],
        },
    }


def frame_identity() -> dict:
    m = 5
    d = s.symbols("d0:5")
    eps = s.Symbol("eps")
    h = s.Matrix(m, m, lambda i, j: s.Symbol(f"h{min(i,j)}{max(i,j)}"))
    lhs = (s.Matrix(d).T * h * s.Matrix(d))[0]
    rhs = sum((h[i,i] - 2*(m-1)*eps)*d[i]**2 for i in range(m))
    for i in range(m):
        for j in range(i+1,m):
            rhs += (eps+h[i,j]/2)*(d[i]+d[j])**2
            rhs += (eps-h[i,j]/2)*(d[i]-d[j])**2
    require(zero(lhs-rhs), "fixed frame identity")
    eps0 = s.Rational(1,4*m)
    require(1-2*(m-1)*eps0 > 0 and eps0 > 0, "positive frame anchors")
    # A concrete singular-normal example checks the addition's exact
    # invisibility on the diagonal quotient without claiming it proves the
    # general smooth-algebra statement.
    dx, dy, t = s.symbols("dx dy t")
    base = s.Matrix([[2, 3], [3, -1]])
    normal = s.Matrix([[0, 1]])
    corrected = base + t*normal.T*normal
    difference = (s.Matrix([dx,dy]).T*(corrected-base)*s.Matrix([dx,dy]))[0]
    require(zero(difference-t*dy**2), "normal-square addition")
    require(base[0,0] > 0 and corrected.subs(t,7).det() > 0, "concrete positive correction")
    GROUPS.append({"name": "denominator-free frame corroboration", "status": "PASS", "checks": 4})
    return {
        "scope": "General smooth-diagonal theorem in diagonal_proofs.tex, assuming exact J^2 data and positive tangent residue.",
        "generic_frame_dimension_checked": m,
        "frame_count": m*m,
        "epsilon": rational(eps0),
        "anchor_axis_coefficient": rational(1-2*(m-1)*eps0),
        "anchor_pair_coefficient": rational(eps0),
        "mixed_inverses_allowed": False,
        "global_finite_cover_supplied": False,
    }


def b0_module() -> dict:
    U, V = s.symbols("U V")
    z = s.Matrix(s.symbols("z0:4"))
    c = s.symbols("c1:5")
    e = [s.eye(4)[:,i] for i in range(4)]
    omega = s.Matrix([U,-1,0,-U*V])
    matrix = (
        c[0]*(1+U**2)*e[0]*e[0].T
        + c[1]*(1+V**2)*omega*omega.T
        + c[2]*(1+U**2)*(1+V**2)*e[1]*e[1].T
        + c[3]*U**2*e[2]*e[2].T
    )
    features = [
        z[0], U*z[0], U*z[0]-z[1]-U*V*z[3],
        U*V*z[0]-V*z[1]-U*V**2*z[3], z[1], U*z[1],
        V*z[1], U*V*z[1], U*z[2],
    ]
    weights = [c[0],c[0],c[1],c[1],c[2],c[2],c[2],c[2],c[3]]
    require(zero((z.T*matrix*z)[0] - sum(w*f*f for w,f in zip(weights,features))), "B0 nine-square identity")
    tests = [
        (0,0,e[0],c[0]),
        (1,1,e[3],2*c[1]),
        (1,1,s.Matrix([0,1,0,-1]),4*c[2]),
        (1,0,e[2],c[3]),
    ]
    test_records = []
    for u0,v0,vector,expected in tests:
        actual = s.expand((vector.T*matrix.subs({U:u0,V:v0})*vector)[0])
        require(actual == expected, "B0 necessary coefficient test")
        test_records.append({"U":u0,"V":v0,"state":[int(vv) for vv in vector],"value":str(actual)})
    T,L,A,eps0 = s.symbols("T L A eps0")
    p = 1-A-L*T*T
    w = 1-L*T*T
    hh = 1-L*(1+T*T)
    raw = [A*T*(1+T*T)*p/w, A*eps0, eps0*p, A*eps0*(1+T*T)*w*p/hh]
    chamber = {T:1,L:s.Rational(2,5),A:s.Rational(1,5),eps0:1}
    sample = [s.cancel(cc.subs(chamber)) for cc in raw]
    require(sample == [s.Rational(4,15),s.Rational(1,5),s.Rational(2,5),s.Rational(12,25)], "B0 rational chamber")
    # Each feature's coefficient in z is a linear combination of the five
    # fixed monomials. This is the exact premise for source-only congruences.
    monomial_set = {(0,0),(1,0),(0,1),(1,1),(1,2)}
    feature_coefficient_data = []
    for feature in features:
        coefficients = [s.diff(feature,zz) for zz in z]
        require(zero(feature-sum(coefficient*zz for coefficient,zz in zip(coefficients,z))), "B0 state linearity")
        for coefficient in coefficients:
            if coefficient != 0:
                require(set(s.Poly(coefficient,U,V).monoms()).issubset(monomial_set), "B0 congruence envelope")
        feature_coefficient_data.append([polynomial(coefficient,(U,V)) for coefficient in coefficients])
    matrix_records = []
    for i in range(4):
        for j in range(i,4):
            if matrix[i,j] != 0:
                matrix_records.append({"row":i,"column":j,"polynomial":polynomial(matrix[i,j],(U,V,*c))})
    GROUPS.append({"name": "B0 displayed coefficient module", "status": "PASS", "checks": 15,
                   "check_count_convention": "1 square identity + 4 tests + 1 chamber + 9 state/envelope groups"})
    return {
        "scope": "Exactly the displayed four-coefficient B0 module. No restoration to the entire original eight-coordinate quantum slack is supplied.",
        "coefficient_order": [str(cc) for cc in c],
        "matrix_upper_entries": matrix_records,
        "feature_coefficients_in_state_order_z0_z1_z2_z3": feature_coefficient_data,
        "weights": [str(ww) for ww in weights],
        "necessary_tests": test_records,
        "coefficient_cone": "c1>=0, c2>=0, c3>=0, c4>=0",
        "coefficient_cone_pencil": {"order":4,"auxiliary_variables":0,"matrix":"diag(c1,c2,c3,c4)"},
        "sample_chamber_T_L_A_eps0": ["1","2/5","1/5","1"],
        "sample_p_w_h": ["2/5","3/5","1/5"],
        "sample_coefficients": [rational(vv) for vv in sample],
        "congruence_envelope_monomials": [[0,0],[1,0],[0,1],[1,1],[1,2]],
        "congruence_envelope_dimension":20,
        "congruence_condition": "The 4x4 congruence matrix depends only on the coefficient/point factor, not on the polynomial-testing variables U,V in the support factor.",
        "congruence_preserves_order_two_blocks_claimed":False,
        "full_B0_restoration_supplied":False,
    }


def event_contact() -> dict:
    u,v,a,b = s.symbols("u v a b")
    variables = (u,v,a,b)
    Du,Dv = 1+u*u,1+v*v
    y1,y2,y3 = s.symbols("y1 y2 y3")
    psi = s.Matrix([-v*y3,y3,y1,y2])
    rows = [s.Matrix(vector) for vector in (
        (0,1,0,0),(v,-1,u*v,-u),(u,0,-1,0),
        (0,0,1,0),(u,u*v,-1,-v),
    )]
    s2,s3,zstate = s.Matrix([0,1,0,u]),s.Matrix([0,0,v,-1]),s.Matrix([0,0,0,1])
    dot = lambda vector: (vector.T*psi)[0]
    restricted = [
        2*(dot(rows[4])**2+Dv*dot(s2)**2+Du*dot(s3)**2-Du*Dv*dot(zstate)**2),
        4*Du*Dv*dot(rows[0])**2,
        4*dot(rows[1])**2,
        4*Dv*dot(rows[2])**2,
        4*Du*Dv*dot(rows[3])**2,
        4*dot(rows[4])**2,
    ]
    full_forms = [
        2*((Dv+u*u*v*v)*y1*y1-2*u*u*v*y1*y2+u*u*y2*y2+Dv*y3*y3+2*u*Dv*y2*y3),
        4*Du*Dv*y3*y3,
        4*(u*v*y1-u*y2-Dv*y3)**2,
        4*Dv*(y1+u*v*y3)**2,
        4*Du*Dv*y1*y1,
        4*(y1+v*y2)**2,
    ]
    require(all(zero(lhs-rhs) for lhs,rhs in zip(restricted,full_forms)), "six physical event restrictions")
    u0,u1,v0,v1 = s.symbols("u0 u1 v0 v1")
    homogeneous = [
        (u0*u0+u1*u1)*v0*v0*(v0*v0+v1*v1)*y3*y3,
        (u1*v1*v0*y1-u1*v0*v0*y2-u0*(v0*v0+v1*v1)*y3)**2,
        (v0*v0+v1*v1)*(u0*v0*y1+u1*v1*y3)**2,
        (u0*u0+u1*u1)*v0*v0*(v0*v0+v1*v1)*y1*y1,
        u0*u0*v0*v0*(v0*y1+v1*y2)**2,
    ]
    require(all(s.cancel(u0*u0*v0**4*form.subs({u:u1/u0,v:v1/v0})-4*hom) == 0
                for form,hom in zip(full_forms[1:],homogeneous)), "five exact projective event homogenizations")
    G = [
        2*((Dv+u*u*v*v)*a*a-2*u*u*v*a+u*u+Dv*b*b+2*u*Dv*b),
        4*Du*Dv*b*b,
        4*(u*v*a-u-Dv*b)**2,
        4*Dv*(a+u*v*b)**2,
        4*Du*Dv*a*a,
        4*(a+v)**2,
    ]
    G = [s.expand(gg) for gg in G]
    require(all(zero(form.subs({y1:a,y2:1,y3:b})-gg) for form,gg in zip(full_forms,G)), "event state-chart substitution")
    J = s.Matrix([G]+[[s.diff(gg,variable) for gg in G] for variable in variables])
    # Domain-GE keeps all arithmetic in QQ[u,v,a,b] and is significantly
    # quicker than fully expanding six permutation determinants.
    cofactors = [s.expand((-1)**j * J[:,[k for k in range(6) if k != j]].det(method="domain-ge")) for j in range(6)]
    require(any(cofactor != 0 for cofactor in cofactors), "event contact nonzero cofactor vector")
    common = s.Poly(cofactors[0],*variables,domain=s.ZZ)
    for cofactor in cofactors[1:]:
        common = s.gcd(common,s.Poly(cofactor,*variables,domain=s.ZZ))
    common_expr = common.as_expr()
    primitive = [s.div(cofactor,common_expr,*variables)[0] for cofactor in cofactors]
    for cofactor,component in zip(cofactors,primitive):
        require(zero(cofactor-common_expr*component), "event cofactor exact common-factor division")
    residual = J*s.Matrix(primitive)
    require(all(zero(entry) for entry in residual), "event stationary-contact kernel")
    degrees = [int(s.Poly(component,*variables).total_degree()) for component in primitive]
    # This expected degree tuple is a reconstruction target, not a substitute
    # for the exact determinant, division, and five polynomial kernel checks.
    require(degrees == [15,13,13,11,15,17], "event primitive degree tuple")
    GROUPS.append({"name": "event stationary-contact cofactors", "status": "PASS", "checks": 26,
                   "check_count_convention": "6 physical restrictions + 5 homogenizations + 1 chart + 1 nonzero vector + 6 divisions + 5 kernel identities + 1 degree tuple + 1 primitive gcd"})
    primitive_gcd = s.Poly(primitive[0],*variables,domain=s.ZZ)
    for component in primitive[1:]:
        primitive_gcd = s.gcd(primitive_gcd,s.Poly(component,*variables,domain=s.ZZ))
    require(primitive_gcd.as_expr() == 1, "event primitive gcd")
    return {
        "scope": "Stationarity of the reduced event-contact equations only. No source validity, positivity, full state-kernel relation, or global chart coverage is inferred.",
        "variable_order": [str(variable) for variable in variables],
        "source_order": ["d","e1","e2","e3","e4","e5"],
        "form_order": ["Gd","G1","G2","G3","G4","G5"],
        "state_chart": "psi=(-v*b,b,a,1)",
        "full_state_order": ["-v*y3","y3","y1","y2"],
        "physical_restriction_source": "Canonical August30 base equation (4.2), lambda packet killed, alpha=beta=0",
        "projective_homogenization_source": "Canonical August30 base equation (11.0)",
        "homogenizing_multiplier": "u0^2*v0^4",
        "homogenized_forms_are_four_times_equation_11_0": True,
        "row_order": ["G","dG/du","dG/dv","dG/da","dG/db"],
        "forms": [polynomial(gg,variables) for gg in G],
        "jacobian_rows": [[polynomial(J[i,j],variables) for j in range(6)] for i in range(5)],
        "signed_cofactor_convention": "X_j=(-1)^j det(J with column j removed), zero-based j",
        "common_factor": polynomial(common_expr,variables),
        "common_factor_factored": str(s.factor(common_expr)),
        "primitive_cofactor_components": [polynomial(component,variables) for component in primitive],
        "primitive_cofactor_degrees": degrees,
        "primitive_cofactor_term_counts": [len(s.Poly(component,*variables).terms()) for component in primitive],
        "kernel_identity": "J*X=0 exactly over Q[u,v,a,b]",
    }


def main() -> None:
    require(s.__version__ == "1.14.0", "tested SymPy version")
    separated = separated_identity()
    frame = frame_identity()
    b0 = b0_module()
    event = event_contact()
    outputs = {
        "SEPARATED_COEFFICIENT_EXACT.json": separated,
        "DIAGONAL_FRAME_EXACT.json": frame,
        "B0_EXACT.json": b0,
        "EVENT_CONTACT_EXACT.json": event,
    }
    for name,value in outputs.items():
        write_json(name,value)
    receipt = {
        "status": "PASS",
        "sympy_version": s.__version__,
        "groups": GROUPS,
        "counts_are_overlapping_algebraic_corroboration_not_a_proof_count": True,
        "global_finite_lift_extraction_complete": False,
        "full_B0_restoration_complete": False,
        "data_sha256": {name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in outputs},
    }
    write_json("DIAGONAL_VERIFICATION.json",receipt)
    print(json.dumps({"status":"PASS","groups":len(GROUPS),"outputs":list(outputs),"global_extraction_complete":False},sort_keys=True))


if __name__ == "__main__":
    main()
