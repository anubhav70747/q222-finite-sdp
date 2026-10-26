#!/usr/bin/env python3
"""Regenerate and exactly check the tilted-CHSH projection extraction.

Run with Python >=3.9. Standard library only. This verifies the listed finite
algebraic identities, exact matrix coefficients and scope witness. The complete
analytical proof and inequalities are in tilted_family_exact_extraction.tex.
It does not verify an exact lift of the full eight-dimensional quantum set.
"""
import hashlib
import itertools
import json
from fractions import Fraction as F
from pathlib import Path
from exact_poly import P, expectation, kron, mat_add, mat_mul, mat_scale

HERE = Path(__file__).resolve().parent
COORDINATES = ["a0", "a1", "b0", "b1", "E00", "E01", "E10", "E11"]
checks = []


def check(name, statement, detail=None):
    if not statement:
        raise AssertionError(name)
    item = {"name": name, "status": "PASS"}
    if detail is not None:
        item["detail"] = detail
    checks.append(item)


def write(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def reduced_word(word):
    result = []
    for x in word:
        if result and result[-1] == x:
            result.pop()
        else:
            result.append(x)
    return tuple(result)


def nc_add(*items):
    ans = {}
    for item in items:
        for key, coefficient in item.items():
            ans[key] = ans.get(key, P()) + coefficient
    return {key: value for key, value in ans.items() if value}


def nc_scale(scalar, item):
    return {key: scalar * value for key, value in item.items() if scalar * value}


def nc_mul(a, b):
    ans = {}
    for (aw, bw), coefficient in a.items():
        for (av, bv), other in b.items():
            key = (reduced_word(aw + av), reduced_word(bw + bv))
            ans[key] = ans.get(key, P()) + coefficient * other
    return {key: value for key, value in ans.items() if value}


alpha, gamma, rho, eps, t, s, y = map(P.var, ["alpha", "gamma", "rho", "eps", "t", "s", "y"])
I = {((), ()): P(1)}
A0, A1 = {((0,), ()): P(1)}, {((1,), ()): P(1)}
B0, B1 = {((), (0,)): P(1)}, {((), (1,)): P(1)}
S = nc_add(nc_mul(A0, nc_add(B0, B1)), nc_mul(A1, nc_add(B0, nc_scale(-1, B1))))
Sp = nc_add(nc_mul(A0, nc_add(B0, nc_scale(-1, B1))), nc_mul(A1, nc_add(B0, B1)))
T = nc_add(nc_scale(alpha, A0), S)
R = nc_add(nc_scale(alpha, A1), nc_scale(-1, Sp))
check("reflection_identity", nc_add(nc_mul(T, T), nc_mul(R, R)) == nc_scale(8 + 2 * alpha**2, I),
      "Exact word reduction uses A0^2=A1^2=B0^2=B1^2=I and cross-party commutation only.")
beta = P.var("beta")
slack = nc_add(nc_scale(beta, I), nc_scale(-1, T))
residual = nc_add(nc_mul(slack, slack), nc_mul(R, R), nc_scale(-2 * beta, slack))
check("reflection_sos_identity", all(not p.reduce_square("beta", 8 + 2 * alpha**2)
                                      for p in residual.values()))

# Independently compute the attaining behavior using real 4x4 Pauli matrices.
Z, X, Id = [[1, 0], [0, -1]], [[0, 1], [1, 0]], [[1, 0], [0, 1]]
p, z, r = map(P.var, ["p", "z", "r"])
psi = [p, P(), P(), z]
Bplus = mat_add(Z, mat_scale(r, X))       # h B0
Bminus = mat_add(Z, mat_scale(-r, X))     # h B1
ops = [kron(Z, Id), kron(X, Id), kron(Id, Bplus), kron(Id, Bminus),
       kron(Z, Bplus), kron(Z, Bminus), kron(X, Bplus), kron(X, Bminus)]
expected = [p*p-z*z, P(), p*p-z*z, p*p-z*z, p*p+z*z, p*p+z*z, 2*p*z*r, -2*p*z*r]
for coordinate, op, target in zip(COORDINATES, ops, expected):
    check("attainment_pauli_" + coordinate, expectation(psi, op) == target)
check("attainment_parameter_relation", ((p*p-z*z)**2 + (2*p*z)**2 - 1)
      .reduce_square("p", 1-z*z) == 0)
check("attainment_observable_reflection", mat_mul(Bplus, Bplus) == mat_scale(1+r*r, Id)
      and mat_mul(Bminus, Bminus) == mat_scale(1+r*r, Id))
check("attainment_chsh", expected[4] + expected[5] + expected[6] - expected[7]
      == 2*(p*p+z*z) + 4*p*z*r)
check("attainment_endpoint_parameters", (1-t*t + 1) == 2-t*t)

# Source/testing products are kept separated before substitution. The relations
# below are the exact radical identities justified analytically in the proof.
v00, v11, v01 = (gamma-alpha)/2, (gamma+alpha)/2, -eps
source_square_norm = v00*v11 - 1
check("source_radical_product", source_square_norm.reduce_square("gamma", 4+alpha*alpha) == 0)
check("source_square_difference", v00-v11 == -alpha)
check("source_square_sum", v00+v11 == gamma)
expanded = v00*(rho+t) + 2*v01*(s/2) + v11*y + v11*(rho-t-y)
check("separated_tensor_sos", expanded == rho*gamma-alpha*t-eps*s)
check("testing_radicand_clear_denominator",
      ((rho+t)*(rho-t)-s*s/4).reduce_square("rho", 2) == 2-t*t-s*s/4)
check("constant_square_certificate", ((1-t)+(1+t))/2 == 1)
absalpha, sigma = P.var("absalpha"), P.var("sigma")
outer = 4-2*sigma*t-eps*s + (absalpha-2)*(1-sigma*t)
check("outer_tilt_endpoint_completion", outer == 2+absalpha-absalpha*sigma*t-eps*s)
check("support_stationary_boundary", (2*alpha*alpha+8-2*gamma*gamma)
      .reduce_square("gamma", 4+alpha*alpha) == 0)

# Exact COO convention: each entry means the actual symmetric matrix entry,
# so an off-diagonal record is copied to both positions without a factor of 2.
def algebraic(a=0, b=0):
    return {"rational": str(F(a)), "sqrt2": str(F(b))}


def record(row, column, a=0, b=0):
    return {"row": row, "column": column, "value": algebraic(a, b)}


def make_pencil(eliminated):
    names = ["constant"] + COORDINATES + ([] if eliminated else ["y"])
    mats = {name: [] for name in names}
    if eliminated:
        mats["constant"] = [record(0, 0, b=1), record(1, 1, b=1), record(2, 2, 1), record(3, 3, 1)]
        mats["a0"] = [record(0, 0, 1), record(1, 1, -1), record(2, 2, -1), record(3, 3, 1)]
        n, blocks = 4, [2, 1, 1]
    else:
        mats["constant"] = [record(0, 0, b=1), record(2, 2, b=1), record(3, 3, 1), record(4, 4, 1)]
        mats["a0"] = [record(0, 0, 1), record(2, 2, -1), record(3, 3, -1), record(4, 4, 1)]
        mats["y"] = [record(1, 1, 1), record(2, 2, -1)]
        n, blocks = 5, [2, 1, 1, 1]
    for name, coefficient in zip(COORDINATES[4:], [F(1, 2), F(1, 2), F(1, 2), F(-1, 2)]):
        mats[name] = [record(0, 1, coefficient)]
    return {
        "schema": "exact-symmetric-upper-coo-v1", "status": "EXACT_FOR_TILTED_PROJECTION_ONLY",
        "scope": "Inverse image of T under t=a0, s=E00+E01+E10-E11; not the full Q222 set.",
        "coordinate_order": COORDINATES, "matrix_order_N": n,
        "auxiliary_dimension_m": 0 if eliminated else 1, "block_orders": blocks,
        "index_base": 0, "coefficient_field": "Q(sqrt(2))", "sqrt2_root": "positive root of X^2-2",
        "coo_semantics": "Upper records are actual entries. Mirror each off-diagonal entry; do not double it.",
        "upper_record_count": sum(len(v) for v in mats.values()), "matrices": mats,
    }


def pencil_polynomial(data):
    n = data["matrix_order_N"]
    a = [[P() for _ in range(n)] for _ in range(n)]
    for name, records in data["matrices"].items():
        factor = P(1) if name == "constant" else P.var(name)
        for item in records:
            i, j, value = item["row"], item["column"], item["value"]
            coefficient = F(value["rational"]) + F(value["sqrt2"])*rho
            a[i][j] += factor*coefficient
            if i != j:
                a[j][i] += factor*coefficient
    return a


lift, direct = make_pencil(False), make_pencil(True)
check("lift_upper_coo_count", lift["upper_record_count"] == 14)
check("direct_upper_coo_count", direct["upper_record_count"] == 12)
lp, dp = pencil_polynomial(lift), pencil_polynomial(direct)
a0, e00, e01, e10, e11 = map(P.var, ["a0", "E00", "E01", "E10", "E11"])
chsh = e00+e01+e10-e11
desired_lift = [[rho+a0, chsh/2, 0, 0, 0], [chsh/2, y, 0, 0, 0],
                [0, 0, rho-a0-y, 0, 0], [0, 0, 0, 1-a0, 0], [0, 0, 0, 0, 1+a0]]
desired_direct = [[rho+a0, chsh/2, 0, 0], [chsh/2, rho-a0, 0, 0],
                  [0, 0, 1-a0, 0], [0, 0, 0, 1+a0]]
check("lift_original_coordinates", lp == desired_lift)
check("direct_original_coordinates", dp == desired_direct)
check("eliminated_determinant", (dp[0][0]*dp[1][1]-dp[0][1]*dp[1][0]).reduce_square("rho", 2)
      == 2-a0*a0-chsh*chsh/4)
check("schur_complement_elimination", ((rho+t)*(rho-t)-s*s/4).reduce_square("rho", 2)
      == 2-t*t-s*s/4)

# A rational behavior that satisfies every relabeled tilted family and all
# probability inequalities but violates a necessary quantum equality.
qstar = [F(0)]*4 + [F(1), F(1), F(1, 4), F(-1, 4)]
probabilities = []
for x, yy, a, b in itertools.product(range(2), range(2), [-1, 1], [-1, 1]):
    prob = (1+a*qstar[x]+b*qstar[2+yy]+a*b*qstar[4+2*x+yy])/4
    probabilities.append({"x": x, "y": yy, "a": a, "b": b, "value": str(prob)})
check("scope_probabilities_nonnegative", all(F(item["value"]) >= 0 for item in probabilities))
for x, yy in itertools.product(range(2), range(2)):
    check(f"scope_probabilities_normalized_{x}{yy}",
          sum(F(item["value"]) for item in probabilities if item["x"] == x and item["y"] == yy) == 1)
patterns = []
for signs in itertools.product([-1, 1], repeat=4):
    if signs[0]*signs[1]*signs[2]*signs[3] == -1:
        value = sum(sgn*v for sgn, v in zip(signs, qstar[4:]))
        patterns.append({"signs": list(signs), "value": str(value)})
max_chsh = max(abs(F(item["value"])) for item in patterns)
check("scope_all_chsh_sign_patterns", len(patterns) == 8 and max_chsh == F(5, 2))
check("scope_below_tsirelson", max_chsh*max_chsh < 8)
check("scope_all_marginals_zero", all(value == 0 for value in qstar[:4]))
check("scope_nonquantum_equality_failure", qstar[4] == qstar[5] == 1 and qstar[6] != qstar[7],
      "E00=E01=1 forces (B0-B1)|psi>=0, hence E10=E11 in any quantum realization.")
check("scope_lift_feasible", qstar[0]*qstar[0] + (qstar[4]+qstar[5]+qstar[6]-qstar[7])**2/4 < 2)

features = {
    "status": "EXACT_FOR_TILTED_PROJECTION_ONLY", "domain": "|t|<=1 and t^2+s^2/4<=2",
    "fixed_testing_functions": {
        "u1": "sqrt(sqrt(2)+t)", "u2": "s/(2*u1)",
        "u3": "sqrt(sqrt(2)-t-u2^2)", "u4": "sqrt(1-t)", "u5": "sqrt(1+t)"},
    "radical_convention": "Unique nonnegative square root in the testing real closed field",
    "nonvanishing_denominator": "sqrt(2)+t >= sqrt(2)-1 > 0",
    "square_spaces": [["u1", "u2"], ["u3"], ["u4"], ["u5"]],
    "product_coordinate_order": ["1", "t", "s", "y"], "auxiliary_definition": "y=u2^2",
    "products": {"u1^2": "sqrt(2)+t", "u1*u2": "s/2", "u2^2": "y",
                 "u3^2": "sqrt(2)-t-y", "u4^2": "1-t", "u5^2": "1+t"},
    "constant_certificate": "1=(u4^2+u5^2)/2",
}
source = {
    "status": "EXACT_FOR_TILTED_PROJECTION_ONLY", "parameters": "alpha real, epsilon in {-1,+1}",
    "source_functions": {"gamma": "sqrt(4+alpha^2)", "beta": "sqrt(2)*gamma",
                         "v0": "sqrt((gamma-alpha)/2)", "v1": "-epsilon*sqrt((gamma+alpha)/2)"},
    "source_relations": ["v0^2+v1^2=gamma", "v0^2-v1^2=-alpha", "v0*v1=-epsilon"],
    "tensor_identity": "beta tensor 1-alpha tensor t-1 tensor epsilon*s=(v0 tensor u1+v1 tensor u2)^2+(v1 tensor u3)^2",
    "tensor_scope": "Kc tensor_R Kt; arbitrary real closed extensions, separately labelled, no disjointness assumption",
    "optimal_support": {"abs(alpha)<=2": "sqrt(8+2*alpha^2)", "abs(alpha)>=2": "2+abs(alpha)"},
    "outer_completion": "2+abs(alpha)-alpha*t-epsilon*s=(4-2*sign(alpha)*t-epsilon*s)+(abs(alpha)-2)*(1-sign(alpha)*t)",
    "endpoint_status": "alpha=0, alpha=+/-2, |alpha|>2, t=+/-1 all covered",
    "general_affine": "For c2!=0 use lambda=abs(c2), alpha=-c1/lambda, epsilon=-sign(c2), plus nonnegative constant. For c2=0 use marginal squares.",
}
witness = {
    "status": "EXACT_SCOPE_COUNTEREXAMPLE", "coordinate_order": COORDINATES, "q": list(map(str, qstar)),
    "probabilities": probabilities, "odd_chsh_sign_patterns": patterns, "maximum_absolute_chsh": str(max_chsh),
    "all_relabeled_tilted_directions": "All four marginals vanish. Every relabeled CHSH has absolute value <=5/2<sqrt(8), and h(alpha)>=sqrt(8).",
    "nonquantum_proof": ["Purify and projectively dilate any hypothetical quantum realization.",
                         "E00=E01=1 implies (A0-B0)|psi>=(A0-B1)|psi>=0.",
                         "Then (B0-B1)|psi>=0 and E10-E11=<psi|A1(B0-B1)|psi>=0.",
                         "The witness has E10-E11=1/2, a contradiction."],
    "consequence": "Tilted-family extraction, including all standard relabelings and probability inequalities, is not full Q222 extraction.",
}
write("TILTED_FIXED_FEATURES.json", features)
write("TILTED_SOURCE_CERTIFICATES.json", source)
write("TILTED_ORIGINAL_COORDINATE_LIFT.coo.json", lift)
write("TILTED_DIRECT_ELIMINATED.coo.json", direct)
write("TILTED_SCOPE_COUNTEREXAMPLE.json", witness)
receipt = {"status": "PASS", "run_kind": "fresh exact rational symbolic regeneration",
           "arithmetic": "Python integers and fractions.Fraction; sparse polynomials; reduced noncommutative words",
           "scope": "Listed finite algebraic identities, exact matrices and rational scope witness; see analytical text for inequalities and quantum reduction.",
           "not_verified": ["full Q222 finite-lift existence", "full Q222 exact extraction", "quantifier elimination"],
           "check_count": len(checks), "checks": checks}
write("CHECK_TILTED_EXTRACTION.json", receipt)
names = ["exact_poly.py", "build_and_check_tilted_extraction.py", "TILTED_FIXED_FEATURES.json",
         "TILTED_SOURCE_CERTIFICATES.json", "TILTED_ORIGINAL_COORDINATE_LIFT.coo.json",
         "TILTED_DIRECT_ELIMINATED.coo.json", "TILTED_SCOPE_COUNTEREXAMPLE.json", "CHECK_TILTED_EXTRACTION.json"]
write("TILTED_MANIFEST.json", {
    "schema": "sha256-v1", "scope": "This newly regenerated tilted extraction component only",
    "files": {name: {"sha256": hashlib.sha256((HERE/name).read_bytes()).hexdigest(),
                     "bytes": (HERE/name).stat().st_size} for name in names}})
print(json.dumps({"status": "PASS", "check_count": len(checks), "matrix_orders": [5, 4],
                  "auxiliary_dimensions": [1, 0], "upper_coo_records": [14, 12],
                  "scope": "tilted projection only"}, sort_keys=True))
