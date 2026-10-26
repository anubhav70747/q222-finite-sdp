#!/usr/bin/env python3
"""Exact fresh checks for the literal Schur-remainder obstruction.

Requires SymPy 1.14.0. Run from any working directory. The proof of the
infinite-family conclusion is in schur_remainder_proof.tex; finite symbolic
checks are recorded here without being labelled a proof-assistant verification.
"""
from pathlib import Path
import hashlib
import json
import sympy as s

HERE = Path(__file__).resolve().parent
checks = []

def check(name, truth, detail=None):
    if not bool(truth):
        raise AssertionError(name)
    row = {"name": name, "status": "PASS"}
    if detail is not None:
        row["detail"] = detail
    checks.append(row)

def zero(name, expression):
    check(name, s.cancel(expression) == 0)

def terms(expr, gens):
    return [{"powers": list(p), "coefficient": str(c)}
            for p, c in s.Poly(expr, *gens, domain=s.QQ).terms()]

u, v, W, X, Y, Z, t, U, V = s.symbols("u v W X Y Z t U V")
R = s.Rational
Du, Dv = 1 + u**2, 1 + v**2
state = (W, X, Y, Z)
c0 = [R(1), R(7, 8), R(9, 8), R(1, 18), R(11, 8), R(9, 8), R(1, 18)]
d, n, m, C, A, D, alpha = c0

P0 = (n*Du*Dv*X**2 + m*(v*W-X+u*v*Y-u*Z)**2
      + C*Dv*(u*W-Y)**2 + A*Du*Dv*Y**2
      + (D+d/2)*(u*W+u*v*X-Y-v*Z)**2
      + (alpha+d/2)*Du*(v*Y-Z)**2
      + R(1, 2)*Dv*(X+u*Z)**2 - R(1, 2)*Du*Dv*Z**2)
P = s.expand(P0 + t*Du*Dv*X**2)
M = s.hessian(P, state)/2
M0 = M.subs(t, 0)
zero("state_matrix_recovers_full_parent", (s.Matrix(state).T*M*s.Matrix(state))[0] - P)
zero("complete_parent_increment", P-P0-t*Du*Dv*X**2)
check("seed_source_is_strictly_positive", all(x > 0 for x in c0))
check("seed_in_reciprocal_sector", c0[2] == c0[5] and c0[3] == c0[6])

rows = s.Matrix([[u, u*v, -1], [0, 0, v], [0, 1, 0]])
zero("upper_signed_core_independent_rows", rows.det()+u*v)
upper_core = R(1, 2)*((u*W+u*v*X-Y)**2 + Du*v**2*Y**2 + Dv*X**2)
zero("upper_signed_core_at_Z_zero", (P0-(n*Du*Dv*X**2+m*(v*W-X+u*v*Y-u*Z)**2+C*Dv*(u*W-Y)**2+A*Du*Dv*Y**2+D*(u*W+u*v*X-Y-v*Z)**2+alpha*Du*(v*Y-Z)**2)).subs(Z, 0)-upper_core)

Luv = s.cancel(M0.det(method="domain-ge")/(Du**2*Dv**2*u**2*v**2))
L0 = R(121, 331776)*(880*U + 9317*U/V + 720*V - 21994 + 880/V + 9477*V/U + 720/U)
zero("full_two_angle_determinant_Laurent_formula", Luv-L0.subs({U:u**2, V:v**2}))
outer = [880, 9317, 720, 880, 9477, 720]
check("six_outer_Laurent_masses_positive", all(x>0 for x in outer))
zero("Laurent_value_at_contact", L0.subs({U:1,V:1}))
zero("Laurent_log_gradient_U_at_contact", (U*s.diff(L0,U)).subs({U:1,V:1}))
zero("Laurent_log_gradient_V_at_contact", (V*s.diff(L0,V)).subs({U:1,V:1}))
for name, substitutions in [
    ("u_sign_gauge", {u:-u, Y:-Y, Z:-Z}),
    ("v_sign_gauge", {v:-v, X:-X, Z:-Z}),
]:
    zero(name, P.subs(substitutions, simultaneous=True)-P)
check("parent_can_be_bihomogenized_to_each_angle_degree_two", s.degree(P,u)<=2 and s.degree(P,v)<=2)
seed = {t:0,u:1,v:1,W:1,X:R(-1,10),Y:R(1,10),Z:1}
zero("seed_state_contact_value", P.subs(seed))
for x in (u,v,W,X,Y):
    zero("seed_contact_derivative_"+str(x), s.diff(P,x).subs(seed))

Apoly = 26883*U**2 + 105314*U + 29403
Bpoly = 29000*U**2 + 50384*U + 21384
Dpoly = Apoly + t*Bpoly
Npoly = 112167*(U-1)**2 + t*(121000*U**2-189648*U+81576)
Mline = M.subs(v,1)
upper_det = s.factor(Mline[:3,:3].det(method="domain-ge"))
full_det = s.factor(Mline.det(method="domain-ge"))
zero("literal_upper_determinant", upper_det-R(11,20736)*(u**2+1)*Dpoly.subs(U,u**2))
zero("literal_full_determinant", full_det-R(11,82944)*(u**2+1)**2*Npoly.subs(U,u**2))
Schur = (U+1)*Npoly/(4*Dpoly)
zero("literal_Schur_ratio", full_det/upper_det-Schur.subs(U,u**2))
resAB = s.resultant(Apoly,Bpoly,U)
check("fixed_polynomials_coprime", resAB == 1386911105418617856)
resDN = s.factor(s.resultant(Dpoly,Npoly,U))
zero("denominator_numerator_resultant", resDN-10441728*(8*t+11)*(1000*t+927)*(33677*t+55550)**2)
zero("no_factor_U_plus_one_in_denominator", Dpoly.subs(U,-1)+49028)
disc = s.discriminant(Dpoly,U)
zero("denominator_discriminant", disc-448*(129472*t**2+10942118*t+17699275))
check("denominator_coefficients_positive_for_nonnegative_t", all(all(a>0 for a in s.Poly(c,t).all_coeffs()) for c in s.Poly(Dpoly,U).all_coeffs()))

data = {
    "status": "EXACT_LITERAL_REMAINDER_OBSTRUCTION",
    "scope": "The literal scalar Schur remainders have infinite linear span. This does not rule out fixed features or a finite lift for the complete parent or for Q222.",
    "source_order": ["d","n","m","C","A","D","alpha"],
    "seed": [str(x) for x in c0],
    "source_family": "seed + t*(0,1,0,0,0,0,0), t >= 0",
    "state_order": list(map(str,state)),
    "full_parent_variable_order": ["u","v","W","X","Y","Z","t"],
    "full_parent": terms(P,(u,v,W,X,Y,Z,t)),
    "seed_Laurent": str(L0),
    "seed_Laurent_outer_integer_masses": outer,
    "seed_Laurent_prefactor": "121/331776",
    "seed_Laurent_constant_integer_mass": -21994,
    "restricted_angle": "v=1; U=u^2>0",
    "A": terms(Apoly,(U,)), "B": terms(Bpoly,(U,)),
    "D": terms(Dpoly,(U,t)), "N": terms(Npoly,(U,t)),
    "Schur": "(U+1)*N(U,t)/(4*D(U,t))",
    "resultant_A_B": str(resAB),
    "resultant_D_N": str(resDN),
    "D_at_minus_one": "-49028",
    "discriminant_D": str(s.factor(disc)),
    "denominator_pairwise_coprime_argument": "D_s-D_t=(s-t)B; a common divisor for s!=t would divide A and B.",
    "linear_independence_argument": "Multiply any finite rational identity by the product of its D_t, then reduce modulo each D_t. Its numerator and all the other denominators are units in that quotient, forcing its coefficient to vanish.",
    "complete_parent_fixed_features": "Four entries of M_0(u,v)^(1/2)*(W,X,Y,Z)^T and sqrt((1+u^2)*(1+v^2))*X. Source weights are 1,1,1,1,t.",
    "global_extraction_status": "OPEN"
}
data_path=HERE/"schur_remainder_data.json"
data_path.write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")
receipt={"status":"PASS","runtime":"SymPy 1.14.0 exact rational polynomial arithmetic","check_groups":len(checks),"checks":checks,"data_sha256":hashlib.sha256(data_path.read_bytes()).hexdigest(),"scope":"Fresh exact identities only; infinite-family and positivity conclusions use the accompanying analytical proof. No global extraction or historical R21 re-audit is claimed."}
(HERE/"schur_remainder_verification.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
print(json.dumps({"status":"PASS","check_groups":len(checks),"data_sha256":receipt["data_sha256"]},sort_keys=True))
