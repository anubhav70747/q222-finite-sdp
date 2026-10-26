#!/usr/bin/env python3
"""Separate-agent review of the regular-C2 data using direct SymPy derivatives.

Does not import either C2 implementation. This complements, rather than
replaces, the independent full Taylor-allocation checker supplied in c2/.
"""
from __future__ import annotations
import hashlib
import json
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
C2 = ROOT / "c2"


def poly_dict(records):
    return {tuple(powers):F(coefficient) for powers,coefficient in records}


def add_to(output,key,value):
    output[key] = output.get(key,F(0))+value
    if not output[key]:
        del output[key]


def multiply(a,b):
    output = {}
    for e,c in a.items():
        for f,d in b.items():
            add_to(output,tuple(x+y for x,y in zip(e,f)),c*d)
    return output


def direct_parent(source,z,Z=1):
    d,n,m,C,A,D,alpha = source
    u,v,W,X,Y = z
    Du,Dv = 1+u*u,1+v*v
    return (n*Du*Dv*X*X+m*(v*W-X+u*v*Y-u*Z)**2
            +C*Dv*(u*W-Y)**2+A*Du*Dv*Y*Y
            +(D+d/2)*(u*W+u*v*X-Y-v*Z)**2
            +(alpha+d/2)*Du*(v*Y-Z)**2
            +(d/2)*Dv*(X+u*Z)**2-(d/2)*Du*Dv*Z*Z)


def sympy_records(expression,variables):
    return {powers:F(str(value)) for powers,value in s.Poly(s.expand(expression),*variables).terms() if value}


def main():
    frame = json.loads((C2/"regular_c2_frames.json").read_text())
    features = json.loads((C2/"regular_c2_feature_list.json").read_text())
    products = json.loads((C2/"regular_c2_product_map.json").read_text())
    rho,eps = F(frame["radius"]),F(frame["epsilon"])
    bound = F(frame["uniform_positive_remainder_lower_bound"])
    c = s.symbols("n m C A D alpha")
    z = s.symbols("u v W X Y")
    ca = list(map(s.Rational,frame["source_anchor"][1:]))
    za = list(map(s.Rational,frame["contact_anchor"]))
    anchor = dict(zip(c+z,ca+za))
    P = direct_parent([s.Integer(1),*c],z)
    contact = s.Matrix([P,*[s.diff(P,variable) for variable in z]])
    assert all(s.cancel(entry.subs(anchor)) == 0 for entry in contact)
    jac = contact.jacobian(c).subs(anchor)
    assert jac == s.Matrix([[s.Rational(value) for value in row] for row in frame["contact_jacobian_at_anchor"]])
    assert jac.det() == s.Rational(114048,78125)
    hessian = s.hessian(P,z).subs(anchor)/2
    assert hessian == s.Matrix([[s.Rational(value) for value in row] for row in frame["H0"]])
    L = s.Matrix([[s.Rational(value) for value in row] for row in frame["LDL_L"]])
    diagonal = s.diag(*map(s.Rational,frame["LDL_D"]))
    assert hessian == L*diagonal*L.T
    assert all(value > 0 for value in diagonal.diagonal())
    # Independently expand the source Jacobian on the contact box.
    h = s.symbols("h0:5")
    J_h = contact.jacobian(c).subs(dict(zip(z,[value+variable for value,variable in zip(za,h)])))
    E = jac.inv()*J_h-s.eye(6)
    row_l1 = [sum(abs(F(str(value))) for j in range(6) for powers,value in s.Poly(s.expand(E[i,j]),*h).terms() if any(powers)) for i in range(6)]
    assert max(row_l1) == F(114199,160)
    assert max(row_l1)*rho == F(11650604933999,11225797109143808) < 1
    # Check every separated coefficient and its independently computed bound.
    all_norms = []
    mixed_count = 0
    min_pure_bound = None
    for ff in frame["frames"]:
        p = poly_dict(ff["positive_source_remainder"])
        q = poly_dict(ff["positive_testing_remainder"])
        rebuilt = {}
        for powers,value in p.items(): add_to(rebuilt,powers+(0,)*5,value)
        for powers,value in q.items(): add_to(rebuilt,(0,)*11+powers,value)
        for pure in (p,q):
            constant = pure.get((0,)*len(next(iter(pure))),F(0))
            lower = constant-sum(abs(value)*rho**sum(powers) for powers,value in pure.items() if any(powers))
            assert lower >= bound
            min_pure_bound = lower if min_pure_bound is None else min(min_pure_bound,lower)
        for mm in ff["mixed_positive_products"]:
            aa = poly_dict(mm["source_constant_free_polynomial"])
            U = poly_dict(mm["source_polynomial"])
            assert all(any(powers) for powers in aa)
            sigma = sum(map(abs,aa.values()))
            assert sigma == F(mm["source_coefficient_l1"])
            assert F(mm["source_shift"]) == sigma*rho
            beta = tuple(mm["testing_monomial"])
            shift = F(mm["testing_shift"])
            assert sum(beta) > 0 and shift == rho**sum(beta)
            assert U == {**aa,(0,)*11:sigma*rho}
            for powers,value in U.items():
                add_to(rebuilt,powers+beta,value)
                add_to(rebuilt,powers+(0,)*5,value*shift)
            mixed_count += 1
        expected = {tuple(a)+tuple(b):F(value) for a,b,value in ff["gamma_polynomial"]}
        assert rebuilt == expected
        norm = sum(abs(value) for powers,value in expected.items() if any(powers))
        assert norm == F(ff["nonconstant_coefficient_l1"])
        all_norms.append(norm)
    assert mixed_count == 735
    Bmax = max(all_norms)
    assert Bmax == F(frame["maximum_nonconstant_coefficient_l1"])
    assert rho == min(F(1,2),eps/(8*(1+Bmax)))
    assert bound == F(5,16)*eps
    assert len(features["features_with_repetitions"]) == 1570
    # Rebuild all 2,355 within-block product polynomials without importing the
    # builder, and compare all nonzero coefficient records.
    target_products = defaultdict(dict)
    monomials = [tuple(e) for e in products["monomials"]]
    assert len(monomials) == 237 and len(set(monomials)) == 237
    assert max(map(sum,monomials)) == 6
    for block,i,j,k,value in products["symmetric_upper_product_records"]:
        assert monomials[k] not in target_products[(block,i,j)]
        target_products[(block,i,j)][monomials[k]] = F(value)
    for block in frame["blocks"]:
        ell = {tuple(int(i==j) for j in range(5)):F(value)
               for i,value in enumerate(frame["frames"][block["frame"]]["linear_testing_function"]) if F(value)}
        V = poly_dict(block["testing_weight"])
        actual = [V,multiply(V,ell),multiply(multiply(V,ell),ell)]
        for (i,j),value in zip(((0,0),(0,1),(1,1)),actual):
            assert target_products[(block["id"],i,j)] == value
    assert len(frame["blocks"]) == 785
    assert len(products["symmetric_upper_product_records"]) == 27747
    # Cross-check the seven full target columns directly in SymPy.
    t = s.symbols("t0:5")
    test_z = [value+variable for value,variable in zip(za,t)]
    expected_targets = {}
    for k,record in enumerate(products["seven_source_target_polynomials"]):
        source = [s.Integer(int(i==k)) for i in range(7)]
        direct = sympy_records(direct_parent(source,test_z),t)
        assert direct == poly_dict(record["polynomial"])
        for powers,value in direct.items():
            expected_targets[(monomials.index(powers),k)] = value
    actual_targets = {(k,nu):F(value) for k,nu,value in products["seven_source_target_records"]}
    assert actual_targets == expected_targets and len(actual_targets) == 295
    assert products["off_diagonal_pairing_multiplier"] == 2
    assert products["gram_scalar_variable_count"] == 2355
    receipt = {
        "status":"PASS",
        "reviewer":"separate diagonal/extraction agent within the same model system",
        "method":"Direct SymPy derivatives and linear algebra; separate Fraction polynomial data cross-check; no C2 module imports; analytical proof reading.",
        "data_sha256":{name:hashlib.sha256((C2/name).read_bytes()).hexdigest() for name in (
            "regular_c2_frames.json","regular_c2_feature_list.json","regular_c2_product_map.json")},
        "analytical_findings":[
            "The six contact equations include the state value and all three normalized state derivatives; Euler homogeneity restores the fourth derivative.",
            "The contact source Jacobian is invertible throughout the stated box. The source coordinate bounds are retained, so no claim is made that every contact-box point lies in the bounded source graph.",
            "The exact polynomial Taylor identity and correct L inverse orientation retain all parent terms.",
            "All source/test coefficient separation is explicit; positivity holds on independent boxes and over independently labelled real closed fields.",
            "The matrix identity gives the stated local Gram representation. Ambient monomial count is not mislabelled as the rank of the product space.",
            "The seven source columns are full parent columns, but the normalized local contact proof does not prove a global projection onto the original eight quantum coordinates."
        ],
        "verified_counts":{"frames":25,"mixed_products":735,"pure_positive_remainders":50,"PSD2_blocks":785,"listed_features":1570,"ambient_monomials":237,"product_records":27747,"target_records":295},
        "jacobian_infinity_norm_bound":str(max(row_l1)*rho),
        "minimum_fresh_pure_coefficient_bound":str(min_pure_bound),
        "editorial_issue_reported":"Missing backslash before quad in matching equation; author reported it corrected.",
        "scope":"One normalized regular-C2 source/contact graph and independent testing box only.",
        "limitations":["Does not independently regenerate the full multivariate Taylor allocation, which the separate C2 checker tests.","No global source validity or finite all-source/all-test cover is supplied.","Not external peer review or proof-assistant formal verification."],
        "global_E1_E2_E3_complete":False,
    }
    (HERE/"REGULAR_C2_SEPARATE_REVIEW.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":"PASS","scope":"regular-C2 local data and proof","product_records":27747,"global_extraction_complete":False}))


if __name__ == "__main__":
    main()
