#!/usr/bin/env python3
"""Independent same-model review of the sibling Schur-remainder component.

This script does not import the sibling checker or its PASS receipt. It builds
the state matrix directly from rank-one packet matrices, then derives the
determinants and resultants again. Requires SymPy 1.14 (review only).
"""
import hashlib
import json
from pathlib import Path
import sympy as s

HERE = Path(__file__).resolve().parent
C2 = HERE.parent / "c2"
u, v, t, U, V = s.symbols("u v t U V")
Q = s.Rational
Du, Dv = 1+u*u, 1+v*v
tests = []


def check(name, truth):
    if not bool(truth):
        raise AssertionError(name)
    tests.append({"name": name, "status": "PASS"})


def ident(name, polynomial):
    check(name, s.cancel(polynomial) == 0)


def gram(row):
    z = s.Matrix(row)
    return z*z.T


# State order W,X,Y,Z. This construction uses packet rows and weights directly.
M0 = (Q(7,8)*Du*Dv*gram([0,1,0,0])
      +Q(9,8)*gram([v,-1,u*v,-u])
      +Q(1,18)*Dv*gram([u,0,-1,0])
      +Q(11,8)*Du*Dv*gram([0,0,1,0])
      +Q(13,8)*gram([u,u*v,-1,-v])
      +Q(5,9)*Du*gram([0,0,v,-1])
      +Q(1,2)*Dv*gram([0,1,0,u])
      -Q(1,2)*Du*Dv*gram([0,0,0,1]))
Mt = M0+t*Du*Dv*gram([0,1,0,0])
core_rows = s.Matrix([[u,u*v,-1],[0,0,v],[0,1,0]])
ident("strict_upper_core_row_determinant", core_rows.det()+u*v)
masses = [880,9317,720,880,9477,720]
exponents = [(1,0),(1,-1),(0,1),(0,-1),(-1,1),(-1,0)]
check("amgm_total_mass", sum(masses) == 21994)
check("amgm_exponent_balance", all(sum(w*e[k] for w,e in zip(masses,exponents)) == 0 for k in (0,1)))
check("amgm_masses_strictly_positive", all(w>0 for w in masses))
laurent = Q(121,331776)*(sum(w*U**e[0]*V**e[1] for w,e in zip(masses,exponents))-21994)
ident("independent_full_determinant", M0.det(method="domain-ge")
      -Du**2*Dv**2*u*u*v*v*laurent.subs({U:u*u,V:v*v}))
Mline = Mt.subs(v,1)
A = 26883*U*U+105314*U+29403
B = 29000*U*U+50384*U+21384
D = A+t*B
N = 112167*(U-1)**2+t*(121000*U*U-189648*U+81576)
det_upper = Mline[:3,:3].det(method="domain-ge")
det_full = Mline.det(method="domain-ge")
ident("independent_restricted_upper_determinant", det_upper-Q(11,20736)*(1+u*u)*D.subs(U,u*u))
ident("independent_restricted_full_determinant", det_full-Q(11,82944)*(1+u*u)**2*N.subs(U,u*u))
ident("independent_literal_remainder", det_full/det_upper-((U+1)*N/(4*D)).subs(U,u*u))
check("independent_resultant_A_B", s.resultant(A,B,U) == 1386911105418617856)
ident("independent_resultant_D_N", s.resultant(D,N,U)
      -10441728*(8*t+11)*(1000*t+927)*(33677*t+55550)**2)
ident("independent_D_minus_one", D.subs(U,-1)+49028)
ident("independent_D_discriminant", s.discriminant(D,U)
      -448*(129472*t*t+10942118*t+17699275))
check("positive_denominator_coefficients", all(all(c>0 for c in s.Poly(p,t).all_coeffs())
                                              for p in s.Poly(D,U).all_coeffs()))
contact = s.Matrix([1,Q(-1,10),Q(1,10),1])
check("seed_contact_kernel", M0.subs({u:1,v:1})*contact == s.zeros(4,1))
check("matrix_is_biquadratic_in_angles", all(s.degree(entry,u)<=2 and s.degree(entry,v)<=2 for entry in Mt))
for name, G, substitutions in [
    ("u",s.diag(1,1,-1,-1),{u:-u}),
    ("v",s.diag(1,-1,1,-1),{v:-v})]:
    check("independent_sign_gauge_"+name, all(s.expand(entry)==0 for entry in G*Mt.subs(substitutions)*G-Mt))

review = {
    "status": "PASS_WITH_STATED_SCOPE", "review_kind": "independent same-model analytical review and separately constructed exact packet-matrix check",
    "runtime": "SymPy "+s.__version__, "fresh_checks": tests, "fresh_check_count": len(tests),
    "analytical_observations": [
        {"topic": "positivity of complete seed", "finding": "The signed core restricted to Z=0 is strictly positive when uv!=0 because its three packet rows have determinant -uv. Positive event packets preserve this. The nonnegative full determinant from balanced weighted AM-GM then gives a nonnegative scalar Schur complement."},
        {"topic": "axes and projective endpoints", "finding": "The matrix is polynomial and separately of degree at most two in each angle. Continuity after bihomogenization includes axes and projective endpoints. This reasoning concerns the complete parent, with source t>=0."},
        {"topic": "literal remainder linear independence", "finding": "D_t has degree two for every t>=0. The displayed nonzero resultants imply coprimality of its numerator and each other denominator. Clearing a finite relation and reducing modulo each D_t forces its scalar coefficient to vanish."},
        {"topic": "fixed-feature obstruction scope", "finding": "An arbitrary fixed finite function space has a finite-dimensional product span. Thus it cannot contain SOS certificates for every literal S_t. The conclusion is about these literal scalar remainders, not about all possible certificates of the parent."},
        {"topic": "complete-family five features", "finding": "M_t=M_0+t*Du*Dv*e_X e_X^T, and M_0 is PSD. Four components of its unique PSD square root applied to the state, plus sqrt(Du*Dv)*X, give five source-independent testing functions. The square root is semialgebraic through its defining polynomial and PSD conditions."},
        {"topic": "global scope", "finding": "Neither the scalar obstruction nor the five-feature one-parameter family proves or disproves full-Q222 finite-lift existence or supplies a global feature atlas. The proof explicitly states this limitation."}
    ],
    "requested_changes": [],
    "not_claimed": ["external peer review", "proof-assistant verification", "global Q222 extraction", "fresh audit of historical R21"],
    "reviewed_files": {name: hashlib.sha256((C2/name).read_bytes()).hexdigest()
                       for name in ["schur_remainder_proof.tex","schur_remainder_check.py","schur_remainder_data.json"]}
}
(HERE/"INDEPENDENT_SCHUR_REVIEW.json").write_text(json.dumps(review,indent=2,sort_keys=True)+"\n")
print(json.dumps({"status":review["status"],"fresh_checks":len(tests),"requested_changes":0},sort_keys=True))
