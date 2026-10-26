"""Generate/check the inverse-Alice fixed polynomial weighted family.

Requires Python 3 and SymPy.  No solver, network, floating-point arithmetic,
or nonlocal filesystem access is used.  Run with --write once to emit the
JSON; subsequent invocations check the existing JSON independently of its
stored claims.  The certificate is LOCAL, not a lift of the full quantum set.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sympy as s

HERE = Path(__file__).resolve().parent
DATA = HERE / "inverse_alice_fixed_blocks.json"
REPORT = HERE / "inverse_alice_fixed_checks.json"
q, b, m, C, A, D, alpha, beta = s.symbols("q b m C A D alpha beta")
p, w, W, X, Y, Z = s.symbols("p w W X Y Z")
d = s.symbols("d")
z, Q = q**2, 1-b
delta = s.Rational(1, 2**12)
Dp, Dw = 1+p**2, 1+w**2
TEST = (p, w, W, X, Y, Z)
LOCAL = {str(v): v for v in (q,b,m,C,A,D,alpha,beta,p,w,W,X,Y,Z,d)}


def zero(expr):
    return s.Poly(s.expand(s.fraction(s.cancel(expr))[0]),
                  p,w,W,X,Y,Z,q,b,m,C,A,D,alpha,beta,d).is_zero


def matstr(matrix):
    return [[s.sstr(v) for v in row] for row in matrix.tolist()]


def outer(v, weight=1):
    v = s.Matrix(v)
    return weight*v*v.T


def generate():
    blocks = []
    def block(name, vector, multiplier, gram, role):
        blocks.append({"id": name, "testing_vector": [s.sstr(v) for v in vector],
                       "scalar_multiplier": s.sstr(multiplier),
                       "source_gram": matstr(s.Matrix(gram)), "role": role})
    phi = [X,Y,W,w*W,w*Z,p*X]
    L = s.Matrix([[1,0,0,1,0,0], [0,0,0,0,0,1],
                  [0,1/q,0,0,-b/q,0], [0,0,q,0,0,0],
                  [0,0,0,0,1/q,0]])
    channels = []
    for y, entries, r in [
        (p, [(3,2,-3*z/2),(3,4,-3*b*z/2)], 3*z),
        (p*(w/delta)**2, [(3,2,-2*delta**2),(3,4,-2*b*delta**2)],4*delta**2),
        (w/delta, [(1,2,-2*q*delta),(1,4,2*q*Q*delta)],4*q*delta)]:
        Aj = s.zeros(5)
        for i,j,c in entries:
            Aj[i,j] = Aj[j,i] = c
        channels.append((y,Aj,r))
    H0 = (s.Rational(1,9)-sum(ch[2] for ch in channels))*s.eye(5)
    block("main_0",phi,1,L.T*H0*L,"constant main reserve")
    for j,(yj,Aj,rj) in enumerate(channels,1):
        block(f"main_{j}_plus",phi,1+yj,L.T*(rj*s.eye(5)+Aj)*L/2,"bounded testing channel")
        block(f"main_{j}_minus",phi,1-yj,L.T*(rj*s.eye(5)-Aj)*L/2,"bounded testing channel")
    block("base_R",[X+w*W],1,[[s.Rational(17,9)]],"baseline square")
    block("base_S",[p*X],1,[[s.Rational(17,9)]],"baseline square")
    block("base_J",[Y,w*Z],1,outer([1/q,-b/q],s.Rational(8,9)),"baseline square")
    block("base_EG",[W,p*w*Z],1,outer([q,-s.Rational(36,25)/q],s.Rational(25,18)),"baseline square")
    block("base_G",[w*Z],1,[[s.Rational(2,225)/z]],"baseline square")
    block("base_gate_G",[w*Z],1-p**2,[[s.Rational(69,50)/z]],"baseline testing gate")
    tJ = s.Rational(3,2)/(b*Q)
    block("positive_J",[Y,w*Z],1,outer([1/q,-b/q],tJ-1),"selected source cap")
    block("positive_pJ",[Y,w*Z],p**2,outer([1/q,-b/q],tJ),"selected source cap")
    block("positive_pF",[p*Y],1,[[3*z/2]],"selected source cap")
    block("positive_wEF",[W-p*Y],w**2,[[3*z/2]],"selected source cap")
    block("positive_wF",[w*Y],Dp,[[s.Rational(3,2)/(b*z)]],"selected source cap")
    block("positive_wS",[p*w*X],1,[[2]],"selected source cap")
    block("positive_GF",[Y-w*Z],1,[[2]],"selected source cap")
    block("positive_pF_Dw",[p*Y],Dw,[[2]],"selected source cap")
    block("event_m",[p*W-p*w*X+Y-w*Z],1,[[4*m]],"whole event")
    block("event_D",[w*W+X-p*w*Y-p*Z],1,[[4*D]],"whole event")
    block("event_beta",[p*X+Z],Dw,[[4*beta]],"whole event")
    block("excess_C",[W-p*Y],Dw,[[4*(C-3*z/8)]],"whole event excess")
    block("excess_A",[Y],Dp*Dw,[[4*(A-s.Rational(3,8)/(b*z))]],"whole event excess")
    block("excess_alpha",[Y-w*Z],Dp,[[4*(alpha-s.Rational(3,8)/(Q*z))]],"whole event excess")
    features = []
    gates = [("one",s.Integer(1))]
    for j,(yj,_,_) in enumerate(channels,1):
        gates += [(f"plus_{j}",1+yj),(f"minus_{j}",1-yj)]
    for name, gate in gates:
        for i,v in enumerate(phi,1):
            features.append({"id":f"phi{i}_{name}","polynomial":s.sstr(v),"radicand":s.sstr(gate)})
    for i,v in enumerate([p*w*Z,p*Y,p*w*Y,w*Y,p*w*X,p*W,p*Z,Z],1):
        features.append({"id":f"extra_{i}","polynomial":s.sstr(v),"radicand":"1"})
    features.append({"id":"extra_9","polynomial":s.sstr(w*Z),"radicand":s.sstr(1-p**2)})
    target = {
        "d":2*(w*W+X-p*w*Y-p*Z)**2+2*Dp*(Y-w*Z)**2+2*Dw*(p*X+Z)**2-2*Dp*Dw*Z**2,
        "m":4*(p*W-p*w*X+Y-w*Z)**2,
        "C":4*Dw*(W-p*Y)**2,"A":4*Dp*Dw*Y**2,
        "D":4*(w*W+X-p*w*Y-p*Z)**2,
        "alpha":4*Dp*(Y-w*Z)**2,"beta":4*Dw*(p*X+Z)**2}
    return {"schema":"q222.polynomial_weighted_gram.v1",
      "scope":"inverse-Alice local testing chart; not a full-Q222 lift",
      "source_parameters":[str(v) for v in (q,b,m,C,A,D,alpha,beta)],
      "source_domain":["0<q<=1/8","0<b<1","m>=0","D>=0","beta>=0",
                       "C>=3*q**2/8","A>=3/(8*b*q**2)","alpha>=3/(8*(1-b)*q**2)"],
      "testing_variables":[str(v) for v in TEST],
      "testing_domain":["0<=p<=1","0<=w<=1/4096","W,X,Y,Z real"],
      "delta":"1/4096", "canonical_coefficient_order":["d","m","C","A","D","alpha","beta"],
      "target_coefficient_polynomials":{k:s.sstr(v) for k,v in target.items()},
      "normalized_target_d":"1",
      "physical_normalizer":{"equation":"kappa_t**2*(1+p**2)*(1+w**2)=1",
                              "feature_action":"multiply every fixed feature by kappa_t",
                              "unit_state":"W**2+X**2+Y**2+Z**2=1"},
      "blocks":blocks,"scalar_features":features,
      "claimed_counts":{"blocks":27,"total_order":66,"maximum_order":6,"labelled_scalar_features":51}}


def check(data):
    checks = {}
    def ck(name, predicate):
        if not predicate:
            raise AssertionError(name)
        checks[name] = True
    parse = lambda a:s.sympify(a,locals=LOCAL)
    blocks = data["blocks"]
    total = 0
    for item in blocks:
        v = s.Matrix([parse(t) for t in item["testing_vector"]])
        h = s.Matrix([[parse(a) for a in row] for row in item["source_gram"]])
        multiplier = parse(item["scalar_multiplier"])
        ck(item["id"]+"_source_separation",not h.free_symbols.intersection(TEST))
        ck(item["id"]+"_testing_separation",not (set().union(*(a.free_symbols for a in v))|multiplier.free_symbols).intersection((q,b,m,C,A,D,alpha,beta)))
        ck(item["id"]+"_symmetric",h == h.T and h.rows == v.rows)
        total += multiplier*(v.T*h*v)[0]
    target = sum((s.Integer(1) if k=="d" else LOCAL[k])*parse(v)
                 for k,v in data["target_coefficient_polynomials"].items())
    ck("complete_raw_parent_identity",zero(total-target))
    E,F,G,R,S = q*W,Y/q,w*Z/q,X+w*W,p*X
    J=F-b*G
    K=R**2+Dw*S**2+z*(G-F)**2+p**2*Dw*z*F**2-2*p*E*(G+w**2*F)+2*q*w*S*(G-F)
    ck("CE20_complete_core",zero(parse(data["target_coefficient_polynomials"]["d"])-2*K))
    B=2*R**2+2*S**2+J**2+s.Rational(3,2)*E**2+s.Rational(3,2)*(1+p**2)*G**2-4*p*E*G
    reserve=sum(t*t for t in (R,S,J,E,G))/9
    baseline=s.Rational(17,9)*(R**2+S**2)+s.Rational(8,9)*J**2+s.Rational(25,18)*(E-s.Rational(36,25)*p*G)**2+s.Rational(2,225)*G**2+s.Rational(69,50)*(1-p**2)*G**2
    ck("fixed_baseline_identity",zero(B-reserve-baseline))
    ck("reserve_margin",s.Rational(3,64)+4*delta**2+delta/2 < s.Rational(1,18))
    # Verify the full original event expansion, independently of block count.
    def weighted_sum(ids):
        out=0
        for v in blocks:
            if v["id"] in ids:
                f=s.Matrix([parse(a) for a in v["testing_vector"]]); h=s.Matrix([[parse(a) for a in row] for row in v["source_gram"]])
                out+=parse(v["scalar_multiplier"])*(f.T*h*f)[0]
        return out
    err=-(3*p*z+4*p*w**2)*E*(J+b*G)-4*q*w*S*(J-Q*G)
    main=weighted_sum([v["id"] for v in blocks if v["id"].startswith("main_")])
    ck("seven_block_gate_identity",zero(main-reserve-err))
    # All ungated polynomial roots of the twenty-two retained squares belong
    # to the fourteen-monomial base list; only source multipliers are omitted.
    ungated=[parse(t["polynomial"]) for t in data["scalar_features"] if t["radicand"]=="1"]
    monomials={mon for a in ungated for mon,c in s.Poly(a,*TEST).terms() if c}
    roots=[J,p*J,p*F,w*(E-p*z*F),w*F,p*w*F,w*S,G-F,p*F,p*w*F,
           p*W-p*w*X+Y-w*Z,w*W+X-p*w*Y-p*Z,p*X+Z,w*(p*X+Z),
           W-p*Y,w*(W-p*Y),Y,p*Y,w*Y,p*w*Y,Y-w*Z,p*(Y-w*Z)]
    ck("twenty_two_retained_roots",len(roots)==22)
    for i,root in enumerate(roots,1):
        ck(f"retained_{i}_fixed_span",all(mon in monomials for mon,c in s.Poly(s.expand(root),*TEST).terms() if c))
    ck("feature_count",len(data["scalar_features"])==51)
    orders=[len(v["testing_vector"]) for v in blocks]
    ck("block_counts",len(orders)==27 and sum(orders)==66 and max(orders)==6)
    # Native source elimination. Rational substitutions are exact.
    x,sigma,h,zz,bb=s.symbols("x sigma h zz bb")
    ee=x*h/sigma; rr=bb*zz; th=1-h; Ks=1+x; Rs=1+x*rr
    Ms=Rs-x*zz*th; Js=1+x*h-sigma*(1-rr)
    Ps=sigma+Ks+sigma*zz*(1+bb*x+2*x*h)+sigma*bb*x**2*h*zz**2+x**2*h*th*zz
    Ds=(Ks*(1+sigma*zz)*(1+sigma*rr)-Js*Ps)/(2*Js*Ps)
    H0s=1+zz*(1-bb); Bs=Ps+x*zz*(Js-sigma)
    calH=Rs*H0s+zz*((1-rr)*x*h*(2+x*rr)-ee*(2+x*(1+h+rr*th))+Ks*ee**2)
    ck("native_source_elimination",s.cancel(ee*Bs-x*H0s-sigma*calH+2*Ps*Js*Ds/sigma)==0)
    epi=s.symbols("epi",nonnegative=True)
    ep=s.Rational(1,2**30)
    ck("C_pi_cap_rational_margin",(1-6*ep)/(2*(1+5*ep))>s.Rational(3,8))
    ck("alpha_cap_rational_margin",(1-ep)/(2*(1+5*ep))>s.Rational(3,8))
    npi=zz*(1+x*h)*(1-ee)+(1-rr)*(1+x*h*zz)
    ck("opening_numerator_identity",s.cancel(npi-(H0s-zz*ee+x*h*zz*(2-rr-ee)))==0)
    # The symbolic inequalities Ps<=1+5eps and calH<=4 are established by
    # positive term bounds, checked here at their maximal rational endpoints.
    Pmajor=1+2*ep+ep*s.Rational(1,64)*(1+3*ep)+ep**3/s.Integer(4096)+ep**2/s.Integer(64)
    Hmajor=(1+ep/s.Integer(64))*s.Rational(65,64)+s.Rational(1,64)*(ep*(2+ep/s.Integer(64))+1+ep)
    ck("P_termwise_majorant",Pmajor<=1+5*ep)
    ck("calH_termwise_majorant",Hmajor<=4)
    report={"scope":data["scope"],"exact_checks":len(checks),"checks":checks,
            "counts":{"blocks":len(orders),"total_order":sum(orders),"maximum_order":max(orders),"scalar_features":len(data["scalar_features"]),"retained_roots":22},
            "reserve_surplus_over_one_eighteenth":str(s.Rational(1,18)-(s.Rational(3,64)+4*delta**2+delta/2)),
            "verification_limit":"Algebraic identity/support checks; ordered-field source inequalities are proved in the accompanying TeX. No global coverage assertion."}
    REPORT.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:report[k] for k in ("scope","exact_checks","counts","reserve_surplus_over_one_eighteenth")},indent=2))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true",help="generate the adjacent explicit JSON, then check it")
    args=parser.parse_args()
    if args.write:
        DATA.write_text(json.dumps(generate(),indent=2)+"\n",encoding="utf-8")
    check(json.loads(DATA.read_text(encoding="utf-8")))
