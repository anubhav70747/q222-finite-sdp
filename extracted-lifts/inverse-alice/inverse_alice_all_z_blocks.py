"""Exact polynomial weighted data for the all-z inverse-Alice chart.

Python 3 + SymPy only.  --write generates the adjacent JSON for P=2**24
(or --alice-bound P); otherwise the stored data is checked.  No optimizer,
floating-point arithmetic, network, or global coverage assertion is used.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sympy as s

HERE=Path(__file__).resolve().parent
DATA=HERE/"inverse_alice_all_z_blocks.json"
REPORT=HERE/"inverse_alice_all_z_checks.json"
q,b,tB,tH,m,C,A,D,alpha,beta=s.symbols("q b tB tH m C A D alpha beta")
p,w,W,X,Y,Z,d=s.symbols("p w W X Y Z d")
TEST=(p,w,W,X,Y,Z)
SOURCE=(q,b,tB,tH,m,C,A,D,alpha,beta)
LOCAL={str(v):v for v in (*TEST,*SOURCE,d)}
z,Q=q*q,1-b
H0=1+z*Q
B=H0/(b*Q)
r=b*z
Dp,Dw=1+p*p,1+w*w


def out(v,weight=1):
    v=s.Matrix(v)
    return weight*v*v.T


def reduced_numerator(expr):
    numerator=s.Poly(s.fraction(s.cancel(expr))[0],tB,tH)
    result=0
    for (i,j),coef in numerator.terms():
        result+=coef*B**(i//2)*H0**(j//2)*tB**(i%2)*tH**(j%2)
    return s.cancel(result)


def generate(P):
    P=s.Integer(P)
    eta=s.Rational(1,2**30)/P**2
    delta=s.Rational(1,2**12)/P
    f=1-25*eta
    tt=p/P
    C0=z/(2*(1-r)); A0=H0/(2*b*z); alpha0=1/(2*Q*z)
    aa,cc=z/tB,1/tH
    dn=Q*z*(1-r)/H0
    phi=[X,Y,W,w*W,w*Z,p*X]
    # Ordered normal roots (R,S,J*,Ebar,G*) = L phi.
    L=s.Matrix([[1,0,0,1,0,0],[0,0,0,0,0,1],
                [0,tB/q,0,0,-b*tB/q,0],
                [0,0,q/P,0,0,0],[0,0,0,0,tH/q,0]])
    channels=[]
    for y,entries,rj in [
        (tt,[(3,4,-2*(1-f)*P**2/tH)],50*eta*P**2),
        (tt*(w/delta)**2,[(3,2,-2*P**2*delta**2/tB),(3,4,-2*b*P**2*delta**2/tH)],4*P**2*delta**2),
        (w/delta,[(1,2,-2*q*delta/tB),(1,4,2*q*Q*delta/tH)],4*delta)]:
        M=s.zeros(5)
        for i,j,v in entries:
            M[i,j]=M[j,i]=v
        channels.append((y,M,rj))
    blocks=[]
    def block(name,vector,multiplier,gram,proof):
        blocks.append({"id":name,"testing_vector":[s.sstr(v) for v in vector],
                       "scalar_multiplier":s.sstr(multiplier),
                       "source_gram":[[s.sstr(a) for a in row] for row in s.Matrix(gram).tolist()],
                       "positivity_proof":proof})
    reserve=(f/2-sum(ch[2] for ch in channels))*s.eye(5)
    block("main_0",phi,1,L.T*reserve*L,"source congruence of a positive scalar matrix")
    for j,(yj,M,rj) in enumerate(channels,1):
        for name,sgn in [("plus",1),("minus",-1)]:
            block(f"main_{j}_{name}",phi,1+sgn*yj,L.T*(rj*s.eye(5)+sgn*M)*L/2,
                  "source congruence of an absolute-row-sum bound")
    lJ=s.Matrix([tB/q,-b*tB/q]); lG=s.Matrix([0,tH/q])
    h1=s.Matrix([q/P,-z/(P*q),-(1-r)/(P*q)])
    h43=s.Matrix([q/P,-4*z/(3*P*q),-4*(1-r)/(3*P*q)])
    block("positive_RS",[X+w*W,p*X],1,(2-f/2)*s.eye(2),"positive scalar identity matrix")
    block("positive_EL",[W,p*Y,p*w*Z],1,
          out(h43,3*f/2)+out(h1,2*f*(P**2/(1-r)-1)),"sum of two source rank-one positive matrices")
    block("positive_JG",[Y,w*Z],1,
          out(lJ,3*f/2)+out(lG,5*f/6)+out([1,-1],2*(1-f)),"sum of three source rank-one positive matrices")
    block("positive_gate_G",[w*Z],1-tt**2,[[2*f*tH**2/(3*q**2)]],"positive source scalar")
    block("positive_pJG",[Y,w*Z],p**2,
          out(lJ,f*(2-s.Rational(2,3)/P**2))
          +(out(lJ)+out(lG))*2*f*dn/(3*P**2)
          +out(cc*lJ-aa*lG,2*f/(3*P**2))
          +out([1,0],2*(1-f)),"sum of source rank-one positive matrices")
    block("positive_pwXY",[p*w*X,p*w*Y],1,2*s.eye(2),"positive scalar identity matrix")
    block("full_C",[W-p*Y,w*(W-p*Y)],1,
          s.diag(4*(C-f*C0),4*C),"diagonal nonnegative source coefficients")
    block("full_A",[Y,w*Y],Dp,s.diag(4*(A-f*A0),4*A),"diagonal nonnegative source coefficients")
    block("event_m",[p*W-p*w*X+Y-w*Z],1,[[4*m]],"nonnegative source event coefficient")
    block("event_D",[w*W+X-p*w*Y-p*Z],1,[[4*D]],"nonnegative source event coefficient")
    block("event_beta",[p*X+Z],Dw,[[4*beta]],"nonnegative source event coefficient")
    block("event_alpha",[Y-w*Z],Dp,[[4*(alpha-f*alpha0)]],"nonnegative source event excess")
    target={
      "d":2*(w*W+X-p*w*Y-p*Z)**2+2*Dp*(Y-w*Z)**2+2*Dw*(p*X+Z)**2-2*Dp*Dw*Z**2,
      "m":4*(p*W-p*w*X+Y-w*Z)**2,
      "C":4*Dw*(W-p*Y)**2,"A":4*Dp*Dw*Y**2,
      "D":4*(w*W+X-p*w*Y-p*Z)**2,
      "alpha":4*Dp*(Y-w*Z)**2,"beta":4*Dw*(p*X+Z)**2}
    scalar=[]
    gates=[("one",s.Integer(1))]
    for j,(yj,_,_) in enumerate(channels,1):
        gates += [(f"plus_{j}",1+yj),(f"minus_{j}",1-yj)]
    for name,g in gates:
        for i,v in enumerate(phi,1):
            scalar.append({"id":f"phi{i}_{name}","polynomial":s.sstr(v),"radicand":s.sstr(g)})
    for i,v in enumerate([p*w*Z,p*Y,p*w*Y,w*Y,p*w*X,p*W,p*Z,Z],1):
        scalar.append({"id":f"extra_{i}","polynomial":s.sstr(v),"radicand":"1"})
    scalar.append({"id":"extra_9","polynomial":s.sstr(w*Z),"radicand":s.sstr(1-tt**2)})
    return {"schema":"q222.polynomial_weighted_gram.v1",
      "scope":"all-z inverse-Alice local source family, testing p in [0,P]; not full-Q222",
      "constants":{"P":str(P),"eta":str(eta),"delta":str(delta),"f":str(f)},
      "source_parameters":[str(v) for v in SOURCE],
      "source_domain":["q>0","0<b<1","b*q**2<1","tB>0","tH>0","m>=0","D>=0","beta>=0",
                        f"C>=({f})*q**2/(2*(1-b*q**2))",f"A>=({f})*(1+(1-b)*q**2)/(2*b*q**2)",f"alpha>=({f})/(2*(1-b)*q**2)"],
      "source_algebraic_equations":["tB**2=(1+(1-b)*q**2)/(b*(1-b))","tH**2=1+(1-b)*q**2"],
      "native_admission":{"condition":f"0<x,sigma<=({eta})/(1+z); native strict gaps; D,beta>=0",
                          "proof":"inverse_alice_all_z_completion.tex"},
      "testing_variables":[str(v) for v in TEST],
      "testing_domain":[f"0<=p<={P}",f"0<=w<={delta}","W,X,Y,Z real"],
      "canonical_coefficient_order":["d","m","C","A","D","alpha","beta"],
      "target_coefficient_polynomials":{k:s.sstr(v) for k,v in target.items()},"normalized_target_d":"1",
      "physical_normalizer":{"equation":"kappa_t**2*(1+p**2)*(1+w**2)=1",
                             "feature_action":"multiply every fixed feature by kappa_t",
                             "unit_state":"W**2+X**2+Y**2+Z**2=1"},
      "blocks":blocks,"scalar_features":scalar,
      "claimed_counts":{"blocks":19,"total_order":62,"maximum_order":6,"labelled_scalar_features":51}}


def check(data):
    parse=lambda a:s.sympify(a,locals=LOCAL)
    checks={}
    def ck(name,condition):
        if not condition:
            raise AssertionError(name)
        checks[name]=True
    total=0
    for item in data["blocks"]:
        vv=s.Matrix([parse(a) for a in item["testing_vector"]])
        HH=s.Matrix([[parse(a) for a in row] for row in item["source_gram"]])
        gg=parse(item["scalar_multiplier"])
        ck(item["id"]+"_pure_source",not HH.free_symbols.intersection(TEST))
        ck(item["id"]+"_pure_testing",not(gg.free_symbols|vv.free_symbols).intersection(SOURCE))
        ck(item["id"]+"_symmetric",all(s.cancel(a)==0 for a in HH-HH.T) and HH.rows==vv.rows)
        total+=gg*(vv.T*HH*vv)[0]
    target=sum((1 if k=="d" else LOCAL[k])*parse(v) for k,v in data["target_coefficient_polynomials"].items())
    ck("complete_raw_parent_identity",reduced_numerator(total-target)==0)
    P=s.Integer(data["constants"]["P"]); eta=s.Rational(data["constants"]["eta"])
    delta=s.Rational(data["constants"]["delta"]); f=s.Rational(data["constants"]["f"])
    tt=p/P; aa=z/tB; cc=1/tH; dn=Q*z*(1-r)/H0
    ck("angular_norm_defect",reduced_numerator(1-aa*aa-cc*cc-dn)==0)
    E,F,G,R,S=q*W,Y/q,w*Z/q,X+w*W,p*X
    J=F-b*G; Jstar=tB*J; Gstar=tH*G; Ebar=E/P; LL=aa*Jstar+cc*Gstar
    Q0=2*R*R+2*S*S+2*P**2/(1-r)*(Ebar-tt*LL)**2+2*(1+p*p)*Jstar**2+2*Gstar**2
    K0=R*R+S*S+z*(G-F)**2+p*p*z*F*F-2*p*E*G
    leading=2*K0+2*z/(1-r)*(W-p*Y)**2+2*H0/(b*z)*Dp*Y**2+2/(Q*z)*Dp*(Y-w*Z)**2
    ck("exact_leading_completion",reduced_numerator(Q0-leading)==0)
    clipped=2*R*R+2*S*S+2*(Ebar-tt*LL)**2+2*(1+tt**2)*Jstar**2+2*Gstar**2
    norm=R*R+S*S+Ebar**2+Jstar**2+Gstar**2
    sos=s.Rational(3,2)*(R*R+S*S)+s.Rational(3,2)*(Ebar-s.Rational(4,3)*tt*LL)**2+s.Rational(3,2)*Jstar**2+s.Rational(5,6)*Gstar**2+s.Rational(2,3)*(1-tt**2)*Gstar**2+s.Rational(4,3)*tt**2*Jstar**2+s.Rational(2,3)*tt**2*(dn*(Jstar**2+Gstar**2)+(cc*Jstar-aa*Gstar)**2)
    ck("clipped_reserve_SOS",reduced_numerator(clipped-norm/2-sos)==0)
    extra=2*(P**2/(1-r)-1)*(Ebar-tt*LL)**2+2*(P**2-1)*tt**2*Jstar**2
    ck("leading_positive_excess",reduced_numerator(Q0-clipped-extra)==0)
    rsum=50*eta*P**2+4*P**2*delta**2+4*delta
    ck("evaluated_error_bound",rsum<s.Rational(1,512))
    ck("positive_main_reserve",f/2-rsum>s.Rational(1,4))
    # Every ungated term is in the exact fourteen-monomial fixed space.
    base=[parse(t["polynomial"]) for t in data["scalar_features"] if t["radicand"]=="1"]
    mon={e for a in base for e,c in s.Poly(a,*TEST).terms() if c}
    for item in data["blocks"]:
        for i,entry in enumerate(item["testing_vector"]):
            ck(item["id"]+f"_span_{i}",all(e in mon for e,c in s.Poly(parse(entry),*TEST).terms() if c))
    orders=[len(v["testing_vector"]) for v in data["blocks"]]
    ck("literal_block_counts",len(orders)==19 and sum(orders)==62 and max(orders)==6)
    ck("literal_scalar_count",len(data["scalar_features"])==51)
    ck("fixed_monomial_count",len(mon)==14)
    report={"scope":data["scope"],"exact_checks":len(checks),"checks":checks,
            "constants":data["constants"],"error_reserve":str(rsum),"main_diagonal_reserve":str(f/2-rsum),
            "counts":{"blocks":19,"total_order":62,"maximum_order":6,"scalar_features":51,"ungated_monomials":14},
            "verification_limit":"Exact identities/support/counts only; native source inequalities and ordered-field positivity are proved separately. No global claim."}
    REPORT.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:report[k] for k in ("scope","exact_checks","constants","error_reserve","main_diagonal_reserve","counts")},indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    parser.add_argument("--alice-bound",type=int,default=2**24)
    args=parser.parse_args()
    if args.alice_bound<1:
        raise ValueError("P must be an integer >=1")
    if args.write:
        DATA.write_text(json.dumps(generate(args.alice_bound),indent=2)+"\n",encoding="utf-8")
    check(json.loads(DATA.read_text(encoding="utf-8")))
