"""Targeted exact full-state checks for the R6.5--R6.7 endpoint repair.

Run with: py -3.14 check_homogeneous.py
No parameter sweep, source-existence claim, or global proof certification.
"""
import json
from pathlib import Path

import sympy as S

u, v, p, q, X, W, Y, Z, H, T, Pz, Qz = S.symbols(
    "u v p q X W Y Z H T Pz Qz", real=True
)
n, m, C, A, D, alpha, K, g, mu, kappa, ac, cc, dc, mc = S.symbols(
    "n m C A D alpha K g mu kappa ac cc dc mc", real=True
)
half = S.Rational(1, 2)
Du, Dv, Dp, Dq = 1 + u**2, 1 + v**2, 1 + p**2, 1 + q**2
M, J = v*W-X+u*v*Y-u*Z, u*W+u*v*X-Y-v*Z
core = J**2 + Du*(v*Y-Z)**2 + Dv*(X+u*Z)**2 - Du*Dv*Z**2
parent = n*Du*Dv*X**2 + m*M**2 + C*Dv*(u*W-Y)**2
parent += A*Du*Dv*Y**2 + D*J**2 + alpha*Du*(v*Y-Z)**2 + half*core

Mu, Ju = p*v*W-p*X+v*Y-Z, W+v*X-p*Y-p*v*Z
cu = Ju**2 + Dp*(v*Y-Z)**2 + Dv*(p*X+Z)**2 - Dp*Dv*Z**2
Pu = n*Dp*Dv*X**2 + m*Mu**2 + C*Dv*(W-p*Y)**2
Pu += A*Dp*Dv*Y**2 + D*Ju**2 + alpha*Dp*(v*Y-Z)**2 + half*cu

Mv, Jv = W-q*X+u*Y-q*u*Z, q*u*W+u*X-q*Y-Z
cv = Jv**2 + Du*(Y-q*Z)**2 + Dq*(X+u*Z)**2 - Du*Dq*Z**2
Pv = n*Du*Dq*X**2 + m*Mv**2 + C*Dq*(u*W-Y)**2
Pv += A*Du*Dq*Y**2 + D*Jv**2 + alpha*Du*(Y-q*Z)**2 + half*cv

Mb, Jb = p*W-p*q*X+Y-Qz, q*W+X-p*q*Y-Pz
cb = Jb**2 + Dp*(Y-Qz)**2 + Dq*(p**2*X**2+2*X*Pz-Pz**2)
Pb = n*Dp*Dq*X**2 + m*Mb**2 + C*Dq*(W-p*Y)**2
Pb += A*Dp*Dq*Y**2 + D*Jb**2 + alpha*Dp*(Y-Qz)**2 + half*cb

checks = []
def zero(name, expr):
    expanded = S.cancel(S.expand(expr))
    if expanded != 0:
        raise AssertionError(f"{name}: {expanded}")
    checks.append(name)

def gram(expr, roots):
    return S.hessian(expr, roots)/2

zero("single-u full-state polynomial homogenization", p**2*parent.subs(u, 1/p)-Pu)
zero("single-v full-state polynomial homogenization", q**2*parent.subs(v, 1/q)-Pv)
zero("double full five-root homogenization", p**2*q**2*parent.subs({u:1/p,v:1/q})-Pb.subs({Pz:p*Z,Qz:q*Z}))

for name, row, arow in (("u", Mu, p*v*W), ("v", Mv, W), ("uv", Mb, p*W)):
    zero(f"{name} exact centered m extraction", m*row**2-(m-K)*(row-g*arow)**2
         -(K*row**2+2*(m-K)*g*arow*row-(m-K)*g**2*arow**2))

Ru = Pu.subs(m,K)+2*mu*p*v*W*Mu-mu*g*p**2*v**2*W**2
Rv = Pv.subs(m,K)+2*mu*W*Mv-mu*g*W**2
Rb = Pb.subs(m,K)+2*mu*p*W*Mb-mu*g*p**2*W**2
resu = Ru.subs({p:0,g:0,alpha:0,mu:kappa})
resv = Rv.subs({q:0,g:0,alpha:0,mu:kappa})
resb = Rb.subs({p:0,q:0,g:0,alpha:0,mu:kappa})
r65 = n*Dv*X**2+C*Dv*W**2+A*Dv*Y**2+(D+half)*(W+v*X)**2
r66 = (n+half)*Du*X**2+C*Du**2*Y**2+(A+half)*Du*Y**2+D*(u*X-Z)**2
r67 = (n+half)*X**2+C*W**2+A*Y**2+D*(X-Pz)**2+(K+half)*(Y-Qz)**2
zero("R6.5 complete kernel restriction", resu.subs(Z,v*Y)-r65)
zero("R6.6 complete kernel restriction including centered correction", resv.subs(W,-u*Y)-r66)
zero("R6.7 complete five-root residue", resb-r67)
zero("R6.5 full residue restoration", resu-r65-(K+half)*(v*Y-Z)**2)
zero("R6.6 full residue restoration", resv-((n+half)*Du*X**2+C*(u*W-Y)**2+(A+half)*Du*Y**2+D*(u*X-Z)**2+K*(W+u*Y)**2+2*kappa*W*(W+u*Y)))
zero("R6.5 determinant at n0=0", gram(r65.subs(n,0),(X,W,Y)).det()-A*C*(D+half)*v**2*Dv**2)
zero("R6.6 determinant", gram(r66,(X,Y,Z)).det()-D*(n+half)*Du**2*(C*Du+A+half))
zero("R6.7 determinant", gram(r67,(X,W,Y,Pz,Qz)).det()-A*C*D*(n+half)*(K+half))

# Exact exceptional corner allocation after retaining n,m,D,alpha and
# nonnegative C-cc,A-ac excess event squares. No p or v division is used.
Jhat = W+Dp*T-p*Y-p*v*H
chat = Jhat**2 + Dp*v**2*Y**2-2*Dp*v*Y*H+(1-p**2*v**2)*H**2
chat += 2*p*Dp*T*Y+2*p*v*Dp*T*H-p**2*Dp*T**2
exceptional = cc*Dv*(W-p*Y)**2+ac*Dp*Dv*Y**2+half*chat
selected = cc*Dv*(W-p*Y)**2+ac*Dp*Dv*Y**2+half*cu
zero("infinity-zero exact whole-core four-root allocation", selected.subs(Z,H-p*X)-exceptional.subs(T,v*X))
exceptional_residue = exceptional.subs({p:0,v:0})
zero("infinity-zero four-root residue", exceptional_residue-(cc*W**2+ac*Y**2+half*(W+T)**2+half*H**2))
zero("infinity-zero four-root determinant", gram(exceptional_residue,(H,W,Y,T)).det()-ac*cc/4)

# R7: retain n,alpha and actual event excesses in m,A,D, not the
# contact-centered rows, then examine the complete capped polynomial.
subcap = {n:0,alpha:0,m:mc,C:cc,A:ac,D:dc}
r7u = cc*Dv*W**2+ac*Dv*Y**2+(dc+half)*(W+v*X)**2+(mc+half)*(v*Y-Z)**2
r7v = mc*(W+u*Y)**2+cc*(u*W-Y)**2+(ac+half)*Du*Y**2+dc*(u*X-Z)**2+half*Du*X**2
r7b = cc*W**2+ac*Y**2+dc*(X-Pz)**2+half*X**2+(mc+half)*(Y-Qz)**2
zero("R7 RU28 full four-state capped residue", Pu.subs(subcap).subs(p,0)-r7u)
zero("R7 RU29 full four-state capped residue", Pv.subs(subcap).subs(q,0)-r7v)
zero("R7 RU30 full five-root capped residue", Pb.subs(subcap).subs({p:0,q:0})-r7b)
zero("R7 RU28 determinant", gram(r7u,(X,W,Y,Z)).det()-cc*ac*Dv**2*(dc+half)*(mc+half)*v**2)

output = {
    "scope": "Targeted polynomial identities, support and residual determinants only",
    "sympy_version": S.__version__,
    "checks": checks,
    "checked_count": len(checks),
    "limitations": [
        "No source nonemptiness or universal source classification checked",
        "No local contact conductor or finite remote proof checked here",
        "Tensor SOS uses the existing bounded separated-Gram lemma, not determinant positivity alone",
    ],
}
Path(__file__).with_name("exact-checks.json").write_text(json.dumps(output,indent=2)+"\n",encoding="utf-8")
print(f"Verified {len(checks)} targeted exact identities and determinants.")
