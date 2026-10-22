"""Division-free independent reader for the seven original Cramer bridges.
Substitute monomials, clear only displayed h powers, and compare cross-products.
Avoids expensive simplify/cancel of nested rational substitutions.
"""
from pathlib import Path
import json,os
import sympy as s
from sympy.polys.rings import ring
from sympy import QQ
R=Path(__file__).resolve().parents[1];OUT=Path(os.environ.get('Q222_REPORT_DIR',str(R/'verification_output')));OUT.mkdir(parents=True,exist_ok=True)
P,b,h,q,r,E=ring('b,hnew,q,r,E',QQ)
sv=s.symbols('b hnew q r E');ov=s.symbols('u v A Y e');loc={str(z):z for z in sv};oloc={str(z):z for z in ov}
raw=json.loads((R/'sources/ORIGINAL_C1_CONTACT.json').read_text());fix=json.loads((R/'sources/X_BOUNDED_V_SMALL_FULL_POLYNOMIALS.json').read_text())
def poly0(expr):return P.from_dict({m:QQ.convert(c)for m,c in s.Poly(expr,*sv).terms()})
def pulled_poly(expr):
 terms=[]
 for (i,j,k,l,m),c in s.Poly(expr,*ov).terms():terms.append(((i+k+l,i-k-l,i+j+2*m,l,m),QQ.convert(c)))
 hp=max([0]+[-mm[1]for mm,c in terms]);dic={}
 for mm,c in terms:
  nn=(mm[0],mm[1]+hp,*mm[2:]);dic[nn]=dic.get(nn,QQ.zero)+c
 return P.from_dict(dic),hp
Dc=poly0(s.sympify(fix['e3'],locals=loc))/2
checks=[]
for key in ['d','e2','e3','e4','e5','alpha','beta']:
 on,od=s.fraction(s.cancel(s.sympify(raw[key],locals=oloc)));N,nh=pulled_poly(on);D,dh=pulled_poly(od)
 cn,cd=s.fraction(s.sympify(fix[key],locals=loc));C,Q=poly0(cn),poly0(cd)
 assert C*D*h**nh == Q*Dc*N*h**dh,key
 checks.append('entire_original_Cramer_'+key)
Dr=Dc.exquo(r)
# Direct residue at b=q=0 is exactly 1, for every h,r,E.
res={m:c for m,c in Dr.items()if m[0]==0 and m[2]==0};assert res=={(0,0,0,0,0):QQ.one}
checks+=['D_divisible_by_r','D_over_r_residue_one']
report={'status':'PASS','count':len(checks),'checks':checks,'method':'exact monomial substitution and polynomial cross-products; no fraction-field division by a shape parameter'}
(OUT/'FAST_ORIGINAL_BRIDGE_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS',len(checks),'division-free original-coefficient bridge checks')
