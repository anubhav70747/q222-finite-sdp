import sympy as S,pickle,json
from fractions import Fraction as F
from pathlib import Path
O=Path(__file__).resolve().parents[2]/'data'
d=pickle.load((O/'lII_normalized_full.pkl').open('rb'))
Q=d['Q'];variables=d['variables'];bounds=d['bounds'];base=d['base']
u,v,hh,r,tt,zz,kk,db,ds,dc,dg,dal,dA,dm,dn,dq=variables
zero={x:0 for x in [zz,kk,db,ds,dc,dg,dal,dA,dm,dn,dq]}
hmin=F(1,64)
def groups(expr):
 num,den=S.fraction(S.cancel(expr));fac=S.factor_list(den)
 hp=0;lower=F(fac[0]);out={}
 for p,e in fac[1]:
  if p==hh:hp+=e
  elif p==hh**2+1:pass
  elif p in [db+1,ds+1]:lower*=F(1,2)**e
  else:raise ValueError((p,e))
 for mon,c in S.Poly(num,*variables).terms():
  term=abs(F(c))/lower
  for var,pow in zip(variables,mon):
   if var==hh:
    pow-=hp
    term*=hmin**pow if pow<0 else F(1)
   else:term*=bounds[var]**pow
  key=mon[:2];out[key]=out.get(key,F(0))+term
 return out
tot={}
for i in range(6):
 for j in range(i,6):
  lead=Q[i,j].subs(zero)
  for ex in [Q[i,j]-lead,lead-base[i,j]]:
   for key,b in groups(ex).items():tot[key]=tot.get(key,F(0))+(1 if i==j else 2)*b
def dyad(q):
 if q==0:return F(0)
 x=F(1)
 while x<q:x*=2
 while x/2>=q:x/=2
 return x
rows=[]
for key,b in sorted(tot.items()):
 alpha=dyad(b)
 rows.append({'exponents':list(key),'budget':str(alpha),'raw_bound':str(b)})
 print(key,'budget=',alpha,'raw=',float(b))
assert sum(dyad(b) for b in tot.values())<F(1,8)
receipt={'domain':{'H_min':'1/4096','u_max':'1/4096','v_max':'1/33554432','epsilon':'1/18446744073709551616'},'budgets':rows,'total_budget':str(sum(dyad(b) for b in tot.values())),'total_lt_one_eighth':True,'p_monomial':'(4096*u)^i*(33554432*v)^j','exclude_constant_from_gate_list':True}
(O/'LII_FIXED_GATE_BUDGETS.json').write_text(json.dumps(receipt,indent=2))
