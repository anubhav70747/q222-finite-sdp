"""Independent verification of original-generator identities and full Pauli parent.
No Groebner basis, producer labels, saved reduced remainder, or floating point is
used as an acceptance oracle. All products are rational sparse polynomials.
"""
from pathlib import Path
import json,time
import sympy as s
from sympy import QQ
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1]; start=time.time(); count=0
import os
OUT=Path(os.environ.get('Q222_REPORT_DIR', str(ROOT/'verification_output')));OUT.mkdir(parents=True,exist_ok=True)
R,u,W,v,q,b,a,r,E,X,Y,e,rr,ss=ring('u,W,v,q,b,a,r,E,X,Y,e,rr,ss',QQ)
n=13
def ck(x,message):
 global count
 if not x: raise AssertionError(message)
 count+=1
def readterms(ts):
 return R.from_dict({tuple(mm)+(0,)*(n-len(mm)):QQ.convert(s.Rational(cc)) for mm,cc in ts})
g=[u*W-b*b*v,b*u*v,b*b*(v*v-q*q)]; small=[u,W,v,q,b]
raw=json.loads((ROOT/'sources/XSMALL_Z_FULL_SOURCE.json').read_text())
def original_part(key):
 d={}
 for mm,cc in raw[key]['sparse']:
  aa,bb,qq,uu,vv,ww,rrr,ee=mm
  d[(uu,ww,vv,qq,bb,aa,rrr,ee,0,0,0,0,0)]=QQ.convert(s.Rational(cc))
 return R.from_dict(d)
reports=[]
for kind in ['scalar','x','y']:
 d=json.loads((ROOT/'data'/f'C1D_Z_{kind.upper()}_ORIGINAL_MULTIPLIERS.json').read_text()); gens=[]
 for spec,saved in zip(d['generator_specs'],d['generators']):
  if 'root_pair' in spec:
   i,j=spec['root_pair']; z=g[i]*g[j]*(small[spec['small']] if 'small' in spec else R.one)
  else: z=small[spec['small']]*g[spec['root']]*(b if kind=='x' else R.one)
  ck(z==readterms(saved),'original generator mismatch');gens.append(z)
 target=original_part(kind)-(8*sum((z*z for z in g),R.zero) if kind=='scalar' else R.zero)
 actual=sum((f*readterms(m) for f,m in zip(gens,d['multipliers'])),R.zero)
 ck(actual==target,f'{kind}: original-source membership identity')
 for gg,ll in zip(d['basis'],d['basis_original_generator_labels']):
  ck(sum((f*readterms(m) for f,m in zip(gens,ll)),R.zero)==readterms(gg),'label identity')
 for ts in d['multipliers']:
  ck(all(len(mm)==8 and min(mm)>=0 and s.Rational(cc).q>0 for mm,cc in ts),'polynomial multiplier')
 reports.append({'target':kind,'original_target_terms':len(target),'multiplier_terms':sum(map(len,d['multipliers']))})
print('PASS three complete original-generator allocations',time.time()-start,flush=True)
Dg=json.loads((ROOT/'data/C1D_Z_COMPLETE_SOURCE_GRAM.json').read_text()); Q=[[R.zero for _ in range(5)] for _ in range(5)]
for z in Dg['upper_entries']:
 Q[z['row']][z['col']]=Q[z['col']][z['row']]=readterms(z['terms'])
B=X+u-a*a*b*b*q*q*E*v*W;N=Y-r*u*W
# Multiply every root by e*rr*ss, so the complete equation is polynomial.
z=[e*rr*ss*f for f in g]+[b*B*rr*ss,e*N]
poly=sum((z[i]*Q[i][j]*z[j] for i in range(5) for j in range(5)),R.zero)
# Quotient reduction in the three pure coefficient square-root relations.
def radical_reduce(p):
 out=R.zero
 for mm,cc in p.items():
  em,rm,sm=mm[10]//2,mm[11]//2,mm[12]//2
  nn=list(mm);nn[10]%=2;nn[11]%=2;nn[12]%=2;nn[7]+=em;nn[6]+=rm
  out+=R.from_dict({tuple(nn):cc})*(1-r)**sm
 return out
poly=radical_reduce(poly)
# Direct Pauli basis, constructed here independently of the source normal form.
I=[[R.one,R.zero],[R.zero,R.one]];Z=[[R.one,R.zero],[R.zero,-R.one]];PX=[[R.zero,R.one],[R.one,R.zero]]
def add(A,B):return [[a+b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def scale(t,A):return [[t*x for x in row] for row in A]
def kron(A,B):return [[A[i//len(B)][j//len(B[0])]*B[i%len(B)][j%len(B[0])] for j in range(len(A[0])*len(B[0]))]for i in range(len(A)*len(B))]
def lin(items):
 out=[[R.zero]*4 for _ in range(4)]
 for co,mat in items:out=add(out,scale(co,mat))
 return out
Du=1+u*u;Dv=1+v*v
Au=add(scale(1-u*u,Z),scale(2*u,PX));Bv=add(scale(1-v*v,Z),scale(2*v,PX))
Id=scale(Du*Dv,kron(I,I));a0=scale(Du*Dv,kron(Z,I));a1=scale(Dv,kron(Au,I));b0=scale(Du*Dv,kron(I,Z));b1=scale(Du,kron(I,Bv))
c00=scale(Du*Dv,kron(Z,Z));c01=scale(Du,kron(Z,Bv));c10=scale(Dv,kron(Au,Z));c11=kron(Au,Bv)
par={
 'd':lin([(1,Id),(-QQ(1,2),c00),(QQ(1,2),c01),(-QQ(1,2),c10),(-QQ(1,2),c11)]),
 'e2':lin([(1,Id),(1,a1),(-1,b1),(-1,c11)]),
 'e3':lin([(1,Id),(-1,a1),(1,b0),(-1,c10)]),
 'e4':lin([(1,Id),(-1,a0),(1,b0),(-1,c00)]),
 'e5':lin([(1,Id),(-1,a1),(1,b1),(-1,c11)]),
 'alpha':lin([(1,Id),(-1,a0),(-1,b1),(1,c01)]),
 'beta':lin([(1,Id),(1,a1),(-1,b0),(-1,c10)])}
fixture=json.loads((ROOT/'sources/X_BOUNDED_V_SMALL_FULL_POLYNOMIALS.json').read_text())
sv=s.symbols('b hnew q r E');b0s,hs,qs,rs,Es=sv
coeff={}
for key in par:
 expr=s.Poly(s.cancel(Es*s.sympify(fixture[key],locals={str(z):z for z in sv})),*sv)
 acc=R.zero
 for mm,co in expr.terms():
  bb,hh,qq,rrr,ee=mm
  nn=[0]*n;nn[4]=bb+hh;nn[5]=hh;nn[3]=qq;nn[6]=rrr;nn[7]=ee
  acc+=R.from_dict({tuple(nn):QQ.convert(co)})
 coeff[key]=acc
state=[W,X,Y,R.one]
physical=sum((coeff[k]*sum((state[i]*par[k][i][j]*state[j] for i in range(4) for j in range(4)),R.zero)for k in par),R.zero)
ck(poly==physical,'Full five-root identity versus independently reconstructed Pauli operator')
print('PASS complete physical Pauli identity',len(poly),'terms',time.time()-start,flush=True)
# Positive pure coefficient clearing, checked as a polynomial identity.
Dc=coeff['e3'].exquo(2*E)
h=a*b
Dfactor=r*(b*b*r+1)*(b*b*h*h*q*q+1)*(b*b*q*q*r-1)*(-E*h*h*q*q+b*b*h*h*q*q*r+1)*(E*h*h*q**4-h*h*q*q+1)*(E*b*b*h*h*q**4-b*b*h*h*q*q-b*b*r+b*b-1)
ck(Dc==Dfactor,'exact factorization of D')
# Full finite box, with no lower bound on E, r or 1-r.
total=QQ.zero; selected_small=[0,1,2,3,4,8,9]
for i in range(5):
 for j in range(i,5):
  correction=Q[i][j]-(8 if i==j else 0)
  ck(all(any(mm[k] for k in selected_small) for mm in correction),'all correction monomials small')
  total+=sum((abs(co) for co in correction.values()),QQ.zero)
ck(str(total)==Dg['correction_total_l1'],'exact coefficient l1 budget')
eps=QQ.one/(2*(1+total));ck(str(eps)==Dg['uniform_box_epsilon'],'exact epsilon')
ck(0<eps<QQ(1,2),'valid positivity box')
# All coefficients are nonnegative at a concrete positive point of the box.
point=[QQ.zero]*n
point[3]=point[4]=eps/2;point[5]=QQ(1,2);point[6]=QQ(1,2);point[7]=QQ(1,2)
def eval_rational(p,x):
 return sum((co*__import__('functools').reduce(lambda aa,bb:aa*bb,(xx**n for xx,n in zip(x,mm)),QQ.one)for mm,co in p.items()),QQ.zero)
example={}
for key,pp in coeff.items():
 val=eval_rational(pp,point);ck(val>0,'strict physical coefficient example: '+key);example[key]={'positive':True,'numerator_digits':len(str(val.numerator)),'denominator_digits':len(str(val.denominator))}
report={'status':'PASS','exact_assertions':count,'allocations':reports,'physical_polynomial_terms':len(poly),
 'epsilon':str(eps),'coefficient_l1':str(total),'strict_positive_coefficient_example':{'b':str(eps/2),'q':str(eps/2),'a':'1/2','r':'1/2','E':'1/2','coefficients':example},
 'proof_scope':'complete full physical C1D Z=1 identity and uniform local parameter rectangle; no universal master coverage',
 'elapsed_seconds':time.time()-start}
(OUT/'SOURCE_ALLOCATION_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS',count,'exact original-source and full-parent checks',time.time()-start,flush=True)
