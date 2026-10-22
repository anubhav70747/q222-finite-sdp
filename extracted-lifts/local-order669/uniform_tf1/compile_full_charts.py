#!/usr/bin/env python3
"""Evaluate fixed full-parent Gram cones on both state-sign charts.
Uses the proved twelve-root principal identity and literal retained root support.
All pencil entries are rational; no SDP solve or guessed feature extension.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import json,hashlib
OUT=Path(__file__).resolve().parent; IN=OUT.parent/'c1_fixed_family'
base=json.loads((OUT/'C1_FIXED_PRINCIPAL_FEATURES.json').read_text())
ret=json.loads((IN/'C1_RETAINED_ROOTS.json').read_text())
N=8;VARS=base['variables']; ZERO=(0,)*N
DELTA=F(1,4096)
def mono(**kw):return tuple(kw.get(v,0) for v in VARS)
def ep(x,y):return tuple(a+b for a,b in zip(x,y))
def add(d,e,k=F(1)):
 for j,c in e.items():
  d[j]=d.get(j,F(0))+k*c
  if not d[j]:del d[j]
def mul(p,q):
 r={}
 for e,c in p.items():
  for f,d in q.items():add(r,{ep(e,f):c*d})
 return r
def scale(p,c):return {e:F(c)*v for e,v in p.items() if c*v}
def polyadd(*ps):
 r={}
 for p in ps:add(r,p)
 return r
def sq(p):return mul(p,p)
def red(e):
 p={e:F(1)}
 for gate,bas in ((6,0),(7,1)):
  while any(v[gate]>=2 for v in p):
   r={}
   for v,c in p.items():
    if v[gate]<2:add(r,{v:c});continue
    w=list(v);w[gate]-=2;add(r,{tuple(w):c});w[bas]+=4;add(r,{tuple(w):-c})
   p=r
 return p
def fs(x):return str(x.numerator) if x.denominator==1 else f'{x.numerator}/{x.denominator}'
def norm_old(e):
 # quarter_u=2^-6*a², quarter_w=2^-3*btest.
 return (2*e[0],e[1],e[2],e[3],e[4],e[5],0,0),F(1,2**(6*e[0]+3*e[1]))
def solve_columns(rows,n,m):
 piv={}
 for i,row in enumerate(rows):
  v=row.copy();rhs={i:F(1)}
  while v:
   p=min(v)
   if p not in piv:
    s=v[p];piv[p]=({j:c/s for j,c in v.items()},{j:c/s for j,c in rhs.items()});break
   vv,rr=piv[p];s=v[p];add(v,vv,-s);add(rhs,rr,-s)
  else:
   assert not rhs,('dependent target coefficient row',i)
 assert len(piv)==n
 free=[j for j in range(m) if j not in piv]
 exp={j:{n+i:F(1)} for i,j in enumerate(free)}
 for p in sorted(piv,reverse=True):
  rr,rhs=piv[p];v=rhs.copy()
  for j,c in rr.items():
   if j!=p:add(v,exp[j],-c)
  exp[p]=v
 for i,row in enumerate(rows):
  v={}
  for j,c in row.items():add(v,exp[j],c)
  assert v=={i:F(1)}
 return exp,free,piv

def target(sign):
 one={ZERO:F(1)}
 u={mono(a=8):DELTA**2};w={mono(btest=4):DELTA}
 W={mono(omega=4):F(1)};Z={mono(zeta=4):F(sign)}
 X={mono(X=1):F(1)};Y={mono(Y=1):F(1)}
 du=polyadd(one,sq(u));dw=polyadd(one,sq(w));uw=mul(u,w)
 r1=polyadd(W,scale(mul(w,X),-1),mul(u,Y),scale(mul(uw,Z),-1))
 r2=polyadd(mul(u,W),scale(Y,-1))
 r3=polyadd(mul(uw,W),mul(u,X),scale(mul(w,Y),-1),scale(Z,-1))
 r4=polyadd(Y,scale(mul(w,Z),-1))
 r5=polyadd(X,mul(u,Z))
 p3=sq(r3);p4=mul(du,sq(r4));p5=mul(dw,sq(r5))
 d=polyadd(scale(p3,2),scale(p4,2),scale(p5,2),scale(mul(mul(du,dw),sq(Z)),-2))
 return [d,scale(sq(r1),4),scale(mul(dw,sq(r2)),4),scale(mul(mul(du,dw),sq(Y)),4),scale(p3,4),scale(p4,4),scale(p5,4)]

receipts={}
for sign in [1,-1]:
 label='POS' if sign==1 else 'NEG'; chart=ret['charts'][str(sign)]
 blocks=[]
 for b in base['blocks']:
  # Last principal feature is unused when d_i=B_i on WZ<=0.
  fs0=b['features'] if sign==1 else b['features'][:-1]
  blocks.append({'block':len(blocks),'order':len(fs0),'category':'principal','features':fs0,'source_matrix_sum':b['source_matrix_sum'],'multiplier':b['multiplier']})
 qmin='1/2-c*dG/S0*(2/lam+k+1/16777216)'
 qslope='c*dG/(S0*16777216)'
 for root in chart['roots']:
  features=[];coeff=[]
  for fi,c in root['entries']:
   e,sc=norm_old(chart['feature_exponents'][fi]);features.append(e);coeff.append(f'({fs(sc)})*({c})')
  branches=[(ZERO,'1')]
  if root['category']=='circuit_unweighted':
   branches=[(ZERO,qmin),(mono(xi=1),qslope),(mono(a=2,xi=1),qslope)]
  for branch,weight in branches:
   bfeatures=[list(ep(e,branch)) for e in features]
   blocks.append({'block':len(blocks),'order':len(features),'category':root['category'],'name':root['name'],'features':bfeatures,'source_vector':coeff,'source_nonnegative_weight':weight,'multiplier':list(branch)})
 cols=[];entry=[];ambient=set();distinct=set()
 for b in blocks:
  fs0=[tuple(e) for e in b['features']];distinct.update(fs0)
  for i in range(b['order']):
   for j in range(i,b['order']):
    col=len(cols);cols.append([b['block'],i,j])
    for e,c in red(ep(fs0[i],fs0[j])).items():
     entry.append([e,col,c*(1 if i==j else 2)]);ambient.add(e)
 ts=target(sign)
 assert all(e in ambient for p in ts for e in p),'Target outside the derived feature support.'
 ambient=sorted(ambient);rowidx={e:i for i,e in enumerate(ambient)}
 coo=[[rowidx[e],col,fs(c)] for e,col,c in entry]
 n=len(ambient);m=len(cols);rows=[{} for _ in range(n)]
 for i,j,c in coo:rows[i][j]=F(c)
 ex,free,piv=solve_columns(rows,n,m)
 target_rows=[{} for _ in range(n)]
 for c,p in enumerate(ts):
  for e,v in p.items():target_rows[rowidx[e]][c]=v
 # Substitute p=T*c to obtain actual source-variable pencil (7 canonical coefficients).
 source_ex={}
 for col,e in ex.items():
  r={}
  for v,c in e.items():
   if v<n:add(r,target_rows[v],c)
   else:add(r,{7+v-n:c})
  source_ex[col]=r
 # Independent product verification of source pencil: A G(c,y)=T c.
 for i,row in enumerate(rows):
  r={}
  for j,c in row.items():add(r,source_ex[j],c)
  assert r==target_rows[i],(sign,i,r,target_rows[i])
 pcoo=[]
 for col,e in sorted(source_ex.items()):
  b,i,j=cols[col]
  for v,c in sorted(e.items()):pcoo.append([b,i,j,v,fs(c)])
 orders=[b['order'] for b in blocks]
 counts={'blocks':len(blocks),'block_orders':orders,'total_matrix_order':sum(orders),'max_block_order':max(orders),'feature_slots':sum(orders),'distinct_features':len(distinct),'gram_variables':m,'ambient_coefficients':n,'coefficient_map_rank':len(piv),'auxiliary_dimension':len(free),'nonzero_coefficient_map_entries':len(coo),'nonzero_affine_pencil_entries':len(pcoo)}
 features={'schema':'Q222.C1.full_local_fixed_features.v1','state_sign':sign,'variables':VARS,'testing_substitution':{'u':'delta^2*a^8','w':'delta*btest^4','W':'omega^4','Z':f'{sign}*zeta^4'},'delta':'1/4096','testing_equations':base['testing_equations'],'normalizer':base['normalizer'],'blocks':blocks,'counts':counts,'claim':'Fixed full-parent local certificate cone; not a full-Q222 lift.'}
 matrix={'schema':'Q222.gram_coefficient_map.v1','row_monomials':[list(e) for e in ambient],'column_gram_entries':cols,'shape':[n,m],'coo':coo,'target_source_order':['d','m','C','A','D','alpha','beta'],'target_coo':[[row,c,fs(v)] for row,r in enumerate(target_rows) for c,v in sorted(r.items())]}
 pencil={'schema':'Q222.C1.full_local_canonical_source_pencil.v1','state_sign':sign,'claim':'Exact representation of membership in this fixed full-parent Gram cone on the stated local testing chart. Not a full-Q222 lift.','equivalence':'A canonical source belongs to this fixed certificate cone iff the displayed rational pencil has a PSD completion.','source_order':['d','m','C','A','D','alpha','beta'],'auxiliary_order':[f'y_{i}' for i in range(len(free))],'block_orders':orders,'total_matrix_order':sum(orders),'constant_entries':[],'coo':pcoo,'convention':'[block,i,j,variable,rational coefficient]; variable0..6 are canonical source coefficients and remaining variables auxiliaries. i<=j stores actual symmetric entries; mirror off-diagonal.','counts':counts}
 (OUT/f'C1_{label}_FULL_FIXED_FEATURES.json').write_text(json.dumps(features,indent=2)+'\n')
 (OUT/f'C1_{label}_FULL_COEFFICIENT_MAP.coo.json').write_text(json.dumps(matrix,indent=2)+'\n')
 (OUT/f'C1_{label}_FULL_SOURCE_PENCIL.coo.json').write_text(json.dumps(pencil,indent=2)+'\n')
 receipt={'status':'passed','counts':counts,'exact_checks':{'all_coefficient_equations':n,'source_target_equations':n,'free_gram_coordinates':len(free)},'input_retained_sha256':hashlib.sha256((IN/'C1_RETAINED_ROOTS.json').read_bytes()).hexdigest(),'qualification':'Full-parent symbolic source reassembly and uniform source admission are verified in separate mathematical files; this receipt verifies fixed features, matrix compilation, and coefficient elimination.'}
 receipts[label]=receipt
print(json.dumps(receipts,indent=2))
(OUT/'full_charts_compilation_receipt.json').write_text(json.dumps(receipts,indent=2)+'\n')
