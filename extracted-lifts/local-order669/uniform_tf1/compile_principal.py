#!/usr/bin/env python3
"""Exact fixed-family compiler for the complete C1 principal twelve-root Gram.
No numerical SDP or prior archive verifier is used. All coefficients are integers.
The source-side PSD certificate is described in UNIFORM_TF1_PROOF.txt.
"""
from pathlib import Path
from fractions import Fraction
from collections import defaultdict
import json, hashlib
OUT=Path(__file__).resolve().parent
IN=OUT.parent/'c1_fixed_family'
VARS=('a','btest','omega','zeta','X','Y','xi','eta')
N=len(VARS)
def mono(**kw):
    return tuple(kw.get(v,0) for v in VARS)
def addexp(x,y):return tuple(a+b for a,b in zip(x,y))
def times(p,q):
    r=defaultdict(Fraction)
    for e,c in p.items():
      for f,d in q.items():r[addexp(e,f)]+=c*d
    return {e:c for e,c in r.items() if c}
def reduce_gate(e):
    """Canonical reduction by xi²=1-a⁴ and eta²=1-btest⁴."""
    p={e:Fraction(1)}
    for gate,base in ((6,0),(7,1)):
      while any(v[gate]>=2 for v in p):
        out=defaultdict(Fraction)
        for v,c in p.items():
          if v[gate]>=2:
            w=list(v);w[gate]-=2;out[tuple(w)]+=c
            w[base]+=4;out[tuple(w)]-=c
          else:out[v]+=c
        p={v:c for v,c in out.items() if c}
    return p
# Primitive rescaled principal testing vector. φ=Aχ with source diagonal scales.
CHI=(mono(X=1),mono(Y=1),mono(zeta=4),mono(omega=4),
     mono(a=8,zeta=4),mono(btest=4,zeta=4),mono(a=8,omega=4),
     mono(a=8,btest=4,omega=4),mono(a=16,omega=4),
     mono(a=12,omega=4),mono(a=4,omega=4),
     mono(a=4,btest=2,omega=2,zeta=2))
expected={(0,0),(0,1),(0,2),(1,0),(2,0),(2,1),(4,0),(4,1)}
gram=json.loads((IN/'C1_PRINCIPAL_GRAM.json').read_text())
monos={tuple(d['monomial_exponents']) for d in gram['coefficient_matrices']}
assert monos==expected
blocks={mono():{'id':'constant','uses':['residual_source_matrix']}}
channel_map=[]
for p,q in sorted(monos-{(0,0)}):
    plus=mono(a=2*p,btest=2*q)
    minus=[mono(a=2*i,xi=1) for i in range(p)]
    minus += [mono(a=2*p,btest=2*j,eta=1) for j in range(q)]
    blocks[mono()]['uses'].append(f'Hplus_{p}_{q}')
    for e,label in [(plus,f'Hplus_{p}_{q}')]+[(m,f'Hminus_{p}_{q}') for m in minus]:
        blocks.setdefault(e,{'id':None,'uses':[]})['uses'].append(label)
    # Direct exact algebraic checks of the only needed scalar identities.
    tau=mono(a=4*p,btest=4*q)
    P={mono():Fraction(1),tau:Fraction(1)}
    rhs=defaultdict(Fraction)
    for e in [mono(),plus]:
      for ee,c in reduce_gate(addexp(e,e)).items():rhs[ee]+=c
    assert {e:c for e,c in rhs.items() if c}==P
    P={mono():Fraction(1),tau:Fraction(-1)}
    rhs=defaultdict(Fraction)
    for e in minus:
      for ee,c in reduce_gate(addexp(e,e)).items():rhs[ee]+=c
    assert {e:c for e,c in rhs.items() if c}==P
    channel_map.append({'monomial':[p,q], 'source_weight':f'1/{4096**(p+q)}',
                        'plus':[list(mono()),list(plus)],'minus':[list(x) for x in minus]})
# Stable ordering: constant, un-gated monomials, xi gates, eta gates.
mult=sorted(blocks,key=lambda e:(e[6]+e[7],e[7],e[0]+e[1],e))
assert len(mult)==16
feature_blocks=[]; columns=[]; entries=[]; ambient=set(); unique=set()
for bi,m in enumerate(mult):
    fs=[addexp(m,c) for c in CHI];unique.update(fs)
    feature_blocks.append({'block':bi,'order':12,'multiplier':list(m),
                           'source_matrix_sum':blocks[m]['uses'],
                           'features':[list(e) for e in fs]})
    for i in range(12):
      for j in range(i,12):
        col=len(columns);columns.append([bi,i,j])
        p=reduce_gate(addexp(fs[i],fs[j]))
        for e,c in p.items():
          cc=c*(1 if i==j else 2);assert cc.denominator==1
          entries.append([e,col,cc.numerator]);ambient.add(e)
ambient=sorted(ambient);rows={e:i for i,e in enumerate(ambient)}
coo=[[rows[e],col,c] for e,col,c in entries]
# Exact row rank, exploiting small sparse columns.
pivots={}
for row in range(len(ambient)):
    v={j:Fraction(c) for i,j,c in coo if i==row}
    while v:
      p=min(v)
      if p not in pivots:
        scale=v[p];pivots[p]={j:c/scale for j,c in v.items()};break
      fac=v[p]
      for j,c in pivots[p].items():
        v[j]=v.get(j,Fraction(0))-fac*c
        if not v[j]:del v[j]
rank=len(pivots)
# Moment matrices: coefficient of m_row in each 12x12 block. Store i<=j,
# with no doubled off-diagonal convention, i.e. product of two functions.
mom=[]
for bi,b in enumerate(feature_blocks):
  arr=[]
  for i in range(12):
    for j in range(i,12):
      for e,c in reduce_gate(addexp(tuple(b['features'][i]),tuple(b['features'][j]))).items():
        arr.append([rows[e],i,j,c.numerator])
  mom.append({'block':bi,'order':12,'entries':arr})
raw=json.dumps({'variables':VARS,'blocks':feature_blocks,'monomials':ambient,'coo':coo},sort_keys=True).encode()
base={'schema':'Q222.C1.uniform_principal_gram_cone.v1','scope':'Admitted bounded-positive-z C1 sources, lambda>0; W,Z>=0 testing chart; principal twelve-root part only. Full parent retained roots are separate.',
 'delta':'1/4096','testing_substitution':{'u':'delta^2*a^8','w':'delta*btest^4','W':'omega^4','Z':'zeta^4'},
 'testing_equations':['xi^2+a^4-1=0','eta^2+btest^4-1=0','omega^8+zeta^8+X^2+Y^2-1=0'],
 'normalizer':'kappa_test^2*(1+delta^4*a^16)*(1+delta^2*btest^8)=1',
 'variables':VARS,'principal_chi':[list(x) for x in CHI],
 'phi_from_chi_diagonal':['1','1','1','1',f'1/{4096**2}',f'1/{4096}',f'1/{4096**2}',f'1/{4096**3}',f'1/{4096**4}',f'1/{4096**3}',f'1/{4096}',f'1/{2**18}'],
 'blocks':feature_blocks,'channel_map':channel_map,
 'counts':{'blocks':16,'block_order':12,'total_block_order':192,'feature_slots':192,'distinct_features':len(unique),'gram_variables':len(columns),'ambient_coefficients':len(ambient),'coefficient_map_rank':rank,'gram_kernel_dimension':len(columns)-rank,'nonzero_map_entries':len(coo)},
 'input_gram_sha256':hashlib.sha256((IN/'C1_PRINCIPAL_GRAM.json').read_bytes()).hexdigest(),
 'matrix_data_sha256':hashlib.sha256(raw).hexdigest()}
(OUT/'C1_FIXED_PRINCIPAL_FEATURES.json').write_text(json.dumps(base,indent=2)+'\n')
(OUT/'C1_FIXED_PRINCIPAL_COEFFICIENT_MAP.coo.json').write_text(json.dumps({'schema':'Q222.gram_coefficient_map.v1','row_monomials':[list(x) for x in ambient],'column_gram_entries':columns,'shape':[len(ambient),len(columns)],'coo':coo,'off_diagonal_convention':'Gram coefficient map doubles i<j; each variable denotes the actual symmetric matrix entry.'},indent=2)+'\n')
(OUT/'C1_FIXED_PRINCIPAL_MOMENT_PENCIL.coo.json').write_text(json.dumps({'schema':'Q222.block_moment_pencil.v1','variables':[f'm_{i}' for i in range(len(ambient))],'row_monomials':[list(x) for x in ambient],'blocks':mom,'off_diagonal_convention':'Stored i<=j actual symmetric entries; mirror for i>j. Each [k,i,j,c] adds c*m_k.'},indent=2)+'\n')
receipt={'status':'passed','checks':{'source_gram_monomial_support':True,'fourteen_scalar_polynomial_identities':True,'integer_coefficient_map':True,'exact_rational_row_rank':rank},'counts':base['counts'],'limits':['This does not certify every C1/C2 source.','This does not compile or prove a full-Q222 pencil.','Retained circuit/event/angular functions are not included in this principal-only export.']}
(OUT/'principal_compilation_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
