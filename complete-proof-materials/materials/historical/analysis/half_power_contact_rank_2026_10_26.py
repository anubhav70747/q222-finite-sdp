"""Exact coefficient ranks for half-power contact faces; no solver."""
import json
from pathlib import Path
import sympy as S

def allowed(family,i,j,l):
    if l==0 and i+j<2:return False
    if family=='C2' and l==3 and i+j>2:return False
    if family=='C1' and l==1 and i>j:return False
    return True

def product(left,right):
    out={}
    for (i,j,l),a in left:
        for (k,h,m),b in right:
            assert (i+k)%2==0 and (j+h)%2==0
            key=((i+k)//2,(j+h)//2,min(l,m),max(l,m))
            out[key]=out.get(key,0)+a*b
    return {key:value for key,value in out.items() if value}

def compute(family):
    blocks=[]
    for parity in [(0,0),(0,1),(1,0),(1,1)]:
        blocks.append([(i,j,l) for l in range(4) for i in range(3)
                       for j in range(3) if (i%2,j%2)==parity and allowed(family,i,j,l)])
    cols=[]
    for block in blocks:
        roots=[[(feature,1),(block[-1],-1)] for feature in block[:-1]]
        for i,r in enumerate(roots):
            for j in range(i,len(roots)):
                coef=product(r,roots[j])
                cols.append({key:value*(1 if i==j else 2) for key,value in coef.items()})
    rows=[(a,b,l,m) for a in range(3) for b in range(3) for l in range(4) for m in range(l,4)]
    mat=S.Matrix([[col.get(key,0) for col in cols] for key in rows])
    zero_rows=[key for i,key in enumerate(rows) if not any(mat[i,j] for j in range(mat.cols))]
    rank=mat.rank()
    assert rank==72
    assert [len(b) for b in blocks]==[14,6,6,4]
    assert len(cols)==127 and len(zero_rows)==12
    return {'family':family,'raw_support_blocks':[len(b) for b in blocks],
            'contact_blocks':[len(b)-1 for b in blocks],
            'symmetric_gram_variables':len(cols),'coefficient_rows':len(rows),
            'identically_absent_rows':zero_rows,'nonzero_rows':len(rows)-len(zero_rows),
            'exact_rank':rank,'gram_kernel_dimension':len(cols)-rank,
            'scope':'Linear coefficient-map rank after monomial/contact normalization; no PSD feasibility or admission asserted.'}

result=[compute(f) for f in ['C1','C2']]
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
