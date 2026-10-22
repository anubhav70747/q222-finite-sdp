#!/usr/bin/env python3
"""Evaluate an affine rational pencil for the already fixed principal Gram cone."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import json
OUT=Path(__file__).resolve().parent
A=json.loads((OUT/'C1_FIXED_PRINCIPAL_COEFFICIENT_MAP.coo.json').read_text())
n,m=A['shape'];rows=[{} for _ in range(n)]
for i,j,c in A['coo']:rows[i][j]=F(c)
def add(d,e,k=F(1)):
    for j,c in e.items():
        d[j]=d.get(j,F(0))+k*c
        if not d[j]:del d[j]
piv={}
for i,row in enumerate(rows):
    v=row.copy();rhs={i:F(1)}
    while v:
        p=min(v)
        if p not in piv:
            s=v[p];v={j:c/s for j,c in v.items()};rhs={j:c/s for j,c in rhs.items()}
            piv[p]=(v,rhs);break
        s=v[p];vv,rr=piv[p];add(v,vv,-s);add(rhs,rr,-s)
    else:
        assert not rhs,'Full row rank expected.'
assert len(piv)==n
free=[j for j in range(m) if j not in piv]
freeidx={j:i for i,j in enumerate(free)}
exp={j:{n+i:F(1)} for j,i in freeidx.items()}
for p in sorted(piv,reverse=True):
    row,rhs=piv[p];v=rhs.copy()
    for j,c in row.items():
        if j==p:continue
        assert j in exp
        add(v,exp[j],-c)
    exp[p]=v
# Independent direct coefficient multiplication, not an inverse assertion.
for i,row in enumerate(rows):
    r={}
    for j,c in row.items():add(r,exp[j],c)
    assert r=={i:F(1)},(i,r)
# Assert the free columns equal their independent named variables.
assert all(exp[j]=={n+i:F(1)} for j,i in freeidx.items())
def fs(x):return str(x.numerator) if x.denominator==1 else f'{x.numerator}/{x.denominator}'
coo=[]
for col,expr in sorted(exp.items()):
    b,i,j=A['column_gram_entries'][col]
    for v,c in sorted(expr.items()):coo.append([b,i,j,v,fs(c)])
pencil={'schema':'Q222.C1.fixed_principal_affine_pencil.v1',
 'claim':'Exact evaluated affine pencil for the derived fixed principal Gram cone, not for Q222.',
 'equivalence':'p is in the 16-block principal Gram cone iff exists y in R^968 such that every listed 12x12 block is PSD.',
 'block_orders':[12]*16,'total_matrix_order':192,'output_dimension':n,'auxiliary_dimension':len(free),
 'constant_entries':[], 'variable_order':{'outputs':[f'p_{i}' for i in range(n)],'auxiliaries':[f'y_{i}' for i in range(len(free))]},
 'output_monomials':A['row_monomials'],'stored_entries':'[block,i,j,variable_index,rational_coefficient], i<=j, actual matrix entries; mirror off-diagonal.',
 'coo':coo,'pivot_gram_columns':sorted(piv),'free_gram_columns':free}
(OUT/'C1_PRINCIPAL_GRAM_CONE_PENCIL.coo.json').write_text(json.dumps(pencil,indent=2)+'\n')
receipt={'status':'passed','equations_verified_exactly':n,'free_coordinates_verified_exactly':len(free),'nonzero_affine_entries':len(coo),'maximum_coefficient_denominator':max(c.denominator for e in exp.values() for c in e.values()),'total_matrix_order':192,'output_dimension':n,'auxiliary_dimension':len(free)}
(OUT/'principal_affine_pencil_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
