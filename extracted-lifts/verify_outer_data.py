"""Exact rational checks of stored outer data and the two-frame separator.

This does not test quantum membership or prove either open reverse inclusion.
"""
from pathlib import Path
from fractions import Fraction as F
import json

HERE=Path(__file__).resolve().parent/'outer-candidates'
result={}
for name in ['STRICT_OUTER_N8368.coo.json','STRICT_OUTER_N10736.coo.json',
             'OUTER_REVERSE_OPEN_N13744.coo.json','OUTER_REVERSE_OPEN_C160_N63888.coo.json']:
    d=json.loads((HERE/name).read_text());N=d['matrix_order'];m=d['auxiliary_count']
    if 'C160' in name:
        sizes=dict(d['blocks']);assert sum(sizes.values())==N
        for b,i,j,k,p,q in d['entries']:
            assert 0<=i<=j<sizes[b] and 0<=k<9+m
            assert F(p,q)!=0
    else:
        owner=[None]*N
        for b,block in enumerate(d['blocks']):
            off=block['offset'];r=block['order']
            assert off>=0 and off+r<=N
            for i in range(off,off+r):assert owner[i] is None;owner[i]=b
        assert all(v is not None for v in owner)
        for k,i,j,v in d['entries']:
            assert 0<=i<=j<N and 0<=k<9+m and owner[i]==owner[j]
            assert F(v)!=0
    assert d['quantum_reverse_inclusion_proved'] is False
    result[name]={'order':N,'auxiliaries':m,'blocks':len(d['blocks']),
                  'rational_records':len(d['entries']),'scope':d['status']}

sep=json.loads((HERE/'strict_outer_separator.json').read_text())
coo=json.loads((HERE/'C1_certificate_source_pencil.coo.json').read_text())
source=list(map(F,json.loads((HERE/'strict_outer_physical_source.json').read_text())['source']))
ell=list(map(F,sep['common_separator']));assert sum(a*b for a,b in zip(source,ell))==-1
for orientation in sep['orientation_certificates']:
    blocks=[[[F(x) for x in row]for row in mat]for mat in orientation['dual_blocks']]
    for mat in blocks:
        assert all(mat[i][j]==mat[j][i] for i in range(len(mat)) for j in range(len(mat)))
        a=[row[:]for row in mat]
        for k in range(len(a)):
            pivot=a[k][k];assert pivot>0
            for i in range(k+1,len(a)):
                for j in range(i,len(a)):
                    a[i][j]-=a[i][k]*a[k][j]/pivot;a[j][i]=a[i][j]
    actual=[F(0)]*(1+len(coo['variables']))
    for k,i,j,v in coo['entries']:
        for b,block in enumerate(coo['blocks'][:4]):
            off=block['offset'];size=block['order']
            if off<=i<=j<off+size:
                actual[k]+=F(v)*blocks[b][i-off][j-off]*(1 if i==j else 2)
                break
    perm=list(range(7)) if orientation['orientation']==0 else [0,2,1,5,6,3,4]
    expected=[F(0)]+[ell[i]for i in perm]+[F(0)]*coo['auxiliary_count']
    expected[1]-=1
    assert actual==expected
result['strict_outer_separator']={'source_pairing':'-1','positive_dual_blocks':8,
                                 'all_visible_and_auxiliary_trace_coefficients_equal':True}
print(json.dumps(result,indent=2))
