"""Replay one rational half-line witness using Python's standard library.

Reads only the literal witness and the exact affine pencil.  The numerical
diagnostic and the witness-recovery script are not imported or executed.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib
import json

HERE=Path(__file__).resolve().parent
path=HERE/'SIMPLE_SOURCE_HALFLINE_EXACT_WITNESS.json'
w=json.loads(path.read_text())
pencil_path=HERE/'HALFLINE_QUADRATIC_SELECTOR_ORDER24.coo.json'
p=json.loads(pencil_path.read_text())
assert hashlib.sha256(pencil_path.read_bytes()).hexdigest()==w['original_pencil_sha256']

def matrix(rows):return [[Q(a) for a in row] for row in rows]
def zero(n,m=None):return [[Q(0) for _ in range(n if m is None else m)] for _ in range(n)]
def tr(A):return [list(a) for a in zip(*A)]
def add(A,B):return [[a+b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def scale(c,A):return [[c*a for a in row] for row in A]
def mul(A,B):
    return [[sum(a*b for a,b in zip(row,col)) for col in tr(B)] for row in A]
def block(A,B,C,D):return [a+b for a,b in zip(A,B)]+[c+d for c,d in zip(C,D)]
def vecmul(A,b):return [sum(a*v for a,v in zip(row,b)) for row in A]
def extract(A,r,c):return [[A[i][j] for j in c] for i in r]
def diag(a):return [[v if i==j else Q(0) for j in range(len(a))] for i,v in enumerate(a)]

source=list(map(Q,w['source']));aux=list(map(Q,w['auxiliary']))
assert source==list(map(Q,['359','0','0','75/16','60','416','41/2','831/2','151/2']))
assert all(a>=0 for a in source)
assert [Q(a) for a in w['source_d1']]==[a/359 for a in source]
assert len(aux)==44 and p['auxiliary_order']==w['auxiliary_order']
J=[matrix([[aux[4*k+2*i+j] for j in range(2)] for i in range(2)]) for k in range(3)]
assert J==[scale(Q(1,64),matrix(a)) for a in w['selector_integer_coefficients']]
BG=matrix(w['BG']);BH=matrix(w['BH'])
assert BG==[aux[12+4*i:16+4*i] for i in range(4)]
assert BH==[aux[28+4*i:32+4*i] for i in range(4)]
for name,M in [('BG',BG),('BH',BH)]:
    compressed=w[name+'_compact']
    assert M==scale(1/Q(compressed['denominator']),matrix(compressed['integer_numerator']))

blocks={name:matrix(w['blocks'][name]) for name in ['SG','SH','TG','TH']}
full=zero(24)
for k,i,j,a in p['entries']:full[i][j]+=a*(source+aux)[k]
offset=0
for name,n in zip(blocks,[8,8,4,4]):
    assert blocks[name]==extract(full,range(offset,offset+n),range(offset,offset+n))
    offset+=n

# Independent direct P,N,E,K formulas, evaluated at three rational values.
# Every compared polynomial has degree at most two, so this proves all
# coefficient equations, not merely agreement on selected testing angles.
d,la,nn,m,C,A,D,al,be=source
def expected(v):
    dv=1+v*v
    P=[[4*m*v*v,-4*m*v],[-4*m*v,4*m+(4*be+2*d)*dv]]
    N=[[4*C*dv+4*D+2*d,(4*D+2*d)*v],[(4*D+2*d)*v,(4*D+2*d)*v*v]]
    E=[[4*A*dv+(4*al+2*d)*v*v,-(4*al+2*d)*v],[-(4*al+2*d)*v,4*al-2*d*v*v]]
    K=[[4*la,4*la*v],[4*la*v,4*la*v*v+4*nn*dv]]
    j=add(add(J[0],scale(v,J[1])),scale(v*v,J[2]))
    off=add(add(K,scale(-1,E)),scale(-1,j))
    return block(add(P,K),j,tr(j),add(N,K)),block(add(P,E),off,tr(off),add(N,E))
def recovered(S,T,v):
    R=[[(Q(int(i==j))) for j in range(4)] for i in range(4)]
    R=[row+[v*a for a in row] for row in R]
    return add(mul(mul(R,S),tr(R)),scale(v,T))
for v in map(Q,[-1,0,2]):
    G,H=expected(v)
    assert recovered(blocks['SG'],blocks['TG'],v)==G
    assert recovered(blocks['SH'],blocks['TH'],v)==H

# Independent exact LDL elimination, allowing zero pivots only with zero
# remaining row. This is also a complete rational PSD test for these data.
results={}
for name,M in blocks.items():
    assert M==tr(M)
    n=len(M);R=[row[:] for row in M];L=diag([Q(1)]*n);piv=[]
    for i in range(n):
        t=R[i][i];piv.append(t)
        assert t>=0
        if not t:
            assert all(R[i][j]==0 for j in range(i+1,n))
            continue
        for j in range(i+1,n):L[j][i]=R[j][i]/t
        for j in range(i+1,n):
            for k in range(j,n):
                R[j][k]-=R[j][i]*R[i][k]/t
                R[k][j]=R[j][k]
    assert mul(mul(L,diag(piv)),tr(L))==M
    z=w['ldl'][name]
    assert L==matrix(z['L']) and piv==list(map(Q,z['D']))
    assert sum(bool(a) for a in piv)==z['rank']
    assert [a.numerator//a.denominator for a in piv if a>0]==z['positive_pivot_floor']
    results[name]={'order':n,'rank':z['rank'],
                   'positive_pivot_floor':z['positive_pivot_floor'],
                   'zero_pivot_indices':[i for i,a in enumerate(piv) if not a]}

# Endpoint and full contact kernels are exact, not residual tolerances.
g=list(map(Q,[34,4,6,-9]));h=list(map(Q,[17,2,12,-18]))
assert all(a==0 for a in blocks['SG'][0]+blocks['SG'][3])
assert all(a==0 for a in blocks['TG'][0]+blocks['TG'][3])
assert vecmul(blocks['SG'],g+[2*a for a in g])==[0]*8
assert vecmul(blocks['SH'],h+[2*a for a in h])==[0]*8
assert vecmul(blocks['TG'],g)==[0]*4
assert vecmul(blocks['TH'],h)==[0]*4

# The positive scaling to d=1 scales every Gram and selector linearly.
for S,T in [('SG','TG'),('SH','TH')]:
    assert recovered(scale(Q(1,359),blocks[S]),scale(Q(1,359),blocks[T]),Q(2))==scale(Q(1,359),recovered(blocks[S],blocks[T],Q(2)))
out={'scope':'Exact rational witness for this one simple C1 source only.',
     'witness_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
     'exact_coefficient_equalities':96,'exact_blocks':results,
     'all_four_blocks_psd':True,'solver_called':False,
     'uses_numerical_diagnostic':False}
(HERE/'simple_source_exact_witness_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
