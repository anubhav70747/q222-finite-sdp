"""Independent exact checks of the angular cone pencil; standard library only.

This checks the published linear matrices and explicit rational feasible
witnesses. The accompanying proof establishes universal cone equality.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import gcd, isqrt
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
p = json.loads((HERE/'ANGULAR_CONE_ORDER49.coo.json').read_text())
E = [(a,b) for a in [-1,0,1] for b in [-1,0,1] if a or b]
P = [[e,(-e[0],-e[1])] for e in [(1,0),(0,1),(1,1),(1,-1)]]
T = [[(a,b),(-a,0),(0,-b)] for a in [-1,1] for b in [-1,1]]
W = [[(-1,1),(-1,-1),(1,0)],[(1,1),(1,-1),(-1,0)],
     [(1,-1),(-1,-1),(0,1)],[(1,1),(-1,1),(0,-1)]]
assert p['order']==49 and [b['order'] for b in p['blocks']]==[2]*24+[1]
assert len(p['input_variables'])==9 and len(p['auxiliary_variables'])==48
assert len(p['entries'])==109

def cross(a,b): return a[0]*b[1]-a[1]*b[0]
def signature(es,bs): return tuple(sorted(zip(es,bs)))
all_circuits=[]
for a,b in combinations(E,2):
    if a==(-b[0],-b[1]):all_circuits.append(signature([a,b],[1,1]))
for a,b,c in combinations(E,3):
    z=[cross(b,c),cross(c,a),cross(a,b)]
    if all(t<0 for t in z):z=[-t for t in z]
    if all(t>0 for t in z):
        g=gcd(*z);all_circuits.append(signature([a,b,c],[t//g for t in z]))
listed=[signature(es,[1,1]) for es in P]+[signature(es,[1,1,1]) for es in T]+[signature(es,[1,1,2]) for es in W]
assert len(all_circuits)==12 and set(all_circuits)==set(listed)

def direct(v):
    blocks=[]
    def val(n):return v.get(n,Q(0))
    def r(e):
        return val(f'q_{e[0]}_{e[1]}')-sum(val(f'r_{k}_{j}') for k,es in enumerate(T+W) for j,f in enumerate(es) if e==f)
    for k,(e,f) in enumerate(P): blocks.append([2*r(e),val(f's_{k}'),2*r(f)])
    for k in range(4):
        a,b,c=(val(f'r_{k}_{j}') for j in range(3));x=val(f'x_{k}');y=val(f'y_{k}');t=val(f's_{k+4}')/3
        blocks.extend([[a,x,b],[c,y,t],[x,t,y]])
    for k in range(4,8):
        a,b,c=(val(f'r_{k}_{j}') for j in range(3));x=val(f'x_{k}');t=val(f's_{k+4}')/4
        blocks.extend([[a,x,b],[x,t,c/2]])
    blocks.append([val('h')+sum(val(f's_{k}') for k in range(12))])
    return blocks

def from_data(v):
    ans=[[Q(0),Q(0),Q(0)] for _ in range(24)]+[[Q(0)]]
    for i,j,k,a,b in p['entries']:
        assert i<=j
        block=i//2
        assert block==j//2
        pos=0 if i%2==j%2==0 else 2 if i%2==j%2==1 else 1
        ans[block][pos]+=Q(a,b)*v.get(p['variables'][k],Q(0))
    return ans

for name in p['variables']:
    assert from_data({name:Q(1)})==direct({name:Q(1)}),name

def square_root(q):
    assert q>=0
    a,b=isqrt(q.numerator),isqrt(q.denominator)
    assert a*a==q.numerator and b*b==q.denominator
    return Q(a,b)

def make_witness(weights,U=Q(2)**6,V=Q(3)**6,extra=None,flip=False,margin=Q(0)):
    v={name:Q(0) for name in p['variables']}
    def mon(e):return U**e[0]*V**e[1]
    for k,(es,bs) in enumerate(zip(P+T+W,[[1,1]]*4+[[1,1,1]]*4+[[1,1,2]]*4)):
        lam=Q(weights[k]); rr=[lam*b/mon(e) for e,b in zip(es,bs)]
        for e,a in zip(es,rr):v[f'q_{e[0]}_{e[1]}']+=a
        v[f's_{k}']=lam*sum(bs)
        if k>=4:
            j=k-4
            for l,a in enumerate(rr):v[f'r_{j}_{l}']=a
            v[f'x_{j}']=square_root(rr[0]*rr[1])
            if k<8:v[f'y_{j}']=square_root(rr[2]*lam)
        if flip and (k<4 or k>=8):v[f's_{k}']*=-1
    for e,a in (extra or {}).items():v[f'q_{e[0]}_{e[1]}']+=Q(a)
    v['h']=-sum(v[f's_{k}'] for k in range(12))+margin
    blocks=from_data(v)
    assert blocks==direct(v)
    for a,b,c in blocks[:24]:assert a>=0 and c>=0 and a*c-b*b>=0
    assert blocks[24][0]>=0
    return v

fixtures=[]
for k in range(12):
    weights=[0]*12;weights[k]=Q(k+1,k+2)
    fixtures.append(make_witness(weights))
fixtures.append(make_witness([Q(k+1,k+2) for k in range(12)]))
fixtures.append(make_witness([0]*12,extra={e:Q(j+1) for j,e in enumerate(E)},margin=Q(2)))
# An opposite pair on a proper face, with an off-face monomial whose
# contribution tends to zero along V -> 0: no finite global minimizer.
fixtures.append(make_witness([1]+[0]*11,extra={(0,1):Q(7)}))
fixtures.append(make_witness([1]*12,flip=True,margin=Q(1,7)))
fixtures.append(make_witness([0]*12))
fixtures.append(make_witness([Q(1,5)]*12,U=Q(1,2)**6,V=Q(5,2)**6))
for v in fixtures:
    for U,V in [(Q(1),Q(1)),(Q(2),Q(3)),(Q(1,5),Q(7)),(Q(2)**6,Q(3)**6)]:
        value=v['h']+sum(v[f'q_{a}_{b}']*U**a*V**b for a,b in E)
        assert value>=0

report={'scope':'Exact matrix reconstruction and rational witnesses; universal cone equality is proved analytically in angular_cone.tex.',
        'all_primitive_positive_circuits':12,'input_columns_checked':9,
        'auxiliary_columns_checked':48,'rational_coo_records':109,
        'order2_blocks':24,'scalar_blocks':1,'total_order':49,
        'exact_feasible_fixtures':len(fixtures),'boundary_and_zero_weights':True,
        'proper_face_with_minimum_only_at_infinity':True,
        'signed_pair_and_weighted_thresholds':True,'all_checks_passed':True}
(HERE/'angular_pencil_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(report,indent=2))
