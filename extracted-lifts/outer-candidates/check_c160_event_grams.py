"""Exact sixteen-event certificates in the fixed 32-dimensional C160 block.

Only rational arithmetic is used.  These certificates prove probability
positivity, not the quantum reverse inclusion of the C160 candidate.
Usage: python check_c160_event_grams.py [directory containing C160 cache]
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import hashlib, json, sys
import sympy as s

HERE = Path(__file__).resolve().parent
DATA = Path(sys.argv[1]) if len(sys.argv)>1 else HERE
CACHE = DATA/'C160_source_compiler_cache.json'
d = json.loads(CACHE.read_text())
refs = d['tensor_phi_basis']['class_reference_coords']
refmap = {c:(r,F(w)) for r,(c,w) in enumerate(refs)}
cols = {c:{r:w} for c,(r,w) in refmap.items()}
for v in d['nullspace_basis']:
    col = {}
    for c,n,den in v['entries']:
        if c == v['free']: continue
        r,w = refmap[c]
        col[r] = col.get(r,F(0))-F(n,den)*w
    cols[v['free']] = {r:w for r,w in col.items() if w}
T = [[F(0) for j in range(9)] for r in range(243)]
for j,vec in enumerate(d['target_affine_section']):
    for c,n,den in vec:
        for r,w in cols[c].items(): T[r][j] += F(n,den)*w
block = d['block_layout'][0]
assert block['free']==list(range(5)) and block['size']==32 and block['coord_start']==0
bits = list(product(range(2), repeat=5))
z = s.symbols('k a0 b0 a1 b1')
k = z[0]
phi = [s.prod(t if bit else 1-t for t,bit in zip(z,b)) for b in bits]
def f(t,sign): return (1-t,t) if sign==1 else (-t,1-t)
QD = s.eye(16).row_join(s.eye(16)).col_join(s.eye(16).row_join(2*s.eye(16)))
events = []
for x,y in product(range(2),repeat=2):
    total = s.zeros(32)
    ia,ib = (1,3)[x],(2,4)[y]
    other = [i for i in range(1,5) if i not in (ia,ib)]
    for alpha,beta in product((1,-1),repeat=2):
        av,bv = f(z[ia],alpha),f(z[ib],beta)
        P = av[0]*bv[0]+k*av[1]*bv[1]
        roots=[]
        for e in product(range(2),repeat=2):
            h = s.Poly(P*s.prod(z[i] if bit else 1-z[i] for i,bit in zip(other,e)),*z)
            assert all(n<=1 for mon,c in h.terms() for n in mon)
            coeff = [int(h.eval(dict(zip(z,b)))) for b in bits]
            assert s.Poly(sum(c*v for c,v in zip(coeff,phi)),*z)==h
            roots.append(coeff)
        R = s.Matrix(roots)
        G = R.T*R
        got = [F(0) for r in range(243)]
        c=0
        for i in range(32):
            for j in range(i,32):
                if G[i,j]:
                    for r,w in cols[c].items(): got[r]+=int(G[i,j])*w
                c+=1
        target = [(T[r][0]-alpha*T[r][1+x]-beta*T[r][3+y]-alpha*beta*T[r][5+2*x+y])/4 for r in range(243)]
        assert got==target, (x,y,alpha,beta)
        total+=G
        events.append({'x':x,'y':y,'alpha':alpha,'beta':beta,'root_rows':roots,'rank':int(R.rank())})
    # Gram representations are not unique.  The sum and QD have the same
    # polynomial, not the same entries (the former has rank at most 16).
    diff=[F(0) for r in range(243)]
    c=0
    for i in range(32):
        for j in range(i,32):
            for r,w in cols[c].items(): diff[r]+=int(total[i,j]-QD[i,j])*w
            c+=1
    assert not any(diff), (x,y)
out={'scope':'Sixteen exact event Gram certificates; quantum reverse inclusion remains unproved',
     'cache_sha256':hashlib.sha256(CACHE.read_bytes()).hexdigest(),
     'basis_bits':bits,'block':0,'gram_order':32,
     'normalization_gram':'[[I16,I16],[I16,2I16]]',
     'equation':'sum_h root_h(z)^2 = D(z) (1+alpha A_x+beta B_y+alpha beta E_xy)/4',
     'events':events}
(HERE/'C160_EVENT_GRAMS.json').write_text(json.dumps(out,separators=(',',':'))+'\n',encoding='utf-8')
print(json.dumps({'event_certificates':len(events),'ranks':[e['rank'] for e in events],
                  'coefficient_equalities':16*243,'normalization_polynomial_equalities':4,
                  'arithmetic':'exact rational','cache_sha256':out['cache_sha256']},indent=2))
