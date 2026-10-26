"""Independent bounded map/parity audit of the frozen Pro669 package."""
from pathlib import Path
from fractions import Fraction
import json
import sympy as sp

BASE=Path(__file__).resolve().parent/'pro669'/'Q222_Fixed_Family_Source_and_Matrices_2026-10-23'
read=lambda name:json.loads((BASE/name).read_text())
raw=read('c1_fixed_family/C1_RETAINED_ROOTS.json')
source_names='b k q sigma p m D beta C a ell h rho dG S0 nu R0 Rx Bx Ay chi lam q0 L0 KA c a0 kappa r'
source_locals={n:sp.Symbol(n,real=True) for n in source_names.split()}
parse=lambda s:sp.sympify(s,locals=source_locals)
stats={}
for sign,label in ((1,'POS'),(-1,'NEG')):
    fd=read(f'uniform_tf1/C1_{label}_FULL_FIXED_FEATURES.json')
    blocks=fd['blocks']
    for block in blocks:
        for e in block['features']:
            assert e[2]+e[3]+4*(e[4]+e[5])==4
        if block['category']!='principal':
            assert len({tuple(e[i]%2 for i in (0,1,2,3,6,7))
                        for e in block['features']})==1
    # Rebuild every retained fixed block directly from the source export.
    chart=raw['charts'][str(sign)]
    k=16
    for root in chart['roots']:
        fs=[];cs=[]
        for fi,c in root['entries']:
            e=chart['feature_exponents'][fi]
            fs.append((2*e[0],e[1],*e[2:],0,0))
            scale=Fraction(1,2**(6*e[0]+3*e[1]))
            cs.append(parse(c)*sp.Rational(scale.numerator,scale.denominator))
        shifts=[(0,)*8]
        if root['category']=='circuit_unweighted':
            shifts=[(0,)*8,(0,0,0,0,0,0,1,0),(2,0,0,0,0,0,1,0)]
        for shift in shifts:
            bb=blocks[k];k+=1
            assert bb['name']==root['name']
            assert bb['category']==root['category']
            assert bb['features']==[[x+y for x,y in zip(e,shift)] for e in fs]
            assert len(bb['source_vector'])==len(cs)
            assert all(sp.expand(parse(a)-b)==0
                       for a,b in zip(bb['source_vector'],cs))
    assert k==len(blocks)
    stats[label]={'blocks':len(blocks),'retained_raw_roots':len(chart['roots']),
       'retained_compiled_blocks':len(blocks)-16,
       'all_retained_maps_exact':True,'all_retained_cover_parities_uniform':True,
       'all_features_state_weight_one':True}

# Negative-state source congruence: compare every row, not just its target.
am=read('c1_fixed_family/C1_ROOT_MATRIX.json')
names=set()
for row in am['source_matrix_negative_state_sign']:
    for expr in row:names.update(str(s) for s in sp.sympify(expr).free_symbols)
loc={s:sp.Symbol(s,positive=True) for s in names}
AA=sp.Matrix([[sp.sympify(e,locals=loc) for e in row]
              for row in am['source_matrix_negative_state_sign']])
de=sp.Rational(1,4096)
alpha,bt,om,ze,X,Y=sp.symbols('alpha bt om ze X Y',real=True)
chi=sp.Matrix([X,Y,ze**4,om**4,alpha**8*ze**4,bt**4*ze**4,
    alpha**8*om**4,alpha**8*bt**4*om**4,alpha**16*om**4,
    alpha**12*om**4,alpha**4*om**4])
sc=sp.diag(1,1,-1,1,-de**2,-de,de**2,de**3,de**4,de**3,de)
got=AA[:,:11]*sc*chi
k,q,sig,h,b,dG,lam,rho,L0,q0=[loc[n] for n in
    ('k','q','sigma','h','b','dG','lam','rho','L0','q0')]
u=de**2*alpha**8;w=de*bt**4;W=om**4;Z=-ze**4
expected=sp.Matrix([X+u*h*Z,(Y-b*w*Z)/q,w*Z/q-u*q*W,
    sp.sqrt(dG)*(sig*Z-u*w*W),q*(u*u-k*k)*W,k*q*(u-k)*W,
    sp.sqrt(L0)*u*Z,sp.sqrt(rho)*u*w*W,sp.sqrt(lam)*k*k*q*W,
    sp.sqrt(lam)*k*w*Z/q,
    sp.sqrt(q0)*q*de*alpha**4*(u-k)*W,
    sp.sqrt(q0)*sp.sqrt(k)*q*(u-k)*W])
assert all(sp.expand(a-b)==0 for a,b in zip(got,expected))

# Exact circuit scalar split on the entire gate cover.
ac,qc0,A=sp.symbols('ac qc0 A')
xi2=1-A**4
assert sp.expand(qc0-ac*de**2*A**8-
    ((qc0-ac*de**2)+ac*de**2*(xi2+xi2*A**4)))==0
stats['negative_source_congruence_all_12_rows_exact']=True
stats['circuit_gate_split_exact']=True
dest=Path(__file__).with_suffix('.json')
dest.write_text(json.dumps(stats,indent=2)+'\n')
print(json.dumps(stats,indent=2))
