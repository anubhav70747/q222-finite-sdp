"""Independent bounded exact checks for the pro133 analytical argument.

Does not read or certify the supplied order133 serialized pencil.
"""
from itertools import combinations
from functools import reduce
from math import gcd
from pathlib import Path
import json
import sympy as s

E=[(i,j) for i in (-1,0,1) for j in (-1,0,1) if (i,j)!=(0,0)]
N=set((3*i,3*j) for i,j in E)
N.update((a,b) for a in (-1,1) for b in (-1,1))
N.update((a,b) for a in (-1,1) for b in (-2,2))
N.update((a,b) for a in (-2,2) for b in (-1,1)); N.add((0,0))
assert len(N)==21
circuits=[]
for k in (2,3):
    for support in combinations(E,k):
        A=s.Matrix(2,k,lambda i,j:support[j][i])
        ker=A.nullspace()
        if len(ker)!=1: continue
        v=ker[0]
        if all(t<0 for t in v): v=-v
        if not all(t>0 for t in v): continue
        mult=s.ilcm(*[t.q for t in v]); z=[int(t*mult) for t in v]
        div=reduce(gcd,z); z=[t//div for t in z]
        circuits.append((support,z))
assert len(circuits)==12
assert sum(len(c)==2 for c,_ in circuits)==4
assert sum(w==[1,1,1] for _,w in circuits)==4
assert sum(sorted(w)==[1,1,2] for _,w in circuits)==4
for support,weights in circuits:
    if len(support)==2: continue
    if weights==[1,1,1]:
        for h in support:
            for l,m in combinations(support,2):
                a=tuple(h[i]+2*l[i] for i in (0,1))
                b=tuple(h[i]+2*m[i] for i in (0,1))
                assert a in N and b in N
                assert all((a[i]-b[i])%2==0 for i in (0,1))
    else:
        g=support[weights.index(2)]
        e,f=[p for p,w in zip(support,weights) if w==1]
        assert all((e[i]-f[i])%2==0 for i in (0,1))

# Equation14, after p_h=a_h^3 with product a_1*a_2*a_3=1.
a,b,c=s.symbols('a b c',positive=True)
nine=s.Rational(1,2)*(a+b+c)*((a-b)**2+(a-c)**2+(b-c)**2)
assert s.expand(nine-(a**3+b**3+c**3-3*a*b*c))==0

# General inverse-free adjugate completion.
h11,h12,h13,h22,h23,h33=s.symbols('h11 h12 h13 h22 h23 h33')
b1,b2,b3,cc=s.symbols('b1 b2 b3 cc')
W,X,Y,Z=s.symbols('W X Y Z')
H=s.Matrix([[h11,h12,h13],[h12,h22,h23],[h13,h23,h33]])
bb=s.Matrix([b1,b2,b3]); xx=s.Matrix([W,X,Y])
detH=H.det(); adjH=H.adjugate(); q=H*xx+bb*Z
FF=(xx.T*H*xx)[0]+2*Z*(bb.T*xx)[0]+cc*Z*Z
detP=cc*detH-(bb.T*adjH*bb)[0]
assert s.expand(detH*FF-(q.T*adjH*q)[0]-detP*Z*Z)==0

# The positive d-core is a full-rank3 Gram at every interior angle.
a0,a1,b0,b1,d=s.symbols('a0 a1 b0 b1 d',real=True)
core_rows=s.Matrix([[a1*b0,a1*b1,-a0*b0],[0,0,b1],[0,a0,0]])
assert s.expand(core_rows.det()**2-a0*a0*a1*a1*b0*b0*b1*b1)==0

# Explicit inverse entry in the exact-whitening example.
aa,v=s.symbols('aa v',positive=True)
HH=s.Matrix([[8*v*v+22,14*v,-14],[14*v,32*v*v+18,-14*v],
             [-14,-14*v,8*aa*(1+v*v)+20*v*v+14]])
rv=(320*v**4+621*v*v+126)/(2*(1+v*v)*(64*v*v+99))
invYY=(HH[:2,:2].det())/HH.det()
assert s.cancel(invYY-1/(8*(1+v*v)*(aa+rv)))==0
assert rv.subs(v,s.Rational(1,2))==s.Rational(241,230)
assert rv.subs(v,1)==s.Rational(1067,652)

out={'scope':'Bounded independent algebra/circuit checks; serialized pencil not audited.',
 'inverse_free_adjugate_identity':True,'core_frame_squared_minor':'a0^2*a1^2*b0^2*b1^2',
 'primitive_positive_circuit_count':len(circuits),
 'circuit_types':{'opposite_pairs':4,'triangles111':4,'triangles112':4},
 'angular_monomial_count':len(N),'circuit_support_and_parity':True,
 'triangle_nine_square_identity':True,'whitening_example_inverse_YY':True,
 'primitive_circuits':[{'support':c,'weights':w} for c,w in circuits]}
(Path(__file__).parent/'pro133_weighted_independent_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='primitive_circuits'},indent=2))
