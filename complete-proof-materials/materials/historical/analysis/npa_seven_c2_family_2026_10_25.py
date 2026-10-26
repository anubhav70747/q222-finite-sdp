"""One exact native contact family, no sampling and no PSD solver."""
import json
from pathlib import Path
import sympy as s
import npa_seven_gram_elimination_2026_10_25 as base

c = s.symbols("c", positive=True)
z = c*c
u, v = s.symbols("u v", real=True)
I = s.eye(2)
Z = s.diag(1, -1)
X = s.Matrix([[0, 1], [1, 0]])
obs = [s.kronecker_product(Z, I),
       s.kronecker_product(((1-u*u)*Z+2*u*X)/(1+u*u), I),
       s.kronecker_product(I, Z),
       s.kronecker_product(I, ((1-v*v)*Z+2*v*X)/(1+v*v))]
ops = [s.eye(4)]+obs+[obs[0]*obs[2], obs[0]*obs[3],
                     obs[1]*obs[2], obs[1]*obs[3]]
weights = {
    base.d: 1, base.la: 0, base.be: 0,
    base.n: (2-5*z)/(8*(1+2*z)),
    base.m: (2+z)/(8*(1+2*z)),
    base.C: (2+z)/(2*(1+2*z)),
    base.A: 3*(2+z)/(8*(1+2*z)),
    base.D: 3*(2-z)/(8*(1+2*z)),
    base.al: 3*z*z/(2*(2-z)*(1+2*z)),
}
psi = s.Matrix([1, -c/2, c/2, 1])
F = s.zeros(4)
for coefficient, operator in zip(base.target, ops):
    F += coefficient.subs(weights)*operator
Fc = F.subs({u:c, v:c})
kernel = [s.factor(a) for a in Fc*psi]
jets = [s.factor((psi.T*s.diff(F,t).subs({u:c,v:c})*psi)[0]) for t in (u,v)]
assert kernel == [0]*4 and jets == [0]*2
G, free = base.calculate("C2", base.patterns["C2"], True)
assert G[2,3]==0 and G[5,6]==0
assert s.expand(G[0,2]-G[4,6]-G[6,6])==0
assert s.expand(G[1,2]+G[4,5])==0
G = G.subs(weights)
hops = [ops[a]-ops[b] for a,b in base.patterns["C2"]]
T = s.Matrix.hstack(*[h.subs({u:c,v:c})*psi for h in hops])
linear = [s.factor(a) for a in G*T.T]
# Short independent proof: four homogeneous universal-word relations reduce
# any contact-supported Gram to two squares, whose forbidden word coefficients
# then force both coefficients to vanish.
aa=3*(2-z)/(2*(1+z)); bb=-(z+4)/(2*(1+z))
dd=-bb; ee=-(z+2)/(2*(1+z))
N=s.Matrix([[0,aa,dd],[0,bb,ee],[-1,-aa,-ee],[1,0,0],
            [0,0,-1],[0,1,0],[0,0,1]])
assert all(s.factor(t)==0 for t in T*N)
assert s.factor(T[:,[0,1,2,4]].det()) == -8*c**4/(1+c*c)
assert s.factor(aa*ee-bb*dd)==1
s11,s12,s13,s22,s23,s33=s.symbols('s11 s12 s13 s22 s23 s33')
S=s.Matrix([[s11,s12,s13],[s12,s22,s23],[s13,s23,s33]])
Rc=N*S*N.T
relations=[Rc[2,3],Rc[5,6],Rc[0,2]-Rc[4,6]-Rc[6,6],Rc[1,2]+Rc[4,5]]
subs={s23:0,s12:-aa*s22,s13:-ee*s33,s11:aa*aa*s22+ee*ee*s33}
assert all(s.factor(t.subs(subs))==0 for t in relations)
constraint_matrix,_=s.linear_eq_to_matrix(relations,[s11,s12,s13,s23])
assert s.factor(constraint_matrix.det()) != 0
coefficient={}
for i in range(7):
    for j in range(7):
        for wi,ci in [(base.patterns['C2'][i][0],1),(base.patterns['C2'][i][1],-1)]:
            for wj,cj in [(base.patterns['C2'][j][0],1),(base.patterns['C2'][j][1],-1)]:
                word=base.multiply(base.words[wi],base.words[wj])
                coefficient[word]=coefficient.get(word,0)+ci*cj*Rc[i,j]
alice_word=s.factor(coefficient[((0,1),())].subs(subs))
bob_word=s.factor(coefficient[((),(0,1))].subs(subs))
assert s.factor(alice_word+2*aa*s22)==0
assert s.factor(bob_word+2*ee*s33)==0
data={
    'source_kernel_residuals':list(map(str,kernel)),
    'angular_jet_residuals':list(map(str,jets)),
    'contact_rank_minor':str(-8*c**4/(1+c*c)),
    'contact_null_basis':[[str(s.factor(t)) for t in N.col(i)] for i in range(3)],
    'ae_minus_bd':'1',
    'four_constraint_matrix_determinant':str(s.factor(constraint_matrix.det())),
    'forbidden_A0A1_coefficient':str(alice_word),
    'forbidden_B0B1_coefficient':str(bob_word),
    'conclusion':'s22=s33=0, hence entire Gram zero; no 1+AB certificate for 0<c^2<=2/5',
    'rational_example':{str(k):str(s.factor(v.subs(c,s.Rational(1,2)))) if hasattr(v,'subs') else str(v) for k,v in weights.items()},
}
Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
print(json.dumps(data,indent=2))
