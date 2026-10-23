#!/usr/bin/env python3
"""Build ONLY the new order-eight pointwise selector pencil and counterexample.

Needs SymPy.  This does not load or regenerate any previously extracted pencil.
The resulting pencil is parameterized by Bob's test v and is NOT a finite lift.
"""
from pathlib import Path
import json
import hashlib
import itertools
import sympy as s

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1] / 'public_repo'
INPUT = HERE/'ESTABLISHED_VALID_SOURCE.json'
if not INPUT.exists():
    INPUT.write_bytes((REPO/'extracted-lifts/outer-candidates/strict_outer_physical_source.json').read_bytes())
v = s.symbols('v', real=True)
names = ['d', 'lambda', 'n', 'm', 'C', 'A', 'D', 'alpha', 'beta']
d, lam, n, m, C, A, D, al, be = s.symbols('d lambda n m C A D alpha beta')
source = [d, lam, n, m, C, A, D, al, be]
y = s.symbols('y11 y12 y21 y22')
p, r, t = s.symbols('p r s')
Dv = 1 + v*v
Ap = s.Matrix([[4*m*v*v, -4*m*v], [-4*m*v, 4*m+(4*be+2*d)*Dv]])
Am = s.Matrix([[4*C*Dv+4*D+2*d, (4*D+2*d)*v],
               [(4*D+2*d)*v, (4*D+2*d)*v*v]])
C0 = s.Matrix([[4*A*Dv+(4*al+2*d)*v*v, -(4*al+2*d)*v],
               [-(4*al+2*d)*v, 4*al-2*d*v*v]])
K = 4*lam*s.Matrix([[1,v],[v,v*v]])+4*n*Dv*s.diag(0,1)

def pencil(J):
    G = s.BlockMatrix([[Ap+K,J],[J.T,Am+K]]).as_explicit()
    H = s.BlockMatrix([[Ap+C0,K-C0-J],[K-C0-J.T,Am+C0]]).as_explicit()
    return G, H, s.diag(G,H)

def export(mat, variables, labels, filename, scope):
    # Full symmetric COO, zero-based.  [parameter, v_power, row, column, integer].
    entries=[]
    reconstructed=s.zeros(8)
    for k,x in enumerate(variables):
        co=mat.diff(x)
        for j in range(3):
            cj=co.applyfunc(lambda z:s.expand(z).coeff(v,j))
            for row,col in itertools.product(range(8), repeat=2):
                if cj[row,col]:
                    assert cj[row,col].is_Integer
                    entries.append([k,j,row,col,int(cj[row,col])])
                    reconstructed[row,col]+=x*v**j*cj[row,col]
    assert all(s.expand(z)==0 for z in mat-reconstructed)
    obj={
      'format':'q222-new-selector-polynomial-pencil-v1',
      'scope':scope,
      'order':8,'blocks':[4,4], 'parameter_order':labels,
      'nonnegative_parameter_indices':list(range(len(variables)-(4 if len(variables)==13 else 3))),
      'free_parameter_indices':list(range(len(variables)-(4 if len(variables)==13 else 3),len(variables))),
      'state_order_in_each_block':['a1','a2','b1','b2'],
      'testing_variable':'v', 'max_testing_degree':2,
      'indexing':'zero-based; full symmetric COO',
      'entry_order':['parameter_index','v_power','row','column','integer_coefficient'],
      'entries':entries,
      'equation':'L(v)=sum(parameter[k]*v**j*coefficient[k,j]); L(v)=diag(G(v),H(v))',
      'not_a_global_lift':True
    }
    (HERE/filename).write_text(json.dumps(obj,indent=2)+'\n')
    return len(entries)

G,H,full=pencil(s.Matrix(2,2,y))
full_count=export(full,source+list(y),names+[str(x) for x in y],
 'POINTWISE_SELECTOR_ORDER8.coo.json',
 'For each fixed Bob test v, the four selector entries are existential scalars. '
 'Canonical validity is forall v exists selector; the quantifier has NOT been eliminated.')
J=s.Matrix([[p*v*v,0],[r*v,t*v*v]])
G3,H3,c1=pencil(J)
zero={lam:0,n:0}
G3,H3,c1=(z.subs(zero) for z in (G3,H3,c1))
c1_source=[d,m,C,A,D,al,be]
c1_names=['d','m','C','A','D','alpha','beta']
c1_count=export(c1,c1_source+[p,r,t],c1_names+['p','r','s'],
 'C1_POLYNOMIAL_SELECTOR3_ORDER8.coo.json',
 'Derived exact form of any parity-symmetrized polynomial Bob selector on C1. '
 'There are three source-only selector parameters, but positivity for every v '
 'is NOT equivalent to canonical validity; the supplied valid source disproves that implication.')

old=json.loads(INPUT.read_text())
z={k:s.Rational(x) for k,x in zip(old['source_order'],old['source'])}
sub={d:z['d'],m:z['m'],C:z['C'],A:z['A'],D:z['D'],al:z['alpha'],be:z['beta'],lam:0,n:0}
u0=s.Rational(old['contact']['u']);v0=s.Rational(old['contact']['v']);U=u0*u0
W,X,Y,Z=map(s.Rational,old['contact']['state'])
Du=1+U
av=s.Matrix([(W+u0*Y)/Du,(X+u0*Z)/Du])
bv=s.Matrix([(W-Y/u0)/Du,(X-Z/u0)/Du])
Ap0=Ap.subs(sub).subs(v,v0);Am0=Am.subs(sub).subs(v,v0)
C00=C0.subs(sub).subs(v,v0)
PP=s.BlockMatrix([[(1+U)*Ap0+U*C00,-U*C00],
                    [-U*C00,U*(1+U)*Am0+U*C00]]).as_explicit()
assert PP*av.col_join(bv)==s.zeros(4,1)
equations=list(J.subs(v,v0)*bv+Ap0*av/U)+list(J.subs(v,v0).T*av+U*Am0*bv)
M, rhs=s.linear_eq_to_matrix(equations,[p,r,t])
assert M.rank()==3 and M.row_join(rhs).rank()==3
solution=s.solve(equations,[p,r,t])
Gbad=G3.subs(sub).subs(solution)
Hbad=H3.subs(sub).subs(solution)
minor=s.factor(Gbad.subs(v,1).extract([1,2,3],[1,2,3]).det())
assert minor<0
det=s.factor(Gbad.det())
target=v**4*(v-16)*(v+16)*(1949306657149583773*v*v-44177987181070080)
prefactor=s.cancel(det/target)
assert prefactor.is_Rational and prefactor>0
aa,bb,cc,ee=4*z['m'],4*z['beta']+2*z['d'],4*z['C'],4*z['D']+2*z['d']
assert 3<ee<4 and 0<cc<s.Rational(1,100)
assert 0<aa+2*bb<5 and solution[r]<s.Rational(5,4) and solution[t]>2
assert minor==2*ee*cc*(aa+2*bb)-ee*(solution[r]-solution[t])**2-2*cc*solution[t]**2

q=lambda x:str(x)
out={
 'status':'Exact obstruction to every polynomial Bob selector in the new Alice quadratic normal form.',
 'does_not_disprove':'A fixed algebraic testing family or an additional source-independent multiplier.',
 'established_input':{
  'repository_commit':'5fc0a0487bb94d7bdfa9929edb5f4f6f3ff50240',
  'path':'extracted-lifts/outer-candidates/strict_outer_physical_source.json',
  'sha256':hashlib.sha256(INPUT.read_bytes()).hexdigest(),
  'validity':'Accepted established input: independently verified order-13744 rational certificate. No old certificate rebuilt.'
 },
 'source_order':old['source_order'],'source':old['source'],
 'contact':old['contact'],
 'transformed_contact':{'a':list(map(q,av)), 'b':list(map(q,bv)), 'U':q(U)},
 'forced_selector':{'form':'[[p*v^2,0],[r*v,s*v^2]]','p':q(solution[p]),'r':q(solution[r]),'s':q(solution[t])},
 'kernel_equations':{
  'unknown_order':['p','r','s'],
  'coefficient_matrix':[[q(x) for x in row] for row in M.tolist()],
  'right_hand_side':list(map(q,rhs)), 'rank':3,
  'meaning':'J*b = -Aplus*a/U and J^T*a = -U*Aminus*b at the exact contact'
 },
 'negative_principal_minor':{'testing_v':'1','zero_based_rows':[1,2,3],'value':q(minor)},
 'determinant':{'positive_prefactor':q(prefactor),
  'remaining_factor':'v^4*(v-16)*(v+16)*(1949306657149583773*v^2-44177987181070080)'},
 'hand_bounds':{'a':'4m','b':'4beta+2d','c':'4C','e':'4D+2d',
  'bounds':['3<e<4','0<c<1/100','0<a+2b<5','r<5/4','s>2'],
  'minor':'2ec(a+2b)-e(r-s)^2-2cs^2 < 2/5-27/16 = -103/80 < 0'},
 'coo_full_nonzeros':full_count, 'coo_c1_nonzeros':c1_count
}
(HERE/'POLYNOMIAL_SELECTOR_COUNTEREXAMPLE.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS: new order-eight pointwise pencil,',full_count,'full COO entries')
print('PASS: derived C1 three-parameter pencil,',c1_count,'full COO entries')
print('PASS: exact contact kernel, rank-three selector system, and negative principal minor')
print('No established matrix archive was read or regenerated.')
