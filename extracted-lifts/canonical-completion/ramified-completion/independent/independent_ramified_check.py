"""Bounded independent ramified-pencil/compact-feature audit; standard library.

Imports no archive code. Reads frozen data, directly expands the homogeneous
SG2 parent and all forty feature products, and checks the two rational lifts.
Only independent_ramified_results.json beside this script is written.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import ast
import hashlib
import json

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('package', nargs='?', type=Path, default=HERE.parent)
PACKAGE = parser.parse_args().package.resolve()
DATA = PACKAGE/'data'
before = {str(p.relative_to(PACKAGE)): hashlib.sha256(p.read_bytes()).hexdigest()
          for p in DATA.rglob('*') if p.is_file()}
load = lambda name: json.loads((DATA/name).read_text())
pencil = load('ramified_order49_pencil.json')
assert (pencil['order'], pencil['auxiliary_dimension']) == (49, 148)
assert sum(map(len, pencil['source_matrices']+pencil['auxiliary_matrices'])) == 411

# Polynomial dictionary in a0,a1,rho,tau,W,X,Y,Z; no CAS or substitutions
# at numerical angle samples are used in the compact identity check.
NAMES = ['a0','a1','rho','tau','W','X','Y','Z']
ZERO = (0,)*8
one = {ZERO: Q(1)}
def add(*polys):
    result = {}
    for poly in polys:
        for exponent, coefficient in poly.items():
            result[exponent] = result.get(exponent, Q(0))+coefficient
    return {e:c for e,c in result.items() if c}
def scale(c, poly): return {e:c*x for e,x in poly.items() if c*x}
def mul(a, b):
    out = {}
    for e,x in a.items():
        for f,y in b.items():
            g = tuple(i+j for i,j in zip(e,f))
            out[g] = out.get(g, Q(0))+x*y
    return {e:c for e,c in out.items() if c}
def power(a, n):
    result = one
    for _ in range(n): result = mul(result, a)
    return result
var = {name: {tuple(int(i==j) for i in range(8)):Q(1)} for j,name in enumerate(NAMES)}
a0,a1,rho,tau,W,X,Y,Z = [var[name] for name in NAMES]
AliceD = add(power(a0,2), power(a1,2))

def parse_poly(text, epsilon):
    def visit(node):
        if isinstance(node,ast.Expression): return visit(node.body)
        if isinstance(node,ast.Name): return scale(epsilon,one) if node.id=='epsilon' else var[node.id]
        if isinstance(node,ast.Constant): return scale(Q(node.value),one)
        if isinstance(node,ast.UnaryOp):
            return scale(-1 if isinstance(node.op,ast.USub) else 1,visit(node.operand))
        if isinstance(node,ast.BinOp):
            if isinstance(node.op,ast.Pow):
                assert isinstance(node.right,ast.Constant) and isinstance(node.right.value,int)
                return power(visit(node.left),node.right.value)
            left,right=visit(node.left),visit(node.right)
            if isinstance(node.op,ast.Add): return add(left,right)
            if isinstance(node.op,ast.Sub): return add(left,scale(-1,right))
            if isinstance(node.op,ast.Mult): return mul(left,right)
        raise AssertionError(ast.dump(node))
    return visit(ast.parse(text.replace('^','**'),mode='eval'))

cover = load('compact_fixed_functions.json')
all_functions = [text for block in cover['feature_blocks'] for text in block['functions']]
assert len(all_functions)==40
polynomial_basis_checks=0
for epsilon in [-1,1]:
    pg = [add(mul(power(a0,2),mul(a1,W)),mul(a0,mul(power(a1,2),Y))),
          scale(epsilon,add(mul(power(a0,2),mul(a1,X)),mul(a0,mul(power(a1,2),Z)))),
          add(mul(power(a1,3),W),scale(-1,mul(a0,mul(power(a1,2),Y)))),
          scale(epsilon,add(mul(power(a1,3),X),scale(-1,mul(a0,mul(power(a1,2),Z)))))]
    ph = [add(mul(a0,mul(power(a1,2),W)),mul(power(a1,3),Y)),
          scale(epsilon,add(mul(a0,mul(power(a1,2),X)),mul(power(a1,3),Z))),
          add(mul(a0,mul(power(a1,2),W)),scale(-1,mul(power(a0,2),mul(a1,Y)))),
          scale(epsilon,add(mul(a0,mul(power(a1,2),X)),scale(-1,mul(power(a0,2),mul(a1,Z)))))]
    g = [power(tau,4),mul(power(rho,2),power(tau,2)),power(rho,4)]
    r = [mul(rho,power(tau,3)),mul(power(rho,3),tau)]
    features = [mul(bob,alice) for boblist,alicelist in [(g,pg),(r,pg),(g,ph),(r,ph)]
                for bob in boblist for alice in alicelist]
    for actual,text in zip(features,all_functions):
        assert text.endswith('/kappa')
        assert actual==parse_poly(text[:-len('/kappa')],epsilon)
    b0,b1=power(tau,4),scale(16*epsilon,power(rho,4))
    BobD=add(power(b0,2),power(b1,2))
    md=add(mul(a0,mul(b1,W)),scale(-1,mul(a0,mul(b0,X))),
           mul(a1,mul(b1,Y)),scale(-1,mul(a1,mul(b0,Z))))
    de=add(mul(a1,mul(b0,W)),mul(a1,mul(b1,X)),
           scale(-1,mul(a0,mul(b0,Y))),scale(-1,mul(a0,mul(b1,Z))))
    ae=add(mul(b1,Y),scale(-1,mul(b0,Z)))
    be=add(mul(a0,X),mul(a1,Z))
    source_forms = [
        scale(2,add(power(de,2),mul(AliceD,power(ae,2)),mul(BobD,power(be,2)),
                    scale(-1,mul(mul(AliceD,BobD),power(Z,2))))),
        scale(4,mul(AliceD,power(add(mul(b0,W),mul(b1,X)),2))),
        scale(4,mul(mul(AliceD,BobD),power(X,2))),
        scale(4,power(md,2)),
        scale(4,mul(BobD,power(add(mul(a1,W),scale(-1,mul(a0,Y))),2))),
        scale(4,mul(mul(AliceD,BobD),power(Y,2))),
        scale(4,power(de,2)),
        scale(4,mul(AliceD,power(ae,2))),
        scale(4,mul(BobD,power(be,2))),
    ]
    targets=[mul(mul(power(a1,2),AliceD),poly) for poly in source_forms]
    for index,sparse in enumerate(pencil['source_matrices']+pencil['auxiliary_matrices']):
        gram={}
        for i,j,c in sparse:
            if i<9 or j<9: continue
            gram=add(gram,scale(Q(c)*(1 if i==j else 2),mul(features[i-9],features[j-9])))
        assert gram==(targets[index] if index<9 else {}), (epsilon,index)
        polynomial_basis_checks+=1

# Independent exact rational reconstruction of both exported feasible lifts.
def matrix(a):return [[Q(x) for x in row] for row in a]
def transpose(a):return list(map(list,zip(*a)))
def mmul(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def msum(a,b):return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def congruence(v,x):return mmul(mmul(v,x),transpose(v))
def exact_ldl(a):
    n=len(a);l=[[Q(i==j) for j in range(n)] for i in range(n)];d=[]
    assert a==transpose(a)
    for j in range(n):
        pivot=a[j][j]-sum(l[j][k]**2*d[k] for k in range(j))
        assert pivot>Q(1,2**18)
        d.append(pivot)
        for i in range(j+1,n):l[i][j]=(a[i][j]-sum(l[i][k]*l[j][k]*d[k] for k in range(j)))/pivot
    assert congruence(l,[[d[i] if i==j else Q(0) for j in range(n)] for i in range(n)])==a
    return d
cert=load('ramified_exact.json')
full=[]
for index,(v,x) in enumerate(zip(cert['face_matrices'],cert['reduced_Gram_blocks'])):
    x=matrix(x);assert exact_ldl(x)==list(map(Q,cert['LDL_pivots'][index]))
    gram=congruence(matrix(v),x)
    diagonal=[Q(cert['D_G' if index<2 else 'D_H'][i][i]) for i in range(4)]
    assert all(diagonal)
    full.append([[value/diagonal[i%4]/diagonal[j%4] for j,value in enumerate(row)] for i,row in enumerate(gram)])
bar=load('ramified_cbar_exact.json');bargrams=list(map(matrix,bar['Gram_matrices']))
source=list(map(Q,cert['source']));sourcebar=list(map(Q,bar['source']))
public_path=PACKAGE.parents[1]/'outer-candidates/strict_outer_physical_source.json'
public_raw=public_path.read_bytes()
assert public_raw==(PACKAGE/'reused/ESTABLISHED_VALID_SOURCE.json').read_bytes()
old=json.loads(public_raw)
assert source==[Q(old['source'][0]),Q(0),Q(0)]+list(map(Q,old['source'][1:]))
assert sourcebar==[Q(1),Q(1,2**34),Q(1,2**14),Q(1,2**24),Q(1,576),Q(161),Q(1,2),Q(76),Q(1,250)]
assert min(sourcebar)>0
delta=[b-a for a,b in zip(source,sourcebar)]
assert delta[0]==0 and min(delta)>=0
for index,key in [(0,'positive_event_additions_to_Q_G'),(2,'positive_event_additions_to_Q_H')]:
    gram=full[index]
    for term in bar[key]:
        weight=Q(term['weight']);z=list(map(Q,term['vector']));assert weight>=0
        gram=msum(gram,[[weight*x*y for y in z] for x in z])
    assert gram==bargrams[index]
assert full[1]==bargrams[1] and full[3]==bargrams[3]
for filename,expected in [('ramified_order49_source_auxiliaries.json',full),('ramified_order49_cbar_auxiliaries.json',bargrams)]:
    values=load(filename);assert values['auxiliary_order']==pencil['auxiliary_names']
    parameters=list(map(Q,values['source']+values['auxiliary_values']))
    assert len(parameters)==157
    lifted=[[Q(0) for _ in range(49)] for _ in range(49)]
    for value,sparse in zip(parameters,pencil['source_matrices']+pencil['auxiliary_matrices']):
        for i,j,x in sparse:
            lifted[i][j]+=value*Q(x)
            if i!=j:lifted[j][i]+=value*Q(x)
    for i in range(9):assert lifted[i][i]==parameters[i]>=0
    for offset,size,gram in zip([9,21,29,41],[12,8,12,8],expected):
        assert [row[offset:offset+size] for row in lifted[offset:offset+size]]==gram
assert before=={str(p.relative_to(PACKAGE)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in DATA.rglob('*') if p.is_file()}
result={
 'scope':'Ramified construction, both exact rescues, and compact testing family only; no order33 separator or universal adequacy audit.',
 'pencil_order':49,'auxiliaries':148,'upper_triangular_records':411,
 'independent_compact_polynomial_basis_checks':polynomial_basis_checks,
 'identity_checked':'kappa^2*featureGram = a1^2*(a0^2+a1^2)*kappa^2*Fh, before imposing a0^2+a1^2=1',
 'all_40_literal_compact_functions_match':True,
 'both_epsilon_signs_checked':True,'four_reduced_Gram_LDLs_exact_above_2neg18':True,
 'cstar_matches_public_source_bytes':True,'cbar_matches_stated_all_positive_source':True,
 'both_148_auxiliary_vectors_reconstruct_actual_Grams':True,
 'cbar_Gram_additions_nonnegative_event_outer_products':True,
 'inputs_unchanged':True,'input_sha256':before,
}
(HERE/'independent_ramified_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({key:value for key,value in result.items() if key!='input_sha256'},indent=2))
