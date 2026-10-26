"""Fixed symbolic a: exact quadratic identities over Q(a), no parameter scan."""
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix

directory=Path(__file__).parent
prior=json.loads((directory/'c2_maxentangled_general_angles_2026_10_25.json').read_text())
a,r,z=s.symbols('a r z')
local={'a':a,'r':r,'z':z}
source=[s.sympify(t,locals=local) for t in prior['clean_source']]
assert s.cancel(a*(source[0]+2*source[1])-2*source[4])==0
scale=s.sympify(prior['positive_common_polynomial_scale'],locals=local)
indices=[0,1,2,3,5,6]
names=['d','n','m','C','D','alpha']
polys=[s.Poly(s.cancel(source[i]*scale),r,z) for i in indices]
pairs=[(i,j) for i in range(6) for j in range(i,6)]
products=[polys[i]*polys[j] for i,j in pairs]
dicts=[p.as_dict() for p in products]
monomials=sorted(set().union(*(set(t) for t in dicts)))
M=s.Matrix([[p.get(mon,0) for p in dicts] for mon in monomials])
dm=DomainMatrix.from_Matrix(M).convert_to(s.QQ.poly_ring(a)).to_field()
red,pivots=dm.rref()
rr=red.to_Matrix()
null_rows=[]
for free in [j for j in range(len(pairs)) if j not in pivots]:
    row=[s.Integer(0)]*len(pairs);row[free]=s.Integer(1)
    for i,pivot in enumerate(pivots):row[pivot]=-rr[i,free]
    null_rows.append(row)
null=s.Matrix(null_rows)
print('Generic coefficient rank',len(pivots),'quadratic dimension',len(null_rows),flush=True)
coords=s.symbols('d n m C D alpha')
identities=[]
vectors=[]
for i in range(null.rows):
    vector=null.row(i).applyfunc(s.cancel)
    denom=s.lcm([s.denom(t) for t in vector])
    entries=[s.Poly(s.cancel(t*denom),a) for t in vector]
    common=s.gcd_list([t.as_expr() for t in entries])
    primitive=[s.cancel(t.as_expr()/common) for t in entries]
    polynomial=s.factor(sum(t*coords[j]*coords[k] for t,(j,k) in zip(primitive,pairs)))
    assert all(s.cancel(t)==0 for t in M*s.Matrix(primitive))
    identities.append(str(polynomial));vectors.append(list(map(str,primitive)))
out={'source_coordinates':names,'identity_eliminated':'2*A=a*(d+2*n)',
     'coefficient_matrix_shape':list(M.shape),'generic_rank':len(pivots),
     'quadratic_identity_dimension':null.rows,'primitive_identities':identities,
     'pairs':pairs,'identity_vectors':vectors}
(directory/'c2_fixed_a_quadrics_2026_10_26.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k not in ['pairs','identity_vectors']},indent=2),flush=True)
# One exact nonzero maximal minor bounds ALL exceptional parameter values.
rowpiv=dm.transpose().rref()[1]
square=M[list(rowpiv),list(pivots)]
print('Factoring one generic maximal minor...',flush=True)
maxminor=s.factor(DomainMatrix.from_Matrix(square).convert_to(s.QQ.poly_ring(a)).det().as_expr())
out['maximal_minor_rows']=list(rowpiv)
out['maximal_minor_columns']=list(pivots)
out['factored_maximal_minor']=str(maxminor)
(directory/'c2_fixed_a_quadrics_2026_10_26.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print('Maximal minor:',maxminor,flush=True)
# A factored basis of the generic identities, with A already eliminated.
dd,nn,mm,CC,DD,al=coords
AA=a*(dd+2*nn)/2
bb=AA+nn; PP=CC+dd/2; QQ=al+dd/2
gg=a*a*(bb+PP+QQ)-bb
factored=[mm*PP-DD*QQ,mm*gg-bb*QQ,DD*gg-bb*PP]
def coefficient_vector(expr):
    poly=s.Poly(s.expand(expr),*coords)
    result=[]
    for j,k in pairs:
        exponents=[0]*6;exponents[j]+=1;exponents[k]+=1
        result.append(poly.coeff_monomial(tuple(exponents)))
    return result
fv=s.Matrix([coefficient_vector(t) for t in factored])
assert DomainMatrix.from_Matrix(fv).to_field().rank()==3
assert all(s.cancel(t)==0 for t in M*fv.T)
exceptional={}
for value in [s.Integer(0),s.Integer(1),s.Integer(-1)]:
    specialized=DomainMatrix.from_Matrix(M.subs(a,value)).convert_to(s.QQ)
    exceptional[str(value)]={'rank':specialized.rank(),
                              'quadratic_identity_dimension':21-specialized.rank()}
# At a=1 an additional linear identity and two quadratic generators give
# every quadratic identity; no normalization by d is used here.
linear=mm+DD-(dd/2+2*nn)
special=[linear*t for t in coords]
special += [mm*CC-(dd/2+nn-mm)**2,
            DD*al-(dd/2+nn-DD)**2]
special_vectors=s.Matrix([coefficient_vector(t) for t in special])
assert special_vectors.rank()==8
assert M.subs(a,1)*special_vectors.T==s.zeros(M.rows,8)
assert exceptional['1']=={'rank':13,'quadratic_identity_dimension':8}
out.update({'generic_factored_basis':
    ['m*P-D*Q','m*Gamma-b*Q','D*Gamma-b*P'],
    'linear_forms':{'A':'a*(d+2*n)/2','b':'A+n','P':'C+d/2',
                    'Q':'alpha+d/2','Gamma':'a^2*(b+P+Q)-b'},
    'exceptional_formal_ranks':exceptional,
    'positive_a_exception':'a=1 only',
    'a1_linear_identity':str(linear),
    'a1_quadratic_basis':list(map(str,special)),
    'quadratically_invisible_direction':
       ['1','-a/(2*(a+1))','0','-1/2','0','-1/2']})
(directory/'c2_fixed_a_quadrics_2026_10_26.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print('Exceptional formal ranks',exceptional,flush=True)
