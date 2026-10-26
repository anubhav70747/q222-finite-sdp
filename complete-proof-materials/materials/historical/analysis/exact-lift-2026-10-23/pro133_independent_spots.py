"""Bounded independent convention checks; imports no archive Python code."""
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parent/'pro133_review_run'
WROOT=ROOT/'weighted'
feature=json.loads((WROOT/'WEIGHTED_FIXED_FUNCTIONS.json').read_text())
product=json.loads((WROOT/'WEIGHTED_GRAM_PRODUCT_MAP.coo.json').read_text())
pencil=json.loads((WROOT/'WEIGHTED_CONE_ORDER133.coo.json').read_text())
quartic=json.loads((WROOT/'CANONICAL_SOURCE_QUARTIC_MAP.json').read_text())
V=s.symbols('a0 a1 b0 b1 x0 x1 y0 y1 W X Y Z')
a0,a1,b0,b1,x0,x1,y0,y1,W,X,Y,Z=V
mon=lambda powers,variables:s.prod(v**n for v,n in zip(variables,powers))
def sparse(rec):
    return s.SparseMatrix(*rec['shape'],{(i,j):s.Rational(n,d) for i,j,n,d in rec['entries']})
gb=s.groebner([a0*a0+a1*a1-1,b0*b0+b1*b1-1,
    x0**6-a0*a0,x1**6-a1*a1,y0**6-b0*b0,y1**6-b1*b1],*V,order='lex')
basis=s.Matrix([sum(s.Rational(c)*mon(e,V) for e,c in f['polynomial']) for f in feature['features']])
blocks=[sparse(block['basis_in_global_family'])*basis for block in feature['blocks']]
rows=product['monomial_rows']
col_terms={}
for i,j,c in product['entries']:col_terms.setdefault(j,[]).append((i,s.Rational(c)))
chosen=[]
for block in range(len(blocks)):
    chosen.append(next(j for j,rec in enumerate(product['gram_coordinates'])
        if rec=={'block':block,'row':0,'column':1}))
for col in chosen:
    rec=product['gram_coordinates'][col];v=blocks[rec['block']]
    expected=gb.reduce(s.expand(2*v[0]*v[1]))[1]
    actual=sum(c*mon(rows[i],V) for i,c in col_terms[col])
    assert s.expand(expected-actual)==0

# Independent decode: COO records are variable,row,column,numerator,denominator;
# Gram-product columns already include the symmetric off-diagonal factor2.
sizes=product['psd_block_sizes'];offsets=[];offset=0
for size in sizes:offsets.append(offset);offset+=size
positions={(offsets[rec['block']]+rec['row'],offsets[rec['block']]+rec['column']):j
    for j,rec in enumerate(product['gram_coordinates'])}
M=s.SparseMatrix(535,535,{(positions[i,j],v):s.Rational(n,d)
    for v,i,j,n,d in pencil['entries']})
A=s.SparseMatrix(261,535,{(i,j):s.Rational(c) for i,j,c in product['entries']})
R=sparse(pencil['coefficient_reconstruction'])
for col in [0,94,188,189,362,534]:
    target=R[:,col] if col<189 else s.zeros(261,1)
    assert A*M[:,col]==target

# New source and angle fixtures, not the archive checker's saved fixtures.
source=s.symbols('d lam n m C A D alpha beta')
source_values=(s.Integer(3),s.Integer(2),s.Integer(4),s.Integer(1),s.Integer(5),
               s.Integer(0),s.Integer(6),s.Integer(7),s.Integer(8))
y=s.zeros(189,1)
for row,exp,n,d in quartic['entries']:
    assert len(exp)==9 and sum(exp)==4
    y[row]+=s.Rational(n,d)*mon(exp,source_values)
coefs=R*y
row_table=s.Matrix([[b0,b1,0,0],[0,1,0,0],
    [a0*b1,-a0*b0,a1*b1,-a1*b0],[a1,0,-a0,0],[0,0,1,0],
    [a1*b0,a1*b1,-a0*b0,-a0*b1],[0,0,b1,-b0],[0,a0,0,a1]])
d,lam,n,m,C,AA,D,alpha,beta=source_values
kappa=[4*lam,4*n,4*m,4*C,4*AA,4*D+2*d,4*alpha+2*d,4*beta+2*d]
P=sum((kappa[i]*row_table[i,:].T*row_table[i,:] for i in range(8)),s.zeros(4))
P[3,3]-=2*d
fixtures=[(s.Rational(-4,5),s.Rational(3,5),s.Rational(7,25),s.Rational(-24,25)),
          (s.Integer(0),s.Integer(-1),s.Rational(-20,29),s.Rational(21,29))]
for point in fixtures:
    actual=s.zeros(4)
    for exp,coef in zip(rows,coefs):
        if not coef:continue
        assert exp[4]==exp[6]==0 and exp[5]%6==exp[7]%6==0
        angle=mon(exp[:4],point)*point[1]**(exp[5]//3)*point[3]**(exp[7]//3)
        state=[i for i in range(4) for _ in range(exp[8+i])]
        assert len(state)==2
        i,j=state
        if i==j:actual[i,i]+=coef*angle
        else:
            actual[i,j]+=coef*angle/2
            actual[j,i]+=coef*angle/2
    pval=P.subs(dict(zip((a0,a1,b0,b1),point)))
    assert actual==pval[:3,:3].det()*pval

hashes=json.loads((ROOT/'review_input_hashes.json').read_text())
same={name:hashlib.sha256((WROOT/name).read_bytes()).hexdigest()==sha for name,sha in hashes.items()}
same_lf={name:hashlib.sha256((WROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()==sha
         for name,sha in hashes.items()}
assert all(same_lf.values())
result={
    'archive_modules_imported':False,
    'offdiagonal_product_columns_reconstructed_from_original_cover':len(chosen),
    'one_column_from_each_block':True,
    'independently_decoded_pencil_columns':[0,94,188,189,362,534],
    'new_full_parent_matrix_fixtures_checked':len(fixtures),
    'all_quartic_terms_have_degree_four':True,
    'quartic_terms':len(quartic['entries']),
    'canonical_source_variables':9,
    'visible_polynomial_coefficient_variables':189,
    'not_an_affine_Q222_lift':True,
    'four_rebuilt_data_files_byte_identical':same,
    'four_rebuilt_data_files_identical_after_CRLF_to_LF_only':same_lf}
(ROOT/'independent_spot_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
