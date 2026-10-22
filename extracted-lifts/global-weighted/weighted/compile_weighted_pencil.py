#!/usr/bin/env python3
"""Eliminate the evaluated product equations and export a rational pencil.

The pencil represents a cone of weighted polynomial coefficients.  Its pullback
under the exported quartic source map is exact for canonical valid slacks.
The quartic substitution is essential: this is NOT an affine Q222 lift.
"""
from pathlib import Path
import json
import sympy as s

HERE=Path(__file__).resolve().parent
F=json.loads((HERE/'WEIGHTED_FIXED_FUNCTIONS.json').read_text())
P=json.loads((HERE/'WEIGHTED_GRAM_PRODUCT_MAP.coo.json').read_text())
D=json.loads((HERE/'DETERMINANT_SOURCE_COEFFICIENTS.json').read_text())

def matrix(data):
    return s.SparseMatrix(*data['shape'],{(i,j):s.Rational(n,d) for i,j,n,d in data['entries']})

A=s.SparseMatrix(*P['matrix_shape'],{(i,j):c for i,j,c in P['entries']})
_,pivots=A.rref()
_,rows=A.T.rref()
pivots,rows=list(pivots),list(rows)
free=[j for j in range(A.cols) if j not in pivots]
rank=len(pivots)
M=A[rows,pivots].inv()
K=-M*A[rows,free]
assert A[:,pivots]*M*A[rows,:]==A

# x = X*(y,t) reconstructs every symmetric Gram coordinate.
X=s.zeros(A.cols,rank+len(free))
for i,p in enumerate(pivots):
    for j in range(rank):
        X[p,j]=M[i,j]
    for j in range(len(free)):
        X[p,rank+j]=K[i,j]
for j,f in enumerate(free):
    X[f,rank+j]=1
assert A[rows,:]*X==s.eye(rank).row_join(s.zeros(rank,len(free)))

offsets=[]
total=0
for order in P['psd_block_sizes']:
    offsets.append(total)
    total+=order
entries=[]
for i,coord in enumerate(P['gram_coordinates']):
    row=offsets[coord['block']]+coord['row']
    col=offsets[coord['block']]+coord['column']
    for j in range(X.cols):
        value=X[i,j]
        if value:
            entries.append([j,row,col,int(value.p),int(value.q)])

source=s.symbols(' '.join(D['variables']))
d,lam,n,m,C,Ac,Dc,alpha,beta=source
kappa=s.Matrix([4*lam,4*n,4*m,4*C,4*Ac,4*Dc+2*d,4*alpha+2*d,4*beta+2*d])
q={}
for rec in D['coefficients']:
    q[tuple(rec['ij'])]=sum(s.Rational(n,d)*s.prod(x**e for x,e in zip(source,p))
                             for p,n,d in rec['terms'])
target_grams=[]
for block in F['blocks']:
    T=matrix(block['original_to_block_basis'])
    if block['name']=='determinant':
        G=s.zeros(21)
        for j,nu in enumerate(block['n_indices']):
            if nu==[0,0]:
                G[j,j]=q[1,1]
            elif nu[0]%3==0 and nu[1]%3==0:
                G[j,j]=q[nu[0]//3+1,nu[1]//3+1]
        target_grams.append(T*G*T.T)
    else:
        _,i,j=block['name'].split('_')
        v=T*s.Matrix([kappa[k-1] for k in block['k_indices']])
        target_grams.append(kappa[int(i)-1]*kappa[int(j)-1]*(v*v.T))
baseline=s.Matrix([s.expand(target_grams[c['block']][c['row'],c['column']])
                   for c in P['gram_coordinates']])
visible=A[rows,:]*baseline
source_entries=[]
for i,expr in enumerate(visible):
    pp=s.Poly(expr,*source)
    for powers,value in pp.terms():
        if value:
            assert sum(powers)==4
            source_entries.append([i,list(powers),int(value.p),int(value.q)])

status='Exact rational pencil for the weighted Gram cone; canonical source input is quartic. NOT an affine lift of Q222.'
coefficient_reconstruction=A[:,pivots]*M
pencil={
    'schema':'Q222.weighted-coefficient-cone-rational-pencil.v1',
    'status':status,
    'matrix_order':total,
    'block_sizes':P['psd_block_sizes'],
    'block_offsets':offsets,
    'visible_coefficient_variables':rank,
    'free_gram_variables':len(free),
    'variable_order':[f'y{i:03d}' for i in range(rank)]+[f't{i:03d}' for i in range(len(free))],
    'constant_matrix':'zero',
    'entries':entries,
    'entry_format':'[variable_index, global_row, global_column, rational_numerator, positive_denominator]; upper triangle only, mirror off-diagonal entries.',
    'constraint':'L(y,t) is positive semidefinite; all 27 diagonal blocks are independent.',
    'selected_coefficient_rows':rows,
    'selected_monomials':[P['monomial_rows'][i] for i in rows],
    'pivot_gram_coordinates':pivots,
    'free_gram_coordinates':free,
    'coefficient_reconstruction':{
        'shape':[A.rows,rank],
        'entries':[[i,j,int(v.p),int(v.q)]
             for i in range(A.rows) for j in range(rank)
             if (v:=coefficient_reconstruction[i,j])]},
}
source_map={
    'schema':'Q222.canonical-source-quartic-map.v1',
    'status':status,
    'variables':[str(x) for x in source],
    'source_domain':'d,lambda,n,m,C,A,D,alpha,beta are all nonnegative.',
    'output_dimension':rank,
    'homogeneous_degree':4,
    'entries':source_entries,
    'entry_format':'[output_index, nine-source-variable exponent vector, rational_numerator, positive_denominator]',
    'equivalence':'The canonical full Bell slack is valid at all quantum tests iff there exist real t with L(y(source),t) positive semidefinite.',
    'warning':'Treating the quartic monomials as independent free variables loses this equivalence and does not produce a lift of Q222.',
}
receipt={
    'matrix_order':total,'psd_blocks':len(P['psd_block_sizes']),
    'independent_testing_functions':len(F['features']),
    'coefficient_rows_before_elimination':A.rows,
    'coefficient_rank':rank,'redundant_product_equations':A.rows-rank,
    'gram_coordinates':A.cols,'free_directions':len(free),
    'pencil_nonzero_upper_triangle_entries':len(entries),
    'quartic_source_map_terms':len(source_entries),
    'exact_elimination_verified':True,'old_order669_recomputed':False,
    'is_full_Q222_affine_lift':False,
}
for name,data in [('WEIGHTED_CONE_ORDER133.coo.json',pencil),
                  ('CANONICAL_SOURCE_QUARTIC_MAP.json',source_map),
                  ('pencil_compilation_receipt.json',receipt)]:
    (HERE/name).write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
