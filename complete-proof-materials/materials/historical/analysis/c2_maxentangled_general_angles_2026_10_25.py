"""One symbolic three-parameter stationary-source calculation, no sweep."""
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix

a,c,r,u,v=s.symbols('a c r u v',positive=True)
dc,nc,mc,Cc,Ac,Dc,alc=s.symbols('d n m C A D alpha')
active=[dc,nc,mc,Cc,Ac,Dc,alc]
I=s.eye(2); Z=s.diag(1,-1); X=s.Matrix([[0,1],[1,0]])
a0=s.kronecker_product(Z,I)
a1=s.kronecker_product(((1-u*u)*Z+2*u*X)/(1+u*u),I)
b0=s.kronecker_product(I,Z)
b1=s.kronecker_product(I,((1-v*v)*Z+2*v*X)/(1+v*v))
ops=[s.eye(4),a0,a1,b0,b1,a0*b0,a0*b1,a1*b0,a1*b1]
coeff=[dc+nc+mc+Cc+Ac+Dc+alc,nc-Ac-alc,mc-Cc-Dc,
       -nc+Cc+Ac,-mc+Dc-alc,-dc/2-nc-Ac,dc/2+alc,
       -dc/2-Cc,-dc/2-mc-Dc]
F=s.zeros(4)
for k,o in zip(coeff,ops):F+=k*o
psi=s.Matrix([1,-c*r,c*r,1]); sub={u:c,v:a*c}
equations=list(F.subs(sub)*psi)+[(psi.T*s.diff(F,t).subs(sub)*psi)[0] for t in (u,v)]
M,_=s.linear_eq_to_matrix([s.factor(k) for k in equations],active)
N=M[:,1:]; rhs=-M[:,0]
aug=N.row_join(rhs).applyfunc(s.cancel)
red,den,piv=DomainMatrix.from_Matrix(aug).to_field().rref_den()
R=red.to_Matrix()/den.as_expr()
assert piv==tuple(range(6))
answer=s.Matrix([1]+[s.factor(R[i,6]) for i in range(6)])
assert all(s.factor(k)==0 for k in M*answer)
z=s.symbols('z',positive=True)
source=[s.factor(k.subs(c*c,z)) for k in answer]
assert all(c not in k.free_symbols for k in source)
data={'active_source_order':list(map(str,active)),
      'source':list(map(str,source)),
      'all_six_contact_equations_verified':True}
Q0=1+(a*(1-r)+r)*z
L0=a+r-1+a*r*z
J0=1-a*r*z
R0=L0+r*Q0
W0=1+(a*a+a+2*r)*z+(a*a+a*r*(2-r)+r*r)*z*z
N0=(1-r)+(a-r*(a*a+3*a)+r*(2-3*r))*z-r*(a*a+a*r*(3-r)+r*r)*z*z
clean=[s.Integer(1),a*N0/(2*r*(a+1)*W0),
       (1-r)*(1+a*a*z)*(1+r*z)*Q0/(2*(a+1)*L0*W0),
       (1+r*z)*R0/(2*(1-r)*W0),
       a*(a+r)*(1+r*z)*Q0/(2*r*(a+1)*W0),
       (a+r)*(1+z)*J0*Q0/(2*(a+1)*L0*W0),
       a*z*z*(a+r)*R0/(2*J0*W0)]
assert all(s.factor(k-l)==0 for k,l in zip(source,clean))
assert s.factor(N0.subs(z,1/(a*r)))==-(a+1)**2*(a+r)/a**2
scale=2*r*(a+1)*(1-r)*L0*J0*W0
polynomials=[s.Poly(s.cancel(k*scale),a,r,z) for k in clean]
pairs=[(i,j) for i in range(7) for j in range(i,7)]
products=[polynomials[i]*polynomials[j] for i,j in pairs]
coefficients=[p.as_dict() for p in products]
monomials=sorted(set().union(*(set(k) for k in coefficients)))
prime=1000003
assert s.isprime(prime)
basis={}; selected=[]; selected_rows=[]
for monomial in monomials:
    row=[int(k.get(monomial,0))%prime for k in coefficients]
    original=row[:]
    for pivot,b in sorted(basis.items()):
        value=row[pivot]
        if value:row=[(x-value*y)%prime for x,y in zip(row,b)]
    pivot=next((i for i,k in enumerate(row) if k),None)
    if pivot is not None:
        inv=pow(row[pivot],-1,prime)
        basis[pivot]=[(x*inv)%prime for x in row]
        selected.append(monomial);selected_rows.append(original)
assert len(basis)==27
columns=sorted(basis)
square=[[row[j] for j in columns] for row in selected_rows]
detmod=1
for j in range(27):
    pivot=next(i for i in range(j,27) if square[i][j])
    if pivot!=j:
        square[j],square[pivot]=square[pivot],square[j]
        detmod=-detmod
    value=square[j][j]
    detmod=detmod*value%prime
    inv=pow(value,-1,prime)
    for i in range(j+1,27):
        ratio=square[i][j]*inv%prime
        square[i]=[(x-ratio*y)%prime for x,y in zip(square[i],square[j])]
assert detmod!=0
relation=[0]*28
for pair,value in [((0,2),-1),((0,5),1),((2,3),-2),((5,6),2)]:
    relation[pairs.index(pair)]=value
assert sum((value*p for value,p in zip(relation,products)),s.Poly(0,a,r,z)).is_zero
relation_polynomial=dc*(Dc-mc)-2*mc*Cc+2*Dc*alc
assert (s.hessian(relation_polynomial,active)/2).rank()==4
old_gap=s.factor(clean[5]*clean[6]-(clean[4]-clean[5])**2)
old_gap_limit=s.limit(old_gap.subs({a:2,r:s.Rational(1,2)}),z,0,dir='+')
assert old_gap_limit==-s.Rational(625,324)
minor=s.factor(DomainMatrix.from_Matrix(N.applyfunc(s.cancel)).det().as_expr())
data.update({
    'positive_helpers':{k:str(v) for k,v in [('Q',Q0),('L',L0),('J',J0),('R',R0),('W',W0),('N',N0)]},
    'clean_source':list(map(str,clean)),
    'exact_native_domain':'a>0,0<r<1,z>0,L>0,N>=0; N>=0 implies J>0',
    'N_at_J_zero':str(s.factor(N0.subs(z,1/(a*r)))),
    'positive_common_polynomial_scale':str(scale),
    'scaled_polynomial_term_counts':[len(p.terms()) for p in polynomials],
    'scaled_polynomial_degrees':[[int(p.degree(t)) for t in (a,r,z)] for p in polynomials],
    'quadratic_source_monomial_pairs':pairs,
    'unique_quadratic_relation_coefficients':relation,
    'quadratic_relation':'d*(D-m)-2*m*C+2*D*alpha=0',
    'quadratic_form_rank':4,
    'coefficient_matrix_rows':len(monomials),
    'coefficient_matrix_rank_lower_bound_mod_prime':27,
    'prime':prime,
    'independent_coefficient_row_monomials':selected,
    'independent_coefficient_columns':columns,
    'nonzero_minor_mod_prime':detmod,
    'coefficient_rows_mod_prime':selected_rows,
    'contact_minor_delete_d':str(minor),
    'old_D_alpha_block_det_limit_a2_r_half':str(old_gap_limit),
})
Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
print(json.dumps({key:data[key] for key in
    ['all_six_contact_equations_verified','exact_native_domain','quadratic_relation',
     'quadratic_form_rank','coefficient_matrix_rows','coefficient_matrix_rank_lower_bound_mod_prime',
     'prime','nonzero_minor_mod_prime','contact_minor_delete_d',
     'old_D_alpha_block_det_limit_a2_r_half']},indent=2))
