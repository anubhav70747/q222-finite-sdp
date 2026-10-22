"""Exact repaired C1 ray and a rational open family; no SDP solver.

Run with Python >=3.11 and SymPy from this directory or the project root.
The input JSON files are copied byte for byte from the historical scientific
construction; the current N13744 canonical source pencil is checked separately.
"""
from __future__ import annotations
from collections import defaultdict
from fractions import Fraction as Q
import hashlib, json, pathlib, sys, time, zipfile
import sympy as sp

sys.set_int_max_str_digits(100000)
HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT = HERE / 'c1_contact_selector_data'
OUT.mkdir(exist_ok=True)
ARCHIVE = ROOT / 'source/recovered/Q222_Extraction_Continuation_Package_2026-10-22.zip'
CURRENT = ROOT / 'analysis/human-organization-2026-10-23/extracted-lifts/outer-candidates'
PACKAGED_PENCIL = HERE.parent / 'outer-candidates/C1_enlarged_certificate_source_pencil.coo.json'
CURRENT_PENCIL = (PACKAGED_PENCIL if PACKAGED_PENCIL.exists()
                 else CURRENT/'C1_enlarged_certificate_source_pencil.coo.json')
NAMES = ['C1_REPAIRED_TWO_FRAME_COUNTEREXAMPLE.json',
         'C1_ENLARGED_REGULAR_CONTACT_SYSTEM.json',
         'C1_REGULAR_CONTACT_DOMAIN.json',
         'C1_REPAIRED_CHAMBER_AUXILIARY_MAP.json']
checks = {}
start = time.monotonic()
if ARCHIVE.exists():
    with zipfile.ZipFile(ARCHIVE) as archive:
        for name in NAMES:
            matches = [n for n in archive.namelist() if n.endswith('/'+name)]
            assert len(matches) == 1
            (OUT/name).write_bytes(archive.read(matches[0]))
        raw = archive.read(next(n for n in archive.namelist()
                           if n.endswith('/inputs/C1_ENLARGED_PRIMAL.coo.json')))
        (OUT/'C1_ENLARGED_PRIMAL.coo.json').write_bytes(raw)
else:
    raw=(OUT/'C1_ENLARGED_PRIMAL.coo.json').read_bytes()
if CURRENT_PENCIL.exists():
    assert raw == CURRENT_PENCIL.read_bytes()
assert hashlib.sha256(raw).hexdigest()=='1e1eb340077940be9c466f176efc4816d95c5271d448b7b52e94aa6a6542e344'
checks['current_N13744_canonical_source_pencil_SHA256'] = hashlib.sha256(raw).hexdigest()
literal = json.loads(raw)
seed, system, domain, chamber = [json.loads((OUT/n).read_text()) for n in NAMES]
c = list(map(Q, seed['coefficients']))
y = list(map(Q, seed['auxiliary_values']))
theta = list(map(Q, seed['contact_parameters']))
p,q,a,b,t = theta
phi = [a,b,t,Q(1)]
features = system['features']
weights = system['weights']
orders = list(map(len, features))

def zeros(n): return [[Q(0) for j in range(n)] for i in range(n)]
def evaluate_pencil(c,y):
    v = [Q(1)] + c + y
    M = zeros(literal['matrix_order'])
    for k,i,j,value in literal['entries']:
        M[i][j] += v[k]*Q(value)
        if i != j: M[j][i] += v[k]*Q(value)
    return [ [row[bl['offset']:bl['offset']+bl['order']]
               for row in M[bl['offset']:bl['offset']+bl['order']]]
             for bl in literal['blocks'] ]

G = evaluate_pencil(c,y)
assert G[-1] == [[Q(1)]]
for g in range(4):
    assert G[g] == [[Q(x) for x in row] for row in seed['grams'][g]]
    factor = seed['rational_PSD_factorizations'][g]
    F = [[Q(x) for x in row] for row in factor['factor']]
    D = list(map(Q, factor['positive_diagonal']))
    assert len(D) == orders[g]-1 and min(D)>0
    assert all(G[g][i][j] == sum(F[i][k]*D[k]*F[j][k]
               for k in range(len(D))) for i in range(orders[g]) for j in range(orders[g]))
    assert all(F[i][j] == 0 for i in range(len(D)) for j in range(i+1,len(D)))
    assert all(F[i][i] != 0 for i in range(len(D)))
    contact = [phi[i]*p**r*q**s for i,r,s in features[g]]
    assert all(sum(G[g][i][j]*contact[j] for j in range(orders[g]))==0
               for i in range(orders[g]))
checks['four_actual_pencil_Grams_and_positive_LDL_ranks'] = [n-1 for n in orders]

# Scalar polynomials in u,v,W,X,Y,Z are indexed by (r,s,i,j), i<=j.
def add_square(dst, form, multiplier, du=False, dv=False, uvweight=(0,0)):
    for i,r,s,x in form:
        for j,R,S,z in form:
            ii,jj = sorted((i,j))
            for dr in ([0,2] if du else [0]):
                for ds in ([0,2] if dv else [0]):
                    dst[(r+R+dr+uvweight[0],s+S+ds+uvweight[1],ii,jj)] += multiplier*x*z

def physical_polynomial(source):
    d,m,C,A,D,al,be = source
    z=defaultdict(Q)
    add_square(z,[(0,0,1,1),(1,0,0,-1),(2,1,1,1),(3,1,0,-1)],4*m)
    add_square(z,[(0,1,0,1),(2,0,0,-1)],4*C,dv=True)
    add_square(z,[(2,0,0,1)],4*A,du=True,dv=True)
    add_square(z,[(0,1,0,1),(1,1,1,1),(2,0,0,-1),(3,0,1,-1)],4*D+2*d)
    add_square(z,[(2,0,1,1),(3,0,0,-1)],4*al+2*d,du=True)
    add_square(z,[(1,0,0,1),(3,1,0,1)],4*be+2*d,dv=True)
    add_square(z,[(3,0,0,1)],-2*d,du=True,dv=True)
    return {k:v for k,v in z.items() if v}

def multiplied_physical(source):
    z=defaultdict(Q)
    for (r,s,i,j),value in physical_polynomial(source).items():
        for dr in [0,1]:
            for ds in [0,1]: z[r+dr,s+ds,i,j]+=value
    return {k:v for k,v in z.items() if v}

monomials = [tuple(x) for x in system['monomials']]
rowid = {x:i for i,x in enumerate(monomials)}
E=defaultdict(Q)
for r,k,v in system['unscaled_physical_target_map']: E[r,k]+=Q(v)
for k in range(7):
    src=[Q(int(i==k)) for i in range(7)]
    target=multiplied_physical(src)
    assert all(key in rowid for key in target)
    assert [target.get(mon,Q(0)) for mon in monomials] == [E[i,k] for i in range(148)]
checks['full_parent_source_map_all_seven_columns'] = True

def gram_polynomial(grams):
    z=defaultdict(Q)
    for g,features_g in enumerate(features):
        for k,(i,r,s) in enumerate(features_g):
            for l,(j,R,S) in enumerate(features_g):
                ii,jj=sorted((i,j)); al,be=weights[g]
                z[r+R+al,s+S+be,ii,jj]+=grams[g][k][l]
    return {k:v for k,v in z.items() if v}
assert gram_polynomial(G[:4]) == multiplied_physical(c)
checks['seed_full_polynomial_identity'] = True

# Six independent contact equations, including the full four-state kernel.
phys=physical_polynomial(c)
state=[Q(0)]*4; angular=[Q(0)]*2
for (r,s,i,j),value in phys.items():
    term=value*p**r*q**s
    state[i]+=term*phi[j]
    state[j]+=term*phi[i]
    angular[0]+=value*r*p**(r-1)*q**s*phi[i]*phi[j]
    angular[1]+=value*s*p**r*q**(s-1)*phi[i]*phi[j]
assert state == [0]*4 and angular == [0]*2
checks['seed_full_contact_four_state_two_angle_equations'] = True

# Verify the contact coefficient map from the actual difference features.
R=defaultdict(Q)
variables=[tuple(x) for x in system['gram_variables']]
for k,(g,i,j) in enumerate(variables):
    fi,fl=features[g][i],features[g][-1]
    fj=features[g][j]
    for (state_i,r,s),si in [(fi,1),(fl,-1)]:
        for (state_j,Rexp,Sexp),sj in [(fj,1),(fl,-1)]:
            ii,jj=sorted((state_i,state_j)); al,be=weights[g]
            R[rowid[r+Rexp+al,s+Sexp+be,ii,jj],k] += si*sj*(1 if i==j else 2)
R={k:v for k,v in R.items() if v}
assert R == {(i,j):Q(v) for i,j,v in system['coefficient_matrix']}

C=sp.zeros(6,148)
for col,(r,s,i,j) in enumerate(monomials):
    C[i,col]+=1; C[j,col]+=1; C[4,col]=r; C[5,col]=s
assert C.rank()==6
Rcols=defaultdict(list)
for (i,j),v in R.items(): Rcols[j].append((i,v))
assert all(sum(Q(C[k,i])*v for i,v in Rcols[j])==0
           for k in range(6) for j in range(392))
I=system['independent_rows']; J=system['pivot_columns']; free=system['free_columns']
P=defaultdict(Q)
for i,j,v in system['pivot_rhs_map']: P[i,j]+=Q(v)
Rij=defaultdict(list)
Iinv={r:i for i,r in enumerate(I)}
for k,col in enumerate(J):
    for r,v in Rcols[col]:
        if r in Iinv: Rij[Iinv[r]].append((k,v))
Prows=defaultdict(list)
for (i,j),v in P.items(): Prows[i].append((j,v))
for i in range(142):
    product=defaultdict(Q)
    for k,v in Rij[i]:
        for j,w in Prows[k]: product[j]+=v*w
    assert {j:v for j,v in product.items() if v} == {i:Q(1)}
checks['contact_map_derived_rank_and_full_compatibility'] = {'rank':142,'constraints':6}

normalized=[]
H=[]
for g in range(4):
    al,be=weights[g]
    diag=[p**r*q**s*phi[i] for i,r,s in features[g]]
    Gt=[[p**al*q**be*diag[i]*G[g][i][j]*diag[j]
         for j in range(orders[g])] for i in range(orders[g])]
    assert all(sum(row)==0 for row in Gt)
    normalized.append(Gt)
    H.append([row[:-1] for row in Gt[:-1]])
h=[H[g][i][j] for g,i,j in variables]
b0=[p**r*q**s*phi[i]*phi[j]*sum(E[k,l]*c[l] for l in range(7))
    for k,(r,s,i,j) in enumerate(monomials)]
assert all(sum(v*h[j] for (r,j),v in R.items() if r==i)==b0[i] for i in range(148))

# Rational strict lower bound from exact inverse traces, without large radicals.
def upper_binary_exponent(x):
    x=abs(x)
    if not x: return 0
    e=x.numerator.bit_length()-x.denominator.bit_length()+1
    while x>Q(2)**e: e+=1
    while x<=Q(2)**(e-1): e-=1
    return e
exponents=[]; inverse_traces=[]
for Hg in H:
    M=sp.Matrix([[sp.Rational(x.numerator,x.denominator) for x in row] for row in Hg])
    inv=M.inv(method='DM')
    trace=Q(str(sp.trace(inv)))
    assert trace>0
    inverse_traces.append(str(trace))
    exponents.append(max(0,upper_binary_exponent(trace)))
kappa=max(exponents)
mu=Q(2)**(-kappa)

# If max|b(theta)-b0|<mu/(2L), all four compressed matrices stay positive.
rowbounds=[[Q(0)]*(n-1) for n in orders]
for pivot_index,col in enumerate(J):
    g,i,j=variables[col]
    bound=sum(abs(v) for k,v in Prows[pivot_index])
    rowbounds[g][i]+=bound
    if i!=j: rowbounds[g][j]+=bound
L=max(x for row in rowbounds for x in row)
epsilon=mu/(2*L)
checks['strict_Gram_open_family']={'compressed_orders':[n-1 for n in orders],
    'inverse_trace_upper_binary_exponents':exponents,'kappa':kappa,
    'mu':str(mu),'L':str(L),'block_row_norm_bounds':[str(max(r)) for r in rowbounds],
    'epsilon':str(epsilon)}

# The complete quarter-parent contact matrix in the data is rederived.
u,v,aa,bb,tt=sp.symbols('u v a b t'); ph=sp.Matrix([aa,bb,tt,1])
Jcontact=sp.zeros(6,7)
for col in range(7):
    unit=[Q(int(i==col)) for i in range(7)]
    M=sp.zeros(4,4)
    for (r,s,i,j),value in physical_polynomial(unit).items():
        value=sp.Rational(value.numerator,value.denominator)*u**r*v**s/4
        if i==j: M[i,j]+=value
        else: M[i,j]+=value/2; M[j,i]+=value/2
    ker=M*ph
    for i in range(4): Jcontact[i,col]=sp.expand(ker[i])
    Jcontact[4,col]=sp.expand((ph.T*M.diff(u)*ph)[0])
    Jcontact[5,col]=sp.expand((ph.T*M.diff(v)*ph)[0])
for i in range(6):
    for j in range(7): assert sp.expand(Jcontact[i,j]-sp.sympify(domain['contact_matrix'][i][j]))==0
Jseed=Jcontact.subs({u:sp.Rational(str(p)),v:sp.Rational(str(q)),
                    aa:sp.Rational(str(a)),bb:sp.Rational(str(b)),tt:sp.Rational(str(t))})
assert Jseed*sp.Matrix([sp.Rational(str(x)) for x in c])==sp.zeros(6,1)
minor=sp.factor(Jseed[:,1:].det(method='domain-ge'))
assert minor!=0
checks['contact_source_chart_all_42_entries_and_seed_minor']=str(minor)
derivatives=sp.zeros(6,5)
shape_vars=[u,v,aa,bb,tt]
seed_sub={u:sp.Rational(str(p)),v:sp.Rational(str(q)),
          aa:sp.Rational(str(a)),bb:sp.Rational(str(b)),tt:sp.Rational(str(t))}
csym=sp.Matrix([sp.Rational(str(x)) for x in c])
for k,var in enumerate(shape_vars):
    derivatives[:,k]=-Jseed[:,1:].inv(method='DM')*(Jcontact.diff(var)*csym).subs(seed_sub)
checks['source_chart_Jacobian_rank_at_seed']=derivatives.rank()
assert derivatives.rank()==5
minor_rows=list(derivatives.T.rref()[1])
checks['source_chart_Jacobian_rank_certificate']={
    'rows_zero_based':minor_rows,
    'five_by_five_minor':str(derivatives[minor_rows,:].det(method='domain-ge'))}

# A short one-parameter subfamily, retaining exact contact and all signs.
family_J=Jcontact.subs({u:sp.Rational(str(p)),v:sp.Rational(str(q)),
                      bb:sp.Rational(str(b)),tt:sp.Rational(str(t))})
family_c=sp.Matrix([1]).col_join(-family_J[:,1:].inv(method='DM')*family_J[:,0])
family_c=[sp.factor(x) for x in family_c]
family_minor=sp.factor(family_J[:,1:].det(method='domain-ge'))
expected_minor=sp.Rational(6165476003,344646229622500000000)*(25*aa-5219)*(50*aa-13)*(166400*aa-219)*(83200*aa**2+17476053*aa+119)
assert sp.expand(family_minor-expected_minor)==0
shift=sp.Symbol('s')
def coeff_bounds(poly):
    pp=sp.Poly(poly,shift)
    if pp.is_zero: return sp.S.Zero,sp.S.Zero,sp.S.Zero,sp.S.Zero
    const=pp.nth(0)
    tail=sum(abs(pp.nth(k)) for k in range(1,pp.degree()+1))
    norm=abs(const)+tail
    deriv=sum(k*abs(pp.nth(k)) for k in range(1,pp.degree()+1))
    return const,tail,norm,deriv
family_sign_margins=[]
for value in family_c[1:]:
    N,D=map(sp.expand,sp.fraction(value.subs(aa,140+shift)))
    N0,Ntail,_,_=coeff_bounds(N); D0,Dtail,_,_=coeff_bounds(D)
    assert abs(N0)>Ntail and abs(D0)>Dtail and N0*D0>0
    family_sign_margins.append([str(abs(N0)-Ntail),str(abs(D0)-Dtail)])
family_phi=[aa,sp.Rational(str(b)),sp.Rational(str(t)),1]
family_target=[]
derivative_bound=Q(0)
for k,(r,s,i,j) in enumerate(monomials):
    value=sp.cancel(sp.Rational(str(p**r*q**s))*family_phi[i]*family_phi[j]
             *sum(sp.Rational(str(E[k,l]))*family_c[l] for l in range(7)))
    family_target.append(value)
    N,D=map(sp.expand,sp.fraction(value.subs(aa,140+shift)))
    N0,Ntail,Nnorm,Nder=coeff_bounds(N); D0,Dtail,Dnorm,Dder=coeff_bounds(D)
    Dlow=abs(D0)-Dtail
    assert Dlow>0
    derivative_bound=max(derivative_bound,Q(str((Nder*Dnorm+Nnorm*Dder)/Dlow**2)))
radius_exponent=max(1,upper_binary_exponent(2*max(Q(1),derivative_bound)/epsilon))
radius=Q(2)**(-radius_exponent)
assert derivative_bound*radius<=epsilon/2

def first_nonpositive_pivot(M):
    M=[row[:] for row in M]
    for i in range(len(M)):
        pivot=M[i][i]
        if pivot<=0: return i,pivot
        for j in range(i+1,len(M)):
            for k in range(j,len(M)):
                M[j][k]-=M[j][i]*M[i][k]/pivot; M[k][j]=M[j][k]
    return None
family_far_checks=[]
for far_shape in [141,70]:
    far_source=[Q(str(x.subs(aa,far_shape))) for x in family_c]
    assert min(far_source)>0
    assert family_J[:,1:].subs(aa,far_shape).det(method='domain-ge')!=0
    far_target=[Q(str(x.subs(aa,far_shape))) for x in family_target]
    far_h=h[:]
    for i,col in enumerate(J):
        far_h[col]+=sum(v*(far_target[I[k]]-b0[I[k]]) for k,v in Prows[i])
    far_H=[zeros(n-1) for n in orders]
    for value,(g,i,j) in zip(far_h,variables): far_H[g][i][j]=far_H[g][j][i]=value
    failures=[]
    for g,M in enumerate(far_H):
        bad=first_nonpositive_pivot(M)
        if bad is not None: failures.append([g,bad[0],str(bad[1])])
    family_far_checks.append({'state_a':far_shape,'source_all_positive':True,
                              'nonpositive_compressed_LDL_pivots':failures})
checks['one_parameter_family']={
 'shape':'(1/13,16,A,-1/32,1/50)',
 'six_positive_source_coefficients_on':'139<=A<=141',
 'certified_selector_interval':f'abs(A-140)<=2^(-{radius_exponent})',
 'target_derivative_bound':str(derivative_bound),
 'larger_parameter_selector_tests':family_far_checks}
family_data={'scope':'Exact rational C1 stationary source family; the explicit constant-free selector is proved on the displayed small interval only.',
 'shape':'(1/13,16,A,-1/32,1/50)',
 'source_order':['d','m','C','A','D','alpha','beta'],
 'source_formulas':[str(x) for x in family_c],
 'positive_source_interval':'139<=A<=141',
 'contact_event_minor':str(family_minor),
 'numerator_denominator_centered_sign_margins':family_sign_margins,
 'target_derivative_upper_bound':str(derivative_bound),
 'certified_radius_power_of_two':radius_exponent,
 'larger_parameter_selector_tests':family_far_checks}
(OUT/'ONE_PARAMETER_SOURCE_FAMILY.json').write_text(json.dumps(family_data,indent=2)+'\n')

# A visible-coordinate inequality, with no convention inferred from status flags.
event_specs=[(1,-1,1,1),(-1,1,1,0),(-1,1,0,0),
             (-1,1,1,1),(-1,-1,0,1),(1,-1,1,0)]
visible=[Q(1),Q(0),Q(0),Q(0),Q(0),Q(-1,2),Q(1,2),Q(-1,2),Q(-1,2)]
for weight,(sa,sb,x,yy) in zip(c[1:],event_specs):
    visible[0]+=weight; visible[1+x]+=sa*weight
    visible[3+yy]+=sb*weight; visible[5+2*x+yy]+=sa*sb*weight
checks['seed_inequality_order']=['constant','a0','a1','b0','b1','E00','E01','E10','E11']
checks['seed_inequality_coefficients']=[str(x) for x in visible]

# Restore the Bell operator in the original observable coordinates.
Id=sp.eye(2); Zpauli=sp.diag(1,-1); Xpauli=sp.Matrix([[0,1],[1,0]])
Du,Dv=1+u*u,1+v*v
Anum=[Du*Zpauli,(1-u*u)*Zpauli+2*u*Xpauli]
Bnum=[Dv*Zpauli,(1-v*v)*Zpauli+2*v*Xpauli]
core=Du*Dv*sp.eye(4)-(sp.kronecker_product(Anum[0],Bnum[0])
       -sp.kronecker_product(Anum[0],Bnum[1])
       +sp.kronecker_product(Anum[1],Bnum[0])
       +sp.kronecker_product(Anum[1],Bnum[1]))/2
bell_columns=[core]+[sp.kronecker_product(Du*Id+sa*Anum[x],Dv*Id+sb*Bnum[yy])
                       for sa,sb,x,yy in event_specs]
for col in range(7):
    unit=[Q(int(i==col)) for i in range(7)]; direct=sp.zeros(4,4)
    for (r,s,i,j),value in physical_polynomial(unit).items():
        vv=sp.Rational(str(value))*u**r*v**s
        if i==j: direct[i,j]+=vv
        else: direct[i,j]+=vv/2; direct[j,i]+=vv/2
    assert all(sp.expand(x)==0 for x in direct-bell_columns[col])
checks['Bell_operator_full_parent_restoration_all_seven_columns']=True
AA=[M.subs(u,sp.Rational(str(p)))/(1+sp.Rational(str(p))**2) for M in Anum]
BB=[M.subs(v,sp.Rational(str(q)))/(1+sp.Rational(str(q))**2) for M in Bnum]
state=sp.Matrix([sp.Rational(str(x)) for x in phi]); state_norm=(state.T*state)[0]
observables=[sp.kronecker_product(M,Id) for M in AA]
observables += [sp.kronecker_product(Id,M) for M in BB]
observables += [sp.kronecker_product(M,N) for M in AA for N in BB]
contact_behavior=[sp.factor((state.T*M*state)[0]/state_norm) for M in observables]
assert sum(sp.Rational(str(x))*q for x,q in zip(visible,[1]+contact_behavior))==0
checks['seed_quantum_contact_behavior_order']=['a0','a1','b0','b1','E00','E01','E10','E11']
checks['seed_quantum_contact_behavior']=[str(x) for x in contact_behavior]

# The seven chamber generators are checked in the same actual pencil.
Y=[[Q(x) for x in row] for row in chamber['auxiliary_matrix']]
eventy=[[Q(x) for x in row] for row in chamber['event_auxiliary_vectors']]
assert len(eventy)==6 and all(len(v)==302 for v in eventy)
assert all(sum(Y[i][k]*c[k] for k in range(7))==y[i] for i in range(302))
def psd_rational(M):
    M=[row[:] for row in M]; rank=0
    for i in range(len(M)):
        pivot=M[i][i]
        if pivot==0:
            assert all(M[i][j]==0 for j in range(i,len(M))); continue
        assert pivot>0; rank+=1
        for j in range(i+1,len(M)):
            for k in range(j,len(M)):
                M[j][k]-=M[j][i]*M[i][k]/pivot; M[k][j]=M[j][k]
    return rank
eventranks=[]
for k in range(6):
    src=[Q(0)]*7; src[k+1]=Q(1)
    yi=[Y[i][k+1] for i in range(302)]
    assert yi==eventy[k]
    grams=evaluate_pencil(src,yi)
    assert gram_polynomial(grams[:4])==multiplied_physical(src)
    eventranks.append([psd_rational(M) for M in grams])
checks['six_event_generator_full_identities_and_PSD_ranks']=eventranks

newdata={'scope':'Exact rational selector on the explicitly bounded open regular C1 contact domain; not a universal C1 or full Q222 lift.',
 'seed_shape':[str(x) for x in theta], 'seed_source':[str(x) for x in c],
 'seed_compressed_gram_entries':[str(x) for x in h],
 'seed_target':[str(x) for x in b0], 'fixed_free_values':[str(h[i]) for i in free],
 'selector':'h_F=h_star_F; h_J=h_star_J+P*(b(theta)_I-b_star_I)',
 'domain':'p,q>0; a*b*t!=0; det(J_e)!=0; -inverse(J_e)*J_d>0; max(abs(b(theta)-b_star))<epsilon',
 'inverse_trace_binary_bound_exponents':exponents,'kappa':kappa,'L':str(L),'epsilon':str(epsilon),
 'exact_inverse_traces_at_seed':inverse_traces,
 'source_chart_Jacobian_rank_certificate':checks['source_chart_Jacobian_rank_certificate'],
 'visible_inequality_order':checks['seed_inequality_order'],
 'visible_seed_inequality_coefficients':checks['seed_inequality_coefficients'],
 'quantum_contact_behavior_order':checks['seed_quantum_contact_behavior_order'],
 'quantum_contact_behavior':checks['seed_quantum_contact_behavior'],
 'supporting_system_file':NAMES[1],'contact_source_file':NAMES[2]}
(OUT/'RATIONAL_OPEN_CONTACT_SELECTOR.json').write_text(json.dumps(newdata,indent=2)+'\n')
checks['elapsed_seconds']=round(time.monotonic()-start,3)
(HERE/'c1_contact_selector_exact_results.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
