"""Exact canonical maps and conditional full-polar assembly.

Python 3 + SymPy. No solver or source sampling. --check compares every
exported rational map to its defining formula without rewriting the data.
"""
from pathlib import Path
from itertools import product, combinations
from fractions import Fraction
import argparse
import hashlib
import json
import sympy as s

HERE=Path(__file__).resolve().parent
DEFAULT=HERE.parent/'halfline'/'HALFLINE_QUADRATIC_SELECTOR_ORDER33.coo.json'
ap=argparse.ArgumentParser()
ap.add_argument('--pencil',type=Path,default=DEFAULT)
ap.add_argument('--check',action='store_true')
args=ap.parse_args()
pencil=json.loads(args.pencil.read_text())
source_order=['d','lambda','n','m','C','A','D','alpha','beta']
bell_order=['f0','fA0','fA1','fB0','fB1','fE00','fE01','fE10','fE11']
assert pencil['source_order']==source_order
event_data=[('lambda',1,1,0,1),('n',1,-1,0,0),('m',1,-1,1,1),
            ('C',-1,1,1,0),('A',-1,1,0,0),('D',-1,1,1,1),
            ('alpha',-1,-1,0,1),('beta',1,-1,1,0)]
H=s.Matrix([1,0,0,0,0,-s.Rational(1,2),s.Rational(1,2),-s.Rational(1,2),-s.Rational(1,2)])
base=[1,-1,1,1]
def event(a,b,x,y):
    v=s.zeros(9,1);v[0]=1;v[1+x]=a;v[3+y]=b;v[5+2*x+y]=a*b
    return v                         # Coefficients of 4p(a,b|x,y).
events=[event(a,b,x,y) for _,a,b,x,y in event_data]
M=s.Matrix.hstack(H,*events)
assert M.det()!=0
locals=[]
for a0,a1,b0,b1 in product((-1,1),repeat=4):
    locals.append(s.Matrix([1,a0,a1,b0,b1,a0*b0,a0*b1,a1*b0,a1*b1]))
neighbors=[];neighbor_indices=[]
for i,e in enumerate(events):
    found=[j for j,v in enumerate(locals) if (H.T*v)[0]==0 and (e.T*v)[0]==4]
    assert len(found)==1
    j=found[0];neighbor_indices.append(j);neighbors.append(locals[j])
    assert [(a.T*locals[j])[0] for a in events]==[4*int(i==k) for k in range(8)]
pr0=s.Matrix([1,0,0,0,0,*base])
V=s.Matrix.vstack(pr0.T,*(a.T for a in neighbors))
assert V*M==s.diag(-1,*([4]*8))
Minv=s.diag(-1,*([s.Rational(1,4)]*8))*V
assert Minv*M==M*Minv==s.eye(9)

relabels=[];pr=[]
for r0,r1,s0 in product((-1,1),repeat=3):
    diagonal=[1,r0,r1,s0,1,r0*s0,r0,r1*s0,r1]
    R=s.diag(*diagonal)
    rg=R*pr0
    assert list(rg[5:])==[r0*s0,-r0,r1*s0,r1]
    assert s.prod(rg[5:])==-1 and R*R==s.eye(9)
    pr.append(rg)
    relabels.append({'A_output_signs':[r0,r1],'B_output_signs':[s0,1],
        'coordinate_and_coefficient_diagonal':diagonal,
        'PR_expectations':list(map(int,rg[1:])),
        'source_to_bell':[[str(a) for a in row] for row in (R*M).tolist()],
        'bell_to_source':[[str(a) for a in row] for row in (Minv*R).tolist()]})
assert len({tuple(a[5:]) for a in pr})==8

# Every pairwise PR midpoint is an explicit convex combination of local
# deterministic vertices. Nonantipodal pairs have two shared edges;
# their four compatible deterministic assignments have the required mean.
midpoints=[]
for i,j in combinations(range(8),2):
    shared=[k for k in range(4) if pr[i][5+k]==pr[j][5+k]]
    assert len(shared) in (0,2)
    compatible=[k for k,v in enumerate(locals) if all(v[5+e]==pr[i][5+e] for e in shared)]
    assert len(compatible)==(16 if not shared else 4)
    avg=sum((locals[k] for k in compatible),s.zeros(9,1))/len(compatible)
    assert avg==(pr[i]+pr[j])/2
    midpoints.append({'PR_pair':[i,j],'local_indices':compatible,'common_weight':str(s.Rational(1,len(compatible)))})

# Two opposite PR orientations cover all 16 event generators on d=0.
all_events=[(a,b,x,y,event(a,b,x,y)) for x,y in product(range(2),repeat=2) for a,b in product((-1,1),repeat=2)]
id_index=next(i for i,r in enumerate(relabels) if r['A_output_signs']==[1,1] and r['B_output_signs']==[1,1])
op_index=next(i for i,r in enumerate(relabels) if r['A_output_signs']==[-1,-1] and r['B_output_signs']==[1,1])
assignments=[]
for a,b,x,y,e in all_events:
    hits=[]
    for g in (id_index,op_index):
        MG=s.diag(*relabels[g]['coordinate_and_coefficient_diagonal'])*M
        for i in range(8):
            if MG[:,1+i]==e:hits.append((g,i))
    assert len(hits)==1
    g,i=hits[0]
    assignments.append({'event':[a,b,x,y],'chart':g,'source_coordinate':source_order[i+1],
                        'coefficient_normalization':'source coefficient1 gives 4p; coefficient1/4 gives p'})

# Adjoint coordinates are the independent upper-triangular entries of
# the four core blocks. Full symmetric COO entries automatically supply
# the required factor2 for an off-diagonal coordinate in trace(A Y).
orders=[8,8,4,4,1,1,1,1,1,1,1,1,1]
assert pencil['blocks']==orders and pencil['auxiliary_count']==44
coords=[];lookup={};offset=0
for block_index,n in enumerate(orders):
    for i in range(n):
        for j in range(i,n):
            lookup[(offset+i,offset+j)]=len(coords)
            coords.append([block_index,i,j])
    offset+=n
adj=s.zeros(53,101)
for k,i,j,a in pencil['entries']:
    adj[k,lookup[(min(i,j),max(i,j))]]+=a
assert adj.rank()==53
assert adj[:9,92:]==s.eye(9) and adj[9:,92:]==s.zeros(44,9)
B=adj[9:,:92];Ac=adj[:9,:92]
assert B.rank()==44
null=B.nullspace();N=s.Matrix.hstack(*null)
assert N.shape==(92,48) and B*N==s.zeros(44,48) and N.rank()==48
source_aux=Ac*N
assert source_aux.rank()==9
def sparse(A):return [[i,j,str(A[i,j])] for i in range(A.rows) for j in range(A.cols) if A[i,j]]
def mat(A):return [[str(a) for a in row] for row in A.tolist()]
data={'format':'q222-conditional-full-polar-assembly-v1',
 'scope':'Conditional wrapper only: exact global equality requires universal canonical quadratic half-line adequacy. Unconditionally this wrapper is an outer approximation of Q222 contained in the nonsignalling polytope.',
 'bell_coefficient_order':bell_order,'expectation_order':[a[1:] for a in bell_order[1:]],
 'source_order':source_order,'event_order':[list(a) for a in event_data],
 'source_to_bell':mat(M),'bell_to_source':mat(Minv),'source_map_determinant':str(M.det()),
 'base_PR_homogeneous':list(map(int,pr0)),'local_deterministic_homogeneous':[list(map(int,a)) for a in locals],
 'canonical_neighbor_local_indices':neighbor_indices,
 'relabelings':relabels,'PR_midpoint_local_decompositions':midpoints,
 'event_generator_assignments_to_two_charts':assignments,
 'base_chart_index':id_index,'antipodal_chart_index':op_index,
 'pencil_sha256':hashlib.sha256(args.pencil.read_bytes()).hexdigest(),
 'dual':{'inner_product':'trace; off-diagonal upper-triangular coordinates have coefficient2',
   'full_independent_Gram_coordinates':101,'core_independent_Gram_coordinates':92,
   'full_adjoint_rank':53,'auxiliary_adjoint_rank':44,'source_adjoint_on_kernel_rank':9,'source_sign_adjoint_identity':True,
   'core_coordinates':coords[:92],'core_auxiliary_adjoint_sparse':sparse(B),
   'core_source_adjoint_sparse':sparse(Ac),'core_nullspace_basis_shape':[92,48],
   'core_nullspace_basis_sparse':sparse(N),'source_auxiliary_product_sparse':sparse(source_aux),
   'chart_constraint':'core Y=N*t PSD in blocks8,8,4,4; source_signs=M_g^T*z-Acore*N*t >=0; z=(1,q)',
   'per_chart_block_orders':orders,'charts':8,'total_block_order':264,'free_auxiliaries_per_chart':48,'total_free_auxiliaries':384},
 'exact_checks':{'source_inverse':True,'nine_vertex_source_identity':True,'all8_PR_orientations':True,
   'all28_PR_pair_midpoints_local':True,'all16_event_generators_in_two_zero_d_charts':True,
   'full_adjoint_rank53':True,'core_auxiliary_rank44':True,'core_kernel_dimension48':True}}
out=HERE/'full_polar_maps.json'
if args.check:assert json.loads(out.read_text())==data
else:out.write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps({'mode':'check' if args.check else 'build','source_map_determinant':str(M.det()),
 'relabelings':8,'PR_midpoint_decompositions':28,'events':16,'full_adjoint_rank':53,'core_auxiliary_rank':44,
 'dual_core_kernel_dimension':48,'conditional_total_order':264,'conditional_auxiliary_count':384,
 'map_file_bytes':out.stat().st_size},indent=2))
