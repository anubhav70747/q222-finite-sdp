"""Exact primal/dual regularity of the supplied ramified order-49 cone.

Python 3 + SymPy (rank only). No solver, source sweep, or adequacy test.
Reads the actual upper-triangular affine pencil and the frozen order-33
Slater witness. Writes only local witness/check data beside this script.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import json
import sympy as sp

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
ap=argparse.ArgumentParser()
ap.add_argument('--pencil',type=Path,default=BASE/'data'/'ramified_order49_pencil.json')
ap.add_argument('--old-slater',type=Path,default=BASE.parent/'full-polar'/'slater_witness.json')
args=ap.parse_args()
p=json.loads(args.pencil.read_text());old=json.loads(args.old_slater.read_text())
assert p['order']==49 and p['auxiliary_dimension']==148
assert [b['order'] for b in p['blocks']]==[1]*9+[12,8,12,8]
assert p['source_names']==old['source_order']

def zero(n,m=None):return [[Q(0) for _ in range(n if m is None else m)] for _ in range(n)]
def eye(n):return [[Q(int(i==j)) for j in range(n)] for i in range(n)]
def mat(A):return [[Q(a) for a in row] for row in A]
def tr(A):return [list(a) for a in zip(*A)]
def add(A,B):return [[a+b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def scale(k,A):return [[k*a for a in row] for row in A]
def extract(A,i,j,n,m=None):return [row[j:j+(n if m is None else m)] for row in A[i:i+n]]
def place(A,B,i,j):
    for r,row in enumerate(B):
        for s,a in enumerate(row):A[i+r][j+s]=a
def sparse_matrix(entries,n=49):
    A=zero(n)
    for i,j,a in entries:
        assert i<=j
        A[i][j]+=Q(a)
        if i!=j:A[j][i]+=Q(a)
    return A
def dot(A,B):return sum(a*b for ar,br in zip(A,B) for a,b in zip(ar,br))
def serialize(A):return [[str(a) for a in row] for row in A]
def psd_pivots(A):
    R=[row[:] for row in A];out=[]
    assert R==tr(R)
    for i in range(len(R)):
        d=R[i][i];out.append(d);assert d>=0
        if not d:
            assert all(R[i][j]==0 for j in range(i+1,len(R)))
            continue
        for j in range(i+1,len(R)):
            for k in range(j,len(R)):
                R[j][k]-=R[j][i]*R[i][k]/d
                R[k][j]=R[j][k]
    return out
def coefficients(A,B):
    return [[[sum(A[4*a+i][4*b+j] for a in range(3) for b in range(3) if a+b==k)
             +sum(B[4*a+i][4*b+j] for a in range(2) for b in range(2) if a+b+1==k)
             for j in range(4)] for i in range(4)] for k in range(5)]

# Moment functional for integral_0^1 trace(Q(s)) ds.  The nine scalar
# source-sign slots are first in the actual pencil.
H3=[[Q(1,i+j+1) for j in range(3)] for i in range(3)]
H2=[[Q(1,i+j+2) for j in range(2)] for i in range(2)]
assert all(a>0 for a in psd_pivots(H3)+psd_pivots(H2))
assert [sum(sp.Matrix(H).inv()[i,i] for i in range(len(H))) for H in (H3,H2)]==[381,54]
W=zero(49)
for b,H in zip(p['blocks'][9:],[H3,H2,H3,H2]):
    for i,row in enumerate(H):
        for j,a in enumerate(row):
            for k in range(4):W[b['offset']+4*i+k][b['offset']+4*j+k]=a
source_mats=[sparse_matrix(a) for a in p['source_matrices']]
aux_mats=[sparse_matrix(a) for a in p['auxiliary_matrices']]
assert [dot(W,A) for A in source_mats]==[Q(2088,5)]*9
assert [dot(W,B) for B in aux_mats]==[0]*148
Wfull=[row[:] for row in W]
for i in range(9):Wfull[i][i]=1
assert [dot(Wfull,A) for A in source_mats]==[Q(2093,5)]*9
assert [dot(Wfull,B) for B in aux_mats]==[0]*148

# Literal strict primal point obtained by composing v=16s^2 with the
# old Q(v)=[I,vI]S[I,vI]^T+vT certificate.
source=list(map(Q,old['source']))
old_blocks={a['name']:mat(a['matrix']) for a in old['blocks']}
Jold=[mat(a) for a in old['J_coefficient_matrices']]
Jnew=[Jold[0],zero(2),scale(16,Jold[1]),zero(2),scale(256,Jold[2])]
values={f'J{k}_{i}{j}':Jnew[k][i][j] for k in range(5) for i in range(2) for j in range(2)}
new_blocks=[];baseline_blocks=[];proof={}
for side in ['G','H']:
    S=old_blocks['S'+side];T=old_blocks['T'+side];Bold=mat(old['B'+side])
    Q0=extract(S,0,0,4);Q2=extract(S,4,4,4)
    Q1=add(add(Bold,tr(Bold)),T)
    A0=zero(12)
    place(A0,Q0,0,0);place(A0,scale(16,T),4,4);place(A0,scale(256,Q2),8,8)
    place(A0,scale(16,Bold),0,8);place(A0,scale(16,tr(Bold)),8,0)
    A=[row[:] for row in A0];B=scale(Q(1,100),eye(8))
    for i in range(4):
        A[i][i+4]-=Q(1,200);A[i+4][i]-=Q(1,200)
        A[i+4][i+8]-=Q(1,200);A[i+8][i+4]-=Q(1,200)
    assert coefficients(A,B)==[Q0,zero(4),scale(16,Q1),zero(4),scale(256,Q2)]
    pp0=psd_pivots(add(A0,scale(-Q(8,25),eye(12))))
    pp=psd_pivots(add(A,scale(-Q(31,100),eye(12))))
    assert max(sum(abs(a-b) for a,b in zip(ar,br)) for ar,br in zip(A,A0))==Q(1,100)
    for i in range(4):
        for j in range(i,4):
            values[f'{side}_B00_{i}{j}']=B[i][j]
            values[f'{side}_B11_{i}{j}']=B[4+i][4+j]
        for j in range(4):
            values[f'{side}_B01_{i}{j}']=B[i][4+j]
            values[f'{side}_A02_{i}{j}']=A[i][8+j]
        for j in range(i+1,4):
            values[f'{side}_K01_{i}{j}']=(A[i][4+j]-A[j][4+i])/2
            values[f'{side}_K12_{i}{j}']=(A[4+i][8+j]-A[4+j][8+i])/2
    assert all(values[f'{side}_{a}_{i}{j}']==0 for a in ['K01','K12'] for i in range(4) for j in range(i+1,4))
    baseline_blocks.append(A0);new_blocks.extend([A,B])
    proof[side]={'baseline_minus_8over25I_pivots':list(map(str,pp0)),
                 'final_minus_31over100I_pivots':list(map(str,pp)),
                 'perturbation_max_abs_row_sum':'1/100'}
assert len(values)==148 and set(values)==set(p['auxiliary_names'])
aux=[values[a] for a in p['auxiliary_names']]
L=zero(49)
for a,M in zip(source+aux,source_mats+aux_mats):
    for i in range(49):
        for j in range(49):L[i][j]+=a*M[i][j]
for i in range(9):assert L[i][i]==source[i]>=Q(1,100)
for b,M in zip(p['blocks'][9:],new_blocks):assert extract(L,b['offset'],b['offset'],b['order'])==M
assert all(a>=0 for a in psd_pivots(add(L,scale(-Q(1,100),eye(49)))))

# Actual trace adjoint ranks, with coefficient2 on off-diagonal core
# coordinates. This confirms that no auxiliary direction is invisible.
coords=[]
for b in p['blocks'][9:]:
    coords.extend((b['offset']+i,b['offset']+j) for i in range(b['order']) for j in range(i,b['order']))
assert len(coords)==228
adj=sp.Matrix([[sp.Rational(M[i][j])*(1 if i==j else 2) for i,j in coords] for M in source_mats+aux_mats])
assert adj[9:,:].rank()==148
full=adj.row_join(sp.Matrix.vstack(sp.eye(9),sp.zeros(148,9)))
assert full.rank()==157

out={'scope':'Regularity and exact dimensions of the supplied sufficient order49 cone only; no universal adequacy claim.',
 'pencil_sha256':hashlib.sha256(args.pencil.read_bytes()).hexdigest(),
 'old_slater_sha256':hashlib.sha256(args.old_slater.read_bytes()).hexdigest(),
 'source_order':p['source_names'],'source':list(map(str,source)),
 'auxiliary_order':p['auxiliary_names'],'auxiliary_values':list(map(str,aux)),
 'strict_primal_core_blocks':{b['name']:serialize(M) for b,M in zip(p['blocks'][9:],new_blocks)},
 'strict_primal_bounds':{'A_G':'31/100','B_G':'1/100','A_H':'31/100','B_H':'1/100','source_signs':'1/100','full_order49':'1/100'},
 'strict_primal_bound_certificates':proof,
 'moment_dual':{'H3':serialize(H3),'H2prime':serialize(H2),
   'core_block_order':['H3 tensor I4','H2prime tensor I4','H3 tensor I4','H2prime tensor I4'],
   'H3_inverse_trace':'381','H2prime_inverse_trace':'54','core_lower_bound':'1/381',
   'core_source_contractions':['2088/5']*9,'full_source_contractions_with_sign_I9':['2093/5']*9,
   'all148_auxiliary_contractions_zero':True},
 'ranks':{'independent_core_Gram_coordinates':228,'auxiliary_adjoint_rank':148,
   'independent_all_Gram_coordinates':237,'full_adjoint_rank':157,'per_chart_free_dual_coordinates':80},
 'eight_chart_outer_wrapper':{'total_order':392,'free_auxiliaries':640,'maximum_block_order':12,'block_count':104,
   'scope':'Outer approximation only; equality requires universal canonical adequacy.'}}
(HERE/'ramified_regularity_witness.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'all148_auxiliaries_reconstruct_actual_pencil':True,
 'all148_auxiliary_moment_contractions_zero':True,'all9_source_core_moments':'2088/5',
 'full_order49_primal_margin':'1/100','core_auxiliary_adjoint_rank':148,
 'full_adjoint_rank':157,'per_chart_dual_auxiliaries':80,
 'conditional_eight_chart_order':392,'eight_chart_outer_auxiliaries':640},indent=2))
