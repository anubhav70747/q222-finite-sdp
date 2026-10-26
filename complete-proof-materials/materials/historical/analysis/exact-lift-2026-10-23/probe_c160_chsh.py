"""One bounded diagnostic of the evaluated C160 candidate; not a proof."""
from pathlib import Path
from fractions import Fraction as F
import json, time, sys
import numpy as np
from scipy.sparse import coo_matrix
import cvxpy as cp

HERE=Path(__file__).resolve().parent
DATA=HERE.parents[1]/'github-q222-finite-sdp/extracted-lifts/outer-candidates'
d=json.loads((DATA/'C160_source_compiler_cache.json').read_text())
refs=d['tensor_phi_basis']['class_reference_coords']
refmap={c:(r,F(w)) for r,(c,w) in enumerate(refs)}
cols={c:{r:w} for c,(r,w) in refmap.items()}
for v in d['nullspace_basis']:
    col={}
    for c,n,den in v['entries']:
        if c==v['free']:continue
        r,w=refmap[c];col[r]=col.get(r,F(0))-F(n,den)*w
    cols[v['free']]={r:w for r,w in col.items() if w}
L=np.zeros((9,243))
for j,vec in enumerate(d['target_affine_section']):
    for c,n,den in vec:
        for r,w in cols[c].items():L[j,r]+=float(F(n,den)*w)*(1 if j==0 else -1)
Y=cp.Variable(243)
constraints=[L[0]@Y==1]
maps=[]
for b in d['block_layout']:
    size=b['size'];rows=[];columns=[];values=[];c=b['coord_start']
    for i in range(size):
        for j in range(i,size):
            for r,w in cols[c].items():
                value=float(w)/(1 if i==j else 2)
                rows.append(i+j*size);columns.append(r);values.append(value)
                if i!=j:rows.append(j+i*size);columns.append(r);values.append(value)
            c+=1
    A=coo_matrix((values,(rows,columns)),shape=(size*size,243)).tocsr()
    maps.append((size,A))
    constraints.append(cp.reshape(A@Y,(size,size),order='F')>>0)
chsh=np.array([0,0,0,0,0,1,-1,1,1])
problem=cp.Problem(cp.Maximize(chsh@L@Y),constraints)
start=time.time()
strict='--strict' in sys.argv
suffix='_strict' if strict else ''
try:
    problem.solve(solver='MOSEK',verbose=False,mosek_params={
        'MSK_DPAR_OPTIMIZER_MAX_TIME':30.0,'MSK_IPAR_NUM_THREADS':2,
        'MSK_DPAR_INTPNT_CO_TOL_REL_GAP':1e-11 if strict else 1e-9,
        'MSK_DPAR_INTPNT_CO_TOL_PFEAS':1e-11 if strict else 1e-8,
        'MSK_DPAR_INTPNT_CO_TOL_DFEAS':1e-11 if strict else 1e-8})
    report={'scope':'Single numerical diagnostic; not a reverse-inclusion proof or exact certificate',
            'status':problem.status,'seconds':time.time()-start,'optimum':problem.value,
            'quantum_chsh_bound':float(2*np.sqrt(2))}
    if Y.value is not None:
        report['behavior']=(L[1:]@Y.value).tolist()
        report['normalization_residual']=float(L[0]@Y.value-1)
        report['minimum_primal_eigenvalue']=float(min(np.linalg.eigvalsh((A@Y.value).reshape((n,n),order='F')).min() for n,A in maps))
        report['normalization_dual']=float(constraints[0].dual_value)
        np.savez(HERE/f'c160_chsh_diagnostic_arrays{suffix}.npz',moments=Y.value,behavior=L[1:]@Y.value)
except Exception as e:
    report={'status':'not_solved','error':str(e),'seconds':time.time()-start}
(HERE/f'c160_chsh_diagnostic{suffix}.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
