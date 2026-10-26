"""One small diagnostic solve; numerical output is not an adequacy proof."""
from pathlib import Path
from fractions import Fraction
import json
import numpy as np
import cvxpy as cp

here=Path(__file__).resolve().parent
d=json.loads((here/'HALFLINE_QUADRATIC_SELECTOR_ORDER24.coo.json').read_text())
c=np.array([359,0,0,75/16,60,416,41/2,831/2,151/2],dtype=float)/359
arrays=np.zeros((53,24,24))
for k,i,j,a in d['entries']:arrays[k,i,j]+=a
A=np.tensordot(c,arrays[:9],axes=1)
y=cp.Variable(44)
t=cp.Variable()
M=cp.Constant(A)
for k in range(44):M=M+y[k]*arrays[k+9]
constraints=[];offset=0
for n in d['blocks']:
    constraints.append(M[offset:offset+n,offset:offset+n]+t*np.eye(n)>>0)
    offset+=n
problem=cp.Problem(cp.Minimize(t),constraints)
problem.solve(solver='CLARABEL',max_iter=150,tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9)
out={'scope':'One numerical diagnostic for the fixed simple source; no exact feasibility or global-adequacy claim.',
     'source_order':d['source_order'],'source_d1':c.tolist(),
     'status':problem.status,'t':None if t.value is None else float(t.value),
     'solve_seconds':problem.solver_stats.solve_time,
     'auxiliary_order':d['auxiliary_order'],'auxiliary':None if y.value is None else y.value.tolist(),
     'dual_blocks':[None if z.dual_value is None else z.dual_value.tolist() for z in constraints]}
if y.value is not None:
    R=A+np.tensordot(y.value,arrays[9:],axes=1)
    offset=0;out['unshifted_min_eigenvalues']=[]
    for n in d['blocks']:
        out['unshifted_min_eigenvalues'].append(float(np.linalg.eigvalsh(R[offset:offset+n,offset:offset+n])[0]));offset+=n
(here/'simple_source_numerical_diagnostic.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['auxiliary_order','auxiliary','dual_blocks']},indent=2))
