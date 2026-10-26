"""Reproduce the three Boolean NC14 support comparisons (separate from 340)."""
from pathlib import Path
from itertools import product
import json
import sympy as s

root=Path(__file__).resolve()
receipt=root.with_name('root_boundary_identity_audit.json')
if not receipt.exists():
    receipt=root.parent.parent/'receipts/root_boundary_identity_audit.json'
data=json.loads(receipt.read_text())['constrained_outer_coefficients']
letters=s.symbols('n m c A D alpha beta')
n,m,c,A,D,alpha,beta=letters
polys={k:{tuple(map(int,key.strip('()').split(','))):s.sympify(val,locals=dict(zip(
    ['n','m','c','A','D','alpha','beta'],letters))) for key,val in row.items()} for k,row in data.items()}
result={k:True for k in data}
for values in product([0,1],repeat=7):
    nv,mv,cv,Av,Dv,av,bv=values
    expected={
      'lambda,m':{(-1,-1):av,(-1,0):Dv or av,(-1,1):Dv,
                  (0,-1):nv or Av*av,(0,1):bv or Dv*cv,(1,-1):Av*nv,
                  (1,0):Av*cv or bv*nv,(1,1):bv*cv},
      'lambda,c':{(-1,-1):av,(-1,0):Dv or av*mv,(-1,1):Dv*mv,
                  (0,-1):nv,(0,1):Dv or bv*mv,(1,-1):0,
                  (1,0):Av,(1,1):bv},
      'lambda,alpha':{(-1,-1):1,(-1,0):Dv or mv,(-1,1):Dv*mv,
                      (0,-1):Av or nv,(0,1):Dv*cv or bv*mv,(1,-1):Av*nv,
                      (1,0):Av*cv or bv*nv,(1,1):bv*cv}}
    substitutions=dict(zip(letters,values))
    for key,row in expected.items():
        for exponent,wanted in row.items():
            actual=bool(s.sympify(polys[key].get(exponent,0)).subs(substitutions)>0)
            if actual!=bool(wanted):raise AssertionError((key,values,exponent,actual,wanted))
out={'tables':result,'assignments_per_table':128,'all_passed':True,
     'scope':'Compares the printed NC14 Boolean conditions with independently expanded determinant coefficients.'}
root.with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print('All three support tables agree on all 128 assignments each.')
