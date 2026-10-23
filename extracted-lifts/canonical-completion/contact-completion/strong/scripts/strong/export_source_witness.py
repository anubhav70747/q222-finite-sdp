"""Export the evaluated six-Gram source-only entries as plain text JSON."""
import json,pickle,sympy as S
from pathlib import Path
O=Path(__file__).resolve().parents[2]/'data'
D=pickle.load((O/'lII_normalized_full.pkl').open('rb'))
Q=D['Q'];u,v=D['variables'][:2]
records=[];support=set()
for i in range(6):
    for j in range(i,6):
        P=S.Poly(Q[i,j],u,v)
        for exponents,coefficient in P.terms():
            if coefficient==0:continue
            support.add(exponents)
            records.append([i,j,list(exponents),str(S.factor(coefficient))])
assert len(records)==51
assert len(support)==14
data={'test_variables':['u','v'],'matrix_order':6,
      'source_parameterization':'see source-box proof and extract_lII.py',
      'upper_entries_by_test_monomial':records,
      'testing_support':[list(x) for x in sorted(support)]}
(O/'LII_NORMALIZED_SOURCE_WITNESS.json').write_text(json.dumps(data,indent=2)+'\n')
print('Exported',len(records),'upper-triangular coefficient records.')
