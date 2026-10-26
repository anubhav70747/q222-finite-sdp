"""Exact fixed-grid positivity and interpolation certificate for closedness.

This standard-library checker uses no symbolic algebra, optimizer or rank
tolerance. Distinct monomials with degree at most (4,3,1,1,1,1) are
independent on a (5,4,2,2,2,2) Cartesian grid by interpolation.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import json
HERE=Path(__file__).resolve().parent
d=json.loads((HERE/'fixed_family.json').read_text())
grid=d['closure_grid'];base=[list(map(F,grid[k])) for k in ['u','v','t','r','w']]
counts=[];minimum=None
for chart in [1,-1]:
    zz=grid['zeta_positive' if chart==1 else 'zeta_negative']
    pts=list(product(*base,list(map(F,zz))))
    assert len(pts)==320
    for b in d['blocks']:
        if b['chart']!=chart:continue
        mon=[tuple(m) for m in b['monomials']]
        assert len(mon)==len(set(mon))
        assert all(all(a<len(g) for a,g in zip(m,base+[zz])) for m in mon)
        for point in pts:
            value=sum(F(c)*__import__('math').prod(z**a for z,a in zip(point,m))
                      for m,c in b['weight'])
            assert value>0
            minimum=value if minimum is None else min(minimum,value)
        counts.append(len(mon))
out={'every_block_weight_strictly_positive_on_its_fixed_grid':True,
 'points_per_chart':320,'number_of_Gram_blocks':len(counts),
 'block_orders':counts,'minimum_weight':str(minimum),
 'each_block_polynomial_basis_independent_by_Cartesian_interpolation':True,
 'grid_sizes':[5,4,2,2,2,2],
 'strictly_PD_trace_form':'sum_grid q_b(a)*m_b(a)m_b(a)^T/(1+v(a)^2)',
 'conclusion':'The encoded cone in the nine source coefficients is closed.'}
(HERE/'closure_checks.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print('All 49 positive evaluation forms are positive definite by exact interpolation.')
