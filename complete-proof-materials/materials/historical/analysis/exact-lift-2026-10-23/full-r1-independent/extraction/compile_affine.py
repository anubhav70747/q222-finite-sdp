"""Evaluate the rational affine system of the locally extracted fixed family.

Standard-library sparse rational products, with SymPy only for the nine small
source-column polynomials. No feasibility or degree search is performed.
"""
from pathlib import Path
from collections import defaultdict
from fractions import Fraction as F
import json
import sympy as s
HERE=Path(__file__).resolve().parent
d=json.loads((HERE/'fixed_family.json').read_text())
u,v,t,r,w,z=s.symbols('u v t r w zeta');variables=(u,v,t,r,w,z)
W=w-u*r;X=t-u*z;Y=r+u*w;Z=z+u*t;Du=1+u*u;Dv=1+v*v
P={}
P['lambda']=(W+v*X)**2;P['n']=Dv*X*X;P['m']=Du*(v*w-t)**2
P['C']=Du*Dv*r*r;P['A']=Dv*Y*Y;P['D']=Du*(r+v*z)**2
P['alpha']=(v*Y-Z)**2;P['beta']=Du*Dv*t*t
P['d']=(P['D']+P['alpha']+P['beta']-Dv*Z*Z)/2
source_order=['d','lambda','n','m','C','A','D','alpha','beta']
source={name:{m:F(co) for m,co in s.Poly(s.expand(P[name]),*variables).terms() if co}
        for name in source_order}
# Exact event-ray roots, independently expanded in the correct new numerator.
events={'lambda':[W+v*X],'n':[X,v*X],'m':[v*w-t,u*(v*w-t)],
        'C':[r,u*r,v*r,u*v*r],'A':[Y,v*Y],
        'D':[r+v*z,u*(r+v*z)],'alpha':[v*Y-Z],'beta':[t,u*t,v*t,u*v*t]}
assert all(s.expand(sum(f*f for f in fs)-P[name])==0 for name,fs in events.items())
source_records=[];gram_records=[];monomials=[]
chart_rows={};allcolumns=[];offset=0
for chart in [1,-1]:
    columns=[];support=set().union(*[set(p) for p in source.values()])
    for block,b in enumerate(d['blocks']):
        if b['chart']!=chart:continue
        q={tuple(m):F(c) for m,c in b['weight']};roots=b['monomials']
        for i,a in enumerate(roots):
            for j in range(i,len(roots)):
                f=roots[j];factor=1 if i==j else 2
                col={tuple(ai+fi+mi for ai,fi,mi in zip(a,f,m)):factor*c for m,c in q.items()}
                columns.append((block,i,j,col));support.update(col)
    ordered=sorted(support);rows={m:offset+i for i,m in enumerate(ordered)}
    chart_rows[chart]=len(ordered)
    monomials.extend([[chart,list(m)] for m in ordered]);offset+=len(ordered)
    for block,i,j,col in columns:
        gram_records.extend([[rows[m],block,i,j,str(c)] for m,c in sorted(col.items())])
    for name,p in source.items():source_records.extend([[rows[m],name,str(c)] for m,c in sorted(p.items())])
    allcolumns.extend(columns)
out={'conventions':'Zero-based rows, blocks and actual symmetric Gram indices; source records are RHS.',
     'source_order':source_order,'test_variables':list(map(str,variables)),
     'block_orders':[len(b['monomials']) for b in d['blocks']],
     'block_charts':[b['chart'] for b in d['blocks']],
     'monomial_rows':monomials,'gram_records':gram_records,'source_records':source_records,
     'source_nonnegativity':'Nine separate scalar blocks, one per source coordinate.',
     'target':'(1+v^2)*S_c; equal to the original quarter-parent divided by Du AFTER state rotation.'}
(HERE/'affine_system.json').write_text(json.dumps(out,separators=(',',':'))+'\n',encoding='utf8')
receipt={'total_order_with_nine_signs':d['total_PSD_order_with_nine_source_signs'],
 'gram_variables':len(allcolumns),'blocks':len(d['blocks']),'equations':len(monomials),
 'rows_per_chart':chart_rows,'gram_records':len(gram_records),'source_records':len(source_records),
 'all8_event_rays_exact':True,'all_columns_state_homogeneous_degree2':all(sum(m[2:])==2 for _,m in monomials),
 'no_unreceived_299_data_used':True}
assert receipt['gram_variables']==d['gram_variables']
assert receipt['all_columns_state_homogeneous_degree2']
(HERE/'affine_checks.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
print(json.dumps(receipt,indent=2))
