"""Read the square registry independently, reconstruct every scalar and matrix
coefficient, check all rectangular sign guards, and export the finite local
TESTING-side algebra and spanning list. Does not assert global adequacy.
"""
from pathlib import Path
from fractions import Fraction as F
import json,collections
ROOT=Path(__file__).resolve().parents[1]
import os
OUT=Path(os.environ.get('Q222_REPORT_DIR', str(ROOT/'verification_output')));OUT.mkdir(parents=True,exist_ok=True)
g=json.loads((ROOT/'data/C1D_Z_COMPLETE_SOURCE_GRAM.json').read_text())
d=json.loads((ROOT/'data/C1D_Z_EXPLICIT_SEPARATED_SQUARE_LIBRARY.json').read_text()); checks=0

def ck(p,msg):
 global checks
 if not p:raise AssertionError(msg)
 checks+=1
inds_c=[g['variable_order'].index(t)for t in d['coefficient_variables']];inds_t=[g['variable_order'].index(t)for t in d['testing_variables']]
scalars=[collections.defaultdict(F)for _ in d['scalar_channels']];zero=(0,)*13;eps=F(d['epsilon'])
# The base separated expansion is checked independently from its displayed five roots.
import sympy as s
vs=s.symbols(' '.join(g['variable_order']));loc={str(v):v for v in vs}
for i,terms in enumerate(d['base_roots']):
 expr=sum(s.sympify(t['coefficient'],locals=loc)*s.sympify(t['testing'],locals=loc)for t in terms)
 target=s.sympify(g['root_vector'][i],locals=loc)
 delta=s.cancel(expr-target).subs(loc['e']**2,loc['E'])
 ck(s.cancel(delta)==0,'separated base root')
 for term in terms:
  ck({str(v)for v in s.sympify(term['coefficient'],locals=loc).free_symbols}<=set(d['coefficient_variables']),'pure coefficient ownership')
  ck({str(v)for v in s.sympify(term['testing'],locals=loc).free_symbols}<=set(d['testing_variables']),'pure testing ownership')
feature_pairs=set();atom_list=['1','u*W','v','u*v','v**2','X+u','v*W','Y'];atom_ids={x:i for i,x in enumerate(atom_list)}
maxdeg=0;total_square_roots=0
for z in d['pieces']:
 scalar=scalars[z['scalar']]
 if 'constant_square_weight'in z:
  scalar[zero]+=F(z['constant_square_weight']);total_square_roots+=1;continue
 w=F(z['half_weight']);budget=F(z['budget']);sign=z['sign'];mm=[0]*13
 for i,k in enumerate(inds_c):mm[k]=z['coefficient_exponents'][i]
 for i,k in enumerate(inds_t):mm[k]=z['testing_exponents'][i]
 ck(w>0 and sign in[-1,1],'positive square coefficient');ck(budget>eps,'strict small-factor guard')
 if z['small_side']=='coefficient':
  ck(mm[3]+mm[4]>0,'coefficient monomial genuinely small')
  expected=[('1','1'),('1','-1')]
 else:
  ck(sum(mm[k]for k in inds_t)>0,'testing monomial genuinely small')
  expected=[(z['budget'],str(sign)),(z['budget'],str(-sign))]
 for ii,rootid in enumerate(z['testing_root_ids']):
  rad=d['testing_radical_relations'][rootid]
  ck((rad['constant'],rad['monomial_coefficient'])==expected[ii] and rad['testing_exponents']==z['testing_exponents'],'testing radical binding')
  for j,co in enumerate(d['scalar_channels'][z['scalar']]['direction']):
   if co:
    for atom in d['base_roots'][j]:feature_pairs.add((rootid,atom_ids[atom['testing']]))
 scalar[zero]+=2*w*budget;scalar[tuple(mm)]+=2*w*sign;total_square_roots+=2
 maxdeg=max(maxdeg,sum(z['testing_exponents']))
# Assemble the exact diagonal-dominance decomposition.
Q=[[collections.defaultdict(F)for _ in range(5)]for _ in range(5)]
for channel,poly in zip(d['scalar_channels'],scalars):
 ck(poly[zero]==F(channel['positive_reference']),'exact scalar constant')
 for i,x in enumerate(channel['direction']):
  for j,y in enumerate(channel['direction']):
   if x*y:
    for mm,co in poly.items():Q[i][j][mm]+=x*y*co
for v in g['upper_entries']:
 target={tuple(m):F(c)for m,c in v['terms']};got={m:c for m,c in Q[v['row']][v['col']].items()if c}
 ck(got==target,'full source Gram reconstruction')
ck(total_square_roots==d['stats']['individual_tensor_squares'],'square count')
# Algebraic testing cover of the whole local closed testing box.
algebra={'scope':'fixed common local testing algebra and explicit spanning list, valid only on the printed coefficient/testing chart',
 'base_variables':d['testing_variables'],'box_radius':str(eps),
 'box_relations':[f'box_{j}**2+{x}**2-({eps})**2=0'for j,x in enumerate(d['testing_variables'])],
 'normalizer_variable':'kappa','normalizer_relation':'kappa**2*(1+u**2)*(1+v**2)*(1+W**2+X**2+Y**2)-1=0',
 'energy_normalizer_variant':'Replace the normalizer relation by kappa**2*(1+u**2)*(1+v**2)-1=0 for the unnormalized state quadratic form.',
 'radical_relations':d['testing_radical_relations'],
 'testing_atoms':atom_list,'feature_pairs':[list(z)for z in sorted(feature_pairs)],
 'feature_formula':'kappa*zeta[radical_id]*testing_atoms[atom_id]',
 'feature_count':len(feature_pairs),'spanning_list_not_linearly_independent_basis':True,
 'algebra_variable_count':5+5+1+len(d['testing_radical_relations']),
 'relation_count':5+1+len(d['testing_radical_relations']),
 'max_relation_degree':max(8,maxdeg),
 'physical_test_projection':'(u,v,W,X,Y) maps to original normalized state (W,X,Y,1)/sqrt(1+W^2+X^2+Y^2) and rational real qubit measurements. All sign choices of box and radical variables are harmless.',
 'real_image':'the entire closed testing box |u|,|v|,|W|,|X|,|Y|<=epsilon',
 'global_coverage_proved':False}
ck(algebra==json.loads((ROOT/'data/C1D_Z_LOCAL_TESTING_ALGEBRA.json').read_text()),'independent complete local algebra reconstruction')
(OUT/'C1D_Z_LOCAL_TESTING_ALGEBRA.json').write_text(json.dumps(algebra,indent=2)+'\n')
report={'status':'PASS','exact_assertions':checks,'local_squares':total_square_roots,'distinct_testing_radicals':len(d['testing_radical_relations']),
 'testing_spanning_functions':len(feature_pairs),'local_algebra_variables':algebra['algebra_variable_count'],'local_relations':algebra['relation_count'],
 'max_relation_degree':algebra['max_relation_degree'],'scope':'full local identity and all sign guards; not global master coverage'}
(OUT/'LOCAL_SQUARE_LIBRARY_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
