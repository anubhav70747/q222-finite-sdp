"""Construct a finite, typed tensor-SOS library for the complete local C1D Gram.
Each printed record specifies two simple positive tensors. Taking their roots in
separate factors and multiplying the indicated five-root combination is the
explicit SOS identity. No mixed square root is introduced.
"""
from pathlib import Path
from fractions import Fraction as F
import json
ROOT=Path(__file__).resolve().parents[1]
g=json.loads((ROOT/'data/C1D_Z_COMPLETE_SOURCE_GRAM.json').read_text())
coef_indices=[3,4,5,6,7,10,11,12] # q,b,a,r,E,e,rr,ss
coef_names=[g['variable_order'][i] for i in coef_indices]
test_indices=[0,1,2,8,9]; test_names=[g['variable_order'][i] for i in test_indices]
entries={(z['row'],z['col']):[(m,F(c)) for m,c in z['terms']] for z in g['upper_entries']}
scalar_records=[];pieces=[]
for i in range(5):
 terms=[(m,c) for m,c in entries[i,i] if any(m)]
 scalar_records.append({'direction':[int(k==i) for k in range(5)],'rho':'4','terms':terms})
for i in range(5):
 for j in range(i+1,5):
  for sign in [1,-1]:
   vec=[0]*5;vec[i]=1;vec[j]=sign
   scalar_records.append({'direction':vec,'rho':'1/2','terms':[(m,sign*c/2)for m,c in entries[i,j]]})
scalar_out=[]
for idx,z in enumerate(scalar_records):
 norm=sum((abs(c) for m,c in z['terms']),F(0));rho=F(z['rho']);scalar_out.append({'index':idx,'direction':z['direction'],'positive_reference':str(rho),'correction_l1':str(norm),'term_count':len(z['terms'])})
 if not norm:
  pieces.append({'scalar':idx,'constant_square_weight':str(rho)})
  continue
 budget=rho/norm
 for m,c in z['terms']:
  cp=[m[i]for i in coef_indices];tp=[m[i]for i in test_indices]
  side='coefficient' if m[3] or m[4] else 'testing'
  assert (m[3]+m[4]>0) if side=='coefficient' else sum(m[i] for i in test_indices)>0
  # On a scalar with coefficient beta and monomial C*T, record |beta|/2.
  # coefficient-small: (budget +/- sign*C)*(1 +/- T).
  # test-small:        (1 +/- C)*(budget +/- sign*T).
  pieces.append({'scalar':idx,'half_weight':str(abs(c)/2),'sign':1 if c>0 else -1,
    'small_side':side,'budget':str(budget),'coefficient_exponents':cp,'testing_exponents':tp})
# Base roots are finite sums of simple tensors in the two independently labelled fields.
base_roots=[
 [{'coefficient':'1','testing':'u*W'},{'coefficient':'-b**2','testing':'v'}],
 [{'coefficient':'b','testing':'u*v'}],
 [{'coefficient':'b**2','testing':'v**2'},{'coefficient':'-b**2*q**2','testing':'1'}],
 [{'coefficient':'b/e','testing':'X+u'},{'coefficient':'-a**2*b**3*q**2*e','testing':'v*W'}],
 [{'coefficient':'1/(rr*ss)','testing':'Y'},{'coefficient':'-r/(rr*ss)','testing':'u*W'}]]
# Polynomial sqrt radicands needed only on the testing cover, deduplicated exactly.
radicands={}
def add_rad(const,coef,pow):
 key=(str(const),str(coef),tuple(pow));radicands[key]=None
for z in pieces:
 if 'constant_square_weight'in z:continue
 if z['small_side']=='coefficient':
  add_rad(F(1),F(1),z['testing_exponents']);add_rad(F(1),F(-1),z['testing_exponents'])
 else:
  add_rad(F(z['budget']),F(z['sign']),z['testing_exponents']);add_rad(F(z['budget']),F(-z['sign']),z['testing_exponents'])
radlist=[{'constant':k[0],'monomial_coefficient':k[1],'testing_exponents':list(k[2])}for k in sorted(radicands)]
lookup={k:i for i,k in enumerate(sorted(radicands))}
for z in pieces:
 if 'constant_square_weight'in z:continue
 bp=tuple(z['testing_exponents'])
 if z['small_side']=='coefficient':ks=[('1','1',bp),('1','-1',bp)]
 else:ks=[(z['budget'],str(z['sign']),bp),(z['budget'],str(-z['sign']),bp)]
 z['testing_root_ids']=[lookup[k]for k in ks]
out={'scope':'EXPLICIT FINITE LOCAL tensor-SOS roots for the original full C1D parent on the stated rectangle; no global cover assertion',
 'source_gram':'C1D_Z_COMPLETE_SOURCE_GRAM.json','epsilon':g['uniform_box_epsilon'],
 'coefficient_variables':coef_names,'testing_variables':test_names,'base_roots':base_roots,'scalar_channels':scalar_out,'pieces':pieces,
 'testing_radical_relations':[{'root':f'zeta{i}','equation':'zeta^2 = constant + monomial_coefficient * testing_monomial',**z}for i,z in enumerate(radlist)],
 'physical_restoration':{'coefficient_factor':'sqrt(r*(1-r)/D)','testing_factor':'1/sqrt((1+u**2)*(1+v**2))','D':'r*(b**2*r+1)*(a**2*b**4*q**2+1)*(b**2*q**2*r-1)*(1-E*a**2*b**2*q**2+a**2*b**4*q**2*r)*(1-a**2*b**2*q**2+E*a**2*b**2*q**4)*(-1+b**2*(1-r)-a**2*b**4*q**2*(1-E*q**2))'},
 'identity_recipe':{'coefficient_small':'half_weight*((budget+sign*C) tensor (1+T) + (budget-sign*C) tensor (1-T))',
 'testing_small':'half_weight*((1+C) tensor (budget+sign*T) + (1-C) tensor (budget-sign*T))',
 'root_rule':'For each of the two nonnegative simple tensors take its two square roots independently, multiply by sqrt(half_weight) and by direction dot base_roots; sum their squares.'},
 'stats':{'scalar_channels':len(scalar_out),'nonconstant_tensor_pairs':sum('half_weight'in z for z in pieces),'individual_tensor_squares':sum(2 if 'half_weight'in z else 1 for z in pieces),'distinct_testing_radicals':len(radlist),'max_separation_rank_per_resulting_root':4},
 'universal_finite_master_feature_space':False}
(ROOT/'data/C1D_Z_EXPLICIT_SEPARATED_SQUARE_LIBRARY.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
print(json.dumps(out['stats'],indent=2))
