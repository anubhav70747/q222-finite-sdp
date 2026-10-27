"""Tiny independent carrier dictionary and right-inverse norm audit."""
from pathlib import Path
from fractions import Fraction as Q
import json,time
import sympy as s
START=time.perf_counter();OUT=Path(__file__).resolve().parent
data=json.loads((OUT/'private-copy/fixed_map.json').read_text());res=json.loads((OUT/'private-copy/exact_results.json').read_text())
old={tuple(f) for block in data['blocks'] for f in block};opened=old|{(2,0,1)}
root1={(2,0,1):-1,(0,2,1):s.Symbol('q'),(0,2,0):1-s.Symbol('q')}
root2={(1,1,0):1,(1,1,1):-1}
assert set(root1)<=opened and set(root2)<=opened
S,T,W,X,Y,Z,q=s.symbols('S T W X Y Z q')
state=[W,X,Y,Z]
def polynomial(root):return sum(coef*S**i*T**j*state[l] for (i,j,l),coef in root.items())
carrier_squares=polynomial(root1)**2+2*q*polynomial(root2)**2
carrier_formula=(-S**2*X+q*T**2*X+(1-q)*T**2*W)**2+2*q*S**2*T**2*(W-X)**2
assert s.expand(carrier_squares-carrier_formula)==0
weighted=opened|{(i,j+2,l) for i,j,l in opened}
assert (len(old),len(opened),len(weighted))==(27,28,49)
assert all((3*i,3*j,l) in {(a,b,l0) for a in range(7) for b in range(13) for l0 in range(4)} for i,j,l in weighted)
# E's nonzero rows are the stored inverse rows at its selected upper-coordinate
# positions; independently compute each block's maximum absolute row sum.
inv=[[Q(v) for v in row] for row in data['inverse_pivot_minor']]
pc=data['pivot_columns'];positions=data['upper_positions'];bounds=[Q(0)]*4
assert len(inv)==len(pc)==66 and all(len(row)==66 for row in inv)
for column,row in zip(pc,inv):
    block=positions[column][0]
    bounds[block]=max(bounds[block],sum(abs(x) for x in row))
assert bounds==[Q(141),Q(0),Q(0),Q(0)]
assert [str(x) for x in bounds]==res['right_inverse_block_row_bounds']
report={'result':'BOUNDED INDEPENDENT CARRIER/NORM CHECK PASSED','carrier_roots':2,'old_features':len(old),'opened_features':len(opened),'weighted_union_features':len(weighted),'fixed_K364_dictionary_membership':True,'carrier_identity':'first root square + 2q times half-power second root square','positivity_domain':'q>0 and positive U,V; half-power root representation, not integer-polynomial SOS on all real U,V','right_inverse_block_row_bounds':list(map(str,bounds)),'elapsed_seconds':round(time.perf_counter()-START,6)}
(OUT/'bounded_carrier_norm_results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
