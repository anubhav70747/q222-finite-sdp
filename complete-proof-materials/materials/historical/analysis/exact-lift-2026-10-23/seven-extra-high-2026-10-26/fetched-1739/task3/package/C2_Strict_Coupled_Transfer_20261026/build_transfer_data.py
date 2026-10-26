#!/usr/bin/env python3
"""Generate compact complete C2 source/target bridge data from original exact map."""
import json
from pathlib import Path
import c2_fixed_map as c
HERE=Path(__file__).resolve().parent
full=json.loads((HERE/'C2_154_CONTACT_MAP.json').read_text())
trg=json.loads((HERE/'C2_CLEARED_POLYNOMIAL_TARGET.json').read_text())
blocks=[list(map(tuple,b)) for b in full['support_blocks']]
rows=list(map(tuple,full['coefficient_rows']))
assert len(rows)==412 and list(map(len,blocks))==[52,36,36,30]
assert len(trg['nonzero_polynomial_coefficients_by_row'])==24
rmap={key:i for i,key in enumerate(rows)}
def pair_classes(key,ban):
    return [dict(block=bi,left=list(i),right=list(j))
            for bi,b in enumerate(blocks)
            for ii,i in enumerate(b)
            for j in b[ii:]
            if c.product_key(i,j)==key and i!=ban and j!=ban]
D={
 'domain':'Actual strict six-event C2 native Omega as in c2_fixed_map.py',
 'status':'STRICT_C2_Z6_UNIVERSAL_PSD_UNPROVED',
 'fixed_library':'Z6, with 154 support-face rows and four PSD parity blocks',
 'active_target_rows_24':list(sorted(map(int,trg['nonzero_polynomial_coefficients_by_row']))),
 'extreme_diagonal':[
   {'key':[12,0,2,2],'row':rmap[(12,0,2,2)],'unique_feature':[6,0,2],
     'coefficient':'2*c**4*r*(2*n+1)*(b-a**2*t)'},
   {'key':[0,12,1,1],'row':rmap[(0,12,1,1)],'unique_feature':[0,6,1],
     'coefficient':'2*c**4*r*(2*n+1)*a**2*r'}],
 'carrier':{
   'positive_gap_h_ge_q':'K*((h-q)*(U*phi_Y-phi_X)**2 + q*(U*phi_Y-V*phi_X)**2)',
   'negative_gap_h_le_q':'K*((q-h)*(V*phi_X-phi_X)**2 + h*(U*phi_Y-V*phi_X)**2)',
   'K':'2*c**4*r*(2*n+1)', 'h':'b-a**2*t', 'q':'a**2*r',
   'note':'carrier PSD and contact-zero for each branch, but residual PSD unproved'},
 'residual_transfer':[
   {'sign':'h>q','key':[6,0,1,2],'row':rmap[(6,0,1,2)],
    'required_coefficient':'2*K*(h-q)',
    'surviving_pairs':pair_classes((6,0,1,2),(6,0,2))},
   {'sign':'h<q','key':[0,6,1,1],'row':rmap[(0,6,1,1)],
    'required_coefficient':'2*K*(q-h)',
    'surviving_pairs':pair_classes((0,6,1,1),(0,6,1))}
 ],
 'angular_balance':{
    'weights':['w_A=q(-1,0)/U_c','w_B=q(-1,1)*V_c/U_c','w_C=q(0,-1)/V_c',
               'w_D=q(0,1)*V_c','w_E=q(1,-1)*U_c/V_c','w_F=q(1,0)*U_c'],
    'source_only_circuit_coefficients':['k=max(w_A-w_F,0)','p=w_A-k','l=w_F-p',
      'q=w_E-k','z=w_D-k'],
    'circuits':['s+1/s-2','s/t+t/s-2','t+1/t-2','1/s+t+s/t-3','t/s+1/t+s-3'],
    'positive_claim':'p,q,z,k,l>=0 from source balance, including exact zero values'
 },
 'full_parent_multiplier':{
   'multiplier':'det(B(u,v;theta))>0',
   'B':'leading W,X,Y 3x3 block of complete four-state parent',
   'rank_one_R_rows':16,
   'wedge_roots':120,
   'angular_roots':21,
   'complete_identity':'det(B)*F = sum_{i<j}[(R_i cross R_j) dot(B*xi+z*Z)]**2 + det(P)*Z**2',
   'limitation':'det(B) depends jointly on native source and testing angles; cannot cancel in polynomial SOS or infer an affine fixed SDP cone'}
}
assert len(D['residual_transfer'][0]['surviving_pairs'])==6
assert len(D['residual_transfer'][1]['surviving_pairs'])==3
(HERE/'C2_TRANSFER_EXACT_DATA.json').write_text(json.dumps(D,indent=2,sort_keys=True)+'\n')
print('Generated C2_TRANSFER_EXACT_DATA.json, with 6/3 surviving fractional classes')
