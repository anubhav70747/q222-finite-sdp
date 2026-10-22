#!/usr/bin/env python3
"""Exact block-diagonal intersection of the two lambda-positive source cones."""
import json,hashlib
from pathlib import Path
p=Path(__file__).resolve().parent
a=json.loads((p/'C1_POS_FULL_SOURCE_PENCIL.coo.json').read_text())
b=json.loads((p/'C1_NEG_FULL_SOURCE_PENCIL.coo.json').read_text())
nb=len(a['block_orders']);na=a['counts']['auxiliary_dimension']
coo=list(a['coo'])
for block,i,j,v,c in b['coo']:
 coo.append([nb+block,i,j,v if v<7 else v+na,c])
d={'schema':'Q222.C1.both_signs_fixed_local_source_pencil.v1',
 'scope':'lambda>0 uniform local C1 extraction; exact intersection of the two derived fixed full-state Gram certificate cones, not the full Q222 polar or a Q222 lift.',
 'source_order':a['source_order'],'auxiliary_dimension':na+b['counts']['auxiliary_dimension'],
 'block_orders':a['block_orders']+b['block_orders'],
 'total_matrix_order':a['total_matrix_order']+b['total_matrix_order'],
 'constant_entries':[],'coo':coo,'convention':a['convention'],
 'counts':{'blocks':nb+len(b['block_orders']),'max_block_order':12,'total_matrix_order':508,'source_coordinates':7,'auxiliaries':2021,'nonzero_affine_entries':len(coo)},
 'testing_domain':'|u|<=2^-24, |w|<=2^-12; every normalized real or complex state; positive and negative state-product charts plus the two exact angle/state sign gauges.',
 'source_admission':'See uniform_source_bounds.tex and the finite coefficient-guard theorem. No blanket admissibility of all physical C1/C2 sources is asserted.',
 'both_inclusions_for_stated_object':'A source has PSD Grams in both displayed fixed dictionaries iff this block-diagonal pencil has a PSD completion. Every such source slack is nonnegative on the stated local testing domain. Admitted physical sources have explicit PSD Gram assignments.',
 'provenance':[{'file':x,'sha256':hashlib.sha256((p/x).read_bytes()).hexdigest()} for x in ['C1_POS_FULL_SOURCE_PENCIL.coo.json','C1_NEG_FULL_SOURCE_PENCIL.coo.json']]}
(p/'C1_BOTH_SIGNS_LOCAL_SOURCE_PENCIL.coo.json').write_text(json.dumps(d,indent=2)+'\n')
assert len(coo)==5326 and sum(d['block_orders'])==508
print(json.dumps(d['counts'],indent=2))
