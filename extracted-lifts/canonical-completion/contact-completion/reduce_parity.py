"""Exact symmetry reduction of the independently reconstructed contact pencil.

Averaging the eight radical sign changes sets cross-character Gram entries
to zero. This preserves both feasibility directions and all source columns.
It does not establish admission of any additional source.
"""
from pathlib import Path
import json

here = Path(__file__).resolve().parent
d = json.loads((here/'contact_affine_independent.json').read_text())
groups = []
for block, order in enumerate(d['block_orders'], 1):
    if block <= 27:
        parts = [list(range(1, order+1))]
    elif block == 28:
        parts = [list(range(1,27)), [27,28], list(range(29,33)), [33,34]]
    else:
        parts = [list(range(1,15)), [15,16]]
    groups.extend((block, part) for part in parts)

position = {}
for newblock, (oldblock, part) in enumerate(groups, 1):
    for newindex, oldindex in enumerate(part, 1):
        position[oldblock, oldindex] = (newblock, newindex)

evenrows = [i for i,m in enumerate(d['monomials'],1) if m[6:] == [0,0,0]]
rowmap = {old:new for new,old in enumerate(evenrows,1)}
gram=[]
discarded=0
for row, block, i, j, value in d['gram_records']:
    bi, ii = position[block,i]
    bj, jj = position[block,j]
    assert (row in rowmap) == (bi == bj)
    if bi == bj:
        gram.append([rowmap[row],bi,ii,jj,value])
    else:
        discarded += 1
source=[]
for row, name, value in d['source_records']:
    assert row in rowmap
    source.append([rowmap[row],name,value])

out={
    'conventions':d['conventions'],
    'variable_order':d['variable_order'][:6],
    'source_order':d['source_order'],
    'block_orders':[len(part) for _,part in groups],
    'original_block_and_indices':groups,
    'monomials':[d['monomials'][i-1][:6] for i in evenrows],
    'gram_records':gram,'source_records':source,
}
(here/'contact_affine_parity_reduced.json').write_text(
    json.dumps(out,separators=(',',':'))+'\n',encoding='utf8')
report={
    'exact_projection_equivalence':'average all eight radical-sign conjugations',
    'blocks':len(groups),'gram_order':sum(out['block_orders']),
    'total_order_with_source_signs':9+sum(out['block_orders']),
    'gram_variables':sum(n*(n+1)//2 for n in out['block_orders']),
    'equations':len(evenrows),'gram_records':len(gram),
    'source_records':len(source),'discarded_cross_character_records':discarded,
    'no_change_to_admitted_source_cone':True,
}
(here/'parity_reduction_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
