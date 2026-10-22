"""Independent standard-library replay of the structured local source pencil.

Checks all coefficient identities and exact reconstruction of the admitted
rational source fixture. It is not a full-Q222 adequacy checker.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import json

HERE=Path(__file__).resolve().parent
pencil=json.loads((HERE/'six_root_structured_source_pencil.coo.json').read_text())
products=json.loads((HERE/'six_root_structured_product_map.json').read_text())
fixture=json.loads((HERE/'six_root_structured_fixture.coo.json').read_text())
blocks=pencil['blocks'];n=pencil['matrix_order'];m=pencil['auxiliary_count']
assert n==161 and m==666 and len(blocks)==19
assert max(x['order'] for x in blocks)==10
offsets=[x['offset'] for x in blocks]
assert sum(x['order'] for x in blocks)==n
position_to_col={tuple(pos):k for k,pos in enumerate(products['gram_positions'])}
expressions=[{} for _ in position_to_col]
for k,i,j,value in pencil['entries']:
    assert 0<=k<=7+m and 0<=i<=j<n
    owner=next(b for b in blocks if b['offset']<=i<b['offset']+b['order'])
    assert j<owner['offset']+owner['order']
    pos=(owner['index'],i-owner['offset'],j-owner['offset'])
    col=position_to_col[pos]
    expressions[col][k]=expressions[col].get(k,F(0))+F(value)
for row in products['coefficient_rows']:
    value=defaultdict(F)
    for col,mult in row['gram']:
        for k,v in expressions[col].items():value[k]+=F(mult)*v
    expected={1+k:F(v) for k,v in row['target']}
    assert {k:v for k,v in value.items() if v}==expected,row['monomial']

fixture_entries={tuple(record[:3]):F(record[3]) for record in fixture['entries']}
source={'d':F(1)}
source.update({k:F(v) for k,v in fixture['source_coefficients'].items()})
values={1+k:source[name] for k,name in enumerate(pencil['canonical_source_order'])}
for j,pos in enumerate(pencil['free_gram_positions']):values[8+j]=fixture_entries.get(tuple(pos),F(0))
values[0]=F(1)
for col,pos in enumerate(products['gram_positions']):
    actual=sum(v*values[k] for k,v in expressions[col].items())
    assert actual==fixture_entries.get(tuple(pos),F(0)),pos
receipt={'all180coefficient_equations_exact':True,'all666auxiliary_cancellations_exact':True,
         'admitted_fixture_reconstructed_by_actual_affine_pencil':True,
         'order':n,'auxiliaries':m,'blocks':19,'max_block_order':10,
         'full_Q222_reverse_inclusion_asserted':False}
(HERE/'six_root_independent_replay.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
