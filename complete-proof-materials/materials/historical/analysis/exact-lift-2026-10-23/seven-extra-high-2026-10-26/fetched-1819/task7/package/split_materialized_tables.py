#!/usr/bin/env python3
"""Extract optional individually inspectable cards from the exact master data."""
import json, gzip
from pathlib import Path
HERE=Path(__file__).resolve().parent
with open(HERE/'K364_FULL128_FACTORED_SPARSE_PENCIL.json') as f:d=json.load(f)
exports={
 'K364_T_3250x9_SPARSE.json':{
    'matrix_order':'3250x9','number_of_nnz':d['counts']['T_nonzero_terms'],
    'source_order':d['source_order'],'keys':d['polynomial_keys'],
    'sparse_columns':d['T_source_columns']},
 'K364_PHYSICAL_128_ACTIONS_TARGETS.json':{
    'coordinate_order':d['observable_homogeneous_order'],'source_order':d['source_order'],
    'canonical_E':d['canonical_E_source_to_Bell_9x9'],
    'groups':d['groups'],'targets':d['group_targets'],
    'PR_sectors':d['PR_sectors']},
 'K364_SG12_SG13_SOURCE_AND_DUAL_BASES.json':{
    'source_order':d['source_order'],'blocks':d['face_templates']},
 'K364_BLOCK_AND_VARIABLE_OFFSETS.json':{
    'counts':d['counts'],'copy_offsets':d['group_copy_offsets'],
    'entry_convention':d['source_dual_conventions']},
 'K364_ASTAR_66430_UPPER.json':{
    'keys':d['polynomial_keys'],'features':d['fixed_testing_state_features'],
    'Astar_upper_entries':d['Gram_dual_moment_upper_entries_map_A_star'],
    'A_upper_convention':'1 diagonal, 2 off-diagonal for symmetric independent P_ij; Astar[y]_ij=y_key'},
}
for name,record in exports.items():
    path=HERE/name
    with open(path,'w') as f:json.dump(record,f,separators=(',',':'));f.write('\n')
    if name.endswith('ASTAR_66430_UPPER.json'):
        with gzip.open(str(path)+'.gz','wt') as f:json.dump(record,f,separators=(',',':'))
    print(name, path.stat().st_size)
