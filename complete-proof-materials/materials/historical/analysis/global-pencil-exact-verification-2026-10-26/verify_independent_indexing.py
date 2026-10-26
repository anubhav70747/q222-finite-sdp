#!/usr/bin/env python3
"""Independent finite indexing checks. Imports no supplied generator/assembler.
Reads original recovered JSON only. No solver, samples, or large LMI expansion.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement
import json,time,sys
START=time.perf_counter()
ROOT=Path(r'C:/Users/anubh/OneDrive/Documents/ChatGPT/ManuscriptWriting/q222/analytical-paper-2026-10-22')
OUT=Path(__file__).resolve().parent
SRC=ROOT/'analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1819/task7/package'
d=json.loads((SRC/'K364_FULL128_FACTORED_SPARSE_PENCIL.json').read_text())
rows=json.loads((SRC/'K364_GLOBAL_ALL_3456_AFFINE_EQUALITIES.json').read_text())
keys=[tuple(x) for x in d['polynomial_keys']];ki={x:i for i,x in enumerate(keys)}
features=[tuple(x) for x in d['fixed_testing_state_features']]
N=len(features);K=len(keys);assert (N,K)==(364,3250)
assert len(ki)==K
stored_A=d['Gram_variable_upper_entries_map_A']
stored_Astar=d['Gram_dual_moment_upper_entries_map_A_star']
pairs=list(combinations_with_replacement(range(N),2))
assert len(pairs)==len(stored_A)==len(stored_Astar)==66430
seen_keys=set();diagonal=offdiagonal=0
for ordinal,(i,j) in enumerate(pairs):
    ai,bi,si=features[i];aj,bj,sj=features[j]
    derived_key=(ai+aj,bi+bj,min(si,sj),max(si,sj));idx=ki[derived_key]
    # Direct coefficient extraction: z_i z_j occurs once in z^T P z
    # for a diagonal variable and twice for one symmetric offdiagonal variable.
    quadratic_terms=[(i,j)] if i==j else [(i,j),(j,i)]
    poly_coefficient=len(quadratic_terms)
    # Independently expand Tr(M(y) P)=sum_ab M_ab(y) P_ba.
    # One stored symmetric P_ij variable occupies (i,j), and (j,i) if distinct.
    matrix_positions={(i,j),(j,i)}
    trace_coefficient=sum(1 for a,b in matrix_positions if {a,b}=={i,j})
    assert poly_coefficient==trace_coefficient
    assert stored_A[ordinal]=={'entry':[i,j],'coefficient_key':idx,'A_scalar_on_symmetric_variable':poly_coefficient}
    assert stored_Astar[ordinal]==[i,j,idx]
    seen_keys.add(idx)
    if i==j:diagonal+=1
    else:offdiagonal+=1
assert seen_keys==set(range(K))
print(f'PASS independent coefficient and trace extraction: {len(pairs)} upper pairs ({diagonal} diagonal, {offdiagonal} offdiagonal); all {K} keys',flush=True)

# Expected partition sizes are derived from dimensions and upper triangles.
def upper_pairs(n):return list(combinations_with_replacement(range(n),2))
packed={n:{pair:i for i,pair in enumerate(upper_pairs(n))} for n in (3,6)}
varspec=[('K_dual_coeffs_y',K),('K_orthant_r',9),('SG12_dual_upper_Y',len(packed[6])),('SG12_r',7),('SG13_dual_upper_Y1',len(packed[3])),('SG13_dual_upper_Y2',len(packed[3])),('SG13_r',7)]
blockspec=[('K_moment',N),('K_source_orthant',9),('SG12_dual_PSD',6),('SG12_orthant',7),('SG13_dual_PSD_1',3),('SG13_dual_PSD_2',3),('SG13_orthant',7)]
kindcounts={'K_adjoint':9,'SG12_target':7,'SG12_skew':3,'SG13_target':7,'SG13_xi':1}
vcopy=sum(n for _,n in varspec);bcopy=sum(n for _,n in blockspec);ecopy=sum(kindcounts.values())
assert (vcopy,bcopy,ecopy)==(3306,399,27)
copies=d['group_copy_offsets'];assert len(copies)==128
varend=8;blockend=0;eqend=0;expected_vars=[]
for g,copy in enumerate(copies):
    assert copy['group_id']==g and copy['equality_offset']==eqend and copy['equality_counts']==kindcounts
    v={};b={}
    for name,n in varspec:v[name]=[varend,n];varend+=n
    for name,n in blockspec:b[name]=[blockend,n];blockend+=n
    assert copy['variables']==v and copy['PSD_blocks']==b
    expected_vars.append(v);eqend+=ecopy
assert (varend,blockend,eqend)==(423176,51072,3456)
assert rows['variable_count']==varend and rows['q_first_8'] is True
assert d['counts']['total_scalar_variables']==varend
assert d['counts']['auxiliary_variables']==varend-8
assert d['counts']['PSD_total_order']==blockend+16
assert d['counts']['PSD_blocks_total']==len(copies)*len(blockspec)+1==897
assert d['counts']['affine_equality_count']==eqend
print(f'PASS all {len(copies)} copies: {len(copies)*len(varspec)} variable partitions, {len(copies)*len(blockspec)} block offsets, all equality offsets; event block completes order {blockend+16}',flush=True)

# No helper from the supplied assembly is used below. Reconstruct target
# E^T Q_g by reading the physical signed permutation directly.
E=[[Q(x) for x in row] for row in d['canonical_E_source_to_Bell_9x9']]
assert len(E)==9 and all(len(r)==9 for r in E)
faces=d['face_templates'];f12=faces['SG12'];f13=faces['SG13']
idx12=f12['source_indices'];idx13=f13['source_indices']
assert idx12==[0,2,4,5,6,7,8] and idx13==[0,2,3,4,5,6,8]
for face,indices in [(f12,idx12),(f13,idx13)]:
    assert face['source_injection_9x7']==[[int(i==j) for j in indices] for i in range(9)]
assert len(f12['source_matrix_bases_6'])==7 and len(f12['free_skew_matrix_bases_6'])==3
assert len(f13['source_matrix_bases_3_and_3'])==7 and len(f13['free_xi_matrix_bases_3_and_3'])==2
assert len(d['T_source_columns'])==9
for col in d['T_source_columns']:
    assert len({key for key,x in col})==len(col)
    assert all(0<=key<K and Q(x)!=0 for key,x in col)

# A basis upper entry v corresponds directly to two full-matrix positions
# off diagonal. This implements a literal symmetric trace coefficient.
def trace_terms(upper,n,start,sign=1):
    seen=set();result=[]
    for i,j,text in upper:
        assert 0<=i<=j<n and (i,j) not in seen
        seen.add((i,j));value=Q(text)
        full_positions={(i,j),(j,i)}
        result.append((start+packed[n][(i,j)],sign*value*len(full_positions)))
    return result
actual_rows=rows['rows'];assert len(actual_rows)==eqend
termcount=0;ordinal=0;rowkinds={name:0 for name in kindcounts}
def check_row(g,kind,constant,terms):
    global ordinal,termcount
    expected={}
    for index,value in terms:
        assert 0<=index<varend
        expected[index]=expected.get(index,Q(0))+Q(value)
    expected={index:value for index,value in expected.items() if value}
    actual=actual_rows[ordinal]
    assert actual['id']==ordinal and actual['group_id']==g and actual['kind']==kind
    assert Q(actual['constant'])==Q(constant)
    entries=actual['entries'];indices=[i for i,c in entries]
    assert indices==sorted(set(indices)) and all(Q(c)!=0 for i,c in entries)
    assert {i:Q(c) for i,c in entries}==expected,(ordinal,kind)
    termcount+=len(entries);rowkinds[kind]+=1;ordinal+=1
for g,copy in enumerate(copies):
    group=d['groups'][g];assert group['id']==g
    mapping=group['observable_map']
    assert len(mapping)==9 and sorted(j for j,s in mapping)==list(range(9))
    assert all(sign in (-1,1) for j,sign in mapping)
    target=[[Q(0) for c in range(9)] for r in range(9)]
    for j,(column,sign) in enumerate(mapping):
        for r in range(9):target[r][column]+=E[j][r]*sign
    assert [[Q(x) for x in row] for row in group['canonical_target_9x9']]==target
    gt=d['group_targets'][g]
    assert gt['group_id']==g and [[Q(x) for x in row] for row in gt['K9x9']]==target
    assert [[Q(x) for x in row] for row in gt['SG12_7x9']]==[target[i] for i in idx12]
    assert [[Q(x) for x in row] for row in gt['SG13_7x9']]==[target[i] for i in idx13]
    var=expected_vars[g]
    def target_terms(source_index):return [(c-1,value) for c,value in enumerate(target[source_index]) if c and value]
    for a in range(9):
        terms=target_terms(a)+[(var['K_dual_coeffs_y'][0]+key,-Q(coef)) for key,coef in d['T_source_columns'][a]]+[(var['K_orthant_r'][0]+a,-1)]
        check_row(g,'K_adjoint',target[a][0],terms)
    for a,source_index in enumerate(idx12):
        terms=target_terms(source_index)+trace_terms(f12['source_matrix_bases_6'][a],6,var['SG12_dual_upper_Y'][0],-1)+[(var['SG12_r'][0]+a,-1)]
        check_row(g,'SG12_target',target[source_index][0],terms)
    for a,basis in enumerate(f12['free_skew_matrix_bases_6']):
        check_row(g,'SG12_skew',0,trace_terms(basis,6,var['SG12_dual_upper_Y'][0]))
    for a,source_index in enumerate(idx13):
        terms=target_terms(source_index)+[(var['SG13_r'][0]+a,-1)]
        for k,name in enumerate(('SG13_dual_upper_Y1','SG13_dual_upper_Y2')):
            terms+=trace_terms(f13['source_matrix_bases_3_and_3'][a][k],3,var[name][0],-1)
        check_row(g,'SG13_target',target[source_index][0],terms)
    terms=[]
    for k,name in enumerate(('SG13_dual_upper_Y1','SG13_dual_upper_Y2')):
        terms+=trace_terms(f13['free_xi_matrix_bases_3_and_3'][k],3,var[name][0])
    check_row(g,'SG13_xi',0,terms)
assert ordinal==3456 and rowkinds=={name:128*n for name,n in kindcounts.items()}
print(f'PASS independently reconstructed all {ordinal} affine rows and {termcount} nonzero entries from base matrices, physical actions, trace contractions and checked offsets: {rowkinds}',flush=True)

# Verify event block and fixed copy maps without repeating 66430 local
# positions for every copy: local map was exhaustively checked above; all
# 128 translations have disjoint verified intervals and identical dimensions.
events=d['physical_event_inequalities'];stored_events=json.loads((SRC/'K364_ALL_EVENT_LMI_ENTRIES.json').read_text())
assert len(events)==len(stored_events)==16
for i,event in enumerate(events):
    x,y,a,b=event['x'],event['y'],event['a'],event['b']
    assert x in (0,1) and y in (0,1) and a in (-1,1) and b in (-1,1)
    row=[0]*9;row[0]=1;row[1+x]=a;row[3+y]=b;row[5+2*x+y]=a*b
    assert event['nonnegative_4p']==row
    expected={'global_row':blockend+i,'global_col':blockend+i,'constant':1,'q_variable_terms':[[j-1,c] for j,c in enumerate(row) if j and c]}
    assert stored_events[i]==expected
assert len({(e['x'],e['y'],e['a'],e['b']) for e in events})==16
elapsed=round(time.perf_counter()-START,6)
report={'result':'ALL TARGETED FINITE INDEXING CHECKS PASSED','upper_pairs':len(pairs),'diagonal_pairs':diagonal,'offdiagonal_pairs':offdiagonal,'coefficient_keys':K,'copies':len(copies),'copy_variable_partitions':len(copies)*len(varspec),'copy_PSD_blocks':len(copies)*len(blockspec),'PSD_blocks_including_events':897,'PSD_total_order':blockend+16,'scalar_variables':varend,'affine_rows':ordinal,'affine_nonzero_entries':termcount,'affine_kind_counts':rowkinds,'event_rows':16,'elapsed_seconds':elapsed,'imports_supplied_generation_assembly_code':False,'solvers':False,'expanded_all_128_moment_entries':False,'scope':'Exact indexing and coefficient consistency of the recovered finite pencil only. No universal source admission or reverse inclusion proof.'}
(OUT/'independent_indexing_results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print('PASS all 16 event entries; repeated moment copies covered by one complete local map and all checked disjoint translations',flush=True)
print(json.dumps(report,indent=2),flush=True)
