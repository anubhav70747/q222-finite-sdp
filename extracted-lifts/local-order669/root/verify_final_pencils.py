#!/usr/bin/env python3
"""Read-only exact replay of the final local Q222 source-cone matrices.

Python standard library only. No generators or archived verifier programs are
imported. Rebuilds all polynomial products and the full canonical parent,
checks coefficient maps, free-coordinate parametrizations, and the final sum.
Writes only FINAL_EXACT_REPLAY.json next to this script.
This does not assert a full-Q222 lift or global source coverage.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from bisect import bisect_right
import hashlib
import json

HERE=Path(__file__).resolve().parent
WORK=HERE.parent
UF=WORK/'uniform_tf1'; SIX=WORK/'candidate_boundary'
CANON=['d','m','C','A','D','alpha','beta']
READ=[]
def load(path):
    READ.append(path)
    return json.loads(path.read_text())
def add(d,e,k=F(1)):
    for x,c in e.items():
        d[x]=d.get(x,F(0))+k*c
        if not d[x]:del d[x]
def ep(x,y):return tuple(a+b for a,b in zip(x,y))
def times(p,q):
    out={}
    for e,c in p.items():
        for f,d in q.items():add(out,{ep(e,f):c*d})
    return out
def plus(*ps):
    out={}
    for p in ps:add(out,p)
    return out
def scale(p,c):return {e:F(c)*v for e,v in p.items() if c*v}
def square(p):return times(p,p)
def unit(n,i,k=1):
    e=[0]*n;e[i]=k;return {tuple(e):F(1)}
def p4(u,w,W,X,Y,Z,n):
    one={(0,)*n:F(1)}
    Du=plus(one,square(u));Dw=plus(one,square(w));uw=times(u,w)
    r1=plus(W,scale(times(w,X),-1),times(u,Y),scale(times(uw,Z),-1))
    r2=plus(times(u,W),scale(Y,-1))
    r3=plus(times(uw,W),times(u,X),scale(times(w,Y),-1),scale(Z,-1))
    r4=plus(Y,scale(times(w,Z),-1))
    r5=plus(X,times(u,Z))
    q3=square(r3);q4=times(Du,square(r4));q5=times(Dw,square(r5))
    d=plus(scale(q3,2),scale(q4,2),scale(q5,2),scale(times(times(Du,Dw),square(Z)),-2))
    return [d,scale(square(r1),4),scale(times(Dw,square(r2)),4),
            scale(times(times(Du,Dw),square(Y)),4),scale(q3,4),scale(q4,4),scale(q5,4)]
def gate_reduce(e):
    # Independent finite binomial expansion of xi^(2k),eta^(2l).
    from math import comb
    k,l=e[6]//2,e[7]//2
    base=list(e);base[6]%=2;base[7]%=2
    out={}
    for i in range(k+1):
        for j in range(l+1):
            v=base.copy();v[0]+=4*i;v[1]+=4*j
            add(out,{tuple(v):F((-1)**(i+j)*comb(k,i)*comb(l,j))})
    return out

def positions(orders):
    return [(b,i,j) for b,n in enumerate(orders) for i in range(n) for j in range(i,n)]
def rank(rows):
    # Exact fraction-free-style sparse Gaussian elimination, no saved pivots.
    piv={}
    for row in rows:
        v=row.copy()
        while v:
            p=min(v)
            if p not in piv:
                q=v[p];piv[p]={j:c/q for j,c in v.items()};break
            add(v,piv[p],-v[p])
    return len(piv)
def pencil_entries(data):
    ex=defaultdict(dict)
    if 'entries' in data:
        offs=[b['offset'] for b in data['blocks']]
        for k,i,j,c in data['entries']:
            assert k>0,'Unexpected constant matrix entry'
            bi,bj=bisect_right(offs,i)-1,bisect_right(offs,j)-1
            assert bi==bj
            add(ex[bi,i-offs[bi],j-offs[bi]],{k-1:F(c)})
    else:
        assert not data.get('constant_entries',[])
        for b,i,j,k,c in data['coo']:add(ex[b,i,j],{k:F(c)})
    orders=data['block_orders']
    for b,i,j in ex:assert 0<=b<len(orders) and 0<=i<=j<orders[b]
    return dict(ex)
def coefficient_equations(rows,gram_positions,ex,target_rows,source_map=None):
    for row,target in zip(rows,target_rows):
        result={}
        for g,c in row.items():add(result,ex.get(gram_positions[g],{}),c)
        expected={}
        for v,c in target.items():
            if source_map is None:add(expected,{v:c})
            else:add(expected,source_map[v],c)
        assert result==expected,('coefficient equation mismatch',result,expected)

def target_rows(monomials,targets):
    ri={e:i for i,e in enumerate(monomials)}
    rows=[{} for _ in monomials]
    for c,p in enumerate(targets):
        for e,v in p.items():
            assert e in ri,('parent monomial absent from products',e)
            rows[ri[e]][c]=v
    return rows

def verify_parametrization(rows,gram_positions,ex,aux_count,explicit_free=None):
    rr=rank(rows)
    assert rr==len(rows)
    assert len(gram_positions)-rr==aux_count
    byunit={}
    for pos,expr in ex.items():
        if len(expr)==1:
            k,v=next(iter(expr.items()))
            if k>=7 and v==1:byunit.setdefault(k,pos)
    assert set(byunit)==set(range(7,7+aux_count))
    selected=[byunit[7+i] for i in range(aux_count)]
    assert len(set(selected))==aux_count
    if explicit_free is not None:
        assert len(explicit_free)==aux_count
        for j,pos in enumerate(explicit_free):assert ex[tuple(pos)]=={7+j:F(1)}
    # The coefficient equations imply all auxiliary columns lie in ker(A).
    # These unit coordinates make them independent; the checked rank proves
    # their count equals nullity, hence they form the entire kernel.
    return rr,len(selected)

def verify_signed(label,sign):
    fd=load(UF/f'C1_{label}_FULL_FIXED_FEATURES.json')
    md=load(UF/f'C1_{label}_FULL_COEFFICIENT_MAP.coo.json')
    pd=load(UF/f'C1_{label}_FULL_SOURCE_PENCIL.coo.json')
    assert pd['source_order']==CANON
    orders=pd['block_orders']; pos=positions(orders)
    assert [list(x) for x in pos]==md['column_gram_entries']
    mons=[tuple(e) for e in md['row_monomials']];ri={e:i for i,e in enumerate(mons)}
    rows=[{} for _ in mons]
    for g,(b,i,j) in enumerate(pos):
        features=fd['blocks'][b]['features']
        assert len(features)==orders[b]
        prod=gate_reduce(ep(features[i],features[j]))
        for e,c in prod.items():
            assert e in ri
            rows[ri[e]][g]=c*(1 if i==j else 2)
    exported=[{} for _ in mons]
    for i,j,c in md['coo']:add(exported[i],{j:F(c)})
    assert rows==exported,'Exported coefficient map differs from fixed functions'
    de=F(1,4096)
    args=[scale(unit(8,0,8),de**2),scale(unit(8,1,4),de),unit(8,2,4),
          unit(8,4),unit(8,5),scale(unit(8,3,4),sign)]
    targets=target_rows(mons,p4(*args,8))
    target_export=[{} for _ in mons]
    for r,c,v in md['target_coo']:add(target_export[r],{c:F(v)})
    assert targets==target_export,'Exported target is not the complete P4 parent'
    ex=pencil_entries(pd)
    coefficient_equations(rows,pos,ex,targets)
    rr,nfree=verify_parametrization(rows,pos,ex,pd['counts']['auxiliary_dimension'])
    report={'status':'passed','blocks':len(orders),'matrix_order':sum(orders),
            'rebuilt_product_columns':len(pos),'exact_coefficient_equations':len(rows),
            'full_parent_target_rebuilt':True,'rank':rr,'free_coordinates_verified':nfree}
    return {'data':pd,'positions':pos,'rows':rows,'targets':targets,'entries':ex,'report':report}

def verify_six():
    pd=load(SIX/'six_root_structured_source_pencil.coo.json')
    md=load(SIX/'six_root_structured_product_map.json')
    gd=load(SIX/'six_root_general_formulas.json')
    assert pd['canonical_source_order']==CANON
    assert F(gd['testing_radius'])==F(1,4096)
    orders=pd['block_orders'];pos=positions(orders)
    assert [list(x) for x in pos]==md['gram_positions']
    mons=[tuple(r['monomial']) for r in md['coefficient_rows']];ri={e:i for i,e in enumerate(mons)}
    rows=[{} for _ in mons]
    # Printed SR2 order, checked against the source formula export.
    alphas=[(1,0),(0,1),(1,1),(2,0),(0,2),(2,1),(3,0)]
    assert list(gd['coefficient_matrices'])==['(0, 0)']+[str(a) for a in alphas]
    for g,(b,i,j) in enumerate(pos):
        features=pd['per_block_feature_monomials'][b]
        e=ep(features[i],features[j]);prod={e:F(1)}
        if 1<=b<=14:
            a,c=alphas[(b-1)//2];sign=1 if b%2 else -1
            shift=(a,c,0,0,0,0)
            add(prod,{ep(e,shift):F(sign*4096**(a+c))})
        for e,c in prod.items():
            assert e in ri
            rows[ri[e]][g]=c*(1 if i==j else 2)
    exported=[{g:F(c) for g,c in r['gram']} for r in md['coefficient_rows']]
    assert rows==exported,'Structured six-root products do not match their map'
    targets=target_rows(mons,p4(*(unit(6,i) for i in range(6)),6))
    exported_targets=[{c:F(v) for c,v in r['target']} for r in md['coefficient_rows']]
    assert targets==exported_targets,'Structured target is not the complete P4 parent'
    ex=pencil_entries(pd)
    coefficient_equations(rows,pos,ex,targets)
    rr,nfree=verify_parametrization(rows,pos,ex,pd['auxiliary_count'],pd['free_gram_positions'])
    report={'status':'passed','blocks':len(orders),'matrix_order':sum(orders),
            'rebuilt_product_columns':len(pos),'exact_coefficient_equations':len(rows),
            'full_parent_target_rebuilt':True,'rank':rr,'free_coordinates_verified':nfree}
    return {'data':pd,'positions':pos,'rows':rows,'targets':targets,'entries':ex,'report':report}

def map_entries(ex,blockoff,variable_map):
    out={}
    for (b,i,j),expr in ex.items():
        v={}
        for k,c in expr.items():add(v,variable_map(k),c)
        if v:out[b+blockoff,i,j]=v
    return out

def main():
    plus=verify_signed('POS',1);minus=verify_signed('NEG',-1);six=verify_six()
    both=load(UF/'C1_BOTH_SIGNS_LOCAL_SOURCE_PENCIL.coo.json')
    bo=pencil_entries(both)
    bp=len(plus['data']['block_orders']);ap=plus['report']['free_coordinates_verified']
    expected=dict(plus['entries'])
    expected.update(map_entries(minus['entries'],bp,lambda k:{k if k<7 else k+ap:F(1)}))
    assert bo==expected,'The two-sign intersection pencil is not its exact block concatenation'
    final=load(HERE/'C1_UNIFORM_LOCAL_ORDER669.coo.json');fe=pencil_entries(final)
    nblock_both=len(both['block_orders']);naux_both=both['auxiliary_dimension']
    expected=map_entries(bo,0,lambda k:{k+7 if k<7 else 14+k-7:F(1)})
    def six_vars(k):
        return {k:F(1),7+k:F(-1)} if k<7 else {14+naux_both+k-7:F(1)}
    expected.update(map_entries(six['entries'],nblock_both,six_vars))
    assert fe==expected,'The combined pencil is not the exact source sum split'
    assert final['source_order']==CANON
    assert len(final['auxiliary_order'])==7+naux_both+six['report']['free_coordinates_verified']==2694
    assert len(final['block_orders'])==120 and sum(final['block_orders'])==669
    assert max(final['block_orders'])==12
    # Directly replay the final combined matrix through all three product maps.
    for component,offset,mapping in [
        (plus,0,{i:{7+i:F(1)} for i in range(7)}),
        (minus,bp,{i:{7+i:F(1)} for i in range(7)}),
        (six,nblock_both,{i:{i:F(1),7+i:F(-1)} for i in range(7)})]:
        gpos=[(b+offset,i,j) for b,i,j in component['positions']]
        coefficient_equations(component['rows'],gpos,fe,component['targets'],mapping)
    unique_reads=list(dict.fromkeys(READ))
    result={'status':'passed','verifier':'Python standard-library Fraction arithmetic; no generator imports',
      'scope':'Exact equations and parametrizations for the evaluated local source cones. Not a proof of full-Q222 reverse inclusion.',
      'components':{'POS':plus['report'],'NEG':minus['report'],'SIX_STRUCTURED':six['report']},
      'combined':{'source_coordinates':7,'auxiliaries':2694,'blocks':120,'matrix_order':669,'largest_block':12,
                  'two_sign_intersection_exact':True,'source_sum_split_exact':True,
                  'direct_final_parent_equations_verified':sum(c['report']['exact_coefficient_equations'] for c in (plus,minus,six)),
                  'all_component_free_directions_verified':sum(c['report']['free_coordinates_verified'] for c in (plus,minus,six)),
                  'nonzero_upper_affine_records':len(final['coo'])},
      'inputs':[{'path':str(p.relative_to(WORK)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in unique_reads],
      'mutation_policy':'Pencil, feature and map files were read only; this verifier writes only FINAL_EXACT_REPLAY.json.'}
    (HERE/'FINAL_EXACT_REPLAY.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='inputs'},indent=2))

if __name__=='__main__':main()
