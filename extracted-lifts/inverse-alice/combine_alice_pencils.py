"""Form the exact intersection of the inner- and inverse-Alice source cones.

Usage: python combine_alice_pencils.py INNER.json OUTER.json OUTPUT.json
The visible seven canonical source coefficients are shared; auxiliary
coordinates are independent. This is a local source cone, not full Q222.
"""
from pathlib import Path
import argparse
import hashlib
import json


def combine(inner_path, outer_path, output):
    inner=json.loads(inner_path.read_text())
    outer=json.loads(outer_path.read_text())
    assert inner['source_order']==outer['source_order']==['d','m','C','A','D','alpha','beta']
    inner_aux=len(inner['auxiliary_order'])
    outer_aux=outer['counts']['auxiliary_dimension']
    assert not inner['constant_entries'] and not outer['constant_entries']
    ni=len(inner['block_orders'])
    joined=list(inner['coo'])
    shifted=[[b+ni,i,j,v if v<7 else v+inner_aux,c]
             for b,i,j,v,c in outer['coo']]
    joined.extend(shifted)
    # Literal reverse extraction verifies that neither component has changed.
    assert joined[:len(inner['coo'])]==inner['coo']
    restored=[[b-ni,i,j,v if v<7 else v-inner_aux,c] for b,i,j,v,c in shifted]
    assert restored==outer['coo']
    orders=inner['block_orders']+outer['block_orders']
    result={
      'schema':'q222.intersection_of_local_source_cones.v1',
      'scope':'Exact pencil for the intersection of the specified inner- and inverse-Alice certificate cones; not full Q222.',
      'equivalence':'c belongs to K_inner intersect K_inverse iff blockdiag(L_inner(c,y),L_inverse(c,z)) is PSD for independent auxiliaries y,z.',
      'source_order':inner['source_order'],
      'auxiliary_dimension':inner_aux+outer_aux,
      'block_orders':orders,'total_matrix_order':sum(orders),
      'constant_entries':[],'coo':joined,
      'convention':outer['convention'],
      'native_admission':{
        'source':'d=1; 0<z<=1; 0<x,sigma<=2^-78/(1+z); native strict gaps; D,beta>=0',
        'testing':'Every projective Alice angle u, |w|<=2^-36, all state vectors',
        'proof':'inverse_alice_all_z_certificate.tex; inner proof in ../local-order669/proof/'},
      'input_sha256':{inner_path.name:hashlib.sha256(inner_path.read_bytes()).hexdigest(),
                      outer_path.name:hashlib.sha256(outer_path.read_bytes()).hexdigest()},
      'counts':{'blocks':len(orders),'total_order':sum(orders),'largest_block':max(orders),
                'auxiliaries':inner_aux+outer_aux,'affine_entries':len(joined)}}
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['counts'],indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('inner',type=Path)
    ap.add_argument('outer',type=Path)
    ap.add_argument('output',type=Path)
    args=ap.parse_args()
    combine(args.inner,args.outer,args.output)
