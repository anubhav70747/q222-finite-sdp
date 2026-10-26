#!/usr/bin/env python3
"""Exact verifier for fixed C2 Z6 pruned-face primal/dual witnesses.

A certificate is a JSON file with {"source": ["a","b","c","t","r"],
 "gram_blocks": [matrix of rational strings, ...]}; Gram orders are 51/35/35/29.
A dual witness has {"source": [...],"dual_moments": [412 rational strings]}.
No SDP, floating point, testing denominator, or absent binary artifact.
Run python verify_exact_certificates.py --primal witness.json
 or python verify_exact_certificates.py --dual moments.json.

The tool verifies supplied witnesses; it does NOT produce one for an
arbitrary strict C2 contact.
"""
import argparse,json,sys
import sympy as sp
from c2_fixed_map import (a,b,c,t,r,n,m,C,A,D,al,coefficient_map,
                          native_weights,target_template,is_strict_native,
                          contact_annihilators)

def rat(s):
    z=sp.Rational(str(s));assert z.is_Rational
    return z

def rational_psd(M):
    """Full exact Schur/LDL test, including singular pivots."""
    H=sp.Matrix(M)
    if H!=H.T: return False
    while H.rows:
        idx=next((i for i in range(H.rows) if H[i,i]>0),None)
        if idx is None:
            return H==sp.zeros(H.rows)
        if any(H[i,i]<0 for i in range(H.rows)):return False
        order=[idx]+[i for i in range(H.rows) if i!=idx]
        H=H.extract(order,order)
        pivot=H[0,0];v=H[1:,0]; B=H[1:,1:]
        H=B-v*v.T/pivot
    return True

def target_exact(source,M):
    vals=dict(zip((a,b,c,t,r),map(rat,source)))
    wt,_=native_weights(); wt={k:sp.cancel(v.subs(vals)) for k,v in wt.items()}
    uc=vals[c]; vc=vals[a]*vals[c]
    state=[vals[b]/vals[a],-vals[c]*vals[r],vals[c]*(vals[b]/vals[a]-vals[a]*vals[t]),sp.Integer(1)]
    target=[sp.Rational(0) for _ in M['row_keys']]
    ridx={key:i for i,key in enumerate(M['row_keys'])}
    for key,expr in target_template().items():
        p,q,l,k=key
        target[ridx[key]]=sp.expand(expr).subs(wt, simultaneous=True)*uc**(p//6)*vc**(q//6)*state[l]*state[k]
    return [sp.cancel(x) for x in target]

def source_ok(source):
    assert len(source)==5 and is_strict_native(source), 'not a strict native C2 point'

def verify_primal(obj,M):
    source_ok(obj['source']); gs=obj['gram_blocks']
    assert len(gs)==4
    sizes=[len(b)-1 for b in M['blocks']]
    lhs=[sp.Rational(0)]*len(M['row_keys'])
    colid=0
    for blk,size in enumerate(sizes):
        G=sp.Matrix([[rat(z) for z in row] for row in gs[blk]])
        assert G.shape==(size,size)
        assert rational_psd(G), ('gram block not PSD',blk)
        for j in range(size):
            for i in range(j+1):
                g=G[i,j]
                for row,val in M['sparse_columns'][colid].items(): lhs[row]+=g*val
                colid+=1
    target=target_exact(obj['source'],M)
    assert lhs==target,'exact full-parent coefficient mismatch'
    print('PASS EXACT NATIVE C2 FULL PARENT + 4 PSD CONTACT-REDUCED GRAMS')

def verify_dual(obj,M):
    source_ok(obj['source'])
    yy=[rat(q) for q in obj['dual_moments']]
    assert len(yy)==len(M['row_keys'])
    nn=contact_annihilators(M['row_keys'])
    assert all(sum(row[k]*yy[k] for k in range(len(yy)))==0 for row in nn), 'wrong six-functional dual gauge'
    colid=0; tr=sp.Rational(0)
    for blk in M['blocks']:
        size=len(blk)-1
        B=sp.zeros(size)
        for j in range(size):
            for i in range(j+1):
                col=M['sparse_columns'][colid]
                value=sum(v*yy[row] for row,v in col.items())/(1 if i==j else 2)
                B[i,j]=B[j,i]=value
                colid+=1
        assert rational_psd(B), 'dual contact moment block not PSD'
        tr+=sp.trace(B)
    assert tr==1,('dual trace normalization failed',tr)
    pair=sum(v*w for v,w in zip(yy,target_exact(obj['source'],M)))
    if pair<0:
        print('EXACT: strict exclusion of L6 at given native source; pairing:',pair)
    elif pair==0:
        print('EXACT: no strictly positive reduced Gram at given source (membership undecided); pairing 0')
    else:
        print('EXACT: positive dual pairing only, not a universal admission certificate; pairing:',pair)

if __name__=='__main__':
    p=argparse.ArgumentParser();grp=p.add_mutually_exclusive_group(required=True)
    grp.add_argument('--primal');grp.add_argument('--dual'); z=p.parse_args()
    with open(z.primal or z.dual) as f:ob=json.load(f)
    M=coefficient_map()
    verify_primal(ob,M) if z.primal else verify_dual(ob,M)
