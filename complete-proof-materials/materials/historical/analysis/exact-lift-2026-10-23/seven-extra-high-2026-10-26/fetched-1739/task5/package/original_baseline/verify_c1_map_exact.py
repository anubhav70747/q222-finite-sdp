#!/usr/bin/env python3
"""Exact symbolic coefficient-space theorem for native C1 14+2 fixed features.
Does not claim any PSD completion. Deterministic QQ linear algebra only.
"""
from c1_small_map import *
from pathlib import Path

def verify():
    z=s.symbols('d m C A D alpha beta')
    dd,m,C,A,D,al,be=z
    raw=coefficients(sym_parent((dd,0,0,m,C,A,D,al,be)))
    x,y,w0,x0,y0,z0=s.symbols('uc vc wc xc yc zc',nonzero=True)
    scales=[x**p*y**q*(w0,x0,y0,z0)[a]*(w0,x0,y0,z0)[b] for p,q,a,b in KEYS]
    f=s.matrix_multiply_elementwise(raw,s.Matrix(scales))
    left=M.T.nullspace()
    universal=[i for i,lf in enumerate(left) if s.cancel((lf.T*f)[0])==0]
    assert len(universal)==6,universal
    # at normalized contact U=V=W=X=Y=Z=1, exact stationarity:
    # 4 state derivatives and two testing derivatives
    derivatives=[]
    for k in range(4):
        derivatives.append(s.Matrix([[int(a==k)+int(b==k) for i,j,a,b in KEYS]]))
    derivatives.extend([s.Matrix([[i for i,j,a,b in KEYS]]),s.Matrix([[j for i,j,a,b in KEYS]])])
    E=s.Matrix.vstack(*(left[i].T for i in universal),*derivatives)
    assert E.rank()==12
    assert E*M==s.zeros(12,M.cols)
    # fixed rational right inverse on image
    pivcols=list(M.rref()[1]);pivrows=list(M[:,pivcols].T.rref()[1]); assert len(pivcols)==len(pivrows)==66
    Mp=M[pivrows,pivcols]; assert Mp.det()!=0
    inv=Mp.inv()
    selector=s.zeros(M.cols,66)
    for r,c in enumerate(pivcols):
        for j in range(66):selector[c,j]=inv[r,j]
    assert M*selector*M[pivrows,:]==M
    # Restrict by the combined native angle/sign symmetry: 42 invariant coefficient rows.
    invrows=[i for i,(p,q,a,b) in enumerate(KEYS)
             if (p+q+int(a in (1,2))+int(b in (1,2)))%2==0]
    A42=M[invrows,BLKCOL]
    pc42=list(A42.rref()[1]);pr42=list(A42[:,pc42].T.rref()[1]);
    assert A42.shape==(42,44) and A42.rank()==36
    assert A42[pr42,pc42].det()==-2**24
    jets=s.Matrix([[int(a==j)+int(b==j) for i,k,a,b in (KEYS[r] for r in invrows)] for j in range(4)] +
          [[KEYS[r][0] for r in invrows],[KEYS[r][1] for r in invrows]])
    assert jets.rank()==6 and jets*A42==s.zeros(6,44)
    assert jets[:,[0,1,3,4,6,14]].det()==4
    assert len(A42.nullspace())==8
    print('PASS: invariant 42x44 reduced C1 map rank 36, 8 exact gauge null vectors',flush=True)
    print('PASS: fixed 36x36 integer pivot determinant -2^24; 6x6 contact-jet minor 4',flush=True)
    print('Reduced pivot columns',pc42,flush=True)
    print('Reduced pivot rows',pr42,flush=True)
    print('PASS: 16 candidate features, 12 contact roots, rank 66, 6 universal zero equations',universal,flush=True)
    print('PASS: all 6 C1 state/angle jets + 6 universal equations span cokernel exactly',flush=True)
    print('PASS: lexicographic fixed rational 66 pivot right inverse on complete coefficient image',flush=True)
    print('Pivot columns:',pivcols,flush=True)
    print('Pivot rows:',pivrows,flush=True)
    return universal, pivcols,pivrows
if __name__=='__main__':verify()
