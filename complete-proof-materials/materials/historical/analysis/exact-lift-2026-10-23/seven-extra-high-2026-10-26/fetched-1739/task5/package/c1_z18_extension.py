#!/usr/bin/env python3
"""Exact parameter-free enlargement of native C1 Z16 by sqrt(uv) X,Y.
No SDP optimization, numerical sources, parameter sweeps or unjustified PSD assertions.
"""
import itertools
import sys
from pathlib import Path
import sympy as S
sys.path.insert(0,str(Path(__file__).resolve().parent/'original_baseline'))
from c1_small_map import (FEATS,ROWS,PAIR,KEYS,KEY_INDEX,M,
                          key,column,coefficients,u,v,W,X,Y,Z)
from eight_gauge_constructor import (INVARIANT_ROWS,BLKCOL,PIV_COLS,
                                     PIV_ROWS,FIXED_INVERSE,GAUGES)

HALF=[('sqrtuv'+name,S.Rational(1,2),S.Rational(1,2),idx)
      for idx,name in enumerate(('W','X','Y','Z'))]
HALF_ROWS=[[(i,S.S.One),(3,-S.S.One)] for i in range(3)]
HALF_PAIRS=list(itertools.combinations_with_replacement(range(3),2))
INT_PAIRS=[(i,j) for i,j in PAIR if (i<7)==(j<7)]

def half_column(i,j):
    accum={}
    for p,a in HALF_ROWS[i]:
        for q,b in HALF_ROWS[j]:
            k=key(HALF[p][1]+HALF[q][1],HALF[p][2]+HALF[q][2],
                  HALF[p][3],HALF[q][3])
            accum[k]=accum.get(k,0)+(1 if i==j else 2)*a*b
    return {k:x for k,x in accum.items() if x}

COLS18=[column(i,j) for i,j in INT_PAIRS]+[half_column(i,j) for i,j in HALF_PAIRS]
KEYS18=sorted({kk for col in COLS18 for kk in col})
INDEX18={key:i for i,key in enumerate(KEYS18)}
M18=S.MutableSparseMatrix(len(KEYS18),len(COLS18),
   {(INDEX18[k],j):val for j,col in enumerate(COLS18) for k,val in col.items()})
OLD_INVARIANT=[kk for kk in KEYS if (kk[0]+kk[1]+int(kk[2] in (1,2))
                                        +int(kk[3] in (1,2)))%2==0]
assert len(KEYS18)==46 and len(OLD_INVARIANT)==42
NEW_KEYS=sorted(set(KEYS18)-set(OLD_INVARIANT))
assert NEW_KEYS==[(1,1,0,1),(1,1,0,2),(1,1,1,3),(1,1,2,3)]
FULL_PARITY=sorted({key(i,j,l,m) for i in range(3) for j in range(3)
                    for l in range(4) for m in range(l,4)
                    if (i+j+int(l in (1,2))+int(m in (1,2)))%2==0})
MISSING_PARITY=sorted(set(FULL_PARITY)-set(OLD_INVARIANT))
assert len(FULL_PARITY)==46 and MISSING_PARITY==[(0,0,0,0),(0,0,0,3),(2,0,1,1),(2,0,1,2)]

def reduce_half_block():
    a,b,c,d,e,f=S.symbols('a b c d e f')
    # odd difference rows h0=sqrtuv(W-Z), h1=sqrtuv(X-Z), h2=sqrtuv(Y-Z)
    coeffs=[a,b,c,d,e,f]
    hmat=S.Matrix([[a,b,c],[b,d,e],[c,e,f]])
    oddcoeff=S.zeros(len(KEYS18),1)
    for z,coef in zip(HALF_PAIRS,coeffs):
       for k,v in half_column(*z).items(): oddcoeff[INDEX18[k]]+=coef*v
    forb=S.Matrix([oddcoeff[INDEX18[k]] for k in NEW_KEYS])
    out=S.solve(list(forb),[b,c,d,f],dict=True)
    # Solving exactly gives b=c=0, d=f=-e. Reparameterize delta=-e.
    assert len(out)==1
    reduced=hmat.subs(out[0]); delta=S.symbols('delta',real=True)
    reduced=S.simplify(reduced.subs(e,-delta))
    assert reduced==S.Matrix([[a,0,0],[0,delta,-delta],[0,-delta,delta]])
    return reduced,forb

# Fixed ninth gauge in first seven integer contact differences.
A42=M[INVARIANT_ROWS,BLKCOL]
KDELTA=S.zeros(12)
for (i,j),val in { (2,3):S.Rational(-1,2),
                   (2,5):S.Rational(1,2),
                   (3,4):S.Rational(1,2),
                   (4,5):S.Rational(-1,2)}.items():
    KDELTA[i,j]=KDELTA[j,i]=val


def gram18_oriented(params, gauges=(0,)*9):
    """Exact native coefficient matching for rational/real-algebraic inputs.
    Uses the fixed rational 36-pivot inverse and full native SD1 source.
    Gauges: s0..s6, gamma, delta; PSD must be checked separately.
    Source field and testing feature list remain independently defined.
    """
    from c1_native import (k,f,g,a,t,NATIVE,CONTACT,Lc,Hc,Qc,R,T,M as GapM,J,
                            P4,Pa,Pd,c3,c2,cb,cal,cd,D as De,be)
    from c1_small_map import sym_parent
    if len(params)!=5 or len(gauges)!=9:
        raise ValueError('Five native parameters and nine source-only gauges required')
    kv=S.sympify(params[0])
    if (kv.is_positive is False or (1-kv).is_nonnegative is False):
        raise ValueError('Source outside the oriented native domain 0<k<=1')
    subs=dict(zip((k,f,g,a,t),map(S.sympify,params)))
    checks={name:S.cancel(expr.subs(subs)) for name,expr in (
        ('k',k),('f',f),('g',g),('a',a),('t',t),('1-k',1-k),
        ('Lc',Lc),('Hc',Hc),('Qc',Qc),('M',GapM),('J',J),
        ('1-tg2',1-t*g*g),('E3',1-k*k*f*f+t*k*k*f*f*g*g),
        ('P4',P4),('Pa',Pa),('Pd',Pd),('c2',c2),('cb',cb),
        ('c3',c3),('cd',cd),('cal',cal),('De',De),('beta',be))}
    for name,val in checks.items():
        if name=='1-k':
            if val.is_nonnegative is not True:
                raise ValueError('Orientation not certified: '+str(val))
        elif val.is_positive is not True:
            raise ValueError('Strict source sign cannot be certified: '+name+'='+str(val))
    weights=[S.cancel(x.subs(subs)) for x in NATIVE]
    contact=[S.cancel(x.subs(subs)) for x in CONTACT]
    R=coefficients(sym_parent(weights))
    R=S.Matrix([R[j]*contact[0]**i*contact[1]**h*contact[2+x]*contact[2+y]
                for j,(i,h,x,y) in enumerate(KEYS)])
    rhs=R[INVARIANT_ROWS,:]
    sol=S.zeros(44,1)
    piv=FIXED_INVERSE*rhs[PIV_ROWS,:]
    for j,x in zip(PIV_COLS,piv):sol[j]=x
    for j,gauge in enumerate(gauges[:8]):sol+=S.sympify(gauge)*GAUGES[j]
    if any(S.cancel(z)!=0 for z in A42*sol-rhs):
        raise ValueError('Unexpected failure of exact native coefficient identity')
    H=S.zeros(12)
    for j,col in enumerate(BLKCOL):
        if col==78:continue
        p,q=PAIR[col]; H[p,q]=H[q,p]=sol[j]
    delta=S.sympify(gauges[8])
    H+=delta*KDELTA
    return H[:7,:7],H[7:,7:],sol[-1],delta,weights,contact


def residual_schur(H7,H5):
    """Exact 4x4+4x4 Schur matrices. Caller first proves forced pivots positive.
    All source divisions are by P7, H7[4,4], H5[4,4], which are strict
    on the genuine strict C1 native domain. NO tests/source contacts inverted.
    """
    rowsP=[1,2]; rowsT=[0,3,4,5,6]; rowsR=[0,3,5,6]
    P=H7.extract(rowsP,rowsP)
    T=H7.extract(rowsT,rowsT)-H7.extract(rowsT,rowsP)*P.inv()*H7.extract(rowsP,rowsT)
    idx4=2
    sch7=T.extract([0,1,3,4],[0,1,3,4])-T.extract([0,1,3,4],[idx4])*T.extract([idx4],[0,1,3,4])/T[idx4,idx4]
    sch5=H5[:4,:4]-H5[:4,4]*H5[4,:4]/H5[4,4]
    return sch7,sch5,P,T[idx4,idx4],H5[4,4]


def exact_psd(A):
    """Exact rational/algebraic PSD test including zero eigenvalues.
    For a symmetric real matrix PSD is equivalent to ALL principal minors >=0.
    No approximate eigenvalues or numerical solvers used.
    """
    for n in range(1,A.rows+1):
        for rows in itertools.combinations(range(A.rows),n):
            minor=S.cancel(A.extract(rows,rows).det())
            if minor.is_nonnegative is not True:
                if not minor.is_number:
                    raise ValueError('Cannot certify symbolic sign: '+str(minor))
                if minor<0:return False
    return True
