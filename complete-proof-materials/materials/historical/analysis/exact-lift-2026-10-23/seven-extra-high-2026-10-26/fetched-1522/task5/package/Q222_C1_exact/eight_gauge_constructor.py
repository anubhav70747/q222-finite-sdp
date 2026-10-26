#!/usr/bin/env python3
"""Exact *coefficient-matching* 7+5+1 affine C1 Gram, NOT a PSD solver.

Accepts exact rational native C1 parameters only where original source is defined.
Source-only normalization and eight gauge coordinates. Constructing a matrix
is not certification; PSD must be independently proved over the whole domain.
"""
from c1_small_map import *
from c1_native import k,f,g,a,t,NATIVE,CONTACT,Lc,Hc,Qc,R,T,M as GapM,J,P4,Pa,Pd,c3,c2,cb,cal,cd,D as De,be
INVARIANT_ROWS=[j for j,(i,h,x,y) in enumerate(KEYS) if (i+h+int(x in (1,2))+int(y in (1,2)))%2==0]
A42=M[INVARIANT_ROWS,BLKCOL]
PIV_COLS=list(A42.rref()[1]);PIV_ROWS=list(A42[:,PIV_COLS].T.rref()[1]);FIXED_INVERSE=A42[PIV_ROWS,PIV_COLS].inv()
GAUGES=A42.nullspace()
assert len(GAUGES)==8

def gram_native(params, gauges=(0,)*8):
    """Return rational 7x7, 5x5, scalar odd coefficient for exact params.
    No optimization, positivity assertion or sampling. Use `M`, `ROWS`, `ODD`
    to reconstruct the original feature polynomial identity.
    """
    if len(params)!=5 or len(gauges)!=8: raise ValueError('Need 5 native parameters and 8 gauge entries')
    vals=dict(zip((k,f,g,a,t),map(s.Rational,params)))
    strict_checks={name:s.cancel(expr.subs(vals)) for name,expr in (
      ('k',k),('f',f),('g',g),('a',a),('t',t),('Lc',Lc),('Hc',Hc),('Qc',Qc),
      ('R',R),('T',T),('M',GapM),('J',J),('1-tg2',1-t*g*g),
      ('E3',1-k*k*f*f+t*k*k*f*f*g*g),('P4',P4),('Pa',Pa),('Pd',Pd),
      ('c3',c3),('c2',c2),('cb',cb),('cal',cal),('cd',cd),('D_event',De),('beta',be))}
    if any(value<=0 for value in strict_checks.values()):
      raise ValueError('Source violates strict native conditions: '+str({x:y for x,y in strict_checks.items() if y<=0}))
    weights=[s.cancel(x.subs(vals)) for x in NATIVE]
    contact=[s.cancel(x.subs(vals)) for x in CONTACT]
    if any(x==0 for x in contact):raise ValueError('Contact state/angle must be nonzero')
    rhs=coefficients(sym_parent(weights))
    rhs=s.Matrix([rhs[j]*contact[0]**i*contact[1]**h*contact[2+x]*contact[2+y]
         for j,(i,h,x,y) in enumerate(KEYS)])
    r=rhs[INVARIANT_ROWS,:]
    solution=s.zeros(44,1);chosen=FIXED_INVERSE*r[PIV_ROWS,:]
    for p,i in enumerate(PIV_COLS):solution[i]=chosen[p]
    for bcol,v in zip(GAUGES,gauges):solution+=s.Rational(v)*bcol
    # The exact six stationarity identities guarantee the remaining rows.
    if A42*solution!=r:raise ValueError('Native stationarity/source premise failed')
    H=s.zeros(12)
    for j,col in enumerate(BLKCOL):
        if col==78:continue
        p,q=PAIR[col];H[p,q]=H[q,p]=solution[j]
    return H[:7,:7],H[7:,7:],solution[-1], weights, contact

if __name__=='__main__':
 print('C1 exact constructor available: gram_native((k,f,g,a,t), eight_gauge_values)')
 print('No witness search or PSD claim is performed.')
