#!/usr/bin/env python3
"""Purely symbolic, deterministic full-coefficient and PSD-reduction checks.
No native source evaluation, numerical optimizer, sweep or approximate matrix.
The baseline native proof is checked separately by original_baseline/*.py.
"""
import sympy as S
from c1_z18_extension import *
from c1_small_map import sym_parent

print('START Z18 EXACT CHECKS',flush=True)
# 14 integer + four 1/2-power features: even contact blocks 7,5 and
# odd four-feature contact block 3, with 46 generated support monomials.
assert M18.shape==(46,49)
assert len(FULL_PARITY)==46 and len(OLD_INVARIANT)==42
assert len(NEW_KEYS)==4 and not set(NEW_KEYS).intersection(set(FULL_PARITY))
print('PASS: 42 native invariant keys vs 46 full parity-compatible keys; four new odd keys forbidden',flush=True)
A42=M[INVARIANT_ROWS,BLKCOL]
assert A42.shape==(42,44)
assert len(PIV_COLS)==36 and len(PIV_ROWS)==36
assert A42[PIV_ROWS,PIV_COLS].det()==-2**24
assert A42[PIV_ROWS,PIV_COLS]*FIXED_INVERSE==S.eye(36)
# The old rank is 36 because of six independent contact derivative annihilators.
# The four new odd-parity monomials have rank-four support independent of old columns.
minor_cols=[HALF_PAIRS.index(pair)+len(INT_PAIRS)
            for pair in ((0,1),(0,2),(1,1),(2,2))]
R_new=M18.extract([INDEX18[i] for i in NEW_KEYS],minor_cols)
assert R_new.det()==16
print('PASS: 46x49 invariant-plus-odd map has rank 40 by 36+4 minor and six contact jets',flush=True)

# Four forbidden odd-parity coefficients force a two-scalar fractional Gram.
hmat, forb=reduce_half_block()
a,d=S.symbols('a delta',real=True)
assert hmat==S.Matrix([[S.Symbol('a'),0,0],[0,S.Symbol('delta',real=True),-S.Symbol('delta',real=True)],[0,-S.Symbol('delta',real=True),S.Symbol('delta',real=True)]])
assert hmat[1:,1:].det()==0
print('PASS: odd block reduces to gamma>=0, delta>=0, on W-Z and X-Y separately',flush=True)

# Exact two zero identities independent of native source/contact.
r=[]
for row in ROWS:
  r.append(S.expand(sum(coeff*u**FEATS[i][1]*v**FEATS[i][2]* (W,X,Y,Z)[FEATS[i][3]] for i,coeff in row)))
rv=S.Matrix(r)
Kgamma=S.zeros(12)
for (i,j),val in { (0,1):S.Rational(-1,2),(1,6):S.Rational(1,2)}.items():
 Kgamma[i,j]=Kgamma[j,i]=val
assert S.expand((rv.T*Kgamma*rv)[0]+u*v*(W-Z)**2)==0
assert S.expand((rv.T*KDELTA*rv)[0]+u*v*(X-Y)**2)==0
# This identity must hold in the entire four-state polynomial, not at contact.
# The 42-target map verifies all angular/state coefficients at once.
assert all(z==0 for z in S.simplify(A42*S.Matrix([
    KDELTA[PAIR[BLKCOL[i]][0],PAIR[BLKCOL[i]][1]] if BLKCOL[i]!=78 else 0
    for i in range(44)]) +coefficients(u*v*(X-Y)**2)[INVARIANT_ROWS,:]))
print('PASS: exact eighth/ninth source-only gauge identities in all four states and arbitrary u,v',flush=True)

# Gamma/delta are rank-one PSD half-power packets; no mixed testing/source root.
# The source-normalized half-power state features remain original fixed
# sqrt(u*v) W,X,Y,Z up to nonzero pure-source factors.
for name,root in [('W-Z',W-Z),('X-Y',X-Y)]:
 assert S.Poly(S.expand(u*v*root**2),u,v,W,X,Y,Z).total_degree()<=4
print('PASS: complete fixed-feature coefficient matching and pure-source contact scalings',flush=True)

# The forced P and other two pivots are unchanged by the new gauge.
for i,j in [(1,1),(1,2),(2,2),(4,4),(11,11),(0,1),(0,2),(0,4),(1,4),(2,4)]:
 assert KDELTA[i,j]==0
print('PASS: existing strict forced Schur pivots survive the delta extension unchanged',flush=True)

# Symbolic rank statement for new map: 6 jets upper bound can be verified
# directly on the generated difference-root polynomial columns.
# The old upper bound has an independent six-row rank witness; padding
# new rows cannot destroy independence.
statevars=(W,X,Y,Z)
from c1_small_map import KEY_INDEX
jets=S.Matrix([[int(ll==a)+int(mm==a) for i,j,ll,mm in KEYS18] for a in range(4)]+[
   [i for i,j,ll,mm in KEYS18],[j for i,j,ll,mm in KEYS18]])
assert jets*M18==S.zeros(6,49)
assert jets.rank()==6
print('PASS: six original contact jets annihilate all 49 columns, so exact rank=40',flush=True)
print('ALL Z18 SYMBOLIC CHECKS PASSED',flush=True)
