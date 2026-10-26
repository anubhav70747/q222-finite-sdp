#!/usr/bin/env python3
"""Exact, source-independent algebraic Schur congruence and forced pivots.
Uses no optimizers, sampling, or assumptions that the residual blocks are PSD.
"""
import sympy as S
from c1_z18_extension import KDELTA
p,q,r = S.symbols('p q r',real=True)
P=S.Matrix([[p,q],[q,r]])
B=S.Matrix(2,5, lambda i,j: S.Symbol(f'b{i}{j}'))
C=S.Matrix(5,5,lambda i,j:S.Symbol(f'c{min(i,j)}{max(i,j)}'))
Full=P.row_join(B).col_join(B.T.row_join(C))
T=C-B.T*P.inv()*B
E=S.eye(7); E[:2,2:]=P.inv()*B
Diag=S.diag(P,T)
assert all(S.cancel(x)==0 for x in (E.T*Diag*E-Full))
print('PASS: exact H7 2x2 strict Schur congruence',flush=True)
# second stage, delete the 5-block's forced strictly positive entry at position 2
T2=S.Matrix(5,5,lambda i,j:S.Symbol(f't{min(i,j)}{max(i,j)}'))
perm=[2,0,1,3,4]
Tp=T2.extract(perm,perm);f=Tp[0,0];R=Tp[0,1:];Tail=Tp[1:,1:]-Tp[1:,:1]*R/f
Sunit=S.eye(5);Sunit[0,1:]=R/f
assert all(S.cancel(x)==0 for x in (Sunit.T*S.diag(S.Matrix([[f]]),Tail)*Sunit-Tp))
print('PASS: exact H7 positive diagonal pivot and 4x4 residual congruence',flush=True)
D=S.Symbol('d11',positive=True);R=S.Matrix([[S.Symbol(f'b{i}') for i in range(4)]])
C=S.Matrix(4,4,lambda i,j:S.Symbol(f'v{min(i,j)}{max(i,j)}'))
Full=S.Matrix([[D]]).row_join(R).col_join(R.T.row_join(C))
Tail=C-R.T*R/D
E=S.eye(5);E[0,1:]=R/D
assert all(S.cancel(x)==0 for x in E.T*S.diag(S.Matrix([[D]]),Tail)*E-Full)
print('PASS: exact H5 last forced positive pivot and 4x4 residual congruence',flush=True)
# The delta gauge touches none of these strict pivots or the forced gamma Schur scalar
for i,j in [(1,1),(1,2),(2,2),(4,4),(11,11),(0,1),(0,2),(0,4),(1,4),(2,4)]:
 assert KDELTA[i,j]==0
print('PASS: delta gauge preserves all source-strict forced pivots and gamma bound',flush=True)
print('No universal PSD conclusion is asserted for either 4x4 tail.',flush=True)
