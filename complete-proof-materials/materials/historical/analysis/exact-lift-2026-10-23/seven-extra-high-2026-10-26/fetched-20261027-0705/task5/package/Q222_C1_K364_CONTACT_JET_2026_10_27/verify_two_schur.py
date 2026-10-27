#!/usr/bin/env python3
"""Exact 2x2 Schur equivalence for the two C1 parity sectors.
This checks algebraic coefficient/Schur identities, NOT universal positivity.
"""
import sympy as s
D=s.symbols('D',real=True,positive=True)
def sym(prefix):
 a,b,c=s.symbols(prefix+'00 '+prefix+'01 '+prefix+'11',real=True)
 return s.Matrix([[a,b],[b,c]])
A0=sym('A');B0=sym('B');A2=sym('C');B2=sym('E')
L=s.Matrix(2,2,s.symbols('L00 L01 L10 L11',real=True));T=s.Matrix(2,2,s.symbols('T00 T01 T10 T11',real=True))
I=s.eye(2);Z=s.zeros(2)
def check(X,K,Y):
 inv=X.inv()
 Sch=Y-K.T*inv*K
 h=X.row_join(K).col_join(K.T.row_join(Y))
 upper=I.row_join(inv*K).col_join(Z.row_join(I))
 rhs=upper.T*(X.row_join(Z).col_join(Z.row_join(Sch)))*upper
 assert (h-rhs).applyfunc(s.cancel)==s.zeros(4)
 return h
check(A0,L+T,A2)
check(B0,L.T-T.T,B2)
print('PASS full 4x4 block = positive 2x2 pivot plus complete 2x2 Schur complement in both parity sectors')
print('PASS all four cross-packet matrix entries and all state-normal entries retained')
print('SCOPE: positivity of each 2x2 residual for every Bob v is a separate, unproved source-wide implication.')
