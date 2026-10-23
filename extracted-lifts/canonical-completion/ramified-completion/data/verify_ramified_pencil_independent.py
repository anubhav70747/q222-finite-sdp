"""Independent standard-library check of the exported affine order-49 pencil.

Checks every source/auxiliary coefficient against the native P,N,E,K
definitions after v=16*s^2, without using the pencil compiler.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("pencil",nargs="?",type=Path,
                    default=Path(__file__).with_name("ramified_order49_pencil.json"))
data=json.loads(parser.parse_args().pencil.read_text())

def zero(n): return [[F(0) for _ in range(n)] for _ in range(n)]
def add(A,B): return [[x+y for x,y in zip(r,s)] for r,s in zip(A,B)]
def sub(A,B): return [[x-y for x,y in zip(r,s)] for r,s in zip(A,B)]
def corner(A,i,j): return [r[j:j+2] for r in A[i:i+2]]
def gram_coeff(A,B,k):
    return [[sum(A[4*p+i][4*q+j] for p in range(3) for q in range(3) if p+q==k)
             +sum(B[4*p+i][4*q+j] for p in range(2) for q in range(2) if p+q+1==k)
             for j in range(4)] for i in range(4)]

def source_coeff(c,k):
    if k%2: return [zero(2) for _ in range(4)]
    d,lam,n,m,C,A,D,al,be=c
    P=[[[0,0],[0,4*m+4*be+2*d]],
       [[0,-4*m],[-4*m,0]],[[4*m,0],[0,4*be+2*d]]]
    N=[[[4*C+4*D+2*d,0],[0,0]],
       [[0,4*D+2*d],[4*D+2*d,0]],[[4*C,0],[0,4*D+2*d]]]
    E=[[[4*A,0],[0,4*al]],
       [[0,-4*al-2*d],[-4*al-2*d,0]],[[4*A+4*al+2*d,0],[0,-2*d]]]
    K=[[[4*lam,0],[0,4*n]],
       [[0,4*lam],[4*lam,0]],[[0,0],[0,4*lam+4*n]]]
    j=k//2
    return [[[16**j*x for x in row] for row in M[j]] for M in (P,N,E,K)]

assert data["order"]==49 and data["source_dimension"]==9
assert data["auxiliary_dimension"]==148
assert data["source_names"]==["d","lambda","n","m","C","A","D","alpha","beta"]
matrices=data["source_matrices"]+data["auxiliary_matrices"]
assert len(matrices)==157
for index,sparse in enumerate(matrices):
    L=zero(49)
    for i,j,x in sparse:
        assert 0<=i<=j<49
        L[i][j]+=F(x)
        if i!=j:L[j][i]+=F(x)
    c=[F(index==i) for i in range(9)]
    for i in range(9):
        for j in range(49): assert L[i][j]==(c[i] if i==j else 0)
    blocks=[]
    for offset,size in [(9,12),(21,8),(29,12),(41,8)]:
        blocks.append([r[offset:offset+size] for r in L[offset:offset+size]])
        for i in range(offset,offset+size):
            for j in range(49):
                if not offset<=j<offset+size:assert L[i][j]==0
    for k in range(5):
        G=gram_coeff(*blocks[:2],k)
        H=gram_coeff(*blocks[2:],k)
        P,N,E,K=source_coeff(c,k)
        assert corner(G,0,0)==add(P,K)
        assert corner(G,2,2)==add(N,K)
        assert corner(H,0,0)==add(P,E)
        assert corner(H,2,2)==add(N,E)
        assert add(corner(G,0,2),corner(H,0,2))==sub(K,E)
print("PASS: all 157 exported affine matrices have the prescribed source signs, block layout, and five polynomial coefficients; every auxiliary cancels from the source equations.")
