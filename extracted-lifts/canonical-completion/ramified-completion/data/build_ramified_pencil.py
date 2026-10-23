"""Compile the explicit ramified half-line affine pencil of order 49.

The coordinate is v=16s^2, s>=0. Variables: nine source coordinates and
148 auxiliaries. Blocks: 9 scalar signs, 12,8,12,8. All matrix entries
are rational and affine; no coefficient equations remain.

This is a sufficient canonical cone. Universal adequacy is NOT asserted.
Standard-library Python only; run python3 build_ramified_pencil.py.
"""
from fractions import Fraction as F
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
names=['d','lambda','n','m','C','A','D','alpha','beta']

def lin(*terms):
    out={}
    for a,p in terms:
        for i,x in p.items():
            out[i]=out.get(i,F(0))+F(a)*x
    return {i:x for i,x in out.items() if x}

def var(name):
    i=len(names);names.append(name);return {i:F(1)}

def zmat(n,m=None):
    return [[{} for j in range(n if m is None else m)] for i in range(n)]

def symvar(prefix,n=4):
    A=zmat(n)
    for i in range(n):
        for j in range(i,n):A[i][j]=A[j][i]=var(f'{prefix}_{i}{j}')
    return A

def fullvar(prefix,n=4):
    return [[var(f'{prefix}_{i}{j}') for j in range(n)] for i in range(n)]

def skewvar(prefix,n=4):
    A=zmat(n)
    for i in range(n):
        for j in range(i+1,n):
            A[i][j]=var(f'{prefix}_{i}{j}');A[j][i]=lin((-1,A[i][j]))
    return A

def setc(A,k,i,j,*terms):
    A[k][i][j]=A[k][j][i]=lin(*terms)

c=[{i:F(1)} for i in range(9)]
d,lam,n,m,C,AA,D,alpha,beta=c
P,N,E,K=([zmat(2) for k in range(5)] for j in range(4))
setc(P,4,0,0,(1024,m));setc(P,2,0,1,(-64,m))
setc(P,0,1,1,(4,m),(4,beta),(2,d));setc(P,4,1,1,(1024,beta),(512,d))
setc(N,0,0,0,(4,C),(4,D),(2,d));setc(N,4,0,0,(1024,C))
setc(N,2,0,1,(64,D),(32,d));setc(N,4,1,1,(1024,D),(512,d))
setc(E,0,0,0,(4,AA));setc(E,4,0,0,(1024,AA),(1024,alpha),(512,d))
setc(E,2,0,1,(-64,alpha),(-32,d));setc(E,0,1,1,(4,alpha));setc(E,4,1,1,(-512,d))
setc(K,0,0,0,(4,lam));setc(K,2,0,1,(64,lam))
setc(K,0,1,1,(4,n));setc(K,4,1,1,(1024,lam),(1024,n))
J=[fullvar(f'J{k}',2) for k in range(5)]
Qs=[]
for side in ['G','H']:
    Q=[zmat(4) for k in range(5)]
    for k in range(5):
        for i in range(2):
            for j in range(2):
                extra=K if side=='G' else E
                Q[k][i][j]=lin((1,P[k][i][j]),(1,extra[k][i][j]))
                Q[k][i+2][j+2]=lin((1,N[k][i][j]),(1,extra[k][i][j]))
                cross=J[k][i][j] if side=='G' else lin((1,K[k][i][j]),(-1,E[k][i][j]),(-1,J[k][i][j]))
                Q[k][i][j+2]=Q[k][j+2][i]=cross
    Qs.append(Q)

blocks=[]
for side,Q in zip(['G','H'],Qs):
    B00=symvar(f'{side}_B00');B11=symvar(f'{side}_B11')
    B01=fullvar(f'{side}_B01');A02=fullvar(f'{side}_A02')
    K01=skewvar(f'{side}_K01');K12=skewvar(f'{side}_K12')
    A=zmat(12);B=zmat(8)
    for i in range(4):
        for j in range(4):
            A[i][j]=Q[0][i][j];A[i+8][j+8]=Q[4][i][j]
            A[i+4][j+4]=lin((1,Q[2][i][j]),(-1,A02[i][j]),(-1,A02[j][i]),(-1,B01[i][j]),(-1,B01[j][i]))
            a01=lin((F(1,2),Q[1][i][j]),(-F(1,2),B00[i][j]),(1,K01[i][j]))
            a12=lin((F(1,2),Q[3][i][j]),(-F(1,2),B11[i][j]),(1,K12[i][j]))
            A[i][j+4]=A[j+4][i]=a01
            A[i+4][j+8]=A[j+8][i+4]=a12
            A[i][j+8]=A[j+8][i]=A02[i][j]
            B[i][j]=B00[i][j];B[i+4][j+4]=B11[i][j]
            B[i][j+4]=B[j+4][i]=B01[i][j]
    # Exact symbolic coefficient check: 2 * 5 * 10 = 100 entries including J.
    for k in range(5):
        for i in range(4):
            for j in range(4):
                terms=[(1,A[4*a+i][4*b+j]) for a in range(3) for b in range(3) if a+b==k]
                terms += [(1,B[4*a+i][4*b+j]) for a in range(2) for b in range(2) if a+b+1==k]
                assert lin(*terms)==Q[k][i][j]
    blocks.extend([A,B])
assert len(names)==157
L=zmat(49)
for i in range(9):L[i][i]=c[i]
for offset,A in zip([9,21,29,41],blocks):
    for i,row in enumerate(A):
        for j,x in enumerate(row):L[offset+i][offset+j]=x
assert all(L[i][j]==L[j][i] for i in range(49) for j in range(49))
entries=[[] for name in names]
for i in range(49):
    for j in range(i,49):
        for k,x in L[i][j].items():entries[k].append([i,j,str(x)])
out={
 'description':'Ramified quadratic-selector enlargement: explicit affine sufficient canonical cone; global adequacy not claimed.',
 'order':49,'source_dimension':9,'auxiliary_dimension':148,
 'source_names':names[:9], 'auxiliary_names':names[9:],
 'blocks':[{'name':names[i],'offset':i,'order':1} for i in range(9)]+[
    {'name':n,'offset':o,'order':z} for n,o,z in zip(['A_G','B_G','A_H','B_H'],[9,21,29,41],[12,8,12,8])],
 'coordinate':'v=16*s^2, s>=0',
 'matrix_convention':'Entries list upper-triangular [row,column,rational] with zero-based indices; reflect off-diagonal entries.',
 'pencil':'sum(c_i*source_matrices[i])+sum(a_j*auxiliary_matrices[j]) is positive semidefinite',
 'source_matrices':entries[:9], 'auxiliary_matrices':entries[9:],
 'feature_blocks':{'A':['I4','s*I4','s^2*I4'],'B':['sqrt(s)*I4','s^(3/2)*I4']},
 'product_exponents_A':[[i+j for j in range(3)] for i in range(3)],
 'product_exponents_B':[[1+i+j for j in range(2)] for i in range(2)],
 'constraints_before_elimination':80,'independent_Gram_entries':228
}
p=HERE/'ramified_order49_pencil.json';p.write_text(json.dumps(out,indent=2)+'\n')
print('Built',p.name,': order49, 9 sources,148 auxiliaries,',sum(map(len,entries)),'upper-triangular nonzero entries.')
print('All polynomial product coefficients checked exactly.')
