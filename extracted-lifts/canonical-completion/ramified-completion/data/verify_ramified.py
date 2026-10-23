#!/usr/bin/env python3
"""Standard-library exact verifier; no numerical solver or CAS is used.

It checks the rational ramified certificate for c*, constructs the explicitly
stated d=0 event addition for cbar, checks every polynomial coefficient, and
optionally exports cbar's actual rational Gram matrices.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent

def zero(n,m=None):return [[F(0) for _ in range(n if m is None else m)] for _ in range(n)]
def eye(n):
    a=zero(n)
    for i in range(n):a[i][i]=F(1)
    return a
def transpose(a):return list(map(list,zip(*a)))
def add(a,b):return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def mul(a,b):
    bt=transpose(b)
    return [[sum((x*y for x,y in zip(r,c)),F(0)) for c in bt] for r in a]
def outer(a):return [[x*y for y in a] for x in a]
def scaled(a,w):return [[w*x for x in r] for r in a]
def decode(a):return [[F(x) for x in r] for r in a]
def encode(a):return [[str(x) for x in r] for r in a]
def congruence(f,x):return mul(mul(f,x),transpose(f))
def block(a,b,c,d):return [r+s for r,s in zip(a,b)]+[r+s for r,s in zip(c,d)]
def sub(a,rows,cols):return [[a[i][j] for j in cols] for i in rows]

def positive_ldl(a):
    n=len(a);l=eye(n);d=[]
    assert a==transpose(a)
    for j in range(n):
        p=a[j][j]-sum((l[j][k]**2*d[k] for k in range(j)),F(0))
        assert p>F(1,2**18),('pivot below certified lower bound',j,p)
        d.append(p)
        for i in range(j+1,n):
            l[i][j]=(a[i][j]-sum((l[i][k]*l[j][k]*d[k] for k in range(j)),F(0)))/p
    dd=zero(n)
    for i,p in enumerate(d):dd[i][i]=p
    assert congruence(l,dd)==a
    return d

def poly_coeff(q,r,k):
    a=zero(4)
    for i in range(3):
        for j in range(3):
            if i+j==k:a=add(a,sub(q,range(4*i,4*i+4),range(4*j,4*j+4)))
    for i in range(2):
        for j in range(2):
            if i+j+1==k:a=add(a,sub(r,range(4*i,4*i+4),range(4*j,4*j+4)))
    return a

def source_coefficients(c):
    d,lam,n,m,C,A,D,al,be=c
    P=[zero(2) for _ in range(5)];N=[zero(2) for _ in range(5)]
    E=[zero(2) for _ in range(5)];K=[zero(2) for _ in range(5)]
    P[0][1][1]=4*m+4*be+2*d
    P[2][0][1]=P[2][1][0]=-64*m
    P[4][0][0]=1024*m;P[4][1][1]=256*(4*be+2*d)
    N[0][0][0]=4*C+4*D+2*d
    N[2][0][1]=N[2][1][0]=16*(4*D+2*d)
    N[4][0][0]=1024*C;N[4][1][1]=256*(4*D+2*d)
    E[0][0][0]=4*A;E[0][1][1]=4*al
    E[2][0][1]=E[2][1][0]=-16*(4*al+2*d)
    E[4][0][0]=256*(4*A+4*al+2*d);E[4][1][1]=-512*d
    K[0][0][0]=4*lam;K[0][1][1]=4*n
    K[2][0][1]=K[2][1][0]=64*lam
    K[4][1][1]=1024*(lam+n)
    return P,N,E,K

def check_coefficients(c,grams,jcoeff=None):
    P,N,E,K=source_coefficients(c)
    gg=[poly_coeff(grams[0],grams[1],k) for k in range(5)]
    hh=[poly_coeff(grams[2],grams[3],k) for k in range(5)]
    jj=[sub(g,range(2),range(2,4)) for g in gg]
    if jcoeff is not None:assert jj==jcoeff
    for k in range(5):
        assert gg[k]==block(add(P[k],K[k]),jj[k],transpose(jj[k]),add(N[k],K[k]))
        cross=add(add(K[k],scaled(E[k],-1)),scaled(jj[k],-1))
        assert hh[k]==block(add(P[k],E[k]),cross,transpose(cross),add(N[k],E[k]))
    return jj

def event_additions(delta):
    """Nonnegative rational weights and 12-vectors for QG and QH."""
    d,lam,n,m,C,A,D,al,be=delta
    assert d==0 and all(x>=0 for x in delta)
    def vector(a,b):return list(map(F,a))+[F(0)]*4+list(map(F,b))
    shared=[
        ('m',4*m,vector([0,-1,0,0],[16,0,0,0])),
        ('beta_constant',4*be,vector([0,1,0,0],[0,0,0,0])),
        ('beta_v',4*be,vector([0,0,0,0],[0,16,0,0])),
        ('C_constant',4*C,vector([0,0,1,0],[0,0,0,0])),
        ('C_v',4*C,vector([0,0,0,0],[0,0,16,0])),
        ('D',4*D,vector([0,0,1,0],[0,0,0,16]))]
    gs=shared+[
        ('lambda',4*lam,vector([1,0,1,0],[0,16,0,16])),
        ('n_constant',4*n,vector([0,1,0,1],[0,0,0,0])),
        ('n_v',4*n,vector([0,0,0,0],[0,16,0,16]))]
    hs=shared+[
        ('A_constant',4*A,vector([1,0,-1,0],[0,0,0,0])),
        ('A_v',4*A,vector([0,0,0,0],[16,0,-16,0])),
        ('alpha',4*al,vector([0,-1,0,1],[16,0,-16,0]))]
    return gs,hs

def run(export=False):
    cert=json.loads((HERE/'ramified_exact.json').read_text())
    c=list(map(F,cert['source']))
    faces=list(map(decode,cert['face_matrices']))
    small=list(map(decode,cert['reduced_Gram_blocks']))
    for i,a in enumerate(small):
        assert positive_ldl(a)==list(map(F,cert['LDL_pivots'][i]))
    full=[congruence(f,x) for f,x in zip(faces,small)]
    # Undo the two explicitly stored diagonal congruences.
    dg=[F(cert['D_G'][i][i]) for i in range(4)]
    dh=[F(cert['D_H'][i][i]) for i in range(4)]
    for index,diag in enumerate([dg,dg,dh,dh]):
        full[index]=[[a/(diag[i%4]*diag[j%4]) for j,a in enumerate(r)] for i,r in enumerate(full[index])]
    j=check_coefficients(c,full,list(map(decode,cert['J_coefficients'])))
    # The repair uses genuinely odd powers of s, so it does not merely
    # reproduce a quadratic polynomial in v.
    assert any(x for k in [1,3] for r in j[k] for x in r)

    cbar=[F(1),F(1,2**34),F(1,2**14),F(1,2**24),F(1,576),F(161),F(1,2),F(76),F(1,250)]
    delta=[b-a for a,b in zip(c,cbar)]
    gs,hs=event_additions(delta)
    repaired=[[list(r) for r in a] for a in full]
    for index,terms in [(0,gs),(2,hs)]:
        for _,w,z in terms:
            assert w>=0
            repaired[index]=add(repaired[index],scaled(outer(z),w))
    jb=check_coefficients(cbar,repaired)
    kd=source_coefficients(delta)[3]
    assert jb==[add(a,b) for a,b in zip(j,kd)]

    if export:
        def terms(t):return [{'event':label,'weight':str(w),'vector':list(map(str,z))} for label,w,z in t]
        result={
            'description':'Actual rational quartic half-line Gram matrices for the all-positive separated source. v=16s^2.',
            'source_order':cert['source_order'],'source':list(map(str,cbar)),
            'source_delta_from_cstar':list(map(str,delta)),
            'block_order':['Q_G','R_G','Q_H','R_H'],
            'block_sizes':[12,8,12,8],
            'Gram_matrices':list(map(encode,repaired)),
            'J_coefficients':list(map(encode,jb)),
            'positive_event_additions_to_Q_G':terms(gs),
            'positive_event_additions_to_Q_H':terms(hs),
            'PSD_proof':'ramified_exact.json positive reduced blocks, their stored congruences, and the nonnegative rank-one event additions above.'}
        (HERE/'ramified_cbar_exact.json').write_text(json.dumps(result,indent=2))
    target=HERE/'ramified_cbar_exact.json'
    if target.exists():
        supplied=json.loads(target.read_text())
        assert list(map(decode,supplied['Gram_matrices']))==repaired
        assert list(map(F,supplied['source']))==cbar
    print('PASS: four rational LDL factorizations; all old and all-positive source coefficients; every event-addition weight; actual exported Gram matrices.')
    return full,repaired

if __name__=='__main__':run('--export' in sys.argv)
