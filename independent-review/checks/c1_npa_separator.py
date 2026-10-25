"""Independent rational universal-word reconstruction of the C1 separator.

Uses only Python's standard library. No supplied affine Gram is imported.
"""
from fractions import Fraction as F
from pathlib import Path
import json

def reduce_letters(seq):
    out=[]
    for letter in seq:
        if out and out[-1]==letter: out.pop()
        else: out.append(letter)
    return tuple(out)

def multiply(x,y):
    return (reduce_letters(x[0]+y[0]),reduce_letters(x[1]+y[1]))

def adjoint(x):
    return (tuple(reversed(x[0])),tuple(reversed(x[1])))

def product(p,q):
    out={}
    for x,c in p.items():
        for y,d in q.items():
            w=multiply(adjoint(x),y)
            out[w]=out.get(w,0)+c*d
    return {w:c for w,c in out.items() if c}

def rank(A):
    A=[[F(x) for x in row] for row in A]
    nr,nc=len(A),len(A[0]); row=0
    for col in range(nc):
        pivot=next((i for i in range(row,nr) if A[i][col]),None)
        if pivot is None: continue
        A[row],A[pivot]=A[pivot],A[row]
        q=A[row][col]; A[row]=[x/q for x in A[row]]
        for i in range(row+1,nr):
            q=A[i][col]
            if q: A[i]=[x-q*y for x,y in zip(A[i],A[row])]
        row+=1
        if row==nr: break
    return row

def solve(A,b):
    nr,nc=len(A),len(A[0])
    R=[[F(x) for x in row]+[F(y)] for row,y in zip(A,b)]
    row=0; pivots=[]
    for col in range(nc):
        pivot=next((i for i in range(row,nr) if R[i][col]),None)
        if pivot is None: continue
        R[row],R[pivot]=R[pivot],R[row]
        q=R[row][col]; R[row]=[x/q for x in R[row]]
        for i in range(nr):
            q=R[i][col]
            if i!=row and q: R[i]=[x-q*y for x,y in zip(R[i],R[row])]
        pivots.append(col); row+=1
        if row==nr: break
    assert not any(all(x==0 for x in line[:-1]) and line[-1] for line in R)
    x=[F(0)]*nc
    for i,col in enumerate(pivots): x[col]=R[i][-1]
    assert all(sum(F(z)*w for z,w in zip(line,x))==F(rhs) for line,rhs in zip(A,b))
    return x

I=((),()); A0=((0,),()); A1=((1,),()); B0=((),(0,)); B1=((),(1,))
A0B0=multiply(A0,B0); A0B1=multiply(A0,B1)
A1B0=multiply(A1,B0); A1B1=multiply(A1,B1)
roots=[{A0:1,I:-1},{A1B0:1,I:-1},{A1B1:1,I:-1},
       {B0:1,A1:-1},{B1:1,A1:-1},{A0B0:1,A1:-1},{A0B1:1,A1:-1}]
upper=[(i,j) for i in range(7) for j in range(i,7)]
columns=[]
for i,j in upper:
    p=product(roots[i],roots[j])
    if i!=j:
        for w,c in product(roots[j],roots[i]).items(): p[w]=p.get(w,0)+c
    columns.append({w:c for w,c in p.items() if c})
words=sorted(set().union(*(p.keys() for p in columns)))
assert len(words)==25
M=[[p.get(w,0) for p in columns] for w in words]
assert rank(M)==15

# Source coefficients of the nine allowed Bell words, not contact evaluations.
# Active source order (d,m,C,A,D,alpha,beta); lambda=n=0.
target={
 I:[1,1,1,1,1,1,1], A0:[0,0,0,-1,0,-1,0],
 A1:[0,1,-1,0,-1,0,1], B0:[0,0,1,1,0,0,-1],
 B1:[0,-1,0,0,1,-1,0], A0B0:[F(-1,2),0,0,-1,0,0,0],
 A0B1:[F(1,2),0,0,0,0,1,0],
 A1B0:[F(-1,2),0,-1,0,0,0,-1],
 A1B1:[F(-1,2),-1,0,0,-1,0,0]}
T=[target.get(w,[0]*7) for w in words]
assert rank([x+y for x,y in zip(M,T)])==15

K=[
 [389,135,2,203,-21,-98,126],
 [135,158,-12,42,0,-161,-84],
 [2,-12,11,0,7,14,28],
 [203,42,0,158,-12,23,142],
 [-21,0,7,-12,11,9,9],
 [-98,-161,14,23,9,277,240],
 [126,-84,28,142,9,240,396]]
L=[[F(i==j) for j in range(7)] for i in range(7)]
diagonal=[]; minors=[]; leading=F(1)
for i in range(7):
    pivot=F(K[i][i])-sum(L[i][j]**2*diagonal[j] for j in range(i))
    assert pivot>0
    diagonal.append(pivot); leading*=pivot; minors.append(leading)
    for j in range(i+1,7):
        L[j][i]=(F(K[j][i])-sum(L[j][h]*L[i][h]*diagonal[h] for h in range(i)))/pivot
assert all(sum(L[i][h]*diagonal[h]*L[j][h] for h in range(7))==K[i][j]
           for i in range(7) for j in range(7))

# Independently find a universal word-coefficient functional whose adjoint is K.
trace_vector=[K[i][j]*(1 if i==j else 2) for i,j in upper]
MT=[list(row) for row in zip(*M)]
functional=solve(MT,trace_vector)
index={w:i for i,w in enumerate(words)}
functional=[(functional[index[w]]+functional[index[adjoint(w)]])/2 for w in words]
assert all(sum(F(c)*functional[h] for h,c in enumerate(row))==rhs
           for row,rhs in zip(MT,trace_vector))
source_trace=[sum(functional[i]*T[i][j] for i in range(25)) for j in range(7)]
assert source_trace==[-61,9,58,44,2,121,100]

cint=[F(1),F(20,29),F(12,29),F(13,58),F(51,58),F(3,29),F(3,58)]
cbdy=[F(1),F(5,6),F(3,10),F(1,3),F(2),F(1,13),F(0)]
violations=[sum(c*b for c,b in zip(source_trace,row)) for row in (cint,cbdy)]
assert violations==[-F(43,29),-F(3169,390)]

# Full nine-word reduction: the seven roots span the two deterministic kernels.
W=[I,A0,A1,B0,B1,A0B0,A0B1,A1B0,A1B1]
C=[[root.get(w,0) for root in roots] for w in W]
def evaluate(word,point):
    result=1
    for x in word[0]: result*=point[x]
    for x in word[1]: result*=point[2+x]
    return result
evaluations=[[evaluate(w,p) for w in W] for p in [(1,1,1,1),(1,-1,-1,-1)]]
assert rank(C)==7 and rank(evaluations)==2
assert all(sum(row[i]*C[i][j] for i in range(9))==0 for row in evaluations for j in range(7))
assert all(sum(F(evaluate(w,point))*T[i][j] for i,w in enumerate(words))==0
           for point in [(1,1,1,1),(1,-1,-1,-1)] for j in range(7))

def label(w):
    return ''.join('A'+str(i) for i in w[0])+''.join('B'+str(i) for i in w[1]) or 'I'
data={
 'word_count':25,'symmetric_Gram_variables':28,'universal_coefficient_map_rank':15,
 'free_real_Gram_parameters':13,'word_order':list(map(label,words)),
 'Gram_upper_triangle_order':[list(x) for x in upper],
 'universal_coefficient_map':M,'source_coefficient_map':[[str(v) for v in row] for row in T],
 'K':K,'K_LDL_positive_pivots':list(map(str,diagonal)),
 'K_leading_principal_minors':list(map(str,minors)),
 'dual_universal_word_functional':{label(w):str(v) for w,v in zip(words,functional) if v},
 'exact_source_trace_order':['d','m','C','A','D','alpha','beta'],
 'exact_source_trace_coefficients':list(map(str,source_trace)),
 'interior_source_violation':str(violations[0]),'beta_zero_source_violation':str(violations[1]),
 'nine_word_root_basis_rank':7,'deterministic_evaluation_rank':2,
 'scope':'Excludes both stated C1 slacks from the full nine-word 1+AB SOS cone; physical validity and extremality audited separately.',
 'no_numerical_optimizer_or_parameter_sweep':True}
Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in data.items() if k not in ['universal_coefficient_map','source_coefficient_map']},indent=2))
