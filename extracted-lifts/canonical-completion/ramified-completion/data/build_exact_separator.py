"""Build the exact order-33 separator from six fixed rational 4x4 matrices.

Standard-library Python only. No numerical optimization, numerical rounding,
or external data are used. Run: python3 build_exact_separator.py
"""
from fractions import Fraction as F
from pathlib import Path
import json

Q = 2**30
G = [
    [[3376576,-4595,'8890485/1048576','-141059/16384'],
     [-4595,82702,'374101/128',-2966],
     ['8890485/1048576','374101/128',48657,-1435456],
     ['-141059/16384',-2966,-1435456,159900862]],
    [[2570,373,'100461333/8388608','-796959/65536'],
     [373,5619,'132105/32',-4192],
     ['100461333/8388608','132105/32',4654,-50824],
     ['-796959/65536',-4192,-50824,1442242]],
    [[20,26,'156439171/8388608','-18103/1024'],
     [26,822,'92585/128',-4442],
     ['156439171/8388608','92585/128',697,-4948],
     ['-18103/1024',-4442,-4948,53164]],
]
H = [
    [[142545776,374885,142247760,-564236],
     [374885,1002,374101,-1483],
     [142247760,374101,141950400,-563057],
     [-564236,-1483,-563057,2251]],
    [[201343592,529527,200922666,-796959],
     [529527,1401,528420,-2096],
     [200922666,528420,200502635,-795292],
     [-796959,-2096,-795292,3163]],
    [[313472891,94149,312878342,-1158592],
     [94149,16684,92585,-2221],
     [312878342,92585,312285067,-1156239],
     [-1158592,-2221,-1156239,4509]],
]
DIAG = [[F(128),F(1,32),F(1,2),F(1,32)],
        [F(1,512),F(1,16),F(1,512),F(1,32)]]
SOURCE = list(map(F, ['1','0','0','1166523000/31661820415017653',
    '322803000/186199784281','218288946200/1359122513',
    '90943434293/222507771414','248802430236/3300726103',
    '205476973/51646655494']))
Y_INTEGER = [-23803618,906393494300170,2744379904,1378590360912,
             5935589516,20425,2600028,91248,278104928]


def encode(A):
    return [[str(x) for x in row] for row in A]


def block(U,R,W):
    return [U[i]+R[i] for i in range(4)]+[R[i]+W[i] for i in range(4)]


def ldl(A):
    n=len(A)
    L=[[F(i==j) for j in range(n)] for i in range(n)]
    p=[]
    for j in range(n):
        pivot=A[j][j]-sum(L[j][k]**2*p[k] for k in range(j))
        assert pivot>F(1,2**27)
        p.append(pivot)
        for i in range(j+1,n):
            L[i][j]=(A[i][j]-sum(L[i][k]*L[j][k]*p[k] for k in range(j)))/pivot
    assert A==[[sum(L[i][k]*p[k]*L[j][k] for k in range(n))
                for j in range(n)] for i in range(n)]
    return {'L':encode(L),'diagonal':[str(x) for x in p],
            'pivot_floor_at_2pow27':[int(x*2**27) for x in p]}


hatted=[[[[F(x)/Q for x in row] for row in A] for A in side] for side in (G,H)]
original=[[[[16**k*d[i]*A[i][j]*d[j] for j in range(4)]
            for i in range(4)] for k,A in enumerate(side)]
          for d,side in zip(DIAG,hatted)]
for k in range(3):
    assert all(original[0][k][i][j]==original[1][k][i][j]
               for i in range(2) for j in range(2,4))
proofs=[{'M':ldl(block(*side)), 'R':ldl(side[1])} for side in hatted]
y=[F(x,2**42) for x in Y_INTEGER]
pairing=sum(a*b for a,b in zip(y,SOURCE))
assert pairing<0
out={
 'description':'Exact strict dual separator for the quadratic half-line completion cone; source-sign dual entries are zero.',
 'source_order':['d','lambda','n','m','C','A','D','alpha','beta'],
 'source':[str(x) for x in SOURCE],
 'source_contact':{'u':'1/13','v':'16','state':['140','-1/32','1/50','1']},
 'bob_scale':'16', 'state_congruences':[[str(x) for x in row] for row in DIAG],
 'hatted_moments':[[encode(A) for A in side] for side in hatted],
 'original_moments':[[encode(A) for A in side] for side in original],
 'hatted_LDL':proofs,
 'dual_source_coefficients':[str(x) for x in y],
 'integer_functional_scale':str(2**42), 'integer_functional':Y_INTEGER,
 'source_pairing':str(pairing), 'integer_source_pairing':str(pairing*2**42)
}
path=Path(__file__).with_name('dyadic_order33_separator.json')
path.write_text(json.dumps(out,indent=2)+'\n')
print('Built:',path.name)
print('All four LDL decompositions have every pivot > 2^-27.')
print('Exact negative pairing:',pairing)
