"""Materialize the 32 exact projection constraints, with rational coefficients.

Standard library only. No optimizer or numerical search.
The common planar membership lift is proved in proof/main.tex.
"""
from itertools import product
from fractions import Fraction as F
from pathlib import Path
import json

root=Path(__file__).resolve().parent
qnames=['a0','a1','b0','b1','E00','E01','E10','E11']
# Each entry is (matrix row, matrix column, local z index, coefficient).
U=[(0,0,0,1),(0,1,1,F(1,2)),(0,2,4,F(1,2)),(0,3,6,1),
   (0,4,7,F(1,2)),(1,1,2,1),(1,2,5,F(1,2)),(1,3,7,F(1,2)),
   (1,4,8,1),(2,2,9,1),(2,3,10,F(1,2)),(2,4,11,F(1,2)),
   (3,3,12,1),(3,4,13,F(1,2)),(4,4,14,1)]
V=[(0,0,3,1),(0,1,4,F(1,2)),(0,2,5,F(1,2)),(0,3,9,1),
   (1,1,6,1),(1,2,7,F(1,2)),(1,3,10,F(1,2)),
   (2,2,8,1),(2,3,11,F(1,2)),(3,3,15,1)]
ell={0:1,2:1,3:1,6:2,8:2,9:2,12:1,14:1,15:1}
M={1:2,7:4,13:2}
S={0:2,2:2,3:-2,8:8,9:-8,11:8,12:-2,14:-2,15:2}
projections=[]
for signs in product((-1,1),repeat=4):
    if signs[0]*signs[1]*signs[2]*signs[3]!=-1: continue
    for x,y in product((0,1),repeat=2):
        m=[0]*8; m[x]=1; m[2+y]=signs[2*x+y]
        ss=[0]*4+list(signs)
        projections.append({'CHSH_signs':list(signs),'selected_pair':[x,y],
                            'M_row':m,'S_row':ss})
assert len(projections)==32

# Independently enumerate all128 physical relabelings of the original pair.
physical=set()
base=((1,1),(1,-1))
for swap,ia,ib,*flips in product((0,1),repeat=7):
    ea,eb=flips[:2],flips[2:]
    m=[0]*8; ss=[0]*8
    m[2*swap+ia]=(-1)**ea[0]
    m[2*(1-swap)+ib]=(-1)**eb[0]
    for i,j in product((0,1),repeat=2):
        ai,bj=i^ia,j^ib
        x,y=(bj,ai) if swap else (ai,bj)
        ss[4+2*x+y]+=base[i][j]*(-1)**(ea[i]+eb[j])
    # W is invariant under M -> -M, by flipping every output.
    first=next(e for e in m if e)
    physical.add((tuple(first*e for e in m),tuple(ss)))
assert physical=={(tuple(p['M_row']),tuple(p['S_row'])) for p in projections}

blocks=[]; entries=[]; eqs=[]
for k,p in enumerate(projections):
    first_z=9+16*k   # variable0=constant1, variables1..8=observed coordinates
    for block,order,template in ((2*k,5,U),(2*k+1,4,V)):
        blocks.append({'block':block,'order':order})
        for i,j,z,coefficient in template:
            value=F(coefficient)
            entries.append([block,i,j,first_z+z,value.numerator,value.denominator])
    norm={0:-1,**{first_z+z:v for z,v in ell.items()}}
    em={first_z+z:v for z,v in M.items()}
    es={first_z+z:v for z,v in S.items()}
    for i,value in enumerate(p['M_row'],1):
        if value: em[i]=-value
    for i,value in enumerate(p['S_row'],1):
        if value: es[i]=-value
    eqs.extend([sorted(norm.items()),sorted(em.items()),sorted(es.items())])

assert len(blocks)==64 and sum(b['order'] for b in blocks)==288
assert len(eqs)==96 and len(entries)==800
assert max(e[3] for e in entries)==520
assert all(0<=i<=j<blocks[b]['order'] and den>0
           for b,i,j,var,num,den in entries)
data={
 'claim':'Globally sound outer constraints; exact support on every listed plane.',
 'not_claimed':'Exact membership of the complete eight-coordinate quantum set.',
 'observed_order':qnames,
 'variable_order':'0 is constant1;1..8 observed;9..520 are32 consecutive blocks of16 auxiliaries.',
 'matrix_format':'[block,row,col,variable,numerator,denominator], upper triangle; mirror off-diagonals once.',
 'equality_format':'Each list contains [variable,coefficient]; its affine sum equals zero.',
 'projection_count':32,'auxiliary_count':512,'equality_count':96,
 'PSD_block_count':64,'total_PSD_order':288,
 'projections':projections,'blocks':blocks,'matrix_entries':entries,'equalities':eqs}
path=root/'ALL32_EXACT_PROJECTION_CONSTRAINTS.json'
path.write_text(json.dumps(data,indent=2)+'\n')
loaded=json.loads(path.read_text())
assert loaded['projections']==projections and len(loaded['matrix_entries'])==800
print('Verified all128 physical actions reduce to32 distinct projection constraints.')
print('Materialized64 PSD blocks, total order288,512 auxiliaries and96 affine equations.')
print(path)
