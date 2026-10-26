#!/usr/bin/env python3
"""Generate a fully indexed rational-data 128-copy Q222 outer SDP pencil.

No optimization, network, searches, or floating point. Requires SymPy solely
for a fixed 9x9 rational basis inversion and small exact sanity checks.
Source and target conventions: SG1/SG2, SG12/SG13.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import gzip, hashlib, json
import sympy as sp

HERE = Path(__file__).resolve().parent
SRC = ('d','lambda','n','m','C','A','D','alpha','beta')
OBS = ('1','A0','A1','B0','B1','E00','E01','E10','E11')
STATES = ('W','X','Y','Z')
PAIRS = [(a,b) for a in range(4) for b in range(a,4)]
PAIR_INDEX={p:i for i,p in enumerate(PAIRS)}
KEYS=[(I,J,a,b) for I in range(13) for J in range(25) for a,b in PAIRS]
KINDEX={k:i for i,k in enumerate(KEYS)}
FEATURES=[(i,j,s) for i in range(7) for j in range(13) for s in range(4)]
ZINDEX={z:i for i,z in enumerate(FEATURES)}
assert len(KEYS)==3250 and len(FEATURES)==364

def frac(x):
    x=F(x)
    return str(x.numerator) if x.denominator==1 else f'{x.numerator}/{x.denominator}'
def tupkey(t): return (t[0],t[1],min(t[2],t[3]),max(t[2],t[3]))
def form(terms):
    d=defaultdict(F)
    for eu,ev,s,co in terms: d[(eu,ev,s)]+=F(co)
    return {k:v for k,v in d.items() if v}
def polynomial_prod(p,q):
    d=defaultdict(F)
    for (x,y),v in p.items():
        for (a,b),w in q.items(): d[(x+a,y+b)]+=v*w
    return {k:v for k,v in d.items() if v}
def squared(root,prefactor):
    """Returns coefficient form indexed by integer u/v degree and state pair."""
    d=defaultdict(F); items=sorted(root.items())
    for i,((eu,ev,s),a) in enumerate(items):
        for (fu,fv,t),b in [t for t in items[i:]]:
            mult=F(a)*b*(1 if (eu,ev,s)==(fu,fv,t) else 2)
            for (mu,mv),w in prefactor.items():
                key=(eu+fu+mu,ev+fv+mv,min(s,t),max(s,t))
                d[key]+=mult*w
    return {k:v for k,v in d.items() if v}
ONE={(0,0):F(1)}
DU={(0,0):F(1),(2,0):F(1)}
DV={(0,0):F(1),(0,2):F(1)}
DUDV=polynomial_prod(DU,DV)

# The original SG2 parent, reconstructed term by term, not from a determinant.
ROOTS={
    'lambda':form([(0,0,0,1),(0,1,1,1)]),
    'n':form([(0,0,1,1)]),
    'm':form([(0,1,0,1),(0,0,1,-1),(1,1,2,1),(1,0,3,-1)]),
    'C':form([(1,0,0,1),(0,0,2,-1)]),
    'A':form([(0,0,2,1)]),
    'D':form([(1,0,0,1),(1,1,1,1),(0,0,2,-1),(0,1,3,-1)]),
    'alpha':form([(0,1,2,1),(0,0,3,-1)]),
    'beta':form([(0,0,1,1),(1,0,3,1)]),
    'Z':form([(0,0,3,1)])
}

ORIG={
    'lambda':[(F(4),ROOTS['lambda'],DU)],
    'n':[(F(4),ROOTS['n'],DUDV)],
    'm':[(F(4),ROOTS['m'],ONE)],
    'C':[(F(4),ROOTS['C'],DV)],
    'A':[(F(4),ROOTS['A'],DUDV)],
    'D':[(F(4),ROOTS['D'],ONE)],
    'alpha':[(F(4),ROOTS['alpha'],DU)],
    'beta':[(F(4),ROOTS['beta'],DV)],
    'd':[(F(2),ROOTS['D'],ONE),
         (F(2),ROOTS['alpha'],DU),
         (F(2),ROOTS['beta'],DV),
         (F(-2),ROOTS['Z'],DUDV)]
}

def original_column(src):
    d=defaultdict(F)
    for coef,root,mult in ORIG[src]:
        for key,v in squared(root,mult).items():d[key]+=coef*v
    return {k:v for k,v in d.items() if v}

def T_columns():
    cols=[]
    for src in SRC:
        c=original_column(src);dv=defaultdict(F)
        for (i,j,s,t),v in c.items():
            for (vpow,w) in [(0,F(1)),(2,F(1))]:dv[(6*i,6*(j+vpow),s,t)]+=v*w
        cols.append(sorted([[KINDEX[key],frac(x)] for key,x in dv.items() if x]))
    return cols

def sympair(a,b): return (a,b) if a<=b else (b,a)

def A_upper():
    """Gram coefficient columns: A(P)_key = sum_{i<=j} mult P_ij.
    A*(y)[i,j]=y_key in all upper entries, so this is a materialized
    exact base moment pencil; no unknown linear operator remains.
    """
    out=[]; hist=[0]*len(KEYS)
    for i,(iu,iv,sa) in enumerate(FEATURES):
        for j in range(i,len(FEATURES)):
            ju,jv,sb=FEATURES[j]
            key=KINDEX[(iu+ju,iv+jv,*sympair(sa,sb))]
            out.append([i,j,key]);hist[key]+=1
    assert len(out)==364*365//2 and all(hist)
    return out,hist

def canonical_E():
    M=sp.zeros(9,9)
    M[0,0]=1
    for idx,x in [(5,-1),(6,1),(7,-1),(8,-1)]:M[idx,0]=sp.Rational(x,2)
    events=[(0,1,1,1),(0,0,1,-1),(1,1,1,-1),
            (1,0,-1,1),(0,0,-1,1),(1,1,-1,1),
            (0,1,-1,-1),(1,0,1,-1)]
    for col,(x,y,a,b) in enumerate(events,1):
        M[0,col]=1;M[1+x,col]=a;M[3+y,col]=b;M[5+2*x+y,col]=a*b
    assert M.det()!=0
    return M

def q_action(swap,sa,sb,af0,af1,bf0,bf1):
    af=(af0,af1);bf=(bf0,bf1)
    M=sp.zeros(9,9);M[0,0]=1
    for x in range(2):
        M[1+x, (3 if swap else 1)+(x^sa)]=af[x]
        M[3+x, (1 if swap else 3)+(x^sb)]=bf[x]
    for x in range(2):
        for y in range(2):
            xx=x^sa;yy=y^sb
            orig=5+(2*yy+xx if swap else 2*xx+yy)
            M[5+2*x+y,orig]=af[x]*bf[y]
    assert all(sum(abs(M[i,j]) for j in range(9))==1 for i in range(9))
    return M

def group_data(E):
    Einv=E.inv();entries=[];hashes=set();prsectors={}
    PR0=sp.Matrix([1,0,0,0,0,1,-1,1,1])
    SIGMA=sp.zeros(9,9)
    for i,j in enumerate([0,2,1,4,3,7,8,5,6]):SIGMA[i,j]=1
    sigma_index=None
    for g,(swap,sa,sb,af0,af1,bf0,bf1) in enumerate(product([0,1],[0,1],[0,1],[-1,1],[-1,1],[-1,1],[-1,1])):
        act=q_action(swap,sa,sb,af0,af1,bf0,bf1)
        raw=tuple(int(t) for t in act)
        assert raw not in hashes
        hashes.add(raw)
        src=Einv*act.T*E
        assert all(x.q in (1,2,4,8) for x in src)
        target=E.T*act
        assert src.T*E.T==target
        pr=tuple(int(x) for x in (act*PR0)[5:9])
        prsectors.setdefault(pr,[]).append(g)
        if src==SIGMA:sigma_index=g
        entries.append({
            'id':g,
            'parameters':[swap,sa,sb,af0,af1,bf0,bf1],
            'observable_map':[[int(next(j for j in range(9) if act[i,j])),int(act[i,next(j for j in range(9) if act[i,j])])] for i in range(9)],
            'pr_sector_correlations':list(pr),
            'source_map_9x9':[[frac(F(x)) for x in row] for row in src.tolist()],
            'canonical_target_9x9':[[frac(F(x)) for x in row] for row in target.tolist()],
        })
    assert len(entries)==128 and len(prsectors)==8 and all(len(ids)==16 for ids in prsectors.values())
    assert sigma_index is not None
    return entries,prsectors,sigma_index

F12_INDEX=[0,2,4,5,6,7,8]
F13_INDEX=[0,2,3,4,5,6,8]

def face12_matrix(x,sk):
    # x: map of all nine source coefficients. sk: free skew entries (01,02,12)
    d,n,C,A,D,alpha,beta=(x[k] for k in ('d','n','C','A','D','alpha','beta'))
    K0=sp.diag(4*C+4*D+2*d,4*A,4*alpha)
    K2=sp.diag(4*C,4*A+4*alpha+2*d,4*D)
    K1=sp.zeros(3,3)
    K1[0,2]=K1[2,0]=-(4*D+2*d)
    K1[1,2]=K1[2,1]=-(4*alpha+2*d)
    S=sp.zeros(3,3)
    for val,(i,j) in zip(sk,[(0,1),(0,2),(1,2)]):S[i,j]=val;S[j,i]=-val
    return (K0.row_join(K1/2+S)).col_join((K1/2-S).row_join(K2))

def face13_matrix(x,xi):
    d,n,m,C,A,D,beta=(x[k] for k in ('d','n','m','C','A','D','beta'))
    p=4*C+4*D+2*d;q=4*m+4*beta+2*d;r=4*D+2*d
    return (sp.Matrix([[4*m,-4*m,xi],[-4*m,q+4*n,q],[xi,q,q]]),
            sp.Matrix([[p+4*A,4*D,-p],[4*D,4*D,-(4*m+r)-xi],[-p,-(4*m+r)-xi,p]]))

def matrices_upper(M):
    return [[i,j,frac(F(M[i,j]))] for i in range(M.rows) for j in range(i,M.cols) if M[i,j]]

def face_templates():
    zer={k:sp.Integer(0) for k in SRC}
    f12=[];f13=[]
    for idx in F12_INDEX:
        x=dict(zer);x[SRC[idx]]=sp.Integer(1)
        f12.append(matrices_upper(face12_matrix(x,[0,0,0])))
    skew=[]
    for k in range(3):
        sk=[0,0,0];sk[k]=1
        skew.append(matrices_upper(face12_matrix(zer,sk)))
    for idx in F13_INDEX:
        x=dict(zer);x[SRC[idx]]=sp.Integer(1)
        f13.append([matrices_upper(mat) for mat in face13_matrix(x,0)])
    xi=[matrices_upper(mat) for mat in face13_matrix(zer,1)]
    return {'SG12':{'source_indices':F12_INDEX,'source_injection_9x7':[[int(i==idx) for idx in F12_INDEX] for i in range(9)],
                    'source_matrix_bases_6':f12,'free_skew_matrix_bases_6':skew,
                    'free_skew_order':[[0,1],[0,2],[1,2]]},
            'SG13':{'source_indices':F13_INDEX,'source_injection_9x7':[[int(i==idx) for idx in F13_INDEX] for i in range(9)],
                    'source_matrix_bases_3_and_3':f13,'free_xi_matrix_bases_3_and_3':xi}}

def event_matrix():
    rows=[]
    for x,y,a,b in product(range(2),range(2),(-1,1),(-1,1)):
        row=[0]*9;row[0]=1;row[1+x]=a;row[3+y]=b;row[5+2*x+y]=a*b
        rows.append({'x':x,'y':y,'a':a,'b':b,'nonnegative_4p':row})
    return rows

def block_layout():
    # An evaluated global LMI is the direct sum of these genuine sparse matrices.
    names=[('K_moment',364),('K_source_orthant',9),('SG12_dual_PSD',6),('SG12_orthant',7),
           ('SG13_dual_PSD_1',3),('SG13_dual_PSD_2',3),('SG13_orthant',7)]
    varnames=[('K_dual_coeffs_y',3250),('K_orthant_r',9),('SG12_dual_upper_Y',21),('SG12_r',7),
              ('SG13_dual_upper_Y1',6),('SG13_dual_upper_Y2',6),('SG13_r',7)]
    copy=[]
    for g in range(128):
        v=8+3306*g;b=399*g
        vars={};blks={}
        for name,n in varnames:vars[name]=[v,n];v+=n
        for name,n in names:blks[name]=[b,n];b+=n
        assert v==8+3306*(g+1) and b==399*(g+1)
        copy.append({'group_id':g,'variables':vars,'PSD_blocks':blks,'equality_offset':27*g,
                     'equality_counts':{'K_adjoint':9,'SG12_target':7,'SG12_skew':3,'SG13_target':7,'SG13_xi':1}})
    return copy

def build_data():
    E=canonical_E();groups,sectors,sigmaidx=group_data(E)
    T=T_columns();A,hist=A_upper();faces=face_templates();events=event_matrix();copy=block_layout()
    scts=[{'correlators':list(k),'group_ids':v} for k,v in sorted(sectors.items())]
    # Explicit canonical and face targets per group; these are the physical q adjoints.
    targets=[]
    for g in groups:
        K=g['canonical_target_9x9']
        targets.append({'group_id':g['id'],'K9x9':K,
                        'SG12_7x9':[K[i] for i in F12_INDEX],
                        'SG13_7x9':[K[i] for i in F13_INDEX]})
    counts={'features':364,'polynomial_state_coefficient_keys':3250,
            'Gram_upper_independent_coordinates':364*365//2,
            'Gram_map_nonzero_columns':len(A),'Gram_moment_keys_with_nonzero_entries':sum(n>0 for n in hist),
            'T_nonzero_terms':sum(len(col) for col in T),'physical_symmetry_actions':128,'PR_sectors':8,
            'source_cone_K_symmetric_PSD_blocks':128,'SG12_dual_PSD_blocks':128,'SG13_dual_PSD_blocks':256,
            'K_orthant_diagonal_blocks_9':128,'SG12_orthant_diagonal_blocks_7':128,
            'SG13_orthant_diagonal_blocks_7':128,'event_diagonal_block_16':1,
            'PSD_blocks_total':897,'PSD_total_order':51088,
            'auxiliary_variables':423168,'physical_q_variables':8,'total_scalar_variables':423176,
            'affine_equality_count':3456,'event_inequality_count':16,
            'equality_per_group':27,'auxiliary_variables_per_group':3306,
            'primal_certificate_Gram_unknowns_not_in_global_dual_pencil':364*365//2,
            'sigma_group_id':sigmaidx,
            'identity_group_id':next(g['id'] for g in groups if g['observable_map']==[[i,1] for i in range(9)])}
    assert sum(k['Gram_moment_keys_with_nonzero_entries'] for k in [counts])==3250
    meta={'kind':'FULLY_MATERIALIZED_FACTORED_RATIONAL_LMI','version':'1.0',
          'main_ref':'3d193cc25ae2fcdfd3817c0bf5ebd54dfada7b55',
          'source_order':SRC,'observable_homogeneous_order':OBS,'state_order':STATES,
          'polynomial_keys':KEYS,'fixed_testing_state_features':FEATURES,
          'Gram_variable_upper_entries_map_A':[{'entry':[a,b],'coefficient_key':k,'A_scalar_on_symmetric_variable':1 if a==b else 2} for a,b,k in A],
          'Gram_dual_moment_upper_entries_map_A_star':A,
          'T_source_columns':T,
          'canonical_E_source_to_Bell_9x9':[[frac(F(x)) for x in row] for row in E.tolist()],
          'groups':groups,'group_targets':targets,'PR_sectors':scts,
          'face_templates':faces,'physical_event_inequalities':events,
          'group_copy_offsets':copy,
          'physical_q_variable_positions':[f'q_{s}' for s in OBS[1:]],
          'source_dual_conventions':{
              'K_target':'E^T Q_g (1,q)=T^T y_g+r_g',
              'K_LMI':'A^* y_g >= 0 with A^* y_g[i,j]=y_g[key(i,j)]',
              'K_orthant':'diag(r_g)>=0',
              'SG12_target':'I12^T E^T Q_g (1,q) = (Tr(M12_k Y12_g))_k+r12_g',
              'SG12_aux':'Tr(S12_t Y12_g)=0 for t=0,1,2',
              'SG13_target':'I13^T E^T Q_g (1,q) = (Tr(B13_k,1 Y1_g)+Tr(B13_k,2 Y2_g))_k+r13_g',
              'SG13_aux':'Tr(B13_xi,1 Y1_g)+Tr(B13_xi,2 Y2_g)=0',
              'faces_oriented_source_orthant':'diag(r_face)>=0',
              'event':'diag(1 +a A_x+b B_y+ab E_xy)>=0, all 16'},
          'counts':counts}
    return meta

def writejson(path,data):
    with open(path,'w',encoding='utf-8',newline='\n') as f:json.dump(data,f,separators=(',',':'),ensure_ascii=True);f.write('\n')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default=str(HERE));args=ap.parse_args()
    p=Path(args.output);p.mkdir(exist_ok=True,parents=True)
    data=build_data();out=p/'K364_FULL128_FACTORED_SPARSE_PENCIL.json'
    writejson(out,data)
    with gzip.open(p/'K364_FULL128_FACTORED_SPARSE_PENCIL.json.gz','wt',encoding='utf-8',compresslevel=9) as f:json.dump(data,f,separators=(',',':'))
    print('WROTE',out,'bytes',out.stat().st_size)
    print(json.dumps(data['counts'],indent=2))

if __name__=='__main__':main()
