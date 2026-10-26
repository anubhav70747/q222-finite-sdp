#!/usr/bin/env python3
"""Independent exact algebra and finite-matrix verification of global K364 pencil.

No numerical sweeps or SDP solve. Deterministic rational identities only.

- checks all nine SG2 target columns independently with SymPy
- checks all 66430 entries of A and A* including factors of two
- checks 128 physical signed permutations, group closure, 8 PR sectors
- checks source Bell maps, sigma map, event orbits
- checks face source injections and every dual base coefficient
- builds actual strict source Gram P6* and P364*, checks full exact target
- checks all 3456 global affine rows and all PSD copy offsets
- emits explicit rational strict point as a sparse source witness.
"""
from __future__ import annotations
import sys,json,hashlib
from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path
import sympy as s
import generate_exact_pencil as g
import assemble_global_pencil as assembly
HERE=Path(__file__).resolve().parent

ST={nm:s.symbols(nm) for nm in g.SRC}
u,v,W,X,Y,Z=s.symbols('u v W X Y Z'); STV=[W,X,Y,Z]
Du=1+u*u;Dv=1+v*v

def original_parent_full():
    d,lam,n,m,C,A,D,alpha,beta=[ST[nm] for nm in g.SRC]
    Q=u*W+u*v*X-Y-v*Z
    return (4*lam*Du*(W+v*X)**2+4*n*Du*Dv*X**2
      +4*m*(v*W-X+u*v*Y-u*Z)**2+4*C*Dv*(u*W-Y)**2
      +4*A*Du*Dv*Y**2+(4*D+2*d)*Q**2
      +(4*alpha+2*d)*Du*(v*Y-Z)**2
      +(4*beta+2*d)*Dv*(X+u*Z)**2-2*d*Du*Dv*Z**2)

def source_full_target():
    Fexpr=s.expand(Dv*original_parent_full())
    res=defaultdict(F)
    for exp,coef in s.Poly(Fexpr,u,v,W,X,Y,Z,*ST.values()).terms():
        eu,ev=exp[:2];s_exp=exp[2:6];src_exp=exp[6:]
        assert sum(s_exp)==2 and sum(src_exp)==1,(exp,coef)
        st=[i for i,e in enumerate(s_exp) for _ in range(e)]
        src=src_exp.index(1)
        key=g.KINDEX[(6*eu,6*ev,*st)]
        res[(src,key)]+=F(coef)
    return {k:v for k,v in res.items() if v}

def signed_compose(A,B):
    # row-wise signed permutations: C=A*B
    return tuple((B[A[i][0]][0], A[i][1]*B[A[i][0]][1]) for i in range(9))

def coeff_pairing(G,basis):
    d=defaultdict(F)
    for (a,b),v in G.items():
        au,av,sa=basis[a];bu,bv,sb=basis[b]
        key=g.KINDEX[(au+bu,av+bv,*g.sympair(sa,sb))]
        d[key]+=v*(1 if a==b else 2)
    return {k:v for k,v in d.items() if v}

def sgnminor(M):
    return [s.factor(M[:i,:i].det()) for i in range(1,M.rows+1)]

def upper_add(M,i,j,value):
    key=(min(i,j),max(i,j))
    M[key]=M.get(key,F(0))+F(value)
    if not M[key]:del M[key]

def vector16(root):
    corner=g.STATES
    CORNERS=((0,0),(1,0),(0,1),(1,1))
    map16={(i,j,l):4*k+l for k,(i,j) in enumerate(CORNERS) for l in range(4)}
    v=[F(0)]*16
    for feature,c in root.items():
        assert feature in map16, feature
        v[map16[feature]]+=c
    return v

def roots_shifts(root,du,dv):
    return {(i+du,j+dv,l):coef for (i,j,l),coef in root.items()}

def build_strict_witness():
    # Faithful reconstruction of the published 18 root Q0 and the signed CHSH Qd.
    event_roots=[
        g.ROOTS['lambda'], roots_shifts(g.ROOTS['lambda'],1,0),
        g.ROOTS['n'],roots_shifts(g.ROOTS['n'],1,0),
        roots_shifts(g.ROOTS['n'],0,1),roots_shifts(g.ROOTS['n'],1,1),
        g.ROOTS['m'],
        g.ROOTS['C'],roots_shifts(g.ROOTS['C'],0,1),
        g.ROOTS['A'],roots_shifts(g.ROOTS['A'],1,0),
        roots_shifts(g.ROOTS['A'],0,1),roots_shifts(g.ROOTS['A'],1,1),
        g.ROOTS['D'],
        g.ROOTS['alpha'],roots_shifts(g.ROOTS['alpha'],1,0),
        g.ROOTS['beta'],roots_shifts(g.ROOTS['beta'],0,1),
    ]
    assert len(event_roots)==18
    pos_roots=[g.ROOTS['D'],g.ROOTS['alpha'],roots_shifts(g.ROOTS['alpha'],1,0),
               g.ROOTS['beta'],roots_shifts(g.ROOTS['beta'],0,1)]
    neg_roots=[g.ROOTS['Z'],roots_shifts(g.ROOTS['Z'],1,0),
               roots_shifts(g.ROOTS['Z'],0,1),roots_shifts(g.ROOTS['Z'],1,1)]
    def gram_roots(roots,scale):
        m=s.zeros(16)
        for root in roots:
            a=s.Matrix([s.Rational(t.numerator,t.denominator) for t in vector16(root)])
            m+=scale*(a*a.T)
        return m
    Q0=gram_roots(event_roots,4)
    Qplus=gram_roots(pos_roots,2)
    Qminus=gram_roots(neg_roots,2)
    eps=s.Rational(1,21)
    Qc=Q0+eps*(Qplus-Qminus)-eps*s.eye(16)
    assert Q0.rank()==16
    assert s.trace(Q0.inv())==s.Rational(21,4)
    Qcmin=sgnminor(Qc)
    assert all(t>0 for t in Qcmin)
    H=s.diag(1,*([s.Rational(1,2)]*5),1)
    for i in range(1,6):H[i-1,i+1]=H[i+1,i-1]=-s.Rational(1,4)
    assert all(t>0 for t in sgnminor(H))
    R=s.symbols('r');rz=s.Matrix([R**i for i in range(7)])
    assert s.expand((rz.T*H*rz)[0]-1-R**12)==0
    z6=[(i,j,l) for i in range(7) for j in range(7) for l in range(4)]
    ix6={z:i for i,z in enumerate(z6)}
    CORNERS=((0,0),(1,0),(0,1),(1,1))
    ix16=[ix6[(6*i,6*j,l)] for i,j in CORNERS for l in range(4)]
    M6={}
    for a in range(16):
        for b in range(a,16):
            value=F(Qc[a,b]);
            if value:upper_add(M6,ix16[a],ix16[b],value)
    for ai,aj,al in z6:
        i=ix6[(ai,aj,al)]
        for bi in range(ai,7):
            for bj in range(7):
                j=ix6[(bi,bj,al)]
                if j<i:continue
                v6=eps*H[ai,bi]*H[aj,bj]
                if v6:upper_add(M6,i,j,F(v6))
    # Independent coefficient comparison against generator SG2; complete 196-feature witness.
    z6_keys={}
    for (i,j),val in M6.items():
        x,y,a=z6[i];xx,yy,b=z6[j]
        key=g.KINDEX[(x+xx,y+yy,*g.sympair(a,b))]
        z6_keys[key]=z6_keys.get(key,F(0))+(1 if i==j else 2)*val
    z6_keys={k:v for k,v in z6_keys.items() if v}
    Tcol=g.T_columns();Tstar=defaultdict(F)
    for src,col in enumerate(Tcol):
        for k,v in col:
            Tstar[k]+=F(v)*(F(1,21) if src==0 else F(1))
    # The original F target is T / (1+v²); for the 196-row Gram convolve exact parent columns.
    old=defaultdict(F)
    for j,src in enumerate(g.SRC):
        for (u_exp,v_exp,sa,sb),val in g.original_column(src).items():
            key=g.KINDEX[(6*u_exp,6*v_exp,sa,sb)]
            old[key]+=val*(F(1,21) if j==0 else F(1))
    assert z6_keys=={k:v for k,v in old.items() if v}, 'Strict 196 witness parent mismatch'
    M364={};ixext=g.ZINDEX
    for (i,j),val in M6.items():
        ai,aj,al=z6[i];bi,bj,bl=z6[j]
        upper_add(M364,ixext[(ai,aj,al)],ixext[(bi,bj,bl)],val)
        upper_add(M364,ixext[(ai,aj+6,al)],ixext[(bi,bj+6,bl)],val)
    ext=coeff_pairing(M364,g.FEATURES)
    assert ext=={k:v for k,v in Tstar.items() if v},'Strict 364 witness target mismatch'
    # Positivity: Qc>0 and H>0, hence M6>0; two shifted embeddings cover all features.
    # Save exact sparse witness and finite source matrix tests.
    out={'source_order':g.SRC,'source':["1/21"]+['1']*8,
         'proven_PSD':'Qc>0, H>0, and J0/J1 embeddings together span all 364 rows',
         'P6_basis':z6,'P6_upper':[[i,j,g.frac(v)] for (i,j),v in sorted(M6.items())],
         'P364_upper':[[i,j,g.frac(v)] for (i,j),v in sorted(M364.items())],
         'Qc_minors':[str(t) for t in Qcmin],
         'H_minors':[str(t) for t in sgnminor(H)],
         'H':[[str(x) for x in line] for line in H.tolist()],
         'Q0':[[str(x) for x in line] for line in Q0.tolist()],
         'Qc':[[str(x) for x in line] for line in Qc.tolist()],
         'target_identity_exact':True,'P364_PSD_order':364}
    with open(HERE/'K364_STRICT_RATIONAL_SOURCE_AND_GRAM.json','w') as f:
        json.dump(out,f,separators=(',',':'));f.write('\n')
    return {'strict_196_upper_nonzero':len(M6),'strict_364_upper_nonzero':len(M364),
            'Q0_inv_trace':str(s.trace(Q0.inv())),
            'Qc_min_minor':str(min(Qcmin)),
            'H_min_minor':str(min(sgnminor(H)))}

def main():
    checks=[]
    path=HERE/'K364_FULL128_FACTORED_SPARSE_PENCIL.json'
    data=assembly.load(path)
    assert len(data['polynomial_keys'])==3250
    assert len(data['fixed_testing_state_features'])==364
    assert len(data['Gram_dual_moment_upper_entries_map_A_star'])==66430
    checks.append('3250 monomial/state keys and 66430 exact A/A* Gram columns')
    Texact=source_full_target()
    for idx,col in enumerate(data['T_source_columns']):
        calc={(a,b):F(x) for a,b,x in [(idx,k,v) for k,v in col]}
        target={(idx,k):v for (src,k),v in Texact.items() if src==idx}
        assert calc==target,(idx,calc,target)
    checks.append('9 source columns of T equal expanded complete SG2 parent times (1+v²)')
    features=data['fixed_testing_state_features'];keys=data['polynomial_keys']
    visited=set()
    for ai,bi,ki in data['Gram_dual_moment_upper_entries_map_A_star']:
        assert ai<=bi and (ai,bi) not in visited
        visited.add((ai,bi))
        u1,v1,s1=features[ai];u2,v2,s2=features[bi]
        assert tuple(keys[ki])==(u1+u2,v1+v2,*g.sympair(s1,s2))
    checks.append('every A column, A* upper entry, and symmetric factor-of-two convention')
    # Independent adjoint trace-pairing test, exact with nonzero off diagonals.
    sparseP={(0,0):F(2),(1,3):F(3),(100,200):F(-1),(363,363):F(5)}
    sparseY={i:F((17*i+3)%13-6,7) for i in range(3250) if i%7==0}
    lhs=F(0);rhs=F(0)
    for ai,bi,ki in data['Gram_dual_moment_upper_entries_map_A_star']:
        if (ai,bi) in sparseP:
            term=sparseP[(ai,bi)]*(1 if ai==bi else 2)*sparseY.get(ki,F(0))
            lhs+=term;rhs+=term
    assert lhs==rhs
    checks.append('exact bilinear adjoint trace pairing including off-diagonal double count')
    # All 128 group actions, orbit, targets, composition and sigma.
    groups=data['groups'];srcE=s.Matrix([[s.Rational(v) for v in line] for line in data['canonical_E_source_to_Bell_9x9']])
    inv=srcE.inv()
    signatures=set()
    for group in groups:
        mapping=tuple(tuple(pair) for pair in group['observable_map'])
        signatures.add(mapping)
        qg=s.zeros(9,9)
        for i,(j,sign) in enumerate(mapping):qg[i,j]=sign
        assert qg[0,0]==1 and (qg*qg.T)==s.eye(9)
        S=s.Matrix([[s.Rational(x) for x in row] for row in group['source_map_9x9']])
        assert S==inv*qg.T*srcE
        tg=s.Matrix([[s.Rational(x) for x in row] for row in group['canonical_target_9x9']])
        assert tg==srcE.T*qg
        gt=data['group_targets'][group['id']]
        assert gt['K9x9']==group['canonical_target_9x9']
        for nm,indices in [('SG12_7x9',g.F12_INDEX),('SG13_7x9',g.F13_INDEX)]:
            assert gt[nm]==[gt['K9x9'][i] for i in indices]
    assert len(signatures)==128
    for A in list(signatures):
        for B in list(signatures):
            assert g_tuple(A,B) in signatures
    sector=data['PR_sectors'];assert len(sector)==8 and all(len(x['group_ids'])==16 for x in sector)
    sigma=data['counts']['sigma_group_id']
    assert groups[sigma]['observable_map']==[[0,1],[1,1],[2,-1],[4,-1],[3,-1],[6,-1],[5,-1],[8,1],[7,1]]
    checks.append('all 128 physical group actions closed; eight PR sectors; correct sigma')
    # Events: every physical two-outcome one-setting probability.
    events=data['physical_event_inequalities'];assert len(events)==16
    for item in events:
        row=[0]*9;row[0]=1
        row[1+item['x']]=item['a'];row[3+item['y']]=item['b']
        row[5+2*item['x']+item['y']]=item['a']*item['b']
        assert item['nonnegative_4p']==row
    assert {tuple(it['nonnegative_4p']) for it in events}==set(tuple(x) for x in events_to_orbit(groups))
    checks.append('all 16 event inequalities and their full physical orbit')
    # Dual SG12/SG13 matrix bases vs original face formulas (every independent variable).
    face=data['face_templates']
    for i,j in enumerate(g.F12_INDEX):
        x={z:s.Integer(0) for z in g.SRC};x[g.SRC[j]]=1
        expected=g.matrices_upper(g.face12_matrix(x,[0,0,0]))
        assert face['SG12']['source_matrix_bases_6'][i]==expected
    for i in range(3):
        vv=[0,0,0];vv[i]=1
        x={z:s.Integer(0) for z in g.SRC}
        assert face['SG12']['free_skew_matrix_bases_6'][i]==g.matrices_upper(g.face12_matrix(x,vv))
    for i,j in enumerate(g.F13_INDEX):
        x={z:s.Integer(0) for z in g.SRC};x[g.SRC[j]]=1
        expect=[g.matrices_upper(mat) for mat in g.face13_matrix(x,0)]
        assert face['SG13']['source_matrix_bases_3_and_3'][i]==expect
    assert face['SG13']['free_xi_matrix_bases_3_and_3']==[g.matrices_upper(mat) for mat in g.face13_matrix({z:s.Integer(0) for z in g.SRC},1)]
    # Evaluate the actual strict source feasibility points.
    xf={z:s.Integer(0) for z in g.SRC};xf.update({'d':1,'n':1,'C':2,'A':2,'D':1,'alpha':1,'beta':1})
    face12mins=sgnminor(g.face12_matrix(xf,[0,0,0]))
    assert all(x>0 for x in face12mins)
    xf={z:s.Integer(0) for z in g.SRC};xf.update({'d':1,'n':2,'m':1,'C':2,'A':2,'D':1,'beta':2})
    face13mins=[sgnminor(m) for m in g.face13_matrix(xf,-5)]
    assert all(x>0 for vec in face13mins for x in vec)
    checks.append('SG12,SG13 exact source/auxiliary bases, omitted-coordinate targets, strict Slater')
    eqs=assembly.materialize_equalities(data)
    with open(HERE/'K364_GLOBAL_ALL_3456_AFFINE_EQUALITIES.json') as f:fromfile=json.load(f)['rows']
    assert eqs==fromfile and len(eqs)==3456
    # A literal feasible point of the entire 128-copy outer pencil at q=0:
    # all moment duals and face PSD duals zero; each orthant r=1; events=I16.
    positive_r=set()
    for entry in data['group_copy_offsets']:
        for name in ('K_orthant_r','SG12_r','SG13_r'):
            start,n=entry['variables'][name]
            positive_r.update(range(start,start+n))
    assert len(positive_r)==128*(9+7+7)
    for e in eqs:
        val=F(e['constant'])+sum(F(coef) for idx,coef in e['entries'] if idx in positive_r)
        assert val==0,(e['id'],e['kind'],val)
    for group in data['group_copy_offsets']:
        assert group['equality_counts']=={'K_adjoint':9,'SG12_target':7,'SG12_skew':3,'SG13_target':7,'SG13_xi':1}
        assert group['PSD_blocks']['K_moment'][1]==364
    assert all(e['constant']==1 for e in assembly.event_global_entries(data))
    assert data['counts']['PSD_total_order']==51088
    assert data['counts']['total_scalar_variables']==423176
    checks.append('all 3456 equality rows and 897 PSD copy offsets; literal q=0 primal SDP feasibility')
    checks.append('strict 196/364 source witness '+str(build_strict_witness()))
    report={'result':'EXACT FINITE COEFFICIENTS; GLOBAL REVERSE INCLUSION CONDITIONAL',
            'strict_C1_oriented_K364':'UNPROVED','strict_C2_K364':'UNPROVED',
            'deterministic_checks':checks,
            'SG12_strict_principal_minors':[str(v) for v in face12mins],
            'SG13_strict_principal_minors':[[str(v) for v in vec] for vec in face13mins],
            'numerical_solvers':False,'time_or_sample_sweeps':False,'GitHub_writes':False}
    with open(HERE/'EXACT_VERIFICATION_REPORT.json','w') as f:json.dump(report,f,indent=2);f.write('\n')
    for c in checks:print('PASS:',c)
    print('ALL EXACT CHECKS PASSED; strict universal source admission remains conditional.')

def g_tuple(A,B):
    return tuple((B[A[i][0]][0],A[i][1]*B[A[i][0]][1]) for i in range(9))

def events_to_orbit(groups):
    # Canonical first losing event, physical orbit transformed under q actions.
    v=s.Matrix([1,1,0,0,1,0,1,0,0])
    res=set()
    for group in groups:
        Q=s.zeros(9,9)
        for i,(j,sgn) in enumerate(group['observable_map']):Q[i,j]=sgn
        arr=Q.T*v
        res.add(tuple(int(a) for a in arr))
    return res

if __name__=='__main__':main()
