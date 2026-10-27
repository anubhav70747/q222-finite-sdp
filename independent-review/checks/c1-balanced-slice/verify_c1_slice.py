#!/usr/bin/env python3
"""Exact global C1 slice and a quantitative five-parameter admission tube.

SymPy only. Reconstructs the full parent, eight constant rational Grams,
all coefficient identities, the order-five source/behaviour pencils and
the contact right inverse. No optimizer, numerical sweep or network.
Universal H1/H2 outside the displayed slice/tube are not claimed.
"""
from pathlib import Path
import argparse,itertools,json
import sympy as s
import native_c1 as nc

BASE=Path(__file__).resolve().parent
r,x,y,W,X,Y,Z=s.symbols('r x y W X Y Z');PSI=(W,X,Y,Z)
u,v=s.symbols('u v')
FEATURES=[u*W,v*W,X,u*v*X,Y,u*v*Y,u*Z,v*Z,
          u*v*W,v*X,u*Y,v*Y,Z,u*v*Z]
SPLIT=[FEATURES[:8],FEATURES[8:],[s.sqrt(u*v)*W,s.sqrt(u*v)*Z],
       [s.sqrt(u*v)*X,s.sqrt(u*v)*Y]]

def parent(weights):
    d,la,n,m,C,A,D,al,be=weights;du,dv=1+u*u,1+v*v
    return (4*la*du*(W+v*X)**2+4*n*du*dv*X**2
            +4*m*(v*W-X+u*v*Y-u*Z)**2+4*C*dv*(u*W-Y)**2
            +4*A*du*dv*Y**2+(4*D+2*d)*(u*W+u*v*X-Y-v*Z)**2
            +(4*al+2*d)*du*(v*Y-Z)**2+(4*be+2*d)*dv*(X+u*Z)**2
            -2*d*du*dv*Z**2)

def root_blocks():
    uu,vv=x*x,y*y
    h=[vv*W-uu*W,uu*Z-uu*W,X+Y,uu*vv*X+uu*vv*Y,X-uu*vv*X]
    g=[uu*vv*W-Z,uu*Y-vv*Y,vv*X+uu*Y]
    hw=x*y*(W-Z);hx=x*y*(X+Y)
    roots=[h+[X+r*uu*W,(r-s.Rational(1,2))*(vv*Z-uu*W)]+[r*f for f in h],
           g+[uu*Y-r*Z,(r-s.Rational(1,3))*(uu*vv*Z-Z)]+[r*f for f in g],
           [hw,r*hw],[hx,r*hx],h+[vv*Z-uu*W],g+[uu*vv*Z-Z],[hw],[hx]]
    omega=(r-s.Rational(1,3))*(s.Rational(1,2)-r)
    return roots,[s.S.One]*4+[omega]*4

def ldl(A):
    M=A.copy();piv=[]
    for k in range(M.rows):
        p=s.cancel(M[k,k]);assert p>0,(k,p);piv.append(p)
        for i in range(k+1,M.rows):
            for j in range(i,M.rows):M[i,j]=M[j,i]=s.cancel(M[i,j]-M[i,k]*M[k,j]/p)
    return piv

def prototype(t):
    return [8,0,0,2,2,2*(1-t)/t,2*(1-2*t)/t,4*t/(1-t),2*(3*t-1)/(1-t)]

def coefficients(expr):
    out={}
    for exp,val in s.Poly(s.expand(expr),u,v,W,X,Y,Z).terms():
        st=[i for i,n in enumerate(exp[2:]) for _ in range(n)]
        out[(exp[0],exp[1],*st)]=val
    return out

def contact_map():
    blocks=[s.Matrix([f-b[-1] for f in b[:-1]]) for b in SPLIT[:2]]
    gamma=u*v*(W-Z)**2;delta=u*v*(X-Y)**2
    columns=[];positions=[]
    for bi,b in enumerate(blocks):
        for i in range(len(b)):
            for j in range(i,len(b)):
                columns.append(coefficients(b[i]*b[j]*(1 if i==j else 2)));positions.append((bi,i,j))
    columns += [coefficients(gamma),coefficients(delta)]
    keys=sorted(set().union(*(set(c) for c in columns)));ki={k:i for i,k in enumerate(keys)}
    M=s.zeros(len(keys),len(columns))
    for j,c in enumerate(columns):
        for k,val in c.items():M[ki[k],j]=val
    ct={u:1,v:1,W:1,X:1,Y:1,Z:1}
    J=s.Matrix([[s.diff(u**i*v**j*PSI[l]*PSI[k],z).subs(ct)
                 for i,j,l,k in keys] for z in (u,v,W,X,Y,Z)])
    pc=list(M.rref()[1]);pr=list(M[:,pc].T.rref()[1]);minor=M.extract(pr,pc)
    assert M.shape==(42,45) and len(pc)==36 and minor.det()==-2**24
    assert J.rank()==6 and J*M==s.zeros(6,45)
    E=s.zeros(45,42);inv=minor.inv()
    for i,p in enumerate(pc):
        for j,row in enumerate(pr):E[p,row]=inv[i,j]
    assert M*E*M==M and E[43:,:]==s.zeros(2,42)
    bounds=[max(sum(abs(z) for z in E[i,:]) for i in range(a,b)) for a,b in [(0,28),(28,43)]]
    assert bounds==[40,s.Rational(17,2)]
    return M,J,E,keys,positions,pc,pr

def root_coefficient_matrices(roots):
    test=[b if bi<2 else [x*y*W,x*y*Z] if bi==2 else [x*y*X,x*y*Y] for bi,b in enumerate(SPLIT)]
    test=[ [s.expand(f.subs({u:x*x,v:y*y})) for f in b] for b in test]
    out=[]
    for bi,root in enumerate(roots):
        b=test[bi%4];M=s.zeros(len(b),len(root))
        for j,f in enumerate(root):
            for exp,c in s.Poly(s.expand(f),x,y,W,X,Y,Z).terms():
                term=x**exp[0]*y**exp[1]*W**exp[2]*X**exp[3]*Y**exp[4]*Z**exp[5]
                i=b.index(term);M[i,j]=c
        out.append(M)
    return out

def source_pencil_data():
    # Source parameter order (m,A,alpha); canonical order d,lambda,n,m,C,A,D,alpha,beta.
    S=s.Matrix([[4,0,0],[0,0,0],[0,0,0],[1,0,0],[1,0,0],[0,1,0],[-1,1,0],[0,0,1],[-1,0,1]])
    E=s.zeros(9);E[0,0]=1
    for i,c in [(5,-1),(6,1),(7,-1),(8,-1)]:E[i,0]=s.Rational(c,2)
    events=[(0,1,1,1),(0,0,1,-1),(1,1,1,-1),(1,0,-1,1),
            (0,0,-1,1),(1,1,-1,1),(0,1,-1,-1),(1,0,1,-1)]
    for j,(x0,y0,a,b) in enumerate(events,1):
        E[0,j]=1;E[1+x0,j]=a;E[3+y0,j]=b;E[5+2*x0+y0,j]=a*b
    return S,E

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=BASE/'exact_results.json');args=ap.parse_args()
    data=json.loads((BASE/'certificate.json').read_text());roots,weights=root_blocks()
    assert [len(b) for b in roots]==[12,8,2,2,6,4,1,1]
    G=[s.zeros(len(b)) for b in roots]
    for val,(bi,i,j) in zip(data['gram_upper_entries'],data['positions']):G[bi][i,j]=G[bi][j,i]=s.Rational(val)
    pivots=[ldl(g-s.eye(g.rows)/4) for g in G]
    target=s.cancel(r*(1-r)*parent(prototype(r))).subs({u:x*x,v:y*y})
    certified=sum(w*(s.Matrix(b).T*g*s.Matrix(b))[0] for b,w,g in zip(roots,weights,G))
    assert s.Poly(s.expand(certified-target),r,x,y,W,X,Y,Z).is_zero
    print('PASS whole closed interval: eight constant rational Grams > I/4; complete four-state polynomial identity',flush=True)
    # Check the original native source dictionary, not an independently free-weight slice.
    sub={nc.k:1,nc.f:1,nc.g:1,nc.a:1-r,nc.t:1-r}
    assert all(s.cancel(a.subs(sub)-b)==0 for a,b in zip(nc.NATIVE,prototype(r)))
    assert all(s.cancel(a.subs(sub)-b)==0 for a,b in zip(nc.CONTACT,[1,1,1,-r,r,1]))
    z=s.symbols('z',positive=True);rpos=(3*z+2)/(6*(z+1))
    checks=[nc.Lc,nc.Hc,nc.Qc,nc.M,nc.J,1-nc.t*nc.g**2,1-nc.k**2*nc.f**2+nc.t*nc.vc,nc.P4,nc.Pa,nc.Pd,nc.c2,nc.cb,nc.c3,nc.cd,nc.cal,nc.D,nc.be]
    for expr in checks:
        num,den=s.fraction(s.cancel(expr.subs(sub).subs(r,rpos)))
        assert all(c>0 for c in s.Poly(num,z).coeffs()) and all(c>0 for c in s.Poly(den,z).coeffs())
    print('PASS genuine native family k=f=g=1, a=t=1-r, 1/3<r<1/2; all native signs and contact identities',flush=True)
    # Necessity of the global slice condition has a two-dimensional maximally
    # entangled test witness at orthogonal measurements, not only a determinant.
    m,A,al,z0,z1=s.symbols('m A al z0 z1',real=True)
    sw=[4*m,0,0,m,m,A,A-m,al,al-m]
    restricted=s.expand(parent(sw).subs({u:1,v:1,W:z0,X:-z1,Y:z1,Z:z0}))
    Q=s.Matrix([[al,-al-2*m],[-al-2*m,2*A+al+4*m]])
    assert s.expand(restricted-16*(s.Matrix([z0,z1]).T*Q*s.Matrix([z0,z1]))[0])==0
    assert s.expand(Q.det()-2*(A*al-2*m*m))==0
    S,Ebell=source_pencil_data();assert S*s.Matrix([m,A,al])==s.Matrix(sw)
    A0=2*m*m/A;assert all(s.cancel(a-b)==0 for a,b in zip([4*m,0,0,m,m,A,A-m,A0,A0-m],[m/s.Integer(2)*v for v in prototype(m/(A+m))]))
    mM=s.Matrix([[0,2],[2,0]]);AM=s.Matrix([[1,0],[0,0]]);alM=s.Matrix([[0,0],[0,2]])
    assert s.expand((m*mM+A*AM+al*alM).det()-2*(A*al-2*m*m))==0
    assert ldl(mM+2*AM+2*alM)==[s.Integer(2),s.Integer(2)]
    print('PASS exact global source slice: A>=m>=0, alpha>=m, A*alpha>=2m^2; rational order-five pencil and dual',flush=True)
    # Prototype Gram in the source-normalized C1 contact frame.
    coefmat=root_coefficient_matrices(roots)
    rawgrams=[coefmat[b]*G[b]*coefmat[b].T+weights[b+4]*coefmat[b+4]*G[b+4]*coefmat[b+4].T for b in range(4)]
    ct=[1,-r,r,1]
    diags=[]
    for bi,b in enumerate(SPLIT):
        vals=[]
        for f in b:
            vals.append(s.cancel(f.subs({u:1,v:1,W:ct[0],X:ct[1],Y:ct[2],Z:ct[3]})))
        diags.append(s.diag(*vals))
    normalized=[(D*g*D/(r*(1-r))).applyfunc(s.cancel) for D,g in zip(diags,rawgrams)]
    x0=s.Matrix([normalized[bi][i,j] for bi,n in enumerate((7,5)) for i in range(n) for j in range(i,n)]+[normalized[2][0,0],normalized[3][0,0]])
    M,J,E,keys,pos,pc,pr=contact_map()
    fn=parent(prototype(r)).subs({X:-r*X,Y:r*Y},simultaneous=True);co=coefficients(fn);b0=s.Matrix([co.get(k,0) for k in keys])
    assert (M*x0-b0).applyfunc(s.cancel)==s.zeros(42,1)
    # Explicit source-evaluation minors: only the distances to the two walls
    # occur in denominators. This checks the claimed inverse-entry list.
    inverse_values=[]
    for bi,n in enumerate((7,5)):
        B=coefmat[bi][:n,:n];inv=B.inv().applyfunc(s.factor);expected=-(2*r-1)/2 if bi==0 else (3*r-1)/3
        assert s.expand(B.det()-expected)==0
        wall=s.Rational(1,2) if bi==0 else s.Rational(1,3)
        values=set(map(str,inv));allowed={str(s.factor(f)) for f in [s.S.Zero,s.S.One,-s.S.One,1/(r-wall),-1/(r-wall),r/(r-wall),-r/(r-wall)]}
        assert values<=allowed,(bi,values);inverse_values.append(sorted(values))
    eta=s.Rational(1,24);radius=eta**2/s.Integer(246960);r0=s.Rational(2,5)
    # A strict native point moving off the one-parameter family in all five
    # parameters. The proved tube, not a numerical eigenvalue, admits it.
    eps=s.Rational(1,2**40);theta=[1-eps,1+eps,1-eps,s.Rational(3,5)+eps,s.Rational(3,5)-eps]
    sub1=dict(zip((nc.k,nc.f,nc.g,nc.a,nc.t),theta))
    assert all(s.cancel(e.subs(sub1))>0 for e in checks)
    w1=[s.cancel(e.subs(sub1)) for e in nc.NATIVE];ct1=[s.cancel(e.subs(sub1)) for e in nc.CONTACT]
    full=parent(w1).subs({u:ct1[0]*u,v:ct1[1]*v,W:ct1[2]*W,X:ct1[3]*X,Y:ct1[4]*Y,Z:ct1[5]*Z},simultaneous=True)
    co1=coefficients(full);assert not set(co1)-set(keys);b1=s.Matrix([co1.get(k,0) for k in keys]);assert J*b1==s.zeros(6,1)
    delta=b1-b0.subs(r,r0);distance=max(abs(v) for v in delta);assert distance<radius
    solved=x0.subs(r,r0)+E*delta;assert M*solved==b1
    h=[s.zeros(7),s.zeros(5)]
    for val,(bi,i,j) in zip(solved,pos):h[bi][i,j]=h[bi][j,i]=val
    hp=[ldl(g) for g in h];assert solved[-1]>0 and solved[-2]>0
    print('PASS fixed 42x45 contact map: rank36, image=six-jet kernel; correction bounds40 and17/2',flush=True)
    print('PASS uniform coefficient tube eta^2/246960; exact off-family five-parameter source and PSD Gram',flush=True)
    def serial(M):return [[str(v) for v in row] for row in M.tolist()]
    pencil={'source_order':['m','A','alpha'],'canonical_source_order':['d','lambda','n','m','C','A','D','alpha','beta'],
            'source_injection':serial(S),'PSD_2x2_coefficient_bases':[serial(a) for a in (mM,AM,alM)],
            'scalar_sign_rows':[[1,0,0],[-1,1,0],[-1,0,1]],'PSD_block_orders':[2,1,1,1],'total_PSD_order':5,
            'Bell_coordinate_order':['1','A0','A1','B0','B1','E00','E01','E10','E11'],
            'behaviour_target_3x9':serial(S.T*Ebell.T),
            'behaviour_dual':'Y>=0, z_m,z_D,z_beta>=0; b_m=4Y01+z_m-z_D-z_beta, b_A=Y00+z_D, b_alpha=2Y11+z_beta',
            'scope':'Exact for this three-dimensional source slice and its source dual, not the complete Q222 set.'}
    args.output.with_name('source_and_behaviour_pencil.json').write_text(json.dumps(pencil,indent=2)+'\n')
    out={'status':'EXACT GLOBAL C1 SLICE AND QUANTIFIED FIVE-PARAMETER TUBE; UNIVERSAL H1/H2 OPEN',
         'interval_r':['1/3','1/2'],'constant_Gram_orders':[12,8,2,2,6,4,1,1],'constant_Gram_lower_bound':'1/4',
         'test_library_size':18,'source_slice_PSD_order':5,'source_slice_criterion':'A>=m>=0, alpha>=m, A*alpha>=2m^2',
         'contact_map_shape':[42,45],'contact_map_rank':36,'pivot_minor':'-2^24','correction_bounds':['40','17/2'],
         'uniform_tube_coefficient_radius':'eta^2/246960','eta_range':'0<eta<1/12','prototype_r_range':'[1/3+eta,1/2-eta]',
         'integer_Gram_lower_bound':'eta^2/441','gamma_lower_bound':'1','delta_lower_bound':'1/9',
         'example_eta':str(eta),'example_r':str(r0),'example_radius':str(radius),'example_distance':str(distance),
         'off_family_native_parameters':list(map(str,theta)),'off_family_Gram_upper_vector':list(map(str,solved)),
         'off_family_integer_min_pivots':list(map(lambda p:str(min(p)),hp)),
         'constant_shifted_Gram_min_pivots':[str(min(p)) for p in pivots],
         'independent_review':False,'scope':'Universal admission is proved on the displayed fixed slice and explicit coefficient tube only. No source-wide H1 or H2 completion.'}
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    args.output.with_name('contact_map.json').write_text(json.dumps({'coefficient_keys':keys,'pivot_rows':pr,'pivot_columns':pc,'right_inverse':serial(E)},indent=2)+'\n')
    print('ALL EXACT GLOBAL SLICE AND TUBE CHECKS PASSED',flush=True)

if __name__=='__main__':main()
