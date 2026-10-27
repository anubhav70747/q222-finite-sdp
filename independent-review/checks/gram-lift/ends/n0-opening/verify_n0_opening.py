#!/usr/bin/env python3
"""Exact C2 n=0 wall-to-strict opening. SymPy only; no optimizer or sweep.

Reconstructs a fixed contact-difference map, proves image = six-jet kernel,
repairs an archived wall proposal exactly, and certifies one strict source
by the resulting source-rational continuation. Universal H2 is not claimed.
"""
from pathlib import Path
import argparse,json
import sympy as s

BASE=Path(__file__).resolve().parent
U,V,W,X,Y,Z=s.symbols('U V W X Y Z')
PSI=(W,X,Y,Z)
# Exponent pairs refer to sqrt(U),sqrt(V). Every normalized feature is 1
# at U=V=1, PSI=(1,1,1,1). Preserve the original larger wall support.
BLOCKS=[[(0,2,0),(2,0,0),(2,2,0),(0,0,1),(0,2,1),(2,2,1),
         (0,0,2),(0,2,2),(2,0,2),(2,2,2),(0,0,3),(0,2,3),(2,0,3)],
        [(2,1,0),(0,1,1),(0,1,2),(2,1,2),(0,1,3)],
        [(1,2,0),(1,2,1),(1,0,2),(1,2,2),(1,0,3)],
        [(1,1,0),(1,1,1),(1,1,2),(1,1,3)]]
SIZES=[len(b)-1 for b in BLOCKS]
CT={U:1,V:1,W:1,X:1,Y:1,Z:1}

def product_column(block,i,j):
    pivot=block[0];ei,ej=block[i+1],block[j+1];out={}
    for x,y,c in [(ei,ej,1),(ei,pivot,-1),(pivot,ej,-1),(pivot,pivot,1)]:
        assert (x[0]+y[0])%2==0 and (x[1]+y[1])%2==0
        key=((x[0]+y[0])//2,(x[1]+y[1])//2,min(x[2],y[2]),max(x[2],y[2]))
        out[key]=out.get(key,0)+c
    return {k:v*(1 if i==j else 2) for k,v in out.items() if v}

def build_map():
    cols=[];positions=[]
    for bi,block in enumerate(BLOCKS):
        for i in range(len(block)-1):
            for j in range(i,len(block)-1):
                cols.append(product_column(block,i,j));positions.append((bi,i,j))
    keys=sorted(set().union(*(set(c) for c in cols)));ki={k:i for i,k in enumerate(keys)}
    M=s.zeros(len(keys),len(cols))
    for j,col in enumerate(cols):
        for k,c in col.items():M[ki[k],j]=c
    J=s.Matrix([[s.diff(U**i*V**j*PSI[l]*PSI[k],z).subs(CT)
                 for i,j,l,k in keys] for z in (U,V,W,X,Y,Z)])
    pc=list(M.rref()[1]);pr=list(M[:,pc].T.rref()[1]);minor=M.extract(pr,pc)
    assert M.shape==(72,104) and len(pc)==len(pr)==66
    assert minor.det()==2**54 and J*M==s.zeros(6,104)
    jp=list(J.rref()[1]);assert len(jp)==6 and J[:,jp].det()==2
    inv=minor.inv();E=s.zeros(104,72)
    for i,ci in enumerate(pc):
        for j,rj in enumerate(pr):E[ci,rj]=inv[i,j]
    assert M*E*M==M
    assert E[78:,:]==s.zeros(26,72)  # correction uses the integer block only
    return M,J,E,keys,positions,pc,pr,inv,jp

def coefficients(expr):
    out={}
    for exp,c in s.Poly(s.expand(expr),U,V,W,X,Y,Z).terms():
        states=[i for i,n in enumerate(exp[2:]) for _ in range(n)]
        assert len(states)==2
        out[(exp[0],exp[1],*states)]=c
    return out

def parent(u,v,psi,weights):
    w,x,y,z=psi;n,m,C,A,D,al=weights;du,dv=1+u*u,1+v*v
    return (4*n*du*dv*x*x+4*m*(v*w-x+u*v*y-u*z)**2
            +4*C*dv*(u*w-y)**2+4*A*du*dv*y*y
            +(4*D+2)*(u*w+u*v*x-y-v*z)**2
            +(4*al+2)*du*(v*y-z)**2+2*dv*(x+u*z)**2-2*du*dv*z*z)

def native(theta,wall=False):
    a,b,c,t,r=theta;z=c*c;B=1+z*r;h=b-a*a*t
    g=(1+z)*b-(1-r+a*a*z*t);j=B-t
    de=1+(1+z)*b+z*(1+z)*b*b+2*z*(1+z)*b*r+z*z*r*r+2*z*r+a*a*z*((1+z)*b+1+z*r*(2-r))
    N=(1+z)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z)*(1-r)*t
    Nn=(b+r)*B*(B+a*a*z*t)-r*de
    gaps=(a,b,c,t,r,1-r,h,g,j,1-z*h,b-a*a*z*r,de,N)
    assert all(v>0 for v in gaps)
    weights=[Nn/(2*r*de),(1+a*a*z)*(1-r)*B*(B+a*a*z*t)/(2*g*de),
             B*N/(2*a*a*t*de),a*a*B*(b+r)*(B+a*a*z*t)/(2*h*de),
             (1+z)*(b+r)*(b-a*a*z*r)*(B+a*a*z*t)/(2*a*a*j*de),
             z*z*(b+r)*N/(2*(1-z*h)*de)]
    assert weights[0]==0 if wall else weights[0]>0
    assert all(v>0 for v in weights[1:])
    return list(map(s.cancel,weights)),[b/a,-c*r,c*h/a,s.Integer(1)]

def target(theta,keys,wall=False):
    a,b,c,t,r=theta;weights,contact=native(theta,wall);n=weights[0]
    q=a/(1+a)
    full=parent(c*U,a*c*V,[contact[i]*PSI[i] for i in range(4)],weights)
    carrier=4*n*c**4*r*r*((-U*X+q*V*X+(1-q)*V*W)**2+2*q*U*V*(W-X)**2)
    co=coefficients(full-carrier);assert not set(co)-set(keys)
    return s.Matrix([co.get(k,0) for k in keys]),contact,weights

def blocks(vector,positions):
    out=[s.zeros(n) for n in SIZES]
    for val,(bi,i,j) in zip(vector,positions):out[bi][i,j]=out[bi][j,i]=val
    return out

def vector(block_matrices,positions):
    return s.Matrix([block_matrices[b][i,j] for b,i,j in positions])

def ldl_pivots(matrix):
    m=matrix.copy();piv=[]
    for k in range(m.rows):
        p=m[k,k]
        if p<=0:return None
        piv.append(p)
        for i in range(k+1,m.rows):
            for j in range(i,m.rows):
                m[i,j]=m[j,i]=s.cancel(m[i,j]-m[i,k]*m[k,j]/p)
    return piv

def normalized_proposal(pt,proposal,contact):
    c,a,b,r=map(s.Rational,pt);ys=[];out=[];cursor=0
    for block in BLOCKS:
        y=[s.sqrt(c)**i*s.sqrt(a*c)**j*contact[l] for i,j,l in block]
        assert all(v.is_Rational and v!=0 for v in y);ys.append(y)
        anchor=max(range(len(y)),key=lambda i:abs(y[i]));cols=[i for i in range(len(y)) if i!=anchor]
        N=s.zeros(len(y),len(cols))
        for j,i in enumerate(cols):N[i,j]=1;N[anchor,j]=-y[i]/y[anchor]
        S=s.zeros(len(cols))
        for i in range(S.rows):
            for j in range(i+1):S[i,j]=S[j,i]=s.Rational(proposal[cursor]);cursor+=1
        G=s.diag(*y)*N*S*N.T*s.diag(*y)
        assert G*s.ones(len(y),1)==s.zeros(len(y),1)
        out.append(G[1:,1:])
    assert cursor==len(proposal)==104
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=BASE/'exact_results.json');args=ap.parse_args()
    M,J,E,keys,pos,pc,pr,inv,jp=build_map()
    print('PASS fixed wall map: 72x104, rank66, minor2^54; six contact jets, minor2; image=ker J',flush=True)
    # All-source symbolic support and zero contact jets of the carrier. No native
    # coefficient identity is inferred from these formal weights alone.
    a,b,c,h,r,q=s.symbols('a b c h r q',positive=True)
    weights=s.symbols('n m C A D alpha',positive=True);n=weights[0]
    fn=parent(c*U,a*c*V,[b/a*W,-c*r*X,c*h/a*Y,Z],weights)
    ca=4*n*c**4*r*r*((-U*X+q*V*X+(1-q)*V*W)**2+2*q*U*V*(W-X)**2)
    cd=coefficients(ca);rd=coefficients(fn-ca)
    assert not set(rd)-set(keys)
    assert all(s.expand(s.diff(ca,z).subs(CT))==0 for z in (U,V,W,X,Y,Z))
    assert s.expand(rd.get((2,0,1,1),0))==0
    print('PASS all-source residual support; carrier six contact jets vanish; forced u^2 X^2 is removed exactly',flush=True)
    seed=json.loads((BASE/'seed_wall_proposal.json').read_text());pt,proposal=seed['point_c_a_b_r'],seed['proposal']
    cv,av,bv,rv=map(s.Rational,pt);t=s.symbols('t')
    z=cv*cv;B=1+z*rv
    de=1+(1+z)*bv+z*(1+z)*bv*bv+2*z*(1+z)*bv*rv+z*z*rv*rv+2*z*rv+av*av*z*((1+z)*bv+1+z*rv*(2-rv))
    num=(bv+rv)*B*(B+av*av*z*t)-rv*de;tv=s.solve(num,t)[0]
    theta0=[av,bv,cv,tv,rv];r0,contact0,weights0=target(theta0,keys,True)
    assert J*r0==s.zeros(6,1)
    p0=vector(normalized_proposal(pt,proposal,contact0),pos)
    x0=p0+E*(r0-M*p0);assert M*x0==r0
    Q0=blocks(x0,pos);piv0=[ldl_pivots(b) for b in Q0];assert all(p is not None for p in piv0)
    mus=[1/max(sum(abs(v) for v in row) for row in b.inv().tolist()) for b in Q0]
    Cb=[max([sum(abs(v) for v in E[j,:]) for j,(bi,i,k) in enumerate(pos) if bi==b]+[s.Integer(0)]) for b in range(4)]
    radius=min(mu/(2*size*bound) for mu,size,bound in zip(mus,SIZES,Cb) if bound>0)
    assert radius>0
    simple_radius=s.Rational(9,1000000)
    assert radius>simple_radius
    print('PASS archived wall seed repaired exactly; all four reduced blocks PD; positive rational continuation radius',flush=True)
    # This is a finite exact certificate, not a native-domain sweep. Decrease one
    # rational step until the explicit coefficient-radius hypothesis is met.
    found=None
    for sign in (-1,1):
        for power in range(1,90):
            theta1=theta0.copy();theta1[3]+=sign*s.Rational(1,2**power)
            try:r1,contact1,weights1=target(theta1,keys)
            except AssertionError:continue
            assert J*r1==s.zeros(6,1)
            delta=r1-r0
            if max(abs(v) for v in delta)>=simple_radius:continue
            x1=x0+E*delta;assert M*x1==r1
            Q1=blocks(x1,pos);piv1=[ldl_pivots(b) for b in Q1];assert all(p is not None for p in piv1)
            found=(theta1,weights1,x1,piv1,power,sign);break
        if found:break
    assert found is not None
    theta1,weights1,x1,piv1,power,sign=found
    old_features={tuple(f) for block in BLOCKS for f in block};opened=old_features|{(2,0,1)}
    weighted=opened|{(i,j+2,l) for i,j,l in opened}
    assert len(opened)==28 and len(weighted)==49
    assert all(0<=3*i<=6 and 0<=3*j<=12 for i,j,l in weighted)
    print('PASS strict finite-n certificate; all coefficients and all four residual blocks checked exactly',flush=True)
    print('PASS 28 fixed half-power features; multiplier embeds in49 fixed K364 features',flush=True)
    out={'status':'EXACT LOCAL WALL-TO-STRICT OPENING; UNIVERSAL H1/H2 REMAIN OPEN',
         'map_shape':[72,104],'rank':66,'pivot_minor':'2^54','jet_minor':'2',
         'unweighted_feature_count':28,'weighted_feature_union_count':49,
         'wall_theta_a_b_c_t_r':list(map(str,theta0)),'wall_weights_n_m_C_A_D_alpha':list(map(str,weights0)),
         'wall_reduced_upper_vector':list(map(str,x0)),
         'wall_reduced_min_pivots':[str(min(p)) for p in piv0],
         'rational_eigenvalue_lower_bounds':list(map(str,mus)),
         'right_inverse_block_row_bounds':list(map(str,Cb)),
         'coefficient_radius':str(radius),'simple_coefficient_radius':str(simple_radius),'strict_step_t':str(sign*s.Rational(1,2**power)),
         'strict_theta_a_b_c_t_r':list(map(str,theta1)),
         'strict_weights_n_m_C_A_D_alpha':list(map(str,weights1)),
         'strict_reduced_upper_vector':list(map(str,x1)),
         'strict_reduced_min_pivots':[str(min(p)) for p in piv1],
         'scope':'Every wall point with a PD Gram on all four larger wall contact blocks has an explicit source-rational finite-n neighborhood. Whole-wall PSD and arbitrary strict-source admission are unproved.',
         'independent_review':False}
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    mapout={'blocks':BLOCKS,'coefficient_keys':keys,'upper_positions':pos,'pivot_rows':pr,'pivot_columns':pc,
            'inverse_pivot_minor':[[str(v) for v in row] for row in inv.tolist()],'jet_pivot_columns':jp,
            'convention':'offdiagonal symmetric variables contribute twice; difference roots f_i-f_0'}
    args.output.with_name('fixed_map.json').write_text(json.dumps(mapout,indent=2)+'\n')
    print('ALL EXACT OPENING CHECKS PASSED; no universal source-admission claim',flush=True)

if __name__=='__main__':main()
