#!/usr/bin/env python3
"""Fixed-map opening at the C2 n=0 / C1 beta=0 and C1 D=0 walls.

Symbolic identities and rational rank/right-inverse checks only. Conditional
on a wall Gram PD on all four larger contact complements. This is not a
proof of whole-wall PSD selection or universal strict-family admission.
"""
import argparse,itertools,json
from pathlib import Path
import sympy as s

BASE=Path(__file__).resolve().parent
U,V,W,X,Y,Z=s.symbols('U V W X Y Z');PSI=[W,X,Y,Z]
uc,vc,wc,xc,yc,zc=s.symbols('uc vc wc xc yc zc',nonzero=True)
d,n,m,C,A,D,al,beta=s.symbols('d n m C A D al beta',nonnegative=True)

def coefficients(expr):
    out={}
    for e,c in s.Poly(s.expand(expr),U,V,W,X,Y,Z).terms():
        states=[i for i,count in enumerate(e[2:]) for _ in range(count)]
        out[(*e[:2],*states)]=c
    return out

def wall_map(kind):
    features=[(i,j,l) for i in range(3) for j in range(3) for l in range(4)
              if (l!=0 or i+j>=2) and (l!=1 or j>=i)
              and (l!=3 or (i+j<=2 if kind=='n_beta' else i>=j))]
    blocks=[[f for f in features if (f[0]%2,f[1]%2)==p] for p in itertools.product((0,1),repeat=2)]
    columns=[];positions=[]
    for bi,b in enumerate(blocks):
        pivot=b[0]
        for i in range(1,len(b)):
            for j in range(i,len(b)):
                col={}
                for a,b0,c in ((b[i],b[j],1),(b[i],pivot,-1),(pivot,b[j],-1),(pivot,pivot,1)):
                    key=((a[0]+b0[0])//2,(a[1]+b0[1])//2,min(a[2],b0[2]),max(a[2],b0[2]))
                    col[key]=col.get(key,0)+c*(1 if i==j else 2)
                columns.append({key:c for key,c in col.items() if c});positions.append((bi,i-1,j-1))
    keys=sorted(set().union(*(set(c) for c in columns)));ki={k:i for i,k in enumerate(keys)}
    M=s.zeros(len(keys),len(columns))
    for j,col in enumerate(columns):
        for k,c in col.items():M[ki[k],j]=c
    ct={U:1,V:1,W:1,X:1,Y:1,Z:1}
    J=s.Matrix([[s.diff(U**i*V**j*PSI[l]*PSI[k],z).subs(ct) for i,j,l,k in keys]
                for z in (U,V,W,X,Y,Z)])
    pc=list(M.rref()[1]);pr=list(M[:,pc].T.rref()[1]);minor=M.extract(pr,pc)
    assert M.shape==(72,104) and len(pc)==66 and abs(minor.det())==2**54
    assert J.rank()==6 and J*M==s.zeros(6,104)
    E=s.zeros(104,72);inv=minor.inv()
    for i,p in enumerate(pc):
        for j,row in enumerate(pr):E[p,row]=inv[i,j]
    assert M*E*M==M and E[78:,:]==s.zeros(26,72)
    norm=max(sum(abs(v) for v in E[i,:]) for i in range(78))
    return features,blocks,keys,M,J,E,positions,pc,pr,minor.det(),norm

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=BASE/'exact_results.json');args=ap.parse_args()
    u,v=s.symbols('u v');du,dv=1+u*u,1+v*v
    F=(4*n*du*dv*X**2+4*m*(v*W-X+u*v*Y-u*Z)**2+4*C*dv*(u*W-Y)**2
       +4*A*du*dv*Y**2+(4*D+2*d)*(u*W+u*v*X-Y-v*Z)**2
       +(4*al+2*d)*du*(v*Y-Z)**2+(4*beta+2*d)*dv*(X+u*Z)**2-2*d*du*dv*Z**2)
    full=s.expand(F.subs({u:uc*U,v:vc*V,W:wc*W,X:xc*X,Y:yc*Y,Z:zc*Z},simultaneous=True))
    cn=4*n*uc**2*xc**2*((U*X-V*(X+W)/2)**2+U*V*(X-W)**2)
    cb=4*beta*uc**2*vc**2*zc**2*((U*V*Z-(U+V)*Z/2)**2+U**2*V*(Z-W)**2+U*V**2*(Z-W)**2)
    cd=4*D*vc**2*zc**2*((V*Z-U*(Z+W)/2)**2+U*V*(Z-W)**2)
    ct={U:1,V:1,W:1,X:1,Y:1,Z:1}
    for carrier in (cn,cb,cd):
        assert carrier.subs(ct)==0
        assert all(s.expand(s.diff(carrier,z).subs(ct))==0 for z in (U,V,W,X,Y,Z))
    assert coefficients(cn)[(2,0,1,1)]==4*n*uc**2*xc**2
    assert s.expand(coefficients(cb)[(2,2,3,3)]-4*beta*uc**2*vc**2*zc**2)==0
    assert s.expand(coefficients(cd)[(0,2,3,3)]-4*D*vc**2*zc**2)==0
    output={};maps={}
    for kind,carrier in [('n_beta',cn+cb),('D',cn+cd)]:
        feats,blocks,keys,M,J,E,pos,pc,pr,det,norm=wall_map(kind)
        co=coefficients(full-carrier);assert not set(co)-set(keys)
        rhs=s.Matrix([co.get(k,0) for k in keys])
        jets=s.Matrix([s.diff(full,z).subs(ct) for z in (U,V,W,X,Y,Z)])
        assert all(s.expand(v)==0 for v in J*rhs-jets)
        added={(2,0,1)}|({(2,2,3),(2,1,3),(1,2,3)} if kind=='n_beta' else {(0,2,3)})
        joint=set(feats)|added
        weighted=joint|{(i,j+2,l) for i,j,l in joint}
        assert all(0<=3*i<=6 and 0<=3*j<=12 for i,j,l in weighted)
        output[kind]={'wall_feature_count':27,'block_orders':list(map(len,blocks)),'map_shape':[72,104],'rank':66,
                      'minor':str(det),'right_inverse_integer_row_bound':str(norm),
                      'joint_opening_feature_count':len(joint),'weighted_joint_feature_count':len(weighted),
                      'scope':'All residual coefficients are reachable from the six stationary jets. PSD is proved locally only from a wall Gram PD on all four contact complements.'}
        maps[kind]={'features':feats,'blocks':blocks,'coefficient_keys':keys,'upper_positions':pos,
                    'pivot_rows':pr,'pivot_columns':pc,'right_inverse':[[str(v) for v in row] for row in E.tolist()]}
        print('PASS',kind,'fixed map72x104 rank66; imM=kerJ; complete symbolic residual support and contact jets',flush=True)
    # Individual openings need fewer extra rows than the joint construction.
    output['individual_feature_counts']={'C2_n':28,'C1_beta':30,'C1_D':28}
    output['status']='EXACT CONDITIONAL LOCAL ZERO-WEIGHT OPENINGS; UNIVERSAL H1/H2 OPEN'
    output['independent_review']=False
    args.output.write_text(json.dumps(output,indent=2)+'\n')
    args.output.with_name('fixed_maps.json').write_text(json.dumps(maps,indent=2)+'\n')
    print('ALL EXACT GENERAL ZERO-WEIGHT OPENING CHECKS PASSED',flush=True)

if __name__=='__main__':main()
