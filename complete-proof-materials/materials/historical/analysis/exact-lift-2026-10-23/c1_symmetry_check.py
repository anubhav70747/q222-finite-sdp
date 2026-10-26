"""Finite exact check of the physical C1 stabilizer and native transport."""
from pathlib import Path
from itertools import product,permutations
import json
import sympy as s

H=((1,-1),(1,1))
names=['lambda','D','C','A','alpha','m','beta','n']
events=[(1,1,0,1),(-1,1,1,1),(-1,1,1,0),(-1,1,0,0),
        (-1,-1,0,1),(1,-1,1,1),(1,-1,1,0),(1,-1,0,0)]
preservers=[]
raw_edge_stabilizers=[]
for swap,pa,pb,aa,bb in product((False,True),permutations((0,1)),permutations((0,1)),product((-1,1),repeat=2),product((-1,1),repeat=2)):
    edge_image=[]
    for a,b,x0,y0 in (events[0],events[7]):
        edge_image.append((bb[y0]*b,aa[x0]*a,pb[y0],pa[x0]) if swap else (aa[x0]*a,bb[y0]*b,pa[x0],pb[y0]))
    if set(edge_image)=={events[0],events[7]}:
        raw_edge_stabilizers.append({'party_swap':swap,'A_input_permutation':pa,'B_input_permutation':pb,'A_signs':aa,'B_signs':bb})
    coeff=[[0,0],[0,0]]
    for i,j in product(range(2),repeat=2):
        oi,oj=(pb[j],pa[i]) if swap else (pa[i],pb[j])
        coeff[oi][oj]+=H[i][j]*aa[i]*bb[j]
    if coeff!=[list(row) for row in H]:
        continue
    perm=[]
    for a,b,x,y in events:
        old=(bb[y]*b,aa[x]*a,pb[y],pa[x]) if swap else (aa[x]*a,bb[y]*b,pa[x],pb[y])
        perm.append(events.index(old))
    preservers.append({'party_swap':swap,'A_input_permutation':pa,'B_input_permutation':pb,
                       'A_signs':aa,'B_signs':bb,'new_event_to_old':dict(zip(names,[names[i] for i in perm])),
                       'event_index_permutation':perm})
assert len(preservers)==16
assert len(raw_edge_stabilizers)==4
dihedral={tuple((eps*i+j)%8 for i in range(8)) for eps in(-1,1) for j in range(8)}
assert {tuple(v['event_index_permutation']) for v in preservers}==dihedral
stabilizers=[v for v in preservers if {v['event_index_permutation'][0],v['event_index_permutation'][7]}=={0,7}]
assert len(stabilizers)==2

checks={}
def ck(name,expr):
    rem=s.cancel(expr)
    if rem!=0:
        raise AssertionError(name+': '+str(rem))
    checks[name]=True

u,v,W,X,Y,Z=s.symbols('u v W X Y Z')
d,lam,n,m,C,A,D,al,be=s.symbols('d lam n m C A D al be')
def parent(cc,uu,vv,psi):
    dd,ll,nn,mm,cc0,aa,dd0,al0,bb=cc
    ww,xx,yy,zz=psi; Du=1+uu*uu;Dv=1+vv*vv
    return (4*ll*Du*(ww+vv*xx)**2+4*nn*Du*Dv*xx*xx
            +4*mm*(vv*ww-xx+uu*vv*yy-uu*zz)**2+4*cc0*Dv*(uu*ww-yy)**2
            +4*aa*Du*Dv*yy*yy+(4*dd0+2*dd)*(uu*ww+uu*vv*xx-yy-vv*zz)**2
            +(4*al0+2*dd)*Du*(vv*yy-zz)**2+(4*bb+2*dd)*Dv*(xx+uu*zz)**2-2*dd*Du*Dv*zz*zz)
psi1=(v*W-X,-W-v*X,Z-v*Y,Y+v*Z)
old=parent((d,lam,n,m,C,A,D,al,be),u,v,(W,X,Y,Z))
new=parent((d,n,lam,C,m,al,be,A,D),1/u,v,psi1)
ck('AS4 full nine-coefficient parent',s.expand(u*u*new-(1+v*v)*old))
stateM=s.Matrix([[v,-1,0,0],[-1,-v,0,0],[0,0,-v,1],[0,0,1,v]])
ck('state congruence norm',sum(stateM.T*stateM-(1+v*v)*s.eye(4)))
assert stateM.T*stateM==(1+v*v)*s.eye(4)
ck('state determinant',stateM.det()-(1+v*v)**2)

x,eta,zeta,r,h=s.symbols('x eta zeta r h')
rho=r+eta; G=h+zeta;V=1+x*eta*zeta
aa=(1-r)/eta;tt=(1-h)/zeta
M=1+x*r-tt;J=1+x*h-aa;E=1-x*eta*h;Q=1-x*zeta*r
trans={x:1/x,eta:x*zeta*rho/G,zeta:x*eta*G/rho,r:Q/G,h:E/rho}
tr=lambda expr:expr.subs(trans,simultaneous=True)
ck('native a transport',tr(aa)-M/(x*rho))
ck('native t transport',tr(tt)-J/(x*G))
ck('native M transport',tr(M)-aa*V/(x*G))
ck('native J transport',tr(J)-tt*V/(x*rho))
ck('native E transport',tr(E)-h*V/G)
ck('native Q transport',tr(Q)-r*V/rho)
ck('native rho transport',tr(rho)-V/G)
ck('native G transport',tr(G)-V/rho)
ck('actual Bob squared angle invariant',tr(x*eta*zeta)-x*eta*zeta)
for symbol in trans:
    ck('native involution '+str(symbol),tr(tr(symbol))-symbol)

# Independent d-normalized native coefficient formulas NC2--NC4.
K=1+x;Rs=1+x*r;Ts=1+x*h
P4=Ts*E+x*(1-r)*G
Pal=Rs*h*V+zeta*M*E
Pd=1+2*x*h+r*x*x*G+x*zeta+eta*x*(K*zeta+x*h*(1-h))
B5=K*V*rho-(eta*Ts-1+r)*Pd
Bbe=K*x*rho*G-(1-h)*Pd
co={'C':K*rho*E/(2*(1-r)*Pd),'A':eta*P4/(2*r*Pd),
    'm':K*x*rho*h*V/(2*M*zeta*Pd),'alpha':x*x*zeta*rho*Pal/(2*Q*Pd),
    'D':B5/(2*J*eta*Pd),'beta':Bbe/(2*(1-h)*Pd)}
pair={'C':'m','m':'C','A':'alpha','alpha':'A','D':'beta','beta':'D'}
for key,target in pair.items():
    ck('native coefficient '+key+' to '+target,tr(co[key])-co[target])

data={'physical_relabellings_total':128,'CHSH_preservers':16,'event_action':'i -> eps*i+j mod8 in cycle lambda,D,C,A,alpha,m,beta,n',
      'raw_zero_edge_stabilizers_without_fixing_CHSH':raw_edge_stabilizers,
      'preserving_relabellings':preservers,'C1_zero_edge_stabilizers':stabilizers,
      'exact_checks':len(checks),'checks':checks,
      'scope':'Physical signed-CHSH relabellings, fixed zero edge {lambda,n}; no claim for recanonicalizing a different zero edge on extra-zero boundary faces.'}
Path(__file__).with_name('c1_symmetry_check.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps({'CHSH_preservers':16,'C1_zero_edge_stabilizers':stabilizers,'exact_checks':len(checks)},indent=2))
