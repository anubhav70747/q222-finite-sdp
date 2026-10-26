"""Exact audit of the sixth-root fractional-channel necessity argument.
Only the 90 integer moments are read from the published half-power checker.
No optimization, floating point, or execution of that checker is used.
"""
from fractions import Fraction as F
from pathlib import Path
import ast,json
ROOT=Path(__file__).resolve().parents[1]
oldpath=ROOT/'github-q222-finite-sdp/independent-review/checks/c2_half_power_obstruction.py'
tree=ast.parse(oldpath.read_text(encoding='utf-8'))
rows=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(q,ast.Name) and q.id=='rows' for q in n.targets))
pairs=[(l,m) for l in range(4) for m in range(l,4)]
Y={(p,q,l,m):F(rows[3*p+q][i]) for p in range(3) for q in range(3) for i,(l,m) in enumerate(pairs)}
T=24000000000000
B=24002021367851
checks={}

def moment(e,f,add=True):
    i,j,l=e;ii,jj,m=f
    if (i+ii)%6 or (j+jj)%6:return F(0)
    return Y[((i+ii)//6,(j+jj)//6,min(l,m),max(l,m))]+(T if add else 0)

def ldl(Q):
    n=len(Q);L=[[F(i==j) for j in range(n)] for i in range(n)];D=[]
    for j in range(n):
        d=Q[j][j]-sum(L[j][k]**2*D[k] for k in range(j))
        assert d>0
        D.append(d)
        for i in range(j+1,n):L[i][j]=(Q[i][j]-sum(L[i][k]*L[j][k]*D[k] for k in range(j)))/d
    assert all(Q[i][j]==sum(L[i][k]*D[k]*L[j][k] for k in range(n)) for i in range(n) for j in range(n))
    return L,D

def quadratic_inverse(Q,b):
    L,D=ldl(Q);z=[]
    for i in range(len(b)):z.append(b[i]-sum(L[i][k]*z[k] for k in range(i)))
    return sum(z[i]**2/D[i] for i in range(len(b)))

parities=[(0,0),(0,1),(1,0),(1,1)]
allblocks=[];oldblocks=[];newblocks=[]
for parity in parities:
    full=[(i,j,l) for l in range(4) for i in range(7) for j in range(7) if (i%2,j%2)==parity]
    keep=[e for e in full if (e[2]!=0 or e[0]+e[1]>=6) and (e[2]!=3 or e[0]+e[1]<=6)]
    old=[e for e in keep if e[0]%3==0 and e[1]%3==0]
    new=[e for e in keep if e not in old]
    allblocks.append(keep);oldblocks.append(old);newblocks.append(new)
checks['unpruned_orders']=[len([1 for l in range(4) for i in range(7) for j in range(7) if (i%2,j%2)==pa]) for pa in parities]==[64,48,48,36]
checks['Newton_retained_orders']=[len(b) for b in allblocks]==[52,36,36,30]
checks['contact_reduced_orders']=[len(b)-1 for b in allblocks]==[51,35,35,29]
checks['old_new_support_counts']=[len(b) for b in oldblocks]==[14,6,6,4] and sum(map(len,newblocks))==124
old=sum(oldblocks,[]);new=sum(newblocks,[])
checks['exact_old_new_adjoint_zero']=all(moment(e,f)==0 for e in old for f in new)
checks['every_new_diagonal_zero']=all(moment(e,e)==0 for e in new)
checks['all_nonzero_new_pairs_are_complementary']=all(((e[0]+f[0])%6==0 and (e[1]+f[1])%6==0 and any((e[k],f[k]) in [(1,5),(5,1),(2,4),(4,2)] for k in [0,1])) for e in new for f in new if moment(e,f))
checks['maximum_entry_exact']=max(abs(moment(e,f)) for e in old+new for f in old+new)==B
schur=[];mins=[]
for block in oldblocks:
    p0=block[0];rest=block[1:]
    aa=moment(p0,p0);bb=[moment(p0,e)-aa for e in rest]
    Q=[[moment(e,f)-moment(e,p0)-moment(p0,f)+aa for f in rest] for e in rest]
    L,diag=ldl(Q);mins.append(min(diag))
    schur.append(aa-quadratic_inverse(Q,bb))
checks['compressed_old_blocks_positive_definite']=True
schur_bounds=[400000000000,15000000000000,19000000000000,23000000000000]
checks['claimed_exact_Schur_gaps']=all(x>y for x,y in zip(schur,schur_bounds))

# Reconstruct the independent native source and full SG2 parent.
a,b,c,t,r=F(13,600),F(101,1000),F(47,15),F(27,35),F(13,154)
z=c*c;Ac=1+z*r;h=b-a*a*t;g=(1+z)*b-(1-r+a*a*z*t);j=1+z*r-t
Delta=((1+z)*z*b*b+2*z*(1+z)*b*r+(1+z)*(1+a*a*z)*b+z*z*(1-a*a)*r*r+2*z*(1+a*a*z)*r+1+a*a*z)
N=(1+z)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z)*(1-r)*t
source=[F(1),F(0),((b+r)*Ac*(Ac+a*a*z*t)-r*Delta)/(2*r*Delta),(1+a*a*z)*(1-r)*Ac*(Ac+a*a*z*t)/(2*g*Delta),Ac*N/(2*a*a*t*Delta),a*a*Ac*(b+r)*(Ac+a*a*z*t)/(2*h*Delta),(1+z)*(b+r)*(b-a*a*z*r)*(Ac+a*a*z*t)/(2*a*a*j*Delta),z*z*(b+r)*N/(2*(1-z*h)*Delta),F(0)]
d,la,n,m,C,A,D,al,be=source
uc,vc=c,a*c;psi=[b/a,-c*r,c*h/a,F(1)]
P={}
def packet(weight,terms,shifts=((0,0),)):
    for su,sv in shifts:
        for p,q,l,x in terms:
            for pp,qq,k,y in terms:
                key=(p+pp+su,q+qq+sv,min(l,k),max(l,k))
                P[key]=P.get(key,F(0))+weight*x*y
DU=((0,0),(2,0));DV=((0,0),(0,2));BOTH=((0,0),(2,0),(0,2),(2,2))
packet(4*la,[(0,0,0,1),(0,1,1,1)],DU)
packet(4*n,[(0,0,1,1)],BOTH)
packet(4*m,[(0,1,0,1),(0,0,1,-1),(1,1,2,1),(1,0,3,-1)])
packet(4*C,[(1,0,0,1),(0,0,2,-1)],DV)
packet(4*A,[(0,0,2,1)],BOTH)
packet(4*D+2*d,[(1,0,0,1),(1,1,1,1),(0,0,2,-1),(0,1,3,-1)])
packet(4*al+2*d,[(0,1,2,1),(0,0,3,-1)],DU)
packet(4*be+2*d,[(0,0,1,1),(1,0,3,1)],DV)
packet(-2*d,[(0,0,3,1)],BOTH)
P={key:value for key,value in P.items() if value}
normP={(p,q,l,k):val*uc**p*vc**q*psi[l]*psi[k] for (p,q,l,k),val in P.items()}
checks['independent_native_source_strict']=all(q>0 for q in [a,b,c,t,r,1-r,h,g,j,1-z*h,b-a*a*z*r,n,m,C,A,D,al])
checks['contact_value_exact_zero']=sum(normP.values())==0
checks['native_WW_Newton_support']=set((p,q) for p,q,l,k in P if l==k==0)=={(0,2),(2,0),(2,2)}
checks['native_ZZ_Newton_support']=set((p,q) for p,q,l,k in P if l==k==3)=={(0,0),(2,0),(0,2)}
pairing=sum(val*(Y[key]+T) for key,val in normP.items())
delta=F(31917121340091771507320448459210586965990611717744093996544823253,1129383645613305483559298619423160947865086012060000000000)
checks['same_exact_negative_pairing']=pairing==-delta
tr_bound=delta/(124*B);off_bound=delta/(2*B)
checks['trace_bound_claim']=tr_bound>F(949,10**11)
checks['offdiagonal_bound_claim']=off_bound>F(588,10**9)
# Entry bound yields ||M_N||_op <= maximum absolute row sum <=123B <124B.
checks['new_absolute_row_sum_bound']=max(sum(abs(moment(e,f)) for f in new) for e in new)<=123*B
out={'checks':checks,'all_pass':all(checks.values()),'T':T,'B':B,'old_orders':list(map(len,oldblocks)),'new_orders':list(map(len,newblocks)),'Schur_gaps':list(map(str,schur)),'strict_Schur_lower_bounds':schur_bounds,'compressed_old_min_LDL':list(map(str,mins)),'delta':str(delta),'trace_lower_bound':str(tr_bound),'offdiagonal_absolute_sum_lower_bound':str(off_bound),'old_blocks_features':oldblocks,'new_blocks_features':newblocks,'source':list(map(str,source)),'scope':'Conditional necessity for every admitting retained sixth-root Gram. No sixth-root admission or exclusion theorem.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'all_pass':out['all_pass'],'checks':checks,'Schur_gaps_approx_only':[float(q) for q in schur],'trace_bound_approx_only':float(tr_bound),'off_bound_approx_only':float(off_bound)},indent=2))
assert out['all_pass']
