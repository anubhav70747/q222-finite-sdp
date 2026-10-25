"""Exact native C2 half-power exclusion, independently reconstructed from SG2."""
from fractions import Fraction as F
from pathlib import Path
import json

a,b,c,t,r=F(13,600),F(101,1000),F(47,15),F(27,35),F(13,154)
z=c*c;Ac=1+z*r;h=b-a*a*t
g=(1+z)*b-(1-r+a*a*z*t);j=1+z*r-t
Delta=((1+z)*z*b*b+2*z*(1+z)*b*r+(1+z)*(1+a*a*z)*b
       +z*z*(1-a*a)*r*r+2*z*(1+a*a*z)*r+1+a*a*z)
N=(1+z)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z)*(1-r)*t
source=[F(1),F(0),((b+r)*Ac*(Ac+a*a*z*t)-r*Delta)/(2*r*Delta),
 (1+a*a*z)*(1-r)*Ac*(Ac+a*a*z*t)/(2*g*Delta),
 Ac*N/(2*a*a*t*Delta),a*a*Ac*(b+r)*(Ac+a*a*z*t)/(2*h*Delta),
 (1+z)*(b+r)*(b-a*a*z*r)*(Ac+a*a*z*t)/(2*a*a*j*Delta),
 z*z*(b+r)*N/(2*(1-z*h)*Delta),F(0)]
d,lam,n,m,C,A,D,al,be=source
uc,vc=c,a*c
psi=[b/a,-c*r,c*h/a,F(1)]
gaps=[h,g,j,1-z*h,b-a*a*z*r]
claimed_gaps=[F(140893,1400000),F(601024007,3465000000),
 F(36637,34650),F(3767363,315000000),F(1255020827,12474000000)]
checks={'all_native_gaps_exact':gaps==claimed_gaps,
 'strict_native_domain':all(q>0 for q in [a,b,c,t,r,1-r,Delta,N]+gaps),
 'all_six_original_retained_events_positive':all(q>0 for q in [n,m,C,A,D,al]),
 'all_four_contact_amplitudes_nonzero':all(q!=0 for q in psi)}

# key=(u exponent,v exponent,first state index,second state index).
# Coefficients are those of the full quadratic polynomial, not the state matrix.
P={}
def mon(key,value):
    P[key]=P.get(key,F(0))+value
def packet(weight,terms,shifts=((0,0),)):
    for su,sv in shifts:
        for p,q,l,coef in terms:
            for pp,qq,k,co in terms:
                mon((p+pp+su,q+qq+sv,min(l,k),max(l,k)),weight*coef*co)
dv=((0,0),(0,2));du=((0,0),(2,0));both=((0,0),(2,0),(0,2),(2,2))
packet(4*lam,[(0,0,0,1),(0,1,1,1)],du)
packet(4*n,[(0,0,1,1)],both)
packet(4*m,[(0,1,0,1),(0,0,1,-1),(1,1,2,1),(1,0,3,-1)])
packet(4*C,[(1,0,0,1),(0,0,2,-1)],dv)
packet(4*A,[(0,0,2,1)],both)
packet(4*D+2*d,[(1,0,0,1),(1,1,1,1),(0,0,2,-1),(0,1,3,-1)])
packet(4*al+2*d,[(0,1,2,1),(0,0,3,-1)],du)
packet(4*be+2*d,[(0,0,1,1),(1,0,3,1)],dv)
packet(-2*d,[(0,0,3,1)],both)
P={key:value for key,value in P.items() if value}
def evaluate(poly,uu,vv,ps):
    return sum(value*uu**p*vv**q*ps[l]*ps[k] for (p,q,l,k),value in poly.items())
checks['complete_parent_value_zero']=evaluate(P,uc,vc,psi)==0
grad=[]
for angle in [0,1]:
    total=F(0)
    for (p,q,l,k),value in P.items():
        exponent=(p,q)[angle]
        if exponent:
            total+=exponent*value*uc**(p-(angle==0))*vc**(q-(angle==1))*psi[l]*psi[k]
    grad.append(total)
for ell in range(4):
    total=F(0)
    for (p,q,l,k),value in P.items():
        if l==ell:total+=value*uc**p*vc**q*psi[k]
        if k==ell:total+=value*uc**p*vc**q*psi[l]
    grad.append(total)
checks['all_six_stationary_equations']=all(q==0 for q in grad)
checks['WW_Newton_diagonal']={(p,q):v for (p,q,l,k),v in P.items() if l==k==0}=={
 (0,2):4*m,(2,0):4*C+4*D+2*d,(2,2):4*C}
checks['ZZ_Newton_diagonal']={(p,q):v for (p,q,l,k),v in P.items() if l==k==3}=={
 (0,0):4*al,(2,0):4*al+2*d+4*m,(0,2):4*D}

rows=[
[161019792,0,0,0,-121999126,-51597381,-92880749,198288996,63902153,19654676],
[187856424,-46850180,181090964,93807406,104512504,-51884026,-170947993,104746412,60911824,56523337],
[198410875,-202408890,134972214,153970197,2021367851,-131764329,-929567199,102173701,98542620,600959153],
[187856424,-20723372,228808703,94672500,-8876925,-32301218,-55239329,117094417,62928564,26333033],
[211714865,-20971443,135517649,91673664,-66738426,-13877461,-150724699,104061701,67646572,392762321],
[165685625,-100092162,132934774,129085407,872170933,-110726908,-814152246,108962482,438514151,328748742],
[259329330,-1150608,147851759,93736392,966590370,-12339807,-94532325,112427108,80542993,535327210],
[166289097,16862629,134869776,98549977,120025352,-54544777,-527190917,121915289,582567055,328748742],
[163696817,-79616153,139866971,0,703323270,-491353386,0,629878509,0,493123113]]
statepairs=[(l,k) for l in range(4) for k in range(l,4)]
Y={(p,q,l,k):rows[3*p+q][i] for p in range(3) for q in range(3) for i,(l,k) in enumerate(statepairs)}
def moment(e,f):
    i,j,l=e;ii,jj,k=f
    assert (i+ii)%2==0 and (j+jj)%2==0
    return F(Y[((i+ii)//2,(j+jj)//2,min(l,k),max(l,k))])
blocks=[]
for parity in [(0,0),(0,1),(1,0),(1,1)]:
    # Explicit order: W,X,Y,Z outermost, followed by i,j lexicographic.
    blocks.append([(i,j,l) for l in range(4) for i in range(3) for j in range(3)
      if (i%2,j%2)==parity and (l!=0 or i+j>=2) and (l!=3 or i+j<=2)])
checks['support_block_sizes']=[len(blk) for blk in blocks]==[14,6,6,4]
dual=[]
for blk in blocks:
    p0=blk[0];free=blk[1:]
    dual.append([[moment(e,f)-moment(e,p0)-moment(p0,f)+moment(p0,p0)
                  for f in free] for e in free])
checks['compressed_block_sizes']=[len(Q) for Q in dual]==[13,5,5,3]
def ldl(Q):
    nn=len(Q);L=[[F(i==j) for j in range(nn)] for i in range(nn)];diag=[]
    for j in range(nn):
        pivot=Q[j][j]-sum(L[j][k]*L[j][k]*diag[k] for k in range(j))
        if pivot<=0:raise AssertionError(('nonpositive pivot',j,pivot))
        diag.append(pivot)
        for i in range(j+1,nn):
            L[i][j]=(Q[i][j]-sum(L[i][k]*L[j][k]*diag[k] for k in range(j)))/pivot
    assert all(Q[i][j]==sum(L[i][k]*diag[k]*L[j][k] for k in range(nn))
               for i in range(nn) for j in range(nn))
    return diag
pivots=[ldl(Q) for Q in dual]
proved_bounds=[89,100,360,2000000]
checks['all_exact_LDL_reconstructions_and_positive_pivots']=True
checks['conservative_exact_LDL_pivot_bounds']=all(min(ps)>bound for ps,bound in zip(pivots,proved_bounds))

# Invertible contact normalization: original u=uc U, v=vc V,
# and original state component psi_l=psi_contact_l * phi_l.
normP={(p,q,l,k):val*uc**p*vc**q*psi[l]*psi[k] for (p,q,l,k),val in P.items()}
pairing=sum(val*Y[key] for key,val in normP.items())
expected=-F(31917121340091771507320448459210586965990611717744093996544823253,
            1129383645613305483559298619423160947865086012060000000000)
checks['exact_negative_full_parent_pairing']=pairing==expected and pairing<0
result={'checks':checks,'all_pass':all(checks.values()),
 'native_tuple':dict(zip(['a','b','c','t','r'],map(str,[a,b,c,t,r]))),
 'source':dict(zip(['d','lambda','n','m','C','A','D','alpha','beta'],map(str,source))),
 'contact':{'u':str(uc),'v':str(vc),'state_Z_equals_1':list(map(str,psi))},
 'native_gaps':dict(zip(['h','g','j','1-z*h','b-a^2*z*r'],map(str,gaps))),
 'feature_order':'state W,X,Y,Z outermost; then i,j lexicographic; first retained feature is pivot',
 'compressed_block_min_LDL_pivots':list(map(lambda q:str(min(q)),pivots)),
 'proved_LDL_lower_bounds':proved_bounds,
 'all_LDL_pivots':[[str(v) for v in q] for q in pivots],
 'negative_pairing':str(pairing),
 'scope':'Full36-feature half-power cone exclusion by necessary support and contact compression. Y is a compressed-face dual; no sixth-root obstruction.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'checks':checks,'all_pass':result['all_pass'],
 'min_LDL_pivots':result['compressed_block_min_LDL_pivots'],'negative_pairing':str(pairing)},indent=2))
if not result['all_pass']:raise SystemExit(1)
